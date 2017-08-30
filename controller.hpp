#ifndef __CONTROLLER_HPP__
#define __CONTROLLER_HPP__

#include <yaml-cpp/yaml.h>
#include <armadillo>
#include <chrono>
#include <functional>
#include <list>

namespace YAML {
template <>
struct convert<arma::vec> {
  static Node encode(arma::vec const& vector) {
    Node node;
    for (const auto& value : vector) {
      node.push_back(value);
    }
    return node;
  }

  static bool decode(Node const& node, arma::vec& vector) {
    if (!node.IsSequence()) {
      return false;
    }
    std::vector<double> values;
    for (auto i = 0; i < node.size(); ++i) {
      values.push_back(node[i].as<double>());
    }
    vector = arma::vec(values);
    return true;
  }
};
}

namespace controller {

auto deg2rad = [](double deg) { return deg * (atan(1.0) / 45.0); };
auto rad2deg = [](double rad) { return rad * (45.0 / atan(1.0)); };

using namespace arma;
using namespace std::chrono;

struct Trajectory {
  struct Parameters {
    double center_offset;
    double lx;
    double ly;
    double kx;
    double ky;
    std::function<vec(double)> trajectory;
    std::function<vec(double)> trajectory_speed;
    vec dynamic_parameters;
    double lv;
    double lw;
    double kv;
    double kw;
    uint64_t period_ms;

    Parameters() {}
    Parameters(std::string const& filename) {
      YAML::Node node = YAML::LoadFile(filename);
      this->center_offset = node["center_offset"].as<double>();
      this->lx = node["lx"].as<double>();
      this->ly = node["ly"].as<double>();
      this->kx = node["kx"].as<double>();
      this->ky = node["ky"].as<double>();
      this->dynamic_parameters = node["teta_scale"].as<double>() * node["teta"].as<arma::vec>();
      this->lv = node["lv"].as<double>();
      this->lw = node["lw"].as<double>();
      this->kv = node["kv"].as<double>();
      this->kw = node["kw"].as<double>();
      this->period_ms = node["period_ms"].as<uint64_t>();
    }
  } parameters;

  bool first;

  Trajectory(Parameters const& parameters) {
    this->parameters = parameters;
    first = true;
  }

  void set_parameters(Parameters const& parameters) { this->parameters = parameters; }

  mat kinematics(mat const& current_pose) {  // [x,y,psi]
    mat trajectory_pose(0, 0, fill::zeros);
    return kinematics(current_pose, trajectory_pose);
  }

  mat kinematics(mat const& current_pose, mat& trajectory_pose) {  // [x,y,psi]
    static auto t0 = high_resolution_clock::now();
    auto t = duration_cast<milliseconds>(high_resolution_clock::now() - t0).count() / 1000.0;

    vec current_pose_vec(current_pose);  // [x,y,psi]'
    auto heading = current_pose_vec(2);
    auto a = this->parameters.center_offset;
    current_pose_vec = current_pose_vec + vec({a * cos(heading), a * sin(heading), 0.0});
    vec current_point = current_pose_vec.subvec(0, 1);  // [x,y]'
    vec trajectory_point = this->parameters.trajectory(t);
    trajectory_pose = mat(join_vert(trajectory_point, vec({0.0})));
    vec error = trajectory_point - current_point;

    vec trajectory_speed = this->parameters.trajectory_speed(t);

    mat invA = {{cos(heading), sin(heading)}, {-(1.0 / a) * sin(heading), (1.0 / a) * cos(heading)}};

    vec L = {this->parameters.lx, this->parameters.ly};
    vec K = {this->parameters.kx, this->parameters.ky};
    mat C = trajectory_speed + L % tanh((K / L) % error);

    return invA * C;  // [v,w]'
  }

  mat dynamics_1(mat const& current_telemetry, mat const& k_speed) {  // [x,y,psi,v,w]', [v_k,w_k]'
    auto speed = current_telemetry.rows(3, 4);                        // [v,w]'
    auto period_ms = this->parameters.period_ms;
    vec dk_speed = this->first ? vec(2, fill::zeros) : (1000.0 / period_ms) * (k_speed - this->last_k_speed);
    this->last_k_speed = k_speed;
    if (this->first)
      this->first = false;
    vec speed_tilde = k_speed - speed;                  // [v_k,w_k]' - [v,w]'
    vec k({this->parameters.kv, this->parameters.kw});  // [kv,kw]'
    vec sigma = dk_speed + k % speed_tilde;             // [sigma1, sigma2]'
                                                        // % -> element-wise multiplication
    double v = speed(0);
    double w = speed(1);
    auto teta = this->parameters.dynamic_parameters;
    mat G = {{sigma(0), 0.0, -w * w, v, 0.0, 0.0}, {0.0, sigma(1), 0.0, 0.0, v * w, w}};
    return G * teta;
  }

  mat dynamics_2(mat const& current_telemetry, mat const& k_speed,
                 bool adaptative = false) {     // [x,y,psi,v,w]', [v_k,w_k]'
    auto speed = current_telemetry.rows(3, 4);  // [v,w]'
    auto period_ms = this->parameters.period_ms;
    vec dk_speed = this->first ? vec(2, fill::zeros) : (1000.0 / period_ms) * (k_speed - this->last_k_speed);
    this->last_k_speed = k_speed;
    if (this->first)
      this->first = false;
    vec speed_tilde = k_speed - speed;                   // [v_k,w_k]' - [v,w]'
    vec L = {this->parameters.lv, this->parameters.lw};  // [lv,lw]'
    vec K = {this->parameters.kv, this->parameters.kw};  // [kv,kw]'
    vec sigma = dk_speed + L % tanh((K / L) % speed_tilde);

    double v = speed(0);
    double w = speed(1);
    double vd = k_speed(0);
    double wd = k_speed(1);

    mat G = {{sigma(0), 0.0, -wd * w, vd, 0.0, 0.0}, {0.0, sigma(1), vd * w - v * wd, 0.0, v * wd, wd}};
    if (adaptative) {
      mat gama = 0.005 * diagmat(vec({0.005, 0.05, 0.0005, 0.01, 0.005, 0.01}));
      mat GAMA = 0.005 * diagmat(vec({0.005, 0.05, 0.0005, 0.01, 0.005, 0.01}));
      vec teta_dot = gama * G.t() * (k_speed - speed) - gama * GAMA * this->parameters.dynamic_parameters;
      teta_dot.print("teta_dot");
      this->parameters.dynamic_parameters += teta_dot * (period_ms / 1000.0);
    }
    auto teta = this->parameters.dynamic_parameters;
    return G * teta;
  }

 private:
  vec last_k_speed;
};

struct Formation {
  struct Parameters {
    double center_offset;
    vec L;
    vec K;

    Parameters() {}
    Parameters(std::string const& filename) {
      YAML::Node node = YAML::LoadFile(filename);
      this->center_offset = node["center_offset"].as<double>();
      this->L = node["L"].as<vec>();
      this->L(3) = deg2rad(this->L(3));
      this->K = node["K"].as<vec>();
      this->K(3) = deg2rad(this->K(3));
    }
  } formation_parameters;
  Trajectory::Parameters trajectory_parameters;
  Trajectory controller_1;
  Trajectory controller_2;

  Formation(Formation::Parameters const& formation_parameters, Trajectory::Parameters const& trajectory_parameters)
      : controller_1(trajectory_parameters), controller_2(trajectory_parameters) {
    this->formation_parameters = formation_parameters;
    this->trajectory_parameters = trajectory_parameters;
  }

  std::pair<mat, mat> kinematics(mat const& current_pose_1, mat const& current_pose_2) {
    mat formation_pose(0, 0, fill::zeros);
    mat trajectory_pose(0, 0, fill::zeros);
    return kinematics(current_pose_1, current_pose_2, formation_pose, trajectory_pose);
  }

  std::pair<mat, mat> kinematics(mat const& current_pose_1, mat const& current_pose_2, mat& formation_pose,
                                 mat& trajectory_pose) {
    static auto t0 = high_resolution_clock::now();
    auto t = duration_cast<milliseconds>(high_resolution_clock::now() - t0).count() / 1000.0;

    vec x = join_vert(vec(current_pose_1.rows(0, 1)), vec(current_pose_2.rows(0, 1)));
    vec q = f(x);
    vec q_des = this->trajectory_parameters.trajectory(t);
    vec dq_des = this->trajectory_parameters.trajectory_speed(t);
    vec q_tilde = q_des - q;

    formation_pose = q;
    trajectory_pose = q_des;

    vec L = this->formation_parameters.L;
    vec K = this->formation_parameters.K;

    vec dq_ref = dq_des + L % tanh((K / L) % q_tilde);
    vec dx_ref = invJ(q) * dq_ref;
    double psi_1 = current_pose_1(2, 0);
    double psi_2 = current_pose_2(2, 0);
    vec v_ref = invK(psi_1, psi_2) * dx_ref;

    return std::make_pair(mat(v_ref.subvec(0, 1)), mat(v_ref.subvec(2, 3)));
  }

  std::pair<mat, mat> dynamics_1(mat const& current_telemetry_1, mat const& current_telemetry_2,
                                 std::pair<mat, mat> const& k_speeds) {
    return dynamics_1(current_telemetry_1, current_telemetry_2, std::get<0>(k_speeds), std::get<1>(k_speeds));
  }

  std::pair<mat, mat> dynamics_1(mat const& current_telemetry_1, mat const& current_telemetry_2, mat const& k_speed_1,
                                 mat const& k_speed_2) {
    return std::make_pair(this->controller_1.dynamics_1(current_telemetry_1, k_speed_1),
                          this->controller_2.dynamics_1(current_telemetry_2, k_speed_2));
  }

  std::pair<mat, mat> dynamics_2(mat const& current_telemetry_1, mat const& current_telemetry_2,
                                 std::pair<mat, mat> const& k_speeds) {
    return dynamics_2(current_telemetry_1, current_telemetry_2, std::get<0>(k_speeds), std::get<1>(k_speeds));
  }

  std::pair<mat, mat> dynamics_2(mat const& current_telemetry_1, mat const& current_telemetry_2, mat const& k_speed_1,
                                 mat const& k_speed_2) {
    return std::make_pair(this->controller_1.dynamics_2(current_telemetry_1, k_speed_1),
                          this->controller_2.dynamics_2(current_telemetry_2, k_speed_2));
  }

 private:
  vec f(vec const& x) {
    return vec({
        (x(0) + x(2)) / 2.0,                                                  // x_f = (x1+x2)/2
        (x(1) + x(3)) / 2.0,                                                  // y_f = (y1+y2)/2
        sqrt((x(2) - x(0)) * (x(2) - x(0)) + (x(3) - x(1)) * (x(3) - x(1))),  // rho_f = sqrt((x2-x1)^2 + (y2-y1)^2)
        atan2(x(3) - x(1), x(2) - x(0))                                       // alpha_f = atan((y2-y1)/(x2-x1))
    });
  }

  mat invJ(vec const& q) {
    mat invJ = zeros<mat>(4, 4);
    invJ(0, 0) = 1.0;
    invJ(1, 1) = 1.0;
    invJ(2, 0) = 1.0;
    invJ(3, 1) = 1.0;

    invJ(0, 2) = -cos(q(3)) / 2.0;         // J^-1(1,3)= -cos(alpha_f)/2
    invJ(0, 3) = q(2) * sin(q(3)) / 2.0;   // J^-1(1,4)= rho_f*sin(alpha_f)/2
    invJ(1, 2) = -sin(q(3)) / 2.0;         // J^-1(2,3)= -sin(alpha_f)/2
    invJ(1, 3) = -q(2) * cos(q(3)) / 2.0;  // J^-1(2,4)= -rho_f*cos(alpha_f)/2

    invJ(2, 2) = -invJ(0, 2);  // J^-1(3,3)= cos(alpha_f)/2
    invJ(2, 3) = -invJ(0, 3);  // J^-1(3,4)= -rho_f*sin(alpha_f)/2
    invJ(3, 2) = -invJ(1, 2);  // J^-1(4,3)= sin(alpha_f)/2
    invJ(3, 3) = -invJ(1, 3);  // J^-1(4,4)= rho_f*cos(alpha_f)/2

    return invJ;
  }

  mat invK(double const& psi1, double const& psi2) {
    mat invK = zeros<mat>(4, 4);
    double a = this->formation_parameters.center_offset;
    invK(0, 0) = cos(psi1);
    invK(0, 1) = sin(psi1);
    invK(1, 0) = -(1.0 / a) * sin(psi1);
    invK(1, 1) = (1.0 / a) * cos(psi1);

    invK(2, 2) = cos(psi2);
    invK(2, 3) = sin(psi2);
    invK(3, 2) = -(1.0 / a) * sin(psi2);
    invK(3, 3) = (1.0 / a) * cos(psi2);

    return invK;
  }
};

struct Path {
  struct Parameters {
    std::list<vec> path;
    double next_point_distance;
    double final_distance;

    Parameters() {}
    Parameters(std::string const& filename) {
      YAML::Node node = YAML::LoadFile(filename);
      for (auto i = 0; i < node["path"].size(); ++i) {
        vec point = node["path"][i].as<vec>();
        if (point.n_elem != 2) {
          throw std::runtime_error("Path point in " + filename + " must be have only two coordinates.");
        }
        this->path.push_back(point);
        node["path"][i].as<vec>().print();
      }
      this->next_point_distance = node["next_point_distance"].as<double>();
      this->final_distance = node["final_distance"].as<double>();
    }
  } path_parameters;

  std::unique_ptr<Trajectory> controller;
  Trajectory::Parameters controller_parameters;
  bool path_done;

  Path() {}
  Path(Path::Parameters const& path_parameters, Trajectory::Parameters const& controller_parameters) {
    this->path_parameters = path_parameters;
    this->controller_parameters = controller_parameters;
    this->controller_parameters.trajectory = [this](double) { return this->path_parameters.path.front(); };
    this->controller_parameters.trajectory_speed = [](double) { return vec(2, fill::zeros); };
    this->controller = std::move(std::unique_ptr<Trajectory>(new Trajectory(this->controller_parameters)));
    this->path_done = false;
  }

  mat kinematics(mat const& current_pose) {
    mat trajectory_pose(0, 0, fill::zeros);
    return kinematics(current_pose, trajectory_pose);
  }

  mat kinematics(mat const& current_pose, mat& trajectory_pose) {
    mat speeds = this->controller->kinematics(current_pose, trajectory_pose);
    auto dist_error = this->path_parameters.path.size() == 1 ? this->path_parameters.final_distance
                                                             : this->path_parameters.next_point_distance;
    if (error(current_pose, trajectory_pose) < dist_error) {
      this->path_parameters.path.pop_front();
      if (this->path_parameters.path.empty()) {
        this->path_done = true;
        return mat(2, 1, fill::zeros);
      }
    }
    return speeds;
  }

  bool on_path() { return !this->path_done; }

 private:
  double error(mat const& current_pose, mat const& desired_pose) {
    return norm(vec(desired_pose.rows(0, 1)) - vec(current_pose.rows(0, 1)));
  }
};

}  // ::controller

#endif  // __CONTROLLER_HPP__