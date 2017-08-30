#include <iostream>
#include <limits>
#include "controller.hpp"
#include "pioneer.hpp"

using namespace std;
using namespace arma;
// using namespace robot::local;
using namespace robot::remote;

int main(int argc, char* argv[]) {
  std::string results_folder = argc == 2 ? argv[1] : "./";
  // Pioneer robot_1;
  // Pioneer robot_2;

  Pioneer robot_1("192.168.1.180", 1110);
  Pioneer robot_2("192.168.1.117", 1111);

  mat initial_pose_1;
  mat initial_pose_2;
  initial_pose_1 << 0.0 << endr << 0.0 << endr << deg2rad(0.0);
  initial_pose_2 << 1000.0 << endr << -2000.0 << endr << deg2rad(0.0);
  robot_1.set_pose(initial_pose_1);
  robot_2.set_pose(initial_pose_2);

  controller::Trajectory::Parameters parameters("parameters.yaml");
  parameters.trajectory = [&](double) { return vec({2000.0, 3000.0, 1000.0, 0.0}); };
  parameters.trajectory_speed = [&](double) { return vec(4, fill::zeros); };

  controller::Formation::Parameters formation_parameters("formation.yaml");

  controller::Formation controller(formation_parameters, parameters);

  mat telemetry_data;
  mat speed_data;
  mat trajectory_data;
  double error = std::numeric_limits<double>::max();
  Loop loop;
  do {
    mat telemetry_1 = robot_1.get_telemetry();
    mat telemetry_2 = robot_2.get_telemetry();  // returns a 5x1 matrix {x,y,psi,v,w} [mm,mm,rad,mm/s,rad/s]
    mat pose_1 = telemetry_1.rows(0, 2);
    mat pose_2 = telemetry_2.rows(0, 2);  // {x,y,psi}
    mat trajectory_pose;
    mat formation_pose;
    auto k_speed = controller.kinematics(pose_1, pose_2, formation_pose, trajectory_pose);
    auto d_speed = controller.dynamics_2(telemetry_1, telemetry_2, k_speed);
    robot_1.set_speed(std::get<0>(d_speed));
    robot_2.set_speed(std::get<1>(d_speed));

    telemetry_data = join_horiz(telemetry_data, join_vert(telemetry_1, telemetry_2));
    trajectory_data = join_horiz(trajectory_data, join_vert(formation_pose,trajectory_pose));
    speed_data = join_horiz(speed_data, join_vert(join_vert(std::get<0>(k_speed), std::get<1>(k_speed)),
                                                  join_vert(std::get<0>(d_speed), std::get<1>(d_speed))));

    error = norm(vec(trajectory_pose.rows(0, 1)) - vec(formation_pose.rows(0, 1)));
    loop.wait();
  } while (error > 100.0);

  // to load on Matlab: load('file.mat','-ASCII');
  telemetry_data.save(results_folder + "telemetry.mat", raw_ascii);
  trajectory_data.save(results_folder + "trajectory.mat", raw_ascii);
  speed_data.save(results_folder + "speed.mat", raw_ascii);

  robot_1.disconnect();
  robot_2.disconnect();
  return 0;
}