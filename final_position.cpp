#include <iostream>
#include "controller.hpp"
#include "pioneer.hpp"

using namespace std;
using namespace arma;
using namespace robot::local;
// using namespace robot::remote;

auto error = [](mat const& current_pose, mat const& desired_pose) {
  return norm(vec(desired_pose.rows(0, 1)) - vec(current_pose.rows(0, 1)));
};

int main(int argc, char* argv[]) {
  std::string results_folder = argc == 2 ? argv[1] : "./";
  Pioneer robot;
  // Pioneer robot("192.168.1.180", 1111);

	mat initial_pose;
  initial_pose << 0.0 << endr << 0.0 << endr << deg2rad(0.0);
  robot.set_pose(initial_pose);

  controller::Trajectory::Parameters parameters("parameters_final_position.yaml");
  const double R = 1000.0;  // [mm]
  const double w = parameters.lx / R; // lx => max_vel
  parameters.trajectory = [&](double) { return vec({3000.0, 2500.0}); };
  parameters.trajectory_speed = [&](double) { return vec({0.0, 0.0}); };

  controller::Trajectory controller(parameters);

  mat telemetry_data;
	mat speed_data;
  mat trajectory_data;
  mat teta;
  
  mat pose, trajectory_pose;
  Loop loop;
  do {
    mat telemetry = robot.get_telemetry();  // returns a 5x1 matrix {x,y,psi,v,w} [mm,mm,rad,mm/s,rad/s]
    pose = telemetry.rows(0, 2);        // {x,y,psi}
    auto k_speed = controller.kinematics(pose, trajectory_pose);
		auto d_speed = controller.dynamics_2(telemetry, k_speed, false);
		robot.set_speed(d_speed);
		
    telemetry_data = join_horiz(telemetry_data, telemetry);
    trajectory_data = join_horiz(trajectory_data, trajectory_pose);
		speed_data = join_horiz(speed_data, join_vert(k_speed, d_speed));
    teta = join_horiz(teta, controller.parameters.dynamic_parameters);
    loop.wait();                
  } while (error(pose, trajectory_pose) > 250.0);

  // to load on Matlab: load('file.mat','-ASCII');
  telemetry_data.save(results_folder + "telemetry.mat", raw_ascii);  
  trajectory_data.save(results_folder + "trajectory.mat", raw_ascii);
  speed_data.save(results_folder + "speed.mat", raw_ascii);  
  teta.save(results_folder + "teta.mat", raw_ascii);  

  robot.disconnect();
  return 0;
}