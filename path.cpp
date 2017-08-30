#include <iostream>
#include "controller.hpp"
#include "pioneer.hpp"

using namespace std;
using namespace arma;
using namespace robot::local;
// using namespace robot::remote;

int main(int argc, char** argv) {
  std::string results_folder = argc == 2 ? argv[1] : "./";
  Pioneer robot;
  // Pioneer robot("192.168.1.143", 1110);

  mat initial_pose;
  initial_pose << 0.0 << endr << 0.0 << endr << deg2rad(0.0);
  robot.set_pose(initial_pose);

  controller::Path::Parameters path_parameters("path.yaml");
  controller::Trajectory::Parameters controller_parameters("path_controller.yaml");
  controller::Path path_controller(path_parameters, controller_parameters);

  mat telemetry_data;
  mat speed_data;
  mat trajectory_data;
  Loop loop;
  do {
    mat telemetry = robot.get_telemetry();  // returns a 5x1 matrix {x,y,psi,v,w} [mm,mm,rad,mm/s,rad/s]
    mat pose = telemetry.rows(0, 2);        // {x,y,psi}
    mat trajectory_pose;
    auto k_speed = path_controller.kinematics(pose, trajectory_pose);
    //auto d_speed = controller.dynamics_2(telemetry, k_speed);
    robot.set_speed(k_speed);

    pose.print("Pose: ");

    telemetry_data = join_horiz(telemetry_data, telemetry);
    trajectory_data = join_horiz(trajectory_data, trajectory_pose);
    //speed_data = join_horiz(speed_data, join_vert(k_speed, d_speed));
    speed_data = join_horiz(speed_data, k_speed);
    loop.wait();
  } while (path_controller.on_path());

  // to load on Matlab: load('file.mat','-ASCII');
  telemetry_data.save(results_folder + "telemetry.mat", raw_ascii);
  trajectory_data.save(results_folder + "trajectory.mat", raw_ascii);
  speed_data.save(results_folder + "speed.mat", raw_ascii);

  robot.disconnect();
  return 0;
}