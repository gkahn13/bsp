#pragma once

#include <iostream>

#include <armadillo>
using namespace arma;

#include <openrave-core.h>
namespace rave = OpenRAVE;

namespace rave_utils {

void cart_to_joint(rave::RobotBase::ManipulatorPtr manip, const mat &matrix4,
		const std::string &ref_frame, const std::string &targ_frame, std::vector<double> &joint_values,
		const int filter_options=0);

mat transform_relative_pose_for_ik(rave::RobotBase::ManipulatorPtr manip,
		const mat &matrix4, const std::string &ref_frame, const std::string &targ_frame);

mat rave_transform_to_mat(rave::Transform rt);

rave::Transform mat_to_rave_transform(mat m);

}
