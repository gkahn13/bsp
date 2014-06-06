#ifndef _RAVE_UTILS_H__
#define _RAVE_UTILS_H__

#include <iostream>

#include <armadillo>
using namespace arma;

#include <openrave-core.h>
namespace rave = OpenRAVE;

namespace rave_utils {

void cart_to_joint(rave::RobotBase::ManipulatorPtr manip, const rave::Transform &matrix4,
		const std::string &ref_frame, const std::string &targ_frame, std::vector<double> &joint_values,
		const int filter_options=0);

rave::Transform transform_relative_pose_for_ik(rave::RobotBase::ManipulatorPtr manip,
		const rave::Transform &matrix4, const std::string &ref_frame, const std::string &targ_frame);

rave::GraphHandlePtr plot_point(rave::EnvironmentBasePtr env, const rave::Vector &pos, rave::Vector &color, float size=.01);

rave::GraphHandlePtr plot_point(rave::EnvironmentBasePtr env, const mat &pos, mat &color, float size=.01);

void plot_transform(rave::EnvironmentBasePtr env, rave::Transform T, std::vector<rave::GraphHandlePtr> &handles);

void save_view(rave::ViewerBasePtr viewer, std::string file_name);

mat rave_transform_to_mat(rave::Transform rt);

rave::Transform mat_to_rave_transform(mat m);

mat rave_vec_to_mat(rave::Vector v, bool is_4d=false);

rave::Vector mat_to_rave_vec(mat m);

}

#endif
