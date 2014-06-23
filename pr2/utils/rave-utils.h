#ifndef _RAVE_UTILS_H__
#define _RAVE_UTILS_H__

#include <iostream>

#include <Eigen/Eigen>
using namespace Eigen;

#include <openrave-core.h>
namespace rave = OpenRAVE;

namespace rave_utils {

void cart_to_joint(rave::RobotBase::ManipulatorPtr manip, const rave::Transform &matrix4,
		const std::string &ref_frame, const std::string &targ_frame, std::vector<double> &joint_values,
		const int filter_options=0);

rave::Transform transform_relative_pose_for_ik(rave::RobotBase::ManipulatorPtr manip,
		const rave::Transform &matrix4, const std::string &ref_frame, const std::string &targ_frame);

void clear_plots(int num=-1);

void plot_point(rave::EnvironmentBasePtr env, rave::Vector pos, rave::Vector color, float size=.01);

void plot_point(rave::EnvironmentBasePtr env, Vector3d pos, Vector3d color, float size=.01);

void plot_segment(rave::EnvironmentBasePtr env, const Vector3d& p0, const Vector3d& p1, Vector3d& color);

void plot_transform(rave::EnvironmentBasePtr env, rave::Transform T);

void save_view(rave::ViewerBasePtr viewer, std::string file_name);

inline Matrix4d rave_to_eigen(const rave::Transform& rt) {
	rave::TransformMatrix rtm(rt);

	Matrix4d m = Matrix4d::Identity();
	int index = 0;
	for(int i=0; i < 3; ++i) {
		for(int j=0; j < 4; ++j) {
			m(i,j) = rtm.m[index++];
		}
	}

	for(int i=0; i < 3; ++i) {
		m(i, 3) = rtm.trans[i];
	}

	return m;
}

inline rave::Transform eigen_to_rave(const Matrix4d& m) {
	rave::TransformMatrix rtm;

	int index = 0;
	for(int i=0; i < 3; ++i) {
		for(int j=0; j < 4; ++j) {
			rtm.m[index++] = m(i,j);
		}
	}

	for(int i=0; i < 3; ++i) {
		rtm.trans[i] = m(i,3);
	}

	rave::Transform rt(rtm);
	return rt;
}

inline Vector3d rave_to_eigen(const rave::Vector& v) {
	return Vector3d(v.x, v.y, v.z);
}

inline rave::Vector eigen_to_rave(const Vector3d& m) {
	return rave::Vector(m(0), m(1), m(2));
}

}

#endif
