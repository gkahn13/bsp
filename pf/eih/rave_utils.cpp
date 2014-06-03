#include "rave_utils.h"

inline double smaller_ang(double x) {
	return fmod(x + M_PI, 2*M_PI) - M_PI;
}

inline double closer_ang(double x, double a) {
	return a + smaller_ang(x-a);
}

namespace rave_utils {

void cart_to_joint(rave::RobotBase::ManipulatorPtr manip, const rave::Transform &matrix4,
		const std::string &ref_frame, const std::string &targ_frame, std::vector<double> &joint_values,
		const int filter_options) {
	std::vector<double> initial_joint_values;
	manip->GetArmDOFValues(initial_joint_values);

	rave::RobotBasePtr robot = manip->GetRobot();
	rave::Transform world_from_EE = transform_relative_pose_for_ik(manip, matrix4, ref_frame, targ_frame);
	bool success = manip->FindIKSolution(rave::IkParameterization(world_from_EE), joint_values, filter_options);

	if (!success) {
		RAVELOG_ERROR("IK Failed. Returning current joint values\n");
		joint_values = initial_joint_values;
	} else {
		for(int i=2; i <= 6; i += 2) {
			joint_values[i] = closer_ang(joint_values[i], initial_joint_values[i]);
		}
	}

}

rave::Transform transform_relative_pose_for_ik(rave::RobotBase::ManipulatorPtr manip,
		const rave::Transform &matrix4, const std::string &ref_frame, const std::string &targ_frame) {
	rave::RobotBasePtr robot = manip->GetRobot();

	rave::Transform world_from_ref, targ_from_EE;
	if (ref_frame == "world") {
		world_from_ref.identity();
	} else {
		rave::RobotBase::LinkPtr ref_link = robot->GetLink(ref_frame);
		world_from_ref = ref_link->GetTransform();
	}

	if (targ_frame == "end_effector") {
		targ_from_EE.identity();
	} else {
		rave::RobotBase::LinkPtr targ_link = robot->GetLink(targ_frame);
		rave::Transform world_from_targ = targ_link->GetTransform();
		rave::Transform  world_from_EE = manip->GetEndEffectorTransform();
		targ_from_EE = world_from_targ.inverse()*world_from_EE;
	}

	rave::Transform ref_from_targ_new = matrix4;
	rave::Transform world_from_EE_new = (world_from_ref*ref_from_targ_new)*targ_from_EE;

	return world_from_EE_new;
}

rave::GraphHandlePtr plot_point(rave::EnvironmentBasePtr env, const rave::Vector &pos, rave::Vector &color, float size) {
	float *ppoints = new float[3];
	ppoints[0] = pos.x;
	ppoints[1] = pos.y;
	ppoints[2] = pos.z;

	color.w = 1;

	rave::GraphHandlePtr h = env->plot3(ppoints, 1, sizeof(float), size, color, 1);
	return h;
}


rave::GraphHandlePtr plot_point(rave::EnvironmentBasePtr env, const mat &pos, mat &color, float size) {
	rave::Vector pos_vec = mat_to_rave_vec(pos);
	rave::Vector color_vec = mat_to_rave_vec(color);
	return plot_point(env, pos_vec, color_vec, size);
}

mat rave_transform_to_mat(rave::Transform rt) {
	rave::TransformMatrix rtm(rt);

	mat m;
	m.eye(4,4);
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

rave::Transform mat_to_rave_transform(mat m) {
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

mat rave_vec_to_mat(rave::Vector v) {
	mat m;
	m << v.x << v.y << v.z << v.w;
	return m;
}

rave::Vector mat_to_rave_vec(mat m) {
	rave::Vector v(m(0), m(1), m(2));
	if (m.n_rows*m.n_cols == 4) {
		v.w = m(3);
	}
	return v;
}


}
