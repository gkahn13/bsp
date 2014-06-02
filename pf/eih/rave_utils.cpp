#include "rave_utils.h"

inline double smaller_ang(double x) {
	return fmod(x + M_PI, 2*M_PI) - M_PI;
}

inline double closer_ang(double x, double a) {
	return a + smaller_ang(x-a);
}

namespace rave_utils {

void cart_to_joint(rave::RobotBase::ManipulatorPtr manip, const mat &matrix4,
		const std::string &ref_frame, const std::string &targ_frame, std::vector<double> &joint_values,
		const int filter_options) {
	std::vector<double> initial_joint_values;
	manip->GetArmDOFValues(initial_joint_values);

	std::cout << "ref_frame: " << ref_frame << "\n";
	std::cout << "targ_frame: " << targ_frame << "\n";
	rave::RobotBasePtr robot = manip->GetRobot();
	std::cout << "matrix4\n" << matrix4 << "\n";
	mat world_from_EE = transform_relative_pose_for_ik(manip, matrix4, ref_frame, targ_frame);
	std::cout << "world_from_EE\n" << world_from_EE << "\n";
	bool success = manip->FindIKSolution(rave::IkParameterization(mat_to_rave_transform(world_from_EE)),
			joint_values, filter_options);

	if (!success) {
		std::cout << "IK Failed. Returning original joint values\n";
		joint_values = initial_joint_values;
	} else {
		std::cout << "IK Success\n";
		for(int i=0; i < joint_values.size(); ++i) { std::cout << joint_values[i] << " "; }
		std::cout << "\n";
		for(int i=2; i <= 6; i += 2) {
			joint_values[i] = closer_ang(joint_values[i], initial_joint_values[i]);
		}
	}

}

mat transform_relative_pose_for_ik(rave::RobotBase::ManipulatorPtr manip,
		const mat &matrix4, const std::string &ref_frame, const std::string &targ_frame) {
	mat world_from_EE_new;
	rave::RobotBasePtr robot = manip->GetRobot();

	mat world_from_ref, targ_from_EE;
	if (ref_frame == "world") {
		world_from_ref.eye(4,4);
		std::cout << "ref_frame == world\n";
		std::cout << world_from_ref << "\n";
	} else {
		rave::RobotBase::LinkPtr ref_link = robot->GetLink(ref_frame);
		world_from_ref = rave_transform_to_mat(ref_link->GetTransform());
	}

	if (targ_frame == "end_effector") {
		targ_from_EE.eye(4,4);
		std::cout << "targ_frame == end_effector\n";
		std::cout << targ_from_EE << "\n";
	} else {
		rave::RobotBase::LinkPtr targ_link = robot->GetLink(targ_frame);
		mat world_from_targ = rave_transform_to_mat(targ_link->GetTransform());
		mat world_from_EE = rave_transform_to_mat(manip->GetEndEffectorTransform());
		targ_from_EE = inv(world_from_targ)*world_from_EE;
	}

	mat ref_from_targ_new = matrix4;
	std::cout << "world_from_ref*ref_from_targ_new\n" << world_from_ref*ref_from_targ_new << "\n";
	world_from_EE_new = (world_from_ref*ref_from_targ_new)*targ_from_EE;

	return world_from_EE_new;
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

}
