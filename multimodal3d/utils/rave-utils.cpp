#include "rave-utils.h"

std::vector<rave::GraphHandlePtr> handles;

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

void clear_plots(int num) {
	if ((num < 0) || (num > handles.size())) {
		handles.clear();
	} else {
		handles = std::vector<rave::GraphHandlePtr>(handles.begin(), handles.end()-num);
	}
}

void plot_point(rave::EnvironmentBasePtr env, rave::Vector pos, rave::Vector color, float size) {
	float *ppoints = new float[3];
	ppoints[0] = pos.x;
	ppoints[1] = pos.y;
	ppoints[2] = pos.z;

	color.w = 1;

	handles.push_back(env->plot3(ppoints, 1, sizeof(float), size, color, 1));
}


void plot_point(rave::EnvironmentBasePtr env, Vector3d pos, Vector3d color, float size) {
	rave::Vector pos_vec = eigen_to_rave(pos);
	rave::Vector color_vec = eigen_to_rave(color);
	plot_point(env, pos_vec, color_vec, size);
}

void plot_segment(rave::EnvironmentBasePtr env, const Vector3d& p0, const Vector3d& p1, Vector3d& color) {
	float points[6];
	points[0] = p0(0); points[1] = p0(1); points[2] = p0(2);
	points[3] = p1(0); points[4] = p1(1); points[5] = p1(2);
	float c[6];
	for(int i=0; i < 3; ++i) { c[i] = color(i); }
	for(int i=0; i < 3; ++i) { c[i+3] = color(i); }
	handles.push_back(env->drawlinestrip(points, 2, sizeof(float)*3, 1.0f, c));
}

//void plot_transform(rave::EnvironmentBasePtr env, rave::Transform T, std::vector<rave::GraphHandlePtr> &handles) {
//	fmat M = conv_to<fmat>::from(rave_transform_to_mat(T));
//
//	fmat o = M.submat(span(0, 2), span(3, 3));
//	fmat I = join_vert(eye<fmat>(3,3), eye<fmat>(3,3));
//	fmat ppoints(6,1);
//	ppoints.submat(span(0,2), span(0,0)) = o;
//	for(int i=0; i < 3; ++i) {
//		ppoints.submat(span(3,5), span(0, 0)) = o + 0.1 * M.submat(span(0, 2), span(i, i));
//		fmat color = I.submat(span(0, 5), span(i, i));
//		handles.push_back(env->drawlinestrip(ppoints.colptr(0), 2, sizeof(float)*3, 1.0f, color.colptr(0)));
//	}
//}

void save_view(rave::ViewerBasePtr viewer, std::string file_name) {
	std::string input = "SetFiguresInCamera 1", output = "";
	viewer->SendCommand(input, output);

	int width = 640, height = 480;
	std::vector<uint8_t> image(3*width*height);
	rave::SensorBase::CameraIntrinsics intrinsics(640, 640, 320, 240);
	viewer->GetCameraImage(image, width, height, viewer->GetCameraTransform(), intrinsics);

//	FILE *f = fopen("out.ppm", "wb");
//	fprintf(f, "P6\n%i %i 255\n", width, height);
//	int index = 0;
	for (int y=0; y < height; ++y) {
		for (int x=0; x < width; ++x) {
			for(int i=0; i < 3; ++i) {
//				fputc(image[index++], f);   // 0 .. 255
			}
		}
	}
//	fclose(f);
}


}
