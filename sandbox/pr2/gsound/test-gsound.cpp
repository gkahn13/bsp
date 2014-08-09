#include <iostream>
#include <map>
#include "../../util/Timer.h"
namespace my_util = util;

#include <openrave-core.h>
namespace rave = OpenRAVE;

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include <Eigen/Eigen>
using namespace Eigen;

#include "gsound/GSound.h"
using namespace gsound;

#include "../utils/rave-utils.h"
#include "../utils/utils.h"

inline Vector3 eigen_to_gsound(const Vector3d& v) {
	return Vector3(v(0), v(1), v(2));
}

inline Matrix3 eigen_to_gsound(const Matrix3d& m) {
	Matrix3 mg;
	for(int i=0; i < 3; ++i) {
		for(int j=0; j < 3; ++j) {
			mg(i,j) = m(i,j);
		}
	}
	return mg;
}

inline Vector3d gsound_to_eigen(const Vector3& v) {
	return Vector3d(v(0), v(1), v(2));
}

inline math::Vector3d rave_to_gsound(const rave::Vector& v) {
	return math::Vector3d(v.x, v.y, v.z);
}

void SetViewer(rave::EnvironmentBasePtr penv, const std::string& viewername) {
    rave::ViewerBasePtr viewer = RaveCreateViewer(penv,viewername);
    BOOST_ASSERT(!!viewer);

    // attach it to the environment:
    penv->Add(viewer);

    // finally call the viewer's infinite loop (this is why a separate thread is needed)
    bool showgui = true;
    viewer->main(showgui);
}

void init_rave(rave::EnvironmentBasePtr& env, boost::shared_ptr<boost::thread>& viewer_thread) {
	rave::RaveInitialize(true, rave::Level_Info);
	env = rave::RaveCreateEnvironment();

	viewer_thread = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(SetViewer, env, "qtcoin")));
	sleep(1);
}


SoundMesh* loadBox(AABB3 box, SoundMaterial& material) {
	ArrayList<SoundVertex> vertices;
	ArrayList<SoundTriangle> triangles;
	ArrayList<SoundMaterial> materials;

	double x_min = box.min(0), y_min = box.min(1), z_min = box.min(2);
	double x_max = box.max(0), y_max = box.max(1), z_max = box.max(2);

	vertices.add(math::Vector3d(x_min, y_min, z_min)); // 0
	vertices.add(math::Vector3d(x_min, y_min, z_max)); // 1
	vertices.add(math::Vector3d(x_min, y_max, z_min)); // 2
	vertices.add(math::Vector3d(x_min, y_max, z_max)); // 3
	vertices.add(math::Vector3d(x_max, y_min, z_min)); // 4
	vertices.add(math::Vector3d(x_max, y_min, z_max)); // 5
	vertices.add(math::Vector3d(x_max, y_max, z_min)); // 6
	vertices.add(math::Vector3d(x_max, y_max, z_max)); // 7

	// x end 0
	triangles.add(SoundTriangle(0, 1, 2, 0));
	triangles.add(SoundTriangle(3, 1, 2, 0));
	// x end 1
	triangles.add(SoundTriangle(4, 5, 6, 0));
	triangles.add(SoundTriangle(7, 5, 6, 0));
	// y end 0
	triangles.add(SoundTriangle(0, 1, 4, 0));
	triangles.add(SoundTriangle(5, 1, 4, 0));
	// y end 1
	triangles.add(SoundTriangle(2, 3, 6, 0));
	triangles.add(SoundTriangle(7, 3, 6, 0));
	// z end 0
	triangles.add(SoundTriangle(0, 2, 4, 0));
	triangles.add(SoundTriangle(6, 2, 4, 0));
	// z end 1
	triangles.add(SoundTriangle(1, 3, 5, 0));
	triangles.add(SoundTriangle(7, 3, 4, 0));

	materials.add(material);

	return new SoundMesh(vertices, triangles, materials);
}

void addRaveToGSound(const rave::EnvironmentBasePtr env, SoundScene& scene) {
	// Create a material object which will be the material for the box's surface.
	// The material is specified in terms of FrequencyResponse objects that dictate
	// how sound is affected when it undergoes the associated interaction.
	SoundMaterial material(
			// The reflection attenuation for the material.
			FrequencyResponse::getLinearHighRolloff(1000)*
			FrequencyResponse::getLinearLowRolloff(200)*0.9,
			// The transmission attenuation per world unit for the material.
			FrequencyResponse::getQuadraticHighRolloff(800)*0.9,
			// The absorption attenuation for the material.
			FrequencyResponse::getLinearHighRolloff()*0.5 );


	std::vector<rave::KinBodyPtr> bodies;
	env->GetBodies(bodies);

	for(int i=0; i < bodies.size(); ++i) {
		ArrayList<SoundVertex> vertices;
		ArrayList<SoundTriangle> triangles;
		ArrayList<SoundMaterial> materials;

		Vector3d body_pos(0, 1.5, 0); // TODO: tmp, hardcoded
		std::cout << "body_pos: " << body_pos.transpose() << "\n";
		std::vector<rave::KinBody::LinkPtr> links = bodies[i]->GetLinks();
		for(int j=0; j < links.size(); ++j) {
			std::vector<rave::KinBody::Link::GeometryPtr> geometries = links[j]->GetGeometries();
			for(int k=0; k < geometries.size(); ++k) {
				rave::TriMesh mesh = geometries[k]->GetCollisionMesh();
				std::cout << "num indices: " << mesh.indices.size() << "\n";
				std::cout << "num vertices: " << mesh.vertices.size() << "\n";
				for(int l = 0; l < mesh.vertices.size(); ++l) {
					vertices.add(rave_to_gsound(mesh.vertices[l]));
					rave_utils::plot_point(env, body_pos+rave_utils::rave_to_eigen(mesh.vertices[l]), Vector3d(1,0,0), .5);
				}
				for(int m=0; m < mesh.indices.size(); m+=3) {
					triangles.add(SoundTriangle(mesh.indices[m], mesh.indices[m+1], mesh.indices[m+2], 0));
				}
			}
		}

		materials.add(material);

		SoundObject* sound_object = new SoundObject(new SoundMesh(vertices, triangles, materials));
		sound_object->setPosition(Vector3(body_pos(0), body_pos(1), body_pos(2)));
		scene.addObject(sound_object);
	}

}

int main ( int argc, char** argv ) {
	// initialize Openrave for viewing
	rave::EnvironmentBasePtr env;
	boost::shared_ptr<boost::thread> viewer_thread;
	init_rave(env, viewer_thread);
	env->Load("/home/gkahn/Research/bsp/pr2/envs/table-env.env.xml");

	// The object which performs sound propagation
	SoundPropagator propagator;

	// Enable all types of propagation paths.
	propagator.setDirectSoundIsEnabled(true);
	propagator.setTransmissionIsEnabled(true);
	propagator.setReflectionIsEnabled(false);
	propagator.setDiffractionIsEnabled(true);
	propagator.setReverbIsEnabled(false);

	//***********************************************************************

	// The object which contains a collection of all sound objects,
	// sound sources, and listeners in the scene.
	SoundScene scene;

	// The object which specifies the location and orientation of
	// the sound receiver in the scene.
	SoundListener listener;

	// Set the listener's starting position, at head height and centered at the XZ origin.
	Vector3d listener_pos(0, 1.5, 0);
	Matrix3d listener_ori = Matrix3d::Identity();

	listener.setPosition(eigen_to_gsound(listener_pos));
	listener.setOrientation(eigen_to_gsound(listener_ori));

	Matrix4d listener_pose = Matrix4d::Identity();
	listener_pose.block<3,1>(0,3) = listener_pos;
	listener_pose.block<3,3>(0,0) = listener_ori;

	//***********************************************************************

//	// Create a material object which will be the material for the box's surface.
//	// The material is specified in terms of FrequencyResponse objects that dictate
//	// how sound is affected when it undergoes the associated interaction.
//	SoundMaterial defaultMaterial(
//			// The reflection attenuation for the material.
//			FrequencyResponse::getLinearHighRolloff(1000)*
//			FrequencyResponse::getLinearLowRolloff(200)*0.9,
//			// The transmission attenuation per world unit for the material.
//			FrequencyResponse::getQuadraticHighRolloff(800)*0.9,
//			// The absorption attenuation for the material.
//			FrequencyResponse::getLinearHighRolloff()*0.5 );
//
//
//	// Create an axis-aligned box that is 4x8x3 meters, centered
//	// at the origin that uses the default material.
//	AABB3 a(-2, 2, -1.5, 1.5, -4, 4);
//	SoundMesh* box = loadBox(AABB3(-2, 2, -1.5, 1.5, -4, 4), defaultMaterial);
//
//	SoundObject* boxObject = new SoundObject(box);
//
//	// Set the position of the box so that it is now 1.5 units higher along the Y-axis.
//	boxObject->setPosition(Vector3(0, 1.5, 0));
//
//	// Add the box to the scene
//	scene.addObject(boxObject);

	addRaveToGSound(env, scene);

	//***********************************************************************

	SoundSource* source = new SoundSource();

	// Set the position of the source so that it is on one side of the box.
	Vector3d source_pos(0, 1, -3);
	source->setPosition(eigen_to_gsound(source_pos));
	rave_utils::plot_point(env, source_pos, Vector3d(0, 0, 1), .5);

	// Set the source to have an intensity of 1.
	// This is the gain applied to the source's audio when there is no
	// distance attenuation.
	source->setIntensity(1);

	// Create a distance attenuation object which specifies how the source's
	// audio decreases in intensity with distance.
	SoundDistanceAttenuation attenuation(1, // Constant attenuation of 1.
			1, // Linear attenuation of 1.
			0); // Quadratic attenuation of 0.

	// Set the distance attenuation for the source.
	source->setDistanceAttenuation(attenuation);

	// Set the reverb distance attenuation for the source to slightly less
	// than the normal distance attenuation.
	source->setReverbDistanceAttenuation(SoundDistanceAttenuation(1, 0.5, 0));

	// Add the sound source to the scene.
	scene.addSource( source );

	// Create a buffer to hold the output of the propagation system.
	SoundPropagationPathBuffer pathBuffer;

	double pos_step = .05;
	std::map<int,rave::Vector> delta_position =
	{
			{'a' , rave::Vector(0, pos_step, 0)},
			{'d' , rave::Vector(0, -pos_step, 0)},
			{'w' , rave::Vector(pos_step, 0, 0)},
			{'x' , rave::Vector(-pos_step, 0, 0)},
			{'+' , rave::Vector(0, 0, pos_step)},
			{'-' , rave::Vector(0, 0, -pos_step)},
	};

	double angle_step = 2.0*(M_PI/180);
	std::map<int,rave::Vector> delta_angle =
	{
			{'p', rave::Vector(angle_step, 0, 0)},
			{'o', rave::Vector(-angle_step, 0, 0)},
			{'k', rave::Vector(0, angle_step, 0)},
			{'l', rave::Vector(0, -angle_step, 0)},
			{'n', rave::Vector(0, 0, angle_step)},
			{'m', rave::Vector(0, 0, -angle_step)},
	};


	my_util::Timer timer;
	double elapsed = 0;
	my_util::Timer_tic(&timer);
	int iters = 0;

	rave_utils::plot_transform(env, rave_utils::eigen_to_rave(listener_pose), 1.5);
	int c;
	while ((c = utils::getch()) != 'q') {
		// update based on keyboard
		rave::Transform pose = rave_utils::eigen_to_rave(listener_pose);
		if (delta_position.count(c) > 0) {
			pose.trans += delta_position[c];
		} else if (delta_angle.count(c) > 0) {
			pose.rot = rave::geometry::quatFromAxisAngle(rave::geometry::axisAngleFromQuat(pose.rot) + delta_angle[c]);
		}
		listener_pose = rave_utils::rave_to_eigen(pose);

		Vector3d pos = listener_pose.block<3,1>(0,3);
		Matrix3d ori = listener_pose.block<3,3>(0,0);

		listener.setPosition(eigen_to_gsound(pos));
		listener.setOrientation(eigen_to_gsound(ori));

		rave_utils::clear_plots(3);
		rave_utils::plot_transform(env, rave_utils::eigen_to_rave(listener_pose), 1.5);

		// Perform sound propagation in the scene.
		propagator.propagateSound( scene, // The scene in which to perform propagation.
				listener, // The listener to use as the sound receiver.
				4, // The maximum depth of the rays shot from the listener.
				1000, // The number of rays to shoot from the listener,
				// influences the quality of the early reflection paths.
				4, // The maximum depth of the rays shot from each sound source.
				100, // The number of rays to shoot from each sound source,
				// influences reverb estimation quality.
				pathBuffer ); // The buffer in which to put the propagation paths.

//		std::cout << "Number of sources: " << pathBuffer.getNumberOfSources() << "\n";
//		std::cout << "Number of propagation paths: " << pathBuffer.getTotalNumberOfPropagationPaths() << "\n";

		SoundSourcePropagationPathBuffer sourcePathBuffer = pathBuffer.getSourceBuffer(0);
		int num_prop_paths = sourcePathBuffer.getNumberOfPropagationPaths();
		std::cout << "Number of propagation paths: " << num_prop_paths << "\n";
		double total_avg_gain = 0;
		double max_gain = 0, min_gain = INFINITY;
		for(int i=0; i < num_prop_paths; ++i) {
			PropagationPath propPath = sourcePathBuffer.getPropagationPath(i);
			FrequencyResponse fr = propPath.getFrequencyAttenuation();
			double avg_gain = fr.getAverageGain();
//			std::cout << "gain " << i << ": " << avg_gain << "\n";
			total_avg_gain += avg_gain;
			max_gain = (max_gain > avg_gain) ? max_gain : avg_gain;
			min_gain = (min_gain < avg_gain) ? min_gain : avg_gain;
		}
		double all_avg_gain = total_avg_gain / double(num_prop_paths);
		std::cout << "all_avg_gain: " << all_avg_gain << "\n";
		std::cout << "max_gain: " << max_gain << "\n";
		std::cout << "min_gain: " << min_gain << "\n";

		elapsed += my_util::Timer_toc(&timer);
		my_util::Timer_tic(&timer);
		iters++;

		sleep(.05);

//		std::cout << "Press enter\n";
//		std::cin.ignore();
	}


//	// Perform sound propagation for 10 seconds.
//	while (elapsed < 10.0) {
//		// Perform sound propagation in the scene.
//		propagator.propagateSound( scene, // The scene in which to perform propagation.
//				listener, // The listener to use as the sound receiver.
//				4, // The maximum depth of the rays shot from the listener.
//				1000, // The number of rays to shoot from the listener,
//				// influences the quality of the early reflection paths.
//				4, // The maximum depth of the rays shot from each sound source.
//				100, // The number of rays to shoot from each sound source,
//				// influences reverb estimation quality.
//				pathBuffer ); // The buffer in which to put the propagation paths.
//
////		std::cout << "Number of sources: " << pathBuffer.getNumberOfSources() << "\n";
////		std::cout << "Number of propagation paths: " << pathBuffer.getTotalNumberOfPropagationPaths() << "\n";
//
//		SoundSourcePropagationPathBuffer sourcePathBuffer = pathBuffer.getSourceBuffer(0);
//		int num_prop_paths = sourcePathBuffer.getNumberOfPropagationPaths();
//		std::cout << "Number of propagation paths: " << num_prop_paths << "\n";
//		double total_avg_gain = 0;
//		for(int i=0; i < num_prop_paths; ++i) {
//			PropagationPath propPath = sourcePathBuffer.getPropagationPath(i);
//			FrequencyResponse fr = propPath.getFrequencyAttenuation();
//			double avg_gain = fr.getAverageGain();
////			std::cout << "gain " << i << ": " << avg_gain << "\n";
//			total_avg_gain += avg_gain;
//		}
//		double all_avg_gain = total_avg_gain / double(num_prop_paths);
//		std::cout << "all_avg_gain: " << all_avg_gain << "\n";
//
//		elapsed += my_util::Timer_toc(&timer);
//		my_util::Timer_tic(&timer);
//		iters++;
//
////		std::cout << "Press enter\n";
////		std::cin.ignore();
//	}
//
//	double avg_time = elapsed / double(iters);
//	std::cout << "avg time: " << 1e3*avg_time << " ms\n";


	std::cout << "Exiting\n";
}
