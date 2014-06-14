#include "fadbad-planar-system.h"

/**
 * Constructors and initializers
 */

FadbadPlanarSystem::FadbadPlanarSystem(const vecb<2>& camera_origin, const vecb<2>& object, bool is_static) {
	init(camera_origin, object, is_static);
}

void FadbadPlanarSystem::init(const vecb<2>& camera_origin, const vecb<2>& object, bool is_static) {
	this->camera_origin = camera_origin;
	camera_fov = M_PI/4;
	camera_max_dist = 15;
	this->object = object;
	this->is_static = is_static;

	robot_origin = vecb<2>::Zero();
	link_lengths << 5, 5, 4;

	Q = (M_PI/4)*matb<U_DIM,U_DIM>::Identity();
//	R = 10*eye<matb>(Z_DIM, Z_DIM);
	vecb<R_DIM> R_diag;
	R_diag << (M_PI/4)*vecb<J_DIM>::Ones(),
			5*vecb<C_DIM>::Ones();
	R = R_diag.asDiagonal();

	// x bound for object is sum of link_lengths
	bdouble max_link_length = link_lengths.sum();
	x_min << -M_PI/2, -M_PI/2, -M_PI/2, -M_PI/2, -max_link_length, -max_link_length;
	x_max << M_PI/2, M_PI/2, M_PI/2, M_PI/2, max_link_length, max_link_length;

	bdouble max_input = M_PI/12;
	u_min << -max_input, -max_input, -max_input, (is_static) ? 0 : -max_input;
	u_max << max_input, max_input, max_input, (is_static) ? 0 : max_input;
}


/**
 * Public methods
 */

vecb<X_DIM> FadbadPlanarSystem::dynfunc(const vecb<X_DIM>& x, const vecb<U_DIM>& u, const vecb<Q_DIM>& q, bool enforce_limits) {
	vecb<X_DIM> x_new = x;
	x_new.segment<J_DIM>(0) += DT*(u + q);
	return x_new;
}

vecb<Z_DIM> FadbadPlanarSystem::obsfunc(const vecb<X_DIM>& x, const vecb<2>& object, const vecb<R_DIM>& r) {
	vecb<Z_DIM> z;
	z(0) = x(0); // joint 0
	z(1) = x(1); // joint 1
	z(2) = x(2); // joint 2
	z(3) = x(3); // camera angle
	z(4) = object(0) - camera_origin(0); // delta x to object
	z(5) = object(1) - camera_origin(1); // delta y to object
	return z + r;
}


matb<Z_DIM,Z_DIM> FadbadPlanarSystem::delta_matrix(const vecb<X_DIM>& x, const vecb<2>& object, const bdouble alpha) {
	matb<Z_DIM,Z_DIM> delta = matb<Z_DIM,Z_DIM>::Zero();

	for(int i=0; i < J_DIM; ++i) {
		delta(i, i) = 1; // TODO: should this depend on SD of link segments?
	}

	std::vector<BeamB> fov = get_fov(x);
	bdouble sd = fadbad_geometry2d::signed_distance(object, fov);
	bdouble sd_sigmoid = 1.0 - 1.0/(1.0 + exp(-alpha*sd));
	for(int i=J_DIM; i < X_DIM; ++i) {
		delta(i, i) = sd_sigmoid;
	}

	return delta;
}

/**
 * \brief Propagates belief through EKF with max-likelihood noise assumption
 */
void FadbadPlanarSystem::belief_dynamics(const vecb<X_DIM>& x_t, const matb<X_DIM,X_DIM>& sigma_t, const vecb<U_DIM>& u_t, const bdouble alpha,
		vecb<X_DIM>& x_tp1, matb<X_DIM,X_DIM>& sigma_tp1) {
	x_tp1.Zero();
	sigma_tp1.Zero();

	// propagate dynamics
	x_tp1 = dynfunc(x_t, u_t, vecb<Q_DIM>::Zero());

	// propagate belief through dynamics
	matb<X_DIM,X_DIM> A = matb<X_DIM,X_DIM>::Zero();
	matb<X_DIM,Q_DIM> M = matb<X_DIM,Q_DIM>::Zero();
	linearize_dynfunc(x_t, u_t, vecb<Q_DIM>::Zero(), A, M);

	matb<X_DIM,X_DIM> sigma_tp1_bar = A*sigma_t*A.transpose() + M*Q*M.transpose();

	// propagate belief through observation
	matb<Z_DIM,X_DIM> H = matb<Z_DIM,X_DIM>::Zero();
	matb<Z_DIM,R_DIM> N = matb<Z_DIM,R_DIM>::Zero();
	linearize_obsfunc(x_tp1, vecb<R_DIM>::Zero(), H, N);

	matb<Z_DIM,Z_DIM> delta = delta_matrix(x_tp1, x_tp1.segment<C_DIM>(J_DIM), alpha);
	matb<X_DIM,Z_DIM> K = sigma_tp1_bar*H.transpose()*delta*(delta*H*sigma_tp1_bar*H.transpose()*delta + R).fullPivLu().solve(matb<Z_DIM,Z_DIM>::Identity())*delta;
	sigma_tp1 = (matb<X_DIM,X_DIM>::Identity() - K*H)*sigma_tp1_bar;
}

std::vector<BeamB> FadbadPlanarSystem::get_fov(const vecb<X_DIM>& x) {
	std::vector<BeamB> beams, new_beams;

	// start out with full field of view
	bdouble angle = x(J_DIM-1);
	vecb<2> right, left;
	right << camera_max_dist*sin(angle+camera_fov/2.0),
			camera_max_dist*cos(angle+camera_fov/2.0);
	left << camera_max_dist*sin(angle-camera_fov/2.0),
			camera_max_dist*cos(angle-camera_fov/2.0);
	beams.push_back(BeamB(camera_origin,
			camera_origin + right,
			camera_origin + left));

	std::vector<SegmentB> link_segments = get_link_segments(x);

	for(int l=0; l < link_segments.size(); ++l) {
		new_beams.clear();
		for(int i=0; i < beams.size(); ++i) {
			std::vector<BeamB> new_beams_i = beams[i].truncate(link_segments[l]);
			for(int j=0; j < new_beams_i.size(); ++j) {
				// prevent small beams from being added
				if (new_beams_i[j].top_length() > bepsilon) {
					new_beams.push_back(new_beams_i[j]);
				}
			}
		}
		beams = new_beams;
	}

	return beams;
}

std::vector<SegmentB> FadbadPlanarSystem::get_link_segments(const vecb<X_DIM>& x) {
	std::vector<SegmentB> link_segments;

	vecb<2> p0 = robot_origin, p1;
	bdouble joint0 = 0, joint1;

	for(int i=0; i < link_lengths.rows(); ++i) {
		joint1 = joint0 + x[i];
		p1 << p0(0) + sin(joint1)*link_lengths(i),
				p0(1) + cos(joint1)*link_lengths(i);
		link_segments.push_back(SegmentB(p0, p1));

		joint0 = joint1;
		p0 = p1;
	}

	return link_segments;
}


void FadbadPlanarSystem::get_limits(vecb<X_DIM>& x_min, vecb<X_DIM>& x_max, vecb<U_DIM>& u_min, vecb<U_DIM>& u_max) {
	x_min = this->x_min;
	x_max = this->x_max;
	u_min = this->u_min;
	u_max = this->u_max;
}

bdouble FadbadPlanarSystem::cost(const std::vector<vecb<X_DIM>, aligned_allocator<vecb<X_DIM>>>& X, const matb<X_DIM,X_DIM>& sigma0,
		const std::vector<vecb<U_DIM>, aligned_allocator<vecb<U_DIM>>>& U, const bdouble alpha) {
	bdouble cost = 0;
	int T = X.size();

	vecb<X_DIM> x_tp1 = vecb<X_DIM>::Zero();
	matb<X_DIM,X_DIM> sigma_t = sigma0, sigma_tp1 = matb<X_DIM,X_DIM>::Zero();
	for(int t=0; t < T-1; ++t) {
		belief_dynamics(X[t], sigma_t, U[t], alpha, x_tp1, sigma_tp1);
		cost += alpha_control*U[t].dot(U[t]);
		if (t < T-2) {
			cost += alpha_belief*sigma_tp1.trace();
		} else {
			cost += alpha_final_belief*sigma_tp1.trace();
		}
		sigma_t = sigma_tp1;
	}

	vecb<2> final_ee_position = get_link_segments(x_tp1).back().p1;
	vecb<2> goal_error = final_ee_position - x_tp1.segment<C_DIM>(J_DIM);
	cost += alpha_goal*goal_error.dot(goal_error);

	return cost;
}


/**
 * Private methods
 */

void FadbadPlanarSystem::linearize_dynfunc(const vecb<X_DIM>& x, const vecb<U_DIM>& u, const vecb<Q_DIM>& q, matb<X_DIM,X_DIM>& A, matb<X_DIM,Q_DIM>& M) {
	A = matb<X_DIM,X_DIM>::Identity();

	for(int j=0; j < Q_DIM; ++j) {
		for(int i=0; i < j+1; ++i) {
			M(i, j) = 1;
		}
	}
}

void FadbadPlanarSystem::linearize_obsfunc(const vecb<X_DIM>& x, const vecb<R_DIM>& r, matb<Z_DIM,X_DIM>& H, matb<Z_DIM,R_DIM>& N) {
	H = matb<Z_DIM,X_DIM>::Identity();
	N = matb<Z_DIM,R_DIM>::Identity();
}
