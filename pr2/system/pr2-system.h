#ifndef __PR2_SYSTEM_H__
#define __PR2_SYSTEM_H__

#include "pr2-sim.h"
#include "../utils/pr2-utils.h"
#include "../utils/utils.h"
#include "figtree.h"

#include <Eigen/Eigen>
#include <Eigen/StdVector>
using namespace Eigen;

#include "../../util/logging.h"

#define TIMESTEPS 5
#define DT 1.0 // Note: if you change this, must change the FORCES matlab file

#define X_DIM (ARM_DIM+3)	// arm + object
#define J_DIM ARM_DIM
#define U_DIM ARM_DIM
#define Q_DIM ARM_DIM
#define Z_DIM (ARM_DIM+3)	// relative object position
#define R_DIM (ARM_DIM+3)

#define TOTAL_VARS (TIMESTEPS*J_DIM + (TIMESTEPS-1)*U_DIM)

#define M_DIM 1000 // number of particles

typedef Matrix<double,X_DIM,1> VectorX;
typedef Matrix<double,J_DIM,1> VectorJ;
typedef Matrix<double,U_DIM,1> VectorU;
typedef Matrix<double,Q_DIM,1> VectorQ;
typedef Matrix<double,Z_DIM,1> VectorZ;
typedef Matrix<double,R_DIM,1> VectorR;
typedef Matrix<double,TOTAL_VARS,1> VectorTOTAL;
typedef Matrix<double,M_DIM,1> VectorM;

typedef Matrix<double,X_DIM,X_DIM> MatrixX;
typedef Matrix<double,J_DIM,J_DIM> MatrixJ;
typedef Matrix<double,U_DIM,U_DIM> MatrixU;
typedef Matrix<double,Q_DIM,Q_DIM> MatrixQ;
typedef Matrix<double,Z_DIM,Z_DIM> MatrixZ;
typedef Matrix<double,R_DIM,R_DIM> MatrixR;
typedef Matrix<double,TOTAL_VARS,TOTAL_VARS> MatrixTOTAL;

typedef Matrix<double,3,M_DIM> MatrixP;

typedef std::vector<VectorX, aligned_allocator<VectorX>> StdVectorX;
typedef std::vector<VectorJ, aligned_allocator<VectorJ>> StdVectorJ;
typedef std::vector<VectorU, aligned_allocator<VectorU>> StdVectorU;
typedef std::vector<Matrix4d, aligned_allocator<Matrix4d>> StdMatrix4d;

#include "voxel-grid.h" // needs to be here, not sure why

class ParticleGaussian;
class PR2System;

struct ParticleGaussian {
	Vector3d mean;
	Matrix3d cov;
	MatrixXd particles;
	double pct;
	Cube ODF;

	ParticleGaussian(const Vector3d& m, const Matrix3d& c, const MatrixXd& P, double p) :
		mean(m), cov(c), particles(P), pct(p) { };
};

/**
 * NOTE: all coordinates are with respect to OpenRAVE 'world' frame
 */
class PR2System {
	const double step = 0.0078125*0.0078125;
	const double INFTY = 1e10;

	const double alpha_control = 0; // .01
	const double alpha_belief = 1; // 1
	const double alpha_final_belief = 1; // 1
	const double alpha_goal = .5; // .5

public:
	PR2System(Vector3d& object);
	PR2System(Vector3d& object, Arm::ArmType arm_type, bool view);
	PR2System(Vector3d& object, Arm::ArmType arm_type, std::string env_file, std::string robot_name, bool view);

	VectorJ dynfunc(const VectorJ& j, const VectorU& u, const VectorQ& q, bool enforce_limits=false);
	VectorZ obsfunc(const VectorJ& j, const Vector3d& object, const VectorR& r);

	MatrixZ delta_matrix(const VectorJ& j, const Vector3d& object, const double alpha, const Cube& ODF);

	void belief_dynamics(const VectorX& x_t, const MatrixX& sigma_t, const VectorU& u_t, const double alpha, const Cube& ODF,
			VectorX& x_tp1, MatrixX& sigma_tp1);
	void execute_control_step(const VectorJ& j_t_real, const VectorJ& j_t, const VectorU& u_t, const MatrixP& P_t, const Cube& ODF,
				VectorJ& j_tp1_real, VectorJ& j_tp1, MatrixP& P_tp1);

	void get_limits(VectorJ& j_min, VectorJ& j_max, VectorU& u_min, VectorU& u_max);

	double cost(const StdVectorJ& J, const Vector3d& obj, const MatrixX& sigma0, const StdVectorU& U, const double alpha, const Cube& ODF);
	double cost_gmm(const StdVectorJ& J, const MatrixJ& j_sigma0, const StdVectorU& U,
				const std::vector<ParticleGaussian>& particle_gmm, const double alpha);
	void cost_gmm_and_grad(StdVectorJ& J, const MatrixJ& j_sigma0, StdVectorU& U,
				const std::vector<ParticleGaussian>& particle_gmm, const double alpha,
				double& cost, VectorTOTAL& grad);

//	VectorTOTAL cost_grad(StdVectorJ& J, const Vector3d& obj, const MatrixX& sigma0, StdVectorU& U, const double alpha);
	VectorTOTAL cost_gmm_grad(StdVectorJ& J, const MatrixJ& j_sigma0, StdVectorU& U,
			const std::vector<ParticleGaussian>& particle_gmm, const double alpha);

	/**
	 * RIPPED APART VERSIONS
	 */

	MatrixZ delta_matrix_ripped(const VectorJ& j, const Vector3d& object, const double alpha, const Matrix4d& sd_vec, const bool& obj_in_fov);
	void belief_dynamics_ripped(const VectorX& x_t, const MatrixX& sigma_t, const VectorU& u_t, const double alpha, const Matrix4d& sd_vec, const bool& obj_in_fov,
			VectorX& x_tp1, MatrixX& sigma_tp1);

	double cost_ripped(const StdVectorJ& J, const Vector3d& obj, const MatrixX& sigma0, const StdVectorU& U,
			const double alpha, const StdMatrix4d& sd_vecs, const std::vector<bool>& obj_in_fovs);
	double cost_gmm_ripped(const StdVectorJ& J, const MatrixJ& j_sigma0, const StdVectorU& U,
					const std::vector<ParticleGaussian>& particle_gmm, const double alpha,
					const std::vector<StdMatrix4d>& objs_sd_vecs, const std::vector<std::vector<bool> >& objs_in_fovs);
	VectorTOTAL cost_gmm_grad_ripped(StdVectorJ& J, const MatrixJ& j_sigma0, StdVectorU& U,
			const std::vector<ParticleGaussian>& particle_gmm, const double alpha);

	/**
	 * END RIPPED APART VERSIONS
	 */

	// use figtree
	void fit_gaussians_to_pf(const MatrixP& P, std::vector<ParticleGaussian>& particle_gmm);

	void update(const StdVector3d& new_pc, const Matrix<double,HEIGHT_FULL,WIDTH_FULL>& zbuffer, const Matrix4d& cam_pose);
	Cube get_ODF(const Vector3d& obj) { return vgrid->get_ODF(obj); }

	void display(const VectorJ& j, bool pause=true);
	void display(const StdVectorJ& J, bool pause=true);
	void display(const VectorJ& j, const std::vector<ParticleGaussian>& particle_gmm, bool pause=true);
	void display(const StdVectorJ& J, const std::vector<ParticleGaussian>& particle_gmm, bool pause=true);

	PR2* get_brett() { return brett; }
	Arm* get_arm() { return arm; }
	Camera* get_camera() { return cam; }
	VoxelGrid* get_voxel_grid() { return vgrid; }
	StdVector3d get_pc() { return pc; }

private:
	PR2* brett;
	Arm* arm;
	Camera* cam;
	Arm::ArmType arm_type;
	VoxelGrid* vgrid;

	Vector3d object;
	VectorJ j_min, j_max, u_min, u_max;
	MatrixQ Q;
	MatrixR R;

	StdVector3d pc;

	void init();

	void linearize_dynfunc(const VectorX& x, const VectorU& u, const VectorQ& q,
			Matrix<double,X_DIM,X_DIM>& A, Matrix<double,X_DIM,Q_DIM>& M);
	void linearize_obsfunc(const VectorX& x, const VectorR& r,
			Matrix<double,Z_DIM,X_DIM>& H);

	void update_particles(const VectorJ& j_tp1_t, const double delta_fov_real, const VectorZ& z_tp1_real, const MatrixP& P_t,
			MatrixP& P_tp1);
	double gauss_likelihood(const Vector3d& v, const Matrix3d& S);
	void low_variance_sampler(const MatrixP& P, const VectorM& W, MatrixP& P_sampled);

};

#endif
