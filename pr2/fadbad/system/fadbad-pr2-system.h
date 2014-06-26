#ifndef __FADBAD_PR2_SYSTEM_H__
#define __FADBAD_PR2_SYSTEM_H__

#include "../fadbad-utils.h"
#include "fadbad-pr2-sim.h"

#include "../../system/pr2-system.h"

#include <Eigen/Eigen>
#include <Eigen/StdVector>
using namespace Eigen;

//#include "../../../util/logging.h"

typedef Matrix<bdouble,X_DIM,1> VectorXb;
typedef Matrix<bdouble,J_DIM,1> VectorJb;
typedef Matrix<bdouble,U_DIM,1> VectorUb;
typedef Matrix<bdouble,Q_DIM,1> VectorQb;
typedef Matrix<bdouble,Z_DIM,1> VectorZb;
typedef Matrix<bdouble,R_DIM,1> VectorRb;
typedef Matrix<bdouble,TOTAL_VARS,1> VectorTOTALb;
typedef Matrix<bdouble,M_DIM,1> VectorMb;

typedef Matrix<bdouble,X_DIM,X_DIM> MatrixXb;
typedef Matrix<bdouble,J_DIM,J_DIM> MatrixJb;
typedef Matrix<bdouble,U_DIM,U_DIM> MatrixUb;
typedef Matrix<bdouble,Q_DIM,Q_DIM> MatrixQb;
typedef Matrix<bdouble,Z_DIM,Z_DIM> MatrixZb;
typedef Matrix<bdouble,R_DIM,R_DIM> MatrixRb;

typedef Matrix<bdouble,3,M_DIM> MatrixPb;

typedef std::vector<VectorXb, aligned_allocator<VectorXb>> StdVectorXb;
typedef std::vector<VectorJb, aligned_allocator<VectorJb>> StdVectorJb;
typedef std::vector<VectorUb, aligned_allocator<VectorUb>> StdVectorUb;

struct FadbadParticleGaussian {
	Vector3b mean;
	Matrix3b cov;
	bdouble pct;

	FadbadParticleGaussian(const Vector3b& m, const Matrix3b& c, bdouble p) :
		mean(m), cov(c), pct(p) { };
};

class FadbadPR2System {
	const bdouble step = 0.0078125*0.0078125;
	const bdouble INFTY = 1e10;

	const bdouble alpha_control = 0; // .01
	const bdouble alpha_belief = 1; // 1
	const bdouble alpha_final_belief = 1; // 1
	const bdouble alpha_goal = 0; // .5

public:
	FadbadPR2System(PR2* brett, Arm* arm, Camera* cam);

	VectorJb dynfunc(const VectorJb& j, const VectorUb& u, const VectorQb& q, bool enforce_limits=false);
	VectorZb obsfunc(const VectorJb& j, const Vector3b& object, const VectorRb& r);

	MatrixZb delta_matrix(const VectorJb& j, const Vector3b& object, const bdouble alpha);

	void belief_dynamics(const VectorXb& x_t, const MatrixXb& sigma_t, const VectorUb& u_t, const bdouble alpha,
			VectorXb& x_tp1, MatrixXb& sigma_tp1);

	void get_limits(VectorJb& j_min, VectorJb& j_max, VectorUb& u_min, VectorUb& u_max);

	bdouble cost(const StdVectorJb& J, const Vector3b& obj, const MatrixXb& sigma0, const StdVectorUb& U, const bdouble alpha);
	bdouble cost_gmm(const StdVectorJb& J, const MatrixJb& j_sigma0, const StdVectorUb& U,
				const std::vector<ParticleGaussian>& particle_gmm, const bdouble alpha);


private:
	PR2* brett;
	Arm* arm;
	Camera* cam;

	Vector3b object;
	VectorJb j_min, j_max, u_min, u_max;
	MatrixQb Q;
	MatrixRb R;

	void init();

	void linearize_dynfunc(const VectorXb& x, const VectorUb& u, const VectorQb& q,
			Matrix<bdouble,X_DIM,X_DIM>& A, Matrix<bdouble,X_DIM,Q_DIM>& M);
	void linearize_obsfunc(const VectorXb& x, const VectorRb& r,
			Matrix<bdouble,Z_DIM,X_DIM>& H);

};

#endif
