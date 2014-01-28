#include <sstream>
#include <string>

typedef double* experiment_t;
typedef std::vector< experiment_t > experiments_t;
experiments_t experiments;


template <size_t _dim, size_t _xDim, size_t _simXDim, size_t _uDim, size_t _zDim, size_t _kDim, size_t _sDim, size_t _bDim>
struct DynamicsBase {
	DynamicsBase(std::string _name) : name(_name), horizon(0), goalr(0.1), xParamReal(zeros<_kDim,1>()), Rint(identity<_uDim>()), Qint(identity<_xDim>()), QintVariance(identity<_xDim>()), QGoal(identity<_xDim>()), QGoalVariance(identity<_xDim>()), x0(zeros<_xDim,1>()), Sigma0(identity<_xDim>()), control_bounds(zeros<_uDim,2>()) {}

	/*
	virtual void saturateControl(Matrix<U_DIM>& u) const {
	if (u(0,0) > control_bounds(0,0) u(0,0) = control_bounds(0,0);
	if (u(0,0) < control_bounds(0,1) u(0,0) = control_bounds(0,1);
	}

	virtual void saturateControl(Matrix<U_DIM>& u) const {
	if (u(0,0) > control_bounds(0,0) u(0,0) = control_bounds(0,0);
	if (u(0,0) < control_bounds(0,1) u(0,0) = control_bounds(0,1);
	if (u(0,1) > control_bounds(1,0) u(0,1) = control_bounds(1,0);
	if (u(0,1) < control_bounds(1,1) u(0,1) = control_bounds(1,1);
	}
	*/

	virtual inline void linearizeDynamics(double step, const Matrix<_xDim>& xBar, const Matrix<_uDim>& uBar, Matrix<_xDim>& c, Matrix<_xDim, _xDim>& A, Matrix<_xDim, _uDim>& B, SymmetricMatrix<_xDim>& M, unsigned int flag) const = 0;

	virtual inline void linearizeObservation(double step, const Matrix<_xDim>& xBar, const Matrix<_uDim>& uBar, Matrix<_zDim, _xDim>& H, SymmetricMatrix<_zDim>& N) const = 0;

	virtual std::string getDirectoryName() const = 0;
	virtual void initExperiment(const experiment_t& e) = 0;
	virtual void saturateControl(Matrix<U_DIM>& u) const = 0;

	virtual inline SymmetricMatrix<Z_DIM> varN(const Matrix<_xDim>& x, const Matrix<_uDim>& u) const = 0;
	virtual inline SymmetricMatrix<Z_DIM> varN(const Matrix<_simXDim>& x, const Matrix<_uDim>& u) const = 0;

	virtual void setupExperiments(experiments_t& experiments) = 0;
	virtual void initControls(std::vector< Matrix<_uDim> >& u) const = 0;
	virtual inline Matrix<_xDim> f(double step, const Matrix<_xDim>& x, const Matrix<_uDim>& u) const = 0;
	virtual inline SymmetricMatrix<_xDim> varM(const Matrix<_xDim>& x, const Matrix<_uDim>& u) const = 0;

	virtual inline void quadratizeFinalCost(const Matrix<_xDim>& xBar, const SymmetricMatrix<_xDim>& SigmaBar, double& s, SymmetricMatrix<_xDim>& S, Matrix<1,_xDim>& sT, Matrix<1,S_DIM>& tT, unsigned int flag) const = 0;
	virtual inline bool quadratizeCost(double step, const Matrix<_xDim>& xBar, const SymmetricMatrix<_xDim>& SigmaBar, const Matrix<_uDim>& uBar, double& q, SymmetricMatrix<_xDim>& Q, 
		SymmetricMatrix<_uDim>& R, Matrix<_uDim, _xDim>& P, Matrix<1,_xDim>& qT, Matrix<1,_uDim>& rT, Matrix<1,S_DIM>& pT, unsigned int flag) const = 0;
	
	//virtual inline Matrix<_simXDim> g(const Matrix<_simXDim>& x, const Matrix<_uDim>& u, const float length) const = 0;
	virtual inline Matrix<_simXDim> freal(double step, const Matrix<_simXDim>& x,const Matrix<_uDim>& u) const = 0;
	virtual inline SymmetricMatrix<_simXDim> varMReal(const Matrix<_simXDim>& x, const Matrix<_uDim>& u) const = 0;

	std::string name;

	int horizon;

	SymmetricMatrix<_uDim> Rint;
	SymmetricMatrix<_xDim> Qint;
	SymmetricMatrix<_xDim> QintVariance;
	SymmetricMatrix<_xDim> QGoal;
	SymmetricMatrix<_xDim> QGoalVariance;

	Matrix<_xDim> x0;
	SymmetricMatrix<_xDim> Sigma0;

	Matrix<_kDim,1> xParamReal;

	Matrix<_xDim> xGoal;
	double goalr;

	Matrix<_uDim,2> control_bounds;

	experiments_t experiments;
};
