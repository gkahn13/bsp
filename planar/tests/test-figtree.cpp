#include <iostream>
#include <stdio.h>
#include <string.h> // for memset
#include <math.h>   // for abs
#include "figtree.h"

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;
namespace np = boost::numpy;

#include <Eigen/Eigen>
#include <Eigen/StdVector>
using namespace Eigen;

template <size_t _dim0, size_t _dim1>
using mat = Matrix<double, _dim0, _dim1>;

template <size_t _dim>
using vec = Matrix<double, _dim, 1>;

#include "../include/planar-utils.h"

#define C_DIM 2
#define M_DIM 1000

inline double uniform(double low, double high) {
	return (high - low)*(rand() / double(RAND_MAX)) + low;
}

void find_mode_iter(const mat<C_DIM,M_DIM>& means, const mat<C_DIM,M_DIM>& particles, mat<C_DIM,M_DIM>& new_means) {
  // The dimensionality of each sample vector.
  int d = C_DIM;

  // The number of targets (vectors at which gauss transform is evaluated).
  int M = M_DIM;

  // The number of sources which will be used for the gauss transform.
  int N = M_DIM;

  // The bandwidth.  NOTE: this is not the same as standard deviation since
  // the Gauss Transform sums terms exp( -||x_i - y_j||^2 / h^2 ) as opposed
  // to  exp( -||x_i - y_j||^2 / (2*sigma^2) ).  Thus, if sigma is known,
  // bandwidth can be set to h = sqrt(2)*sigma.
  double h = sqrt(2)*3;

  // Desired maximum absolute error after normalizing output by sum of weights.
  // If the weights, q_i (see below), add up to 1, then this is will be the
  // maximum absolute error.
  // The smaller epsilon is, the more accurate the results will be, at the
  // expense of increased computational complexity.
  double epsilon = 1e-2; // 1e-2

  // The source array.  It is a contiguous array, where
  // ( x[i*d], x[i*d+1], ..., x[i*d+d-1] ) is the ith d-dimensional sample.
  // For example, below N = 20 and d = 7, so there are 20 rows, each
  // a 7-dimensional sample.
  double *x = new double[d*N];

  // The target array.  It is a contiguous array, where
  // ( y[j*d], y[j*d+1], ..., y[j*d+d-1] ) is the jth d-dimensional sample.
  // For example, below M = 10 and d = 7, so there are 10 rows, each
  // a 7-dimensional sample.
  double *y = new double[d*M];

  // Number of weights.  For each set of weights a different Gauss Transform is computed,
  // but by giving multiple sets of weights at once some overhead can be shared.
  int W = 3;
  // The weight array.  The ith weight is associated with the ith source sample.
  // To evaluate the Gauss Transform with the same sources and targets, but
  // different sets of weights, add another row of weights and set W = 2.
  double *q = new double[W*N];


  // source array is the particles themselves
//  std::cout << "source array:\n";
  for(int j=0; j < N; ++j) {
	  for(int i=0; i < d; ++i) {
		  x[j*d+i] = particles(i,j);
//		  std::cout << particles(i,j) << " ";
	  }
//	  std::cout << "\n";
  }
//  std::cout << "\n";

  // weights are (x,y) values of particles
  // plus last row of 1's to get the normalization
//  std::cout << "weight array:\n";
  for(int i=0; i < d; ++i) {
	  for(int j=0; j < N; ++j) {
		  q[i*N+j] = particles(i,j);
//		  std::cout << particles(i,j) << " ";
	  }
//	  std::cout << "\n";
  }
  for(int j=0; j < N; ++j) {
	  q[d*N+j] = 1;
//	  std::cout << q[d*N+j] << " ";
  }
//  std::cout << "\n\n";


  // target array is the current means
//  std::cout << "target array:\n";
  for(int j=0; j < M; ++j) {
	  for(int i=0; i < d; ++i) {
		  y[j*d+i] = means(i,j);
//		  std::cout << means(i,j) << " ";
	  }
//	  std::cout << "\n";
  }
//  std::cout << "\n";

  // allocate array into which the result of the Gauss Transform will be stored for each
  // target sample.  The first M elements will correspond to the Gauss Transform computed
  // with the first set of weights, second M elements will correspond to the G.T. computed
  // with the second set of weights, etc.
  double * output = new double[W*M];

  // initialize all output arrays to zero
  memset( output        , 0, sizeof(double)*W*M );
  //
  // RECOMMENDED way to call figtree().
  //
  // Evaluates the Gauss transform using automatic method selection (the automatic
  // method selection function analyzes the inputs -- including source/target
  // points, weights, bandwidth, and error tolerance -- to automatically choose
  // between FIGTREE_EVAL_DIRECT, FIGTREE_EVAL_DIRECT_TREE, FIGTREE_EVAL_IFGT,
  // FIGTREE_EVAL_IFGT_TREE.
  // This function call makes use of the default parameters for the eval method
  // and param method, and is equivalent to
  // figtree( d, N, M, W, x, h, q, y, epsilon, g_auto, FIGTREE_EVAL_AUTO, FIGTREE_PARAM_NON_UNIFORM ).
  figtree( d, N, M, W, x, h, q, y, epsilon, output );

  for(int j=0; j < M; ++j) {
	  for(int i=0; i < d; ++i) {
		  new_means(i,j) = output[i*M+j] / output[d*M+j];
	  }
//	  std::cout << output[d*M+j] << " ";
  }
//  std::cout << "\n";



  // deallocate memory
  delete [] x;
  delete [] y;
  delete [] q;
  delete [] output;
}

void display(const mat<C_DIM,M_DIM>& P, const std::vector<vec<C_DIM>,aligned_allocator<vec<C_DIM>>>& modes) {
	Py_Initialize();
	np::initialize();
	py::numeric::array::set_module_and_type("numpy", "ndarray");

	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string planar_dir = bsp_dir + "bsp/planar/tests";

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+planar_dir+"')"), main_namespace);
	py::object plot_mod = py::import("plot_mean_shift");

	py::object plot_mean_shift = plot_mod.attr("plot_mean_shift");

	np::ndarray P_ndarray = planar_utils::eigen_to_ndarray(P);

	py::list modes_pylist;
	for(int i=0; i < modes.size(); ++i) {
		modes_pylist.append(planar_utils::eigen_to_ndarray(modes[i]));
	}

	plot_mean_shift(P_ndarray, modes_pylist);
}

int main() {
	srand(time(0));

	mat<C_DIM,M_DIM> particles, means, new_means;

	boost::mt19937 igen;
	boost::normal_distribution<> nd(0.0, 1.0);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > gen(igen, nd);


	for(int j=0; j < M_DIM; ++j) {
//		if (j < M_DIM/2) {
//			particles(0, j) = gen() - 4;
//			particles(1, j) = gen();
//		} else {
//			particles(0, j) = gen() + 4;
//			particles(1, j) = gen();
//		}

		particles(0,j) = uniform(-10,10);
		particles(1,j) = uniform(0,10);

//		if (j < M_DIM/2) {
//			particles(0,j) = uniform(-10,0);
//			particles(1,j) = uniform(0,10);
//		} else {
//			particles(0,j) = uniform(0,10);
//			particles(1,j) = uniform(0,10);
//		}
	}

	means = particles;

//	std::cout << "initial means:\n" << means.block<C_DIM,10>(0,0) << "\n\n";

	int d = C_DIM;
	int M = M_DIM;
	int N = M_DIM;
	double h = sqrt(2)*4;
	double epsilon = .1; // 1e-2
	double *x = new double[d*N];
	double *y = new double[d*M];
	int W = 3;
	double *q = new double[W*N];
	double * output = new double[W*M];

	for(int j=0; j < N; ++j) {
		for(int i=0; i < d; ++i) {
			x[j*d+i] = particles(i,j);
			y[j*d+i] = particles(i,j);
		}
	}

	for(int i=0; i < d; ++i) {
		for(int j=0; j < N; ++j) {
			q[i*N+j] = particles(i,j);
		}
	}
	for(int j=0; j < N; ++j) {
		q[d*N+j] = 1;
	}

//	int eval_method;
//	figtreeChooseEvaluationMethod(d, N, M, W, x, h, y, epsilon, FIGTREE_PARAM_NON_UNIFORM, 0, &eval_method);
//	std::cout << "eval_method: " << eval_method << "\n";

	util::Timer total_timer, figtree_timer;
	double figtree_time = 0;
	util::Timer_tic(&total_timer);

	double max_diff = INFINITY;
	while(max_diff > 1e-3) {
		  for(int j=0; j < M; ++j) {
			  for(int i=0; i < d; ++i) {
				  y[j*d+i] = means(i,j);
			  }
		  }

		  memset( output        , 0, sizeof(double)*W*M );

		  util::Timer_tic(&figtree_timer);
		  figtree( d, N, M, W, x, h, q, y, epsilon, output ); // .45

//		  figtree( d, N, M, W, x, h, q, y, epsilon, output, FIGTREE_EVAL_DIRECT ); // > 5
//		  figtree( d, N, M, W, x, h, q, y, epsilon, output, FIGTREE_EVAL_DIRECT_TREE ); // > 10
//		  figtree( d, N, M, W, x, h, q, y, epsilon, output, FIGTREE_EVAL_IFGT, FIGTREE_PARAM_UNIFORM, 1 ); // .5
//		  figtree( d, N, M, W, x, h, q, y, epsilon, output, FIGTREE_EVAL_IFGT, FIGTREE_PARAM_NON_UNIFORM, 1 ); // .4
//		  figtree( d, N, M, W, x, h, q, y, epsilon, output, FIGTREE_EVAL_IFGT_TREE, FIGTREE_PARAM_UNIFORM, 1 ); // 1
//		  figtree( d, N, M, W, x, h, q, y, epsilon, output, FIGTREE_EVAL_IFGT_TREE, FIGTREE_PARAM_NON_UNIFORM, 1 ); // .45
		  figtree_time += util::Timer_toc(&figtree_timer);

		  for(int j=0; j < M; ++j) {
			  for(int i=0; i < d; ++i) {
				  new_means(i,j) = output[i*M+j] / output[d*M+j];
			  }
		  }


//		find_mode_iter(means, particles, new_means);

		max_diff = (means - new_means).colwise().norm().maxCoeff();
//		std::cout << "max_diff: " << max_diff << "\n";
		means = new_means;

//		std::cout << "means[" << t << "]:\n" << means.block<C_DIM,10>(0,0) << "\n\n";
//		std::cin.ignore();
	}

//	std::cout << "final means:\n" << means.block<C_DIM,10>(0,0) << "\n\n";

	std::vector<vec<C_DIM>, aligned_allocator<vec<C_DIM>>> modes;
	for(int i=0; i < M_DIM; ++i) {
		bool is_new_mode = true;
		for(int j=0; j < modes.size(); ++j) {
			if ((means.col(i)-modes[j]).norm() < .05) {
				is_new_mode = false;
				break;
			}
		}

		if (is_new_mode) {
			modes.push_back(means.col(i));
		}
	}

	double total_time = util::Timer_toc(&total_timer);

//	std::cout << "total time: " << total_time << "\n";
	std::cout << "figtree time: " << figtree_time << "\n\n";

//	std::cout << "number of modes: " << modes.size() << "\n";
//	for(int i=0; i < modes.size(); ++i) {
//		std::cout << modes[i].transpose() << "\n";
//	}
//	std::cout << "\n";

	display(particles, modes);
}
