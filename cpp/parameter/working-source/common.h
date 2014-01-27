#ifndef __COMMON_H__
#define __COMMON_H__

#include "callisto.h"

#include "matrix.h"
#include "utils.h"

#include "POMDP-On4.h"

#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <typeinfo>
#include <sstream>
#include <direct.h> // Needed for _mkdir
#include <process.h>
#include <limits>

#define DT 1.0/100.0

#define SIMULATE
#define SYSTEM SIMULATE

/* The following is related to displaying. */
struct Obstacle {
	std::vector< Matrix<2> > vlist;
};

std::vector<Obstacle> obslist;

int *cal_path, *cal_beliefs;
int cal_environment, cal_obstacles, cal_box, cal_goal, cal_sensor, cal_cylinder, cal_point, cal_rrt;
//double goalr;

int cal_line, cal_test, cal_poly;

// Forward declarations
struct Planner;
/*
typedef std::vector<Planner*> planners_t;
planners_t planners;

typedef std::vector< experiment_t > experiments_t;
experiments_t experiments;
*/
struct bar_t {
	std::vector< Matrix<X_DIM> > x;
	std::vector< SymmetricMatrix<X_DIM> > sigma;
	std::vector< Matrix<U_DIM, X_DIM> > L;
	std::vector< Matrix<U_DIM> > u;

	void resize(const int steps, const Matrix<X_DIM>& xDefault, const SymmetricMatrix<X_DIM>& sigmaDefault, const Matrix<U_DIM>& uDefault, const Matrix<U_DIM, X_DIM>& LDefault) {
		x.clear();
		sigma.clear();
		u.clear();
		L.clear();

		x.resize(steps+1, xDefault);
		sigma.resize(steps+1, sigmaDefault);
		u.resize(steps, uDefault);
		L.resize(steps, LDefault);
	}
};

struct hat_t {
	std::vector< Matrix<X_DIM> > x;
	std::vector< SymmetricMatrix<X_DIM> > sigma;
	std::vector< Matrix<U_DIM> > u;

	void resize(const int steps, const Matrix<X_DIM>& xDefault, const SymmetricMatrix<X_DIM>& sigmaDefault, const Matrix<U_DIM>& uDefault) {
		x.clear();
		sigma.clear();
		u.clear();

		x.resize(steps+1, xDefault);
		sigma.resize(steps+1, sigmaDefault);
		u.resize(steps, uDefault);
	}
};

struct metrics_t {
	metrics_t() : totalCost(0.0) {
		likelihoods = std::vector< Matrix<K_DIM,K_DIM> >();
	}

	double totalCost;
	std::vector< Matrix<K_DIM,K_DIM> > likelihoods;
	double runtime;
};

#include "systems.hpp"
#include "planners.hpp"

void saveMetrics(const planners_t& planners, const std::string& directory, const std::string& name) {
	std::string metrics_file = "stats/" + directory;

	int err = _mkdir(metrics_file.c_str());
	if (err && errno != EEXIST) {
		std::cout << "Could not make directory: " << metrics_file << std::endl;
		exit(1);
	}

	// Name and open file
	metrics_file += "/" + name + "_metrics.txt";
	std::ofstream metrics_fh;
	metrics_fh.open(metrics_file, std::ios::out);

	// Write header
	metrics_fh << "Runtime\tTotal Cost\tFinal Likelihood\t" << std::endl;

	// Write metrics
	size_t num_plans = planners.size();
	for (size_t i = 0; i < num_plans; ++i) {
		metrics_fh << planners[i]->metrics.runtime << '\t' << planners[i]->metrics.totalCost << '\t';

		for (int j = 0; j < K_DIM; ++j) {
			metrics_fh << planners[i]->metrics.likelihoods[planners[i]->metrics.likelihoods.size()-1][j];
		}

		metrics_fh << std::endl;
	}

	// Save and close file
	metrics_fh.flush();
	metrics_fh.close();
}

void saveData(const planners_t& planners, const std::string& directory, const std::string& name) {
	std::string data_file = "stats/" + directory;

	int err = _mkdir(data_file.c_str());
	if (err && errno != EEXIST) {
		std::cout << "Could not make directory: " << data_file << std::endl;
		exit(1);
	}

	// Name and open file
	data_file += "/" + name + "_data.txt";
	std::ofstream data_fh;
	data_fh.open(data_file, std::ios::out);

	// Write header
	size_t num_plans = planners.size();
	size_t iters = planners[0]->hat.u.size();
	for (size_t k = 0; k < iters+1; ++k) {
		for (int i = 0; i < X_DIM; ++i) {
			data_fh << "x[" << k << ':' << i << "]" << '\t';
		}

		for (int i = 0; i < X_DIM; ++i) {
			for (int j = 0; j < X_DIM; ++j) {
				data_fh << "sigma(" << k << ':' << i << ',' << j << ")" << '\t';
			}
		}

		for (int i = 0; i < K_DIM; ++i) {
			data_fh << "likelihood[" << k << ':' << i << "]\t";
		}

		for (int i = 0; i < U_DIM-1; ++i) {
			data_fh << "u[" << k << ':' << i << "]" << '\t';
		}
		data_fh << "u[" << k << ':' << (U_DIM-1) << "]";
	}

	data_fh << std::endl;

	// Iterate over data and write to file
	for (size_t plan = 0; plan < num_plans; ++plan) {
		for(size_t i = 0; i < iters; ++i) {
			for (int j = 0; j < X_DIM; ++j) {
				data_fh << planners[plan]->hat.x[i][j] << '\t';
			}

			for (int j = 0; j < X_DIM; ++j) {
				for (int k = 0; k < X_DIM; ++k) {
					data_fh << planners[plan]->hat.sigma[i](j,k) << '\t';
				}
			}

			for (int j = 0; j < K_DIM; ++j) {
				data_fh << planners[plan]->metrics.likelihoods[i][j] << '\t';
			}

			for (int j = 0; j < U_DIM; ++j) {
				data_fh << planners[plan]->hat.u[i][j] << '\t';
			}
		}

		// There is one more belief than there are controls so write the final belief to file
		for (int j = 0; j < X_DIM; ++j) {
			data_fh << planners[plan]->hat.x[iters][j] << '\t';
		}

		for (int j = 0; j < X_DIM; ++j) {
			for (int k = 0; k < X_DIM; ++k) {
				data_fh << planners[plan]->hat.sigma[iters](j,k) << '\t';
			}
		}

		for (int j = 0; j < K_DIM; ++j) {
			data_fh << planners[plan]->metrics.likelihoods[iters][j] << '\t';
		}

		for (int j = 0; j < U_DIM-1; ++j) {
			data_fh << "0\t";
		}

		data_fh << '0';

		data_fh << std::endl;
	}

	// Save and close file
	data_fh.flush();
	data_fh.close();
}

void gatherData(Dynamics& experiment, System& sys) {
	srand(time(NULL));

	// Setup experiments
	experiments_t experiments;

	experiment.setupExperiments(experiments);

	// Setup planners
	planners_t planners;
	/*
	Random random(experiment);
	random.setName("random");
	planners.push_back(&random);
	*/
	NoReplan noreplan(experiment);
	noreplan.setName("noreplan");
	planners.push_back(&noreplan);
	/*
	Replan replanwoiters(experiment);
	replanwoiters.subIters = 1;
	replanwoiters.setName("replanwoiters");
	planners.push_back(&replanwoiters);
	
	Replan replanwiters(experiment);
	replanwiters.setName("replanwiters");
	planners.push_back(&replanwiters);
	*/
	std::vector< Matrix<U_DIM> > u(experiment.horizon, zeros<U_DIM,1>());
	experiment.initControls(u);

	for (experiments_t::iterator e = experiments.begin(); e != experiments.end(); ++e) {
		experiment.initExperiment(*e);

		for (planners_t::iterator p = planners.begin(); p != planners.end(); ++p) {
			Planner* planner = *p;
			planner->plan(sys, u);
			planner->saveLogs();
		}
	}
}

void runMPC(const Dynamics& e, System& sys) 
{
	srand(0);
	std::vector< Matrix<U_DIM> > u(e.horizon, zeros<U_DIM,1>());

	Planner* planner = new ReplanMPC(e,-1,-1);
	//Planner* planner = new Replan(e,-1,1);
	//Planner* planner = new Replan(e,-1,-1);
	e.initControls(u);
	planner->plan(sys, u);
}


// TODO DWEBB -- this will need to be generalized to handle different numbers of parameters
void analyzePlanner(const Dynamics& e, System& sys, size_t iters = 30) {
	std::vector< Matrix<U_DIM> > u(e.horizon, zeros<U_DIM,1>());

	time_t t;
	time(&t);
	for (int k = 3; k > -1; k--) {
		//srand(t); // Seed the random number generator so that all noise vectors are the same between experiments.
		srand(0);

		// Setup stats
		size_t last_state = e.horizon;
		size_t num_data_points = last_state + 1;
		typedef std::pair<double, double> stats_t; // Mean and variance
		typedef std::pair< stats_t, stats_t > stat_set_t; // stats for mass estimate and variance estimate
		typedef std::vector< stat_set_t > stat_sets_t;
		stat_sets_t sss(num_data_points, std::make_pair(std::make_pair(0,0),std::make_pair(0,0)));

		//stat_set_t temp;

		// Setup planners
		std::string name;
		switch (k) {
		case 0:
			name = "random";
			break;
		case 1:
			name = "noreplan";
			break;
		case 2:
			name = "replanwoiters";
			break;
		case 3:
			name = "replanwiters";
			break;
		}
		planners_t planners(iters);
		for (size_t i = 0; i < iters; ++i) {
			switch(k) {
			case 0:
				planners[i] = new Random(e);
				std::cout << "Random instantiated" << std::endl;
				break;
			case 1:
				planners[i] = new NoReplan(e);
				std::cout << "NoReplan instantiated" << std::endl;
				break;
			case 2:
				planners[i] = new Replan(e,-1,1);
				std::cout << "Replan instantiated" << std::endl;
				break;
			case 3:
				planners[i] = new Replan(e,-1,-1);
				break;
			}
			e.initControls(u);
			planners[i]->plan(sys, u);
		}

		// Save estimates of mass and variance
		float curr_val = 0.0f;
		for (size_t i = 0; i < num_data_points; ++i) {
			for (size_t j = 0; j < iters; ++j) {
				curr_val = planners[j]->hat.x[i][X_DIM-1];
				sss[i].first.first += curr_val;
				sss[i].first.second += curr_val*curr_val;

				curr_val = planners[j]->hat.sigma[i](X_DIM-1,X_DIM-1);
				sss[i].second.first += curr_val;
				sss[i].second.second += curr_val*curr_val;
			}

			sss[i].first.first /= iters;
			sss[i].first.second /= iters;
			sss[i].first.second -= sss[i].first.first*sss[i].first.first;

			sss[i].second.first /= iters;
			sss[i].second.second /= iters;
			sss[i].second.second -= sss[i].second.first*sss[i].second.first;
		}

		std::string dir = e.getDirectoryName();
		saveMetrics(planners, dir, name);
		saveData(planners, dir, name);
	}
}

#endif //__COMMON_H__