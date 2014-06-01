typedef std::vector<Planner*> planners_t;
planners_t planners;


struct Planner {
	Planner(const Dynamics& e) : step(DT), experiment(e) {
		name = "unknown";
	}

	void setName(const std::string& n) {
		name = n;
	}

	void plan(System& sys, const std::vector< Matrix<U_DIM> >& u) {
		clock_t start_t = clock();

		_plan(sys, u);

		clock_t end_t = clock();

		metrics.runtime = ((double) (end_t - start_t) / CLOCKS_PER_SEC);

		std::cout << " Time: " << metrics.runtime << " " << CLOCKS_PER_SEC;
		std::cout << std::endl;
	}

	virtual void _plan(System& sys, const std::vector< Matrix<U_DIM> >& u) = 0;

	/*
	This function is supposed to take an initial belief state and a control policy (uBar, L)
	to determine the inputs for solvePOMDP. solvePOMDP assumes
	*/
	void evalControls(double step, const Matrix<X_DIM>& xHat, bar_t& bar) {
		int steps = bar.u.size();

		Matrix<X_DIM> xCurrent, xNext;
		xCurrent = xHat;

		for (int t = 0; t < steps; ++t) {
			bar.u[t] = bar.L[t]*(xCurrent - bar.x[t]) + bar.u[t];
			xNext = experiment.f(step, xCurrent, bar.u[t]);
			xCurrent = xNext;
		}
	}

	void calculateMetrics(const hat_t& hat, metrics_t& metrics) {
		int numSteps = hat.u.size();
		int numBeliefs = numSteps + 1;

		// Add up intermediate costs
		for (int i = 0; i < numSteps; ++i) {
			metrics.totalCost += scalar((~(hat.x[i]-experiment.xGoal))*experiment.Qint*(hat.x[i]-experiment.xGoal)) + scalar((~hat.u[i])*experiment.Rint*hat.u[i]) + scalar(vecTh(experiment.QintVariance)*vec(hat.sigma[i]));
		}

		// Add on final cost
		metrics.totalCost += scalar((~(hat.x[numSteps]-experiment.xGoal))*experiment.QGoal*(hat.x[numSteps]-experiment.xGoal)) + scalar(vecTh(experiment.QGoalVariance)*vec(hat.sigma[numSteps]));

		// Calculate likelihood for each belief
		metrics.likelihoods.resize(numBeliefs);
		for (int i = 0; i < numBeliefs; ++i) {
			Matrix<K_DIM,1> xParamEst = hat.x[i].subMatrix<K_DIM>(X_DIM-K_DIM,0);
			SymmetricMatrix<K_DIM> SigmaParamEst = hat.sigma[i].subSymmetricMatrix<K_DIM>(X_DIM-K_DIM);

			Matrix<K_DIM,1> xDiff = xParamEst - experiment.xParamReal;
			Matrix<1,K_DIM> xDiffT = ~xDiff;
			// TODO DWEBB fix this for K_DIM > 1
			//metrics.likelihoods[i] = exp(-((xDiffT/SigmaParamEst)*xDiff)/2)/(sqrt(det(SigmaParamEst)*pow(2*M_PI,K_DIM)));
		}
	}

	virtual void saveLogs() {
		_saveMetrics();
		_saveData();
	}

	std::string directory() {
		std::string dir = "logs/" + experiment.getDirectoryName();

		int err = _mkdir(dir.c_str());
		if (err && errno != EEXIST) {
			std::cout << "Could not make directory: " << dir << std::endl;
			exit(1);
		}

		return dir;
	}

	virtual void _saveMetrics() {
		// Name and open file
		std::string metrics_file = directory();
		metrics_file += "/" + name;
		metrics_file += "_metrics.txt";
		std::ofstream metrics_fh;
		metrics_fh.open(metrics_file, std::ios::out);

		// Write metrics
		metrics_fh << "Runtime\tTotal Cost\tFinal Likelihood" << std::endl;
		metrics_fh << metrics.runtime << '\t' << metrics.totalCost << '\t';

		for (int i = 0; i < K_DIM; ++i) {
			metrics_fh << metrics.likelihoods[metrics.likelihoods.size()-1][i];
		}

		// Save and close file
		metrics_fh.flush();
		metrics_fh.close();
	}

	virtual void _saveData() {
		// Name and open file
		std::string data_file = directory();
		data_file += "/" + name;
		data_file += "_data.txt";
		std::ofstream data_fh;
		data_fh.open(data_file, std::ios::out);

		// Write header
		for (int i = 0; i < X_DIM; ++i) {
			data_fh << "x[" << i << "]" << '\t';
		}

		for (int i = 0; i < X_DIM; ++i) {
			for (int j = 0; j < X_DIM; ++j) {
				data_fh << "sigma(" << i << ',' << j << ")" << '\t';
			}
		}

		for (int i = 0; i < K_DIM; ++i) {
			data_fh << "likelihood[" << i << "]\t";
		}

		for (int i = 0; i < U_DIM; ++i) {
			data_fh << "u[" << i << "]" << '\t';
		}

		data_fh << std::endl;

		// Iterate over data and write to file
		int iters = (int)hat.u.size();
		for(int i = 0; i < iters; ++i) {
			for (int j = 0; j < X_DIM; ++j) {
				data_fh << hat.x[i][j] << '\t';
			}

			for (int j = 0; j < X_DIM; ++j) {
				for (int k = 0; k < X_DIM; ++k) {
					data_fh << hat.sigma[i](j,k) << '\t';
				}
			}

			for (int j = 0; j < K_DIM; ++j) {
				data_fh << metrics.likelihoods[i][j] << '\t';
			}

			for (int j = 0; j < U_DIM; ++j) {
				data_fh << hat.u[i][j] << '\t';
			}

			data_fh << std::endl;
		}

		// There is one more belief than there are controls so write the final belief to file
		for (int j = 0; j < X_DIM; ++j) {
			data_fh << hat.x[iters][j] << '\t';
		}

		for (int j = 0; j < X_DIM; ++j) {
			for (int k = 0; k < X_DIM; ++k) {
				data_fh << hat.sigma[iters](j,k) << '\t';
			}
		}

		for (int j = 0; j < K_DIM; ++j) {
			data_fh << metrics.likelihoods[iters][j] << '\t';
		}

		for (int j = 0; j < U_DIM; ++j) {
			data_fh << "0\t";
		}

		// Save and close file
		data_fh.flush();
		data_fh.close();
	}

	std::string name;
	bar_t bar;
	hat_t hat;
	metrics_t metrics;
	double step;

	const Dynamics& experiment;
};

struct Random : Planner {
	Random(const Dynamics& e) : Planner(e) {}

	void _plan(System& sys, const std::vector< Matrix<U_DIM> >& u) {
		int steps = u.size();
		//bar.resize(steps, zeros<X_DIM,1>(), zeros<X_DIM>(), zeros<U_DIM,1>(), zeros<U_DIM,X_DIM>());
		bar.resize(steps, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>(), zeros<U_DIM,X_DIM>());
		bar.u = u;
		//hat.resize(steps, zeros<X_DIM,1>(), zeros<X_DIM>(), zeros<U_DIM,1>());
		hat.resize(steps, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>());

		sys.executeControls(bar, hat);

		calculateMetrics(hat, metrics);
	}
};

struct NoReplan : Planner {
	NoReplan(const Dynamics& e) : Planner(e) {}

	void _plan(System& sys, const std::vector< Matrix<U_DIM> >& u) {
		int steps = u.size();
		//bar.resize(steps, zeros<X_DIM,1>(), zeros<X_DIM>(), zeros<U_DIM,1>(), zeros<U_DIM,X_DIM>());
		bar.resize(steps, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>(), zeros<U_DIM,X_DIM>());
		bar.u = u;
		//hat.resize(steps, zeros<X_DIM,1>(), zeros<X_DIM>(), zeros<U_DIM,1>());
		hat.resize(steps, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>());

		//solvePOMDP(experiment, step, linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, bar.x, bar.sigma, bar.u, bar.L);
		solvePOMDP<X_DIM, U_DIM, Z_DIM>(experiment, step, bar.x, bar.sigma, bar.u, bar.L);

		sys.executeControls(bar, hat);

		for (int i = 0; i < steps; ++i) {
			bar.L[i] = zeros<U_DIM,X_DIM>();
		}

		calculateMetrics(hat, metrics);
	}
};

struct Replan : Planner {
	Replan(const Dynamics& e) : Planner(e), initIters(-1), subIters(-1) {}
	Replan(const Dynamics& e, int iI, int sI) : Planner(e), initIters(iI), subIters(sI) {}

	void _plan(System& sys, const std::vector< Matrix<U_DIM> >& u) {
		int steps = u.size();

		// Resize all the buckets - make the bar buckets one larger than necessary to accomadate loop structure below.
		// This prevents the need to check whether the buckets are too small through each pass of the loop.
		bar.resize(steps+1, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>(), zeros<U_DIM,X_DIM>());
		std::vector< Matrix<SIM_X_DIM> > xReal(steps+1, experiment.x0.subMatrix<SIM_X_DIM>(0,0));
		hat.resize(steps, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>());

		bar.u = u;
		bar.u.resize(steps+1, zeros<U_DIM,1>());

		int iterations = initIters;
		for (int t = 0; t < steps; t++) {
			std::cout << "Step: " << t << std::endl;

			// Convert (uBar, L) tuple to uBar by evaluating the current control policy
			evalControls(step, hat.x[t], bar);

			// Shorten all bar buckets to accomodate shortening horizon
			bar.x.resize(bar.x.size()-1);
			bar.x[0] = hat.x[t];

			bar.sigma.resize(bar.sigma.size()-1);
			bar.sigma[0] = hat.sigma[t];

			reverse(bar.u.begin(),bar.u.end());
			bar.u.resize(bar.u.size()-1);
			reverse(bar.u.begin(),bar.u.end());

			bar.L.resize(bar.L.size()-1);

			// Solve the POMDP
			//solvePOMDP(experiment, step, linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, bar.x, bar.sigma, bar.u, bar.L, iterations);
			solvePOMDP<X_DIM, U_DIM, Z_DIM>(experiment, step, bar.x, bar.sigma, bar.u, bar.L, iterations);

			// Simulate one step to get new belief
			sys.executeControlStep(bar.x[0], bar.u[0], bar.L[0], xReal[t], xReal[t+1], hat.x[t], hat.sigma[t], hat.x[t+1], hat.sigma[t+1], hat.u[t]);

			iterations = subIters;
		}

		evalControls(step, hat.x[steps], bar);
		bar.L.resize(steps, zeros<U_DIM,X_DIM>());

		calculateMetrics(hat, metrics);
	}

	int initIters;
	int subIters;
};

struct ReplanMPC : Planner {
	ReplanMPC(const Dynamics& e) : Planner(e), initIters(-1), subIters(-1) {}
	ReplanMPC(const Dynamics& e, int iI, int sI) : Planner(e), initIters(iI), subIters(sI) {}

	void _plan(System& sys, const std::vector< Matrix<U_DIM> >& u) {
		int steps = u.size();
		int simsteps = 1000;

		// Resize all the buckets - make the bar buckets one larger than necessary to accomadate loop structure below.
		// This prevents the need to check whether the buckets are too small through each pass of the loop.
		bar.resize(steps, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>(), zeros<U_DIM,X_DIM>());
		
		std::vector< Matrix<SIM_X_DIM> > xReal(simsteps+1, experiment.x0.subMatrix<SIM_X_DIM>(0,0));
		
		hat.resize(simsteps, experiment.x0, experiment.Sigma0, zeros<U_DIM,1>());

		bar.u = u;
		bar.u.resize(steps, zeros<U_DIM,1>());

		int iterations = initIters;
		for (int t = 0; t < simsteps; t++) {
			
			std::cout << "Step: " << t << std::endl;

			// Convert (uBar, L) tuple to uBar by evaluating the current control policy
			evalControls(step, hat.x[t], bar);

			bar.x[0] = hat.x[t];
			bar.sigma[0] = hat.sigma[t];

			// Solve the POMDP
			//solvePOMDP(experiment, step, linearizeDynamics, linearizeObservation, quadratizeFinalCost, quadratizeCost, bar.x, bar.sigma, bar.u, bar.L, iterations);
			solvePOMDP<X_DIM, U_DIM, Z_DIM>(experiment, step, bar.x, bar.sigma, bar.u, bar.L, iterations);

			// Simulate one step to get new belief
			sys.executeControlStep(bar.x[0], bar.u[0], bar.L[0], xReal[t], xReal[t+1], hat.x[t], hat.sigma[t], hat.x[t+1], hat.sigma[t+1], hat.u[t]);

			iterations = subIters;
		}
	}

	int initIters;
	int subIters;
};

