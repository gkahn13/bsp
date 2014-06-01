class System {
public:
	System(const Dynamics& e, std::string _name = "simulation") : name(_name), step(DT), experiment(e) {}

	virtual void executeControlStep(const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, const Matrix<U_DIM, X_DIM>& L,
				const Matrix<SIM_X_DIM>& xReal, Matrix<SIM_X_DIM>& xRealNext,
				const Matrix<X_DIM>& xHat, const SymmetricMatrix<X_DIM>& SigmaHat,
				Matrix<X_DIM>& xHatNext, SymmetricMatrix<X_DIM>& SigmaHatNext,
				Matrix<U_DIM>& uHat) = 0;
	virtual void executeControls(const bar_t& bar, hat_t& hat) = 0;

protected:
	std::string name;

	const double step;
	const Dynamics& experiment;
};

class Simulation : public System {
public:
	Simulation(const Dynamics& e) : System(e, "simulation") {}

	virtual void executeControlStep(const Matrix<X_DIM>& xBar, const Matrix<U_DIM>& uBar, const Matrix<U_DIM, X_DIM>& L,
					const Matrix<SIM_X_DIM>& xReal, Matrix<SIM_X_DIM>& xRealNext,
					const Matrix<X_DIM>& xHat, const SymmetricMatrix<X_DIM>& SigmaHat,
					Matrix<X_DIM>& xHatNext, SymmetricMatrix<X_DIM>& SigmaHatNext,
					Matrix<U_DIM>& uHat) {
		Matrix<X_DIM, X_DIM> A;
		Matrix<X_DIM, U_DIM> B;
		SymmetricMatrix<X_DIM> M;
		Matrix<Z_DIM, X_DIM> H;
		SymmetricMatrix<Z_DIM> N;
		Matrix<X_DIM> xHatPrime;
		SymmetricMatrix<X_DIM> Lambdat;
		Matrix<X_DIM> linPt;
		SymmetricMatrix<X_DIM> W;

		// Calculate control
		uHat = L*(xHat - xBar) + uBar;
		experiment.saturateControl(uHat);

		// Update real state
#if defined(MAX_LIKELIHOOD)
		xRealNext = experiment.freal(step, xReal, uHat);
#else
		xRealNext = sampleGaussian(experiment.freal(step, xReal, uHat), experiment.varMReal(xReal, uHat));
#endif

		// Sense real state
#if defined(MAX_LIKELIHOOD)
		Matrix<Z_DIM> z = experiment.h(xRealNext);
#else
		Matrix<Z_DIM> z = sampleGaussian(experiment.h(xRealNext), experiment.varN(xRealNext,uHat));
#endif

	#if ITERATE_KALMAN > 0
		int iterations = ITERATE_KALMAN;
		linPt = xHat;

		for (int i = 0; i < iterations; ++i) {
			// Advance state using the dynamics
			experiment.linearizeDynamics(step, linPt, uHat, xHatPrime, A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);

			xHatPrime += A*(xHat - linPt);

			// Calculate Lambda/SigmaHatPrime
			SigmaHatNext = Lambdat = SymProd(A,SigmaHat*(~A)) + M;

			// Calculate observation matrix
			experiment.linearizeObservation(step, xHatPrime, uHat, H, N);

			// Calculate the Kalman gain
			Matrix<X_DIM, Z_DIM> HT = ~H;
			Matrix<Z_DIM, X_DIM> HSigma = H*Lambdat;
			//Matrix<X_DIM, Z_DIM> Kt = ~HSigma*(SymProd(HSigma,~H) + N);
			W = SymProd(~HSigma,(SymProd(HSigma,~H) + N) % HSigma);
			Matrix<X_DIM, Z_DIM> Kt = (~HSigma)*(!(HSigma*HT + N));

			// Correct the new state using the Kalman gain and the observation
			xHatNext = xHatPrime + Kt*(z - experiment.h(xHatPrime));

			// Calculate new covariance
			SigmaHatNext -= W;

			linPt = xHat - (SigmaHat*((~A)/Lambdat)*(xHatNext - xHatPrime));
		}
	#else

		// Advance state using the dynamics
		experiment.linearizeDynamics(xHat, uHat, xHatNext, A, B, M, COMPUTE_c|COMPUTE_A|COMPUTE_M);

		// Calculate observation matrix
		experiment.linearizeObservation(xHatNext, uHat, H, N);

		// Calculate Lambda and the Kalman gain
		SigmaHatNext = SymProd(A,SigmaHat*(~A)) + M;
		Matrix<X_DIM, Z_DIM> HT = ~H;
		Matrix<Z_DIM, X_DIM> HSigma = H*SigmaHatNext;
		//Matrix<X_DIM, Z_DIM> Kt = ~HSigma*(SymProd(HSigma,~H) + N);
		W = SymProd(~HSigma,(SymProd(HSigma,~H) + N) % HSigma);
		Matrix<X_DIM, Z_DIM> Kt = (~HSigma)*(!(HSigma*HT + N));

		// Calculate new covariance
		SigmaHatNext -= W;

		// Correct the new state using the Kalman gain and the observation
		xHatNext = xHatNext + Kt*(z - h(xHatNext));

		/*
		// Calculate Lambda and the Kalman gain
		Matrix<X_DIM, X_DIM> Lambdat = A*SigmaHat*(~A) + M;
		Matrix<X_DIM, Z_DIM> HT = ~H;
		Matrix<X_DIM, Z_DIM> Kt = Lambdat*HT*(!(H*Lambdat*HT + N));

		// Sense
		Matrix<Z_DIM> z = sampleGaussian(h(xRealNext), varN(xRealNext,uHat));

		// Correct the new state using the Kalman gain and the observation
		xHatNext = xHatNext + Kt*(z - h(xHatNext));

		// Calculate new covariance
		SigmaHatNext = SymSum(0.5*(Lambdat - Kt*H*Lambdat));
		*/
	#endif
	}

	virtual void executeControls(const bar_t& bar, hat_t& hat) {
		// Setup space for data
		int numSteps = bar.u.size();
		std::vector< Matrix<SIM_X_DIM> > xReal(numSteps+1);

		// Initialize belief
		xReal[0] = experiment.x0.subMatrix<SIM_X_DIM>(0,0);
		hat.x[0] = experiment.x0;
		hat.sigma[0] = experiment.Sigma0;

		// Iterate over the designated number of controls
		for (int t = 0; t < numSteps; ++t) {
			executeControlStep(bar.x[t], bar.u[t], bar.L[t], xReal[t], xReal[t+1], hat.x[t], hat.sigma[t], hat.x[t+1], hat.sigma[t+1], hat.u[t]);
		}
	}
};
