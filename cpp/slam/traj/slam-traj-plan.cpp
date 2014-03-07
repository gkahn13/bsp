#include "../slam.h"
#include "slam-traj.h"

#include "util/Timer.h"

int main(int argc, char* argv[])
{
	std::vector<std::vector<Matrix<P_DIM> > > l_list = landmarks_list();

	std::ofstream f;
	logDataHandle("slam/data/slam-traj", f);

	for(int i=0; i < l_list.size(); ++i) {
		initProblemParams(l_list[i]);

		std::vector<Matrix<B_DIM> > B_total(T*NUM_WAYPOINTS);
		std::vector<Matrix<U_DIM> > U_total((T-1)*NUM_WAYPOINTS);
		int B_total_idx = 0, U_total_idx = 0;

		util::Timer solveTimer, trajTimer;
		double totalSolveTime = 0, trajTime = 0;

		double totalTrajCost = 0;

		//x0[2] = nearestAngleFromTo(0, x0[2]); // need to remod back to near 0
		for(int i=0; i < NUM_WAYPOINTS; ++i) {
			LOG_INFO("Going to waypoint %d",i);
			// goal is waypoint position + direct angle + landmarks
			xGoal.insert(0, 0, waypoints[i]);

			// want to be facing the next waypoint
			if (i < NUM_WAYPOINTS - 1) {
				xGoal[2] = atan2(waypoints[i+1][1] - waypoints[i][1], waypoints[i+1][0] - waypoints[i][0]);
			} else {
				xGoal[2] = atan2(xGoal[1] - x0[1], xGoal[0] - x0[0]);
			}


			xGoal.insert(C_DIM, 0, x0.subMatrix<L_DIM,1>(C_DIM,0));

			util::Timer_tic(&trajTimer);
			util::Timer_tic(&solveTimer);

			std::vector<Matrix<U_DIM> > U(T-1);
			bool initTrajSuccess = initTraj(x0.subMatrix<C_DIM,1>(0,0), xGoal.subMatrix<C_DIM,1>(0,0), U, T);
			if (!initTrajSuccess) {
				LOG_ERROR("Failed to initialize trajectory, exiting slam-state");
				exit(-1);
			}

			double initTrajTime = util::Timer_toc(&trajTimer);
			double solveTime = util::Timer_toc(&solveTimer);

			trajTime += initTrajTime;
			totalSolveTime += solveTime;

			vec(x0, SqrtSigma0, B_total[B_total_idx++]);
			for(int t=0; t < T-1; ++t) {
				B_total[B_total_idx] = beliefDynamics(B_total[B_total_idx-1], U[t]);
				B_total_idx++;
				U_total[U_total_idx++] = U[t];
			}

			unVec(B_total[B_total_idx-1], x0, SqrtSigma0);

		}

		LOG_INFO("Total trajectory cost: %4.10f", totalTrajCost);
		LOG_INFO("Total trajectory solve time: %5.3f ms", trajTime*1000);
		LOG_INFO("Total solve time: %5.3f ms", totalSolveTime*1000);

		logDataToFile(f, B_total, totalSolveTime*1000, trajTime*1000);

		//pythonDisplayTrajectory(B_total, U_total, waypoints, landmarks, T*NUM_WAYPOINTS, true);

	}

	return 0;
}
