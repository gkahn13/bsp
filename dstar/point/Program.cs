using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Expressions;
using System.Text;
using System.Diagnostics;
using DifferentiableFunction;
using matrix;
using utilities;
//using lapack;
using System.IO;

// Figure out how to pass in matrix arguments to a function and evaluate a runtime function -- DONE
// -- order variables to pass in to function? -- SORT OF
// test ekf implementation -- DONE
// belief dynamics function compile -- DONE
// cost function based on belief eval function compile -- DONE
// cleanup code -- DONE
// derivative of the cost function w.r.t state and control inputs -- Jacobian and Hessian -- DONE

namespace Example_CreatingRuntimeFunction
{
    public class Program
    {
        double DT;
        int T;
        int XDIM, UDIM, ZDIM, QDIM, RDIM;
        Function[,] Q, R;
        Function[][] X, U;
        Function[] q, r;
        Function[,] Sigma_0;

        double[][] Xvals, Uvals;
        double[,] Sigma_0val;

        const uint COMPUTE_COST = (1 << 0);
        const uint COMPUTE_JACOBIAN = (1 << 1);
        const uint COMPUTE_HESSIAN = (1 << 2);
        const uint COMPUTE_DIAGONAL_HESSIAN = (1 << 3);
        
        public Program()
        {
            DT = 1;
            XDIM = 2; UDIM = 2; ZDIM = 2; QDIM = 2; RDIM = 2;
            Q = VM.identity<Function>(QDIM, 0, 1);
            R = VM.identity<Function>(RDIM, 0, 1);
        }

        public static Variable[] initializeInputVariables(Function func, Variable[] testInput, out bool[] ivi)
        {
            int i, j = 0;
            Variable[] tDom = func.trueDomain();
            ivi = new bool[testInput.GetLength(0)];
            Variable[] inputVar = new Variable[tDom.GetLength(0)];
            for (i = 0; i < testInput.GetLength(0); i++)
            {
                if (tDom.Contains<Variable>(testInput[i]))
                {
                    ivi[i] = true;
                    inputVar[j++] = testInput[i];
                }
                else
                    ivi[i] = false;
            }
            return inputVar;
        }

        Function[,] mdivide(Function[,] A, Function[,] B)
        {
            // Cholesky factorization A = L*~L
            // check if symmetric
            bool check = VM.checkSymmetric<Function>(A);
            int size = A.GetLength(0);
            Function[,] L = VM.constant<Function>(size, size, 0);
            for (int i = 0; i < size; ++i)
            {
                for (int j = i; j < size; ++j)
                {
                    Function sum = A[j, i];
                    for (int k = 0; k < i; ++k)
                    {
                        sum -= L[j, k] * L[i, k];
                    }
                    if (i == j)
                    {
                        L[i, i] = Function.sqrt(sum);
                    }
                    else
                    {
                        L[j, i] = sum / L[i, i];
                    }
                }
            }

            int ncols = B.GetLength(1);
            // Backward and forward substitution
            Function[,] M = new Function[size, ncols];
            for (int i = 0; i < size; ++i)
            {
                for (int k = 0; k < ncols; ++k)
                {
                    Function sum = B[i, k];
                    for (int j = 0; j < i; ++j)
                    {
                        sum -= L[i, j] * M[j, k];
                    }
                    M[i, k] = sum / L[i, i];
                }
            }
            for (int i = size - 1; i != -1; --i)
            {
                for (int k = 0; k < ncols; ++k)
                {
                    Function sum = M[i, k];
                    for (int j = i + 1; j < size; ++j)
                    {
                        sum -= L[j, i] * M[j, k];
                    }
                    M[i, k] = sum / L[i, i];
                }
            }
            return M;
        }

        double[] dynfunc(double[] x_t, double[] u_t)
        {
            double[] x_tp1 = new double[2];
            x_tp1[0] = x_t[0] + u_t[0] * DT;
            x_tp1[1] = x_t[1] + u_t[1] * DT;
            return x_tp1;
        }

        Function[] dynfunc(Function[] x_t, Function[] u_t, Function[] q_t)
        {
            Function[] x_tp1 = VM.plus(VM.plus(x_t, VM.mult(u_t, DT)), VM.mult(q_t, 0.01));
            return x_tp1;
        }

        Function[] obsfunc(Function[] x_t, Function[] r_t)
        {
            Function intensity = Function.sqrt(x_t[0] * x_t[0] * 0.5 * 0.5 + 1e-6);
            Function[] z_t = VM.plus(x_t, VM.mult(r_t, intensity));
            return z_t;
        }

        void computeDynJacobians(Function[] x_t, Function[] u_t, Function[] q_t, out Function[,] A, out Function[,] M)
        {
            // A = d(dynfunc)/dx
            A = new Function[XDIM, XDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                Function[] Aicol = Function.D(dynfunc(x_t, u_t, q_t), x_t[i]);
                for (int j = 0; j < XDIM; ++j)
                {
                    A[i, j] = Aicol[j];
                }
            }

            // M = d(dynfunc)/dq
            M = new Function[XDIM, QDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                Function[] Micol = Function.D(dynfunc(x_t, u_t, q_t), q_t[i]);
                for (int j = 0; j < QDIM; ++j)
                {
                    M[i, j] = Micol[j];
                }
            }
        }

        void computeObsJacobians(Function[] x_tp1, Function[] r_t, out Function[,] H, out Function[,] N)
        {
            // H = d(obsfunc)/dx
            H = new Function[ZDIM, XDIM];
            for (int i = 0; i < ZDIM; ++i)
            {
                Function[] Hicol = Function.D(obsfunc(x_tp1, r_t), x_tp1[i]);
                for (int j = 0; j < XDIM; ++j)
                {
                    H[i, j] = Hicol[j];
                }
            }

            // N = d(obsfunc)/dr
            N = new Function[ZDIM, RDIM];
            for (int i = 0; i < ZDIM; ++i)
            {
                Function[] Nicol = Function.D(obsfunc(x_tp1, r_t), r_t[i]);
                for (int j = 0; j < RDIM; ++j)
                {
                    N[i, j] = Nicol[j];
                }
            }

            // testing
            //VM.print<Function>(N);
            //Function[,] Nd = Function.derivative(N);
            //VM.print<Function>(Nd);
        }

        void EKF(Function[] x_t, Function[] u_t, Function[] q_t, Function[] r_t, Function[,] Sigma_t, out Function[] x_tp1, out Function[,] Sigma_tp1)
        {
            VM.checkSize<Function>(Sigma_t, XDIM, XDIM);

            Function[,] A, M;
            computeDynJacobians(x_t, u_t, q_t, out A, out M);

            Sigma_tp1 = VM.plus(VM.mult(VM.mult(A, Sigma_t), VM.transpose<Function>(A)), VM.mult(VM.mult(M, Q), VM.transpose<Function>(M)));
            //VM.print<Function>(Sigma_tp1);

            for (int i = 0; i < QDIM; ++i) { q_t[i] = 0; }

            x_tp1 = dynfunc(x_t, u_t, q_t);
            //VM.print<Function>(x_tp1);

            Function[,] H, N;
            computeObsJacobians(x_tp1, r_t, out H, out N);

            Function[,] K1 = VM.mult(Sigma_tp1, VM.transpose<Function>(H));
            Function[,] K2 = VM.plus(VM.mult(VM.mult(H, Sigma_tp1), VM.transpose<Function>(H)), VM.mult(VM.mult(N, R), VM.transpose<Function>(N)));

            Function[,] K = VM.transpose(mdivide(VM.transpose(K2), VM.transpose(K1)));

            Sigma_tp1 = VM.mult(VM.minus(VM.identity<Function>(XDIM, 0, 1), VM.mult(K, H)), Sigma_tp1);
            //VM.print<Function>(Sigma_tp1);
        }

        void beliefUpdate(Function[][] X, Function[][] U, Function[] q, Function[] r, Function[,] Sigma_0, out Function[] x_T, out Function[,] Sigma_T)
        {
            int T = X.GetLength(0);

            Function[] x_tp1 = X[0];
            Function[,] Sigma_t = Sigma_0, Sigma_tp1;

            for (int t = 0; t < T - 1; ++t)
            {
                EKF(X[t], U[t], q, r, Sigma_t, out x_tp1, out Sigma_tp1);
                Sigma_t = Sigma_tp1;
            }

            x_T = x_tp1;
            Sigma_T = Sigma_t;
        }

        Function costfunc(Function[][] X, Function[][] U, Function[] q, Function[] r, Function[,] Sigma_0, Function alpha_belief, Function alpha_control, Function alpha_final_belief)
        {
            int T = X.GetLength(0);

            Function cost = Constant.newConstant(0.0);
            
            Function[] x_tp1 = X[0];
            Function[,] Sigma_t = Sigma_0, Sigma_tp1;

            for (int t = 0; t < T - 1; ++t)
            {
                cost = cost + alpha_belief * VM.trace(Sigma_t);
                cost = cost + alpha_control * VM.dot(U[t], U[t]);

                EKF(X[t], U[t], q, r, Sigma_t, out x_tp1, out Sigma_tp1);
                Sigma_t = Sigma_tp1;
            }

            cost = cost + alpha_final_belief * VM.trace(Sigma_t);
            return cost;
        }

        // one EKF time step check
        void testEKFStep()
        {
            int nvars = 12;
            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }
            Function[] x = new Function[] { vars[0], vars[1] };
            Function[] u = new Function[] { vars[2], vars[3] };

            Function[] q = new Function[] { vars[4], vars[5] };
            Function[] r = new Function[] { vars[6], vars[7] };

            Function[,] Sigma = new Function[,] { { vars[8], vars[9] }, { vars[10], vars[11] } };

            Function[] x_tp1;
            Function[,] Sigma_tp1;

            EKF(x, u, q, r, Sigma, out x_tp1, out Sigma_tp1);

            Function[] packed = VM.concatenate<Function>(x_tp1, VM.matrixToVector<Function>(Sigma_tp1));
            Function beliefvec = Function.derivative(packed);
            beliefvec.printOperatorCounts();

            bool[] inputVarIndices;
            vars = initializeInputVariables(beliefvec, vars, out inputVarIndices);
            beliefvec.orderVariablesInDomain(vars);

            RuntimeFunction brun = beliefvec.compile();

            double[] result = new double[6], varvals = { -4, 2, 0.214285714285714, -0.285714285714286, 0, 0, 1, 0, 0, 1 };
            brun.eval(result, varvals);
            for (int i = 0; i < 6; ++i)
            {
                Console.Write("b[" + i + "]: " + result[i] + " ");
            }
            Console.WriteLine();
        }

        int initVars(int T, Variable[] vars)
        {
            X = new Function[T][];
            U = new Function[T - 1][];

            int idx = 0;
            for (int t = 0; t < T; ++t)
            {
                X[t] = new Function[XDIM];
                for (int i = 0; i < XDIM; ++i)
                {
                    X[t][i] = vars[idx++];
                }
            }
            for (int t = 0; t < (T - 1); ++t)
            {
                U[t] = new Function[UDIM];
                for (int i = 0; i < UDIM; ++i)
                {
                    U[t][i] = vars[idx++];
                }
            }
            q = new Function[] { vars[idx++], vars[idx++] };
            r = new Function[] { vars[idx++], vars[idx++] };

            Sigma_0 = new Function[XDIM, XDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                for (int j = 0; j < XDIM; ++j)
                {
                    Sigma_0[i, j] = vars[idx++];
                }
            }
            return idx;
        }

        int initVarVals(int T, double[] varvals)
        {
            // For evaluation
            double[] u = new double[UDIM];
            double[] start = new double[XDIM];
            start[0] = -3.5; start[1] = 2;
            double[] goal = new double[XDIM];
            goal[0] = -3.5; goal[1] = -2;

            u[0] = (goal[0] - start[0]) / (T - 1);
            u[1] = (goal[1] - start[1]) / (T - 1);

            //Console.WriteLine("u[0]: " + u[0] + " u[1]: " + u[1]);

            // initialize varvals
            int idx = 0;

            Xvals = new double[T][];
            Uvals = new double[T - 1][];

            for (int t = 0; t < (T - 1); ++t)
            {
                Uvals[t] = new double[UDIM];
                Uvals[t][0] = u[0];
                Uvals[t][1] = u[1];
            }

            Xvals[0] = new double[XDIM];
            Xvals[0][0] = start[0]; Xvals[0][1] = start[1];
            for (int t = 1; t < T; ++t)
            {
                Xvals[t] = new double[XDIM];
                Xvals[t] = dynfunc(Xvals[t - 1], Uvals[t - 1]);
            }
            //Console.WriteLine("X[T]: " + Xvals[T - 1][0] + " " + Xvals[T - 1][1]);

            Sigma_0val = VM.identity<double>(XDIM, 0, 1);
            
            idx = 0;

            for (int t = 0; t < T; ++t)
            {
                for (int i = 0; i < XDIM; ++i)
                {
                    varvals[idx++] = Xvals[t][i];
                }
            }
            for (int t = 0; t < (T - 1); ++t)
            {
                for (int i = 0; i < UDIM; ++i)
                {
                    varvals[idx++] = Uvals[t][i];
                }
            }
            for (int i = 0; i < 4; ++i)
            {
                varvals[idx++] = 0;
            }
            varvals[idx++] = Sigma_0val[0, 0]; varvals[idx++] = Sigma_0val[0, 1]; varvals[idx++] = Sigma_0val[1, 0]; varvals[idx++] = Sigma_0val[1, 1];
            return idx;
        }

        void testBeliefUpdate()
        {
            // num timesteps
            int T = 15;

            // variable instantiations
            int nparams = 3;
            int nvars = T * XDIM + (T - 1) * UDIM + QDIM + RDIM + (XDIM * XDIM) + nparams;
            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            int idx = initVars(T, vars);

            Function alpha_belief = vars[idx++];
            Function alpha_control = vars[idx++];
            Function alpha_final_belief = vars[idx++];

            Function[] x_T;
            Function[,] Sigma_T;
            beliefUpdate(X, U, q, r, Sigma_0, out x_T, out Sigma_T);

            Function[] packed = VM.concatenate<Function>(x_T, VM.matrixToVector<Function>(Sigma_T));
            Function beliefvec = Function.derivative(packed);
            beliefvec.printOperatorCounts();

            bool[] inputVarIndices;
            vars = initializeInputVariables(beliefvec, vars, out inputVarIndices);
            beliefvec.orderVariablesInDomain(vars);

            RuntimeFunction beliefrun = beliefvec.compile();

            //beliefrun.compileCCodeToFile("test.c");

            // For evaluation
            double[] varvals = new double[nvars];
            double[] inputvarvals = new double[beliefrun.domainDimension];

            idx = initVarVals(T, varvals);

            int inputidx = 0;
            for (int i = 0; i < nvars; ++i)
            {
                if (inputVarIndices[i])
                {
                    inputvarvals[inputidx++] = varvals[i];
                }
            }

            double[] result = new double[6];
            beliefrun.eval(result, inputvarvals);
            for (int i = 0; i < 6; ++i)
            {
                Console.Write("b[" + i + "]: " + result[i] + " ");
            }
            Console.WriteLine();
        }

        void testCostFunc(uint flag)
        {
            // num timesteps
            T = 15;

            // variable instantiations
            int nparams = 3;
            int nvars = T * XDIM + (T - 1) * UDIM + QDIM + RDIM + (XDIM * XDIM) + nparams;

            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            int idx = initVars(T, vars);

            Function alpha_belief = vars[idx++];
            Function alpha_control = vars[idx++];
            Function alpha_final_belief = vars[idx++];

            Function cost, costJac, costHess, costDiagHess;
            RuntimeFunction costRun, costJacRun, costHessRun, costDiagHessRun;
            Function[] costJacFunc, costDiagHessFunc;
            Function[,] costHessFunc;

            Function[] XU = VM.concatenate<Function>(VM.jaggedToLinear<Function>(X), VM.jaggedToLinear<Function>(U));
            
            // init variable values
            double[] varvals = new double[nvars];

            idx = initVarVals(T, varvals);

            varvals[idx++] = 10;
            varvals[idx++] = 1;
            varvals[idx++] = 10;

            cost = costfunc(X, U, q, r, Sigma_0, alpha_belief, alpha_control, alpha_final_belief);
            
            if ((flag & COMPUTE_COST) != 0)
            {
                cost = Function.derivative(cost);
                cost.printOperatorCounts();

                bool[] inputVarIndices;
                Variable[] costvars = initializeInputVariables(cost, vars, out inputVarIndices);
                cost.orderVariablesInDomain(costvars);

                System.IO.StreamWriter fh = new System.IO.StreamWriter("cost-mask.txt");
                fh.Write(nvars + " ");
                for (int i = 0; i < nvars; ++i)
                {
                    if (inputVarIndices[i])
                    {
                        fh.Write("1 ");
                    }
                    else
                    {
                        fh.Write("0 ");
                    }
                }
                fh.WriteLine();
                fh.Close();        

                cost.compileCCodeToFile("cost.c");

                costRun = cost.compile();

                double[] inputvarvals = new double[costRun.domainDimension];

                int inputidx = 0;
                for (int i = 0; i < nvars; ++i)
                {
                    if (inputVarIndices[i])
                    {
                        inputvarvals[inputidx++] = varvals[i];
                    }
                }
                for (int i = 0; i < costRun.domainDimension; ++i)
                {
                    Console.Write(inputvarvals[i] + " ");
                }
                Console.WriteLine();

                double[] result = new double[1];
                costRun.eval(result, inputvarvals);
                Console.Write("cost: " + result[0]);
                Console.WriteLine();
            }
 
            if ((flag & COMPUTE_JACOBIAN) != 0) 
            {
                costJacFunc = new Function[XU.Length];
                for (int i = 0; i < XU.Length; ++i)
                {
                    costJacFunc[i] = Function.D(cost, XU[i]);
                }

                costJac = Function.derivative(costJacFunc);
                costJac.printOperatorCounts();
                
                bool[] inputVarIndices;
                Variable[] costJacVars = initializeInputVariables(costJac, vars, out inputVarIndices);
                costJac.orderVariablesInDomain(costJacVars);

                //costJac.compileCCodeToFile("costJac.c");

                costJacRun = costJac.compile();

                double[] inputvarvals = new double[costJacRun.domainDimension];

                int inputidx = 0;
                for (int i = 0; i < nvars; ++i)
                {
                    if (inputVarIndices[i])
                    {
                        inputvarvals[inputidx++] = varvals[i];
                    }
                }

                double[] result = new double[XU.Length];
                costJacRun.eval(result, inputvarvals); 
                for (int i = 0; i < XU.Length; ++i)
                {
                    Console.Write("J[" + i + "]: " + result[i] + " ");
                }
                Console.WriteLine();
            }

            if ((flag & COMPUTE_DIAGONAL_HESSIAN) != 0)
            {
                costDiagHessFunc = new Function[XU.Length];
                for (int i = 0; i < XU.Length; ++i)
                {
                    costDiagHessFunc[i] = Function.D(cost, XU[i], XU[i]);
                }

                costDiagHess = Function.derivative(costDiagHessFunc);
                costDiagHess.printOperatorCounts();

                bool[] inputVarIndices;
                Variable[] costDiagHessVars = initializeInputVariables(costDiagHess, vars, out inputVarIndices);
                costDiagHess.orderVariablesInDomain(costDiagHessVars);

                //costDiagHess.compileCCodeToFile("costDiagHess.c");

                costDiagHessRun = costDiagHess.compile();

                double[] inputvarvals = new double[costDiagHessRun.domainDimension];

                int inputidx = 0;
                for (int i = 0; i < nvars; ++i)
                {
                    if (inputVarIndices[i])
                    {
                        inputvarvals[inputidx++] = varvals[i];
                    }
                }

                double[] result = new double[XU.Length];
                costDiagHessRun.eval(result, inputvarvals);
                for (int i = 0; i < XU.Length; ++i)
                {
                    Console.Write("DH[" + i + "]: " + result[i] + " ");
                }
                Console.WriteLine();
            }
            
            if ((flag & COMPUTE_HESSIAN) != 0)
            {
                costHessFunc = new Function[XU.Length, XU.Length];

                // symmetric Hessian
                for (int j = 0; j < XU.Length; ++j)
                {
                    for (int i = j; i < XU.Length; ++i)
                    {
                        //costHessFunc[i, j] = Function.D(costJacFunc[j], XU[i]);
                        costHessFunc[i, j] = Function.D(cost, XU[i], XU[j]);
                        costHessFunc[j, i] = costHessFunc[i, j];
                    }
                }

                costHess = Function.derivative(VM.matrixToVector<Function>(costHessFunc));
                costHess.printOperatorCounts();

                bool[] inputVarIndices;
                Variable[] costHessVars = initializeInputVariables(costHess, vars, out inputVarIndices);
                costHess.orderVariablesInDomain(costHessVars);

                //costHess.compileCCodeToFile("costHess.c");

                //costHessRun = costHess.compile();

                /*
                double[] inputvarvals = new double[costHessRun.domainDimension];

                int inputidx = 0;
                for (int i = 0; i < nvars; ++i)
                {
                    if (inputVarIndices[i])
                    {
                        inputvarvals[inputidx++] = varvals[i];
                    }
                }

                double[] result = new double[XU.Length*XU.Length];
                costJacRun.eval(result, inputvarvals);
                for (int i = 0; i < XU.Length; ++i)
                {
                    Console.Write("J[" + i + "]: " + result[i] + " ");
                }
                Console.WriteLine();
                */
            }
        }

        void computeCostGradDiagHess(int T)
        {
            // num timesteps
            //T = 10;

            // variable instantiations
            int nparams = 3;
            int nvars = T * XDIM + (T - 1) * UDIM + QDIM + RDIM + (XDIM * XDIM) + nparams;

            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            int idx = initVars(T, vars);

            Function alpha_belief = vars[idx++];
            Function alpha_control = vars[idx++];
            Function alpha_final_belief = vars[idx++];
 
            Function[] XU = VM.concatenate<Function>(VM.jaggedToLinear<Function>(X), VM.jaggedToLinear<Function>(U));

            Function cost = costfunc(X, U, q, r, Sigma_0, alpha_belief, alpha_control, alpha_final_belief);

            Function[] costJacDiagHessFunc = new Function[2*XU.Length+1];
            //Function[] costJacDiagHessFunc = new Function[1];
            costJacDiagHessFunc[0] = cost;
            idx = 1;
            for (int i = 0; i < XU.Length; ++i)
            {
                costJacDiagHessFunc[idx++] = Function.D(cost, XU[i]);
            }
            for (int i = 0; i < XU.Length; ++i)
            {
                costJacDiagHessFunc[idx++] = Function.D(cost, XU[i], XU[i]);
            }
            
            Function costJacDiagHess = Function.derivative(costJacDiagHessFunc);
            costJacDiagHess.printOperatorCounts();

            bool[] inputVarIndices;
            Variable[] costJacDiagHessVars = initializeInputVariables(costJacDiagHess, vars, out inputVarIndices);
            costJacDiagHess.orderVariablesInDomain(costJacDiagHessVars);

            costJacDiagHess.compileCCodeToFile("costStateJacDiagHess"+T+".c");

            string fileName = "state-masks-"+T+".txt";
            Console.WriteLine("Writing " + fileName);
            System.IO.StreamWriter fh = new System.IO.StreamWriter(fileName);
            fh.Write(nvars + " ");
            for (int i = 0; i < nvars; ++i)
            {
                if (inputVarIndices[i])
                {
                    fh.Write("1 ");
                }
                else
                {
                    fh.Write("0 ");
                }
            }
            fh.WriteLine();
            fh.Close();        
        }

        Function controlCostFunc(Function[][] U, Function[] q, Function[] r, Function[] x_0, Function[] x_goal, Function[,] Sigma_0, Function alpha_belief, Function alpha_control, Function alpha_final_belief, Function alpha_goal_state)
        {
            int T = U.GetLength(0) + 1;

            Function cost = Constant.newConstant(0.0);

            Function[] x_t = x_0, x_tp1;
            Function[,] Sigma_t = Sigma_0, Sigma_tp1;

            for (int t = 0; t < T - 1; ++t)
            {
                cost = cost + alpha_belief * VM.trace(Sigma_t);
                cost = cost + alpha_control * VM.dot(U[t], U[t]);

                EKF(x_t, U[t], q, r, Sigma_t, out x_tp1, out Sigma_tp1);
                x_t = x_tp1;
                Sigma_t = Sigma_tp1;
            }

            cost = cost + alpha_final_belief * VM.trace(Sigma_t) + alpha_goal_state*VM.dot(VM.minus(x_t,x_goal), VM.minus(x_t,x_goal));
            return cost;
        }

        void computeControlCostGradDiagHess()
        {
            // num timesteps
            T = 40;

            // variable instantiations
            int nparams = 4;
            int nvars = (T - 1) * UDIM + QDIM + RDIM + XDIM + XDIM + (XDIM * XDIM) + nparams;

            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            U = new Function[T - 1][];

            int idx = 0;
            for (int t = 0; t < (T - 1); ++t)
            {
                U[t] = new Function[UDIM];
                for (int i = 0; i < UDIM; ++i)
                {
                    U[t][i] = vars[idx++];
                }
            }
            q = new Function[] { vars[idx++], vars[idx++] };
            r = new Function[] { vars[idx++], vars[idx++] };

            Function[] x_0 = new Function[XDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                x_0[i] = vars[idx++];
            }
            Function[] x_goal = new Function[XDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                x_goal[i] = vars[idx++];
            }
            Sigma_0 = new Function[XDIM, XDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                for (int j = 0; j < XDIM; ++j)
                {
                    Sigma_0[i, j] = vars[idx++];
                }
            }
            
            Function alpha_belief = vars[idx++];
            Function alpha_control = vars[idx++];
            Function alpha_final_belief = vars[idx++];
            Function alpha_goal_state = vars[idx++];

            Function[] Uvec = VM.jaggedToLinear<Function>(U);

            Function cost = controlCostFunc(U, q, r, x_0, x_goal, Sigma_0, alpha_belief, alpha_control, alpha_final_belief, alpha_goal_state);

            Function[] costJacDiagHessFunc = new Function[2*Uvec.Length + 1];
            //Function[] costJacDiagHessFunc = new Function[1];
            costJacDiagHessFunc[0] = cost;
            idx = 1;
            for (int i = 0; i < Uvec.Length; ++i)
            {
                costJacDiagHessFunc[idx++] = Function.D(cost, Uvec[i]);
            }
            for (int i = 0; i < Uvec.Length; ++i)
            {
                costJacDiagHessFunc[idx++] = Function.D(cost, Uvec[i], Uvec[i]);
            }
            Function costJacDiagHess = Function.derivative(costJacDiagHessFunc);
            costJacDiagHess.printOperatorCounts();

            bool[] inputVarIndices;
            Variable[] costJacDiagHessVars = initializeInputVariables(costJacDiagHess, vars, out inputVarIndices);
            costJacDiagHess.orderVariablesInDomain(costJacDiagHessVars);

            costJacDiagHess.compileCCodeToFile("costControlJacDiagHess40.c");
            //costJacDiagHess.compileCCodeToFile("costControl25.c");

            System.IO.StreamWriter fh = new System.IO.StreamWriter("control-masks-40.txt");
            fh.Write(nvars + " ");
            for (int i = 0; i < nvars; ++i)
            {
                Console.WriteLine(inputVarIndices[i]);
                if (inputVarIndices[i])
                {
                    fh.Write("1 ");
                }
                else
                {
                    fh.Write("0 ");
                }
            }
            fh.WriteLine();
            fh.Close();
        }

        static void Main(string[] args)
        { 
            // T = args[0]
            int T;
           
            if (args.Length == 0) {
                T = 15;
            } else {
                T = int.Parse(args[0]);
            }
           
            Console.WriteLine("Creating files for T = "+T);
           
            Function.newContext();
            Function.printCompilerSource = false;

            Program prog = new Program();
            prog.T = T;

            var stopwatch = new Stopwatch();
            stopwatch.Start();
            //prog.testEKFStep();

            //prog.testBeliefUpdate();

            //prog.testCostFunc(Program.COMPUTE_COST);

            prog.computeCostGradDiagHess(T);

            //prog.computeControlCostGradDiagHess();

            stopwatch.Stop();
            
            Console.WriteLine("Finished in " + (stopwatch.ElapsedMilliseconds/1000.0) + " s, Ctrl-C to exit");
            //Console.Read();
        }
    }
}