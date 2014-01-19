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

namespace Point_Dstar
{
    public class Point
    {
        protected double DT;
        protected int T;
        protected int XDIM, UDIM, ZDIM, QDIM, RDIM;
        protected Function[,] Q, R;
        protected Function[][] X, U;
        protected Function[] q, r;
        protected Function[,] Sigma_0;

        public Point() : this(15) {
        }
        
        public Point(int timesteps)
        {
            T = timesteps;
            DT = 1;
            XDIM = 2;
            UDIM = 2;
            ZDIM = 2;
            QDIM = 2;
            RDIM = 2;
            Q = VM.identity<Function>(QDIM, 0, 1);
            R = VM.identity<Function>(RDIM, 0, 1);
        }
        
        protected static Variable[] initializeInputVariables(Function func, Variable[] testInput, out bool[] ivi)
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
        }

        protected void EKF(Function[] x_t, Function[] u_t, Function[] q_t, Function[] r_t, Function[,] Sigma_t, out Function[] x_tp1, out Function[,] Sigma_tp1)
        {
            VM.checkSize<Function>(Sigma_t, XDIM, XDIM);

            Function[,] A, M;
            computeDynJacobians(x_t, u_t, q_t, out A, out M);

            Sigma_tp1 = VM.plus(VM.mult(VM.mult(A, Sigma_t), VM.transpose<Function>(A)), VM.mult(VM.mult(M, Q), VM.transpose<Function>(M)));

            for (int i = 0; i < QDIM; ++i) { q_t[i] = 0; }

            x_tp1 = dynfunc(x_t, u_t, q_t);
            
            Function[,] H, N;
            computeObsJacobians(x_tp1, r_t, out H, out N);

            Function[,] K1 = VM.mult(Sigma_tp1, VM.transpose<Function>(H));
            Function[,] K2 = VM.plus(VM.mult(VM.mult(H, Sigma_tp1), VM.transpose<Function>(H)), VM.mult(VM.mult(N, R), VM.transpose<Function>(N)));

            Function[,] K = VM.transpose(mdivide(VM.transpose(K2), VM.transpose(K1)));

            Sigma_tp1 = VM.mult(VM.minus(VM.identity<Function>(XDIM, 0, 1), VM.mult(K, H)), Sigma_tp1);
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


    }
}