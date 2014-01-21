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

using System.Reflection;

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
        int DIM, XDIM, UDIM, ZDIM, QDIM, RDIM;
        Function[,] Q, R;
        Function[][] X, U;
        Function[] q, r;
        Function[,] Sigma_0;
        double l1, l2, l3, l4;
        double[] cam0, cam1;

        double[][] Xvals, Uvals;
        double[,] Sigma_0val;

        const uint COMPUTE_COST = (1 << 0);
        const uint COMPUTE_JACOBIAN = (1 << 1);
        const uint COMPUTE_HESSIAN = (1 << 2);
        const uint COMPUTE_DIAGONAL_HESSIAN = (1 << 3);
        
        public Program()
        {
            DT = 0.5;
            DIM = 3; XDIM = 6; UDIM = 6; ZDIM = 4; QDIM = 6; RDIM = 4;
            Q = VM.identity<Function>(QDIM, 0, 1);
            R = VM.identity<Function>(RDIM, 0, 1);
            l4 = 2.375; l3 = 10.375; l2 = 8; l1 = 7.25;
            cam0 = new double[DIM]; cam1 = new double[DIM];
            cam0[0] = -4;  cam0[1] = 23; cam0[2] = 0;
            cam1[0] = 4;  cam1[1] = 23; cam1[2] = 0;
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
            double[] x_tp1 = VM.plus(x_t, VM.mult(u_t, DT));
            return x_tp1;
        }

        Function[] dynfunc(Function[] x_t, Function[] u_t, Function[] q_t)
        {
            Function[] x_tp1 = VM.plus(x_t, VM.mult(VM.plus(u_t, q_t), DT));
            return x_tp1;
        }
        
        // joint angles -> end effector position
        double[] g(double[] x)
        {
            double a0 = x[0], a1 = x[1], a2 = x[2], a3 = x[3], a4 = x[4], a5 = x[5];
            
            double[] p = new double[3];
            p[0] = Math.Sin(a0)*(Math.Cos(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)+Math.Sin(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2))+Math.Cos(a0)*Math.Sin(a3)*Math.Sin(a4)*l4;
            p[1] = -Math.Sin(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)+Math.Cos(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2)+l1;
            p[2] = Math.Cos(a0)*(Math.Cos(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)+Math.Sin(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2))-Math.Sin(a0)*Math.Sin(a3)*Math.Sin(a4)*l4;
            
            return p;   
        }
        
        // joint angles -> end effector position
        Function[] g(Function[] x)
        {
            Function a0 = x[0], a1 = x[1], a2 = x[2], a3 = x[3], a4 = x[4], a5 = x[5];
            
            Function[] p = new Function[3];
            p[0] = Function.sin(a0)*(Function.cos(a1)*(Function.sin(a2)*(Function.cos(a4)*l4+l3)+Function.cos(a2)*Function.cos(a3)*Function.sin(a4)*l4)+Function.sin(a1)*(Function.cos(a2)*(Function.cos(a4)*l4+l3)-Function.sin(a2)*Function.cos(a3)*Function.sin(a4)*l4+l2))+Function.cos(a0)*Function.sin(a3)*Function.sin(a4)*l4;
            p[1] = -Function.sin(a1)*(Function.sin(a2)*(Function.cos(a4)*l4+l3)+Function.cos(a2)*Function.cos(a3)*Function.sin(a4)*l4)+Function.cos(a1)*(Function.cos(a2)*(Function.cos(a4)*l4+l3)-Function.sin(a2)*Function.cos(a3)*Function.sin(a4)*l4+l2)+l1;
            p[2] = Function.cos(a0)*(Function.cos(a1)*(Function.sin(a2)*(Function.cos(a4)*l4+l3)+Function.cos(a2)*Function.cos(a3)*Function.sin(a4)*l4)+Function.sin(a1)*(Function.cos(a2)*(Function.cos(a4)*l4+l3)-Function.sin(a2)*Function.cos(a3)*Function.sin(a4)*l4+l2))-Function.sin(a0)*Function.sin(a3)*Function.sin(a4)*l4;
            
            return p;   
        }
        
        
        void analyticalJacobian(double[] x)
        {
            double[,] J = new double[DIM,XDIM];
            double a0 = x[0];
            double a1 = x[1];
            double a2 = x[2];
            double a3 = x[3];
            double a4 = x[4];
            double a5 = x[5];
            
            // J: Jacobian of end-effector position w.r.t joint angles (3x6 matrix)
            J[0,0] = Math.Cos(a0)*(Math.Cos(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)+Math.Sin(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2))-Math.Sin(a0)*Math.Sin(a3)*Math.Sin(a4)*l4;
            J[0,1] = Math.Sin(a0)*(Math.Cos(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2)-Math.Sin(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4));
            J[0,2] = Math.Sin(a0)*(Math.Sin(a1)*(-Math.Sin(a2)*(Math.Cos(a4)*l4+l3)-Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)+Math.Cos(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4));
            J[0,3] = Math.Sin(a0)*(Math.Sin(a1)*Math.Sin(a2)*Math.Sin(a3)*Math.Sin(a4)*l4-Math.Cos(a1)*Math.Cos(a2)*Math.Sin(a3)*Math.Sin(a4)*l4)+Math.Cos(a0)*Math.Cos(a3)*Math.Sin(a4)*l4;
            J[0,4] = Math.Sin(a0)*(Math.Cos(a1)*(Math.Cos(a2)*Math.Cos(a3)*Math.Cos(a4)*l4-Math.Sin(a2)*Math.Sin(a4)*l4)+Math.Sin(a1)*(-Math.Cos(a2)*Math.Sin(a4)*l4-Math.Sin(a2)*Math.Cos(a3)*Math.Cos(a4)*l4))+Math.Cos(a0)*Math.Sin(a3)*Math.Cos(a4)*l4;
            J[0,5] = 0;
            
            J[1,0] = 0;
            J[1,1] = -Math.Cos(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)-Math.Sin(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2);
            J[1,2] = Math.Cos(a1)*(-Math.Sin(a2)*(Math.Cos(a4)*l4+l3)-Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)-Math.Sin(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4);
            J[1,3] = Math.Cos(a1)*Math.Sin(a2)*Math.Sin(a3)*Math.Sin(a4)*l4+Math.Sin(a1)*Math.Cos(a2)*Math.Sin(a3)*Math.Sin(a4)*l4;
            J[1,4] = Math.Cos(a1)*(-Math.Cos(a2)*Math.Sin(a4)*l4-Math.Sin(a2)*Math.Cos(a3)*Math.Cos(a4)*l4)-Math.Sin(a1)*(Math.Cos(a2)*Math.Cos(a3)*Math.Cos(a4)*l4-Math.Sin(a2)*Math.Sin(a4)*l4);
            J[1,5] = 0;
            
            J[2,0] = -Math.Sin(a0)*(Math.Cos(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)+Math.Sin(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2))-Math.Cos(a0)*Math.Sin(a3)*Math.Sin(a4)*l4;
            J[2,1] = Math.Cos(a0)*(Math.Cos(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4+l2)-Math.Sin(a1)*(Math.Sin(a2)*(Math.Cos(a4)*l4+l3)+Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4));
            J[2,2] = Math.Cos(a0)*(Math.Sin(a1)*(-Math.Sin(a2)*(Math.Cos(a4)*l4+l3)-Math.Cos(a2)*Math.Cos(a3)*Math.Sin(a4)*l4)+Math.Cos(a1)*(Math.Cos(a2)*(Math.Cos(a4)*l4+l3)-Math.Sin(a2)*Math.Cos(a3)*Math.Sin(a4)*l4));
            J[2,3] = Math.Cos(a0)*(Math.Sin(a1)*Math.Sin(a2)*Math.Sin(a3)*Math.Sin(a4)*l4-Math.Cos(a1)*Math.Cos(a2)*Math.Sin(a3)*Math.Sin(a4)*l4)-Math.Sin(a0)*Math.Cos(a3)*Math.Sin(a4)*l4;
            J[2,4] = Math.Cos(a0)*(Math.Cos(a1)*(Math.Cos(a2)*Math.Cos(a3)*Math.Cos(a4)*l4-Math.Sin(a2)*Math.Sin(a4)*l4)+Math.Sin(a1)*(-Math.Cos(a2)*Math.Sin(a4)*l4-Math.Sin(a2)*Math.Cos(a3)*Math.Cos(a4)*l4))-Math.Sin(a0)*Math.Sin(a3)*Math.Cos(a4)*l4;
            J[2,5] = 0;
            
            for(int i = 0; i < DIM; ++i) {
                for(int j = 0; j < XDIM; ++j) {
                    //Console.WriteLine(J[i,j]);
                }
            }
            
            double[] pos = g(x);
            
            double[,] H = new double[ZDIM,XDIM];
            // H = dh/dx, jacobian of observation function (ZDIM by XDIM)
		    for (int i = 0; i < XDIM; ++i) {
		        H[0, i] = (J[0, i] * (pos[1] - cam0[1]) - (pos[0] - cam0[0]) * J[1, i]) / ((pos[1] - cam0[1]) * (pos[1] - cam0[1]));
                H[1, i] = (J[2, i] * (pos[1] - cam0[1]) - (pos[2] - cam0[2]) * J[1, i]) / ((pos[1] - cam0[1]) * (pos[1] - cam0[1]));
                H[2, i] = (J[0, i] * (pos[1] - cam1[1]) - (pos[0] - cam1[0]) * J[1, i]) / ((pos[1] - cam1[1]) * (pos[1] - cam1[1]));
                H[3, i] = (J[2, i] * (pos[1] - cam1[1]) - (pos[2] - cam1[2]) * J[1, i]) / ((pos[1] - cam1[1]) * (pos[1] - cam1[1]));        
		    }
		    
            for(int i = 0; i < ZDIM; ++i) {
                for(int j = 0; j < XDIM; ++j) {
                    Console.WriteLine(H[i,j]);
                }
            }
		    
        }
        
        // J = T from the paper
        void symbolicJacobian(Function[] x, out Function[,] H)
        {
            /*
            // J = dg/dx
            J = new Function[DIM,XDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                Function[] Jicol = Function.D(g(x), x[i]);
                for (int j = 0; j < DIM; ++j)
                {
                    J[j,i] = Jicol[j];
                }
            }
            */
            
            // H = d(obsfunc)/dx
            H = new Function[ZDIM, XDIM];
            Function[] r = new Function[RDIM];
            for (int i = 0; i < RDIM; ++i) { r[i] = 0; }
            for (int i = 0; i < XDIM; ++i)
            {
                Function[] Hicol = Function.D(obsfunc(x, r), x[i]);
                for (int j = 0; j < ZDIM; ++j)
                {
                    H[j, i] = Hicol[j];
                }
            }
        }
        
        void testKinematics()
        {
            double[] xtest = new double[XDIM];
            xtest[0] = Math.PI/2; xtest[1] = -1.5431281995798991; xtest[2] = -0.047595544887998331;
            xtest[3] = 1.4423058659586809; xtest[4] = 1.5334368368992011; xtest[5] = -1.1431255223182604;
            
            analyticalJacobian(xtest);
            
            int nvars = 6;
            Variable[] vars = new Variable[nvars];
            Function[] x = new Function[XDIM];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); x[i] = vars[i];}
            
            //Function[,] J, H;
            Function[,] H;
            //symbolicJacobian(x, out J, out H);
            symbolicJacobian(x, out H);
            
            Function Hvec = Function.derivative(VM.jaggedToLinear<Function>(VM.toJaggedArray<Function>(H)));
            Hvec.printOperatorCounts();


            bool[] inputVarIndices;
            vars = initializeInputVariables(Hvec, vars, out inputVarIndices);
            Hvec.orderVariablesInDomain(vars);
            
            RuntimeFunction Hrun = Hvec.compile();
            Console.WriteLine(Hrun.domainDimension + " " + Hrun.rangeDimension);

            double[] result = new double[Hrun.rangeDimension];
            Hrun.eval(result, xtest);
            for (int j = 0; j < Hrun.rangeDimension; ++j)
            {
                Console.WriteLine(result[j]);
            }
            Console.WriteLine();

            /*
            Function Jvec = Function.derivative(VM.jaggedToLinear<Function>(VM.toJaggedArray<Function>(J)));
            //Jvec.printOperatorCounts();
            
            
            bool[] inputVarIndices;
            vars = initializeInputVariables(Jvec, vars, out inputVarIndices);
            for (int i = 0; i < nvars; ++i)
                {
                    if (inputVarIndices[i])
                    {
                        Console.WriteLine("1 ");
                    }
                    else
                    {
                        Console.WriteLine("0 ");
                    }
                }
            
              
            Jvec.orderVariablesInDomain(vars);
            
            RuntimeFunction Jrun = Jvec.compile();
            Console.WriteLine(Jrun.domainDimension + " " + Jrun.rangeDimension);

            double[] result = new double[Jrun.rangeDimension];
            Jrun.eval(result, xtest);
            for (int j = 0; j < Jrun.rangeDimension; ++j)
            {
                Console.WriteLine(result[j]);
            }
            Console.WriteLine();
            */
        }
        
        // observation function (also called h)
        /*double[] obsfunc(double[] x_t)
        {
            double[] ee_pos = g(x_t);
            
            double[] obs = new double[ZDIM];
            obs[0] = (ee_pos[0] - cam0[0])/(ee_pos[2] - cam0[2]);
            obs[1] = (ee_pos[1] - cam0[1])/(ee_pos[2] - cam0[2]);
            obs[2] = (ee_pos[0] - cam1[0])/(ee_pos[2] - cam1[2]);
            obs[3] = (ee_pos[1] - cam1[1])/(ee_pos[2] - cam1[2]);
            
            return obs;
        }*/
        
        // observation function (also called h)
        Function[] obsfunc(Function[] x_t, Function[] r_t)
        {
            Function[] ee_pos = g(x_t);
            
            Function[] obs = new Function[ZDIM];
            
            obs[0] = (ee_pos[0] - cam0[0])/(ee_pos[1] - cam0[1]) + r_t[0];
            obs[1] = (ee_pos[2] - cam0[2])/(ee_pos[1] - cam0[1]) + r_t[1];
            obs[2] = (ee_pos[0] - cam1[0])/(ee_pos[1] - cam1[1]) + r_t[2];
            obs[3] = (ee_pos[2] - cam1[2])/(ee_pos[1] - cam1[1]) + r_t[3];
            
            return obs;
        }

        void computeObsJacobians(Function[] x_tp1, Function[] r_t, out Function[,] H, out Function[,] N)
        {
            // H = d(obsfunc)/dx
            H = new Function[ZDIM, XDIM];
            for (int i = 0; i < XDIM; ++i)
            {
                Function[] Hicol = Function.D(obsfunc(x_tp1, r_t), x_tp1[i]);
                for (int j = 0; j < ZDIM; ++j)
                {
                    H[j, i] = Hicol[j];
                }
            }

            // N = d(obsfunc)/dr
            N = new Function[ZDIM, RDIM];
            for (int i = 0; i < RDIM; ++i)
            {
                Function[] Nicol = Function.D(obsfunc(x_tp1, r_t), r_t[i]);
                for (int j = 0; j < ZDIM; ++j)
                {
                    N[j, i] = Nicol[j];
                }
            }

        }

        void EKF(Function[] x_t, Function[] u_t, Function[] q_t, Function[] r_t, Function[,] Sigma_t, out Function[] x_tp1, out Function[,] Sigma_tp1)
        {
            VM.checkSize<Function>(Sigma_t, XDIM, XDIM);

            Function[,] A, M;
            //computeDynJacobians(x_t, u_t, q_t, out A, out M);
            A = VM.identity<Function>(XDIM, 0, 1); // d(dynfunc)/dx
            M = VM.mult(VM.identity<Function>(UDIM, 0, 1), DT); // d(dynfunc)/dm (noise)

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
            int nvars = XDIM + UDIM + QDIM + RDIM + XDIM*XDIM;
            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }
            Function [] x = new Function[XDIM];
            Function [] u = new Function[UDIM];
            Function [] q = new Function[QDIM];
            Function [] r = new Function[RDIM];
            Function [,] Sigma = new Function[XDIM,XDIM];
            
            int idx = 0;
            for(int i=0; i < XDIM; ++i) { x[i] = vars[idx++]; }
            for(int i=0; i < UDIM; ++i) { u[i] = vars[idx++]; }
            for(int i=0; i < QDIM; ++i) { q[i] = vars[idx++]; }
            for(int i=0; i < RDIM; ++i) { r[i] = vars[idx++]; }
            for(int i=0; i < XDIM; ++i)
            { 
               for(int j=0; j < XDIM; ++j)
               {
                   Sigma[i,j] = vars[idx++];
               }
            }
 
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

            double[] result = new double[nvars];
            double[] varvals = new double[nvars];
            
            double[] joints = {Math.PI/2, -1.5431281995798991, -0.047595544887998331,
                               1.4423058659586809, 1.5334368368992011, -1.1431255223182604};
            double[] input = {-.1, .1, -.1, -.1, -.1, .1};
            double[] x_noise = {0, 0, 0, 0, 0, 0};
            double[] z_noise = {0, 0, 0, 0};
            double[] S = new double[XDIM*XDIM];
            for(int i=0; i < XDIM; ++i)
            {
                for(int j=0; j < XDIM; ++j)
                {
                    if (i == j)
                    {
                        S[i+j*XDIM] = 1;   
                    } else
                    {
                        S[i+j*XDIM] = 0;
                    }
                }    
            }
            
            Console.WriteLine("Before copying");

            joints.CopyTo(varvals, 0);
            input.CopyTo(varvals, XDIM);
            x_noise.CopyTo(varvals, XDIM+UDIM);
            z_noise.CopyTo(varvals, XDIM+UDIM+QDIM);
            S.CopyTo(varvals, XDIM+UDIM+QDIM+RDIM);
            
            Console.WriteLine("After copying");
    
            brun.eval(result, varvals);
            printBelief(result);
        }
        
        void printBelief(double[] b)
        {
            int idx = 0;
            string[] names = {"X","U","Q","R","S"};
            int[] dims = {XDIM,UDIM,QDIM,RDIM,XDIM*XDIM};
            
            for(int i=0; i < dims.Length; ++i)
            {
                Console.WriteLine(names[i]);
                for(int j=0; j < dims[i]; ++j)
                {
                    Console.WriteLine(b[idx++]);
                }
                Console.WriteLine();
            }
        }

        int initVars(Variable[] vars)
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
            
            Function [] q = new Function[QDIM];
            for(int i=0; i < QDIM; ++i) { q[i] = vars[idx++]; }
            
            Function [] r = new Function[RDIM];
            for(int i=0; i < RDIM; ++i) { r[i] = vars[idx++]; }
            
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

        // TODO: (part way done)
        int initVarVals(int T, double[] varvals)
        {
            // For evaluation
            double[] start = {Math.PI/2, -1.5431281995798991, -0.047595544887998331,
                               1.4423058659586809, 1.5334368368992011, -1.1431255223182604};
            double[] goal = {11.5, 11.5, 0}; // TODO: find joint angles!
            
            double[] u = new double[UDIM];
            for(int i=0; i < UDIM; ++i) { u[i] = 0; } // zero initialization for now
            
            //Console.WriteLine("u[0]: " + u[0] + " u[1]: " + u[1]);

            // initialize varvals
            int idx = 0;

            Xvals = new double[T][];
            Uvals = new double[T - 1][];

            for (int t = 0; t < (T - 1); ++t)
            {
                Uvals[t] = new double[UDIM];
                for(int i=0; i < UDIM; ++i) { Uvals[t][i] = u[i]; }
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

        // TODO!
        void testBeliefUpdate()
        {
            // num timesteps
            int T = 15;

            // variable instantiations
            int nparams = 3;
            int nvars = T * XDIM + (T - 1) * UDIM + QDIM + RDIM + (XDIM * XDIM) + nparams;
            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            int idx = initVars(vars);

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

        // TODO!
        void testCostFunc(uint flag)
        {
            // num timesteps
            T = 15;

            // variable instantiations
            int nparams = 3;
            int nvars = T * XDIM + (T - 1) * UDIM + QDIM + RDIM + (XDIM * XDIM) + nparams;

            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            int idx = initVars(vars);

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
        
        void computeCostGradDiagHess(string eval_name)
        {
            uint flag = 0;
            flag |= ((uint)(Convert.ToUInt16(eval_name.Contains("Cost"))) << 0);
            flag |= ((uint)(Convert.ToUInt16(eval_name.Contains("Grad"))) << 1);
            flag |= ((uint)(Convert.ToUInt16(eval_name.Contains("DiagHess"))) << 3);
            
            // variable instantiations
            int nparams = 3;
            int nvars = T * XDIM + (T - 1) * UDIM + QDIM + RDIM + (XDIM * XDIM) + nparams;

            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            int idx = initVars(vars);

            Function alpha_belief = vars[idx++];
            Function alpha_control = vars[idx++];
            Function alpha_final_belief = vars[idx++];
 
            Function[] XU = VM.concatenate<Function>(VM.jaggedToLinear<Function>(X), VM.jaggedToLinear<Function>(U));

            Function cost = costfunc(X, U, q, r, Sigma_0, alpha_belief, alpha_control, alpha_final_belief);

            List<Function> costJacDiagHessFuncList = new List<Function>();
            if ((flag & COMPUTE_COST) != 0)
            {
                costJacDiagHessFuncList.Add(cost);
            }
            
            if ((flag & COMPUTE_JACOBIAN) != 0)
            {
	            for (int i = 0; i < XU.Length; ++i)
	            {
	                costJacDiagHessFuncList.Add(Function.D(cost, XU[i]));
	            }
            }
            
            if ((flag & COMPUTE_DIAGONAL_HESSIAN) != 0)
            {
	            for (int i = 0; i < XU.Length; ++i)
	            {
	                costJacDiagHessFuncList.Add(Function.D(cost, XU[i], XU[i]));
	            }
            }
            
            Function costJacDiagHess = Function.derivative(costJacDiagHessFuncList.ToArray());
            costJacDiagHess.printOperatorCounts();

            bool[] inputVarIndices;
            Variable[] costJacDiagHessVars = initializeInputVariables(costJacDiagHess, vars, out inputVarIndices);
            costJacDiagHess.orderVariablesInDomain(costJacDiagHessVars);

            costJacDiagHess.compileCCodeToFile("state-symeval-"+T+".c");

            string fileName = "state-masks-"+T+".txt";
            Console.WriteLine("Writing " + fileName);
            System.IO.StreamWriter fh = new System.IO.StreamWriter(fileName);
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

        static void Main(string[] args)
        { 
            int T = 15;
            string eval_name = "CostGradDiagHess"; // default to comput cost, gradient, and diagonal hessian
           
            if (args.Length >= 1) {
                T = int.Parse(args[0]);
            }
            if (args.Length >= 2) {
                eval_name = args[1];
            }
           
            Console.WriteLine("Creating files for "+eval_name+" for T = "+T);
           
            Function.newContext();
            Function.printCompilerSource = false;

            Program p = new Program();
            p.T = 15;

            var stopwatch = new Stopwatch();
            stopwatch.Start();

            p.computeCostGradDiagHess(eval_name);

           
            stopwatch.Stop();
            
            Console.WriteLine("Finished in " + (stopwatch.ElapsedMilliseconds/1000.0) + " s");
        }
        
        /*
        static void Main(string[] args)
        { 
            Function.newContext();
            Function.printCompilerSource = false;

            Program prog = new Program();
            prog.T = 15;

            var stopwatch = new Stopwatch();
            stopwatch.Start();
            
            //prog.testKinematics();
            
            prog.testEKFStep();

            //prog.testBeliefUpdate();

            //prog.testCostFunc(Program.COMPUTE_COST);

            //prog.computeCostGradDiagHess(T);

            //prog.computeControlCostGradDiagHess();

            stopwatch.Stop();
            
            Console.WriteLine("Finished in " + (stopwatch.ElapsedMilliseconds/1000.0) + " s, Ctrl-C to exit");
            //Console.Read();
        }
        */
    }
}