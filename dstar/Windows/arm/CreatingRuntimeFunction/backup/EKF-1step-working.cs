using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Expressions;
using System.Text;
using System.Diagnostics;
using DifferentiableFunction;
using matrix;
using utilities;
using lapack;
using System.IO;

// Figure out how to pass in matrix arguments to a function and evaluate a runtime function
// -- order variables to pass in to function?
// test ekf implementation
// beliefeval function compile
// cost function based on belief eval function compile
// derivative of the cost function w.r.t state and control inputs

namespace Example_CreatingRuntimeFunction
{
    public class Program
    {
        double DT;
        int XDIM, UDIM, ZDIM, QDIM, RDIM;
        Function[,] Q, R;

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

        Function[,] mdivide(Function[,] A, Function[,] B) {
	        // Cholesky factorization A = L*~L
            // check if symmetric
            bool check = VM.checkSymmetric<Function>(A);
            int size = A.GetLength(0);
	        Function[,] L = VM.constant<Function>(size, size, 0);
	        for (int i = 0; i < size; ++i) {
		        for (int j = i; j < size; ++j) {
			        Function sum = A[j,i];
			        for (int k = 0; k < i; ++k) {
				        sum -= L[j,k]*L[i,k];
			        }
			        if (i == j) {
                        L[i, i] = Function.sqrt(sum);
			        } else {
				        L[j,i] = sum / L[i,i];
			        }
		        }
	        }

            int ncols = B.GetLength(1);
	        // Backward and forward substitution
	        Function[,] M = new Function[size,ncols];
	        for (int i = 0; i < size; ++i) {
		        for (int k = 0; k < ncols; ++k) {
			        Function sum = B[i,k];
			        for (int j = 0; j < i; ++j) {
				        sum -= L[i,j]*M[j,k];
			        }
			        M[i,k] = sum / L[i,i];
		        }
	        }
	        for (int i = size - 1; i != -1; --i) {
		        for (int k = 0; k < ncols; ++k) {
			        Function sum = M[i,k];
			        for (int j = i + 1; j < size; ++j) {
				        sum -= L[j,i]*M[j,k];
			        }
			        M[i,k] = sum / L[i,i];
		        }
	        }
	        return M;
        }

        public Program()
        {
            DT = 1;
            XDIM = 2; UDIM = 2; ZDIM = 2; QDIM = 2; RDIM = 2;
            Q = VM.identity<Function>(QDIM, 0, 1);
            R = VM.identity<Function>(RDIM, 0, 1);
        }

        Function[] dynfunc(Function[] x_t, Function[] u_t, Function[] q_t)
        {
            Function[] x_tp1 = VM.plus(VM.plus(x_t, VM.mult(u_t, DT)), VM.mult(q_t, 0.01));
            return x_tp1;
        }

        Function[] obsfunc(Function[] x_t, Function[] r_t)
        {
            Function intensity = Function.sqrt(x_t[0]*x_t[0]*0.5*0.5 + 1e-6);
            Function[] z_t = VM.plus(x_t, VM.mult(r_t, intensity));
            return z_t;
        }

        void EKF(Function[] x_t, Function[,] Sigma_t, Function[] u_t, Function[] q_t, Function[] r_t, out Function[] x_tp1, out Function[,] Sigma_tp1)
        {
            VM.checkSize<Function>(Sigma_t, XDIM, XDIM);

            Function[,] A = new Function[XDIM, XDIM];
            for (int i = 0; i < XDIM; ++i) {
                Function[] Aicol = Function.D(dynfunc(x_t, u_t, q_t), x_t[i]);
                for (int j = 0; j < XDIM; ++j) {
                    A[i, j] = Aicol[j];
                }
            }
            //Function[,] Ad = Function.derivative(A);
            //VM.print<Function>(Ad);

            Function[,] M = new Function[XDIM, QDIM];
            for (int i = 0; i < XDIM; ++i) {
                Function[] Micol = Function.D(dynfunc(x_t, u_t, q_t), q_t[i]);
                for (int j = 0; j < QDIM; ++j) {
                    M[i, j] = Micol[j];
                }
            }
            //VM.print<Function>(M);
            //Function[,] Md = Function.derivative(M);
            //VM.print<Function>(Md);

            Sigma_tp1 = VM.plus(VM.mult(VM.mult(A, Sigma_t), VM.transpose<Function>(A)), VM.mult(VM.mult(M, Q), VM.transpose<Function>(M)));
            //VM.print<Function>(Sigma_tp1);

            for (int i = 0; i < QDIM; ++i) { q_t[i] = 0; }
            x_tp1 = dynfunc(x_t, u_t, q_t);
            //VM.print<Function>(x_tp1);
            //Console.WriteLine();

            Function[,] H = new Function[ZDIM, XDIM];
            for (int i = 0; i < ZDIM; ++i) {
                Function[] Hicol = Function.D(obsfunc(x_tp1, r_t), x_tp1[i]);
                for (int j = 0; j < XDIM; ++j) {
                    H[i, j] = Hicol[j];
                }
            }
            //VM.print<Function>(H);
            //Function[,] Hd = Function.derivative(H);
            //VM.print<Function>(Hd);

            Function[,] N = new Function[ZDIM, RDIM];
            for (int i = 0; i < ZDIM; ++i) {
                Function[] Nicol = Function.D(obsfunc(x_tp1, r_t), r_t[i]);
                for (int j = 0; j < RDIM; ++j) {
                    N[i, j] = Nicol[j];
                }
            }
            //VM.print<Function>(N);
            //Function[,] Nd = Function.derivative(N);
            //VM.print<Function>(Nd);

            Function[,] K1 = VM.mult(Sigma_tp1, VM.transpose<Function>(H));
            Function[,] K2 = VM.plus(VM.mult(VM.mult(H, Sigma_tp1), VM.transpose<Function>(H)), VM.mult(VM.mult(N, R), VM.transpose<Function>(N)));
            
            Function[,] K = VM.transpose(mdivide(VM.transpose(K2),VM.transpose(K1)));

            Function[,] I = VM.identity<Function>(XDIM, 0, 1);
            Sigma_tp1 = VM.mult(VM.minus(I, VM.mult(K, H)), Sigma_tp1);
            //VM.print<Function>(Sigma_tp1);
        }

        static void Main(string[] args)
        {
            Function.newContext();
            Function.printCompilerSource = true;
            //Function.turnOnAllSimplification();

            Program prog = new Program();

            //int T = 10;
            //Function[][] X = new Function[T][];
            //Function[][] U = new Function[T-1][];

            //Function[] x = new Function[] {-4, 2};
            //Function[] u = new Function[] {0.214285714285714, -0.285714285714286};
            //Function[,] Sigma = new Function[,] { { 1, 0 }, { 0, 1 } };

            Variable[] vars = new Variable[8];
            for (int i = 0; i < 8; ++i) { vars[i] = new Variable("vars_" + i); }
            Function[] x = new Function[] {vars[0], vars[1]};
            Function[] u = new Function[] {vars[2], vars[3]};
            Function[,] Sigma = VM.identity<Function>(2, 0, 1);

            Function[] q_t = new Function[] {vars[4], vars[5]};
            Function[] r_t = new Function[] {vars[6], vars[7]};

            Function[] x_tp1;
            Function[,] Sigma_tp1;
            
            prog.EKF(x, Sigma, u, q_t, r_t, out x_tp1, out Sigma_tp1);

            Function[] packed = VM.concatenate<Function>(x_tp1, VM.matrixToVector<Function>(Sigma_tp1));
            Function beliefvec = Function.derivative(packed);
            beliefvec.printOperatorCounts();

            bool[] invDynInputVarIndices;
            vars = initializeInputVariables(beliefvec, vars, out invDynInputVarIndices);
            beliefvec.orderVariablesInDomain(vars);

            RuntimeFunction brun = beliefvec.compile();

            double[] result = new double[6], varvals = { -4, 2, 0.214285714285714, -0.285714285714286, 0, 0 };
            brun.eval(result, varvals);
            for (int i = 0; i < 6; ++i) {
                Console.Write("b[" + i + "]: " + result[i] + " ");
            }
            Console.WriteLine();

            Console.Read();
        }
    }
}



        /*
        Function beliefDynamics(Function[,] b_t, Function[,] u_t)
        {

        function b_tp1 = belief_dynamics(b_t, u_t, model)

        dynamics_func = model.dynamics_func;
        obs_func = model.obs_func;
        qDim = model.qDim;
        rDim = model.rDim;

        [x_t, SqrtSigma_t] = decompose_belief(b_t, model);
        Sigma_t = SqrtSigma_t*SqrtSigma_t;

        [x_tp1, Sigma_tp1] = ekf(x_t, Sigma_t, u_t, z_tp1, model);

        % Compute square root for storage
        % Several different choices available -- we use the principal square root
        [V,D] = eigs(Sigma_tp1);
        SqrtSigma_tp1 = V*sqrt(D)*V';

        b_tp1 = compose_belief(x_tp1, SqrtSigma_tp1, model);
        }
        

        

        Function computeCost(Function[][] X, Function[][] U)
        {
            T = model.T;
        belief_cost = zeros(1,T);
        control_cost = zeros(1,T-1);

        b_t = b1;
        for t=1:T-1
            [~, SqrtSigma_t] = decompose_belief(b_t, model);
            belief_cost(t) = model.alpha_belief*sum(sum_square(SqrtSigma_t));
            control_cost(t) = model.alpha_control*sum_square(U(:,t));
            b_t = belief_dynamics(b_t, U(:,t), [], model);
        end
        [~, SqrtSigma_T] = decompose_belief(b_t, model);
        belief_cost(T) = model.alpha_final_belief*sum(sum_square(SqrtSigma_T));

        cost = sum(belief_cost(:)) + sum(control_cost(:));

        }
        */

        
/*

            Function[,] A = new Function[,]{
                    {2.57821520308325,3.00219163401278,-0.423828490485074},
                    {3.00219163401278,8.06195766958555,-1.89719328262036},
                    {-0.423828490485074,-1.89719328262036,0.502255458897483}};
            VM.print<Function>(A);

            Function[,] B = new Function[,]{
                    {1,2,2,1},
                    {1,2,1,2},
                    {2,1,2,1}};
            
            int numvars = 21;
            Variable[] varlist = new Variable[21];
            for (int i = 0; i < numvars; ++i)
            {
                varlist[i] = new Variable("var_" + i);
            }
            
            int idx = 0;
            
            Function[,] A = new Function[3,3];
            for(int i = 0; i < 3; ++i) {
                for(int j = 0; j < 3; ++j) {
                    A[i,j] = varlist[idx++];
                }
            }

            Function[,] B = new Function[3, 4];
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 4; ++j) {
                    B[i, j] = varlist[idx++];
                }
            }

            Program prog = new Program();

            Function[,] F = prog.mdivide(A, B);
            VM.print<Function>(F);
            
            Function.printCompilerSource = true;
            //F[0].compileCCodeToFile("test.c");
          
            int nrows = F.GetLength(0);
            int ncols = F.GetLength(1);
            List<Function> L = new List<Function>();
            for (int i = 0; i < nrows; ++i)
            {
                for (int j = 0; j < ncols; ++j)
                {
                    L.Add(F[i, j]);
                }
            }
            Function G = new Function(L);
            RuntimeFunction grun = G.compile();
           
            Function Fvec = new Function(VM.matrixToVector<Function>(F));
            Fvec.printOperatorCounts();
            //Fvec.orderVariablesInDomain(varlist);
            RuntimeFunction f = Fvec.compile();
            //double[] result = new double[1], vars = {};
            //f.eval(result, vars);
            //Console.WriteLine("result[0]:" + result[0]);
*/

    
//This is example code in Ch 2.7 
//namespace Example_CreatingRuntimeFunction
//{
//    class Program
//    {
//        static void Main(string[] args)
//        {
//            Function.newContext();
//            Variable a = new Variable("a"), b = new Variable("b");
//            Function f = new Function(a * b, Function.sin(b));
//            Function g = Function.D(f, b);
//            Function e = Function.derivative(g);
//            e.lhsName = "e";
//            e.print();

//            Function.printCompilerSource = true;
//            e.compileCCodeToFile("test.c");
//            //RuntimeFunction erun = e.compile();

            //Variable[] a = new Variable[] { new Variable(), new Variable() };
            //Function[] f = new Function[] {a[0] * Function.sqrt(a[1]), Function.sin(a[1]) };
            //Function[,] g1 = Function.D(VM.vectorToMatrix<Function>(f, 1), a);
            ////Function[] g2 = Function.D(f, a[1]);

            //Function[,] e1 = Function.derivative(g1);
            //Console.WriteLine("nrows: " + e1.GetLength(0) + " ncols: " + e1.GetLength(1));
            //VM.print<Function>(e1);
            ////Function e2 = Function.derivative(g2);
            ////e2.print();

//            //double[] result = new double[2], vars = { 1, Math.PI };
//            //erun.eval(result, vars);
//            //Console.WriteLine("e[0]:" + result[0] + " e[1]:" + result[1]);
//            //Console.Read();
//        }
//    }
//}
