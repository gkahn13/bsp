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
    public class PointControl : Point
    {

        public PointControl(int timesteps) : base(timesteps) { }

        Function costFunc(Function[][] U, Function[] q, Function[] r, Function[] x_0, Function[] x_goal, Function[,] Sigma_0, Function alpha_belief, Function alpha_control, Function alpha_final_belief, Function alpha_goal_state)
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
        

        void computeCostGradDiagHess(string eval_name)
        {
            uint flag = 0;
            flag |= ((uint)(Convert.ToUInt16(eval_name.Contains("Cost"))) << 0);
            flag |= ((uint)(Convert.ToUInt16(eval_name.Contains("Grad"))) << 1);
            flag |= ((uint)(Convert.ToUInt16(eval_name.Contains("DiagHess"))) << 3);
            
            // variable instantiations
            int nparams = 4;
            int nvars = (T - 1) * UDIM + QDIM + RDIM + XDIM + XDIM + (XDIM * XDIM) + nparams;

            Variable[] vars = new Variable[nvars];
            for (int i = 0; i < nvars; ++i) { vars[i] = new Variable("vars_" + i); }

            U = new Function[T - 1][];

            // initialize vars
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

            Function cost = costFunc(U, q, r, x_0, x_goal, Sigma_0, alpha_belief, alpha_control, alpha_final_belief, alpha_goal_state);

            /*Function[] costJacDiagHessFunc = new Function[2*Uvec.Length + 1];
            costJacDiagHessFunc[0] = cost;
            idx = 1;
            for (int i = 0; i < Uvec.Length; ++i)
            {
                costJacDiagHessFunc[idx++] = Function.D(cost, Uvec[i]);
            }
            for (int i = 0; i < Uvec.Length; ++i)
            {
                costJacDiagHessFunc[idx++] = Function.D(cost, Uvec[i], Uvec[i]);
            }*/
            List<Function> costJacDiagHessFuncList = new List<Function>();
            if ((flag & COMPUTE_COST) != 0)
            {
                costJacDiagHessFuncList.Add(cost);
            }
            
            if ((flag & COMPUTE_JACOBIAN) != 0)
            {
	            for (int i = 0; i < Uvec.Length; ++i)
	            {
	                costJacDiagHessFuncList.Add(Function.D(cost, Uvec[i]));
	            }
            }
            if ((flag & COMPUTE_DIAGONAL_HESSIAN) != 0)
            {
	            for (int i = 0; i < Uvec.Length; ++i)
	            {
	                costJacDiagHessFuncList.Add(Function.D(cost, Uvec[i], Uvec[i]));
	            }
            }
            
            
            Function costJacDiagHess = Function.derivative(costJacDiagHessFuncList.ToArray());
            costJacDiagHess.printOperatorCounts();

            bool[] inputVarIndices;
            Variable[] costJacDiagHessVars = initializeInputVariables(costJacDiagHess, vars, out inputVarIndices);
            costJacDiagHess.orderVariablesInDomain(costJacDiagHessVars);

            

            costJacDiagHess.compileCCodeToFile("control-symeval-"+T+".c");

            string fileName = "control-masks-"+T+".txt";
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
            // T = args[0]
            int T = 15;
            string eval_name = "CostGradDiagHess"; // default to comput cost, gradient, and diagonal hessian
           
            if (args.Length >= 1) {
                T = int.Parse(args[0]);
            }
            if (args.Length >= 2) {
                eval_name = args[1];
            }
           
            Console.WriteLine("Creating files for T = "+T);
           
            Function.newContext();
            Function.printCompilerSource = false;

            PointControl pc = new PointControl(T);

            var stopwatch = new Stopwatch();
            stopwatch.Start();
            
            pc.computeCostGradDiagHess(eval_name);

            stopwatch.Stop();
            
            Console.WriteLine("Finished in " + (stopwatch.ElapsedMilliseconds/1000.0) + " s");
        }
    }
}