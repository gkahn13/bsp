#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>



// horizon is total lifetime of planning
// timesteps is how far into future accounting for during MPC
#define HORIZON 500
#define TIMESTEPS 15

#define DT 0.1


// for ILQG
//#define TIMESTEPS 500
//#define DT 1.0/100.0

// J_DIM == SIM_X_DIM in original file
#define J_DIM 4 // number of joints in state (2 position and 2 velocity)
#define K_DIM 4 // number of parameters in state (2 masses and 2 lengths)

#define X_DIM 8
#define U_DIM 2
#define Z_DIM 4
#define Q_DIM 8
#define R_DIM 4

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM+S_DIM)
#define XU_DIM (X_DIM*T+U_DIM*(T-1))

namespace dynamics {
const double mass1 = 0.5;
const double mass2 = 0.5;

//coefficient of friction 
const double b1 = 0.0; 
const double b2 = 0.0; 

const double length1 = 0.5;
const double length2 = 0.5;

const double gravity = 9.82;
}


const double diffEps = 0.0078125 / 16;


const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


double *inputVars, *vars;
std::vector<int> maskIndices;



using namespace CasADi;
using namespace std;



SXMatrix jointdynfunc(const SXMatrix& x, const SXMatrix& u, SXMatrix& l1, SXMatrix& l2, SXMatrix& m1, SXMatrix& m2)
{

    SXMatrix I1(1), I2(1);
    I1 = m1*(l1*l1)/12; 
    I2 = m2*(l2*l2)/12;

    SXMatrix A(U_DIM,U_DIM); 

    A(0,0) = l1*l1*(0.25*m2+m1)+I1; 

    A(0,1) = 0.5*m2*l1*l2*cos(x(0)-x(1));

    A(1,0) = 0.5*m2*l1*l2*cos(x(0)-x(1));

    A(1,1) = l2*l2*0.25*m2+I2; 

    SXMatrix B(U_DIM,1); 

    ////cout<<B(0,0)<<B(1,0)<<"\n";

    B(0,0) = dynamics::gravity*l1*sin(x(0))*(0.5*m1+m2) - 0.5*m2*l1*l2*(x(3)*x(3)*sin(x(0)-x(1)))+u(0)-dynamics::b1*x(2);
    //cout<<97<<"\n";
    ////cout<<B(0,0)<<"\n";
    B(1,0) = 0.5*m2*l2*(l1*(x[2]*x[2])*sin(x(0)-x(1))+dynamics::gravity*sin(x(1))) + u(1)-dynamics::b2*x(3);
   
    //cout<<100<<"\n";
    SXMatrix xd(U_DIM);
    xd = solve(A,B);

    //cout<<104<<"\n";
    SXMatrix jNew(J_DIM,1) ;

    jNew(0,0) = x(2);
    jNew(1,0) = x(3);
    jNew(2,0) = xd(0);
    jNew(3,0) = xd(1);

    return jNew;

}


// for both joints and params
SXMatrix dynfunc(const SXMatrix& x, const SXMatrix& u)
{
    // RK4 integration
    SXMatrix k1(J_DIM), k2(J_DIM), k3(J_DIM), k4(J_DIM), jinit(J_DIM);

    jinit = x(Slice(0,J_DIM));

    SXMatrix x4(1),x5(1),x6(1),x7(1);

    x4 = x(4);
    x5 = x(5);
    x6 = x(6);
    x7 = x(7);

    SXMatrix length1(1), length2(1), mass1(1), mass2(1);
    length1 = 1/x4;
    length2 = 1/x5;
    mass1 = 1/x6;
    mass2 = 1/x7;
    //cout<<136<<"\n";
    k1 = jointdynfunc(jinit, u, length1, length2, mass1, mass2);
    //cout<<138<<"\n";
    k2 = jointdynfunc(jinit + 0.5*DT*k1, u, length1, length2, mass1, mass2);
    k3 = jointdynfunc(jinit + 0.5*DT*k2, u, length1, length2, mass1, mass2);
    k4 = jointdynfunc(jinit + DT*k3, u, length1, length2, mass1, mass2);

    SXMatrix xNew(X_DIM,1);
    jinit = jinit + DT*(k1 + 2.0*(k2 + k3) + k4)/6.0;
    xNew(0,0) = jinit(0);
    xNew(1,0) = jinit(1);
    xNew(2,0) = jinit(2);
    xNew(3,0) = jinit(3);
    xNew(4,0) = x4;
    xNew(5,0) = x5;
    xNew(6,0) = x6;
    xNew(7,0) = x7;
    //cout<<152<<"\n";
    return xNew;
}


// Observation model
SXMatrix obsfunc(const SXMatrix & x)
{
    SXMatrix z(Z_DIM,1);



    SXMatrix length1(1), length2(1);
    length1 = 1/x(4);
    length2 = 1/x(5);

    SXMatrix cosx0 = cos(x(0));
    SXMatrix sinx0 = sin(x(0));
    SXMatrix cosx1 = cos(x(1));
    SXMatrix sinx1 = sin(x(1));

    z(0,0) = length1*cosx0 + length2*cosx1;
    z(1,0) = length1*sinx0 + length2*sinx1;
    z(2,0) = x(2);
    z(3,0) = x(3);


    return z;
}


inline SXMatrix  varQ()
{
    SXMatrix  S = SXMatrix(DMatrix::eye(X_DIM));
    S(0,0) = 0.01;
    S(1,1) = 0.01;
    S(2,2) = 0.01;
    S(3,3) = 0.01;
    S(4,4) = 1e-6;
    S(5,5) = 1e-6;
    S(6,6) = 1e-6;
    S(7,7) = 1e-6;
    return S;
}

// for N = dh(x,r)/dr
// this returns N*~N
inline SXMatrix varR()
{
    SXMatrix S = SXMatrix(DMatrix::eye(Z_DIM));
    S(0,0) = 1e-2;
        S(1,1) = 1e-2;
        S(2,2) = 1e-3;
        S(3,3) = 1e-3;
    return S;
}


void EKF(const SXMatrix& x_t, const SXMatrix& u_t, const SXMatrix& Sigma_t, SXMatrix& x_tp1, SXMatrix& Sigma_tp1)
{  
    SXMatrix A(X_DIM,X_DIM), MMT(X_DIM,X_DIM);


    SXMatrix dyn = dynfunc(x_t,u_t);

    A = jacobian(dyn,x_t);


    Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + varQ();
    x_tp1 = dynfunc(x_t, u_t);
  

    SXMatrix H(Z_DIM,X_DIM), NNT(Z_DIM), R(R_DIM,R_DIM);

    //Caclulate H using previous Casadi 
    SXMatrix obs = obsfunc(x_t);

    H = jacobian(obs,x_t);

    SXMatrix K = mul(mul(Sigma_tp1, trans(H)), solve(mul(H, mul(Sigma_tp1, trans(H))) + varR(),SXMatrix(DMatrix::eye(Z_DIM))));
  
    Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));
}


// params(0) = alpha_belief
// params(1) = alpha_control
// params[2] = alpha_final_belief
SXMatrix costfunc(const SXMatrix& X, const SXMatrix& U, const SXMatrix& Sigma_0, const SXMatrix& params)
{
    SXMatrix cost = 0;

    SXMatrix x_tp1(X_DIM,1);
    SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
    SXMatrix x_t(X_DIM,1), u_t(U_DIM,1);

    int offset = 0;

    for (int t = 0; t < (T-1); ++t)
    {

        x_t = X(Slice(t*X_DIM,(t+1)*X_DIM));

        u_t = U(Slice(t*U_DIM,(t+1)*U_DIM));
        offset += U_DIM;

        cost += params(0)*trace(Sigma_t);
        cost += params(1)*inner_prod(u_t, u_t);

        EKF(x_t, u_t, Sigma_t, x_tp1, Sigma_tp1);
        Sigma_t = Sigma_tp1;
    }

    cost += params[2]*trace(Sigma_t);

    return cost;
}




void generateCode(FX fcn, const std::string& name){
    cout << "Generating code for " << name << endl;

    fcn.generateCode(name + ".c");
}


int main(int argc, char* argv[])
{
    
    cout << "Creating casadi file for T = " << T << endl;

    //vector<SXMatrix> X, U;
    int nXU = T*X_DIM+(T-1)*U_DIM;
    SXMatrix X = ssym("X",T*X_DIM,1);
    SXMatrix U = ssym("U",(T-1)*U_DIM,1);
    SXMatrix Sigma_0 = ssym("S0",X_DIM,X_DIM);
    SXMatrix params = ssym("params",3); // alpha_control, alpha_belief, alpha_final_belief
    cout <<"before generate"<<"\n";
    // Objective
    SXMatrix f = costfunc(X, U, Sigma_0, params);

    SXMatrix grad_f = gradient(f,U);
    cout <<"before generate"<<"\n";
  

    // Create functions

    
    vector<SXMatrix> inp;
    inp.push_back(X);
    inp.push_back(U);
    inp.push_back(Sigma_0);
    inp.push_back(params);

    SXFunction f_fcn(inp,f);
    f_fcn.init();

    vector<SXMatrix> out;
    out.push_back(f);
    out.push_back(grad_f);
    SXFunction grad_f_fcn(inp,out);
    grad_f_fcn.init();

    //SXFunction hess_f_fcn(inp,diag_hess_f);
    //hess_f_fcn.init();

    // Generate code
    generateCode(f_fcn,"parameter-controls-cost");
    generateCode(grad_f_fcn,"parameter-controls-grad");
    //generateCode(hess_f_fcn,"slam-state-diag-hess");
    //#define TEST
    #ifdef TEST
    double length1_est = .05, // inverse = 20
            length2_est = .05, // inverse = 20
            mass1_est = .12, // inverse = 9.52
            mass2_est = .13; // inverse = 11.24

    SXMatrix x0(X_DIM,1); 
    SXMatrix u0(U_DIM,1);
    // position, then velocity
    x0(0,0) = -M_PI/2.0; x0(1,0) = -M_PI/2.0; x0(2,0) = 0; x0(3,0) = 0;
    // parameter start estimates (alphabetical, then numerical order)
    x0(4,0) = 1/length1_est; x0(5,0) = 1/length2_est; x0(6,0) = 1/mass1_est; x0(7,0) = 1/mass2_est;

    u0(0,0) = 0.1; 
    u0(1,0) = 0.1; 

    SXMatrix Sigma_t(X_DIM,X_DIM);

    Sigma_t(4,4) = sqrt(0.5);
    Sigma_t(5,5) = sqrt(0.5);
    Sigma_t(6,6) = 1.0;
    Sigma_t(7,7) = 1.0;

    SXMatrix Sigma_tp1(X_DIM,X_DIM);
    SXMatrix xtp1(X_DIM,1); 
    //Fill in X0
    
    double t_x0[nXU];
    double t_x1[X_DIM*X_DIM];
    double t_x2[3];

    //Fill in X and U
    for (int t = 0; t < T-1; ++t) {

        /*for(int i=0; i<X_DIM; i++){
            t_x0(i,0) = x0(i,0);
        }
        for(int i=0; i<U_DIM; i++){
            t_x0(i,0) = u0(i,0);
        }
        */
        cout<<364<<"\n";
        EKF(x0, u0, Sigma_t, xtp1, Sigma_tp1);
        x0=xtp1;
        Sigma_t =Sigma_tp1;
    }

  /*  for(int i=0; i<X_DIM; i++){
            t_x0(i,0) = x0(i,0);
    }

    //Fill in Alpha

   t_x2[0] = 10; t_x2[1]=0; t_x2[2] = 10; 
    */



    #endif



    return 0;
}
