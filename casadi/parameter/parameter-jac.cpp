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
#define RSCALE 1e4
#define QSCALE 1
#define DT 1.0/5.0


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
const double mass1 = 0.1;
const double mass2 = 0.1; 


//coefficient of friction 
const double b1 = 0.0; 
const double b2 = 0.0; 

const double length1 = 0.15; 
const double length2 = 0.15; 

const double gravity = 9.86; 
}


const double diffEps = 0.0078125 / 16;


const int T = TIMESTEPS;
const double INFTY = 1e10;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const double alpha_belief = 10, alpha_final_belief = 10, alpha_control = 0, alpha_goal_state = 1;
const double alpha_joint_belief = 0, alpha_param_belief = 10,
              alpha_final_joint_belief = 0, alpha_final_param_belief = 10, alpha_goal_joint_state = 10, alpha_goal_param_state = 1;

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



// Jacobians: dh(x,r)/dx, dh(x,r)/dr
void linearizeObservation(const SXMatrix& x,SXMatrix& H)
{
    SXMatrix xr(X_DIM), xl(X_DIM);
    xr = x; 
    xl = x; 
    for (int i = 0; i < X_DIM; ++i) {

        xr[i] += diffEps; xl[i] -= diffEps;
        SXMatrix obs = (obsfunc(xr) - obsfunc(xl))/(2.0*diffEps);
        //cout<<196<<"\n";
        for (int j = 0; j< Z_DIM; ++j){
           H(j,i) = obs[j]; 
        }
        
        xr[i] = x[i]; xl[i] = x[i];
    }


}


inline SXMatrix  getMMT()
{
    SXMatrix  S = SXMatrix(DMatrix::eye(X_DIM));
    S(0,0) = 0.25*0.001 + 0.0000000000001;
    S(1,1) = 0.25*0.001 + 0.0000000000001;
    S(2,2) = 1.0*0.001 + 0.0000000000001;
    S(3,3) = 1.0*0.001 + 0.0000000000001;
    S(4,4) = 0.0005;
    S(5,5) = 0.0005;
    S(6,6) = 0.0001;
    S(7,7) = 0.0001;
    S *= QSCALE;
    return S;
}

// for N = dh(x,r)/dr
// this returns N*~N 
inline SXMatrix getNNT()
{
    SXMatrix S = SXMatrix(DMatrix::eye(Z_DIM));
    S(0,0) = 0.0001;
    S(1,1) = 0.0001;
    S(2,2) = 0.00001;
    S(3,3) = 0.00001;
    S *= RSCALE; 
    return S;
}


// Jacobians: df(x,u,q)/dx, df(x,u,q)/dq
void linearizeDynamics(const SXMatrix& x, const SXMatrix& u, SXMatrix& A)
{
    //g is control input steer angle
    SXMatrix xr(X_DIM), xl(X_DIM);
    SXMatrix dyn(X_DIM);

    xr = x; 
    xl = x; 

    for (size_t i = 0; i < X_DIM; ++i) {

        xr[i] += diffEps; xl[i] -= diffEps;
        dyn = (dynfunc(xr, u) - dynfunc(xl, u)) /(2.0*diffEps);

        for (int j = 0; j< X_DIM; ++j){
           A(j,i) = dyn[j]; 
        }

         xr[i] = x[i]; xl[i] = x[i];
    }


    return;
}


void EKF(const SXMatrix& x_t, const SXMatrix& u_t, const SXMatrix& Sigma_t, SXMatrix& x_tp1, SXMatrix& Sigma_tp1)
{
    SXMatrix A(X_DIM,X_DIM), MMT(X_DIM,X_DIM);
 
    linearizeDynamics(x_t, u_t, A);

    
    //cout << "M" << endl << M << endl;

    MMT = getMMT(); 


    Sigma_tp1 = mul(mul(A,Sigma_t),trans(A)) + MMT;

    //cout << "Sigma_tp1 first" << endl << Sigma_tp1 << endl;

    x_tp1 = dynfunc(x_t, u_t);

    //cout << "x_tp1" << endl << x_tp1 << endl;


    SXMatrix H(Z_DIM,X_DIM), NNT(Z_DIM), R(R_DIM,R_DIM);
    //cout<<292<<"\n";
    linearizeObservation(x_tp1, H);

    //cout <<  << endl << H << endl;
    //cout << "N" << endl << N << endl;


    NNT = getNNT(); 

    //K = ((Sigma_tp1*~H)/(H*Sigma_tp1*~H*delta + RC));
    SXMatrix K = mul(mul(Sigma_tp1, trans(H)), solve(mul(H, mul(Sigma_tp1, trans(H))) + NNT,SXMatrix(DMatrix::eye(Z_DIM))));

    Sigma_tp1 = Sigma_tp1 - mul(K,mul(H,Sigma_tp1));
}



// params(0) = alpha_belief
// params(1) = alpha_control
// params[2] = alpha_final_belief
SXMatrix costfunc(const SXMatrix& XU, const SXMatrix& Sigma_0, const SXMatrix& params)
{
    SXMatrix cost = 0;

    SXMatrix x_tp1(X_DIM,1);
    SXMatrix Sigma_t = Sigma_0, Sigma_tp1(X_DIM,X_DIM);
    SXMatrix x_t(X_DIM,1), u_t(U_DIM,1);

    int offset = 0;

    for (int t = 0; t < (T-1); ++t)
    {

        x_t = XU(Slice(offset,offset+X_DIM));

        offset += X_DIM;
        u_t = XU(Slice(offset,offset+U_DIM));
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

SXFunction casadiLinearizeDynamics() {
    SXMatrix x = ssym("x",X_DIM,1);
    SXMatrix u = ssym("u",U_DIM,1);

    SXMatrix dyn = dynfunc(x,u);

    SXMatrix dyndx = jacobian(dyn,x);

    vector<SXMatrix> inp;
    inp.push_back(x);
    inp.push_back(u);
    
    SXFunction dyn_fcn(inp,dyndx);
    dyn_fcn.init();
    
    generateCode(dyn_fcn,"parameter-dyndx");
    
    return dyn_fcn;
}

void casadiLinearizeObservation() {
    SXMatrix x = ssym("x",X_DIM,1);

    SXMatrix dyn = obsfunc(x);

    SXMatrix dyndx = jacobian(dyn,x);

    vector<SXMatrix> inp;
    inp.push_back(x);
    
    SXFunction dyn_fcn(inp,dyndx);
    dyn_fcn.init();
    
    generateCode(dyn_fcn,"parameter-obsdx");
    
}

int main(int argc, char* argv[])
{
    SXFunction dyn_fcn = casadiLinearizeDynamics();
    casadiLinearizeObservation(); 
 
    /*SXMatrix x = ssym("x",X_DIM,1);
    SXMatrix u = ssym("u",U_DIM,1);

    SXMatrix dyn = dynfunc(x,u);

    SXMatrix dyndx = jacobian(dyn,x);

    vector<SXMatrix> inp;
    inp.push_back(x);
    inp.push_back(u);
    
    SXFunction dyn_fcn(inp,dyndx);
    dyn_fcn.init();
    
    generateCode(dyn_fcn,"parameter-dyndx");
   */
    #define TEST

    // test evaluate function
    #ifdef TEST
    double x0[X_DIM];
    double Sigma0[X_DIM][X_DIM];

    x0[0] = 0; x0[1] = 0; x0[2] = 0;
    x0[3] = 0; x0[4] = 1/0.05;
    x0[5] = 1/0.05; x0[6] = 1/0.12;
    x0[7] = 1/0.13; 

 

    for(int i = 0; i < X_DIM; ++i) {
        for(int j = 0; j < X_DIM; ++j) {
            Sigma0[i][j] = 0;
        }
    }
    Sigma0[4][4] = sqrt(0.5); 
    Sigma0[5][5] = sqrt(0.5); 
    Sigma0[6][6] = 1.0; 
    Sigma0[7][7] = 1.0; 

   

    //double t_r0[1];

    double u0[U_DIM];
    for(int i = 0; i < U_DIM; ++i) {
        u0[i] = 0.0;  
      
    }
 
    

    dyn_fcn.setInput(x0,0);
    dyn_fcn.setInput(u0,1);
  
    dyn_fcn.evaluate();

    double A[X_DIM*X_DIM];
    dyn_fcn.getOutput(A,0);

    //cout << "cost: " << setprecision(12) << cost << endl;
    
    #endif


    return 0;
}
