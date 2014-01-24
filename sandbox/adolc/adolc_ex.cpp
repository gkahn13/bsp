

#include <cmath>
#include <float.h>
#include <iostream>
#include <assert.h>
#include <adolc/adolc.h>

#include </home/gkahn/Berkeley/Research/bsp/cpp/util/amatrix.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;

void initMatrix(matrix<adouble>& x, matrix<double>& xp)
{
	for(int i=0; i < x.size1(); ++i) {
		for(int j=0; j < x.size2(); ++j) {
			x(i,j) <<= xp(i,j);
		}
	}
}

void init(adouble* x, double *xp, int n)
{
	for(int i=0; i < n; ++i)
	{
		x[i] <<= xp[i];
	}
}

/*
// converted
// joint angles -> end effector position
matrix<adouble> g(const matrix<adouble>& x)
{
    adouble a0 = x(0,0), a1 = x(1,0), a2 = x(2,0), a3 = x(3,0), a4 = x(4,0), a5 = x(5,0);

    matrix<adouble> p(G_DIM, 1);
    p(0,0) = sin(a0) * (cos(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + sin(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2)) + cos(a0) * sin(a3) * sin(a4) * l4;
    p(1,0) = -sin(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + cos(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2) + l1;
    p(2,0) = cos(a0) * (cos(a1) * (sin(a2) * (cos(a4) * l4 + l3) + cos(a2) * cos(a3) * sin(a4) * l4) + sin(a1) * (cos(a2) * (cos(a4) * l4 + l3) - sin(a2) * cos(a3) * sin(a4) * l4 + l2)) - sin(a0) * sin(a3) * sin(a4) * l4;

    return p;
}*/

adouble matrixFunc(matrix<adouble>& x, matrix<double>& xp)
{
	initMatrix(x, xp);
	matrix<adouble> result = prod(trans(x), x);
	return result(0,0);
}

adouble func(adouble* x, double *xp, int n)
{
	init(x, xp, n);

	adouble total = 0;
	for(int i=0; i < n; ++i)
	{
		total += x[i]*x[i];
	}
	return total;
}

int main(int argc, char* argv[])
{

	Matrix<2,2> M = identity<2>();
	Matrix<2,2> aM = identity<2>();

	aM = aM + M;

	std::cout << aM(0,0).value() << std::endl;
	return 0;

	int n,i,j;
	int tape_stats[STAT_SIZE];


	n = 2;

	matrix<double> xpM(n,1), xpM2(n,1);
	double ypM, ypM2, tmp_double;
	matrix<adouble> xM(n,1), xM2(n,1);
	adouble yM, yM2, tmp_adouble;

	sqrt(yM);

	xpM(0,0) = 2; xpM(1,0) = 3;
	xpM2(0,0) = 2; xpM2(1,0) = 3;

	std::cout << "xp: " << xpM << std::endl;

	trace_on(1);
	yM = matrixFunc(xM, xpM);
	tmp_adouble = yM;
	yM >>= ypM;

	trace_on(2);
	yM2 = matrixFunc(xM2, xpM2);
	yM2 >>= ypM2;
	trace_off(2);

	trace_off(1);

	std::cout <<"tmp_adouble " << ypM << std::endl;

	tapestats(1,tape_stats);             // reading of tape statistics
	cout<<"maxlive "<<tape_stats[NUM_MAX_LIVES]<<"\n";
	// ..... print other tape stats

	//double* g = new double[n];
	//gradient(1,n,xp,g);                  // gradient evaluation

	double* gM = new double[n];
	double* xpM_arr = new double[n];
	std::copy(xpM.begin1(), xpM.end1(), xpM_arr);
	gradient(1,n,xpM_arr,gM);

	double* gM2 = new double[n];
	double* xpM2_arr = new double[n];
	std::copy(xpM2.begin1(), xpM2.end1(), xpM2_arr);
	gradient(1,n,xpM2_arr,gM2);

	//std::cout << "g: ";
	//for(int i=0; i < n; ++i) { std::cout << g[i] << " "; }
	//std::cout << std::endl;

	std::cout << "gM: ";
	for(int i=0; i < n; ++i) { std::cout << gM[i] << " "; }
	std::cout << std::endl;

	std::cout << "gM2: ";
	for(int i=0; i < n; ++i) { std::cout << gM2[i] << " "; }
	std::cout << std::endl;

	double** H   = (double**) malloc(n*sizeof(double*));
	for(i=0; i<n; i++)
		H[i] = (double*)malloc((i+1)*sizeof(double));
	hessian(1,n,xpM_arr,H);                   // H equals (n-1)g since g is
}


/*
int m = 1;
double** J = new double*[n];
for(int i=0; i < m; ++i) { J[i] = new double[1]; }
jacobian(1,m,n,xp,J);

std::cout << "J" << std::endl;
for(int i=0; i < n; ++i) {
	for(int j=0; j < m; ++j) {
		std::cout << J[i][j] << " ";
	}
	std::cout << std::endl;
}
*/
