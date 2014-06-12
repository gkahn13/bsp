#include <iostream>

#include "badiff.h"
using namespace fadbad;

B<double> func(const B<double>& x, const B<double>& y) {
    B<double> z=sqrt(x);
    return y*z+sin(z);
}

int main(int argc, char* argv[]) {
    B<double> x,y,f;    // Declare variables x,y,f
    x=1;                // Initialize variable x
    y=2;                // Initialize variable y
    f=func(x,y);        // Evaluate function and record DAG
    f.diff(0,1);        // Differentiate f (index 0 of 1)
    double fval=f.x();  // Value of function
    double dfdx=x.d(0); // Value of df/dx (index 0 of 1)
    double dfdy=y.d(0); // Value of df/dy (index 0 of 1)

    std::cout << "f(x,y)=" << fval << "\n";
    std::cout << "df/dx(x,y)=" << dfdx << "\n";
    std::cout << "df/dy(x,y)=" << dfdy << "\n";
}
