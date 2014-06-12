#include <symbolic/casadi.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <cstdlib>

namespace AD = CasADi;
using namespace CasADi;

#define X_DIM 2

/**
 * Attempts so far
 *
 * isEqual(x(0), x(1)) doesn't work as expect
 * I believe it compares the expressions (not the values)
 *
 * x(0).__le__(x(1)) will return 1 if x(0) <= x(1), otherwise 0
 *
 * x(0).if_else_zero(val) returns 0 if x(0) == 0 else returns val
 */

SXMatrix f(SXMatrix& x) {
//	if (isEqual(x(0), x(1))) {
//		return 1;
//	} else {
//		return -1;
//	}

//	return (x(0).__le__(x(1)));

//	return x(0).if_else_zero(99);

//	return sign(x(0) - x(1));

	// no idea what hasNZ does
//	if (x(0).hasNZ(0, 0)) {
//		return 2;
//	} else {
//		return -2;
//	}

	// implement x(0) <= x(1)
	// best I can do is if (x(0) <= x(1)) ? 0 : val
	return (sign(x(0) - x(1)) + 1).if_else_zero(99);
}

SXFunction inequality_function() {
	SXMatrix x = ssym("x", X_DIM);

	SXMatrix result = f(x);

	std::vector<SXMatrix> input;
	input.push_back(x);

	SXFunction f_casadi(input, result);
	f_casadi.init();

	return f_casadi;
}

void test_inequality() {
	SXFunction f_casadi = inequality_function();

	// create input arrays
	double x_arr[X_DIM];

	// fill in input arrays
	x_arr[0] = 1.9;
	x_arr[1] = 2;

	// set function inputs
	f_casadi.setInput(x_arr, 0);

	// evaluate
	f_casadi.evaluate();
	double result;
	f_casadi.getOutput(&result, 0);

	std::cout << "result: " << result << "\n";
}

int main(int argc, char* argv[]) {
	test_inequality();
}
