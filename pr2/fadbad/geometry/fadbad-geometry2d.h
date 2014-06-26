#ifndef __FADBAD_GEOMETRY2D_H__
#define __FADBAD_GEOMETRY2D_H__

#include "../fadbad-utils.h"

#include <Eigen/Eigen>
using namespace Eigen;

#include <assert.h>

typedef Matrix<bdouble,2,1> Vector2b;
typedef std::vector<Vector2b> StdVector2b;

const bdouble bepsilon = 1e-5;

inline bdouble fadbad_triangle_area(const Vector2b& a, const Vector2b& b, const Vector2b& c) {
	return abs((c(0)*(a(1) - b(1)) + a(0)*(b(1) - c(1)) + b(0)*(c(1) - a(1))) / 2.0);
}

class FadbadTriangle2d {
public:
	Vector2b a, b, c;

	FadbadTriangle2d(const Vector2b& a_pt, const Vector2b& b_pt, const Vector2b& c_pt) : a(a_pt), b(b_pt), c(c_pt) { };

	Vector2b closest_point_to(const Vector2b& x);
	bool is_inside(const Vector2b& x);
	bdouble area();
};


#endif
