#include "../system/pr2-sim.h"
#include "../system/pr2-system.h"
#include "../geometry/geometry3d.h"
#include "../utils/rave-utils.h"
#include "../utils/pr2-utils.h"
#include "../../util/Timer.h"

#include "../fadbad/geometry/fadbad-geometry3d.h"
#include "../fadbad/fadbad-utils.h"

#include <openrave-core.h>
namespace rave = OpenRAVE;

TimerCollection tc;

void test_cost() {
	Vector3d vec_d(1,2,3);
	Vector3b vec_b = vec_d.cast<bdouble>();

//	bint val;

//	StdVector2b arr(3, Vector2b::Ones());
//	arr[val];

//	bint arr[3];
//	arr[val];

//	int n = 2;
//	val = static_cast<bint>(n >= 0 ? n + 0.5 : n - 0.5);
}

int main() {
	test_cost();
	return 0;
}
