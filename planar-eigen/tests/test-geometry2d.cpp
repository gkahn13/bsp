#include "../include/planar-utils.h"
#include "../include/geometry2d.h"

#include <Python.h>

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;
namespace np = boost::numpy;

#include <Eigen/Eigen>
#include <Eigen/StdVector>
using namespace Eigen;

template <size_t _dim0, size_t _dim1>
using mat = Matrix<double, _dim0, _dim1>;

template <size_t _dim>
using vec = Matrix<double, _dim, 1>;


void test_segment_intersection() {
	vec<2> p0, p1, p0_other, p1_other;
	vec<2> intersection;
	bool is_intersecting;

	{
	p0 << 0 , 0;
	p1 << 0 , 1;
	p0_other << -.5 , .5;
	p1_other << .5 , .5;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "Should be: 0 \t .5\n";
	std::cout << intersection.transpose() << "\n\n";
	}

	{
	p0 << 20 , 10;
	p1 << 10 , 5;
	p0_other << 25 , 15;
	p1_other << 15 , 10;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "is_intersecting should be false\n";
	std::cout << "is_intersecting: " << is_intersecting << "\n\n";
	}

	{
	p0 << 3 , -1;
	p1 << 4 , 2;
	p0_other << 4 , 0;
	p1_other << -1 , 0;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "Should be: 3.33 \t 0\n";
	std::cout << intersection.transpose() << "\n\n";
	}

	{
	p0 << 0 , 0;
	p1 << 1 , 0;
	p0_other << 1 , 0;
	p1_other << 0 , 1;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "Should be: 1 \t 0\n";
	std::cout << intersection.transpose() << "\n\n";
	}

}

std::vector<Segment> test_scope_help() {
	// works because when push_back called
	// creates copy of object that stays with vector scope
	std::vector<Segment> list;
	for(int i=0; i < 1000; ++i) {
		list.push_back(Segment(i*vec<2>::Ones(), i*vec<2>::Ones()));
	}
	return list;
}

void test_scope() {
	std::vector<Segment> list = test_scope_help();

	for(int i=0; i < list.size(); ++i) {
		std::cout << "i: " << i << "\n";
		std::cout << list[i].p0.transpose() << "\n";
		std::cout << list[i].p1.transpose() << "\n\n";
	}
}


void test_beam() {
	vec<2> base, a, b, p0, p1;
	base << 0 , 0;
	a << 1 , 1;
	b << -1 , 1;
	Beam orig_beam(base, a, b);
	std::vector<Beam> beams;

	std::cout << "intersects right and left\n";
	p0 << -1, .5;
	p1 << 1.5, .75;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	std::cout << "intersects right and top\n";
	p0 << .3, 1.2;
	p1 << 1.5, .75;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	std::cout << "intersects top and left\n";
	p0 << .4, 1.5;
	p1 << -.75, .5;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	std::cout << "intersects right\n";
	p0 << .1, .3;
	p1 << 1, .4;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	std::cout << "intersects top\n";
	p0 << -.5, .75;
	p1 << .2, 1.2;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	p0 << .5, .75;
	p1 << -.2, 1.2;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	std::cout << "intersects left\n";
	p0 << -.8, .2;
	p1 << .05, .1;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	std::cout << "fully inside\n";
	p0 << -.3, .4;
	p1 << .7, .75;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	std::cout << "fully outside\n";
	p0 << -1, .25;
	p1 << -.8, .7;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	p0 << .3, 1.25;
	p1 << .5, 1.1;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);

	p0 << .3, .25;
	p1 << .5, .4;
	beams = orig_beam.truncate(Segment(p0, p1));
	geometry2d::plot_beams(beams);
}

void test_segment_distance_to() {
	vec<2> p, p0, p1;

	{
	p0 << 0, 0;
	p1 << 0, 10;
	Segment s(p0, p1);
	p << 5, 5;
	double dist = s.distance_to(p);
	vec<2> c = s.closest_point_to(p);

	std::cout << "dist should be 5\n";
	std::cout << "dist: " << dist << "\n";
	std::cout << "closest point should be (0, 5)\n";
	std::cout << "closest point: " << c.transpose() << "\n\n";
	}

	{
	p0 << 0, 0;
	p1 << 10, 10;
	Segment s(p0, p1);
	p << 13, 14;
	double dist = s.distance_to(p);
	vec<2> c = s.closest_point_to(p);

	std::cout << "dist should be 5\n";
	std::cout << "dist: " << dist << "\n";
	std::cout << "closest point should be (10, 10)\n";
	std::cout << "closest point: " << c.transpose() << "\n\n";
	}
}

void test_segment_closest_point() {
	{
	vec<2> p0, p1, p;
	p0 << -.0112, -1.1158;
	p1 << 1.2724, 1.983;
	Segment s(p0, p1);
	p << 0, 0;
	vec<2> c = s.closest_point_to(p);

	std::cout << "closest point: " << c.transpose() << "\n\n";
	}
}

void test_halfspace() {
	vec<2> p0, p1;

	{
	Halfspace h(vec<2>::Ones(), vec<2>::Zero());
	Segment seg_part(vec<2>::Zero(), vec<2>::Zero());

	{
	p0 << 1, 1;
	p1 << 2, 1;
	Segment s(p0, p1);
	bool c = h.contains_part(s, seg_part);

	std::cout << "should contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part should be (1,1), (2,1)\n";
	std::cout << "seg_part\n" << seg_part.p0.transpose() << "\n" << seg_part.p1.transpose() << "\n\n";
	}

	{
	p0 << 1, 1;
	p1 << -1, -1;
	Segment s(p0, p1);
	bool c = h.contains_part(s, seg_part);

	std::cout << "should contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part should be (0,0), (1,1)\n";
	std::cout << "seg_part\n" << seg_part.p0.transpose() << "\n" << seg_part.p1.transpose() << "\n\n";
	}

	{
	p0 << 0, 0;
	p1 << 1, -1;
	Segment s(p0, p1);
	bool c = h.contains_part(s, seg_part);

	std::cout << "should not contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part\n" << seg_part.p0.transpose() << "\n" << seg_part.p1.transpose() << "\n\n";
	}
	}


	{
	p0 << 2, 0;
	p1 << -2, 0;
	Halfspace h(p0, p1);
	Segment seg_part(vec<2>::Zero(), vec<2>::Zero());

	{
	p0 << -9, 2;
	p1 << -2, 5;
	Segment s(p0, p1);
	bool c = h.contains_part(s, seg_part);

	std::cout << "should not contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part\n" << seg_part.p0.transpose() << "\n" << seg_part.p1.transpose() << "\n\n";
	}
	}

}

int main(int argc, char* argv[]) {
//	test_segment_intersection();
//	test_scope();
//	test_beam();
//	test_segment_distance_to();
	test_segment_closest_point();
//	test_halfspace();
}
