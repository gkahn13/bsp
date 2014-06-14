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

void test_segment_intersection() {
	vec p0, p1, p0_other, p1_other;
	vec intersection;
	bool is_intersecting;

	{
	p0 << 0 << 0;
	p1 << 0 << 1;
	p0_other << -.5 << .5;
	p1_other << .5 << .5;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "Should be: 0 \t .5\n";
	std::cout << intersection.t() << "\n";
	}

	{
	p0 << 20 << 10;
	p1 << 10 << 5;
	p0_other << 25 << 15;
	p1_other << 15 << 10;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "is_intersecting should be false\n";
	std::cout << "is_intersecting: " << is_intersecting << "\n\n";
	}

	{
	p0 << 3 << -1;
	p1 << 4 << 2;
	p0_other << 4 << 0;
	p1_other << -1 << 0;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "Should be: 3.33 \t 0\n";
	std::cout << intersection.t() << "\n";
	}

	{
	p0 << 0 << 0;
	p1 << 1 << 0;
	p0_other << 1 << 0;
	p1_other << 0 << 1;

	Segment seg0(p0, p1), seg1(p0_other, p1_other);
	is_intersecting = seg0.intersection(seg1, intersection);

	std::cout << "Should be: 1 \t 0\n";
	std::cout << intersection.t() << "\n";
	}

}

std::vector<Segment> test_scope_help() {
	// works because when push_back called
	// creates copy of object that stays with vector scope
	std::vector<Segment> list;
	for(int i=0; i < 1000; ++i) {
		list.push_back(Segment(i*ones<vec>(2), i*ones<vec>(2)));
	}
	return list;
}

void test_scope() {
	std::vector<Segment> list = test_scope_help();

	for(int i=0; i < list.size(); ++i) {
		std::cout << "i: " << i << "\n";
		std::cout << list[i].p0.t();
		std::cout << list[i].p1.t() << "\n";
	}
}


void test_beam() {
	Beam orig_beam({0,0}, {1,1}, {-1,1});
	std::vector<Beam> beams;

	std::cout << "intersects right and left\n";
	beams = orig_beam.truncate(Segment({-1,.5}, {1.5, .75}));
	geometry2d::plot_beams(beams);

	std::cout << "intersects right and top\n";
	beams = orig_beam.truncate(Segment({.3, 1.2}, {1.5, .75}));
	geometry2d::plot_beams(beams);

	std::cout << "intersects top and left\n";
	beams = orig_beam.truncate(Segment({.4, 1.5}, {-.75, .5}));
	geometry2d::plot_beams(beams);

	std::cout << "intersects right\n";
	beams = orig_beam.truncate(Segment({.1,.3}, {1, .4}));
	geometry2d::plot_beams(beams);

	std::cout << "intersects top\n";
	beams = orig_beam.truncate(Segment({-.5,.75}, {.2, 1.2}));
	geometry2d::plot_beams(beams);

	beams = orig_beam.truncate(Segment({.5,.75}, {-.2, 1.2}));
	geometry2d::plot_beams(beams);

	std::cout << "intersects left\n";
	beams = orig_beam.truncate(Segment({-.8,.2}, {.05, .1}));
	geometry2d::plot_beams(beams);

	std::cout << "fully inside\n";
	beams = orig_beam.truncate(Segment({-.3,.4}, {.7, .75}));
	geometry2d::plot_beams(beams);

	std::cout << "fully outside\n";
	beams = orig_beam.truncate(Segment({-1, .25}, {-.8, .7}));
	geometry2d::plot_beams(beams);

	beams = orig_beam.truncate(Segment({.3, 1.25}, {.5, 1.1}));
	geometry2d::plot_beams(beams);

	beams = orig_beam.truncate(Segment({.3, .25}, {.5, .4}));
	geometry2d::plot_beams(beams);
}

void test_segment_distance_to() {
	{
	Segment s({0, 0}, {0, 10});
	vec p = {5, 5};
	double dist = s.distance_to(p);
	vec c = s.closest_point_to(p);

	std::cout << "dist should be 5\n";
	std::cout << "dist: " << dist << "\n";
	std::cout << "closest point should be (0, 5)\n";
	std::cout << "closest point: " << c.t() << "\n";
	}

	{
	Segment s({0, 0}, {10, 10});
	vec p = {13, 14};
	double dist = s.distance_to(p);
	vec c = s.closest_point_to(p);

	std::cout << "dist should be 5\n";
	std::cout << "dist: " << dist << "\n";
	std::cout << "closest point should be (10, 10)\n";
	std::cout << "closest point: " << c.t();
	}
}

void test_segment_closest_point() {
	{
	Segment s({-.0112, -1.1158}, {1.2724, 1.983});
	vec p = {0, 0};
	vec c = s.closest_point_to(p);

	std::cout << "closest point: " << c.t();
	}
}

void test_halfspace() {
	{
	Halfspace h({1,1}, {0,0});
	Segment seg_part({0,0},{0,0});

	{
	Segment s({1,1},{2,1});
	bool c = h.contains_part(s, seg_part);

	std::cout << "should contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part should be (1,1), (2,1)\n";
	std::cout << "seg_part\n" << seg_part.p0.t() << seg_part.p1.t() << "\n";
	}

	{
	Segment s({1,1},{-1,-1});
	bool c = h.contains_part(s, seg_part);

	std::cout << "should contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part should be (0,0), (1,1)\n";
	std::cout << "seg_part\n" << seg_part.p0.t() << seg_part.p1.t() << "\n";
	}

	{
	Segment s({0,0},{1,-1});
	bool c = h.contains_part(s, seg_part);

	std::cout << "should not contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part\n" << seg_part.p0.t() << seg_part.p1.t() << "\n";
	}
	}


	{
	Halfspace h({2,0}, {-2,0});
	Segment seg_part({0,0},{0,0});

	{
	Segment s({-9,2},{-2,5});
	bool c = h.contains_part(s, seg_part);

	std::cout << "should not contain part\n";
	std::cout << "contains part: " << c << "\n";
	std::cout << "seg_part\n" << seg_part.p0.t() << seg_part.p1.t() << "\n";
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
