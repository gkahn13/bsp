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
	Py_Initialize();
	np::initialize();

	py::numeric::array::set_module_and_type("numpy", "ndarray");

	std::string working_dir = boost::filesystem::current_path().normalize().string();
	std::string bsp_dir = working_dir.substr(0,working_dir.find("bsp"));
	std::string planar_dir = bsp_dir + "bsp/planar";

	py::object main_module = py::import("__main__");
	py::object main_namespace = main_module.attr("__dict__");
	py::exec("import sys, os", main_namespace);
	py::exec(py::str("sys.path.append('"+planar_dir+"')"), main_namespace);
	py::object plot_mod = py::import("plot_planar");
	py::object plot_beams = plot_mod.attr("plot_beams");


	std::vector<Beam> beams;
	beams.push_back(Beam({0,0}, {1,1}, {-1,1}));
//	beams.push_back(Beam({0,0}, {0, .5}, {-1, .5}));

	beams = beams[0].truncate(Segment({-.25,.5}, {0, .5}));

	py::list beams_pylist;
	for(int i=0; i < beams.size(); ++i) {
		mat m = join_horiz(beams[i].base, beams[i].a);
		m = join_horiz(m, beams[i].b);
		beams_pylist.append(planar_utils::arma_to_ndarray(m));
	}

	try {
		plot_beams(beams_pylist);
	}
	catch(py::error_already_set const &)
	{
		PyErr_Print();
	}
}

int main(int argc, char* argv[]) {
//	test_segment_intersection();
//	test_scope();
	test_beam();
}
