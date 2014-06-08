#ifndef __PLANAR_UTILS_H__
#define __PLANAR_UTILS_H__

#include <Python.h>

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numpy.hpp>
#include <boost/filesystem.hpp>

namespace py = boost::python;
namespace np = boost::numpy;

#include <armadillo>
using namespace arma;

namespace planar_utils {

np::ndarray arma_to_ndarray(mat& m);

}

#endif
