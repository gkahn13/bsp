#include "../include/planar-utils.h"

namespace planar_utils {

np::ndarray arma_to_ndarray(mat& m) {
	if (m.is_colvec()) {
		py::tuple shape = py::make_tuple(m.n_rows);
		np::dtype dtype = np::dtype::get_builtin<float>();
		np::ndarray n = np::zeros(shape, dtype);
		for(int i=0; i < m.n_rows; ++i) {
			n[py::make_tuple(i)] = m(i);
		}

		return n;
	} else {
		py::tuple shape = py::make_tuple(m.n_rows, m.n_cols);
		np::dtype dtype = np::dtype::get_builtin<float>();
		np::ndarray n = np::zeros(shape, dtype);
		for(int i=0; i < m.n_rows; ++i) {
			for(int j=0; j < m.n_cols; ++j) {
				n[py::make_tuple(i,j)] = m(i,j);
			}
		}

		return n;
	}
}

}
