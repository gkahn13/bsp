#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <cmath>
#include <float.h>
#include <iostream>
#include <assert.h>
#include <limits>
#include <iomanip>

template <size_t _numRows, size_t _numColumns = 1> class Matrix {

private:
	double _elems[_numRows * _numColumns];

public:
	// Retrieval
	inline size_t numRows() const { 
		return _numRows; 
	}
	inline size_t numColumns() const { 
		return _numColumns; 
	}

	// Subscript operator
	__forceinline double& operator () (size_t row, size_t column) {
		assert(row < _numRows && column < _numColumns);
		return _elems[row * _numColumns + column]; 
	}
	__forceinline double  operator () (size_t row, size_t column) const {
		assert(row < _numRows && column < _numColumns);
		return _elems[row * _numColumns + column]; 
	}

	__forceinline double& operator [] (size_t elt) {
		assert(elt < _numRows * _numColumns);
		return _elems[elt]; 
	}
	__forceinline double  operator [] (size_t elt) const {
		assert(elt < _numRows * _numColumns);
		return _elems[elt]; 
	}

	// Reset to zeros
	inline void reset() { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] = double(0);
		}
	}

	// Submatrix
	template <size_t numRows, size_t numColumns>
	inline Matrix<numRows, numColumns> subMatrix(size_t row, size_t column) const {
		assert(row + numRows <= _numRows && column + numColumns <= _numColumns);
		Matrix<numRows, numColumns> m;
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				m(i, j) = (*this)(row + i, column + j);
			}
		}
		return m;
	}

	template <size_t numRows>
	inline Matrix<numRows> subMatrix(size_t row, size_t column) const {
		assert(row + numRows <= _numRows && column < _numColumns);
		Matrix<numRows> m;
		for (size_t i = 0; i < numRows; ++i) {
			m[i] = (*this)(row + i, column);
		}
		return m;
	}

	inline Matrix<_numRows> column(size_t columnNr) const {
		assert(columnNr < _numColumns);
		Matrix<_numRows> m;
		for (size_t i = 0; i < _numRows; ++i) {
			m[i] = (*this)(i, columnNr);
		}
		return m;
	}

	inline Matrix<1, _numColumns> row(size_t rowNr) const {
		assert(rowNr < _numRows);
		Matrix<1, _numColumns> m;
		for (size_t i = 0; i < _numColumns; ++i) {
			m[i] = (*this)(rowNr, i);
		}
		return m;
	}

	// Insert
	template <size_t numRows, size_t numColumns>
	inline void insert(size_t row, size_t column, const Matrix<numRows, numColumns>& q) {
		assert(row + numRows <= _numRows && column + numColumns <= _numColumns);
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				(*this)(row + i, column + j) = q(i, j);
			}
		}
	}

	// Matrix addition
	inline const Matrix<_numRows, _numColumns>& operator+=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] += q._elems[i];
		}
		return *this; 
	}

	// Matrix subtraction
	inline const Matrix<_numRows, _numColumns>& operator-=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] -= q._elems[i];
		}
		return *this; 
	}

	// Scalar multiplication
	inline const Matrix<_numRows, _numColumns>& operator*=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] *= a;
		}
		return *this;
	}

	// Scalar division
	inline const Matrix<_numRows, _numColumns>& operator/=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] /= a;
		}
		return *this;
	}
};

template <size_t _size> class SymmetricMatrix {
private:
	double _elems[((_size+1)*_size)/2];

public:
	// Retrieval
	inline size_t numRows() const { 
		return _size; 
	}
	inline size_t numColumns() const { 
		return _size; 
	}

	// Subscript operator
	__forceinline double& operator () (size_t row, size_t column) {
		assert(row < _size && column < _size);
		if (row >= column) {
			return _elems[_size * column + row - ((column + 1)*column) / 2];
		} else {
			return _elems[_size * row + column - ((row + 1)*row) / 2];
		}
	}
	__forceinline double operator () (size_t row, size_t column) const {
		assert(row < _size && column < _size);
		if (row >= column) {
			return _elems[_size * column + row - ((column + 1)*column) / 2];
		} else {
			return _elems[_size * row + column - ((row + 1)*row) / 2];
		}
	}

	__forceinline double& operator [] (size_t elt) {
		assert(elt < ((_size+1)*_size)/2);
		return _elems[elt]; 
	}
	__forceinline double  operator [] (size_t elt) const {
		assert(elt < ((_size+1)*_size)/2);
		return _elems[elt]; 
	}

	inline void reset() { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] = double(0);
		}
	}

	// Matrix addition
	inline const SymmetricMatrix<_size>& operator+=(const SymmetricMatrix<_size>& M) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] += M._elems[i];
		}
		return *this;
	}

	// Matrix subtraction
	inline const SymmetricMatrix<_size>& operator-=(const SymmetricMatrix<_size>& M) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] -= M._elems[i];
		}
		return *this;
	}

	// Scalar multiplication
	inline const SymmetricMatrix<_size>& operator*=(double a) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] *= a;
		}
		return *this;
	}

	// Scalar division
	inline const SymmetricMatrix<_size>& operator/=(double a) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] /= a;
		}
		return *this;
	}

	// Extract symmetric subMatrix starting from a given diagonal element 
	template <size_t _numRows>
	inline SymmetricMatrix<_numRows> subSymmetricMatrix(size_t diag) const {
		assert(diag + _numRows <= _size);
		SymmetricMatrix<_numRows> L;
		size_t index = -1;
		for (size_t i = 0; i < _numRows; ++i) {
			for(size_t j = 0; j <= i; ++j) {
				L[++index] = (*this)(diag + i, diag + j);
			}
		}
		return L;
	}

	// Cast
	inline operator Matrix<_size, _size>() const {
		Matrix<_size,_size> M;
		for (size_t i = 0; i < _size; ++i) {
			for (size_t j = 0; j < _size; ++j) {
				M(i,j) = (*this)(i,j);
			}
		}
		return M;
	}
};

// Unary minus
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator-(const Matrix<_numRows, _numColumns>& M) {
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = -M[i];
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator-(const SymmetricMatrix<_size>& M) {
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = -M[i];
	}
	return L;
}

// Unary plus
template <size_t _numRows, size_t _numColumns>
inline const Matrix<_numRows, _numColumns>& operator+(const Matrix<_numRows, _numColumns>& M) { 
	return M; 
}

template <size_t _size>
inline const SymmetricMatrix<_size>& operator+(const SymmetricMatrix<_size>& M) { 
	return M; 
}

// Transpose
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numColumns, _numRows> operator~(const Matrix<_numRows, _numColumns>& M) {
	Matrix<_numColumns, _numRows> L;
	for (size_t i = 0; i < _numColumns; ++i) {
		for (size_t j = 0; j < _numRows; ++j) {
			L(i, j) = M(j, i);
		}
	}
	return L;
}


// Eigen decomposition of a real symmetric matrix by converting to tridiagonal form followed by QL iteration
// From: Numerical recipes in C++ (3rd edition)
/*
template <size_t _size>
inline void jacobi2(const Matrix<_size, _size>&	q, Matrix<_size, _size>& V, Matrix<_size, _size>& D) 
{
	// Initializations
	Matrix<_size, _size> z = q;
	Matrix<_size> d, e;
	d.reset();
	e.reset();

	int l,k,j,i,m,iter;
	double scale,hh,h,g,f,s,r,p,dd,c,b,absf,absg,pfg;

	//Convert to tridiagonal form
	for (i=_size-1;i>0;i--) 
	{
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++)
				scale += abs(z(i,k));
			if (scale == 0.0)
				e[i]=z(i,l);
			else {
				for (k=0;k<i;k++) {
					z(i,k) /= scale;
					h += z(i,k)*z(i,k);
				}
				f=z(i,l);
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z(i,l)=f-g;
				f=0.0;
				for (j=0;j<i;j++) {
					z(j,i)=z(i,j)/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += z(j,k)*z(i,k);
					for (k=j+1;k<i;k++)
						g += z(k,j)*z(i,k);
					e[j]=g/h;
					f += e[j]*z(i,j);
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=z(i,j);
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						z(j,k) -= (f*e[k]+g*z(i,k));
				}
			}
		} else {
			e[i]=z(i,l);
		}
		d[i]=h;
	}
	d[0]=0.0;

	e[0]=0.0;
	for (i=0;i<_size;i++) {
		if (d[i] != 0.0) {
			for (j=0;j<i;j++) {
				g=0.0;
				for (k=0;k<i;k++)
					g += z(i,k)*z(k,j);
				for (k=0;k<i;k++)
					z(k,j) -= g*z(k,i);
			}
		}
		d[i]=z(i,i);
		z(i,i)=1.0;
		for (j=0;j<i;j++) z(j,i)=z(i,j)=0.0;
	}

	// QL iteration
	for (i=1;i<_size;i++) e[i-1]=e[i];
	e[_size-1]=0.0;
	for (l=0;l<_size;l++) {
		iter=0;
		do {
			for (m=l;m<_size-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if (abs(e[m]) <= std::numeric_limits<double>::epsilon()*dd) break;
			}
			if (m != l) {
				if (iter++ == 30) { 
					std::cerr << "Too many iterations in tqli" << std::endl;
					std::exit(-1);
				}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				absg=abs(g);
				r = ( (absg > 1.0) ? absg*sqrt(1.0+(1.0/absg)*(1.0/absg)) : sqrt(1.0+absg*absg));
				g=d[m]-d[l]+e[l]/(g+((g>=0.0) ? fabs(r):-fabs(r)));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					absf=abs(f);
					absg=abs(g);
					pfg = (absf > absg ? absf*sqrt(1.0+(absg/absf)*(absg/absf)) : (absg == 0.0 ? 0.0 : absg*sqrt(1.0+(absf/absg)*(absf/absg))));
					e[i+1]=(r=pfg); 
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<_size;k++) {
						f=z(k,i+1);
						z(k,i+1)=s*z(k,i)+c*f;
						z(k,i)=c*z(k,i)-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
	
	// Sort eigenvalues and reorder eigenvectors
	//for (int i=0;i<_size-1;i++) {
	//	double p=d[k=i];
	//	for (int j=i;j<_size;j++)
	//		if (d[j] >= p) p=d[k=j];
	//	if (k != i) {
	//		d[k]=d[i];
	//		d[i]=p;
	//		for (int j=0;j<_size;j++) {
	//			p=z(j,i);
	//			z(j,i)=z(j,k);
	//			z(j,k)=p;
	//		}
	//	}
	//}

	//Copy over data
	V = z;
	D.reset();
	for(size_t i=0;i<_size;++i) {
		D(i,i) = d[i];
	}
}
*/

template <size_t _size>
inline const SymmetricMatrix<_size>& operator~(const SymmetricMatrix<_size>& M) { 
	return M; 
}


// Matrix trace
template <size_t _size>
inline double tr(const Matrix<_size, _size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i){
		trace += q(i, i);
	}
	return trace;
}

template <size_t _size>
inline double tr(const SymmetricMatrix<_size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i){
		trace += q(i, i);
	}
	return trace;
}

// Matrix 1-norm
template <size_t _numRows, size_t _numColumns>
inline double norm(const Matrix<_numRows, _numColumns>& q) {
	double norm1 = double(0);
	for (size_t j = 0; j < _numColumns; ++j) {
		double colabssum = double(0);
		for (size_t i = 0; i < _numRows; ++i) {
			colabssum += abs(q(i,j));
		}
		if (colabssum > norm1) {
			norm1 = colabssum;
		}
	}
	return norm1;
}

// Identity matrix
template <size_t _size>
inline SymmetricMatrix<_size> identity() {
	SymmetricMatrix<_size> m;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			m(i,j) = (i == j ? 1.0 : 0.0);
		}
	}
	return m;
}

// Zero matrix
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> zeros() {
	Matrix<_numRows, _numColumns> m;
	m.reset();
	return m;
}

template <size_t _size>
inline SymmetricMatrix<_size> zeros() {
	SymmetricMatrix<_size> L;
	L.reset();
	return L;
}

template <size_t _numRows>
inline Matrix<_numRows> zero() {
	Matrix<_numRows> m;
	m.reset();
	return m;
}

// Matrix determinant
template <size_t _size>
inline double det(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	double D = double(1);

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = std::abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
			D = -D;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
			D = -D;
		}

		D *= m(row_p[k], col_p[k]);
		if (D == double(0)) {
			return double(0);
		}

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
		}  
	}

	return D;
}

template <size_t _size>
inline double det(const SymmetricMatrix<_size>& q) { 
	Matrix<_size, _size> m = q;
	double D = double(1);

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		D *= m(k, k);
		if (D == double(0)) {
			return double(0);
		}

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(i, k) / m(k, k);
			for (size_t j = k + 1; j < _size; ++j) {
				m(i, j) -= factor * m(k, j);
			}
		}  
	}

	return D;
}


// P%Q solves PX = Q for X
template <size_t _size, size_t _numColumns>
inline Matrix<_size, _numColumns> operator%(const Matrix<_size, _size>& p, const Matrix<_size, _numColumns>& q) {
	Matrix<_size, _size> m(p);
	Matrix<_size, _numColumns> inv(q);

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		assert(maximum != double(0));

		// swap rows and columns
		std::swap(row_p[k], row_p[max_row]);
		std::swap(col_p[k], col_p[max_col]);

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
			for (size_t j = 0; j < _numColumns; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		} 
	}

	// Backward substitution
	for (size_t k = _size - 1; k != -1; --k) {
		double quotient = m(row_p[k], col_p[k]);

		for (size_t j = 0; j < _numColumns; ++j) {
			inv(row_p[k], j) /= quotient;
		}

		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			for (size_t j = 0; j < _numColumns; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		}
	}

	// reshuffle result
	size_t invrow_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		invrow_p[row_p[i]] = i;
	}
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			std::swap(inv(col_p[i], j), inv(row_p[i], j));
		}
		row_p[invrow_p[col_p[i]]] = row_p[i];
		invrow_p[row_p[i]] = invrow_p[col_p[i]];
	}

	return inv;
}

template <size_t _size, size_t _numColumns>
inline Matrix<_size, _numColumns> operator%(const SymmetricMatrix<_size>& p, const Matrix<_size, _numColumns>& q) {
	// Cholesky factorization p = L*~L
	SymmetricMatrix<_size> L; // abuse SymmetricMatrix for triangular matrix
	L.reset();
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = i; j < _size; ++j) {
			double sum = p(j,i);
			for (size_t k = 0; k < i; ++k) {
				sum -= L(j,k)*L(i,k);
			}
			if (i == j) {
				assert(sum > 0.0);
				L(i,i) = sqrt(sum);
			} else {
				L(j,i) = sum / L(i,i);
			}
		}
	}

	// Backward and forward substitution
	Matrix<_size, _numColumns> M;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t k = 0; k < _numColumns; ++k) {
			double sum = q(i,k);
			for (size_t j = 0; j < i; ++j) {
				sum -= L(i,j)*M(j,k);
			}
			M(i,k) = sum / L(i,i);
		}
	}
	for (size_t i = _size - 1; i != -1; --i) {
		for (size_t k = 0; k < _numColumns; ++k) {
			double sum = M(i,k);
			for (size_t j = i + 1; j < _size; ++j) {
				sum -= L(j,i)*M(j,k);
			}
			M(i,k) = sum / L(i,i);
		}
	}
	return M;
}

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator/(const Matrix<_numRows, _size>& p, const Matrix<_size, _size>& q) {
	return ~(~q%~p);
}

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator/(const Matrix<_numRows, _size>& p, const SymmetricMatrix<_size>& q) {
	return ~(q%~p);
}

// Matrix pseudo-inverse
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numColumns, _numRows> pseudoInverse(const Matrix<_numRows, _numColumns>& q) { 
	if (_numColumns <= _numRows) {
		Matrix<_numColumns, _numColumns> Vec;
		SymmetricMatrix<_numColumns> Val;
		jacobi(SymProd(~q,q), Vec, Val);

		for (size_t i = 0; i < _numColumns; ++i) {
			if (abs(Val(i,i)) <= sqrt(DBL_EPSILON)) {
				Val(i,i) = 0.0;
			} else {
				Val(i,i) = 1.0 / Val(i,i);
			}
		}
		return SymProd(Vec,Val*~Vec)*~q;
	} else {
		Matrix<_numRows, _numRows> Vec;
		SymmetricMatrix<_numRows> Val;
		jacobi(SymProd(q,~q), Vec, Val);

		for (size_t i = 0; i < _numRows; ++i) {
			if (abs(Val(i,i)) <= sqrt(DBL_EPSILON)) {
				Val(i,i) = 0.0;
			} else {
				Val(i,i) = 1.0 / Val(i,i);
			}
		}
		return ~q*SymProd(Vec,Val*~Vec);
	}
}

// Matrix inverse
template <size_t _size>
inline Matrix<_size, _size> operator!(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	Matrix<_size, _size> inv = identity<_size>();

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
		swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;

		// eliminate column
		assert(maximum != double(0));
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);

			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}

			for (size_t j = 0; j < k; ++j) {
				inv(row_p[i], row_p[j]) -= factor * inv(row_p[k], row_p[j]);
			}
			inv(row_p[i], row_p[k]) = -factor;
		} 
	}

	// Backward substitution
	for (size_t k = _size - 1; k != -1; --k) {
		double quotient = m(row_p[k], col_p[k]);

		for (size_t j = 0; j < _size; ++j) {
			inv(row_p[k], j) /= quotient;
		}

		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			for (size_t j = 0; j < _size; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		}
	}

	// reshuffle result
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(col_p[i], j) = inv(row_p[i], j);
		}
	}

	return m; 
}

// Matrix exponentiation
#define _MATRIX_B0 1729728e1
#define _MATRIX_B1 864864e1
#define _MATRIX_B2 199584e1
#define _MATRIX_B3 2772e2
#define _MATRIX_B4 252e2
#define _MATRIX_B5 1512e0
#define _MATRIX_B6 56e0
#define _MATRIX_B7 1e0
#define _NORMLIM 9.504178996162932e-1

template <size_t _size>
inline Matrix<_size, _size> exp(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> A(q);
	int s = (int) std::max(double(0), ceil(log(norm(A)/_NORMLIM)*M_LOG2E)); 

	A /= pow(2.0,s);
	Matrix<_size, _size> A2(A*A);
	Matrix<_size, _size> A4(A2*A2);
	Matrix<_size, _size> A6(A2*A4);
	Matrix<_size, _size> U( A*(A6*_MATRIX_B7 + A4*_MATRIX_B5 + A2*_MATRIX_B3 + identity<_size>()*_MATRIX_B1) );
	Matrix<_size, _size> V( A6*_MATRIX_B6 + A4*_MATRIX_B4 + A2*_MATRIX_B2 + identity<_size>()*_MATRIX_B0 ); 
	Matrix<_size, _size> R7 = (V - U)%(V + U);

	for (int i = 0; i < s; ++i) {
		R7 *= R7;
	}
	return R7;
}

// Input stream
template <size_t _numRows, size_t _numColumns>
inline std::istream& operator>>(std::istream& is, Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			is >> q(i,j);
		}
	}
	return is;
}

// Output stream
template <size_t _numRows, size_t _numColumns>
inline std::ostream& operator<<(std::ostream& os, const Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			os  << std::left << std::setw(13) << q(i,j) << " ";
		}
		os << std::endl;
	}
	return os;
}

template <size_t _size>
inline std::ostream& operator<<(std::ostream& os, const SymmetricMatrix<_size>& q) {
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			os << std::left << std::setw(13) << q(i,j) << " ";
		}
		os << std::endl;
	}
	return os;
}


// Hilbert matrix
template <size_t _size>
inline Matrix<_size, _size> hilbert() {
	Matrix<_size, _size> m;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(i, j) = double(1) / (double) (i + j + 1);
		}
	}
	return m;
}

// Given a positive-definite symmetric matrix q, construct its Cholesky decomposition, q = z.z'
// From: Numerical recipes in C++ (3rd edition)
template <size_t _size>
inline void chol(const SymmetricMatrix<_size>& q, Matrix<_size, _size>& z)
{
	z = q;
	int n = _size;
	int i,j,k;
	double sum;
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			for (sum = z(i,j),k = i-1; k >= 0; k--) { 
				sum -= z(i,k)*z(j,k);
			}
			if (i == j) {
				// Matrix, with rounding errors, is not positive-definite.
				if (sum <= 0.0) {
					std::cerr << "Cholesky failed, matrix not positive definite" << std::endl;
					std::exit(-1);
				}
				z(i,i) = sqrt(sum);
			} else {
				z(j,i) = sum/z(i,i);
			}
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) { 
			z(j,i) = 0.;
		}
	}
}

// Eigen decomposition of a real symmetric matrix by converting to tridiagonal form followed by QL iteration
// From: Numerical recipes in C++ (3rd edition)
template <size_t _size>
inline void jacobi(const SymmetricMatrix<_size>& q, Matrix<_size, _size>& z, SymmetricMatrix<_size>& D) 
{
	// Initializations
	z = q;
	Matrix<_size> d, e;
	d.reset();
	e.reset();

	int l,k,j,i,m,iter;
	double scale,hh,h,g,f,s,r,p,dd,c,b,absf,absg,pfg;

	//Convert to tridiagonal form
	for (i=_size-1;i>0;i--) 
	{
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++)
				scale += abs(z(i,k));
			if (scale == 0.0)
				e[i]=z(i,l);
			else {
				for (k=0;k<i;k++) {
					z(i,k) /= scale;
					h += z(i,k)*z(i,k);
				}
				f=z(i,l);
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z(i,l)=f-g;
				f=0.0;
				for (j=0;j<i;j++) {
					z(j,i)=z(i,j)/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += z(j,k)*z(i,k);
					for (k=j+1;k<i;k++)
						g += z(k,j)*z(i,k);
					e[j]=g/h;
					f += e[j]*z(i,j);
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=z(i,j);
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						z(j,k) -= (f*e[k]+g*z(i,k));
				}
			}
		} else {
			e[i]=z(i,l);
		}
		d[i]=h;
	}
	d[0]=0.0;

	e[0]=0.0;
	for (i=0;i<_size;i++) {
		if (d[i] != 0.0) {
			for (j=0;j<i;j++) {
				g=0.0;
				for (k=0;k<i;k++)
					g += z(i,k)*z(k,j);
				for (k=0;k<i;k++)
					z(k,j) -= g*z(k,i);
			}
		}
		d[i]=z(i,i);
		z(i,i)=1.0;
		for (j=0;j<i;j++) z(j,i)=z(i,j)=0.0;
	}

	// QL iteration
	for (i=1;i<_size;i++) e[i-1]=e[i];
	e[_size-1]=0.0;
	for (l=0;l<_size;l++) {
		iter=0;
		do {
			for (m=l;m<_size-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if (abs(e[m]) <= std::numeric_limits<double>::epsilon()*dd) break;
			}
			if (m != l) {
				if (iter++ == 30) { 
					std::cerr << "Too many iterations in tqli" << std::endl;
					std::exit(-1);
				}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				absg=abs(g);
				r = ( (absg > 1.0) ? absg*sqrt(1.0+(1.0/absg)*(1.0/absg)) : sqrt(1.0+absg*absg));
				g=d[m]-d[l]+e[l]/(g+((g>=0.0) ? fabs(r):-fabs(r)));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					absf=abs(f);
					absg=abs(g);
					pfg = (absf > absg ? absf*sqrt(1.0+(absg/absf)*(absg/absf)) : (absg == 0.0 ? 0.0 : absg*sqrt(1.0+(absf/absg)*(absf/absg))));
					e[i+1]=(r=pfg); 
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<_size;k++) {
						f=z(k,i+1);
						z(k,i+1)=s*z(k,i)+c*f;
						z(k,i)=c*z(k,i)-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}

	//Copy over data
	D.reset();
	for(size_t i=0;i<_size;++i) {
		D(i,i) = d[i];
	}
}

// Matrix addition
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator+(const Matrix<_numRows, _numColumns>& M, const Matrix<_numRows, _numColumns>& N) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] + N[i];
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator+(const SymmetricMatrix<_size>& M, const Matrix<_size, _size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) + N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator+(const Matrix<_size, _size>& M, const SymmetricMatrix<_size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) + N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator+(const SymmetricMatrix<_size>& M, const SymmetricMatrix<_size>& N) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] + N[i];
	}
	return L;
}

// Matrix subtraction
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator-(const Matrix<_numRows, _numColumns>& M, const Matrix<_numRows, _numColumns>& N) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] - N[i];
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator-(const SymmetricMatrix<_size>& M, const Matrix<_size, _size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) - N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator-(const Matrix<_size, _size>& M, const SymmetricMatrix<_size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) - N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator-(const SymmetricMatrix<_size>& M, const SymmetricMatrix<_size>& N) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] - N[i];
	}
	return L;
}

// Scalar multiplication
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator*(const Matrix<_numRows, _numColumns>& M, double a) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] * a;
	}
	return L;
}

template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator*(double a, const Matrix<_numRows, _numColumns>& M) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = a * M[i];
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator*(const SymmetricMatrix<_size>& M, double a) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] * a;
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator*(double a, const SymmetricMatrix<_size>& M) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = a * M[i];
	}
	return L;
}

// Scalar division
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator/(const Matrix<_numRows, _numColumns>& M, double a) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] / a;
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator/(const SymmetricMatrix<_size>& M, double a) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] / a;
	}
	return L;
}

// Matrix multiplication
template <size_t _numRows, size_t _numColumns, size_t __numColumns>
inline Matrix<_numRows, __numColumns> operator*(const Matrix<_numRows, _numColumns>& M, const Matrix<_numColumns, __numColumns>& N) {
	Matrix<_numRows, __numColumns> L;
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < __numColumns; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _numColumns; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

template <size_t _size, size_t _numColumns>
inline Matrix<_size, _numColumns> operator*(const SymmetricMatrix<_size>& M, const Matrix<_size, _numColumns>& N) {
	Matrix<_size, _numColumns> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator*(const Matrix<_numRows, _size>& M, const SymmetricMatrix<_size>& N) {
	Matrix<_numRows, _size> L;
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator*(const SymmetricMatrix<_size>& M, const SymmetricMatrix<_size>& N) {
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

// Compute product M*N of which one knows that the results is symmetric (and save half the computation)
template <size_t _size, size_t _numRows>
inline SymmetricMatrix<_size> SymProd(const Matrix<_size, _numRows>& M, const Matrix<_numRows, _size>& N) {
	SymmetricMatrix<_size> S;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			double temp = double(0);
			for (size_t k = 0; k < _numRows; ++k) {
				temp += M(i, k) * N(k, j);
			}
			S(i, j) = temp;
		}
	}
	return S;
}

// Compute sum M+M^T
template <size_t _size>
inline SymmetricMatrix<_size> SymSum(const Matrix<_size, _size>& M) {
	SymmetricMatrix<_size> S;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			S(i,j) = M(i,j) + M(j,i);
		}
	}
	return S;
}

inline double scalar(const Matrix<1,1>& M) {
	return M[0];
}


// Principal square root

template <size_t _size>
inline void jacobi(const Matrix<_size, _size>& q, Matrix<_size, _size>& V, Matrix<_size, _size>& D) {
	D = q;
	V = identity<_size>();

	while (true) {
		double maximum = 0; size_t max_row = 0; size_t max_col = 0;
		for (size_t i = 0; i < _size; ++i) {
			for (size_t j = i + 1; j < _size; ++j) {
				if (abs(D(i,j)) > maximum) {
					maximum = abs(D(i,j));
					max_row = i;
					max_col = j;
				}
			}
		}

		if (maximum <= DBL_EPSILON) {
			break;
		}

		double theta = (D(max_col, max_col) - D(max_row, max_row)) / (2 * D(max_row, max_col));
		double t = 1 / (abs(theta) + sqrt(theta*theta+1));
		if (theta < 0) t = -t;
		double c = 1 / sqrt(t*t+1); 
		double s = c*t;

		Matrix<_size, _size> R = identity<_size>();
		R(max_row,max_row) = c;
		R(max_col,max_col) = c;
		R(max_row,max_col) = s;
		R(max_col,max_row) = -s;

		// update D // 
		//std::cout << D << std::endl;
		//D = ~R * D * R;

		double temp1 = c*c*D(max_row, max_row) + s*s*D(max_col, max_col) - 2*c*s*D(max_row, max_col);
		double temp2 = s*s*D(max_row, max_row) + c*c*D(max_col, max_col) + 2*c*s*D(max_row, max_col);
		D(max_row, max_col) = 0;
		D(max_col, max_row) = 0;
		D(max_row, max_row) = temp1;
		D(max_col, max_col) = temp2;
		for (int j = 0; j < _size; ++j) {
			if ((j != max_row) && (j != max_col)) {
				temp1 = c * D(j, max_row) - s * D(j, max_col);
				temp2 = c * D(j, max_col) + s * D(j, max_row);
				D(j, max_row) = (D(max_row, j) = temp1);
				D(j, max_col) = (D(max_col, j) = temp2);
			}
		}
		//std::cout << D << std::endl << std::endl;


		V = V * R;
	} 
}

template <size_t _size>
inline Matrix<_size, _size> sqrtm(const Matrix<_size, _size>& X) {
  Matrix<_size,_size> V;
  Matrix<_size,_size> D;
  jacobi(X, V, D);
  for (size_t i = 0; i < _size; ++i) {
		if (D(i,i) > 0) {
			D(i,i) = sqrt(D(i,i));
		} else {
			D(i,i) = 0;
		}
	}
	return (V*D*~V);
}
#endif