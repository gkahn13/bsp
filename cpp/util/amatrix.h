#ifndef __MATRIX_H__
#define __MATRIX_H__

#define _USE_MATH_DEFINES



#include <cmath>
#include <float.h>
#include <iostream>
#include <assert.h>



template <size_t _numRows, size_t _numColumns = 1> class Matrix {

private:
  adouble _elems[_numRows * _numColumns];
  //adouble *_elems;

public:
  /*Matrix() { _elems = new adouble[_numRows * _numColumns]; }
  Matrix( const Matrix<_numRows, _numColumns>& q ) {
    _elems = new adouble[_numRows * _numColumns];
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      _elems[i] = q._elems[i]; 
    }
  }
  
  ~Matrix() { delete[] _elems; }
  
  Matrix<_numRows, _numColumns>& operator = (const Matrix<_numRows, _numColumns>& q) {
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      _elems[i] = q._elems[i]; 
    }
    return (*this);
  }*/

  // Retrieval
  inline size_t numRows() const { 
    return _numRows; 
  }
  inline size_t numColumns() const { 
    return _numColumns; 
  }
  
  // Subscript operator
  inline adouble& operator () (size_t row, size_t column) {
    assert(row < _numRows && column < _numColumns);
    return _elems[row * _numColumns + column]; 
  }
  inline adouble  operator () (size_t row, size_t column) const {
    assert(row < _numRows && column < _numColumns);
    return _elems[row * _numColumns + column]; 
  }

  inline adouble& operator [] (size_t elt) {
    assert(elt < _numRows * _numColumns);
    return _elems[elt]; 
  }
  inline adouble  operator [] (size_t elt) const {
    assert(elt < _numRows * _numColumns);
    return _elems[elt]; 
  }

  // Reset to zeros
  inline void reset() { 
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      _elems[i] = adouble(0);
    }
  }

  // Submatrix
  template <size_t nRows, size_t nCols>
  inline Matrix<nRows, nCols> subMatrix(size_t row, size_t column) const {
    assert(row + nRows <= _numRows && column + nCols <= _numColumns);
    Matrix<nRows, nCols> m;
    for (size_t i = 0; i < nRows; ++i) {
      for (size_t j = 0; j < nCols; ++j) {
        m(i, j) = (*this)(row + i, column + j);
      }
    }
    return m;
  }

  template <size_t nRows>
  inline Matrix<nRows> subMatrix(size_t row, size_t column) const {
    assert(row + nRows <= _numRows && column <= _numColumns);
    Matrix<nRows> m;
    for (size_t i = 0; i < nRows; ++i) {
      m[i] = (*this)(row + i, column);
    }
    return m;
  }

  // Insert
  template <size_t nRows, size_t nCols>
  inline void insert(size_t row, size_t column, const Matrix<nRows, nCols>& q) {
    assert(row + nRows <= _numRows && column + nCols <= _numColumns);
    for (size_t i = 0; i < nRows; ++i) {
      for (size_t j = 0; j < nCols; ++j) {
        (*this)(row + i, column + j) = q(i, j);
      }
    }
  }
  
  // Unary minus
  inline Matrix<_numRows, _numColumns> operator-() const {
    Matrix<_numRows, _numColumns> m;
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      m._elems[i] = -_elems[i];
    }
    return m;
  }
  
  // Unary plus
  inline const Matrix<_numRows, _numColumns>& operator+() const {
    return *this; 
  }

  // Equality
  inline bool operator==(const Matrix<_numRows, _numColumns>& q) const {
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      if (_elems[i] < q._elems[i] || _elems[i] > q._elems[i]) {
        return false;
      }
    }
    return true;
  }

  // Inequality
  inline bool operator!=(const Matrix<_numRows, _numColumns>& q) const {
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
        if (_elems[i] < q._elems[i] || _elems[i] > q._elems[i]) {
        return true;
      }
    }
    return false;
  }

  // Matrix addition
  inline Matrix<_numRows, _numColumns> operator+(const Matrix<_numRows, _numColumns>& q) const {
    Matrix<_numRows, _numColumns> m;
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      m._elems[i] = _elems[i] + q._elems[i];
    }
    return m;
  }
  inline const Matrix<_numRows, _numColumns>& operator+=(const Matrix<_numRows, _numColumns>& q) {
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      _elems[i] += q._elems[i];
    }
    return *this; 
  }
  
  // Matrix subtraction
  inline Matrix<_numRows, _numColumns> operator-(const Matrix<_numRows, _numColumns>& q) const {
    Matrix<_numRows, _numColumns> m;
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      m._elems[i] = _elems[i] - q._elems[i];
    }
    return m;
  }
  inline const Matrix<_numRows, _numColumns>& operator-=(const Matrix<_numRows, _numColumns>& q) {
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      _elems[i] -= q._elems[i];
    }
    return *this; 
  }

  // Scalar multiplication
  inline Matrix<_numRows, _numColumns> operator*(adouble a) const {
    Matrix<_numRows, _numColumns> m;
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      m._elems[i] = _elems[i] * a;
    }
    return m;
  }
  inline const Matrix<_numRows, _numColumns>& operator*=(adouble a) {
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      _elems[i] *= a;
    }
    return *this;
  }
  
  // Scalar division
  inline Matrix<_numRows, _numColumns> operator/(adouble a) const {
    Matrix<_numRows, _numColumns> m;
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      m._elems[i] = _elems[i] / a;
    }
    return m;
  }
  inline const Matrix<_numRows, _numColumns>& operator/=(adouble a) {
    for (size_t i = 0; i < _numRows * _numColumns; ++i) {
      _elems[i] /= a;
    }
    return *this;
  }

  // Matrix multiplication
  template <size_t nCols>
  inline Matrix<_numRows, nCols> operator*(const Matrix<_numColumns, nCols>& q) const {
    Matrix<_numRows, nCols> m;
    for (size_t i = 0; i < _numRows; ++i) {
      for (size_t j = 0; j < nCols; ++j) {
        adouble temp = adouble(0);
        for (size_t k = 0; k < _numColumns; ++k) {
          temp += (*this)(i, k) * q(k, j);
        }
        m(i, j) = temp;
      }
    }
    return m;
  }
  
  inline const Matrix<_numRows, _numColumns>& operator*=(const Matrix<_numColumns, _numColumns>& q) {
    return ((*this) = (*this) * q); 
  }

  // Matrix transpose
  inline Matrix<_numColumns, _numRows> operator~() const {
    Matrix<_numColumns, _numRows> m;
    for (size_t i = 0; i < _numColumns; ++i) {
      for (size_t j = 0; j < _numRows; ++j) {
        m(i, j) = (*this)(j, i);
      }
    }
    return m;
  }
};

/*template < >
class Matrix<1,1>: public Matrix
{
  // Casting to adouble for 1x1 matrix
  inline operator adouble() const {
    return _elems[0];
  }
};*/


// Scalar multiplication 
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator*(adouble a, const Matrix<_numRows, _numColumns>& q) { return q*a; }

// Matrix trace
template <size_t _size>
inline adouble tr(const Matrix<_size, _size>& q) {
  adouble trace = adouble(0);
  for (size_t i = 0; i < _size; ++i){
    trace += q(i, i);
  }
  return trace;
}

// Matrix 1-norm
template <size_t _numRows, size_t _numColumns>
inline adouble norm(const Matrix<_numRows, _numColumns>& q) {
  adouble norm1 = adouble(0);
  for (size_t j = 0; j < _numColumns; ++j) {
    adouble colabssum = adouble(0);
    for (size_t i = 0; i < _numRows; ++i) {
      colabssum += std::abs(q(i,j));
    }
    if (colabssum > norm1) {
      norm1 = colabssum;
    }
  }
  return norm1;
}


// Identity matrix
template <size_t _size>
inline Matrix<_size, _size> identity() {
  Matrix<_size, _size> m;
  for (size_t i = 0; i < _size; ++i) {
    for (size_t j = 0; j < _size; ++j) {
      m(i, j) = (i == j ? adouble(1) : adouble(0));
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

template <size_t _numRows>
inline Matrix<_numRows> zeros() {
  Matrix<_numRows> m;
  m.reset();
  return m;
}

// Matrix determinant
template <size_t _size>
inline adouble det(const Matrix<_size, _size>& q) {
  Matrix<_size, _size> m(q);
  adouble D = adouble(1);

  size_t row_p[_size];
  size_t col_p[_size];
  for (size_t i = 0; i < _size; ++i) {
    row_p[i] = i; col_p[i] = i;
  }

  // Gaussian elimination
  for (size_t k = 0; k < _size; ++k) {
    // find maximal pivot element
    adouble maximum = adouble(0); size_t max_row = k; size_t max_col = k;
    for (size_t i = k; i < _size; ++i) {
      for (size_t j = k; j < _size; ++j) {
        adouble abs_ij = std::abs(m(row_p[i], col_p[j]));
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
    if (D == adouble(0)) {
      return adouble(0);
    }

    // eliminate column
    for (size_t i = k + 1; i < _size; ++i) {
      adouble factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
      for (size_t j = k + 1; j < _size; ++j) {
        m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
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
    adouble maximum = adouble(0); size_t max_row = k; size_t max_col = k;
    for (size_t i = k; i < _size; ++i) {
      for (size_t j = k; j < _size; ++j) {
        adouble abs_ij = std::abs(m(row_p[i], col_p[j]));
        if (abs_ij > maximum) {
          maximum = abs_ij; max_row = i; max_col = j;
        }
      }
    }

    assert(maximum != adouble(0));

    // swap rows and columns
    std::swap(row_p[k], row_p[max_row]);
    std::swap(col_p[k], col_p[max_col]);
    
    // eliminate column
    for (size_t i = k + 1; i < _size; ++i) {
      adouble factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
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
    adouble quotient = m(row_p[k], col_p[k]);

    for (size_t j = 0; j < _numColumns; ++j) {
      inv(row_p[k], j) /= quotient;
    }
    
    for (size_t i = 0; i < k; ++i) {
      adouble factor = m(row_p[i], col_p[k]);
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

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator/(const Matrix<_numRows, _size>& p, const Matrix<_size, _size>& q) {
  return ~(~q%~p);
}

// Matrix pseudo-inverse
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numColumns, _numRows> pseudoInverse(const Matrix<_numRows, _numColumns>& q) {
  if (_numColumns <= _numRows) {
    Matrix<_numColumns, _numColumns> Vec, Val;
    jacobi(~q*q, Vec, Val);

    for (size_t i = 0; i < _numColumns; ++i) {
      if (std::abs(Val(i,i)) <= std::sqrt(DBL_EPSILON)) {
      //if (std::abs(Val(i,i)) <= DBL_EPSILON) {
        Val(i,i) = 0.0;
      } else {
        Val(i,i) = 1.0 / Val(i,i);
      }
    }
    return (Vec*Val*~Vec)*~q;
  } else {
    Matrix<_numRows, _numRows> Vec, Val;
    jacobi(q*~q, Vec, Val);

    for (size_t i = 0; i < _numRows; ++i) {
      if (std::abs(Val(i,i)) <= std::sqrt(DBL_EPSILON)) {
        Val(i,i) = 0.0;
      } else {
        Val(i,i) = 1.0 / Val(i,i);
      }
    }
    return ~q*(Vec*Val*~Vec);
  }
}

/*template <size_t _numRows, size_t _numColumns>
inline Matrix<_numColumns, _numRows> pseudoInverse(const Matrix<_numRows, _numColumns>& q) {
  Matrix<_numRows, _numColumns> m(q);
  
  size_t row_p[_numRows];
  size_t col_p[_numColumns];
  for (size_t i = 0; i < _numRows; ++i) {
    row_p[i] = i;
  }
  for (size_t i = 0; i < _numColumns; ++i) {
    col_p[i] = i;
  }

  // Gaussian elimination to determine rank
  size_t rank = (_numRows < _numColumns ? _numRows : _numColumns);
  for (size_t k = 0; k < rank; ++k) {
    // find maximal pivot element
    adouble maximum = adouble(0); size_t max_row = k; size_t max_col = k;
    for (size_t i = k; i < _numRows; ++i) {
      for (size_t j = k; j < _numColumns; ++j) {
        adouble abs_ij = std::abs(m(row_p[i], col_p[j]));
        if (abs_ij > maximum) {
          maximum = abs_ij; max_row = i; max_col = j;
        }
      }
    }

    // check if zero
    if (maximum <= DBL_EPSILON) {
      rank = k;
      break;
    }

    // swap rows and columns
    if (k != max_row) {
      size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
    }
    if (k != max_col) {
      size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
    }

    // eliminate column
    for (size_t i = k + 1; i < _numRows; ++i) {
      adouble factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
      for (size_t j = k + 1; j < _numColumns; ++j) {
        m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
      }
    } 
  }

  if (rank == 0) { // zero or empty matrix
    return zeros<_numColumns, _numRows>();
  } else if (rank == _numRows && rank == _numColumns) { // full rank square matrix
    // convert so that compiler accepts number of rows and columns for inverse
    Matrix<_numRows, _numRows> q2;
    for (size_t i = 0; i < _numRows * _numRows; ++i) {
      q2[i] = q[i];
    }
    q2 = !q2;
    Matrix<_numColumns, _numRows> q3;
    for (size_t i = 0; i < _numRows * _numRows; ++i) {
      q3[i] = q2[i];
    }
    return q3;
  } else if (rank == _numRows) { // full row-rank matrix
    return ~q/(q*~q);
  } else if (rank == _numColumns) { // full column-rank matrix
    return (~q*q)%~q;
  } // else: not full rank, perform rank decomposition

  // bring m into reduced row echelon form
  
  // normalize rows such that 1's appear on the diagonal
  for (size_t k = 0; k < rank; ++k) {
    adouble quotient = m(row_p[k], col_p[k]);
    for (size_t j = 0; j < k; ++j) {
      m(row_p[k], col_p[j]) = adouble(0);
    }
    m(row_p[k], col_p[k]) = adouble(1);
    for (size_t j = k + 1; j < _numColumns; ++j) {
      m(row_p[k], col_p[j]) /= quotient;
    }
  }
  
  // subtract rows such that 0's appear above leading row elements
  for (size_t k = rank - 1; k != -1; --k) {
    for (size_t i = 0; i < k; ++i) {
      adouble factor = m(row_p[i], col_p[k]);
      m(row_p[i], col_p[k]) = adouble(0);
      for (size_t j = k + 1; j < _numColumns; ++j) {
        m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
      }
    } 
  }

  // copy m into smaller matrix C and swap columns
  Matrix<(_numRows < _numColumns ? _numRows : _numColumns), _numColumns> C;
  for (size_t k = 0; k < rank; ++k) {
    for (size_t j = 0; j < _numColumns; ++j) {
      C(k,j) = m(row_p[k], j);
    }
  }
  for (size_t k = rank; k < C.numRows(); ++k) {
    for (size_t j = 0; j < _numColumns; ++j) {
      C(k,j) = adouble(0);
    }
  }

  // B is copy of A with columns swapped and non-pivot columns left out
  Matrix<_numRows, (_numRows < _numColumns ? _numRows : _numColumns)> B;
  for (size_t k = 0; k < _numRows; ++k) {
    for (size_t j = 0; j < rank; ++j ) {
      B(k,j) = q(k, col_p[j]);
    }
    for (size_t j = rank; j < B.numColumns(); ++j) {
      B(k,j) = adouble(0);
    }
  }

  // Now, q = B*C and B and C are of full row/column rank
  return ~C*!((~B*B)*(C*~C))*~B;
}*/

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
    adouble maximum = adouble(0); size_t max_row = k; size_t max_col = k;
    for (size_t i = k; i < _size; ++i) {
      for (size_t j = k; j < _size; ++j) {
        adouble abs_ij = std::abs(m(row_p[i], col_p[j]));
        if (abs_ij > maximum) {
          maximum = abs_ij; max_row = i; max_col = j;
        }
      }
    }

    // swap rows and columns
    size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
           swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;

    // eliminate column
    assert(maximum != adouble(0));
    for (size_t i = k + 1; i < _size; ++i) {
      adouble factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
      
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
    adouble quotient = m(row_p[k], col_p[k]);
    
    for (size_t j = 0; j < _size; ++j) {
      inv(row_p[k], j) /= quotient;
    }
    
    for (size_t i = 0; i < k; ++i) {
      adouble factor = m(row_p[i], col_p[k]);
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

template <size_t _size>
inline void jacobi(const Matrix<_size, _size>& m, Matrix<_size, _size>& V, Matrix<_size, _size>& D) {
  D = m;
  V = identity<_size>();

  if (_size <= 1) {
    return;
  }

  size_t pivotRow = 0;
  size_t zeroCount = 0;

  while (true) {
    adouble maximum = 0;
    size_t p, q;
    for (size_t i = 0; i < pivotRow; ++i) {
      if (std::abs(D(i,pivotRow)) > maximum) {
        maximum = std::abs(D(i,pivotRow));
        p = i; q = pivotRow;
      }
    }
    for (size_t j = pivotRow + 1; j < _size; ++j) {
      if (std::abs(D(pivotRow,j)) > maximum) {
        maximum = std::abs(D(pivotRow,j));
        p = pivotRow; q = j;
      }
    }
    
    pivotRow = (pivotRow + 1) % _size;

    if (maximum <= DBL_EPSILON) {
      ++zeroCount;
      if (zeroCount == _size) {
        break;
      } else {
        continue;
      }
    } else {
      zeroCount = 0;
    }

    adouble theta = 0.5*(D(q,q) - D(p,p)) / D(p,q);
    adouble t = 1 / (fabs(theta) + sqrt(theta*theta+1));
    if (theta < 0) t = -t;
    adouble c = 1 / sqrt(t*t+1);
    adouble s = c*t;
    adouble tau = s / (1 + c);

    // update D // 
    
    // update D(r,p) and D(r,q)
    for (size_t r = 0; r < p; ++r) {
      adouble Drp = D(r,p); adouble Drq = D(r,q);
      D(r,p) -= s*(Drq + tau*Drp);
      D(r,q) += s*(Drp - tau*Drq);
    }
    for (size_t r = p + 1; r < q; ++r) {
      adouble Drp = D(p,r); adouble Drq = D(r,q);
      D(p,r) -= s*(Drq + tau*Drp);
      D(r,q) += s*(Drp - tau*Drq);
    }
    for (size_t r = q + 1; r < _size; ++r) {
      adouble Drp = D(p,r); adouble Drq = D(q,r);
      D(p,r) -= s*(Drq + tau*Drp);
      D(q,r) += s*(Drp - tau*Drq);
    }

    // update D(p,p), D(q,q), and D(p,q);
    D(p,p) -= t*D(p,q);
    D(q,q) += t*D(p,q);
    D(p,q) = 0;

    // Update V
    for (size_t r = 0; r < _size; ++r) {
      adouble Vrp = V(r,p); adouble Vrq = V(r,q);
      V(r,p) -= s*(Vrq + tau*Vrp);
      V(r,q) += s*(Vrp - tau*Vrq);
    }
  }

  // clean up D
  for (size_t i = 0; i < _size - 1; ++i) {
    for (size_t j = i+1; j < _size; ++j) {
      D(j,i) = D(i,j) = 0;
    }
  }
}


// Matrix exponentiation
#define _B0 1729728e1
#define _B1 864864e1
#define _B2 199584e1
#define _B3 2772e2
#define _B4 252e2
#define _B5 1512e0
#define _B6 56e0
#define _B7 1e0
#define _NORMLIM 9.504178996162932e-1

template <size_t _size>
inline Matrix<_size, _size> exp(const Matrix<_size, _size>& q) {
  Matrix<_size, _size> A(q);
  int s = (int) std::max(adouble(0), ceil(log(norm(A)/_NORMLIM)*M_LOG2E));

  A /= pow(2.0,s);
  Matrix<_size, _size> A2(A*A);
  Matrix<_size, _size> A4(A2*A2);
  Matrix<_size, _size> A6(A2*A4);
  Matrix<_size, _size> U( A*(A6*_B7 + A4*_B5 + A2*_B3 + identity<_size>()*_B1) );
  Matrix<_size, _size> V( A6*_B6 + A4*_B4 + A2*_B2 + identity<_size>()*_B0 );
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
      os << q(i,j) << "\t";
      //printf_s("%24.24g ", q(i,j));
    }
    //printf_s("\n");
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
      m(i, j) = adouble(1) / (adouble) (i + j + 1);
    }
  }
  return m;
}

// Solution to Continuous-time Algebraic Ricatti Equation ~AX + XA - XGX + Q = 0
template <size_t _size>
inline Matrix<_size, _size> care(const Matrix<_size, _size>& A, const Matrix<_size, _size>& G, const Matrix<_size, _size>& Q) {
  Matrix<2*_size, 2*_size> H, H2, W, SqrtH2;
  H.insert(0,0, A);      H.insert(0, _size, -G);
  H.insert(_size,0, -Q); H.insert(_size, _size, -~A);

  H2 = H*H;

  SqrtH2 = identity<2*_size>();
  for (size_t i = 0; i < 100; ++i) {
    SqrtH2 = 0.5*(SqrtH2 + !SqrtH2*H2);
  }

  W = H - SqrtH2;

  return W.subMatrix<_size, _size>(_size, 0) * !W.subMatrix<_size, _size>(0, 0);
}


// Principal square root
template <size_t _size>
inline Matrix<_size, _size> sqrt(const Matrix<_size, _size>& X) {
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
