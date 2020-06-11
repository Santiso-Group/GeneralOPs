/*
** Copyright 2007-2011 Erik Santiso.
** This file is part of mymath.
** mymath is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** mymath is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with mymath. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** Simple N-dimensional (square) matrix class
**
** Note that elements are stored by row.
*/

/*
** Things that may be good to add:
**
** - Direct compatibility with Vector3D, Vector4D, Matrix3D, Matrix4D and
**   Quaternion (e.g. be able to do MatrixND*Vector3D if N = 3, etc.).
**   Right now it's only possible to copy by constructing the ND-version
**   from the 3D/4D-version (slow).
**
** - Functions (probably in a separate file) to calculate determinant,
**   inverse, LU decomposition, etc. as needed
*/

#ifndef H_MATRIXND
#define H_MATRIXND

#include <iostream>
#include <string>
#include <sstream>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/matrix3D.h"
#include "mymath/include/matrix4D.h"
#include "mymath/include/vectorND.h"

class MatrixND
{
public:

// Constructors

  MatrixND();                                     // Defines an empty matrix
  explicit MatrixND(size_t const numDimensions);  // Defines a null matrix with given
                                                  // number of dimensions
  MatrixND(Matrix3D const &a);                    // Defines a MatrixND from a Matrix3D
  MatrixND(Matrix4D const &a);                    // Defines a MatrixND from a Matrix4D
  MatrixND(VectorND const &v);                    // Defines a diagonal matrix given a vector
                                                  // of diagonal elements
  MatrixND(size_t const numDimensions,            // Defines a matrix with given number
           Real const &r);                        // of dimensions and diagonal elements
                                                  // equal to a given number

// Assignment

  MatrixND &operator=(Real const &r);             // Makes all diagonal elements equal to r
  MatrixND &operator+=(MatrixND const &a);        // Assignment by addition
  MatrixND &operator-=(MatrixND const &a);        // Assignment by subtraction
  MatrixND &operator*=(Real const &r);            // Assignment by multiplication by scalar
  MatrixND &operator/=(Real const &r);            // Assignment by division by scalar

// Access elements/rows by index

  VectorND const &operator[](size_t const i) const;             // Returns the i-th row (const)
  Real &operator()(size_t const i, size_t const j);             // Can modify
  Real const &operator()(size_t const i, size_t const j) const; // Cannot modify

// Access by rows or columns

  VectorND const row(size_t const i) const;               // Returns the i-th row as a vector
  VectorND const column(size_t const i) const;            // Returns the i-th column as a vector
  void setRow(size_t const i, VectorND const &row);       // Replaces the i-th row with the given vector
  void setColumn(size_t const i, VectorND const &column); // Replaces the i-th column with the given vector

// Clear and resize

  void clear();                             // Clears the matrix
  void resize(size_t const numDimensions);  // Set the number of dimensions - extra elements are set to zero

// Size

  size_t const size() const; // Returns the number of dimensions

// Operators

  friend MatrixND operator+(MatrixND const &a1, MatrixND const &a2);      // Addition
  friend MatrixND operator-(MatrixND const &a);                           // Negation
  friend MatrixND operator-(MatrixND const &a1, MatrixND const &a2);      // Subtraction
  friend MatrixND operator*(Real const &r, MatrixND const &a);            // Multiplication by scalar
  friend MatrixND operator*(MatrixND const &a, Real const &r);            // Multiplication by scalar
  friend MatrixND operator/(MatrixND const &a, Real const &r);            // Division by scalar
  friend MatrixND operator*(MatrixND const &a1, MatrixND const &a2);      // Matrix product
  friend VectorND operator*(MatrixND const &a, VectorND const &v);        // Multiplication by vector
  friend VectorND operator*(VectorND const &v, MatrixND const &a);        // Multiplication by vector
  friend bool operator==(MatrixND const &a1, MatrixND const &a2);         // Comparison (equal)
  friend bool operator!=(MatrixND const &a1, MatrixND const &a2);         // Comparison (not equal)
  friend bool operator==(MatrixND const &a, Real const &r);               // Comparison with scalar
  friend bool operator!=(MatrixND const &a, Real const &r);               // Comparison with scalar
  friend std::ostream& operator<<(std::ostream &str, MatrixND const &a);  // Output
  friend std::istream& operator>>(std::istream &str, MatrixND &a);        // Input

// Other operations

  friend MatrixND const outer(VectorND const &v1, VectorND const &v2);  // Outer (Kronecker, dyadic) product
  friend Real const trace(MatrixND const &a);                           // Returns the trace of the matrix
  friend MatrixND const transpose(MatrixND const &a);                   // Returns the transpose of the matrix
  friend VectorND const diagonal(MatrixND const &a);                    // Returns the diagonal of the matrix
  friend MatrixND const diagonal(VectorND const &v);                    // Returns a diagonal matrix with the elements
                                                                        // of v in the diagonal
  static MatrixND const identity(size_t const numDimensions,            // Returns an identity matrix with given
                                 Real const &r = 1.0);                  // number of dimensions (times r)
  friend Real const normFrobenius(MatrixND const &a);                   // Returns the Frobenius norm

private:

  std::vector<VectorND> matrixND_;
};

/*
** End of class MatrixND
*/

// Inlines

inline MatrixND::MatrixND()
:
matrixND_()
{}

inline MatrixND::MatrixND(size_t const numDimensions)
:
matrixND_(numDimensions, VectorND(numDimensions, 0.0))
{}

inline MatrixND::MatrixND(Matrix3D const &a)
:
matrixND_()
{
  for(size_t i = 0; i < 3; i++) matrixND_.push_back(VectorND(a.row(i)));
}

inline MatrixND::MatrixND(Matrix4D const &a)
:
matrixND_()
{
  for(size_t i = 0; i < 4; i++) matrixND_.push_back(VectorND(a.row(i)));
}


inline MatrixND::MatrixND(VectorND const &v)
:
matrixND_()
{
  size_t const numDim = v.size();
  for(size_t i = 0; i < numDim; i++)
  {
    VectorND row(numDim, 0.0);
    row[i] = v[i];
    matrixND_.push_back(row);
  }
}

inline MatrixND::MatrixND(size_t const numDimensions, Real const &r)
:
matrixND_(numDimensions, VectorND(numDimensions, 0.0))
{ 
  for(size_t i = 0; i < numDimensions; i++) matrixND_[i][i] = r;
}

inline MatrixND &MatrixND::operator=(Real const &r)
{
  size_t const numDim = matrixND_.size();
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j < numDim; j++)
      matrixND_[i][j] = (i==j)?r:0.0;
  return *this;
}

inline MatrixND &MatrixND::operator+=(MatrixND const &a)
{
  size_t const numDim = size();
  assert(numDim == a.size());
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j < numDim; j++)
      matrixND_[i][j] += a.matrixND_[i][j];
  return *this;
}

inline MatrixND &MatrixND::operator-=(MatrixND const &a)
{
  size_t const numDim = size();
  assert(numDim == a.size());
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j < numDim; j++)
      matrixND_[i][j] -= a.matrixND_[i][j];
  return *this;
}

inline MatrixND &MatrixND::operator*=(Real const &r)
{
  size_t const numDim = size();
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j< numDim; j++)
      matrixND_[i][j] *= r;
  return *this;
}

inline MatrixND &MatrixND::operator/=(Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  size_t const numDim = size();
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j< numDim; j++)
      matrixND_[i][j] *= invr;
  return *this;
}

inline VectorND const &MatrixND::operator[](size_t const i) const
{
  assert(i < size());
  return matrixND_[i];
}

inline Real &MatrixND::operator()(size_t const i, size_t const j)
{
  size_t const numDim = size();
  assert(i < numDim && j < numDim);
  return matrixND_[i][j];
}

inline Real const &MatrixND::operator()(size_t const i, size_t const j) const
{
  size_t const numDim = size();
  assert(i < numDim && j < numDim);
  return matrixND_[i][j];
}

inline VectorND const MatrixND::row(size_t const i) const
{
  assert(i < size());
  return matrixND_[i];
}

inline VectorND const MatrixND::column(size_t const i) const
{
  size_t const numDim = size();
  assert(i < numDim);
  VectorND column;
  for(size_t j = 0; j < numDim; j++)
  {
    column.push_back(matrixND_[j][i]);
  }
  return column;
}

inline void MatrixND::setRow(size_t const i, VectorND const &row)
{
  size_t const numDim = size();
  assert(numDim == row.size() && i < numDim);
  matrixND_[i] = row;
}

inline void MatrixND::setColumn(size_t const i, VectorND const &column)
{
  size_t const numDim = size();
  assert(numDim == column.size() && i < numDim);
  for(size_t j = 0; j < numDim; j++)
    matrixND_[j][i] = column[j];
}

inline void MatrixND::clear()
{
  matrixND_.clear();
}

inline void MatrixND::resize(size_t const numDimensions)
{
  for(size_t i = 0; i < size(); i++) matrixND_[i].resize(numDimensions);
  matrixND_.resize(numDimensions, VectorND(numDimensions, 0.0));
}

inline size_t const MatrixND::size() const
{
  return matrixND_.size();
}

inline MatrixND operator+(MatrixND const &a1, MatrixND const &a2)
{
  return MatrixND(a1) += a2;
}

inline MatrixND operator-(MatrixND const &a)
{
  return MatrixND(a.size()) -= a;
}

inline MatrixND operator-(MatrixND const &a1, MatrixND const &a2)
{
  return MatrixND(a1) -= a2;
}

inline MatrixND operator*(Real const &r, MatrixND const &a)
{
  return MatrixND(a) *= r;
}

inline MatrixND operator*(MatrixND const &a, Real const &r)
{
  return MatrixND(a) *= r;
}

inline MatrixND operator/(MatrixND const &a, Real const &r)
{
  return MatrixND(a) /= r;
}

inline MatrixND operator*(MatrixND const &a1, MatrixND const &a2)
{
  size_t const numDim = a1.size();
  assert(numDim == a2.size());
  MatrixND product(numDim);
  for(size_t i = 0; i < numDim; i++)
  for(size_t j = 0; j < numDim; j++)
  {
    Real sum = 0.0;
    for(size_t k = 0; k < numDim; k++) sum += a1(i, k)*a2(k, j);
    product(i, j) = sum;
  }
  return product;
}

inline VectorND operator*(MatrixND const &a, VectorND const &v)
{
  size_t const numDim = a.size();
  assert(numDim == v.size());
  VectorND product;
  for(size_t i = 0; i < numDim; i++)
  {
    Real sum = 0.0;
    for(size_t j = 0; j < numDim; j++) sum += a(i, j)*v(j);
    product.push_back(sum);
  }
  return product;
}

inline VectorND operator*(VectorND const &v, MatrixND const &a)
{
  size_t const numDim = a.size();
  assert(numDim == v.size());
  VectorND product;
  for(size_t i = 0; i < numDim; i++)
  {
    Real sum = 0.0;
    for(size_t j = 0; j < numDim; j++) sum += a(j, i)*v(j);
    product.push_back(sum);
  }
  return product;
}

inline bool operator==(MatrixND const &a1, MatrixND const &a2)
{
  size_t const numDim = a1.size();
  if(numDim != a2.size()) return false;
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j < numDim; j++)
      if(a1(i, j) != a2(i, j)) return false;
  return true;
}

inline bool operator!=(MatrixND const &a1, MatrixND const &a2)
{
  size_t const numDim = a1.size();
  if(numDim != a2.size()) return true;
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j < numDim; j++)
      if(a1(i, j) != a2(i, j)) return true;
  return false;
}

inline bool operator==(MatrixND const &a, Real const &r)
{
  size_t const numDim = a.size();
  for(size_t i = 0; i < numDim; i++)
  {
    if(a(i, i) != r) return false;
    for(size_t j = 0; j < i; j++)
      if(a(i, j) != 0.0) return false;
    for(size_t j = i + 1; j < numDim; j++)
      if(a(i, j) != 0.0) return false;
  }
  return true;
}

inline bool operator!=(MatrixND const &a, Real const &r)
{
  size_t const numDim = a.size();
  for(size_t i = 0; i < numDim; i++)
  {
    if(a(i, i) != r) return true;
    for(size_t j = 0; j < i; j++)
      if(a(i, j) != 0.0) return true;
    for(size_t j = i + 1; j < numDim; j++)
      if(a(i, j) != 0.0) return true;
  }
  return false;
}

inline std::ostream& operator<<(std::ostream &str, MatrixND const &a)
{
  size_t const numDim = a.size();
  for(size_t i = 0; i < numDim; i++)
  {
    for(size_t j = 0; j < numDim; j++)
      str << a(i, j) << " ";
    str << std::endl;
  }
  return str;
}

inline std::istream& operator>>(std::istream &str, MatrixND &a)
{
  a.clear();
  size_t currentRow = 0;
  while(!str.eof())
  {
    std::string line;
    std::getline(str, line);
    if(str.eof()) break;
    std::istringstream lineStream(line);
    if(line.find_first_not_of(" \n\t") == line.npos) break; // To handle trailing whitespace
    VectorND row;
    while(!lineStream.eof())
    {
      Real v_i;
      lineStream >> v_i;
      if(!lineStream) break; // To handle trailing characters
      row.push_back(v_i);
    }
    if(currentRow == 0) a.resize(row.size());
    if(currentRow > a.size() || row.size() != a.size())
    {
      std::cerr << "Error reading matrix - matrix not square" << std::endl;
      a.clear();
      return str;
    }
    a.setRow(currentRow, row);
    ++currentRow;
  }
  if(!str && !str.eof())
  {
    std::cerr << "Error reading matrix from stream" << std::endl;
    a.clear();
  }
  return str;
}

inline MatrixND const outer(VectorND const &v1, VectorND const &v2)
{
  size_t const numDim = v1.size();
  assert(numDim == v2.size());
  MatrixND product(numDim, 0.0);
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j < numDim; j++) product(i, j) = v1(i)*v2(j);
  return product;
}

inline Real const trace(MatrixND const &a)
{
  Real tr = 0.0;
  for(size_t i = 0; i < a.size(); i++)
    tr += a(i, i);
  return tr;
}

inline MatrixND const transpose(MatrixND const &a)
{
  size_t const numDim = a.size();
  MatrixND t(numDim, 0.0);
  for(size_t i = 0; i < numDim; i++)
    for(size_t j = 0; j < numDim; j++) t(i, j) = a(j, i);
  return t;
}
 
inline VectorND const diagonal(MatrixND const &a)
{
  size_t const numDim = a.size();
  VectorND d(numDim, 0.0);
  for(size_t i = 0; i < numDim; i++) d[i] = a(i, i);
  return d;
}
 
inline MatrixND const diagonal(VectorND const &v)
{
  return MatrixND(v);
}

inline MatrixND const MatrixND::identity(size_t const numDimensions, Real const &r)
{
  return MatrixND(numDimensions, r);
}

inline Real const normFrobenius(MatrixND const &a)
{
  Real norm2 = 0; size_t const numDim = a.size();
  for(size_t i = 0; i < numDim; ++i)
  for(size_t j = 0; j < numDim; ++j)
    norm2 += a(i, j)*a(i, j);
  return sqrt(norm2);
}

#endif
