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
** Simple 3D matrix class
*/

#ifndef H_MATRIX3D
#define H_MATRIX3D

#include <iostream>
#include <cmath>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/vector3D.h"

class Matrix3D
{
public:

// Constructors

  Matrix3D();                         // Defines a null matrix
  Matrix3D(Real const &r);            // Defines a matrix with all diagonal elements equal to a number
  Matrix3D(Vector3D const &v);        // Defines a diagonal matrix given a vector of diagonal elements
  Matrix3D(Vector3D const &v1,        // Defines a matrix by giving the columns as vectors
           Vector3D const &v2,  
           Vector3D const &v3);
  Matrix3D(Real const a[9]);          // Defines a matrix from an array of reals
  Matrix3D(SmallReal const a[9]);     // Defines a matrix from an array of small reals
  Matrix3D(Real const a[3][3]);       // Defines a matrix from a 2-dimensional array of reals
  Matrix3D(SmallReal const a[3][3]);  // Defines a matrix from a 2-dimensional array of small reals

// Size (for compatibility with templates)

  size_t const size() const; // Returns the size of the matrix

// Assignment

  Matrix3D &operator+=(Matrix3D const &a);  // Assignment by addition
  Matrix3D &operator-=(Matrix3D const &a);  // Assignment by subtraction
  Matrix3D &operator*=(Real const &r);      // Assignment by multiplication by scalar
  Matrix3D &operator/=(Real const &r);      // Assignment by division by scalar

// Access elements by index
  
  Real &operator[](size_t const i);                             // Can modify
  Real const &operator[](size_t const i) const;                 // Cannot modify
  Real &operator()(size_t const i, size_t const j);             // Can modify
  Real const &operator()(size_t const i, size_t const j) const; // Cannot modify

// Access by rows or columns

  Vector3D const row(size_t const i) const;                // Returns the i-th row as a vector
  Vector3D const column(size_t const i) const;             // Returns the i-th column as a vector
  void setRow(size_t const i, Vector3D const &row);        // Replaces the i-th row by the given vector
  void setColumn(size_t const i, Vector3D const &column);  // Replaces the i-th column by the given vector

// Return as array

  void getArray(Real array[9]);         // Copies the matrix into an array of reals
  void getArray(SmallReal array[9]);    // Copies the matrix into an array of small reals
  void getArray(Real array[3][3]);      // Copies the matrix into a 2-dimensional array of reals
  void getArray(SmallReal array[3][3]); // Copies the matrix into a 2-dimensional array of small reals

// Operators

  friend Matrix3D operator+(Matrix3D const &a1, Matrix3D const &a2);      // Addition
  friend Matrix3D operator-(Matrix3D const &a);                           // Negation
  friend Matrix3D operator-(Matrix3D const &a1, Matrix3D const &a2);      // Subtraction
  friend Matrix3D operator*(Real const &r, Matrix3D const &a);            // Multiplication by scalar
  friend Matrix3D operator*(Matrix3D const &a, Real const &r);            // Multiplication by scalar
  friend Matrix3D operator/(Matrix3D const &a, Real const &r);            // Division by scalar
  friend Matrix3D operator*(Matrix3D const &a1, Matrix3D const &a2);      // Matrix product
  friend Vector3D operator*(Matrix3D const &a, Vector3D const &v);        // Multiplication by vector
  friend Vector3D operator*(Vector3D const &v, Matrix3D const &a);        // Multiplication by vector
  friend Matrix3D operator^(Matrix3D const &a, int const n);              // Integer power of a matrix (use only for small n)
  friend bool operator==(Matrix3D const &a1, Matrix3D const &a2);         // Comparison (equal)
  friend bool operator!=(Matrix3D const &a1, Matrix3D const &a2);         // Comparison (not equal)
  friend bool operator==(Matrix3D const &a1, Real const &r);              // Comparison with scalar
  friend bool operator!=(Matrix3D const &a1, Real const &r);              // Comparison with scalar
  friend std::ostream& operator<<(std::ostream &str, Matrix3D const &a);  // Output

// Other operations

  friend Matrix3D const outer(Vector3D const &v1, Vector3D const &v2);        // Outer (Kronecker, dyadic) product
  friend Real const trace(Matrix3D const &a);                                 // Returns the trace of the matrix
  friend Real const determinant(Matrix3D const &a);                           // Returns the determinant of the matrix
  friend Matrix3D const transpose(Matrix3D const &a);                         // Returns the transpose of the matrix
  friend Matrix3D const inverse(Matrix3D const &a);                           // Returns the inverse
  friend Vector3D const diagonal(Matrix3D const &a);                          // Returns the diagonal of the matrix
  friend Matrix3D const diagonal(Vector3D const &v);                          // Returns a diagonal matrix with the 
                                                                              // elements of v in the diagonal
  static Matrix3D const identity(Real const &r = 1.0);                        // Returns an identity matrix (times r)
                                                                              
private:

  Real elements_[9];  // Matrix elements  
};

/*
** End of class Matrix3D
*/

// Inlines

inline Matrix3D::Matrix3D()
{
  for(size_t i = 0; i < 9; i++) elements_[i] = 0.0;
}

inline Matrix3D::Matrix3D(Real const &r)
{
  elements_[0] =   r; elements_[3] = 0.0; elements_[6] = 0.0;
  elements_[1] = 0.0; elements_[4] =   r; elements_[7] = 0.0;
  elements_[2] = 0.0; elements_[5] = 0.0; elements_[8] =   r;
}

inline Matrix3D::Matrix3D(Vector3D const &v)
{
  elements_[0] = v.x; elements_[3] = 0.0; elements_[6] = 0.0;
  elements_[1] = 0.0; elements_[4] = v.y; elements_[7] = 0.0;
  elements_[2] = 0.0; elements_[5] = 0.0; elements_[8] = v.z;
}

inline Matrix3D::Matrix3D(Vector3D const &v1, Vector3D const &v2, Vector3D const &v3)
{
  elements_[0] = v1.x; elements_[3] = v2.x; elements_[6] = v3.x;
  elements_[1] = v1.y; elements_[4] = v2.y; elements_[7] = v3.y;
  elements_[2] = v1.z; elements_[5] = v2.z; elements_[8] = v3.z;
}

inline Matrix3D::Matrix3D(Real const a[9])
{
  for(size_t i = 0; i < 9; i++) elements_[i] = a[i];
}

inline Matrix3D::Matrix3D(SmallReal const a[9])
{
  for(size_t i = 0; i < 9; i++) elements_[i] = (Real)a[i];
}

inline Matrix3D::Matrix3D(Real const a[3][3])
{
  elements_[0] = a[0][0]; elements_[3] = a[0][1]; elements_[6] = a[0][2];
  elements_[1] = a[1][0]; elements_[4] = a[1][1]; elements_[7] = a[1][2];
  elements_[2] = a[2][0]; elements_[5] = a[2][1]; elements_[8] = a[2][2];
}

inline Matrix3D::Matrix3D(SmallReal const a[3][3])
{
  elements_[0] = (Real)a[0][0]; elements_[3] = (Real)a[0][1]; elements_[6] = (Real)a[0][2];
  elements_[1] = (Real)a[1][0]; elements_[4] = (Real)a[1][1]; elements_[7] = (Real)a[1][2];
  elements_[2] = (Real)a[2][0]; elements_[5] = (Real)a[2][1]; elements_[8] = (Real)a[2][2];
}

inline size_t const Matrix3D::size() const
{
  return 3;
}

inline Matrix3D &Matrix3D::operator+=(Matrix3D const &a)
{
  for(size_t i = 0; i < 9; i++) elements_[i] += a.elements_[i];
  return *this;
}

inline Matrix3D &Matrix3D::operator-=(Matrix3D const &a)
{
  for(size_t i = 0; i < 9; i++) elements_[i] -= a.elements_[i];
  return *this;
}

inline Matrix3D &Matrix3D::operator*=(Real const &r)
{
  for(size_t i = 0; i < 9; i++) elements_[i] *= r;
  return *this;
}

inline Matrix3D &Matrix3D::operator/=(Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  for(size_t i = 0; i < 9; i++) elements_[i] *= invr;
  return *this;
}

inline Real &Matrix3D::operator[](size_t const i)
{
  assert(i < 9);
  return elements_[i];
}

inline Real const &Matrix3D::operator[](size_t const i) const
{
  assert(i < 9);
  return elements_[i];
}

inline Real &Matrix3D::operator()(size_t const i, size_t const j)
{
  assert(i < 3 && j < 3);
  return elements_[i + 3*j];
}

inline Real const &Matrix3D::operator()(size_t const i, size_t const j) const
{
  assert(i < 3 && j < 3);
  return elements_[i + 3*j];
}

inline Vector3D const Matrix3D::row(size_t const i) const
{
  assert(i < 3);
  return Vector3D(elements_[i], elements_[i + 3], elements_[i + 6]);
}

inline Vector3D const Matrix3D::column(size_t const i) const
{
  assert(i < 3);
  size_t ti = 3*i;
  return Vector3D(elements_[ti], elements_[ti + 1], elements_[ti + 2]);
}

inline void Matrix3D::setRow(size_t const i, Vector3D const &row)
{
  assert(i < 3);
  elements_[i    ] = row.x;
  elements_[i + 3] = row.y;
  elements_[i + 6] = row.z;
}
 
inline void Matrix3D::setColumn(size_t const i, Vector3D const &column)
{
  assert(i < 3);
  size_t ti = 3*i;
  elements_[ti    ] = column.x;
  elements_[ti + 1] = column.y;
  elements_[ti + 2] = column.z;
}

inline void Matrix3D::getArray(Real array[9])
{
  for(size_t i = 0; i < 9; i++) array[i] = elements_[i];
}

inline void Matrix3D::getArray(SmallReal array[9])
{
  for(size_t i = 0; i < 9; i++) array[i] = (SmallReal)elements_[i];
}

inline void Matrix3D::getArray(Real array[3][3])
{
  for(size_t i = 0; i < 3; i++)
    for(size_t j = 0; j < 3; j++) array[i][j] = elements_[i + 3*j];
}

inline void Matrix3D::getArray(SmallReal array[3][3])
{
  for(size_t i = 0; i < 3; i++)
    for(size_t j = 0; j < 3; j++) array[i][j] = (SmallReal)elements_[i + 3*j];
}

inline Matrix3D operator+(Matrix3D const &a1, Matrix3D const &a2)
{
  return Matrix3D(a1) += a2;
}

inline Matrix3D operator-(const Matrix3D &a)
{
  return Matrix3D() -= a;
}

inline Matrix3D operator-(Matrix3D const &a1, Matrix3D const &a2)
{
  return Matrix3D(a1) -= a2;
}

inline Matrix3D operator*(Real const &r, Matrix3D const &a)
{
  return Matrix3D(a) *= r;
}

inline Matrix3D operator*(Matrix3D const &a, Real const &r)
{
  return Matrix3D(a) *= r;
}

inline Matrix3D operator/(Matrix3D const &a, Real const &r)
{
  return Matrix3D(a) /= r;
}

inline Matrix3D operator*(Matrix3D const &a1, Matrix3D const &a2)
{
  Matrix3D product;
  for(size_t i = 0; i < 3; i++) {
  for(size_t j = 0; j < 3; j++) {
    Real sum = 0.0;
    for(size_t k = 0; k < 3; k++) sum += a1(i, k)*a2(k, j);
    product(i, j) = sum;
  }}
  return product;
}

inline Vector3D operator*(Matrix3D const &a, Vector3D const &v)
{
  Vector3D product;
  for(size_t i = 0; i < 3; i++) {
    Real sum = 0.0;
    for(size_t j = 0; j < 3; j++) sum += a(i, j)*v(j);
    product(i) = sum;
  }
  return product;
}

inline Vector3D operator*(Vector3D const &v, Matrix3D const &a)
{
  Vector3D product;
  for(size_t i = 0; i < 3; i++) {
    Real sum = 0.0;
    for(size_t j = 0; j < 3; j++) sum += a(j, i)*v(j);
    product(i) = sum;
  }
  return product;
}

inline Matrix3D operator^(Matrix3D const &a, int const n)
{
  // Note: This is not very efficient if n can be large.
  // A better option would be to use the prime factor 
  // decomposition of n.

  Matrix3D power(1.0), base((n > 0)?a:inverse(a));
  size_t ex = (n > 0)?n:-n;
  for(size_t i = 0; i < ex; i++) power = power*base;
  return power;
}

inline Matrix3D const outer(Vector3D const &v1, Vector3D const &v2)
{
  Matrix3D product;
  for(size_t i = 0; i < 3; i++) 
    for(size_t j = 0; j < 3; j++) product(i, j) = v1(i)*v2(j);
  return product;
}

inline bool operator==(Matrix3D const &a1, Matrix3D const &a2)
{
  return a1.elements_[0] == a2.elements_[0] &&
         a1.elements_[1] == a2.elements_[1] &&
         a1.elements_[2] == a2.elements_[2] &&
         a1.elements_[3] == a2.elements_[3] &&
         a1.elements_[4] == a2.elements_[4] &&
         a1.elements_[5] == a2.elements_[5] &&
         a1.elements_[6] == a2.elements_[6] &&
         a1.elements_[7] == a2.elements_[7] &&
         a1.elements_[8] == a2.elements_[8];
}

inline bool operator!=(Matrix3D const &a1, Matrix3D const &a2)
{
  return a1.elements_[0] != a2.elements_[0] ||
         a1.elements_[1] != a2.elements_[1] ||
         a1.elements_[2] != a2.elements_[2] ||
         a1.elements_[3] != a2.elements_[3] ||
         a1.elements_[4] != a2.elements_[4] ||
         a1.elements_[5] != a2.elements_[5] ||
         a1.elements_[6] != a2.elements_[6] ||
         a1.elements_[7] != a2.elements_[7] ||
         a1.elements_[8] != a2.elements_[8];
}

inline bool operator==(Matrix3D const &a1, Real const &r)
{
  return (a1 == Matrix3D(r));
}

inline bool operator!=(Matrix3D const &a1, Real const &r)
{
  return (a1 != Matrix3D(r));
}

inline std::ostream& operator<<(std::ostream &str, Matrix3D const &a)
{
  return str << std::endl << a.elements_[0] << " " << a.elements_[3] << " " << a.elements_[6] 
             << std::endl << a.elements_[1] << " " << a.elements_[4] << " " << a.elements_[7] 
             << std::endl << a.elements_[2] << " " << a.elements_[5] << " " << a.elements_[8];
}

inline Real const trace(Matrix3D const &a)
{
  return a.elements_[0] + a.elements_[4] + a.elements_[8];
}

// The #defines in the routines below are not necessary, but they improve speed a lot

inline Real const determinant(Matrix3D const &a)
{
#define a a.elements_

  return a[0]*(a[4]*a[8] - a[7]*a[5]) +
         a[3]*(a[7]*a[2] - a[1]*a[8]) +
         a[6]*(a[1]*a[5] - a[4]*a[2]);
#undef a
}

inline Matrix3D const transpose(Matrix3D const &a)
{
  Matrix3D t;

#define a a.elements_
#define t t.elements_

  t[0] = a[0]; t[3] = a[1]; t[6] = a[2];
  t[1] = a[3]; t[4] = a[4]; t[7] = a[5];
  t[2] = a[6]; t[5] = a[7]; t[8] = a[8];

#undef t
#undef a

  return t;
}

inline Matrix3D const inverse(Matrix3D const &a)
{
#define a a.elements_

  Vector3D col1(a[4]*a[8] - a[7]*a[5],
                a[7]*a[2] - a[1]*a[8],
                a[1]*a[5] - a[4]*a[2]);
  Real det = a[0]*col1.x+ a[3]*col1.y + a[6]*col1.z;
  assert(det);
  Vector3D col2(a[6]*a[5] - a[3]*a[8],
                a[0]*a[8] - a[6]*a[2],
                a[3]*a[2] - a[0]*a[5]);
  Vector3D col3(a[3]*a[7] - a[6]*a[4],
                a[6]*a[1] - a[0]*a[7],
                a[0]*a[4] - a[3]*a[1]);
#undef a

  return Matrix3D(col1, col2, col3)/det;
}

inline Vector3D const diagonal(Matrix3D const &a)
{
  return Vector3D(a.elements_[0], a.elements_[4], a.elements_[8]);
}

inline Matrix3D const diagonal(Vector3D const &v)
{
  return Matrix3D(v);
}

inline Matrix3D const Matrix3D::identity(Real const &r)
{
  return Matrix3D(r);
}

#endif
