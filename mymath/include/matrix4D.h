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
** Simple 4D matrix class
*/

#ifndef H_MATRIX4D
#define H_MATRIX4D

#include <iostream>
#include <cmath>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/vector4D.h"
#include "mymath/include/quaternion.h"

class Matrix4D
{
public:

// Constructors

  Matrix4D();                         // Defines a null matrix
  Matrix4D(Real const &r);            // Defines a matrix with all diagonal elements equal to a number
  Matrix4D(Vector4D const &v);        // Defines a diagonal matrix given a vector of diagonal elements
  Matrix4D(Vector4D const &v1,        // Defines a matrix by giving the columns as vectors
           Vector4D const &v2,  
           Vector4D const &v3,
           Vector4D const &v4);
  Matrix4D(Real const a[16]);         // Defines a matrix from an array of reals
  Matrix4D(SmallReal const a[16]);    // Defines a matrix from an array of small reals
  Matrix4D(Real const a[4][4]);       // Defines a matrix from a 2-dimensional array of reals
  Matrix4D(SmallReal const a[4][4]);  // Defines a matrix from a 2-dimensional array of small reals

// Size (for compatibility with templates)

  size_t const size() const; // Returns the size of the matrix

// Assignment

  Matrix4D &operator+=(Matrix4D const &a);  // Assignment by addition
  Matrix4D &operator-=(Matrix4D const &a);  // Assignment by subtraction
  Matrix4D &operator*=(Real const &r);      // Assignment by multiplication by scalar
  Matrix4D &operator/=(Real const &r);      // Assignment by division by scalar

// Access elements by index
  
  Real &operator[](size_t const i);                             // Can modify
  Real const &operator[](size_t const i) const;                 // Cannot modify
  Real &operator()(size_t const i, size_t const j);             // Can modify
  Real const &operator()(size_t const i, size_t const j) const; // Cannot modify

// Access by rows or columns

  Vector4D const row(size_t const i) const;                // Returns the i-th row as a vector
  Vector4D const column(size_t const i) const;             // Returns the i-th column as a vector
  void setRow(size_t const i, Vector4D const &row);        // Replaces the i-th row by the given vector
  void setColumn(size_t const i, Vector4D const &column);  // Replaces the i-th column by the given vector

// Return as array

  void getArray(Real array[16]);         // Copies the matrix into an array of reals
  void getArray(SmallReal array[16]);    // Copies the matrix into an array of small reals
  void getArray(Real array[4][4]);      // Copies the matrix into a 2-dimensional array of reals
  void getArray(SmallReal array[4][4]); // Copies the matrix into a 2-dimensional array of small reals

// Operators

  friend Matrix4D operator+(Matrix4D const &a1, Matrix4D const &a2);      // Addition
  friend Matrix4D operator-(Matrix4D const &a);                           // Negation
  friend Matrix4D operator-(Matrix4D const &a1, Matrix4D const &a2);      // Subtraction
  friend Matrix4D operator*(Real const &r, Matrix4D const &a);            // Multiplication by scalar
  friend Matrix4D operator*(Matrix4D const &a, Real const &r);            // Multiplication by scalar
  friend Matrix4D operator/(Matrix4D const &a, Real const &r);            // Division by scalar
  friend Matrix4D operator*(Matrix4D const &a1, Matrix4D const &a2);      // Matrix product
  friend Vector4D operator*(Matrix4D const &a, Vector4D const &v);        // Multiplication by vector
  friend Vector4D operator*(Vector4D const &v, Matrix4D const &a);        // Multiplication by vector
  friend Matrix4D operator^(Matrix4D const &a, int const n);              // Integer power of a matrix (use only for small n)
  friend bool operator==(Matrix4D const &a1, Matrix4D const &a2);         // Comparison (equal)
  friend bool operator!=(Matrix4D const &a1, Matrix4D const &a2);         // Comparison (not equal)
  friend bool operator==(Matrix4D const &a1, Real const &r);              // Comparison with scalar
  friend bool operator!=(Matrix4D const &a1, Real const &r);              // Comparison with scalar
  friend std::ostream& operator<<(std::ostream &str, Matrix4D const &a);  // Output

// Other operations

  friend Matrix4D const outer(Vector4D const &v1, Vector4D const &v2);        // Outer (Kronecker, dyadic) product
  friend Matrix4D const outer(Quaternion const &q1, Quaternion const &q2);    // Outer (Kronecker, dyadic) product 
                                                                              // for quaternions
  friend Real const trace(Matrix4D const &a);                                 // Returns the trace of the matrix
  friend Real const determinant(Matrix4D const &a);                           // Returns the determinant of the matrix
  friend Matrix4D const transpose(Matrix4D const &a);                         // Returns the transpose of the matrix
  friend Matrix4D const inverse(Matrix4D const &a);                           // Returns the inverse
  friend Vector4D const diagonal(Matrix4D const &a);                          // Returns the diagonal of the matrix
  friend Matrix4D const diagonal(Vector4D const &v);                          // Returns a diagonal matrix with the 
                                                                              // elements of v in the diagonal
  static Matrix4D const identity(Real const &r = 1.0);                        // Returns an identity matrix (times r)

private:

  Real  elements_[16];  // Matrix elements  
};

/*
** End of class Matrix4D
*/

// Inlines

inline Matrix4D::Matrix4D()
{
  for(size_t i = 0; i < 16; i++) elements_[i] = 0.0;
}

inline Matrix4D::Matrix4D(Real const &r)
{
  elements_[ 0] =   r; elements_[ 4] = 0.0; elements_[ 8] = 0.0; elements_[12] = 0.0;
  elements_[ 1] = 0.0; elements_[ 5] =   r; elements_[ 9] = 0.0; elements_[13] = 0.0;
  elements_[ 2] = 0.0; elements_[ 6] = 0.0; elements_[10] =   r; elements_[14] = 0.0;
  elements_[ 3] = 0.0; elements_[ 7] = 0.0; elements_[11] = 0.0; elements_[15] =   r;
}

inline Matrix4D::Matrix4D(Vector4D const &v)
{
  elements_[ 0] = v.w; elements_[ 4] = 0.0; elements_[ 8] = 0.0; elements_[12] = 0.0;
  elements_[ 1] = 0.0; elements_[ 5] = v.x; elements_[ 9] = 0.0; elements_[13] = 0.0;
  elements_[ 2] = 0.0; elements_[ 6] = 0.0; elements_[10] = v.y; elements_[14] = 0.0;
  elements_[ 3] = 0.0; elements_[ 7] = 0.0; elements_[11] = 0.0; elements_[15] = v.z;
}

inline Matrix4D::Matrix4D(Vector4D const &v1, Vector4D const &v2, Vector4D const &v3, Vector4D const &v4)
{
  elements_[ 0] = v1.w; elements_[ 4] = v2.w; elements_[ 8] = v3.w; elements_[12] = v4.w;
  elements_[ 1] = v1.x; elements_[ 5] = v2.x; elements_[ 9] = v3.x; elements_[13] = v4.x;
  elements_[ 2] = v1.y; elements_[ 6] = v2.y; elements_[10] = v3.y; elements_[14] = v4.y;
  elements_[ 3] = v1.z; elements_[ 7] = v2.z; elements_[11] = v3.z; elements_[15] = v4.z;
}

inline Matrix4D::Matrix4D(Real const a[16])
{
  for(size_t i = 0; i < 16; i++) elements_[i] = a[i];
}

inline Matrix4D::Matrix4D(SmallReal const a[16])
{
  for(size_t i = 0; i < 16; i++) elements_[i] = (Real)a[i];
}

inline Matrix4D::Matrix4D(Real const a[4][4])
{
  elements_[ 0] = a[0][0]; elements_[ 4] = a[0][1]; elements_[ 8] = a[0][2]; elements_[12] = a[0][3];
  elements_[ 1] = a[1][0]; elements_[ 5] = a[1][1]; elements_[ 9] = a[1][2]; elements_[13] = a[1][3];
  elements_[ 2] = a[2][0]; elements_[ 6] = a[2][1]; elements_[10] = a[2][2]; elements_[14] = a[2][3];
  elements_[ 3] = a[3][0]; elements_[ 7] = a[3][1]; elements_[11] = a[3][2]; elements_[15] = a[3][3];
}

inline Matrix4D::Matrix4D(SmallReal const a[4][4])
{
  elements_[ 0] = (Real)a[0][0]; elements_[ 4] = (Real)a[0][1]; elements_[ 8] = (Real)a[0][2]; elements_[12] = (Real)a[0][3];
  elements_[ 1] = (Real)a[1][0]; elements_[ 5] = (Real)a[1][1]; elements_[ 9] = (Real)a[1][2]; elements_[13] = (Real)a[1][3];
  elements_[ 2] = (Real)a[2][0]; elements_[ 6] = (Real)a[2][1]; elements_[10] = (Real)a[2][2]; elements_[14] = (Real)a[2][3];
  elements_[ 3] = (Real)a[3][0]; elements_[ 7] = (Real)a[3][1]; elements_[11] = (Real)a[3][2]; elements_[15] = (Real)a[3][3];
}

inline size_t const Matrix4D::size() const
{
  return 4;
}

inline Matrix4D &Matrix4D::operator+=(Matrix4D const &a)
{
  for(size_t i = 0; i < 16; i++) elements_[i] += a.elements_[i];
  return *this;
}

inline Matrix4D &Matrix4D::operator-=(Matrix4D const &a)
{
  for(size_t i = 0; i < 16; i++) elements_[i] -= a.elements_[i];
  return *this;
}

inline Matrix4D &Matrix4D::operator*=(Real const &r)
{
  for(size_t i = 0; i < 16; i++) elements_[i] *= r;
  return *this;
}

inline Matrix4D &Matrix4D::operator/=(Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  for(size_t i = 0; i < 16; i++) elements_[i] *= invr;
  return *this;
}

inline Real &Matrix4D::operator[](size_t const i)
{
  assert(i < 16);
  return elements_[i];
}

inline Real const &Matrix4D::operator[](size_t const i) const
{
  assert(i < 16);
  return elements_[i];
}

inline Real &Matrix4D::operator()(size_t const i, size_t const j)
{
  assert(i < 4 && j < 4);
  return elements_[i + 4*j];
}

inline Real const &Matrix4D::operator()(size_t const i, size_t const j) const
{
  assert(i < 4 && j < 4);
  return elements_[i + 4*j];
}

inline Vector4D const Matrix4D::row(size_t const i) const
{
  assert(i < 4);
  return Vector4D(elements_[i], elements_[i + 4], elements_[i + 8], elements_[i + 12]);
}

inline Vector4D const Matrix4D::column(size_t const i) const
{
  assert(i < 4);
  size_t fi = 4*i;
  return Vector4D(elements_[fi], elements_[fi + 1], elements_[fi + 2], elements_[fi + 3]);
}

inline void Matrix4D::setRow(size_t const i, Vector4D const &row)
{
  assert(i < 4);
  elements_[i     ] = row[0];
  elements_[i +  4] = row[1];
  elements_[i +  8] = row[2];
  elements_[i + 12] = row[3];
}
 
inline void Matrix4D::setColumn(size_t const i, Vector4D const &column)
{
  assert(i < 4);
  size_t fi = 4*i;
  elements_[fi    ] = column[0];
  elements_[fi + 1] = column[1];
  elements_[fi + 2] = column[2];
  elements_[fi + 3] = column[3];
}

inline void Matrix4D::getArray(Real array[16])
{
  for(size_t i = 0; i < 16; i++) array[i] = elements_[i];
}

inline void Matrix4D::getArray(SmallReal array[16])
{
  for(size_t i = 0; i < 16; i++) array[i] = (SmallReal)elements_[i];
}

inline void Matrix4D::getArray(Real array[4][4])
{
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++) array[i][j] = elements_[i + 4*j];
}

inline void Matrix4D::getArray(SmallReal array[4][4])
{
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++) array[i][j] = (SmallReal)elements_[i + 4*j];
}

inline Matrix4D operator+(Matrix4D const &a1, Matrix4D const &a2)
{
  return Matrix4D(a1) += a2;
}

inline Matrix4D operator-(const Matrix4D &a)
{
  return Matrix4D() -= a;
}

inline Matrix4D operator-(Matrix4D const &a1, Matrix4D const &a2)
{
  return Matrix4D(a1) -= a2;
}

inline Matrix4D operator*(Real const &r, Matrix4D const &a)
{
  return Matrix4D(a) *= r;
}

inline Matrix4D operator*(Matrix4D const &a, Real const &r)
{
  return Matrix4D(a) *= r;
}

inline Matrix4D operator/(Matrix4D const &a, Real const &r)
{
  return Matrix4D(a) /= r;
}

inline Matrix4D operator*(Matrix4D const &a1, Matrix4D const &a2)
{
  Matrix4D product;
  for(size_t i = 0; i < 4; i++) {
  for(size_t j = 0; j < 4; j++) {
    Real sum = 0.0;
    for(size_t k = 0; k < 4; k++) sum += a1(i,k)*a2(k,j);
    product(i, j) = sum;
  }}
  return product;
}

inline Vector4D operator*(Matrix4D const &a, Vector4D const &v)
{
  Vector4D product;
  for(size_t i = 0; i < 4; i++) {
    Real sum = 0.0;
    for(size_t j = 0; j < 4; j++) sum += a(i, j)*v(j);
    product(i) = sum;
  }
  return product;
}

inline Vector4D operator*(Vector4D const &v, Matrix4D const &a)
{
  Vector4D product;
  for(size_t i = 0; i < 4; i++) {
    Real sum = 0.0;
    for(size_t j = 0; j < 4; j++) sum += a(j, i)*v(j);
    product(i) = sum;
  }
  return product;
}

inline Matrix4D operator^(Matrix4D const &a, int const n)
{

  // Note: This is not very efficient if n can be large.
  // A better option would be to use the prime factor 
  // decomposition of n.

  if(n == 0) return Matrix4D::identity();
  Matrix4D power(1.0), base((n > 0)?a:inverse(a));
  size_t ex = (n > 0)?n:-n;
  for(size_t i = 0; i < ex; i++) power = power*base;
  return power;
}

inline Matrix4D const outer(Vector4D const &v1, Vector4D const &v2)
{
  Matrix4D product;
  for(size_t i = 0; i < 4; i++) 
    for(size_t j = 0; j < 4; j++) product(i, j) = v1(i)*v2(j);
  return product;
}

inline Matrix4D const outer(Quaternion const &q1, Quaternion const &q2)
{
  Matrix4D product;
  product[ 0] = q1.w*q2.w; product[ 4] = q1.w*q2.x; product[ 8] = q1.w*q2.y; product[12] = q1.w*q2.z;
  product[ 1] = q1.x*q2.w; product[ 5] = q1.x*q2.x; product[ 9] = q1.x*q2.y; product[13] = q1.x*q2.z;
  product[ 2] = q1.y*q2.w; product[ 6] = q1.y*q2.x; product[10] = q1.y*q2.y; product[14] = q1.y*q2.z;
  product[ 3] = q1.z*q2.w; product[ 7] = q1.z*q2.x; product[11] = q1.z*q2.y; product[15] = q1.z*q2.z;
  return product;
}

inline bool operator==(Matrix4D const &a1, Matrix4D const &a2)
{
  bool result = true;
  for(size_t i = 0; i < 16; i++) result &= (a1.elements_[i] == a2.elements_[i]);
  return result;
}

inline bool operator!=(Matrix4D const &a1, Matrix4D const &a2)
{
  bool result = false;
  for(size_t i = 0; i < 16; i++) result |= (a1.elements_[i] != a2.elements_[i]);
  return result;
}

inline bool operator==(Matrix4D const &a1, Real const &r)
{
  return (a1 == Matrix4D(r));
}

inline bool operator!=(Matrix4D const &a1, Real const &r)
{
  return (a1 != Matrix4D(r));
}

inline std::ostream& operator<<(std::ostream &str, Matrix4D const &a)
{
  return str << std::endl 
    << a.elements_[ 0] << " " << a.elements_[ 4] << " " << a.elements_[ 8] << " " << a.elements_[12] << std::endl
    << a.elements_[ 1] << " " << a.elements_[ 5] << " " << a.elements_[ 9] << " " << a.elements_[13] << std::endl
    << a.elements_[ 2] << " " << a.elements_[ 6] << " " << a.elements_[10] << " " << a.elements_[14] << std::endl
    << a.elements_[ 3] << " " << a.elements_[ 7] << " " << a.elements_[11] << " " << a.elements_[15];
}

inline Real const trace(Matrix4D const &a)
{
  return a.elements_[0] + a.elements_[5] + a.elements_[10] + a.elements_[15];
}

// The #defines in the routines below are not necessary, but they improve speed a lot

inline Real const determinant(Matrix4D const &a)
{
#define a a.elements_

  return a[ 0]*(a[ 5]*(a[10]*a[15] - a[14]*a[11]) + a[ 9]*(a[14]*a[ 7] - a[ 6]*a[15]) + a[13]*(a[ 6]*a[11] - a[10]*a[ 7])) -
         a[ 4]*(a[ 1]*(a[10]*a[15] - a[14]*a[11]) + a[ 9]*(a[14]*a[ 3] - a[ 2]*a[15]) + a[13]*(a[ 2]*a[11] - a[10]*a[ 3])) +
         a[ 8]*(a[ 1]*(a[ 6]*a[15] - a[14]*a[ 7]) + a[ 5]*(a[14]*a[ 3] - a[ 2]*a[15]) + a[13]*(a[ 2]*a[ 7] - a[ 6]*a[ 3])) -
         a[12]*(a[ 1]*(a[ 6]*a[11] - a[10]*a[ 7]) + a[ 5]*(a[10]*a[ 3] - a[ 2]*a[11]) + a[ 9]*(a[ 2]*a[ 7] - a[ 6]*a[ 3]));
#undef a
}

inline Matrix4D const transpose(Matrix4D const &a)
{
  Matrix4D t;

#define a a.elements_
#define t t.elements_

  t[ 0] = a[ 0]; t[ 4] = a[ 1]; t[ 8] = a[ 2]; t[12] = a[ 3]; 
  t[ 1] = a[ 4]; t[ 5] = a[ 5]; t[ 9] = a[ 6]; t[13] = a[ 7];
  t[ 2] = a[ 8]; t[ 6] = a[ 9]; t[10] = a[10]; t[14] = a[11];
  t[ 3] = a[12]; t[ 7] = a[13]; t[11] = a[14]; t[15] = a[15];

#undef t
#undef a

  return t;
}

inline Matrix4D const inverse(Matrix4D const &a)
{
#define a a.elements_

  Vector4D col1(a[ 5]*(a[10]*a[15] - a[14]*a[11]) + a[ 9]*(a[14]*a[ 7] - a[ 6]*a[15]) + a[13]*(a[ 6]*a[11] - a[10]*a[ 7]),
                a[ 1]*(a[14]*a[11] - a[10]*a[15]) + a[ 9]*(a[ 2]*a[15] - a[14]*a[ 3]) + a[13]*(a[10]*a[ 3] - a[ 2]*a[11]),
                a[ 1]*(a[ 6]*a[15] - a[14]*a[ 7]) + a[ 5]*(a[14]*a[ 3] - a[ 2]*a[15]) + a[13]*(a[ 2]*a[ 7] - a[ 6]*a[ 3]),
                a[ 1]*(a[10]*a[ 7] - a[ 6]*a[11]) + a[ 5]*(a[ 2]*a[11] - a[10]*a[ 3]) + a[ 9]*(a[ 6]*a[ 3] - a[ 2]*a[ 7]));

  Real det = a[0]*col1.w + a[4]*col1.x + a[8]*col1.y + a[12]*col1.z;

  assert(det);

  Vector4D col2(a[ 4]*(a[14]*a[11] - a[10]*a[15]) + a[ 8]*(a[ 6]*a[15] - a[14]*a[ 7]) + a[12]*(a[10]*a[ 7] - a[ 6]*a[11]),
                a[ 0]*(a[10]*a[15] - a[14]*a[11]) + a[ 8]*(a[14]*a[ 3] - a[ 2]*a[15]) + a[12]*(a[ 2]*a[11] - a[10]*a[ 3]),
                a[ 0]*(a[14]*a[ 7] - a[ 6]*a[15]) + a[ 4]*(a[ 2]*a[15] - a[14]*a[ 3]) + a[12]*(a[ 6]*a[ 3] - a[ 2]*a[ 7]),
                a[ 0]*(a[ 6]*a[11] - a[10]*a[ 7]) + a[ 4]*(a[10]*a[ 3] - a[ 2]*a[11]) + a[ 8]*(a[ 2]*a[ 7] - a[ 6]*a[ 3]));

  Vector4D col3(a[ 4]*(a[ 9]*a[15] - a[13]*a[11]) + a[ 8]*(a[13]*a[ 7] - a[ 5]*a[15]) + a[12]*(a[ 5]*a[11] - a[ 9]*a[ 7]),
                a[ 0]*(a[13]*a[11] - a[ 9]*a[15]) + a[ 8]*(a[ 1]*a[15] - a[13]*a[ 3]) + a[12]*(a[ 9]*a[ 3] - a[ 1]*a[11]),
                a[ 0]*(a[ 5]*a[15] - a[13]*a[ 7]) + a[ 4]*(a[13]*a[ 3] - a[ 1]*a[15]) + a[12]*(a[ 1]*a[ 7] - a[ 5]*a[ 3]),
                a[ 0]*(a[ 9]*a[ 7] - a[ 5]*a[11]) + a[ 4]*(a[ 1]*a[11] - a[ 9]*a[ 3]) + a[ 8]*(a[ 5]*a[ 3] - a[ 1]*a[ 7]));

  Vector4D col4(a[ 4]*(a[13]*a[10] - a[ 9]*a[14]) + a[ 8]*(a[ 5]*a[14] - a[13]*a[ 6]) + a[12]*(a[ 9]*a[ 6] - a[ 5]*a[10]),
                a[ 0]*(a[ 9]*a[14] - a[13]*a[10]) + a[ 8]*(a[13]*a[ 2] - a[ 1]*a[14]) + a[12]*(a[ 1]*a[10] - a[ 9]*a[ 2]),
                a[ 0]*(a[13]*a[ 6] - a[ 5]*a[14]) + a[ 4]*(a[ 1]*a[14] - a[13]*a[ 2]) + a[12]*(a[ 5]*a[ 2] - a[ 1]*a[ 6]),
                a[ 0]*(a[ 5]*a[10] - a[ 9]*a[ 6]) + a[ 4]*(a[ 9]*a[ 2] - a[ 1]*a[10]) + a[ 8]*(a[ 1]*a[ 6] - a[ 5]*a[ 2]));
#undef a

  return Matrix4D(col1, col2, col3, col4)/det;
}

inline Vector4D const diagonal(Matrix4D const &a)
{
  return Vector4D(a.elements_[0], a.elements_[5], a.elements_[10], a.elements_[15]);
}

inline Matrix4D const diagonal(Vector4D const &v)
{
  return Matrix4D(v);
}

inline Matrix4D const Matrix4D::identity(Real const &r)
{
  return Matrix4D(r);
}

#endif
