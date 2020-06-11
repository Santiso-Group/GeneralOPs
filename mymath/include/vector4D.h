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
** Simple 4D vector class
**
** Note that all data members (w, x, y, z) are public.
*/

#ifndef H_VECTOR4D
#define H_VECTOR4D

#include <iostream>
#include <cmath>
#include "common/include/types.h"
#include "common/include/assert.h"

class Vector4D
{
public:

  Real  w, x, y, z;  // Vector components

// Constructors

  Vector4D();                          // Defines a null vector
  Vector4D(Real const &r);             // Defines a vector with all coordinates equal to a number
  Vector4D(Real const &neww,           // Defines a vector by giving the coordinates
           Real const &newx, 
           Real const &newy, 
           Real const &newz); 
  Vector4D(Real const v[4]);           // Defines a vector from an array of reals
  Vector4D(SmallReal const v[4]);      // Defines a vector from an array of small reals

// Assignment

  Vector4D &operator+=(Vector4D const &v); // Assignment by addition
  Vector4D &operator-=(Vector4D const &v); // Assignment by subtraction
  Vector4D &operator*=(Real const &r);     // Assignment by multiplication by scalar
  Vector4D &operator/=(Real const &r);     // Assignment by division by scalar

// Access elements by index
  
  Real &operator[](size_t const i);             
  Real const &operator[](size_t const i) const; 
  Real &operator()(size_t const i);             
  Real const &operator()(size_t const i) const; 

// Size (for compatibility with templates)

  size_t const size() const; // Returns the size of the vector

// Normalization

  void normalize();      // Normalizes the vector

// Return as array

  void getArray(Real array[4]);       // Copies the vector into an array of reals
  void getArray(SmallReal array[4]);  // Copies the vector into an array of small reals

// Operators

  friend Vector4D operator+(Vector4D const &v1, Vector4D const &v2);     // Addition
  friend Vector4D operator-(Vector4D const &v);                          // Negation
  friend Vector4D operator-(Vector4D const &v1, Vector4D const &v2);     // Subtraction
  friend Vector4D operator*(Real const &r, Vector4D const &v);           // Multiplication by scalar
  friend Vector4D operator*(Vector4D const &v, Real const &r);           // Multiplication by scalar
  friend Vector4D operator/(Vector4D const &v, Real const &r);           // Division by scalar
  friend Real operator*(Vector4D const &v1, Vector4D const &v2);         // Dot product
  friend bool operator==(Vector4D const &v1, Vector4D const &v2);        // Comparison (equal)
  friend bool operator!=(Vector4D const &v1, Vector4D const &v2);        // Comparison (not equal)
  friend bool operator==(Vector4D const &v, Real const &r);              // Comparison with scalar
  friend bool operator!=(Vector4D const &v, Real const &r);              // Comparison with scalar
  friend std::ostream& operator<<(std::ostream &str, Vector4D const &v); // Output

// Other operations

  friend Real const dot(Vector4D const &v1, Vector4D const &v2);        // Dot product
  friend Real const norm(Vector4D const &v);                            // Norm
  friend Real const norm2(Vector4D const &v);                           // Norm squared
  friend Vector4D const unit(Vector4D const &v);                        // Unit vector in the direction of v
};

/*
** End of class Vector4D
*/

// Inlines

inline Vector4D::Vector4D() 
{
  w = x = y = z = 0.0;
}

inline Vector4D::Vector4D(Real const &r)
{
  w = x = y = z = r;
}

inline Vector4D::Vector4D(Real const &neww, Real const &newx, Real const &newy, Real const &newz)
{
  w = neww; x = newx; y = newy; z = newz;
}

inline Vector4D::Vector4D(Real const v[4])
{
  w = v[0]; x = v[1]; y = v[2]; z = v[3];
}

inline Vector4D::Vector4D(SmallReal const v[4])
{
  w = (Real)v[0]; x = (Real)v[1]; y = (Real)v[2]; z = (Real)v[3];
}

inline Vector4D &Vector4D::operator+=(Vector4D const &v)
{
  w += v.w; x += v.x; y += v.y; z += v.z;
  return *this;
}

inline Vector4D &Vector4D::operator-=(Vector4D const &v)
{
  w -= v.w; x -= v.x; y -= v.y; z -= v.z;
  return *this;
}

inline Vector4D &Vector4D::operator*=(Real const &r)
{
  w *= r; x *= r; y *= r; z*= r;
  return *this;
}

inline Vector4D &Vector4D::operator/=(Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  w *= invr; x *= invr; y *= invr; z *= invr;
  return *this;
}

inline Real &Vector4D::operator[](size_t const i)
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline Real const &Vector4D::operator[](size_t const i) const
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline Real &Vector4D::operator()(size_t const i)
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline Real const &Vector4D::operator()(size_t const i) const
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline size_t const Vector4D::size() const
{
  return 4;
}

inline void Vector4D::normalize()
{
  Real length = sqrt(w*w + x*x + y*y + z*z);
  assert(length);
  Real invLength = 1.0/length;
  w *= invLength; x *= invLength; y *= invLength; z *= invLength;
}

inline void Vector4D::getArray(Real array[4])
{
  array[0] = w; array[1] = x; array[2] = y; array[3] = z;
}

inline void Vector4D::getArray(SmallReal array[4])
{
  array[0] = (SmallReal)w;
  array[1] = (SmallReal)x;
  array[2] = (SmallReal)y;
  array[3] = (SmallReal)z;
}

inline Vector4D operator+(Vector4D const &v1, Vector4D const &v2)
{
  return Vector4D(v1.w + v2.w, v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector4D operator-(Vector4D const &v)
{
  return Vector4D(-v.w, -v.x, -v.y, -v.z);
}

inline Vector4D operator-(Vector4D const &v1, Vector4D const &v2)
{
  return Vector4D(v1.w - v2.w, v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector4D operator*(Real const &r, Vector4D const &v)
{
  return Vector4D(r*v.w, r*v.x, r*v.y, r*v.z);
}

inline Vector4D operator*(Vector4D const &v, Real const &r)
{
  return Vector4D(r*v.w, r*v.x, r*v.y, r*v.z);
}

inline Vector4D operator/(Vector4D const &v, Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  return Vector4D(invr*v.w, invr*v.x, invr*v.y, invr*v.z);
}

inline Real operator*(Vector4D const &v1, Vector4D const &v2)
{
  return v1.w*v2.w + v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline bool operator==(Vector4D const &v1, Vector4D const &v2)
{
  return v1.w == v2.w && v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

inline bool operator!=(Vector4D const &v1, Vector4D const &v2)
{
  return v1.w != v2.w || v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
}

inline bool operator==(Vector4D const &v, Real const &r)
{
  return v.w == r && v.x == r && v.y == r && v.z == r;
}

inline bool operator!=(Vector4D const &v, Real const &r)
{
  return v.w != r || v.x != r || v.y != r || v.z != r;
}

inline std::ostream& operator<<(std::ostream &str, Vector4D const &v)
{
  return str << v.w << " " << v.x << " " << v.y << " " << v.z;
}

inline Real const dot(Vector4D const &v1, Vector4D const &v2)
{
  return v1.w*v2.w + v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline Real const norm(Vector4D const &v)
{
  return sqrt(v.w*v.w + v.x*v.x + v.y*v.y + v.z*v.z);
}

inline Real const norm2(Vector4D const &v)
{
  return v.w*v.w + v.x*v.x + v.y*v.y + v.z*v.z;
}

inline Vector4D const unit(Vector4D const &v)
{
  Real length = sqrt(v.w*v.w + v.x*v.x + v.y*v.y + v.z*v.z);
  assert(length);
  Real invLength = 1.0/length;
  return Vector4D(invLength*v.w, invLength*v.x, invLength*v.y, invLength*v.z);
}

#endif
