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
** Simple 3D vector class
**
** Note that all data members (x, y, z) are public.
*/

#ifndef H_VECTOR3D
#define H_VECTOR3D

#include <iostream>
#include <cmath>
#include "common/include/types.h"
#include "common/include/assert.h"

class Vector3D
{
public:

  Real  x, y, z;  // Vector components

// Constructors

  Vector3D();                          // Defines a null vector
  Vector3D(Real const &r);             // Defines a vector with all coordinates equal to a number
  Vector3D(Real const &newx,           // Defines a vector by giving the coordinates
           Real const &newy, 
           Real const &newz); 
  Vector3D(Real const v[3]);           // Defines a vector from an array of reals
  Vector3D(SmallReal const v[3]);      // Defines a vector from an array of small reals

// Assignment

  Vector3D &operator+=(Vector3D const &v); // Assignment by addition
  Vector3D &operator-=(Vector3D const &v); // Assignment by subtraction
  Vector3D &operator*=(Real const &r);     // Assignment by multiplication by scalar
  Vector3D &operator/=(Real const &r);     // Assignment by division by scalar

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

  void getArray(Real array[3]);       // Copies the vector into an array of reals
  void getArray(SmallReal array[3]);  // Copies the vector into an array of small reals

// Operators

  friend Vector3D operator+(Vector3D const &v1, Vector3D const &v2);     // Addition
  friend Vector3D operator-(Vector3D const &v);                          // Negation
  friend Vector3D operator-(Vector3D const &v1, Vector3D const &v2);     // Subtraction
  friend Vector3D operator*(Real const &r, Vector3D const &v);           // Multiplication by scalar
  friend Vector3D operator*(Vector3D const &v, Real const &r);           // Multiplication by scalar
  friend Vector3D operator/(Vector3D const &v, Real const &r);           // Division by scalar
  friend Real operator*(Vector3D const &v1, Vector3D const &v2);         // Dot product
  friend bool operator==(Vector3D const &v1, Vector3D const &v2);        // Comparison (equal)
  friend bool operator!=(Vector3D const &v1, Vector3D const &v2);        // Comparison (not equal)
  friend bool operator==(Vector3D const &v, Real const &r);              // Comparison with scalar
  friend bool operator!=(Vector3D const &v, Real const &r);              // Comparison with scalar
  friend std::ostream& operator<<(std::ostream &str, Vector3D const &v); // Output

// Other operations

  friend Real const dot(Vector3D const &v1, Vector3D const &v2);        // Dot product
  friend Vector3D const cross(Vector3D const &v1, Vector3D const &v2);  // Cross product
  friend Real const norm(Vector3D const &v);                            // Norm
  friend Real const norm2(Vector3D const &v);                           // Norm squared
  friend Vector3D const unit(Vector3D const &v);                        // Unit vector in the direction of v
  friend Real const angle(Vector3D const &v1, Vector3D const &v2);      // Returns the angle (in radians)
                                                                        // between v1 and v2
  friend Real const dihedral(Vector3D const &v1, Vector3D const &v2,    // Returs the dihedral angle (in
                             Vector3D const &v3);                       // radians) defined by v1, v2 and v3
};

/*
** End of class Vector3D
*/

// Inlines

inline Vector3D::Vector3D() 
{
  x = y = z = 0.0;
}

inline Vector3D::Vector3D(Real const &r)
{
  x = y = z = r;
}

inline Vector3D::Vector3D(Real const &newx, Real const &newy, Real const &newz)
{
  x = newx; y = newy; z = newz;
}

inline Vector3D::Vector3D(Real const v[3])
{
  x = v[0]; y = v[1]; z = v[2];
}

inline Vector3D::Vector3D(SmallReal const v[3])
{
  x = (Real)v[0]; y = (Real)v[1]; z = (Real)v[2];
}

inline Vector3D &Vector3D::operator+=(Vector3D const &v)
{
  x += v.x; y += v.y; z += v.z;
  return *this;
}

inline Vector3D &Vector3D::operator-=(Vector3D const &v)
{
  x -= v.x; y -= v.y; z -= v.z;
  return *this;
}

inline Vector3D &Vector3D::operator*=(Real const &r)
{
  x *= r; y *= r; z*= r;
  return *this;
}

inline Vector3D &Vector3D::operator/=(Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  x *= invr; y *= invr; z *= invr;
  return *this;
}

inline Real &Vector3D::operator[](size_t const i)
{
  assert(i < 3);
  return (i == 0)?x:(i == 1)?y:z;
}

inline Real const &Vector3D::operator[](size_t const i) const
{
  assert(i < 3);
  return (i == 0)?x:(i == 1)?y:z;
}

inline Real &Vector3D::operator()(size_t const i)
{
  assert(i < 3);
  return (i == 0)?x:(i == 1)?y:z;
}

inline Real const &Vector3D::operator()(size_t const i) const
{
  assert(i < 3);
  return (i == 0)?x:(i == 1)?y:z;
}

inline size_t const Vector3D::size() const
{
  return 3;
}

inline void Vector3D::normalize()
{
  Real length = sqrt(x*x + y*y + z*z);
  assert(length);
  Real invLength = 1.0/length;
  x *= invLength; y *= invLength; z *= invLength;
}

inline void Vector3D::getArray(Real array[3])
{
  array[0] = x; array[1] = y; array[2] = z;
}

inline void Vector3D::getArray(SmallReal array[3])
{
  array[0] = (SmallReal)x;
  array[1] = (SmallReal)y;
  array[2] = (SmallReal)z;
}

inline Vector3D operator+(Vector3D const &v1, Vector3D const &v2)
{
  return Vector3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector3D operator-(Vector3D const &v)
{
  return Vector3D(-v.x, -v.y, -v.z);
}

inline Vector3D operator-(Vector3D const &v1, Vector3D const &v2)
{
  return Vector3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector3D operator*(Real const &r, Vector3D const &v)
{
  return Vector3D(r*v.x, r*v.y, r*v.z);
}

inline Vector3D operator*(Vector3D const &v, Real const &r)
{
  return Vector3D(r*v.x, r*v.y, r*v.z);
}

inline Vector3D operator/(Vector3D const &v, Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  return Vector3D(invr*v.x, invr*v.y, invr*v.z);
}

inline Real operator*(Vector3D const &v1, Vector3D const &v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline bool operator==(Vector3D const &v1, Vector3D const &v2)
{
  return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

inline bool operator!=(Vector3D const &v1, Vector3D const &v2)
{
  return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
}

inline bool operator==(Vector3D const &v, Real const &r)
{
  return v.x == r && v.y == r && v.z == r;
}

inline bool operator!=(Vector3D const &v, Real const &r)
{
  return v.x != r || v.y != r || v.z != r;
}

inline std::ostream& operator<<(std::ostream &str, Vector3D const &v)
{
  return str << v.x << " " << v.y << " " << v.z;
}

inline Real const dot(Vector3D const &v1, Vector3D const &v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline Vector3D const cross(Vector3D const &v1, Vector3D const &v2)
{
  return Vector3D(v1.y*v2.z - v1.z*v2.y,
                  v1.z*v2.x - v1.x*v2.z,
                  v1.x*v2.y - v1.y*v2.x);
}

inline Real const norm(Vector3D const &v)
{
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline Real const norm2(Vector3D const &v)
{
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

inline Vector3D const unit(Vector3D const &v)
{
  Real length = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  assert(length);
  Real invLength = 1.0/length;
  return Vector3D(invLength*v.x, invLength*v.y, invLength*v.z);
}

inline Real const angle(Vector3D const &v1, Vector3D const &v2)
{
  Real arg = v1*v2/(norm(v1)*norm(v2));
  if(arg > 1.0) arg = 1.0;  // Avoid roundoff problems
  if(arg < -1.0) arg = -1.0;
  return acos(arg);
}

inline Real const dihedral(Vector3D const &v1, Vector3D const &v2, Vector3D const &v3)
{
  Vector3D n1 = cross(v1, v2);
  Vector3D n2 = cross(v2, v3);
  Vector3D n3 = cross(v2, n1);
  Real normn2 = norm(n2);
  return -atan2(n3*n2/(normn2*norm(n3)), n1*n2/(normn2*norm(n1)));
}

#endif
