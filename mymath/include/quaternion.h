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
** Simple quaternion class
**
** Note that all data members (w, x, y, z) are public.
**
** The routines for conversion between quaternions and rotation matrices are based
** on the ones by Martin Baker at:
** http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/index.htm
*/

#ifndef H_QUATERNION
#define H_QUATERNION

#include <iostream>
#include <cmath>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/matrix3D.h"
#include "mymath/include/vector4D.h"

Real const QUAT_TOLERANCE = 1.0E-8; // Tolerance for zero components (see rotationMatrix below)

class Quaternion
{
public:

  Real  w, x, y, z;  // Quaternion components

// Constructors

  Quaternion();                     // Defines a null quaternion
  Quaternion(Real const &r);        // Defines a pure scalar quaternion (w = r, (x,y,z) = 0)
  Quaternion(Vector3D const &v);    // Defines a pure vector quaternion (w = 0, (x,y,z) = v)
  Quaternion(Real const &r,         // Defines a quaternion by giving the scalar and vector parts
             Vector3D const &v);
  Quaternion(Vector4D const &v);    // Defines a quaternion from a 4D vector
  Quaternion(Vector3D const &axis,  // Defines a unit quaternion equivalent to a rotation around
             Real const &angle);    // an axis (not necessarily normalized) by an angle (in radians)
  Quaternion(Matrix3D const &a);    // Defines a quaternion equivalent to the given rotation matrix
                                    // NOTE: This does *not* verify that the matrix is orthogonal!
  Quaternion(Real const &neww,      // Defines a quaternion by giving all the components
             Real const &newx,
             Real const &newy, 
             Real const &newz); 

// Assignment

  Quaternion &operator+=(Quaternion const &q);  // Assignment by addition
  Quaternion &operator+=(Real const &r);        // Assignment by addition of scalar
  Quaternion &operator+=(Vector3D const &v);    // Assignment by addition of vector
  Quaternion &operator-=(Quaternion const &q);  // Assignment by subtraction
  Quaternion &operator-=(Real const &r);        // Assignment by subtraction of scalar
  Quaternion &operator-=(Vector3D const &v);    // Assignment by subtraction of vector
  Quaternion &operator*=(Real const &r);        // Assignment by multiplication by scalar
  Quaternion &operator/=(Real const &r);        // Assignment by division by scalar

// Access elements by index

  Real &operator[](size_t const i);             
  Real const &operator[](size_t const i) const; 
  Real &operator()(size_t const i);             
  Real const &operator()(size_t const i) const; 

// Normalization

  void normalize();      // Normalizes the quaternion

// Operators

  friend Quaternion operator+(Quaternion const &q1, Quaternion const &q2);  // Addition
  friend Quaternion operator+(Quaternion const &q, Real const &r);          // Addition of scalar
  friend Quaternion operator+(Real const &r, Quaternion const &q);          // Addition of scalar
  friend Quaternion operator+(Quaternion const &q, Vector3D const &v);      // Addition of vector
  friend Quaternion operator+(Vector3D const &v, Quaternion const &q);      // Addition of vector
  friend Quaternion operator+(Real const &r, Vector3D const &v);            // Scalar plus vector
  friend Quaternion operator+(Vector3D const &v, Real const &r);            // Scalar plus vector
  friend Quaternion operator-(Quaternion const &q);                         // Negation
  friend Quaternion operator-(Quaternion const &q1, Quaternion const &q2);  // Subtraction
  friend Quaternion operator-(Quaternion const &q, Real const &r);          // Subtraction of scalar
  friend Quaternion operator-(Real const &r, Quaternion const &q);          // Subtraction of scalar
  friend Quaternion operator-(Quaternion const &q, Vector3D const &v);      // Subtraction of vector
  friend Quaternion operator-(Vector3D const &v, Quaternion const &q);      // Subtraction of vector
  friend Quaternion operator-(Real const &r, Vector3D const &v);            // Scalar minus vector
  friend Quaternion operator-(Vector3D const &v, Real const &r);            // Vector minus scalar
  friend Quaternion operator*(Real const &r, Quaternion const &q);          // Multiplication by scalar
  friend Quaternion operator*(Quaternion const &q, Real const &r);          // Multiplication by scalar
  friend Quaternion operator*(Quaternion const &q1, Quaternion const &q2);  // Quaternion product
  friend Quaternion operator*(Quaternion const &q, Vector3D const &v);      // Multiplication by vector
  friend Quaternion operator*(Vector3D const &v, Quaternion const &q);      // Multiplication by vector
  friend Quaternion operator/(Quaternion const &q, Real const &r);          // Division by scalar
  friend Quaternion operator/(Real const &r, Quaternion const &q);          // Division of scalar by quaternion
  friend Quaternion operator~(Quaternion const &q);                         // Quaternion conjugate
  friend Quaternion operator%(Quaternion const &q1, Quaternion const &q2);  // Quaternion rotating from q2 to q1 (= q1*~q2)
  friend bool operator==(Quaternion const &q1, Quaternion const &q2);       // Comparison (equal)
  friend bool operator!=(Quaternion const &q1, Quaternion const &q2);       // Comparison (not equal)
  friend bool operator==(Quaternion const &q, Real const &r);               // Comparison with scalar
  friend bool operator!=(Quaternion const &q, Real const &r);               // Comparison with scalar
  friend bool operator==(Quaternion const &q, Vector3D const &v);           // Comparison with vector
  friend bool operator!=(Quaternion const &q, Vector3D const &v);           // Comparison with vector
  friend std::ostream& operator<<(std::ostream &str, Quaternion const &q);  // Output

// Other operations

  friend Real const dot(Quaternion const &q1,                   // Dot product (same as scalar(q1%q2), but faster)
                        Quaternion const &q2);
  friend Real const norm(Quaternion const &q);                  // Norm
  friend Real const norm2(Quaternion const &q);                 // Norm squared
  friend Quaternion const conjugate(Quaternion const &q);       // Quaternion conjugate
  friend Real const scalar(Quaternion const &q);                // Scalar part of the quaternion
  friend Vector3D const vector(Quaternion const &q);            // Vector part of the quaternion
  friend Quaternion const axisAngle(Vector3D const &axis,       // Returns a quaternion representing a rotation around an
                                    Real const &angle);         // axis (not necessarily normalized) by an angle (in radians)
  friend Matrix3D const rotationMatrix(Quaternion const &q);    // Returns a 3D rotation matrix equivalent to the quaternion
  friend Quaternion const quaternion(Matrix3D const &a);        // Returns a normalized quaternion equivalent to a 3D rotation matrix
                                                                // Note: This does *not* verify that the matrix is orthogonal!
  friend void rotateVector(Quaternion const &q, Vector3D &v);   // Rotates the vector v by the given quaternion
  friend Vector3D const rotatedVector(Quaternion const &q,      // Returns the vector v rotated by the quaternion q
                                      Vector3D const &v);

};

/*
** End of class Quaternion
*/

// Inlines

inline Quaternion::Quaternion()
{
  w = x = y = z = 0.0;
}

inline Quaternion::Quaternion(Real const &r)
{
  w = r; x = y = z = 0.0;
}

inline Quaternion::Quaternion(Vector3D const &v)
{
  w = 0.0; x = v.x; y = v.y; z = v.z;
}

inline Quaternion::Quaternion(Real const &r, Vector3D const &v)
{
  w = r; x = v.x; y = v.y; z = v.z;
}

inline Quaternion::Quaternion(Vector4D const &v)
{
  w = v.w; x = v.x; y = v.y; z = v.z;
}

inline Quaternion::Quaternion(Vector3D const &axis, Real const &angle)
{
  Real fact = sin(angle/2.0)/norm(axis);
  w = cos(angle/2.0); x = fact*axis.x; y = fact*axis.y; z = fact*axis.z;
}
 
inline Quaternion::Quaternion(Matrix3D const &a)
{
  Real tr, s;
  tr = 1.0 + trace(a);
  if( tr > QUAT_TOLERANCE ) // Avoid roundoff error
  {
    s = 0.5/sqrt(tr);
    w = 0.25/s;
    x = s*(a(2,1) - a(1,2));
    y = s*(a(0,2) - a(2,0));
    z = s*(a(1,0) - a(0,1));
  } 
  else if( a(0,0) > a(1,1) && a(0,0) > a(2,2) ) 
  {
    s = 2.0*sqrt(1.0 + a(0,0) - a(1,1) - a(2,2));
    w = (a(2,1) - a(1,2))/s;
    x = 0.25*s;
    y = (a(0,1) + a(1,0))/s;
    z = (a(0,2) + a(2,0))/s;
  } 
  else if( a(1,1) > a(2,2) ) 
  {
    s = 2.0*sqrt(1.0 + a(1,1) - a(0,0) - a(2,2));
    w = (a(0,2) - a(2,0))/s;
    x = (a(0,1) + a(1,0))/s;
    y = 0.25*s;
    z = (a(1,2) + a(2,1))/s;
  }
  else 
  {
    s = 2.0*sqrt(1.0 + a(2,2) - a(0,0) - a(1,1));
    w = (a(1,0) - a(0,1))/s;
    x = (a(0,2) + a(2,0))/s;
    y = (a(1,2) + a(2,1))/s;
    z = 0.25*s;
  }
}

inline Quaternion::Quaternion(Real const &neww, Real const &newx, Real const &newy, Real const &newz)
{
  w = neww; x = newx; y = newy; z = newz;
}

// Assignment

inline Quaternion &Quaternion::operator+=(Quaternion const &q)
{
  w += q.w; x += q.x; y += q.y; z += q.z;
  return *this;
}

inline Quaternion &Quaternion::operator+=(Real const &r)
{
  w += r;
  return *this;
}

inline Quaternion &Quaternion::operator+=(Vector3D const &v)
{
  x += v.x; y += v.y; z += v.z;
  return *this;
}

inline Quaternion &Quaternion::operator-=(Quaternion const &q)
{
  w -= q.w; x -= q.x; y -= q.y; z -= q.z;
  return *this;
}

inline Quaternion &Quaternion::operator-=(Real const &r)
{
  w -= r;
  return *this;
}

inline Quaternion &Quaternion::operator-=(Vector3D const &v)
{
  x -= v.x; y -= v.y; z -= v.z;
  return *this;
}

inline Quaternion &Quaternion::operator*=(Real const &r)
{
  w *= r; x *= r; y *= r; z *= r;
  return *this;
}

inline Quaternion &Quaternion::operator/=(Real const &r)
{
  assert(r);
  Real invr = 1.0/r;
  w*= invr; x *= invr; y *= invr; z *= invr;
  return *this;
}

inline Real &Quaternion::operator[](size_t const i)
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline Real const &Quaternion::operator[](size_t const i) const
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline Real &Quaternion::operator()(size_t const i)
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline Real const &Quaternion::operator()(size_t const i) const
{
  assert(i < 4);
  return (i == 0)?w:(i == 1)?x:(i == 2)?y:z;
}

inline void Quaternion::normalize()
{
  Real length = sqrt(w*w + x*x + y*y + z*z);
  assert(length);
  Real invLength = 1.0/length;
  w*= invLength; x *= invLength; y *= invLength; z *= invLength;
}

inline Quaternion operator+(Quaternion const &q1, Quaternion const &q2)
{
  return Quaternion(q1.w + q2.w, q1.x + q2.x, q1.y + q2.y, q1.z + q2.z);
}

inline Quaternion operator+(Quaternion const &q, Real const &r)
{
  return Quaternion(q.w + r, q.x, q.y, q.z);
}

inline Quaternion operator+(Real const &r, Quaternion const &q)
{
  return Quaternion(q.w + r, q.x, q.y, q.z);
}

inline Quaternion operator+(Quaternion const &q, Vector3D const &v)
{
  return Quaternion(q.w, q.x + v.x, q.y + v.y, q.z + v.z);
}

inline Quaternion operator+(Vector3D const &v, Quaternion const &q)
{
  return Quaternion(q.w, q.x + v.x, q.y + v.y, q.z + v.z);
}

inline Quaternion operator+(Real const &r, Vector3D const &v)
{
  return Quaternion(r, v.x, v.y, v.z);
}

inline Quaternion operator+(Vector3D const &v, Real const &r)
{
  return Quaternion(r, v.x, v.y, v.z);
}

inline Quaternion operator-(Quaternion const &q)
{
  return Quaternion(-q.w, -q.x, -q.y, -q.z);
}

inline Quaternion operator-(Quaternion const &q1, Quaternion const &q2)
{
  return Quaternion(q1.w - q2.w, q1.x - q2.x, q1.y - q2.y, q1.z - q2.z);
}

inline Quaternion operator-(Quaternion const &q, Real const &r)
{
  return Quaternion(q.w - r, q.x, q.y, q.z);
}

inline Quaternion operator-(Real const &r, Quaternion const &q)
{
  return Quaternion(r - q.w, -q.x, -q.y, -q.z);
}

inline Quaternion operator-(Quaternion const &q, Vector3D const &v)
{
  return Quaternion(q.w, q.x - v.x, q.y - v.y, q.z - v.z);
}

inline Quaternion operator-(Vector3D const &v, Quaternion const &q)
{
  return Quaternion(-q.w, v.x - q.x, v.y - q.y, v.z - q.z);
}

inline Quaternion operator-(Real const &r, Vector3D const &v)
{
  return Quaternion(r, -v.x, -v.y, -v.z);
}

inline Quaternion operator-(Vector3D const &v, Real const &r)
{
  return Quaternion(-r, v.x, v.y, v.z);
}

inline Quaternion operator*(Real const &r, Quaternion const &q)
{
  return Quaternion(r*q.w, r*q.x, r*q.y, r*q.z);
}

inline Quaternion operator*(Quaternion const &q, Real const &r)
{
  return Quaternion(r*q.w, r*q.x, r*q.y, r*q.z);
}

inline Quaternion operator*(Quaternion const &q1, Quaternion const &q2)
{
  return Quaternion(q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z, 
                    q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y, 
                    q1.w*q2.y + q1.y*q2.w + q1.z*q2.x - q1.x*q2.z,
                    q1.w*q2.z + q1.z*q2.w + q1.x*q2.y - q1.y*q2.x);
}

inline Quaternion operator*(Quaternion const &q, Vector3D const &v)
{
  return Quaternion(-q.x*v.x - q.y*v.y - q.z*v.z, 
                     q.w*v.x + q.y*v.z - q.z*v.y, 
                     q.w*v.y + q.z*v.x - q.x*v.z,
                     q.w*v.z + q.x*v.y - q.y*v.x);
}

inline Quaternion operator*(Vector3D const &v, Quaternion const &q)
{
  return Quaternion(-v.x*q.x - v.y*q.y - v.z*q.z, 
                     v.x*q.w + v.y*q.z - v.z*q.y, 
                     v.y*q.w + v.z*q.x - v.x*q.z,
                     v.z*q.w + v.x*q.y - v.y*q.x);
}

inline Quaternion operator/(Quaternion const &q, Real const &r)
{
  assert(r);
  Real invr = 1.0/r;
  return Quaternion(invr*q.w, invr*q.x, invr*q.y, invr*q.z);
}

inline Quaternion operator/(Real const &r, Quaternion const &q)
{
  Real fact = r/(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
  return Quaternion(fact*q.w, -fact*q.x, -fact*q.y, -fact*q.z);
}

inline Quaternion operator~(Quaternion const &q)
{
  return Quaternion(q.w, -q.x, -q.y, -q.z);
}

inline Quaternion operator%(Quaternion const &q1, Quaternion const &q2)
{
  return Quaternion(q1.w*q2.w + q1.x*q2.x + q1.y*q2.y + q1.z*q2.z, 
                    q1.x*q2.w - q1.w*q2.x - q1.y*q2.z + q1.z*q2.y, 
                    q1.y*q2.w - q1.w*q2.y - q1.z*q2.x + q1.x*q2.z,
                    q1.z*q2.w - q1.w*q2.z - q1.x*q2.y + q1.y*q2.x);
}

inline bool operator==(Quaternion const &q1, Quaternion const &q2)
{
  return q1.w == q2.w && q1.x == q2.x && q1.y == q2.y && q1.z == q2.z;
}

inline bool operator!=(Quaternion const &q1, Quaternion const &q2)
{
  return q1.w != q2.w || q1.x != q2.x || q1.y != q2.y || q1.z != q2.z;
}

inline bool operator==(Quaternion const &q, Real const &r)
{
  return q.w == r && q.x == 0.0 && q.y == 0.0 && q.z == 0.0;
}

inline bool operator!=(Quaternion const &q, Real const &r)
{
  return q.w != r || q.x != 0.0 || q.y != 0.0 || q.z != 0.0;
}

inline bool operator==(Quaternion const &q, Vector3D const &v)
{
  return q.w == 0.0 && q.x == v.x && q.y == v.y && q.z == v.z;
}

inline bool operator!=(Quaternion const &q, Vector3D const &v)
{
  return q.w != 0.0 || q.x != v.x || q.y != v.y || q.z != v.z;
}

inline std::ostream& operator<<(std::ostream &str, Quaternion const &q)
{
  return str << q.w << " + i*" << q.x << " + j*" << q.y << " + k*" << q.z;
}

inline Real const dot(Quaternion const &q1, Quaternion const &q2)
{
  return q1.w*q2.w + q1.x*q2.x + q1.y*q2.y + q1.z*q2.z;
}

inline Real const norm(Quaternion const &q)
{
  return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
}

inline Real const norm2(Quaternion const &q)
{
  return q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z;
}

inline Quaternion const conjugate(Quaternion const &q)
{
  return Quaternion(q.w, -q.x, -q.y, -q.z);
}

inline Real const scalar(Quaternion const &q)
{
  return q.w;
}

inline Vector3D const vector(Quaternion const &q)
{
  return Vector3D(q.x, q.y, q.z);
}

inline Quaternion const axisAngle(Vector3D const &axis, Real const &angle)
{
  return Quaternion(axis, angle);
}

inline Matrix3D const rotationMatrix(Quaternion const &q)
{
  Real fact = 2.0/norm2(q);
  Matrix3D a;
  a(0,0) = 1.0 - fact*(q.y*q.y + q.z*q.z);
  a(0,1) = fact*(q.x*q.y - q.w*q.z);
  a(0,2) = fact*(q.x*q.z + q.w*q.y);
  a(1,0) = fact*(q.x*q.y + q.w*q.z);
  a(1,1) = 1.0 - fact*(q.x*q.x + q.z*q.z);
  a(1,2) = fact*(q.y*q.z - q.w*q.x);
  a(2,0) = fact*(q.x*q.z - q.w*q.y);
  a(2,1) = fact*(q.y*q.z + q.w*q.x);
  a(2,2) = 1.0 - fact*(q.x*q.x + q.y*q.y);
  return a;
}

inline Quaternion const quaternion(Matrix3D const &a)
{
  return Quaternion(a);
}

inline void rotateVector(Quaternion const &q, Vector3D &v)
{
  Real fact = 2.0/norm2(q);
  v = Vector3D((1.0 - fact*(q.y*q.y + q.z*q.z))*v.x + fact*(q.x*q.y - q.w*q.z)*v.y + fact*(q.x*q.z + q.w*q.y)*v.z,
               fact*(q.x*q.y + q.w*q.z)*v.x + (1.0 - fact*(q.x*q.x + q.z*q.z))*v.y + fact*(q.y*q.z - q.w*q.x)*v.z,
               fact*(q.x*q.z - q.w*q.y)*v.x + fact*(q.y*q.z + q.w*q.x)*v.y + (1.0 - fact*(q.x*q.x + q.y*q.y))*v.z);
}

inline Vector3D const rotatedVector(Quaternion const &q, Vector3D const &v)
{
  Real fact = 2.0/norm2(q);
  return Vector3D((1.0 - fact*(q.y*q.y + q.z*q.z))*v.x + fact*(q.x*q.y - q.w*q.z)*v.y + fact*(q.x*q.z + q.w*q.y)*v.z,
                  fact*(q.x*q.y + q.w*q.z)*v.x + (1.0 - fact*(q.x*q.x + q.z*q.z))*v.y + fact*(q.y*q.z - q.w*q.x)*v.z,
                  fact*(q.x*q.z - q.w*q.y)*v.x + fact*(q.y*q.z + q.w*q.x)*v.y + (1.0 - fact*(q.x*q.x + q.y*q.y))*v.z);
}
 
#endif
