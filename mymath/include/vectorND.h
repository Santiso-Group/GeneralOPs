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
** Simple N-dimensional vector class
*/

/*
** Things that may be good to add:
** 
** - Compatibility with Vector3D and Vector4D (that is, ability to add
**   VectorND + Vector3D if N = 3, etc.)
*/

#ifndef H_VECTORND
#define H_VECTORND

#include <iostream>
#include <cmath>
#include <vector>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/vector4D.h"

class VectorND
{
public:

// Constructors

  VectorND();                                     // Defines an empty vector
  explicit VectorND(size_t const numDimensions);  // Defines a null vector with given 
                                                  // number of dimensions
  VectorND(size_t const numDimensions,            // Defines a vector with given number of
           Real const &r);                        // dimensions and all components equal to r
  VectorND(Vector3D const &v);                    // Defines a VectorND from a Vector3D
  VectorND(Vector4D const &v);                    // Defines a VectorND from a Vector4D
  VectorND(std::vector<Real> const &v);           // Defines a VectorND from a vector of reals

// Assignment

  VectorND &operator=(Real const &r);       // Make all components equal to a number
  VectorND &operator+=(VectorND const &v);  // Assignment by addition
  VectorND &operator-=(VectorND const &v);  // Assignment by subtraction
  VectorND &operator*=(Real const &r);      // Assignment by multiplication by scalar
  VectorND &operator/=(Real const &r);      // Assignment by division by scalar

// Access elements by index
  
  Real &operator[](size_t const i);             // Can modify
  Real const &operator[](size_t const i) const; // Cannot modify
  Real &operator()(size_t const i);             // Can modify
  Real const &operator()(size_t const i) const; // Cannot modify

// Clear, add components, change number of dimensions (similar to std::vector)

  void clear();                             // Clears the vector
  void push_back(Real const &r);            // Add a component to the vector
  void resize(size_t const numDimensions,   // Set the number of dimensions
              Real const &r = 0.0);         // and make extra components equal to r

// Size

  size_t const size() const;  // Returns the number of dimensions

// Normalization

  void normalize();      // Normalizes the vector

// Operators

  friend VectorND operator+(VectorND const &v1, VectorND const &v2);     // Addition
  friend VectorND operator-(VectorND const &v);                          // Negation
  friend VectorND operator-(VectorND const &v1, VectorND const &v2);     // Subtraction
  friend VectorND operator*(Real const &r, VectorND const &v);           // Multiplication by scalar
  friend VectorND operator*(VectorND const &v, Real const &r);           // Multiplication by scalar
  friend VectorND operator/(VectorND const &v, Real const &r);           // Division by scalar
  friend Real operator*(VectorND const &v1, VectorND const &v2);         // Dot product
  friend bool operator==(VectorND const &v1, VectorND const &v2);        // Comparison (equal)
  friend bool operator!=(VectorND const &v1, VectorND const &v2);        // Comparison (not equal)
  friend bool operator==(VectorND const &v, Real const &r);              // Comparison with scalar
  friend bool operator!=(VectorND const &v, Real const &r);              // Comparison with scalar
  friend std::ostream& operator<<(std::ostream &str, VectorND const &v); // Output
  friend std::istream& operator>>(std::istream &str, VectorND &v);       // Input

// Other operations

  friend Real const dot(VectorND const &v1, VectorND const &v2);  // Dot product
  friend Real const norm(VectorND const &v);                      // Euclidean norm
  friend Real const norm2(VectorND const &v);                     // Euclidean norm squared
  friend Real const normInfinity(VectorND const &v);              // Infinity (maximum) norm
  friend VectorND const unit(VectorND const &v);                  // Unit vector in the direction of v

private:

  std::vector<Real> vectorND_;
};

/*
** End of class VectorND
*/

// Inlines

inline VectorND::VectorND() 
:
vectorND_()
{}

inline VectorND::VectorND(size_t const numDimensions)
:
vectorND_(numDimensions, 0.0)
{}

inline VectorND::VectorND(size_t const numDimensions, Real const &r)
:
vectorND_(numDimensions, r)
{}

inline VectorND::VectorND(Vector3D const &v)
:
vectorND_()
{
  for(size_t i = 0; i < 3; i++) vectorND_.push_back(v[i]);
}

inline VectorND::VectorND(Vector4D const &v)
:
vectorND_()
{
  for(size_t i = 0; i < 4; i++) vectorND_.push_back(v[i]);
}

inline VectorND::VectorND(std::vector<Real> const &v)
:
vectorND_(v)
{}

inline VectorND &VectorND::operator=(Real const &r)
{
  for(size_t i = 0; i < vectorND_.size(); i++)
    vectorND_[i] = r;
  return *this;
}

inline VectorND &VectorND::operator+=(VectorND const &v)
{
  size_t const numDim = size();
  assert(numDim == v.size());
  for(size_t i = 0; i < numDim; i++)
    vectorND_[i] += v.vectorND_[i];
  return *this;
}

inline VectorND &VectorND::operator-=(VectorND const &v)
{
  size_t const numDim = size();
  assert(numDim == v.size());
  for(size_t i = 0; i < numDim; i++)
    vectorND_[i] -= v.vectorND_[i];
  return *this;
}

inline VectorND &VectorND::operator*=(Real const &r)
{
  for(size_t i = 0; i < vectorND_.size(); i++)
    vectorND_[i] *= r;
  return *this;
}

inline VectorND &VectorND::operator/=(Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  for(size_t i = 0; i < vectorND_.size(); i++)
    vectorND_[i] *= invr;
  return *this;
}

inline Real &VectorND::operator[](size_t const i)
{
  assert(i < vectorND_.size());
  return vectorND_[i];
}

inline Real const &VectorND::operator[](size_t const i) const
{
  assert(i < vectorND_.size());
  return vectorND_[i];
}

inline Real &VectorND::operator()(size_t const i)
{
  assert(i < vectorND_.size());
  return vectorND_[i];
}

inline Real const &VectorND::operator()(size_t const i) const
{
  assert(i < vectorND_.size());
  return vectorND_[i];
}

inline void VectorND::clear()
{
  vectorND_.clear();
}

inline void VectorND::push_back(Real const &r)
{
  vectorND_.push_back(r);
}

inline void VectorND::resize(size_t const numDimensions, Real const &r)
{
  vectorND_.resize(numDimensions, r);
}

inline size_t const VectorND::size() const
{
  return vectorND_.size();
}

inline void VectorND::normalize()
{
  Real length = 0.0;
  size_t const numDim = vectorND_.size();
  for(size_t i = 0; i < numDim; i++)
    length += vectorND_[i]*vectorND_[i];
  length = sqrt(length);
  assert(length);
  Real invLength = 1.0/length;
  for(size_t i = 0; i < numDim; i++)
    vectorND_[i] *= invLength;
}

inline VectorND operator+(VectorND const &v1, VectorND const &v2)
{
  size_t const numDim = v1.size();
  assert(numDim == v2.size());
  VectorND sum;
  for(size_t i = 0; i < numDim; i++)
    sum.vectorND_.push_back(v1[i] + v2[i]);
  return sum;
}

inline VectorND operator-(VectorND const &v)
{
  size_t const numDim = v.size();
  VectorND neg;
  for(size_t i = 0; i < numDim; i++)
    neg.vectorND_.push_back(-v[i]);
  return neg;
}

inline VectorND operator-(VectorND const &v1, VectorND const &v2)
{
  size_t const numDim = v1.size();
  assert(numDim == v2.size());
  VectorND diff;
  for(size_t i = 0; i < numDim; i++)
    diff.vectorND_.push_back(v1[i] - v2[i]);
  return diff;
}

inline VectorND operator*(Real const &r, VectorND const &v)
{
  size_t const numDim = v.size();
  VectorND prod;
  for(size_t i = 0; i < numDim; i++)
    prod.vectorND_.push_back(r*v[i]);
  return prod;
}

inline VectorND operator*(VectorND const &v, Real const &r)
{
  size_t const numDim = v.size();
  VectorND prod;
  for(size_t i = 0; i < numDim; i++)
    prod.vectorND_.push_back(r*v[i]);
  return prod;
}

inline VectorND operator/(VectorND const &v, Real const &r)
{
  assert(r);
  Real const invr = 1.0/r;
  size_t const numDim = v.size();
  VectorND quot;
  for(size_t i = 0; i < numDim; i++)
    quot.vectorND_.push_back(invr*v[i]);
  return quot;
}

inline Real operator*(VectorND const &v1, VectorND const &v2)
{
  size_t const numDim = v1.size();
  assert(numDim == v2.size());
  Real dot = 0.0;
  for(size_t i = 0; i < numDim; i++)
    dot += v1[i]*v2[i];
  return dot;
}

inline bool operator==(VectorND const &v1, VectorND const &v2)
{
  size_t const numDim = v1.size();
  if(numDim != v2.size()) return false;
  for(size_t i = 0; i < numDim; i++)
    if(v1[i] != v2[i]) return false;
  return true;
}

inline bool operator!=(VectorND const &v1, VectorND const &v2)
{
  size_t const numDim = v1.size();
  if(numDim != v2.size()) return true;
  for(size_t i = 0; i < numDim; i++)
    if(v1[i] != v2[i]) return true;
  return false;
}

inline bool operator==(VectorND const &v, Real const &r)
{
  for(size_t i = 0; i < v.size(); i++)
    if(v[i] != r) return false;
  return true;
}

inline bool operator!=(VectorND const &v, Real const &r)
{
  for(size_t i = 0; i < v.size(); i++)
    if(v[i] != r) return true;
  return false;
}

inline std::ostream& operator<<(std::ostream &str, VectorND const &v)
{
  for(size_t i = 0; i < v.size(); i++)
    str << v[i] << std::endl;
  return str;
}

inline std::istream& operator>>(std::istream &str, VectorND &v)
{
  v.clear();
  while(!str.eof())
  {
    Real v_i;
    str >> v_i;
    if(str.eof()) break;
    v.push_back(v_i);
  }
  if(!str && !str.eof())
  {
    std::cerr << "Error reading vector from stream" << std::endl;
    v.clear();
  }
  return str;
}

inline Real const dot(VectorND const &v1, VectorND const &v2)
{
  size_t const numDim = v1.size();
  assert(numDim == v2.size());
  Real dot = 0.0;
  for(size_t i = 0; i < numDim; i++)
    dot += v1[i]*v2[i];
  return dot;
}

inline Real const norm(VectorND const &v)
{
  Real norm2 = 0.0;
  for(size_t i = 0; i < v.size(); i++)
    norm2 += v[i]*v[i];
  return sqrt(norm2);
}

inline Real const norm2(VectorND const &v)
{
  Real norm2 = 0.0;
  for(size_t i = 0; i < v.size(); i++)
    norm2 += v[i]*v[i];
  return norm2;
}

inline Real const normInfinity(VectorND const &v)
{
  Real normInfinity = 0.0;
  for(size_t i = 0; i < v.size(); i++)
  {
    Real absv = fabs(v[i]);
    if(absv > normInfinity) normInfinity = absv;
  }
  return normInfinity;
}

inline VectorND const unit(VectorND const &v)
{
  size_t const numDim = v.size();
  Real length = 0.0;
  for(size_t i = 0; i < numDim; i++)
    length += v[i]*v[i];
  length = sqrt(length);
  assert(length);
  Real invLength = 1.0/length;
  VectorND unit;
  for(size_t i = 0; i < numDim; i++)
    unit.vectorND_.push_back(invLength*v[i]);
  return unit;
}

#endif
