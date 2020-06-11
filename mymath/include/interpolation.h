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
** Routines for interpolation/approximation - add as needed
**
** Many ideas for the b-spline routines in this code come from the
** open source spline library by John Burkardt. 
** See http://people.sc.fsu.edu/~burkardt/cpp_src/spline/spline.html
** Also see Carl deBoor, "A Practical Guide to Splines", Springer,
** New York (2001); ISBN 0-387-95366-3.
*/

#ifndef H_MYMATH_INTERPOLATION
#define H_MYMATH_INTERPOLATION

#include <utility>
#include <vector>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/vectorND.h"

// Given an ordered set of values and a number, return the index
// of the interval in the data containing the number.
// The routine returns 0 if value < data[0], N-1 if value > data[N-1],
// where N is the number of data points.
// Based on the code by John Burkardt.
// TypeData is a std::vector-like data type with a size() method
// and an operator[] that returns Real.

template<typename TypeData>
size_t const findBracket(TypeData const &data, Real const &value)
{
  size_t numPoints = data.size();
  for(size_t i = 2; i < numPoints; ++i)
    if(value < data[i - 1]) return i - 1;
  return numPoints - 1;
}

// Fit a b-spline to approximate the given (N-dimensional)
// trajectory, and evaluate it at the given value.
// Based on the code by John Burkardt.
// Template types follow the same convention as findBracket, but:
// TypeData must also allow basic vector arithmetic (addition,
// subtraction, multiplication by scalar).

template<typename TypeData>
TypeData const evaluateBSpline(std::vector<Real> const &parameter, 
                               std::vector<TypeData> const &point,
                               Real const &value)
{
  assert(parameter.size() == point.size());
  assert(parameter.size() > 1);
  size_t numPoints = point.size();
  size_t left = 
    findBracket<std::vector<Real> >(parameter, value);
  size_t right = left + 1;
  
  Real scaled = (value - parameter[left - 1])/
                       (parameter[right - 1] - parameter[left - 1]);
  TypeData result;
  Real bFunction; // Basis function
  bFunction = (((-scaled + 3.0)*scaled - 3.0)*scaled + 1.0)/6.0;
  if(left > 1)
    result = point[left - 2]*bFunction;
  else
    result = (2.0*point[0] - point[1])*bFunction;
  bFunction = ((3.0*scaled - 6.0)*scaled*scaled + 4.0)/6.0;
    result += point[left - 1]*bFunction;
  bFunction = (((-3.0*scaled + 3.0)*scaled + 3.0)*scaled + 1.0)/6.0;
    result += point[right - 1]*bFunction;
  bFunction = scaled*scaled*scaled/6.0;
  if(right < numPoints)
    result += point[right]*bFunction;
  else
    result += (2.0*point[numPoints - 1] - point[numPoints - 2])*bFunction;
  return result;
}

// Piecewise linear interpolation for a given trajectory.
// Template types follow the same convention as evaluateBSpline. 

template<typename TypeData>
TypeData const evaluateLinearSpline(std::vector<Real> const &parameter, 
                                    std::vector<TypeData> const &point,
                                    Real const &value)
{
  assert(parameter.size() == point.size());
  assert(parameter.size() > 1);
  size_t left = 
    findBracket<std::vector<Real> >(parameter, value);
  size_t right = left + 1;
  
  Real scaled = (value - parameter[left - 1])/
                (parameter[right - 1] - parameter[left - 1]);
  TypeData result = point[left - 1] + 
                    scaled*(point[right - 1] - point[left - 1]);
  return result;
}

#endif

