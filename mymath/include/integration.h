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
** Routines for numerical integration
*
** All algorithms are implemented as templates. The specific data types
** are explained for each algorithm.
*/

#ifndef H_MYMATH_INTEGRATION
#define H_MYMATH_INTEGRATION

#include <vector>
#include "common/include/types.h"
#include "common/include/assert.h"

// Line integral using trapezoid rule 
// Both TypeTrajectory and TypeFields must support two levels of operator[],
// with each supporting a size() method (e.g. std::vector<std::vector<Real> >).
// and must contain the same number of elements.
// The first argument is a list of points in space
// describing the trajectory, the second argument is a set
// of "vectors" describing the fiels.
// TypeElement should allow for basic arithmetic.

template <typename TypeTrajectory, typename TypeFields, typename TypeElement>
TypeElement const trapezoidLineIntegral(TypeTrajectory const &line,
                                        TypeFields const &field)
{
  // Store number of points
  size_t numPoints = line.size();
  if(field.size() != numPoints)
  {
    std::cerr << "Error in trapezoidLineIntegral - inconsistent line and field sets" << std::endl;
    return 0;
  }
  if(numPoints == 0) return 0;
  size_t numDimensions = line[0].size();
  if(field[0].size() != numDimensions)
  {
    std::cerr << "Error in trapezoidLineIntegral - inconsistent line and field sets" << std::endl;
    return 0;
  }
  TypeElement integral = 0;
  for(size_t i = 0; i < numPoints - 1; ++i)
  {
    TypeElement deltaIntegral = 0;
    for(size_t j = 0; j < numDimensions; ++j)
      deltaIntegral += 0.5*(field[i][j] + field[i + 1][j])*(line[i + 1][j] - line[i][j]);
    integral += deltaIntegral;
  }
  return integral;
}

// Arc length using trapezoid rule
// Same convention as trapezoidLineIntegral, returns the
// arc length for the curve defined by the given set
// of points.

template <typename TypeTrajectory, typename TypeElement>
TypeElement trapezoidArcLength(TypeTrajectory const &line)
{
  // Store number of points
  size_t numPoints = line.size();
  if(numPoints == 0) return 0;
  size_t numDimensions = line[0].size();
  TypeElement arcLength = 0;
  for(size_t i = 0; i < numPoints - 1; ++i)
  {
    TypeElement deltaLength2 = 0;
    for(size_t j = 0; j < numDimensions; ++j)
    {
      TypeElement delta = (line[i + 1][j] - line[i][j]);
      deltaLength2 += delta*delta;
    }
    arcLength += sqrt(deltaLength2);
  }
  return arcLength;
}

// Cumulative integral using trapezoid rule 
// This is the same as trapezoidLineIntegral, but
// it returns a vector of values containing the
// cumulative integral up to each point.

template <typename TypeTrajectory, typename TypeFields, typename TypeElement>
std::vector<TypeElement> const trapezoidCumulativeIntegral(TypeTrajectory const &line,
                                                           TypeFields const &field)
{
  std::vector<TypeElement> cumulativeIntegral(1, 0);
  // Store number of points
  size_t numPoints = line.size();
  if(field.size() != numPoints)
  {
    std::cerr << "Error in trapezoidLineIntegral - inconsistent line and field sets" << std::endl;
    return cumulativeIntegral; 
  }
  if(numPoints == 0) return cumulativeIntegral;
  size_t numDimensions = line[0].size();
  if(field[0].size() != numDimensions)
  {
    std::cerr << "Error in trapezoidLineIntegral - inconsistent line and field sets" << std::endl;
    return cumulativeIntegral; 
  }
  
  for(size_t i = 0; i < numPoints - 1; ++i)
  {
    TypeElement deltaIntegral = 0;
    TypeElement &currentValue = cumulativeIntegral[i];
    for(size_t j = 0; j < numDimensions; ++j)
      deltaIntegral += 0.5*(field[i][j] + field[i + 1][j])*(line[i + 1][j] - line[i][j]);
    cumulativeIntegral.push_back(currentValue + deltaIntegral);
  }
   return cumulativeIntegral;
}

// Cumulative arc length using trapezoid rule
// Same convention as trapezoidCumulativeIntegral

template <typename TypeTrajectory, typename TypeElement>
std::vector<TypeElement> trapezoidCumulativeArcLength(TypeTrajectory const &line)
{
  std::vector<TypeElement> cumulativeLength(1, 0);
  // Store number of points
  size_t numPoints = line.size();
  if(numPoints == 0) return cumulativeLength; 
  size_t numDimensions = line[0].size();
  for(size_t i = 0; i < numPoints - 1; ++i)
  {
    TypeElement deltaLength2 = 0;
    TypeElement &currentValue = cumulativeLength[i];
    for(size_t j = 0; j < numDimensions; ++j)
    {
      TypeElement delta = (line[i + 1][j] - line[i][j]);
      deltaLength2 += delta*delta;
    }
    cumulativeLength.push_back(currentValue + sqrt(deltaLength2));
  }
  return cumulativeLength; 
}

#endif

