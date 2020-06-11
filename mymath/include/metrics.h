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
** Routines for calculation of various metrics
*
** All algorithms are implemented as templates. The specific data types
** are explained for each algorithm.
*/

#ifndef H_MYMATH_METRICS
#define H_MYMATH_METRICS

#include <vector>
#include <algorithm>
#include "common/include/types.h"
#include "common/include/assert.h"

// Auxiliary function to calculate distance for template objects.
// TypePoint must support one level of operator[], with the
// element being TypeElement, and a size() method.
// TypeElement should allow for basic arithmetic

template <typename TypePoint, typename TypeElement>
TypeElement const distance(TypePoint const &point1, TypePoint const &point2)
{
  size_t const numElements = point1.size();
  if(point2.size() != numElements)
  {
    std::cerr << "Error in metrics::distance - inconsistent number of dimensions" << std::endl;
    return 0.0;
  }
  TypeElement distanceSq = 0.0;
  for(size_t i = 0; i < numElements; ++i)
  {
    TypeElement diff = point1[i] - point2[i];
    distanceSq += diff*diff;
  }
  return sqrt(distanceSq);
}

// Auxiliary function for the calculation of Frechet distance.
// See T. Eiter and H. Mannila, "Computing Discrete Frechet Distance",
// Technical Report CD-TR 94/64, Christian Doppler Laboratory for
// Expert Systems, TU Vienna, Austria, 1994.
// TypeTrajectory and TypeElement as defined below in frechetDistance.

template <typename TypeTrajectory, typename TypePoint, typename TypeElement>
TypeElement const frechetAux(size_t const i, size_t const j,
                             std::vector<std::vector<TypeElement> > &frechetMat,
                             TypeTrajectory const &line1,
                             TypeTrajectory const &line2)
{
  if(frechetMat[i][j] >= 0.0) return frechetMat[i][j];
  else if(i == 0 && j == 0) frechetMat[i][j] = distance<TypePoint, TypeElement>(line1[0], line2[0]);
  else if(i > 0 && j == 0) 
  {
    Real dist = distance<TypePoint, TypeElement>(line1[i], line2[0]); 
    Real fAux = frechetAux<TypeTrajectory, TypePoint, TypeElement>(i - 1, 0, frechetMat, line1, line2);
    frechetMat[i][j] = std::max(dist, fAux);
  }
  else if(i == 0 && j > 0)
  {
    Real dist = distance<TypePoint, TypeElement>(line1[0], line2[j]); 
    Real fAux = frechetAux<TypeTrajectory, TypePoint, TypeElement>(0, j - 1, frechetMat, line1, line2);
    frechetMat[i][j] = std::max(dist, fAux);
  }
  else if(i > 0 && j > 0)
  {
    Real dist = distance<TypePoint, TypeElement>(line1[i], line2[j]);
    Real fAux1 = frechetAux<TypeTrajectory, TypePoint, TypeElement>(i - 1, j - 1, frechetMat, line1, line2);
    Real fAux2 = frechetAux<TypeTrajectory, TypePoint, TypeElement>(i - 1, j, frechetMat, line1, line2);
    Real fAux3 = frechetAux<TypeTrajectory, TypePoint, TypeElement>(i, j - 1, frechetMat, line1, line2);
    frechetMat[i][j] = std::max(dist, std::min(fAux1, std::min(fAux2, fAux3)));
  }
  else frechetMat[i][j] = REAL_VERYBIG;
  return frechetMat[i][j];
}

// Frechet distance between two polygonal curves.
// See ref. by Eiter and Mannila above.
// Also see A. Mascret et al., "Coastline Matching Process Based on the
// Discrete Frechet Distance", in "Progress in Spatial Data Handling",
// Part 7, Springer (2006), ISBN 978-3-540-35588-5 (Print) 978-3-540-35589-2 (Online)
// DOI 10.1007/3-540-35589-8_25, Pages 383-400
//
// TypeTrajectory must support two levels of operator[],
// with each supporting a size() method 
// (e.g. std::vector<std::vector<Real> >).
// Both arguments contain list of points in space
// describing the trajectories.
// TypePoint is the type of the elements of TypeTrajectory,
// and TypeElement is the type of the elements of TypePoint.
// TypeElement should allow for basic arithmetic.

template <typename TypeTrajectory, typename TypePoint, typename TypeElement>
TypeElement const frechetDistance(TypeTrajectory const &line1,
                                  TypeTrajectory const &line2)
{
  // Store numbers of points
  size_t const numPoints1 = line1.size();
  size_t const numPoints2 = line2.size();
  std::vector<std::vector<TypeElement> > frechetMat(numPoints1, std::vector<TypeElement>(numPoints2, -1));
  return frechetAux<TypeTrajectory, TypePoint, TypeElement>(numPoints1 - 1, numPoints2 - 1, frechetMat, line1, line2);
}

#endif

