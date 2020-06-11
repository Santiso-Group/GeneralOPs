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
** Optimization routines - add as needed
**
** For an example on how to use BFGS, see the file "optimization_example.h"
*/

#ifndef H_OPTIMIZATION
#define H_OPTIMIZATION

#include "common/include/types.h"
#include "mymath/include/vectorND.h"
#include "mymath/include/matrixND.h"

// Some global constants for the optimization routines

size_t const DEFAULT_LINESEARCH_MAX_ITERATIONS = 100;
size_t const DEFAULT_BFGS_MAX_ITERATIONS = 1000;
Real   const DEFAULT_OPTIMIZATION_X_TOLERANCE = 2.0*REAL_EPSILON;
Real   const DEFAULT_GRADIENT_TOLERANCE = 2.0*REAL_EPSILON;
Real   const DEFAULT_LINESEARCH_MAX_LENGTH = 10.0;
Real   const SQRT_EPSILON = sqrt(REAL_EPSILON);
Real   const WOLFE_C1 = 1.E-4;
Real   const WOLFE_C2_BFGS = 0.9; // Not used yet

// Interpolating line search: Used by BFGS

template <typename Type>
VectorND const lineSearch(Type *myClass,
                          Real const (Type::*target)(VectorND const &),
                          VectorND const &gradient,
                          VectorND const &xInitial,
                          Real &newf,
                          VectorND &searchDirection,
                          Real const &xTolerance = DEFAULT_OPTIMIZATION_X_TOLERANCE,
                          size_t const maxIterations = DEFAULT_LINESEARCH_MAX_ITERATIONS,
                          Real const &maxSearchLength = DEFAULT_LINESEARCH_MAX_LENGTH)
{
  // Initialize variables
  size_t numIterations = 0;
  Real olderLength = 1.0, oldLength = 1.0, newLength = 1.0;
  Real initialf = newf, oldf = newf;
  VectorND newx = xInitial;
  // Reduce search length if necessary
  Real searchDirectionLength = norm(searchDirection);
  if(searchDirectionLength > maxSearchLength) 
    searchDirection *= maxSearchLength/searchDirectionLength;
  Real directionalDerivative = gradient*searchDirection;
  while(numIterations <= maxIterations)
  {
    // Do housekeeping and find new value of target function
    oldf = newf;
    olderLength = oldLength; oldLength = newLength;
    newx = xInitial + newLength*searchDirection;
    newf = (myClass->*target)(newx);
    // Check for x convergence
    if(normInfinity(newx - xInitial) < xTolerance) 
      return xInitial;
    // Armijo rule
    if(newf <= initialf + WOLFE_C1*newLength*directionalDerivative)
      return newx;  // Done
    else
    {
      if(newLength == 1.0)  // First step: use quadratic interpolation
        newLength = -0.5*directionalDerivative/(newf - oldf - directionalDerivative);
      else  // Use cubic interpolation
      {
        Real oldTerm = newf - initialf - directionalDerivative*oldLength;
        Real olderTerm = oldf - initialf - directionalDerivative*olderLength;
        Real den1 = 1.0/(oldLength*oldLength*(oldLength - olderLength));
        Real den2 = 1.0/(olderLength*olderLength*(oldLength - olderLength));
        Real aCubic = den1*oldTerm - den2*olderTerm;
        Real bCubic = -olderLength*den1*oldTerm + oldLength*den2*olderTerm;
        Real disc = bCubic*bCubic - 3.0*aCubic*directionalDerivative;
        if(disc > 0.0) newLength = (-bCubic + sqrt(disc))/(3.0*aCubic);
        else newLength = 0.5*oldLength;
        if(newLength > 0.5*oldLength) newLength = 0.5*oldLength;
      }
      if(newLength < 0.1*oldLength) newLength = 0.1*oldLength;
    }
    ++numIterations;
  }
  std::cerr << "Warning from lineSearch: Maximum number of iterations reached" << std::endl;
  return newx;
}

// BFGS minimization

template <typename Type>
VectorND const bfgs(Type *myClass, 
                    Real const (Type::*target)(VectorND const &),
                    VectorND const (Type::*gradient)(VectorND const &), 
                    VectorND const &xInitial,
                    Real const &xTolerance = DEFAULT_OPTIMIZATION_X_TOLERANCE,
                    Real const &gradientTolerance = DEFAULT_GRADIENT_TOLERANCE,
                    size_t const maxIterations = DEFAULT_BFGS_MAX_ITERATIONS,
                    size_t const maxLinesearchIterations = DEFAULT_LINESEARCH_MAX_ITERATIONS,
                    Real const &maxLinesearchLength = DEFAULT_LINESEARCH_MAX_LENGTH)
{
  // Initialize variables
  VectorND newx, oldx;
  Real newtarget;
  VectorND newGradient, oldGradient;
  size_t numIterations = 0;
  size_t const numDimensions = xInitial.size();
  MatrixND const identity = MatrixND::identity(numDimensions);
  MatrixND inverseHessian = identity;
  newx = xInitial;
  newtarget = (myClass->*target)(xInitial);
  newGradient = (myClass->*gradient)(xInitial);
  // Check for consistency
  if(newGradient.size() != numDimensions)
  {
    std::cerr << "Error in BFGS: Inconsistent initial point and target gradient" << std::endl;
    return xInitial;
  }
  // Main loop
  while(numIterations < maxIterations)
  {
    // Do housekeeping and find search direction
    oldx = newx; oldGradient = newGradient;
    VectorND searchDirection = -inverseHessian*newGradient;
    // Line search
    newx = lineSearch(myClass, target, newGradient, newx, newtarget,
                      searchDirection, xTolerance, maxLinesearchIterations,
                      maxLinesearchLength);
    VectorND deltax = newx - oldx;
    // Test for x convergence
    if(normInfinity(deltax) < xTolerance) return newx;
    newGradient = (myClass->*gradient)(newx);
    // Test for gradient convergence
    if(norm(newGradient) < gradientTolerance) return newx;
    VectorND deltaGradient = newGradient - oldGradient;
    Real deltaDot = deltaGradient*deltax;
    if(deltaDot > SQRT_EPSILON) // Avoid overflow in inverse Hessian
    {
      Real inverseDot = 1.0/deltaDot;
      // The BFGS update
      inverseHessian = (identity - inverseDot*outer(deltax, deltaGradient))*
                       inverseHessian*
                       (identity - inverseDot*outer(deltaGradient, deltax)) +
                       inverseDot*outer(deltax, deltax);
    }
    ++numIterations;
  }
  std::cerr << "Warning from BFGS: Reached maximum number of iterations" << std::endl;
  return newx;
}

#endif
