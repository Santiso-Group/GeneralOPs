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
** Root finding routines - add as needed
**
** For examples, see "rootfinding_examples.h"
*/

#ifndef H_ROOT_FINDING
#define H_ROOT_FINDING

#include <iostream>
#include "common/include/types.h"

// Some global constants for the root finding routines

size_t const DEFAULT_NEWTON_MAX_ITERATIONS = 200;     // Also used for secant method
size_t const DEFAULT_BISECTION_MAX_ITERATIONS = 2000;
Real   const DEFAULT_ROOT_X_TOLERANCE = 2.0*REAL_EPSILON;
Real   const DEFAULT_ROOT_TARGET_TOLERANCE = 2.0*REAL_EPSILON;

// Single-variable Newton-Raphson method

template <typename Type>
Real const newtonRaphson(Type *myClass, 
                         Real const (Type::*target)(Real const &),
                         Real const (Type::*derivative)(Real const &), 
                         Real const &xInitial,
                         Real const &xTolerance = DEFAULT_ROOT_X_TOLERANCE,
                         Real const &targetTolerance = DEFAULT_ROOT_TARGET_TOLERANCE,
                         size_t const maxIterations = DEFAULT_NEWTON_MAX_ITERATIONS)
{
  // Initialize variables
  size_t numIterations = 0;
  Real newx = xInitial;
  // Main loop
  while(numIterations < maxIterations)
  {
    // Store old x
    Real oldx = newx; 
    // Find new x
    Real targetFunction = (myClass->*target)(newx);
    Real targetDerivative = (myClass->*derivative)(newx);
    if(fabs(targetDerivative) < targetTolerance)
    {
      std::cerr << "Zero derivative in newtonRaphson" << std::endl;
      return xInitial;
    }
    Real deltax = targetFunction/targetDerivative;
    newx = oldx - deltax;
    // Test for x convergence
    if(fabs(deltax) < xTolerance) return newx;
    // Test for target function convergence
    if(fabs(targetFunction) < targetTolerance) return newx;
    ++numIterations;
  }
  std::cerr << "Warning from newtonRaphson: Reached maximum number of iterations" << std::endl;
  return newx;
}

// Secant method

template <typename Type>
Real const secant(Type *myClass, 
                  Real const (Type::*target)(Real const &),
                  Real const &xInitial1, Real const &xInitial2,
                  Real const &xTolerance = DEFAULT_ROOT_X_TOLERANCE,
                  Real const &targetTolerance = DEFAULT_ROOT_TARGET_TOLERANCE,
                  size_t const maxIterations = DEFAULT_NEWTON_MAX_ITERATIONS)
{
  // Initialize variables
  size_t numIterations = 0;
  Real newx = xInitial1;
  Real x0 = xInitial1;
  Real x1 = xInitial2;
  Real f0 = (myClass->*target)(x0);
  Real f1 = (myClass->*target)(x1);
  // Main loop
  while(numIterations < maxIterations)
  {
    if(fabs(f1 - f0) < targetTolerance)
    {
      std::cerr << "Horizontal secant line in secant" << std::endl;
      return xInitial1;
    }
    // Find new x
    newx = x0 - f0*(x1 - x0)/(f1 - f0);
    Real newf = (myClass->*target)(newx);
    // Test for target function convergence
    if(fabs(newf) < targetTolerance) return newx;
    // Housekeeping
    if(fabs(f0) > fabs(f1))
    {
      x0 = newx;
      f0 = newf;
    }
    else
    {
      x1 = newx;
      f1 = newf;
    }
    // Test for x convergence
    if(fabs(x1 - x0) < xTolerance) return newx;
    ++numIterations;
  }
  std::cerr << "Warning from secant: Reached maximum number of iterations" << std::endl;
  return newx;
}

// Bisection method

template <typename Type>
Real const bisection(Type *myClass, 
                     Real const (Type::*target)(Real const &),
                     Real const &xInitial1, Real const &xInitial2,
                     Real const &xTolerance = DEFAULT_ROOT_X_TOLERANCE,
                     Real const &targetTolerance = DEFAULT_ROOT_TARGET_TOLERANCE,
                     size_t const maxIterations = DEFAULT_BISECTION_MAX_ITERATIONS)
{
  // Initialize variables
  size_t numIterations = 0;
  Real newx = xInitial1;
  Real x0 = xInitial1;
  Real x1 = xInitial2;
  Real f0 = (myClass->*target)(x0);
  Real f1 = (myClass->*target)(x1);
  if(f0*f1 > 0.0)
  {
    std::cerr << "Error in bisection - initial points do not bracket a root" << std::endl;
    return x0;
  }
  
  // Main loop
  while(numIterations < maxIterations)
  {
    // Find new x
    newx = 0.5*(x0 + x1);
    Real newf = (myClass->*target)(newx);
    // Test for target function convergence
    if(fabs(newf) < targetTolerance) return newx;
    // Housekeeping
    if(newf*f0 > 0.0)
    {
      x0 = newx;
      f0 = newf;
    }
    else
    {
      x1 = newx;
      f1 = newf;
    }
    // Test for x convergence
    if(fabs(x1 - x0) < xTolerance) return newx;
    ++numIterations;
  }
  std::cerr << "Warning from bisection: Reached maximum number of iterations" << std::endl;
  return newx;
}

#endif
