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
** Sample class for optimization routines
**
** The function in this example is the Mueller potential energy surface
** (K. Mueller and L.D. Brown, Theoret. Chim. Acta 53, 75 (1979)).
** It has the three local minima: (-0.558224, 1.44173), (-0.0500108, 0.466694)
** and (0.623499, 0.0280378).
*/

#ifndef H_OPTIMIZATION_EXAMPLE
#define H_OPTIMIZATION_EXAMPLE

#include <cmath>
#include "common/include/types.h"
#include "mymath/include/vector4D.h"
#include "mymath/include/vectorND.h"
#include "mymath/include/optimization.h"

class OptimizationExample
{
public:

// Constructors

  OptimizationExample();

// Target function and its gradient

  Real const targetFunction(VectorND const &x);
  VectorND const targetGradient(VectorND const &x);

// Minimization test

  VectorND const minimize(VectorND const &x0);

private:

  Vector4D A, a, b, c, x_0, y_0;  // Parameters

};

/*
** End of class optimizationExample
*/

// Inlines

inline OptimizationExample::OptimizationExample()
:
A(-200.0, -100.0, -170.0, 15.0), a(-1.0, -1.0, -6.5, 0.7), b(0.0, 0.0, 11.0, 0.6),
c(-10.0, -10.0, -6.5, 0.7), x_0(1.0, 0.0, -0.5, -1.0), y_0(0.0, 0.5, 1.5, 1.0)
{}

inline Real const OptimizationExample::targetFunction(VectorND const &x)
{
  Real energy = 0.0;
  for(size_t i = 0; i < 4; i++)
    energy += A[i]*exp(a[i]*(x[0]-x_0[i])*(x[0]-x_0[i]) + 
                       b[i]*(x[0]-x_0[i])*(x[1]-y_0[i]) +
                       c[i]*(x[1]-y_0[i])*(x[1]-y_0[i]));
  return energy;
}

inline VectorND const OptimizationExample::targetGradient(VectorND const &x)
{
  Real factor = 0.0;
  VectorND gradient(2, 0.0);
  for(size_t i = 0; i < 4; i++)
  {
    factor = A[i]*exp(a[i]*(x[0]-x_0[i])*(x[0]-x_0[i]) + 
                      b[i]*(x[0]-x_0[i])*(x[1]-y_0[i]) +
                      c[i]*(x[1]-y_0[i])*(x[1]-y_0[i]));
    gradient[0] += factor*(2*a[i]*(x[0]-x_0[i]) + b[i]*(x[1]-y_0[i]));
    gradient[1] += factor*(b[i]*(x[0]-x_0[i]) + 2*c[i]*(x[1]-y_0[i]));
  } 
  return gradient;
}

VectorND const OptimizationExample::minimize(VectorND const &x0)
{
  return bfgs<OptimizationExample>(this, &OptimizationExample::targetFunction, &OptimizationExample::targetGradient, x0);
}

#endif
