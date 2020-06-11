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
** Sample classes for root finding routines
**
** The single-variable example solves x^3 = cos(x)
*/

#ifndef H_ROOT_FINDING_EXAMPLE
#define H_ROOT_FINDING_EXAMPLE

#include <cmath>
#include "common/include/types.h"
#include "mymath/include/rootfinding.h"

class RootFinding1DExample
{
public:

// Constructors

  RootFinding1DExample();

// Target function and its derivative

  Real const targetFunction(Real const &x);
  Real const targetDerivative(Real const &x);

// Root finding test

  Real const findRootNewton(Real const &x0);
  Real const findRootSecant(Real const &x1, Real const &x2);
  Real const findRootBisection(Real const &x1, Real const &x2);
};

/*
** End of class RootFinding1DExample
*/

// Inlines

inline RootFinding1DExample::RootFinding1DExample()
{}

inline Real const RootFinding1DExample::targetFunction(Real const &x)
{
  return x*x*x - cos(x);
}

inline Real const RootFinding1DExample::targetDerivative(Real const &x)
{
  return 3.0*x*x + sin(x);
}

Real const RootFinding1DExample::findRootNewton(Real const &x0)
{
  return newtonRaphson<RootFinding1DExample>(this, &RootFinding1DExample::targetFunction, &RootFinding1DExample::targetDerivative, x0);
}

Real const RootFinding1DExample::findRootSecant(Real const &x1, Real const &x2)
{
  return secant<RootFinding1DExample>(this, &RootFinding1DExample::targetFunction, x1, x2);
}

Real const RootFinding1DExample::findRootBisection(Real const &x1, Real const &x2)
{
  return bisection<RootFinding1DExample>(this, &RootFinding1DExample::targetFunction, x1, x2);
}

#endif
