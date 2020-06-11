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
** Special functions - add as needed
*/

/*
** Some notes:
**
** - Many of the formulas used here are really single-precision
**   (error ~ 1E-7). For a later version it would be good to
**   update to double precision - see e.g. Zhang and Jin,
**   "Computation of Special Functions", Wiley (1996).
** - The calculation of the log derivative of 1F1 for large
**   negative x is not implemented. This requires doing the 
**   perturbation series with the first term in the asymptotic
**   expansion in Abramowitz and Stegun, Chapter 13. I do not
**   need this right now, so I am leaving it for a future
**   version.
*/

#ifndef H_SPECIAL_FUNCTIONS
#define H_SPECIAL_FUNCTIONS

#include <cmath>
#include <iostream>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/mpconstants.h"

// Some global constants for the special function routines

size_t const MAXIMUM_CONFLUENT_TERMS = 250;         // Maximum number of terms in series expansion
                                                    // (for 1F1 and related functions)
size_t const MAXIMUM_CONTINUED_FRACTIONS = 1000;    // Maximum number of terms in continued
                                                    // fraction expansion
Real   const FUNCTION_TOLERANCE = 2.0*SMALLREAL_EPSILON; // Tolerance for estimation of some functions
                                                         // (Would be REAL_EPSILON for double precision)

// Modified Bessel function of the first kind and order 0
// (Abramowitz and Stegun, section 9.8)
// Note that this has error ~ 2E-7

Real const bessel_I0(Real const &x)
{
  Real const absx = fabs(x);
  if(absx < 3.75)
  {
    Real const tsq = x*x/14.0625;
    return 1.0 + tsq*(3.5156229 + 
                 tsq*(3.0899424 +
                 tsq*(1.2067492 +
                 tsq*(0.2659732 +
                 tsq*(0.0360768 +
                 tsq*0.0045813)))));
  }
  else
  {
    Real const tm1 = 3.75/absx;
    return (exp(absx)/sqrt(absx))*(0.39894228 + 
                                   tm1*(0.01328592 +
                                   tm1*(0.00225319 +
                                   tm1*(-0.00157565 +
                                   tm1*(0.00916281 +
                                   tm1*(-0.02057706 +
                                   tm1*(0.02635537 +
                                   tm1*(-0.01647633 +
                                   tm1*0.00392377))))))));
  }
}

// Modified Bessel function of the first kind and order 1
// (Abramowitz and Stegun, section 9.8)
// Note that this has error ~ 2E-7

Real const bessel_I1(Real const &x)
{
  Real const absx = fabs(x);
  if(absx < 3.75)
  {
    Real const tsq = x*x/14.0625;
    return x*(0.5 + tsq*(0.87890594 + 
                    tsq*(0.51498869 +
                    tsq*(0.15084934 +
                    tsq*(0.02658733 +
                    tsq*(0.00301532 +
                    tsq*0.00032411))))));
  }
  else
  {
    Real const tm1 = 3.75/absx;
    Real const signx = (x<0)?-1.0:1.0;
    return signx*(exp(absx)/sqrt(absx))*(0.39894228 + 
                                         tm1*(-0.03988024 +
                                         tm1*(-0.00362018 +
                                         tm1*(0.00163801 +
                                         tm1*(-0.01031555 +
                                         tm1*(0.02282967 +
                                         tm1*(-0.02895312 +
                                         tm1*(0.01787654 +
                                         tm1*-0.00420059))))))));
  }
}

// Ratio I1/I0 (useful to avoid errors for large x)

Real const bessel_I1_over_I0(Real const &x)
{
  Real const absx = fabs(x);
  if(absx < 3.75) return bessel_I1(x)/bessel_I0(x);
  else
  {
    Real const tm1 = 3.75/absx;
    Real const signx = (x<0)?-1.0:1.0;
    Real const numerator = signx*(0.39894228 + 
                                  tm1*(-0.03988024 +
                                  tm1*(-0.00362018 +
                                  tm1*(0.00163801 +
                                  tm1*(-0.01031555 +
                                  tm1*(0.02282967 +
                                  tm1*(-0.02895312 +
                                  tm1*(0.01787654 +
                                  tm1*-0.00420059))))))));
    Real const denominator = (0.39894228 + 
                              tm1*(0.01328592 +
                              tm1*(0.00225319 +
                              tm1*(-0.00157565 +
                              tm1*(0.00916281 +
                              tm1*(-0.02057706 +
                              tm1*(0.02635537 +
                              tm1*(-0.01647633 +
                              tm1*0.00392377))))))));
    return numerator/denominator;
  }
}

// Gamma function
// Replaces the idiotic gamma() from the standard library
Real const euler_gamma(Real const &x)
{
  if(x == floor(x) && x <= 0.0) return REAL_VERYBIG;
  else if(x == 1.0) return 1.0;
  else if(x < 1.0) return euler_gamma(x + 1.0)/x;
  else if(x > 1.0 && x <= 2.0)
  {
    // Small x (Abramowitz and Stegun, chapter 6)
    Real y = x - 1.0;
    return 1.0 + 
           y*(-0.577191652 +
           y*(0.988205891 +
           y*(-0.897056937 +
           y*(0.918206857 +
           y*(-0.756704078 +
           y*(0.482199394 +
           y*(-0.193527818 +
           y*0.035868343)))))));
  }
  else if(x > 2.0 && x <= 20.0)  return (x - 1.0)*euler_gamma(x - 1.0);
  else  // x > 20
  {
    // Large x (Abramowitz and Stegun, chapter 6)
    Real invxsq = 1.0/(x*x);
    Real correction = 1.0/x*(1.0/12.0 +
                      invxsq*(-1.0/360.0 +
                      invxsq*(1.0/1260.0 +
                      invxsq*-1.0/1680.0)));
    return exp(-x + (x-0.5)*log(x) + LOG_SQRT_TWO_PI + correction);
  }
}

// Confluent (Kummer's) hypergeometric function of scalar argument 1F1(a, b, x)

Real const confluent(Real const &a, Real const &b, Real const &x)
{
  if(b == 0.0 || (b < 0.0 && b == floor(b))) return REAL_VERYBIG;
  if(x < -700.0) return 0.0;          // Underflow
  if(x > 700.0) return REAL_VERYBIG;  // Overflow
  if(x < 0.0) return exp(x)*confluent(b - a, b, -x);  // Kummer transformation
  if(x < 20.0)
  {
    // Use Taylor expansion
    Real newTerm = a/b*x;   // Term in series expansion
    Real M = 1.0 + newTerm; // Current value of function
    size_t numTerms = 1;
    while(numTerms < MAXIMUM_CONFLUENT_TERMS)
    {
      ++numTerms;
      Real n = (Real)numTerms;
      newTerm *= (a + n - 1.0)*x/(n*(b + n - 1.0));
      M += newTerm;
      if(fabs(newTerm/M) < FUNCTION_TOLERANCE) return M;
    }
    std::cerr << "Warning: Maximum number of terms reached in confluent" << std::endl;
    return M;
  }
  else
  {
    // Use large-x expansion (see Abramowitz and Stegun, chapter 13)
    Real const invx = 1.0/x;
    Real const aminusb = a - b;
    Real const gammab = euler_gamma(b);

    Real newTerm1, newTerm2;  // Terms in first and second series expansions
    Real M;                   // Current value of function
    newTerm1 = gammab*cos(PI*a)*pow(x, -a)/euler_gamma(-aminusb);
    newTerm2 = gammab*exp(x)*pow(x, aminusb)/euler_gamma(a);
    M = newTerm1 + newTerm2;
    size_t numTerms = 0;
    while (numTerms < MAXIMUM_CONFLUENT_TERMS)
    {
      ++numTerms;
      Real n = (Real)numTerms;
      newTerm1 *= -invx*(a + n - 1)*(aminusb + n)/n;
      newTerm2 *= invx*(n - 1 - aminusb)*(n - a)/n;
      Real newTerm = newTerm1 + newTerm2;
      M += newTerm;
      if(fabs(newTerm/M) < FUNCTION_TOLERANCE) return M;
    }
    std::cerr << "Warning: Maximum number of terms reached in confluent" << std::endl;
    return M;
  }
}

// Log derivative of the confluent (Kummer's) hypergeometric function of 
// scalar argument, dln(1F1(a, b, x))/dx

Real const d_ln_confluent(Real const &a, Real const &b, Real const &x)
{
  if(b == 0.0 || (b < 0.0 && b == floor(b))) return REAL_VERYBIG;
  else if(fabs(x) < 20.0) return (a/b)*confluent(a + 1.0, b + 1.0, x)/confluent(a, b, x);
  else if(fabs(x) < 1000.0)
  {
    // Continued fraction expansion from Cuyt et al., "Handbook of Continued
    // Fractions of Special Functions", Springer (2008), chapter 16
    Real aoverb = a/b;
    Real oldNum, oldDen;  // Old numerator and denominator
    Real newNum, newDen;  // New numerator and denominator
    Real dlnM = 1;        // Current value of function
    size_t numTerms = 0;  // Number of terms in continued fraction
    oldDen = 1.0; oldNum = 0.0;
    newDen = 1.0; newNum = 1.0;
    while(numTerms < MAXIMUM_CONFLUENT_TERMS)
    {
      Real olddlnM = dlnM;
      ++numTerms;
      Real n = (Real)numTerms;
      Real aTerm; // Numerator of general term
      if(n == 2.0*floor(n/2.0))
        aTerm = (a + 0.5*n)/((b + n - 1.0)*(b + n));
      else
        aTerm = -(b - a + 0.5*n - 0.5)/((b + n - 1.0)*(b + n));
      Real num = newNum + aTerm*oldNum*x;
      Real den = newDen + aTerm*oldDen*x;
      dlnM = aoverb*num/den;
      if(fabs((dlnM - olddlnM)/dlnM) < FUNCTION_TOLERANCE) return dlnM;
      oldNum = newNum; newNum = num;
      oldDen = newDen; newDen = den;
    }
    std::cerr << "Warning: Maximum number of terms reached in confluent" << std::endl;
    return dlnM;
  }
  else
  {
    // Asymptotic expansion (derived using a perturbation series from the
    // expression in Abramowitz and Stegun, Chapter 13)
    if(x < 0)
    {
      std::cerr << "Warning: d_ln_confluent not implemented for large negative x" << std::endl;
      return 0.0;
    }
    Real invx = 1.0/x;
    Real f0 = 1;
    Real f1 = a - b;
    Real f2 = f1*(1.0 - a);
    Real f3 = f2*(b + 2.0 - a - a);
    return f0 + invx*(f1 + invx*(f2 + invx*f3));
  }
}

#endif
