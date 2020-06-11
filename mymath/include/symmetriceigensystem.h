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
** Find all eigenvalues and eigenvectors of a symmetric matrix.
**
** Notes:
**
** - The eigenvalues are sorted in ascending order
** - Although this class contains functions to return the determinant,
**   trace, and whether the matrix is singular, this is a very
**   inefficient way to find those quantities. If that's all you
**   need, consider using the trace() method in the corresponding
**   matrix class, or class SymmetricEigensystem.
** - See also notes in linearalgebra.h
**
** Adapted from the public domain JAMA library implementation 
** of the EISPACK routines available on the netlib repository.
** JAMA can be obtained at http://math.nist.gov/javanumerics/jama/,
** the netlib repository is at http://www.netlib.org.
*/

#ifndef H_SYMMETRIC_EIGENSYSTEM
#define H_SYMMETRIC_EIGENSYSTEM

#include "common/include/types.h"
#include "common/include/assert.h"
#include <iostream>
#include <cmath>
#include <vector>

size_t const MAX_QL_ITERATIONS = 200; // Maximum number of iterations for eigenvalue decomposition

template<typename Matrix, typename Vector>
class SymmetricEigensystem
{
protected:

  SymmetricEigensystem(); // Prevent construction of empty eigensystem

public:

// Constructor
  
  explicit SymmetricEigensystem(Matrix const &A); // Create the eigensystem for matrix A

// Accessors 

  Vector const &eigenvalues() const;               // Returns the vector of eigenvalues of A
  Matrix const &eigenvectors() const;              // Returns the matrix of eigenvectors of A
  Real const &eigenvalue(size_t const i) const;    // Returns the i-th eigenvalue of A
  Vector const eigenvector(size_t const i) const;  // Returns the i-th eigenvector of A 
  bool const isSingular() const;                   // Returns true if the matrix is singular
  Real const &determinant() const;                 // Returns the determinant of A
  Real const &trace() const;                       // Returns the trace of A

private:

  Vector eigenvalues_;  // The eigenvalues
  Matrix eigenvectors_; // The eigenvector matrix
  Real trace_;          // The trace
  Real determinant_;    // The determinant
};

#define HYPOTH(a,b) ((fabs(a)>fabs(b))?fabs(a)*sqrt(1.0+(b*b/a/a)):fabs(b)*sqrt(1.0+(a*a/b/b)))
// ^^ To avoid underflow 

template<typename Matrix, typename Vector>
inline SymmetricEigensystem<Matrix, Vector>::SymmetricEigensystem(Matrix const &A)
{
  assert(Matrix().size() == Vector().size());

  size_t const n = A.size();

  // Householder reduction to tridiagonal matrix (EISPACK's tred2)

  eigenvectors_ = A;
  eigenvalues_ = A.row(n - 1);
  Vector subdiagonal = eigenvalues_; // After Householder reduction, contains the elements under the diagonal

  for(size_t i = n - 1; i > 0; i--)
  {
    Real h = 0.0;     // Normalization constant
    Real scale = 0.0; // Scaling factor to prevent over/underflow
    for(size_t j = 0; j < i; j++)
      scale += fabs(eigenvalues_(j));
    if(scale == 0.0)
    {
      subdiagonal(i) = eigenvalues_(i - 1);
      for(size_t j = 0; j < i; j++)
      {
        eigenvalues_(j) = eigenvectors_(i - 1, j);
        eigenvectors_(i, j) = eigenvectors_(j, i) = 0.0;
      }
    }
    else
    {
      // Householder vector
      for(size_t j = 0; j < i; j++)
      {
        eigenvalues_(j) /= scale;
        h += eigenvalues_(j)*eigenvalues_(j);
      }
      Real f = eigenvalues_(i - 1);
      Real g = (f > 0)?-sqrt(h):sqrt(h);
      subdiagonal(i) = scale*g;
      h -= f*g;
      eigenvalues_(i - 1) = f - g;
      for(size_t j = 0; j < i; j++)
        subdiagonal(j) = 0.0;
      // Transform remaining columns
      for(size_t j = 0; j < i; j++)
      {
        f = eigenvalues_(j);
        eigenvectors_(j, i) = f;
        g = subdiagonal(j) + eigenvectors_(j, j)*f;
        for(size_t k = j + 1; k <= i - 1; k++)
        {
          g += eigenvectors_(k, j)*eigenvalues_(k);
          subdiagonal(k) += f*eigenvectors_(k, j);
        }
        subdiagonal(j) = g;
      }
      f = 0.0;
      for(size_t j = 0; j < i; j++)
      {
        subdiagonal(j) /= h;
        f += subdiagonal(j)*eigenvalues_(j);
      }
      Real hh = f/(h + h);
      for(size_t j = 0; j < i; j++)
        subdiagonal(j) -= hh*eigenvalues_(j);
      for(size_t j = 0; j < i; j++)
      {
        f = eigenvalues_(j);
        g = subdiagonal(j);
        for(size_t k = j; k <= i-1; k++)
          eigenvectors_(k, j) -= (f*subdiagonal(k) + g*eigenvalues_(k));
        eigenvalues_(j) = eigenvectors_(i - 1, j);
        eigenvectors_(i, j) = 0.0;
      }
    }
    eigenvalues_(i) = h;
  }

  // Accumulate Householder transformations
  for(size_t i = 0; i < n - 1; i++)
  {
    eigenvectors_(n - 1, i) = eigenvectors_(i, i);
    eigenvectors_(i, i) = 1.0;
    Real h = eigenvalues_(i + 1);
    if(h != 0.0)
    {
      for(size_t j = 0; j <= i; j++)
        eigenvalues_(j) = eigenvectors_(j, i+1)/h;
      for(size_t j = 0; j <= i; j++)
      {
        Real g = 0.0;
        for(size_t k = 0; k <= i; k++)
          g += eigenvectors_(k, i+1)*eigenvectors_(k, j);
        for(size_t k = 0; k <= i; k++)
          eigenvectors_(k, j) -= g*eigenvalues_(k);
      }
    }
    for(size_t j = 0; j <= i; j++)
      eigenvectors_(j, i + 1) = 0.0;
  }
  for(size_t i = 0; i < n; i++)
  {
    eigenvalues_(i) = eigenvectors_(n - 1, i);
    eigenvectors_(n - 1, i) = 0.0;
  }
  eigenvectors_(n - 1, n - 1) = 1.0;
  subdiagonal(0) = 0.0;

  // QL algorithm for symmetric tridiagonal matrices (EISPACK's tql2)

  for(size_t i = 1; i < n; i++)
    subdiagonal(i - 1) = subdiagonal(i);
  subdiagonal(n - 1) = 0.0;

  Real f = 0.0;
  Real test1 = 0.0;
  Real epsilon = pow(2.0, -52.0);

  // Find small subdiagonal element
  for(size_t l = 0; l < n; l++)
  {
    Real test2 = fabs(eigenvalues_(l)) + fabs(subdiagonal(l));
    test1 = (test1 > test2)?test1:test2;
    size_t m = l;
    while (m < n)
    {
      if(fabs(subdiagonal(m)) <= epsilon*test1) break;
      if(m == n - 1) break; // To avoid crashes in case of underflow
      m++;
    }
    // If m = l, we have an eigenvalue. Iterate if m > l.
    if(m > l)
    {
      size_t numIterations = 0;
      do
      {
        numIterations++;
        if(numIterations > MAX_QL_ITERATIONS) break;
        // Implicit shift
        Real g = eigenvalues_(l);
        Real p = (eigenvalues_(l + 1) - g)/(subdiagonal(l) + subdiagonal(l));
        //Real r = (p > 0)?sqrt(1.0 + p*p):-sqrt(1.0 + p*p);
        Real r = (p > 0)?HYPOTH(1.0,p):-HYPOTH(1.0,p);
        eigenvalues_(l) = subdiagonal(l)/(p + r);
        eigenvalues_(l + 1) = subdiagonal(l)*(p + r);
        Real nextEigenvalue = eigenvalues_(l + 1);
        Real h = g - eigenvalues_(l);
        for(size_t i = l+2; i < n; i++)
          eigenvalues_(i) -= h;
        f += h;
        // Implicit QL transformation
        p = eigenvalues_(m);
        Real c = 1.0; Real c2 = 1.0; Real c3 = 1.0;
        Real nextSubdiagonal = subdiagonal(l + 1);
        Real s = 0.0; Real s2 = 0.0;
        for(size_t i = m - 1; i >= l; i--)
        {
          c3 = c2; c2 = c; s2 = s;
          g = c*subdiagonal(i);
          h = c*p;
          //r = sqrt(p*p + subdiagonal(i)*subdiagonal(i));
          r = HYPOTH(p,subdiagonal(i));
          subdiagonal(i + 1) = s*r;
          s = subdiagonal(i)/r;
          c = p/r;
          p = c*eigenvalues_(i) - s*g;
          eigenvalues_(i + 1) = h + s*(c*g + s*eigenvalues_(i));
          // Accumulate QL transformation
          for(size_t j = 0; j < n; j++)
          {
            h = eigenvectors_(j, i+1);
            eigenvectors_(j, i+1) = s*eigenvectors_(j, i) + c*h;
            eigenvectors_(j, i) = c*eigenvectors_(j, i) - s*h;
          }
          if(i == 0) break; // Avoid unsigned underflow
        }
        p = -s*s2*c3*nextSubdiagonal*subdiagonal(l)/nextEigenvalue;
        subdiagonal(l) = s*p;
        eigenvalues_(l) = c*p;
      }
      while(fabs(subdiagonal(l)) > epsilon*test1);
      if(numIterations > MAX_QL_ITERATIONS)
      {
        std::cerr << "Too many iterations in symmetricEigensystem" << std::endl;
        return;
      }
    }
    eigenvalues_(l) += f;
    subdiagonal(l) = 0.0;
  }
  // Sort
  for(size_t i = 0; i < n-1; i++)
  {
    size_t imin = i;
    Real minValue = eigenvalues_(i);
    for(size_t j = i+1; j < n; j++)
    {
      if(eigenvalues_(j) < minValue)
      {
        imin = j;
        minValue = eigenvalues_(j);
      }
    }
    if(imin != i)
    {
      eigenvalues_(imin) = eigenvalues_(i);
      eigenvalues_(i) = minValue;
      for(size_t j = 0; j < n; j++)
      {
        minValue = eigenvectors_(j, i);
        eigenvectors_(j, i) = eigenvectors_(j, imin);
        eigenvectors_(j, imin) = minValue;
      }
    }
  }
  // Get determinant and trace
  determinant_ = 1.0;
  trace_ = 0.0;
  for(size_t i = 0; i < n; ++i)
  {
    determinant_ *= eigenvalues_[i];
    trace_ += eigenvalues_[i];
  }
}

#undef HYPOTH

template<typename Matrix, typename Vector>
inline Vector const &SymmetricEigensystem<Matrix, Vector>::eigenvalues() const
{
  return eigenvalues_;
}

template<typename Matrix, typename Vector>
inline Matrix const &SymmetricEigensystem<Matrix, Vector>::eigenvectors() const
{
  return eigenvectors_;
}

template<typename Matrix, typename Vector>
inline Real const &SymmetricEigensystem<Matrix, Vector>::eigenvalue(size_t const i) const
{
  assert(i < eigenvalues_.size());
  return eigenvalues_[i];
}

template<typename Matrix, typename Vector>
inline Vector const SymmetricEigensystem<Matrix, Vector>::eigenvector(size_t const i) const
{
  assert(i < eigenvectors_.size());
  return eigenvectors_.column(i);
}

template<typename Matrix, typename Vector>
inline bool const SymmetricEigensystem<Matrix, Vector>::isSingular() const
{
  return (determinant_ == 0.0);
}

template<typename Matrix, typename Vector>
inline Real const &SymmetricEigensystem<Matrix, Vector>::determinant() const
{
  return determinant_;
}

template<typename Matrix, typename Vector>
inline Real const &SymmetricEigensystem<Matrix, Vector>::trace() const
{
  return trace_;
}

/*
** End of class SymmetricEigensystem
*/

#endif

