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
** LU decomposition (Crout) with pivoting. Useful to solve linear systems
** of equations, inverting matrices, calculating determinants, etc.
**
** See also notes in linearalgebra.h
**
** Adapted from the public domain JAMA library implementation 
** of the EISPACK routines available on the netlib repository.
** JAMA can be obtained at http://math.nist.gov/javanumerics/jama/,
** the netlib repository is at http://www.netlib.org.
*/

#ifndef H_LU_DECOMPOSITION
#define H_LU_DECOMPOSITION

#include <vector>

template<typename Matrix, typename Vector>
class LUDecomposition
{
protected:

  LUDecomposition();           // Prevent construction of empty LU decomposition

public:

// Constructor
  
  explicit LUDecomposition(Matrix const &A); // Create a LU decomposition of the matrix A

// Accessors

  Matrix const getL() const;                    // Returns the lower-triangular matrix
  Matrix const getU() const;                    // Returns the upper-triangular matrix
  std::vector<size_t> const &pivot() const;     // Returns the pivot vector
  bool const isSingular() const;                // Returns true if the matrix is singular
  Real const &determinant() const;              // Returns the determinant of A
  Vector const solve(Vector const &b) const;    // Returns the solution of Ax=b for a vector b
  Matrix const solve(Matrix const &B) const;    // Returns the solution of AX=B for a matrix B
  std::vector<Vector> const 
    solve(std::vector<Vector> const &B) const ; // Returns the solution of Ax=b

private:

  Matrix LU_;                 // Lower and upper triangular matrices in compact form
  std::vector<size_t> pivot_; // The pivot array
  Real determinant_;          // The determinant
  int pivotSign_;             // Sign of the permutation (odd = -1, even = +1)
};

template<typename Matrix, typename Vector>
inline LUDecomposition<Matrix, Vector>::LUDecomposition(Matrix const &A)
{
  assert(Matrix().size() == Vector().size());
  LU_ = A;
  size_t size = A.size();
  for(size_t i = 0; i < size; ++i)
    pivot_.push_back(i);
  pivotSign_ = 1;

  // Loop over columns
  for(size_t j = 0; j < size; ++j)
  {
    Vector LUColj = LU_.column(j);
    for(size_t i = 0; i < size; ++i)
    {
      size_t max = (i < j)?i:j;
      Real dotProduct = 0.0;
      for(size_t k = 0; k < max; ++k)
        dotProduct += LU_(i, k)*LUColj[k];
      LU_(i, j) = LUColj[i] -= dotProduct;
    } // End of LU loop over rows
    // Pivoting
    size_t pivIndx = j;
    for(size_t i = j + 1; i < size; ++i)
      if(fabs(LUColj[i]) > fabs(LUColj[pivIndx])) pivIndx = i;
    if(pivIndx != j)
    {
      for(size_t k = 0; k < size; ++k)
      {
        Real dummy = LU_(pivIndx, k);
        LU_(pivIndx, k) = LU_(j, k);
        LU_(j, k) = dummy;
      }
      size_t dummy = pivot_[pivIndx];
      pivot_[pivIndx] = pivot_[j];
      pivot_[j] = dummy;
      pivotSign_ *= -1;
    } // End of pivoting
    if(j < size && LU_(j, j) != 0.0)
      for(size_t i = j + 1; i < size; ++i)
        LU_(i, j) /= LU_(j, j);
  } // End of loop over columns

  // Calculate the determinant
  determinant_ = 1.0;
  for(size_t i = 0; i < size; ++i)
    determinant_ *= LU_(i, i);
}

template<typename Matrix, typename Vector>
inline Matrix const LUDecomposition<Matrix, Vector>::getL() const
{
  Matrix L = LU_;
  for(size_t i = 0; i < LU_.size(); ++i)
  {
    L(i, i) = 1.0;
    for(size_t j = i + 1; j < LU_.size(); ++j)
      L(i, j) = 0.0;
  }
  return L;
}

template<typename Matrix, typename Vector>
inline Matrix const LUDecomposition<Matrix, Vector>::getU() const
{
  Matrix U = LU_;
  for(size_t i = 0; i < LU_.size(); ++i)
    for(size_t j = 0; j < i; ++j)
      U(i, j) = 0.0;

  return U;
}

template<typename Matrix, typename Vector>
inline std::vector<size_t> const &LUDecomposition<Matrix, Vector>::pivot() const
{
  return pivot_;
}

template<typename Matrix, typename Vector>
inline bool const LUDecomposition<Matrix, Vector>::isSingular() const
{
  return (determinant_ == 0.0);
}

template<typename Matrix, typename Vector>
inline Real const &LUDecomposition<Matrix, Vector>::determinant() const
{
  return determinant_;
}

template<typename Matrix, typename Vector>
inline Vector const LUDecomposition<Matrix, Vector>::solve(Vector const &b) const
{
  assert(LU_.size() == b.size());
  assert(!isSingular());
  // Permute the rows
  Vector x = b; // To copy type
  size_t size = b.size();
  for(size_t i = 0; i < size; ++i)
    x(i) = b(pivot_[i]);
  // Solve for U*x
  for(size_t k = 0; k < size; ++k)
  for(size_t i = k + 1; i < size; ++i)
    x(i) -= x(k)*LU_(i, k);
  // Solve for x
  for(size_t k = size - 1; k >= 0; --k)
  {
    x(k) /= LU_(k, k);
    for(size_t i = 0; i < k; ++i)
      x(i) -= x(k)*LU_(i, k);
    if(k == 0) break;
  }
  return x;
}

template<typename Matrix, typename Vector>
inline Matrix const LUDecomposition<Matrix, Vector>::solve(Matrix const &B) const
{
  assert(LU_.size() == B.size());
  assert(!isSingular());
  // Permute the rows
  Matrix X = B; // To copy the type
  size_t size = B.size();
  for(size_t i = 0; i < size; ++i)
    X.setRow(i, B.row(pivot_[i]));
  // Solve for U*X
  for(size_t k = 0; k < size; ++k)
  for(size_t i = k + 1; i < size; ++i)
  for(size_t j = 0; j < size; ++j)
    X(i, j) -= X(k, j)*LU_(i, k);
  // Solve for X
  for(size_t k = size - 1; k >= 0; --k)
  {
    for(size_t j = 0; j < size; ++j)
      X(k, j) /= LU_(k, k);
    for(size_t i = 0; i < k; ++i)
    for(size_t j = 0; j < size; ++j)
      X(i, j) -= X(k, j)*LU_(i, k);
    if(k == 0) break;
  }
  return X;
}

template<typename Matrix, typename Vector>
inline std::vector<Vector> const LUDecomposition<Matrix, Vector>::solve(std::vector<Vector> const &B) const
{
  if(B.size() == 0) return std::vector<Vector>();
  assert(LU_.size() == B[0].size());
  assert(!isSingular());
  // Permute the rows
  std::vector<Vector> X = B; // To copy the type
  size_t size = B[0].size();
  size_t numVecs = B.size();
  for(size_t i = 0; i < size; ++i)
  for(size_t j = 0; j < numVecs; ++j)
    X[j][i] = B[j][pivot_[i]];
  // Solve for U*X
  for(size_t k = 0; k < size; ++k)
  for(size_t i = k + 1; i < size; ++i)
  for(size_t j = 0; j < numVecs; ++j)
    X[j][i] -= X[j][k]*LU_(i, k);
  // Solve for X
  for(size_t k = size - 1; k >= 0; --k)
  {
    for(size_t j = 0; j < numVecs; ++j)
      X[j][k] /= LU_(k, k);
    for(size_t i = 0; i < k; ++i)
    for(size_t j = 0; j < numVecs; ++j)
      X[j][i] -= X[j][k]*LU_(i, k);
    if(k == 0) break;
  }
  return X;
}

/*
** End of class LUDecomposition
*/

#endif

