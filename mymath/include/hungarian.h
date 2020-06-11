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
** Hungarian algorithm
**
** Given a cost matrix, finds a permutation of the rows such that
** the sum of the diagonal elements will be minimal.
** Many ideas for this code come from the explanation by Prof. Bob
** Pilgrim at http://216.249.163.93/bob.pilgrim/445/munkres.html
*/

#ifndef H_HUNGARIAN
#define H_HUNGARIAN

#include "common/include/types.h"
#include "mymath/include/vectorND.h"
#include "mymath/include/matrixND.h"

enum MaskValue { HUNGARIAN_NONE, HUNGARIAN_STAR, HUNGARIAN_PRIME };
enum HungarianTask { COVER_COLUMNS, COVER_ZEROS };

// The Hungarian algorithm

inline std::vector<size_t> const hungarianMatch(MatrixND const &costMatrix)
{
  size_t size = costMatrix.size();
  
  // Deal with silly cases
  if(size == 0) return std::vector<size_t>();
  if(size == 1) return std::vector<size_t>(1, 0);

  // Copy the cost matrix
  MatrixND cost(costMatrix);

  std::vector<bool> rowCover(size, false);                          // "Covered" rows
  std::vector<bool> colCover(size, false);                          // "Covered" columns
  std::vector<std::vector<MaskValue> >                              // Mask matrix (keeps track
    maskMatrix(size, std::vector<MaskValue>(size, HUNGARIAN_NONE)); // of zeros)
  
  // Subtract the smallest element from each row
  for(size_t i = 0; i < size; i++)
  {
    Real smallest = cost(i, 0);
    for(size_t j = 1; j < size; j++)
      if(cost(i, j) < smallest) smallest = cost(i, j);
    for(size_t j = 0; j < size; j++)
      cost(i, j) -= smallest;
  }
  // Find starred zeros
  for(size_t i = 0; i < size; i++)
    for(size_t j = 0; j < size; j++)
      if(cost(i, j) == 0.0 && !rowCover[i] && !colCover[j])
      {
        maskMatrix[i][j] = HUNGARIAN_STAR;
        rowCover[i] = colCover[j] = true;
      }
  rowCover = colCover = std::vector<bool>(size, false);
  
  // Find optimal cover
  bool done = false;
  HungarianTask nextTask = COVER_COLUMNS;
  while(!done)  // Main loop
  {
    switch(nextTask)
    {
      case(COVER_COLUMNS):
      {
        // Cover columns containing starred zeros
        for(size_t i = 0; i < size; i++)
          for(size_t j = 0; j < size; j++)
            if(maskMatrix[i][j] == HUNGARIAN_STAR) colCover[j] = true;
        // Count covered columns and choose next task
        size_t numCovered = 0;
        for(size_t j = 0; j < size; j++)
          if(colCover[j]) ++numCovered;
        if(numCovered >= size) done = true; 
        else nextTask = COVER_ZEROS;
      } // End of case COVER_COLUMNS
        break;
      case(COVER_ZEROS):
      {
        // Try to cover all zeros
        size_t const max = (size_t)(-1);
        size_t row = max;
        size_t col = max;
        bool doneZeros = false;
        while(!doneZeros)
        {
          row = col = max;
          // Find an uncovered zero
          bool foundZero = false;
          for(size_t i = 0; i < size; i++)
          {
            for(size_t j = 0; j < size; j++)
              if(cost(i, j) == 0.0 && !rowCover[i] && !colCover[j])
              {
                row = i;
                col = j;
                foundZero = true;
              }
            if(foundZero) break;
          }
          if(row == max) doneZeros = true;
          else
          {
            // Prime the uncovered zero and search for star in the same row
            maskMatrix[row][col] = HUNGARIAN_PRIME;
            bool foundStar = false;
            for(size_t j = 0; j < size; j++)
              if(maskMatrix[row][j] == HUNGARIAN_STAR) foundStar = true;
            if(foundStar)
            {
              col = max;
              for(size_t j = 0; j < size; j++)
                if(maskMatrix[row][j] == HUNGARIAN_STAR) col = j;
              rowCover[row] = true;
              colCover[col] = false;
            }
            else doneZeros = true;
          }
        } // Done covering zeros
        if(row == max)
        {
          // Find smallest uncovered element
          Real minVal = REAL_VERYBIG;
          for(size_t i = 0; i < size; i++)
            for(size_t j = 0; j < size; j++)
              if(!rowCover[i] && !colCover[j] && cost(i, j) < minVal)
                minVal = cost(i, j);
          // Update cost matrix
          for(size_t i = 0; i < size; i++)
            for(size_t j = 0; j < size; j++)
            {
              if(rowCover[i]) cost(i, j) += minVal;
              if(!colCover[j]) cost(i, j) -= minVal;
            }
          nextTask = COVER_ZEROS;
        }
        else
        {
          // Update primed and starred zeros
          std::vector<size_t> chainRows(1, row); // Row indices of chain of zeros;
          std::vector<size_t> chainCols(1, col); // Column indices of chain of zeros;
          // Build chain of alternating primed and starred zeros
          bool foundChain = false;
          while(!foundChain)
          {
            // Find a star in current column
            size_t rowStar = max;
            size_t const colIndex = chainCols[chainCols.size() - 1];
            for(size_t i = 0; i < size; i++)
              if(maskMatrix[i][colIndex] == HUNGARIAN_STAR)
                rowStar = i;
            if(rowStar == max) foundChain = true;
            else
            {
              chainRows.push_back(rowStar);
              chainCols.push_back(colIndex);
            }
            if(!foundChain)
            {
              // Find a prime in current row
              size_t const rowIndex = chainRows[chainRows.size() - 1];
              size_t colPrime = 0;
              for(size_t j = 0; j < size; j++)
                if(maskMatrix[rowIndex][j] == HUNGARIAN_PRIME) colPrime = j;
              chainRows.push_back(rowIndex);
              chainCols.push_back(colPrime);
            }            
          } // End while(!foundChain)
          // Unstar starred zeros and star primed zeros
          for(size_t i = 0; i < chainCols.size(); i++)
            if(maskMatrix[chainRows[i]][chainCols[i]] == HUNGARIAN_STAR)
              maskMatrix[chainRows[i]][chainCols[i]] = HUNGARIAN_NONE;
            else maskMatrix[chainRows[i]][chainCols[i]] = HUNGARIAN_STAR;
          // Clear the covers and the primes
          for(size_t i = 0; i < size; i++)
            rowCover[i] = colCover[i] = false;
          for(size_t i = 0; i < size; i++)
            for(size_t j = 0; j < size; j++)
              if(maskMatrix[i][j] == HUNGARIAN_PRIME) maskMatrix[i][j] = HUNGARIAN_NONE;
          nextTask = COVER_COLUMNS;
        } // End if(row == max)
      } // End of case COVER_ZEROS
        break;
    } // End of switch(nextTask)
  } // End of main loop
  // Store matching on output array and return
  std::vector<size_t> matching;
  for(size_t i = 0; i < size; i++)
    for(size_t j = 0; j < size; j++)
      if(maskMatrix[i][j] == HUNGARIAN_STAR)
        matching.push_back(j);
  return matching;
}

#endif
