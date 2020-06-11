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
** Basic statistics routines - add as needed
**
** All functions are implemented as templates. When the input is
** a set of elements, the first type parameter corresponds to the 
** container type of the data to sort, while the second type parameter 
** corresponds to the type of the elements. For example, TypeData can 
** be std::vector<size_t> with the corresponding TypeElement being 
** <size_t>. TypeData must support the operator [] to extract elements, 
** and the size() function returning the number of elements. Note that 
** calling such a function with incompatible TypeData and TypeElement 
** (e.g. std::vector<Real> and int) can give unpredictable results.
*/

#ifndef H_STATISTICS
#define H_STATISTICS

#include <vector>
#include <cmath>
#include <ctime>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "randomc/include/randomc.h"

// Simple average

template <typename TypeData, typename TypeElement>
TypeElement const average(TypeData const &data)
{
  // Get number of elements
  size_t size = data.size();
  // Check consistency
  assert(sizeof(data[0]) == sizeof(TypeElement));
  // Calculate average
  TypeElement sum = 0;
  for(size_t i = 0; i < size; i++)
    sum += data[i];
  return sum/size;
}

// Simple variance

template <typename TypeData, typename TypeElement>
TypeElement const variance(TypeData const &data)
{
  // Get number of elements
  size_t size = data.size();
  // Check consistency
  assert(sizeof(data[0]) == sizeof(TypeElement));
  // Calculate standard deviation
  TypeElement sum = 0;
  for(size_t i = 0; i < size; i++)
    sum += data[i]*data[i];
  TypeElement avg = average<TypeData, TypeElement>(data);
  return (sum/size - avg*avg);
}

// Bootstrapping distribution of the mean
// Uses the randomc library by Agner Fog

template <typename TypeData, typename TypeElement>
std::vector<TypeElement> const bootstrap(TypeData const &data,
                                         size_t const numSamples)
{
  CRandomMersenne ranGen(0);       // Random number generator
  std::vector<TypeElement> result; // The bootstrap distribution
  ranGen.RandomInit((int)time(0)); // Seed random number generator
  size_t size = data.size();

  for(size_t iSample = 0; iSample < numSamples; ++iSample)
  {
    // Add a bootstrap sample of the mean
    Real sum = 0.0;
    for(size_t i = 0; i < size; ++i)
      sum += data[ranGen.IRandomX(0, size - 1)];
    result.push_back(sum/size);
  }
  return result;
}

#endif

