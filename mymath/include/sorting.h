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
** Sorting routines - add as needed
**
** Note that sorting is really implemented as indexing - the original
** array is not changed. Sorting is done in ascending order.
**
** All sorting algorithms are implemented as templates. The first
** type parameter corresponds to the container type of the data
** to sort, while the second type parameter corresponds to the
** type of the elements. For example, TypeData can be std::vector<size_t>
** with the corresponding TypeElement being <size_t>. TypeData must
** support the operator [] to extract elements, and the size() function
** returning the number of elements. TypeElement must, obviously, support
** comparison. Note that calling the sorting routine with incompatible
** TypeData and TypeElement (e.g. std::vector<Real> and int) can give
** unpredictable results.
*/

/*
** Some notes:
**
** - Need to implement a "real" introSort - the one below uses the STL's
**   version, but wastes time copying data.
*/

#ifndef H_SORTING
#define H_SORTING

#include <algorithm>
#include <vector>
#include "common/include/types.h"
#include "common/include/assert.h"

typedef std::vector<size_t> IndexTable;
typedef std::vector<size_t> RankTable;

// Utility to convert an index table to a rank table
// Note that repeated entries have the same rank

template<typename TypeData, typename TypeElement>
RankTable const rankTable(TypeData const &data,
                          IndexTable const &indexTable)
{
  RankTable rankTable(indexTable.size(), 0);
  rankTable[indexTable[0]] = 0;
  for(size_t i = 1; i < indexTable.size(); i++)
  {
    if(data[indexTable[i]] == data[indexTable[i - 1]])
      rankTable[indexTable[i]] = rankTable[indexTable[i - 1]];
    else
     rankTable[indexTable[i]] = i;
  }
  return rankTable;
}

// Insertion sort

template <typename TypeData, typename TypeElement>
IndexTable const insertionSort(TypeData const &data)
{
  // Store data array size
  size_t size = data.size();
  // Initialize index table
  IndexTable index;
  for(size_t i = 0; i < size; i++)
    index.push_back(i);
  // Check consistency
  assert(sizeof(data[0]) == sizeof(TypeElement));
  // The insertion sort algorithm
  for(size_t i = 1; i < size; i++)
  {
    TypeElement element = data[i];
    int j = (int)i - 1;
    while(j >= 0 && data[index[j]] > element)
    {
      index[j + 1] = index[j];
      j--;
    }
    index[j + 1] = i;
  }
  return index;
}

// Comparison (for use with introsort below)

template <typename TypeElement>
bool comparePair(std::pair<TypeElement, size_t> p1,
                 std::pair<TypeElement, size_t> p2)
{
  return p1.first < p2.first;
}

// Introsort
// Note: This is not my implementation, but
// a "translation" of the STL's. 

template <typename TypeData, typename TypeElement>
IndexTable const introSort(TypeData const &data)
{
  std::vector<std::pair<TypeElement, size_t> > v;
  for(size_t i = 0; i < data.size(); i++)
    v.push_back(std::pair<TypeElement, size_t>(data[i], i));
  std::sort(v.begin(), v.end(), comparePair<TypeElement>);
  IndexTable table;
  for(size_t i = 0; i < data.size(); i++)
    table.push_back(v[i].second);
  return table;
}

#endif
