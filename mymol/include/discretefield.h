/*
** Copyright 2012 Erik Santiso.
** This file is part of mymol.
** mymol is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** mymol is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with mymol. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** Discrete field class
*
** This class implements a discrete field defined on a grid in a
** periodic system. It associates with each grid cell a property that
** can be a scalar, vector, or any other type that supports basic
** addition, multiplication by scalar, and construction from scalar.
**
** The class also implements a simple cubic switching function to
** allow smoothing of the data. The switching function width is a
** real value between 0.0 and 0.5, and measures the extent to which
** a point in a given cell contributes to the adjacent cells. A value
** of zero means that no switching function is used.
**
** The switching function only affects the adjacent cells in scaled
** coordinates. If the lattice is very tilted or elongated, this may 
** mean that cells with centers that are closer may not get any weight.
** Ideally one would try to use a grid such that cells have similar
** lengths in all directions.
**
** The switching function is conservative, in the sense that the
** sum of the contributions to all cells equals the original value.
**
** Note that this is a template class but some of the implementation
** is in src/discretefield.cpp. If you get a linker error it probably
** means that you need to instantiate the template with the Type
** you want to use, recompile mymol and relink. See the end of 
** src/discretefield.cpp.
*/

#ifndef H_DISCRETEFIELD
#define H_DISCRETEFIELD

#include <iostream>
#include <vector>
#include <string>
#include "common/include/types.h"
#include "mymol/include/grid.h"

// CellData - Simple container to store the list of cells to which
// a position contributes. Used to implement the switching function.

struct CellData
{
  std::vector<size_t> cellIndices; // Indices of the cells
  std::vector<Real>   cellWeights; // Weight in each cell

  // Clear the indices and weights
  inline void clear()
  {
    cellIndices.clear();
    cellWeights.clear();
  }

  // Return the number of entries
  inline size_t const size() const
  {
    return cellIndices.size();
  }

  // Output (mostly for debugging)
  friend std::ostream &operator<<(std::ostream &outStream, CellData const &data)
  {
    outStream << "Indices: ";
    for(size_t i = 0; i < data.cellIndices.size(); ++i)
      outStream << data.cellIndices[i] << " ";
    outStream << std::endl << "Weights: ";
    for(size_t i = 0; i < data.cellWeights.size(); ++i)
      outStream << data.cellWeights[i] << " ";
    return outStream;
  }
};

// End of CellData

template<typename Type>
class DiscreteField
{
public:

// Constructors

  explicit DiscreteField(std::string const &name = ""); // Defines an empty discrete field with an
                                                        // optionally given name
  DiscreteField(Grid const &grid,                       // Defines a discrete field on a given grid with
                Real const &switchWidth = 0,            // optionally given switching function width
                std::string const &name = "");          // and name

// Interface

  void addValue(Vector3D const &position, // Adds the given value at the given position to the
                Type const &value);       // discrete field (using switching function if given)
  void addValue(size_t const cellIndex,   // Adds the given value directly to the given cell index
                Type const &value);
  void setGrid(Grid const &grid);         // Set the grid (this resets the field values and weights)
  void setSwitchWidth(Real const &width); // Set the switching function width (resets the field values and weights)
  void setName(std::string const &name);  // Set the data label
  void clear();                           // Resets all the field values and weights to zero
  void reset();                           // Clears the field values and the grid 

// Accessors

  Grid const &grid() const;                          // Returns the grid
  Real const &switchWidth() const;                   // Returns the switching function width
  std::string const &name() const;                   // Returns the data label
  size_t const numPoints() const;                    // Returns the total number of points added so far
  Type const &totalValue(size_t const i) const;      // Returns the total field value at the i-th cell
  Real const &totalWeight(size_t const i) const;     // Returns the total weight (or total number of values if no
                                                     // switching function is used) added to the i-th cell
  Type const value(size_t const i,                   // Return the mean field value at the i-th cell (total/weight)
                   bool const silent = true) const;  // - set silent to false to get error messages when weight is
                                                     // zero (otherwise value is silently set to zero)

// Print field information (mostly for debugging)
  template<typename ValueType>
  friend std::ostream &operator<<(std::ostream &outStream, DiscreteField<ValueType> const &field);

private:

// Functions

  void resizeValues(size_t const size);              // Resizes the field values and cell weights
                                                     // to the given size.
  CellData const cellData(Vector3D const &position); // Finds the indices of the cells that the given
                                                     // position contributes to, and the switching
                                                     // function weights.

// Data members

  Real switchWidth_;               // Switching function width
  Grid grid_;                      // The grid
  size_t numPoints_;               // The total number of points added so far
  std::vector<Type> totalValues_;  // The total field values added to each cell
  std::vector<Real> totalWeights_; // The total weights (or total number of values if no switching
                                   // function is used) added to each cell
  std::string name_;               // Label for the data
};

/*
** End of class DiscreteField
*/ 

// Inlines (everything needs to be here for template)

template<typename Type>
inline void DiscreteField<Type>::setGrid(Grid const &grid)
{
  grid_ = grid;
  resizeValues(grid.numCells());
  clear();
}

template<typename Type>
inline void DiscreteField<Type>::setSwitchWidth(Real const &width)
{
  assert(width >= 0.0 && width <= 0.5);
  switchWidth_ = width;
  clear();
}

template<typename Type>
inline void DiscreteField<Type>::setName(std::string const &name)
{
  name_ = name;
}

template<typename Type>
inline void DiscreteField<Type>::clear()
{
  numPoints_ = 0;
  for(size_t i = 0; i < grid_.numCells(); ++i)
    totalValues_[i] = totalWeights_[i] = 0.0;
}

template<typename Type>
inline void DiscreteField<Type>::reset()
{
  grid_.clear();
  resizeValues(grid_.numCells());
  clear();
}

template<typename Type>
inline Grid const &DiscreteField<Type>::grid() const
{
  return grid_;
}

template<typename Type>
inline Real const &DiscreteField<Type>::switchWidth() const
{
  return switchWidth_;
}

template<typename Type>
inline std::string const &DiscreteField<Type>::name() const
{
  return name_;
}

template<typename Type>
inline size_t const DiscreteField<Type>::numPoints() const
{
  return numPoints_;
}

template<typename Type>
inline Type const &DiscreteField<Type>::totalValue(size_t const i) const
{
  assert(i < totalValues_.size());
  return totalValues_[i];
}

template<typename Type>
inline Real const &DiscreteField<Type>::totalWeight(size_t const i) const
{
  assert(i < totalWeights_.size());
  return totalWeights_[i];
}

template<typename Type>
void DiscreteField<Type>::resizeValues(size_t const size)
{
  totalValues_.resize(size);
  totalWeights_.resize(size);
}

#endif

