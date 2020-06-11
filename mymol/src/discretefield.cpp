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
** Discrete Field Class
*/

#include "common/include/assert.h"
#include "mymol/include/lattice.h"
#include "mymol/include/discretefield.h"

// Constructors

template<typename Type>
DiscreteField<Type>::DiscreteField(std::string const &name)
:
switchWidth_(0.0), grid_(), numPoints_(0), totalValues_(1, 0.0), totalWeights_(1, 0.0), name_(name)
{}

template<typename Type>
DiscreteField<Type>::DiscreteField(Grid const &grid, Real const &switchWidth, std::string const &name)
:
switchWidth_(switchWidth), grid_(grid), numPoints_(0),
totalValues_(grid.numCells(), 0.0), totalWeights_(grid.numCells(), 0.0), name_(name)
{
  if(switchWidth < 0.0 || switchWidth > 0.5)
  {
    std::cerr << "Error in DiscreteField: Switching function width must "
              << "be between 0.0 and 0.5" << std::endl;
    reset();
  }
}

// Interface

template<typename Type>
void DiscreteField<Type>::addValue(Vector3D const &position, Type const &value)
{
  // Get the cells to add to and weights
  CellData myCellData = cellData(position);

  // Add the values
  for(size_t i = 0; i < myCellData.size(); ++i)
  {
    size_t const &iCell = myCellData.cellIndices[i];
    Real const &weight = myCellData.cellWeights[i];
    if(iCell >= totalValues_.size() || iCell >= totalWeights_.size())
    {
      std::cerr << "Internal error in DiscreteField - bad cell index" << std::endl;
      reset();
    }
    totalValues_[iCell] += weight*value;
    totalWeights_[iCell] += weight;
 }
 
  // Update accumulator
  ++numPoints_;
}

template<typename Type>
void DiscreteField<Type>::addValue(size_t const cellIndex, Type const &value)
{
  assert(cellIndex < grid_.numCells());
  totalValues_[cellIndex] += value;
  totalWeights_[cellIndex] += 1.0;
 
  // Update accumulator
  ++numPoints_;
}

// Accessors

template<typename Type>
Type const DiscreteField<Type>::value(size_t const i, bool const silent) const
{
  assert(i < totalValues_.size());
  if(totalWeights_[i] == 0.0)
  {
    if(!silent)
      std::cerr << "Error in DiscreteField::value - weight is zero" << std::endl;
    return 0.0;
  }
  return totalValues_[i]/totalWeights_[i];
}

// Print field information (mostly for debugging)

template<typename Type>
std::ostream &operator<<(std::ostream &outStream, DiscreteField<Type> const &field)
{
  Grid const &grid = field.grid_;
  outStream << "Data label: " << field.name_ << std::endl;
  outStream << "Grid lattice information: " << std::endl;
  outStream << grid.lattice() << std::endl;
  outStream << "Grid origin: " << grid.origin() << std::endl;
  outStream << "Number of divisions in each direction: "
            << grid.numDivisions(0) << " "
            << grid.numDivisions(1) << " "
            << grid.numDivisions(2) << std::endl;
  outStream << "Number of cells: " << grid.numCells() << std::endl;
  outStream << "Switching function width: " << field.switchWidth_ << std::endl;
  outStream << "Total values added: " << field.numPoints_ << std::endl;
  outStream << "Field values (cell index, cell center, total value, total weight, mean value):";
  for(size_t i = 0; i < grid.numCells(); ++i)
  {
    std::cout << std::endl << i << " " << grid.cellCenter(i) << " " << field.totalValues_[i] << " "
              << field.totalWeights_[i] << " " << field.value(i, true);
  }
 
  return outStream;
}

// Private functions

template<typename Type>
CellData const DiscreteField<Type>::cellData(Vector3D const &position)
{
  CellData cellData;

  if(switchWidth_ <= 0.0)
  {
    // No switching function - use center of cell
    cellData.cellIndices.push_back(grid_.cellIndex(position));
    cellData.cellWeights.push_back(1.0);
    return cellData;
  }

  // Switching function: This repeats some stuff implemented in grid,
  // but it is faster and simpler to do the math in scaled coordinates

  // Get lattice and grid information  
  Lattice const &lattice = grid_.lattice();
  size_t const numDim = lattice.numDimensions();
  size_t const numDiv[3] = { grid_.numDivisions(0),
                             grid_.numDivisions(1),
                             grid_.numDivisions(2) };
  
  // Get scaled position, indices of nearest cell center, and nearest
  // cell center in scaled coordinates
  Vector3D const pos = position - grid_.origin();
  Vector3D scaled = 0.0;       // Scaled position
  Vector3D center = 0.0;       // Nearest cell center
  size_t iCell[3] = {0, 0, 0}; // Indices of nearest center in each direction

  for(size_t iDim = 0; iDim < numDim; ++iDim)
  {
    scaled[iDim] = lattice.reciprocalVector(iDim)*pos;
    scaled[iDim] -= 1 + floor(scaled[iDim] - 0.5); // Bring back to unit cell
    iCell[iDim] = (size_t)floor(numDiv[iDim]*(scaled[iDim] + 0.5 - floor(scaled[iDim] + 0.5)));
    if(iCell[iDim] > numDiv[iDim] - 1) iCell[iDim] = numDiv[iDim] - 1;
    center[iDim] = ((Real)iCell[iDim] + 0.5)/(Real)numDiv[iDim] - 0.5;
  }

  // Get the indices of the cells that this position contributes to,
  // and the weight for each cell
  std::vector<size_t> indx[3]; // Indices that contribute to weight in each direction
  std::vector<Real> weight[3]; // Contributions to weight in each direction

  Vector3D delta = scaled - center; // Distance from cell center in scaled coordinates
  for(size_t iDim = 0; iDim < 3; ++iDim)
  {
    indx[iDim].push_back(iCell[iDim]);
    Real scaledDiff = (fabs(numDiv[iDim]*delta[iDim]) - 0.5)/switchWidth_;

    if(scaledDiff > -1.0 && numDiv[iDim] > 1)
    {
      if(delta[iDim] > 0.0)
      {
        // Contributes to the cell ahead
        if(iCell[iDim] == numDiv[iDim] - 1) indx[iDim].push_back(0); // PBCs
        else indx[iDim].push_back(iCell[iDim] + 1);
      }
      else
      {
        // Contributes to the cell before
        if(iCell[iDim] == 0) indx[iDim].push_back(numDiv[iDim] - 1); // PBCs
        else indx[iDim].push_back(iCell[iDim] - 1);
      }
      Real const switchWeight = 0.5 + 0.25*scaledDiff*(scaledDiff*scaledDiff - 3);
      weight[iDim].push_back(switchWeight);
      weight[iDim].push_back(1.0 - switchWeight);
    }
    else
      weight[iDim].push_back(1.0);
  }

  // Multiply all the contributions together and return the data
  for(size_t ix = 0; ix < indx[0].size(); ++ix)
  for(size_t iy = 0; iy < indx[1].size(); ++iy)
  for(size_t iz = 0; iz < indx[2].size(); ++iz)
  {
    cellData.cellIndices.push_back(indx[0][ix] + numDiv[0]*(indx[1][iy] + numDiv[1]*indx[2][iz]));
    cellData.cellWeights.push_back(weight[0][ix]*weight[1][iy]*weight[2][iz]);
  }

  return cellData;
}

// Instantiations to be used (this is needed to be able to compile into library, add as needed)

template class DiscreteField<Real>;
template class DiscreteField<Vector3D>;

template std::ostream &operator<<(std::ostream &outStream, DiscreteField<Real> const &field);
template std::ostream &operator<<(std::ostream &outStream, DiscreteField<Vector3D> const &field);

