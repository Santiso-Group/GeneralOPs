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
** Simple grid class
*/

#include "common/include/assert.h"
#include "mymol/include/grid.h"

// Constructors

Grid::Grid()
:
lattice_(), origin_(0.0)
{
  for(size_t i = 0; i < 3; ++i) numDivisions_[i] = 1;
  updateGrid(); 
}

Grid::Grid(size_t const n, Lattice const &lattice, Vector3D const &origin)
:
lattice_(lattice), origin_(origin)
{
  assert(n > 0);
  for(size_t i = 0; i < 3; ++i) numDivisions_[i] = n;
  updateGrid();
}

Grid::Grid(size_t const nx, size_t const ny, size_t const nz, Lattice const &lattice, Vector3D const &origin)
:
lattice_(lattice), origin_(origin)
{
  assert(nx*ny*nz > 0); 
  numDivisions_[0] = nx;
  numDivisions_[1] = ny;
  numDivisions_[2] = nz;
  updateGrid();
}

// Other properties

size_t const Grid::cellIndex(Vector3D const position) const
{  
  size_t const numDim = lattice_.numDimensions();

  if(numDim == 0) return 0; // Trivial case

  size_t indices[3] = {0, 0, 0}; // Grid indices in each direction
  Vector3D const pos = position - origin_;

  for(size_t iDim = 0; iDim < numDim; ++iDim)
  {
    Real scaled = lattice_.reciprocalVector(iDim)*pos;
    indices[iDim] = (size_t)floor(numDivisions_[iDim]*(scaled + 0.5 - floor(scaled + 0.5)));
    if(indices[iDim] > numDivisions_[iDim] - 1) indices[iDim] = numDivisions_[iDim] - 1;
  }
  return indices[0] + numDivisions_[0]*(indices[1] + numDivisions_[1]*indices[2]);  
}

std::vector<size_t> const Grid::adjacentCells(size_t const i) const
{
  std::vector<size_t> adjacentIndices;
  std::vector<size_t> const &indices = cellIndices(i);

  for(size_t iDim = 0; iDim < 3; ++iDim)
  {
    std::vector<size_t> adjacent = indices;
    size_t const iplus = (indices[iDim] == numDivisions_[iDim] - 1)?0:(indices[iDim] + 1);
    adjacent[iDim] = iplus;
    adjacentIndices.push_back(cellIndex(adjacent));
    size_t const iminus = (indices[iDim] == 0)?(numDivisions_[iDim] - 1):(indices[iDim] - 1);
    adjacent[iDim] = iminus;
    adjacentIndices.push_back(cellIndex(adjacent));
  }
  return adjacentIndices; 
}

Vector3D const Grid::point(size_t const i) const
{
  assert(i < numPoints_);

  size_t indices[3] = {0, 0, 0}; // Grid indices in each direction
  size_t indx = i;

  size_t const nx = (numDivisions_[0] + 1);
  size_t const nxy = nx*(numDivisions_[1] + 1);

  indices[2] = indx/nxy;
  indx -= indices[2]*nxy;
  indices[1] = indx/nx;
  indx -= indices[1]*nx;
  indices[0] = indx;

  Vector3D center = origin_;
  size_t numDim = lattice_.numDimensions();

  for(size_t iDim = 0; iDim < numDim; ++iDim)
    center += (Real(indices[iDim])/Real(numDivisions_[iDim]) - 0.5)*
              lattice_.latticeVector(iDim);

  return center;
}

Vector3D const Grid::cellCenter(size_t const i) const
{
  assert(i < numCells_);

  size_t indices[3] = {0, 0, 0}; // Grid indices in each direction
  size_t indx = i;

  size_t const nxy = numDivisions_[0]*numDivisions_[1];

  indices[2] = indx/nxy;
  indx -= indices[2]*nxy;
  indices[1] = indx/numDivisions_[0];
  indx -= indices[1]*numDivisions_[0];
  indices[0] = indx;

  Vector3D center = origin_;
  size_t numDim = lattice_.numDimensions();

  for(size_t iDim = 0; iDim < numDim; ++iDim)
    center += ((Real(indices[iDim]) + 0.5)/Real(numDivisions_[iDim]) - 0.5)*
               lattice_.latticeVector(iDim);

  return center;
}

// Print grid information (mostly for debugging)

std::ostream &operator<<(std::ostream &outStream, Grid const &grid)
{
  outStream << "Grid lattice information: " << std::endl;
  outStream << grid.lattice_ << std::endl;
  outStream << "Grid origin: " << grid.origin_ << std::endl;
  outStream << "Number of divisions in each direction: "
            << grid.numDivisions_[0] << " "
            << grid.numDivisions_[1] << " "
            << grid.numDivisions_[2] << std::endl;
  outStream << "Number of points: " << grid.numPoints_ << std::endl;
  outStream << "Coordinates of all points (index, position): " << std::endl;
  for(size_t i = 0; i < grid.numPoints_; ++i)
   outStream << i << " " << grid.point(i) << std::endl; 
  outStream << "Number of cells: " << grid.numCells_ << std::endl;
  outStream << "Centers of all cells (index, position): ";
  for(size_t i = 0; i < grid.numCells_; ++i)
    outStream << std::endl << i << " " << grid.cellCenter(i);

  return outStream;
}

// Private methods

void Grid::updateGrid()
{
  // Check that number of lattice dimensions is consistent with grid divisions
  bool isConsistent = true;
  size_t numDim = lattice_.numDimensions();

  for(size_t iDim = numDim; iDim < 3; ++iDim)
  {
    if(numDivisions_[iDim] != 1)
    {
      isConsistent = false;
      numDivisions_[iDim] = 1;
    }
  }
  if(!isConsistent)
  {
    std::cerr << "Warning from Grid: Number of grid divisions was not consistent "
              << "with lattice dimension." << std::endl
              << "The extra divisions were removed." << std::endl;
  }

  // Calculate total number of points and cells
  numPoints_ = (numDivisions_[0] + 1)*(numDivisions_[1] + 1)*(numDivisions_[2] + 1);
  numCells_ = numDivisions_[0]*numDivisions_[1]*numDivisions_[2];

  return;
}

std::vector<size_t> const Grid::cellIndices(size_t const ind) const 
{
  std::vector<size_t> indices(3,0);
  size_t indx = ind;
  size_t const nxy = numDivisions_[0]*numDivisions_[1];

  indices[2] = indx/nxy;
  indx -= indices[2]*nxy;
  indices[1] = indx/numDivisions_[0];
  indx -= indices[1]*numDivisions_[0];
  indices[0] = indx;

  return indices;
}

