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
** Grid class
*
** This class implements a regular grid of up to three dimensions
** that can be used to calculate or interpolate properties on grid 
** cells in periodic systems.
**
** By convention, the grid is applied in the directions where the
** lattice has a non-zero lattice vector, the other components of a
** given vector are ignored. The class gives a warning if the number
** of divisions in a non-periodic direction is different from one.
**
** Also, the grid is defined so that the origin is at its center,
** and the cells are counted starting from the corner with the
** smallest scaled coordinates.
**
** For a later version if may be helpful to implement a more general
** grid (rectilinear or not).
*/

#ifndef H_GRID
#define H_GRID

#include <iostream>
#include <vector>
#include "common/include/types.h"
#include "mymath/include/vector3D.h"
#include "mymol/include/lattice.h"

class Grid
{
public:

// Constructors

  Grid();                                                 // Defines a grid with default parameters
  Grid(size_t const n, Lattice const &lattice,            // Defines a grid with equal number of divisions
       Vector3D const &origin = 0.0);                     // in all directions, a lattice and, optionally,
                                                          // an origin.
  Grid(size_t const nx, size_t const ny, size_t const nz, // Defines a grid given the numbers of divisions
       Lattice const &lattice,                            // in all directions, a lattice and, optionally,
       Vector3D const &origin = 0.0);                     // an origin.

// Interface

  void setLattice(Lattice const &lattice);               // Set the supercell
  void setOrigin(Vector3D const &origin);                // Set the origin of the grid
  void setNumDivisions(size_t const n);                  // Set the numbers of divisions in all directions
                                                         // to the given value
  void setNumDivisions(size_t const nx, size_t const ny, // Set the numbers of divisions in each direction
                       size_t const nz);
  void clear();                                          // Resets the grid to default parameters

// Accessors

  Lattice const &lattice() const;                   // Returns the lattice
  size_t const numDivisions(size_t const i) const;  // Returns the number of divisions in the i-th direction
  size_t const numPoints(size_t const i) const;     // Returns the number of points in the i-th direction
  size_t const numPoints() const;                   // Returns the total number of points
  size_t const numCells() const;                    // Returns the total number of cells
  Vector3D const &origin() const;                   // Returns the origin

// Other properties

  size_t const cellIndex(Vector3D const position) const;         // Returns the index of the cell where the given position lies
  Vector3D const point(size_t const i) const;                    // Returns the coordinates of the i-th grid point
  Vector3D const cellCenter(size_t const i) const;               // Returns the coordinates of the center of the i-th cell
  std::vector<size_t> const adjacentCells(size_t const i) const; // Returns the (unsorted) indices of the cells adjacent to the i-th cell

// Print grid information (mostly for debugging, can be long)

  friend std::ostream &operator<<(std::ostream &outStream, Grid const &grid);

private:

// Functions

  void updateGrid();                                                // Updates the numbers of cells and points and checks for consistency
  std::vector<size_t> const cellIndices(size_t const ind) const;    // Convert a global cell index to indices in x, y and z directions
  size_t const cellIndex(std::vector<size_t> const indices) const ; // Convert indices in x, y and z directions to global cell index

// Data members

  size_t numDivisions_[3]; // Number of divisions in each direction
  size_t numPoints_;       // Total number of points
  size_t numCells_;        // Total number of cells
  Lattice lattice_;        // The supercell
  Vector3D origin_;        // Origin of the grid
};

/*
** End of class Grid
*/ 

// Inlines

inline void Grid::setLattice(Lattice const &lattice)
{
  lattice_ = lattice;
  updateGrid();
}

inline void Grid::setOrigin(Vector3D const &origin)
{
  origin_ = origin;
}

inline void Grid::setNumDivisions(size_t const n)
{
  assert(n > 0);
  numDivisions_[0] = numDivisions_[1] = numDivisions_[2] = n;
  updateGrid();
}

inline void Grid::setNumDivisions(size_t const nx, size_t const ny, size_t const nz)
{
  assert(nx*ny*nz > 0);
  numDivisions_[0] = nx; numDivisions_[1] = ny; numDivisions_[2] = nz;
  updateGrid();
}

inline void Grid::clear()
{
  numDivisions_[0] = numDivisions_[1] = numDivisions_[2] = 1;
  lattice_.clear();
  origin_ = 0.0;
  updateGrid();
}

inline Lattice const &Grid::lattice() const
{
  return lattice_;
}

inline size_t const Grid::numDivisions(size_t const i) const
{
  assert(i < 3);
  return numDivisions_[i];
}

inline size_t const Grid::numPoints(size_t const i) const
{
  assert(i < 3);
  return numDivisions_[i] + 1;
}

inline size_t const Grid::numPoints() const
{
  return numPoints_;
}

inline size_t const Grid::numCells() const
{
  return numCells_;
}

inline Vector3D const &Grid::origin() const
{
  return origin_;
}

inline size_t const Grid::cellIndex(std::vector<size_t> const indices) const
{
  assert(indices.size() == 3);
  return indices[0] + numDivisions_[0]*(indices[1] + numDivisions_[1]*indices[2]);
}

#endif

