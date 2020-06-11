/*
** Copyright 2008-2011 Erik Santiso.
** This file is part of crystdist.
** crystdist is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** crystdist is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with crystdist. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** Crystal order parameter containers
**
** Notes: 
**
** - There's no special reason why the order parameter grid and/or the
**   switching function should be here - maybe move to mymol/mymath?
** - Right now the CrystalOrderParameters struct is not safe - it is 
**   the responsibility of the programmer to test that values make sense,
**   e.g. resize does not check that the size is the same as the grid size.
*/

#ifndef H_CRYSTAL_ORDER_PARAMETERS
#define H_CRYSTAL_ORDER_PARAMETERS

#include <vector>
#include "mymath/include/vector3D.h"
#include "mymol/include/lattice.h"

// Grid to calculate order parameters

// Note: If any of the numbers of divisions is zero, no grid is used.
// Instead, the OPs will be calculated per molecule (see opcalculator.cpp)

struct OrderParameterGrid
{
  size_t x;         // Number of divisions in the x-direction
  size_t y;         // Number of divisions in the y-direction
  size_t z;         // Number of divisions in the z-direction
  Vector3D origin;  // Origin of the unit cell (for compatibility
                    // with NAMD)

  OrderParameterGrid()
  :
  x(1), y(1), z(1), origin(0.0)
  {}

  OrderParameterGrid(size_t const x, size_t const y, size_t const z,
                     Vector3D const &origin = 0.0)
  :
  x(x), y(y), z(z), origin(origin)
  {}

  inline void clear()
  {
    x = y = z = 1; origin = 0.0;
  }

  inline size_t const numCells() const
  {
    return x*y*z;
  }

  // Return the grid index corresponding to a given position and unit cell
  inline size_t const getIndex(Vector3D const &position, Lattice const &lattice) const
  {
    Vector3D const pos = position - origin; // For compatibility with NAMD
    Vector3D const scaled(lattice.reciprocalVector(0)*pos,
                          lattice.reciprocalVector(1)*pos,
                          lattice.reciprocalVector(2)*pos);
    size_t ix = (size_t)floor(x*(scaled.x + 0.5 - floor(scaled.x + 0.5)));
    size_t iy = (size_t)floor(y*(scaled.y + 0.5 - floor(scaled.y + 0.5)));
    size_t iz = (size_t)floor(z*(scaled.z + 0.5 - floor(scaled.z + 0.5)));
    if(ix > x - 1) ix = x - 1;
    if(iy > y - 1) iy = y - 1;
    if(iz > z - 1) iz = z - 1;
    return ix + (iy + y*iz)*x;
  }
  
  // Input
  friend std::istream &operator>>(std::istream &inStream, OrderParameterGrid &grid)
  {
    inStream >> grid.x >> grid.y >> grid.z
             >> grid.origin.x >> grid.origin.y >> grid.origin.z;
    return inStream;
  }

  // Output
  friend std::ostream &operator<<(std::ostream &outStream, OrderParameterGrid const &grid)
  {
    outStream << grid.x << " " << grid.y << " " << grid.z << " " 
              << grid.origin.x << " " << grid.origin.y << " " << grid.origin.z;
    return outStream;
  }
};

// Container to store order parameters

struct CrystalOrderParameters
{
  OrderParameterGrid grid;                  // The order parameter grid
  Real switchWidth;                         // Width of switching function
  Real cutoffSq;                            // Square of distance cutoff (0 = none)
  std::vector<Real> distance;               // Distance order parameters
  std::vector<Real> bondOrientation;        // Bond orientation order parameters
  std::vector<Real> relativeOrientation;    // Relative orientation order parameters
  std::vector<std::vector<Real> > internal; // Internal DOF order parameters
  std::vector<Real> localDensity;           // Local densities
  std::vector<Real> total;                  // Total order parameters

  CrystalOrderParameters()
  :
  grid(), switchWidth(0.0), cutoffSq(0.0), distance(), bondOrientation(), 
  relativeOrientation(), internal(), localDensity(), total()
  {}
  
  // Clear the order parameters, but not the grid or the switching width
  void clearOPs()
  {
    distance.clear();
    bondOrientation.clear();
    relativeOrientation.clear();
    internal.clear();
    localDensity.clear();
    total.clear();
  }

  // Clear the order parameters and the grid
  void clear()
  {
    switchWidth = 0.0;
    cutoffSq = 0.0;
    grid.clear();
    distance.clear();
    bondOrientation.clear();
    relativeOrientation.clear();
    internal.clear();
    localDensity.clear();
    total.clear();
  }

  // Resize the order parameter arrays
  void resize(size_t const size, size_t const numInternalDOFs = 0)
  {
    distance.resize(size, 0.0);
    bondOrientation.resize(size, 0.0);
    relativeOrientation.resize(size, 0.0);
    internal.resize(size, std::vector<Real>(numInternalDOFs, 0.0));
    localDensity.resize(size, 0.0);
    total.resize(size, 0.0);
  }

  // Size (useful to treat as a vector, e.g. for the string method)
  size_t const size() const
  {
    return distance.size() + bondOrientation.size() +
           relativeOrientation.size() + internal.size()*internal[0].size() +
           localDensity.size() + total.size();
  }

  // Access by index (useful for string method, etc.)
  Real const &operator[](size_t const i) const
  {
    if(i < distance.size())
      return distance[i];
    else if(i < distance.size() + bondOrientation.size())
      return bondOrientation[i - distance.size()];
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size())
      return relativeOrientation[i - distance.size() - bondOrientation.size()];
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size() +
                internal.size()*internal[0].size())
    {
      size_t iM = i - distance.size() - bondOrientation.size() - relativeOrientation.size();
      size_t nInt = internal[0].size();
      if(iM/nInt >= internal.size())
      {
        std::cerr << "Index out of bounds in CrystalOrderParameters::operator[]" << std::endl;
        return distance[0];
      }
      else return internal[iM/nInt][iM%nInt];
    }
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size() +
                internal.size()*internal[0].size() + localDensity.size())
      return localDensity[i - distance.size() - bondOrientation.size() - relativeOrientation.size() -
                          internal.size()*internal[0].size()];
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size() +
                internal.size()*internal[0].size() + localDensity.size() + total.size())
      return total[i - distance.size() - bondOrientation.size() - relativeOrientation.size() -
                   internal.size()*internal[0].size() - localDensity.size()];
    else
    {
      std::cerr << "Index out of bounds in CrystalOrderParameters::operator[]" << std::endl;
      return distance[0];
    }
  }

  // Access by index (non-const version)
  Real &operator[](size_t const i)
  {
    if(i < distance.size())
      return distance[i];
    else if(i < distance.size() + bondOrientation.size())
      return bondOrientation[i - distance.size()];
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size())
      return relativeOrientation[i - distance.size() - bondOrientation.size()];
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size() +
                internal.size()*internal[0].size())
    {
      size_t iM = i - distance.size() - bondOrientation.size() - relativeOrientation.size();
      size_t nInt = internal[0].size();
      if(iM/nInt >= internal.size())
      {
        std::cerr << "Index out of bounds in CrystalOrderParameters::operator[]" << std::endl;
        return distance[0];
      }
      else return internal[iM/nInt][iM%nInt];
    }
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size() +
                internal.size()*internal[0].size() + localDensity.size())
      return localDensity[i - distance.size() - bondOrientation.size() - relativeOrientation.size() -
                          internal.size()*internal[0].size()];
    else if(i < distance.size() + bondOrientation.size() + relativeOrientation.size() +
                internal.size()*internal[0].size() + localDensity.size() + total.size())
      return total[i - distance.size() - bondOrientation.size() - relativeOrientation.size() -
                   internal.size()*internal[0].size() - localDensity.size()];
    else
    {
      std::cerr << "Index out of bounds in CrystalOrderParameters::operator[]" << std::endl;
      return distance[0];
    }
  }

  // Basic vector arithmetic (useful to treat OPs as vectors, e.g. string method)
  CrystalOrderParameters &operator+=(CrystalOrderParameters const &c)
  {
    size_t const numParameters = size();
    assert(numParameters == c.size());
    for(size_t i = 0; i < numParameters; ++i)
      (*this)[i] += c[i];
    return *this;
  }

  CrystalOrderParameters &operator-=(CrystalOrderParameters const &c)
  {
    size_t const numParameters = size();
    assert(numParameters == c.size());
    for(size_t i = 0; i < numParameters; ++i)
     (*this)[i] -= c[i];
    return *this;
  }

  CrystalOrderParameters &operator*=(Real const &r)
  {
    for(size_t i = 0; i < size(); ++i)
      (*this)[i] *= r;
    return *this;
  }

  CrystalOrderParameters &operator/=(Real const &r)
  {
    assert(r);
    Real const invr = 1.0/r;
    for(size_t i = 0; i < size(); ++i)
     (*this)[i] *= invr;
    return *this;
  }

  friend CrystalOrderParameters operator+(CrystalOrderParameters const &c1, CrystalOrderParameters const &c2)
  {
    CrystalOrderParameters sum = c1;
    sum += c2;
    return sum;
  }

  friend CrystalOrderParameters operator-(CrystalOrderParameters const &c1, CrystalOrderParameters const &c2)
  {
    CrystalOrderParameters diff = c1;
    diff -= c2;
    return diff;
  }

  friend CrystalOrderParameters operator*(Real const &r, CrystalOrderParameters const &c)
  {
    CrystalOrderParameters prod = c;
    prod *= r;
    return prod;
  }

  friend CrystalOrderParameters operator*(CrystalOrderParameters const &c, Real const &r)
  {
    CrystalOrderParameters prod = c;
    prod *= r;
    return prod;
  }

  // Output (mostly for debugging)
  friend std::ostream &operator<<(std::ostream &outStream, CrystalOrderParameters const &ops)
  {
    outStream << "Order parameter grid: " << std::endl;
    outStream << ops.grid << std::endl;
    outStream << "Distance order parameters: " << std::endl;
    for(size_t i = 0; i < ops.distance.size(); ++i)
      outStream << ops.distance[i] << std::endl;
    outStream << "Bond orientation order parameters: " << std::endl;
    for(size_t i = 0; i < ops.bondOrientation.size(); ++i)
      outStream << ops.bondOrientation[i] << std::endl;
    outStream << "Relative orientation order parameters: " << std::endl;
    for(size_t i = 0; i < ops.relativeOrientation.size(); ++i)
      outStream << ops.relativeOrientation[i] << std::endl;
    outStream << "Internal order parameters: " << std::endl;
    for(size_t i = 0; i < ops.internal.size(); ++i)
    {
      for(size_t j = 0; j < ops.internal[i].size(); ++j)
        outStream << ops.internal[i][j] << " ";
      outStream << std::endl;
    }
    outStream << "Local densities: " << std::endl;
    for(size_t i = 0; i < ops.localDensity.size(); ++i)
      outStream << ops.localDensity[i] << std::endl;
    outStream << "Total order parameters: " << std::endl;
    for(size_t i = 0; i < ops.total.size(); ++i)
      outStream << ops.total[i] << std::endl;
    return outStream;
  }
};

#endif
