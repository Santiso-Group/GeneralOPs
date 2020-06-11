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
** Crystal order parameter calculator class
**
** This class calculates crystal order parameters for a given
** System<PointMolecule> and a given set of crystal distribution.
** parameters. The order parameters are calculated for a spatial 
** grid defined by the number of divisions of the unit cell in
** each direction. For example, if grid.x = grid.y = grid.z = 5, 
** there are 125 order parameters (5 divisions in each direction).
** If any of the numbers of divisions is zero, the class calculates
** order parameters per molecule instead of using the spatial grid.
**
** When the spatial grid is used, a switching function is used
** to avoid discontinuities when a molecule moves from one grid
** cell to another. The width of the switching function can be given
** between 0.0 (no switching function) and 0.5.
**
** It is also possible to define a cutoff for the determination of
** a molecule's neighbors. If this cutoff is set to zero, all molecules
** are included in the calculation.
*/

/*
** Some notes:
** 
** - This assumes that all pointmolecules on the given system are
**   of the same kind - the case of hybrid relative configurations
**   is not implemented yet (see note in crystdist.cpp)
*/

#ifndef H_ORDER_PARAMETER_CALCULATOR
#define H_ORDER_PARAMETER_CALCULATOR

#include <vector>
#include <cmath>
#include "mymol/include/lattice.h"
#include "mymol/include/system.h"
#include "crystdist/include/statparameters.h"
#include "crystdist/include/crystalops.h"

// Container to store the list of cells to which a molecule contributes
// order parameters and their weights. As with the grid object, the
// x, y, z directions really refer to crystal axes

struct CellData
{
  std::vector<size_t> cellIndices;  // Indices of the cells
  std::vector<Real> cellWeights;    // Weight in each cell

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

class OrderParameterCalculator
{
public:

// Constructors

  OrderParameterCalculator();                                                       // Define an empty OP calculator
  OrderParameterCalculator(System<PointMolecule> const &system,                     // Define an OP calculator with
                           CrystalDistributionParameters const &parameters,         // given system, set of statistical
                           OrderParameterGrid const &grid = OrderParameterGrid(),   // parameters and, optionally, the
                           Real const &switchWidth = 0.0,                           // grid, switching function width,
                           Real const &cutoff = 0.0,                                // distance cutoff, and whether
                           bool useLogs = false,                                    // to use the log of the OPs
                           Real const &pValue = 0.0);                               // or a p-norm for averaging.

// Interface

  void calculate(System<PointMolecule> const &system,               // Calculates crystal OPs for a given system
                 CrystalDistributionParameters const &parameters,   // and a given set of statistical parameters.
                 bool const useLogs = false,                        // If useLogs = true, stores the log of the OPs.
                 Real const &pValue = 0.0);                         // If pValue is not zero, uses p-norm.
  void setGrid(OrderParameterGrid const &grid);                     // Sets the order parameter grid (clears the ops)
  void setSwitchWidth(Real const &switchWidth);                     // Sets the width of the switching function
                                                                    // (clears the ops)
  void setCutoff(Real const &cutoff);                               // Sets the distance cutoff (clears the ops)
  void clearOPs();                                                  // Clears the order parameters only (not the grid)
  void clear();                                                     // Clears the order parameters and the grid

// Accessors

  OrderParameterGrid const &grid() const;                 // Returns the order parameter grid
  Real const &switchWidth() const;                        // Returns the width of the switching function
  Real const cutoff() const;                              // Returns the distance cutoff
  size_t const numCells() const;                          // Returns the number of grid cells
  CrystalOrderParameters const &orderParameters() const;  // Returns the order parameters

private:

  CrystalOrderParameters ops_;  // The order parameters

// Private functions

  CellData const cellData(Vector3D const &position, Lattice const &lattice);
};

/*
** End of class OrderParameterCalculator
*/

// Inlines

inline void OrderParameterCalculator::setGrid(OrderParameterGrid const &grid)
{
  ops_.grid = grid;
  ops_.clearOPs();
}

inline void OrderParameterCalculator::setSwitchWidth(Real const &switchWidth)
{
  if(switchWidth < 0.0 || switchWidth > 0.5)
  {
    std::cerr << "Error in OrderParameterCalculator::setSwitchWidth: "
              << "switching function width must be between 0 and 0.5" << std::endl;
    return;
  }
  ops_.switchWidth = switchWidth;
  ops_.clearOPs();
}

inline void OrderParameterCalculator::setCutoff(Real const &cutoff)
{
  if(cutoff < 0.0)
  {
    std::cerr << "Error in OrderParameterCalculator::setCutoff: "
              << "cutoff must not be negative" << std::endl;
    return;
  }
  ops_.cutoffSq = cutoff*cutoff;
  ops_.clearOPs();
}

inline void OrderParameterCalculator::clearOPs()
{
  ops_.clearOPs();
}

inline void OrderParameterCalculator::clear()
{
  ops_.clear();
}

inline OrderParameterGrid const &OrderParameterCalculator::grid() const
{
  return ops_.grid;
}

inline Real const &OrderParameterCalculator::switchWidth() const
{
  return ops_.switchWidth;
}

inline Real const OrderParameterCalculator::cutoff() const
{
  return sqrt(ops_.cutoffSq);
}

inline size_t const OrderParameterCalculator::numCells() const
{
  return ops_.grid.numCells();
}

inline CrystalOrderParameters const &OrderParameterCalculator::orderParameters() const
{
  return ops_;
}

#endif
