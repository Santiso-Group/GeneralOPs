/*
** Copyright 2007-2011 Erik Santiso.
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
** Simple lattice class
**
** Note about lattice types:
**
** The latticeType variable controls how the minimum image convention is applied.
** If set to SUPERCELL, the standard minimum image convention is used. This is
** correct only if the maximum distance is less than half of the smallest cell
** vector's length (which, for example, would be the case in a standard MD simulation).
** If set to CRYSTAL, an exhaustive search over cells is done and the closest image
** is chosen. This is, of course, much slower, and should only be used when distances
** can be larger than half of the smallest cell vector's length (which can happen, e.g.
** when analyzing a crystal structure).
*/

#ifndef H_LATTICE
#define H_LATTICE

#include <vector>
#include <string>
#include <iostream>
#include "common/include/assert.h"
#include "common/include/types.h"
#include "mymath/include/vector3D.h"

enum LatticeType { CRYSTAL, SUPERCELL };

class Lattice 
{
public:

// Constructors

  Lattice(LatticeType const &latticeType = SUPERCELL);  // Define a non-periodic lattice. Optionally, set
                                                        // the lattice type
  Lattice(Vector3D const &a,                            // Define a 1D periodic lattice by giving the lattice
          LatticeType const &latticeType = SUPERCELL);  // vector and, optionally, set the lattice type
  Lattice(Vector3D const &a1,                           // Define a 2D periodic lattice by giving the lattice
          Vector3D const &a2,                           // vectors and, optionally, set the lattice type
          LatticeType const &latticeType = SUPERCELL);
  Lattice(Vector3D const &a1,                           // Define a 3D periodic lattice by giving the lattice
          Vector3D const &a2,                           // vectors and, optionally, set the lattice type
          Vector3D const &a3,
          LatticeType const &latticeType = SUPERCELL);
  Lattice(Real const &a,                                // Define a lattice by giving the lengths of the
          Real const &b,                                // lattice vectors, and, optionally the angles between
          Real const &c,                                // them (in degrees) and the lattice type. If any
          Real const &alpha = 90.0,                     // length is set to 0, the cell is not periodic in
          Real const &beta = 90.0,                      // that direction (but angles should be consistent).
          Real const &gamma = 90.0,                     // For 1D, a is aligned with x-axis, b with y-axis,
          LatticeType const &latticeType = SUPERCELL);  // and c with z-axis. For other cases, a is aligned
                                                        // with x-axis and b is in the xy-plane.

// Interface

  void setType(LatticeType const &latticeType);                 // Sets the lattice type only
  void setLattice(Vector3D const &a,                            // Makes the lattice 1D periodic, sets the
                  LatticeType const &latticeType = SUPERCELL);  // lattice vector and, optionally, the lattice type
  void setLattice(Vector3D const &a1,                           // Makes the lattice 2D periodic, sets the
                  Vector3D const &a2,                           // lattice vectors and, optionally, the lattice type
                  LatticeType const &latticeType = SUPERCELL);
  void setLattice(Vector3D const &a1,                           // Makes the lattice 3D periodic, sets the 
                  Vector3D const &a2,                           // lattice vectors and, optionally, the lattice type
                  Vector3D const &a3, 
                  LatticeType const &latticeType = SUPERCELL);
  void setLattice(Real const &a,                                // Sets the lengths of the lattice vectors,
                  Real const &b,                                // and, optionally the angles between
                  Real const &c,                                // them (in degrees) and the lattice type.
                  Real const &alpha = 90.0,                     // If any length is set to 0, the cell is not
                  Real const &beta = 90.0,                      // periodic in that direction (but angles should
                  Real const &gamma = 90.0,                     // be consistent)
                  LatticeType const &latticeType = SUPERCELL);
  void clear(LatticeType const &latticeType = SUPERCELL);       // Clears the lattice (system is not periodic)
                                                                // Optionally, sets the lattice type

// Accessors

  Vector3D const &latticeVector(size_t const i) const;    // Returs the i-th lattice vector
  Vector3D const &reciprocalVector(size_t const i) const; // Returns the i-th reciprocal lattice vector
  LatticeType const &type() const;                        // Returns the lattice type

// Other functions

  Vector3D const difference(Vector3D const &v1,         // Returns the vector pointing from v1 to the
                            Vector3D const &v2) const;  // closest image of v2 (equivalent to v2 - v1)
  size_t const numDimensions() const;                   // Returns the number of dimensions along which the
                                                        // lattice is periodic
  Real const a() const;                                 // Returns the length of the first lattice vector. 
                                                        // If not defined, returns zero
  Real const b() const;                                 // Returns the length of the second cell vector.
                                                        // If not defined, returns zero
  Real const c() const;                                 // Returns the length of the third cell vector.
                                                        // If not defined, returns zero
  Real const alpha() const;                             // Returns the angle (in degrees) between the second and 
                                                        // third lattice vectors. If not defined, returns zero
  Real const beta() const;                              // Returns the angle (in degrees) between the first and 
                                                        // third lattice vectors. If not defined, returns zero
  Real const gamma() const;                             // Returns the angle (in degrees) between the first and 
                                                        // second lattice vectors. If not defined, returns zero
  Real const length() const;                            // Returns the length of the unit cell (1D)
  Real const area() const;                              // Returns the area of the unit cell (2D)
  Real const volume() const;                            // Returns the volume of the unit cell (3D)
  bool const isPeriodic() const;                        // Whether the unit cell is periodic or not

// Print lattice information (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, Lattice const &lattice); 

private:

// Functions

  void calculateLatticeVectors(Real const &a,       // Calculates the cell vectors from
                               Real const &b,       // the crystallographic constants
                               Real const &c,
                               Real const &alpha,
                               Real const &beta,
                               Real const &gamma);  
  void updateLattice();                             // Updates reciprocal lattice vectors and image positions
                                                    // after changing cell parameters or the lattice type.

// Data members

  LatticeType latticeType_;                 // Lattice type (see note at the top of this file)
  std::vector<Vector3D> latticeVectors_;    // Lattice vectors
  std::vector<Vector3D> reciprocalVectors_; // Reciprocal lattice vectors
  std::vector<Vector3D> imagePositions_;    // Positions of the periodic images in adjacent cells
};

/*
** End of class Lattice
*/

#endif

