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
** Relative configuration class
**
** This is the basic object that crystdist does statistics on.
** A relative configuration is a set of two sets of internal degrees
** of freedom (one for each molecule in the pair), a distance,
** a bond orientation (on the coordinate frame of molecule 1), and
** a relative orientation (of molecule 2 with respect to molecule 1).
**
** Note that all data members are public. A method is provided to
** calculate a relative configuration from a pair of point molecules
** (all quantities are relative to the molecule frame of the first
** molecule).
*/

/*
** Some notes:
**
** - When the two objects are vector-like, the bond orientation and
**   relative orientation are simply the angles (in radians) between
**   the molecule-centered "z"-axis and the bond vector/the other
**   molecule's "z"-axis. These are stored in the real part of the
**   corresponding quaternions. 
**
** - Defining the relative configuration for a hybrid pair (e.g.
**   GENERAL - LINEAR_ASYMMETRIC) is not yet implemented. The bond
**   orientation can be defined just like in the symmetric case, i.e. 
**   use q_i_*r_ij*q_i for the first, and an angle for the second.
**   The problem is the relative orientation. Two options are:
**   (a) Use q_i_*u_j*q_i - the problem here is that the antisymmetry
**   between q_ij and q_ji is lost.
**   (b) Choose an arbitrary axis (z?) in the GENERAL object, and
**   do it as the case for two vectors/directions. This has the
**   right mathematical symmetry but loses information.
*/

#ifndef H_RELATIVE_CONFIGURATION
#define H_RELATIVE_CONFIGURATION

#include <string>
#include <vector>
#include <utility>
#include "common/include/assert.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/quaternion.h"
#include "mymol/include/lattice.h"
#include "crystdist/include/pointmolecule.h"

typedef std::pair< std::vector<InternalDOF>, std::vector<InternalDOF> > InternalDOFData;
typedef std::pair< MoleculeType, MoleculeType > TypeData;

class RelativeConfiguration
{
public:

  TypeData types;                 // Types of the original point molecules
  std::string name;               // Name (built from the point molecules' names)
  InternalDOFData internalDOFs;   // Internal degrees of freedom
  Real distance;                  // Distance
  Quaternion bondOrientation;     // Bond orientation
  Quaternion relativeOrientation; // Relative orientation

// Constructors

  RelativeConfiguration();                            // Define an empty relative configuration
  RelativeConfiguration(PointMolecule const &first,   // Define from a pair of point molecules
                        PointMolecule const &second,  // and, optionally, a lattice
                        Lattice const &lattice = Lattice());

// Interface

  void calculate(PointMolecule const &first,    // Calculate the relative configuration from a pair
                 PointMolecule const &second,   // of point molecules and, optionally, a lattice
                 Lattice const &lattice = Lattice());
  void clear();                                 // Clears the relative configuration

// Length squared of difference (compares two relative configurations)

  friend Real const squareDeviation(RelativeConfiguration const &conf1,
                                    RelativeConfiguration const &conf2);

// Comparison (used for sorting)

  friend bool operator<(RelativeConfiguration const &conf1, RelativeConfiguration const &conf2);

// I/O for relative configurations
 
  friend std::istream &operator>>(std::istream &inStream, RelativeConfiguration &conf);
  friend std::ostream &operator<<(std::ostream &outStream, RelativeConfiguration const &conf); 
};

/*
** End of class RelativeConfiguration
*/

// Inlines

inline void RelativeConfiguration::clear()
{
  types = TypeData(GENERAL, GENERAL);
  name.clear();
  internalDOFs = InternalDOFData(std::vector<InternalDOF>(), std::vector<InternalDOF>());
  distance = 0.0;
  bondOrientation = 0.0;
  relativeOrientation = 1.0;
}

#endif
