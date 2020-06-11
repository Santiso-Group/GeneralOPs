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
** Point molecule class
**
** A "point molecule" is a representation of a molecule that is
** useful to study order in molecular crystals/liquid crystals. 
** The molecule is represented by:
**
** - A position - for now, this is the position of the center of mass.
**
** - An orientation - this is stored as a quaternion. There are three
**   ways to represent absolute orientations:
**
**   (1) For linear and symmetric molecules (D_infinity_h, e.g. H2 or
**       CO2), or planar molecules that have in-plane rotational
**       symmetry (e.g. benzene), the orientation is represented by
**       a directionless vector pointing along the molecule's axis
**       (linear) or the normal to the plane (planar).
**
**   (2) For linear but asymmetric molecules (C_infinity_v, e.g. CO),
**       the orientation is given by a vector pointing along the
**       axis of the molecule. This is treated the same as case (1)
**       within this class, but the orientational statistics will
**       be calculated differently in classes that depend on this one.
**
**   (3) For molecules that do not fall in cases (1) or (2), the 
**       orientation is defined by an orthogonal coordinate frame
**       defined by three non-collinear atoms in the molecule.
**
** - A set of internal degrees of freedom representing the internal
**   configuration of the molecule. These can be interatomic distances,
**   angles, or dihedral angles. For dihedrals there is also a
**   symmetry number in case there is rotational symmetry around the
**   central bond. For example, for X-COO (carboxylate) the symmetry
**   number would be 2, for X-CH3 it would be 3. This number has no
**   effect on this class, but can be used to do angle statistics.
**
** Note that all data members are public. A method is provided
** to create a point molecule from a molecule and a molecule map.
*/

#ifndef H_POINT_MOLECULE
#define H_POINT_MOLECULE

#include <vector>
#include <string>
#include "common/include/assert.h"
#include "common/include/types.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/quaternion.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"
#include "crystdist/include/moleculemap.h"

struct InternalDOF
{
  InternalDOFType type;   // Type of internal DOF (see moleculemap.h)
  Real value;             // Value of the internal DOF
  size_t symmetryNumber;  // Symmetry number (only meaningful for dihedrals)

  InternalDOF()
  :
  type(DISTANCE), value(0.0), symmetryNumber(1)
  {}

  InternalDOF(InternalDOFType const &type, size_t const symmetryNumber = 1)
  :
  type(type), value(0.0), symmetryNumber(symmetryNumber)
  {}

  InternalDOF(InternalDOFType const &type, Real const &value, size_t const symmetryNumber = 1)
  :
  type(type), value(value), symmetryNumber(symmetryNumber)
  {}
};

class PointMolecule
{
public:

  std::string name;                      // Name of the point molecule (e.g. residue ID)
  Vector3D position;                     // Position of the point molecule
  Quaternion orientation;                // Absolute orientation of the point molecule
  MoleculeType type;                     // Type of point molecule (see moleculemap.h)
  std::vector<InternalDOF> internalDOFs; // Internal degrees of freedom

// Constructors

  PointMolecule(MoleculeType const &type = GENERAL); // Defines an empty point molecule with an
                                                     // optionally given type
  PointMolecule(Molecule &molecule,                  // Defines a point molecule from a molecule
                MoleculeMap const &map);             // and a molecule map

// Interface

  void set(Molecule &molecule, MoleculeMap const &map); // Defines a point molecule from a molecule
                                                        // and a molecule map
  void clear();                                         // Clears the point molecule

// Accessors

  size_t const numInternalDOFs() const;                     // Returns the number of internal degrees
                                                            // of freedom
  InternalDOF const &internalDOF(size_t const index) const; // Returns the index-th internal degree of freedom

// Build a System<PointMolecule> from a System<Molecule> and a molecule map. Only the molecules with
// the same residue name as the ones in the molecule map vector are converted, others are ignored.

  friend System<PointMolecule> const getPointMolecules(System<Molecule> &system, 
                                                       MoleculeMap const &map);

// Print point molecule information (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, PointMolecule const &pointMolecule); 
};

/*
** End of class PointMolecule
*/

// Inlines

inline void PointMolecule::clear()
{
  name.clear();
  position = 0.0;
  orientation = 1.0;
  type = GENERAL;
  internalDOFs.clear();
}

inline size_t const PointMolecule::numInternalDOFs() const
{
  return internalDOFs.size();
}

inline InternalDOF const &PointMolecule::internalDOF(size_t const index) const
{
  assert(index < internalDOFs.size());
  return internalDOFs[index];
}

// Prototype for getPointMolecules - otherwise gcc pukes
System<PointMolecule> const getPointMolecules(System<Molecule> &system, 
                                              MoleculeMap const &map);

#endif

