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
** Geometry class
**
** Simple container for a set of atoms and bonds with a label (name).
*/

#ifndef H_GEOMETRY
#define H_GEOMETRY

#include <vector>
#include "mymath/include/vector3D.h"
#include "mymol/include/atom.h"
#include "mymol/include/bond.h"

Real const DEFAULT_BOND_FACTOR_SQUARED = 1.44;  // Square of default distance factor
                                                // for automatic determination of bonds

class Geometry 
{
public:

// Constructors

  Geometry(); // Defines an empty geometry

// Destructor

  virtual ~Geometry();  // For virtual derivation

// Interface

  void setName(std::string const &name);                    // Set the geometry label
  virtual void addAtom(Atom const &atom);                   // Adds an atom to the geometry
  virtual void addBond(Bond const &bond);                   // Adds a bond to the geometry
  virtual void setPosition(size_t const index,              // Set the position of the index-th atom
                           Vector3D const &newPosition);    // This does *not* recalculate bonds, if you
                                                            // need that call calculateBonds afterwards
  virtual void setx(size_t const index, Real const &newx);  // Set the x-coordinate of the index-th atom
                                                            // to newx. Does not recalculate bonds.
  virtual void sety(size_t const index, Real const &newy);  // Set the y-coordinate of the index-th atom
                                                            // to newx. Does not recalculate bonds.
  virtual void setz(size_t const index, Real const &newz);  // Set the z-coordinate of the index-th atom
                                                            // to newx. Does not recalculate bonds.
  virtual void clear();                                     // Clears the geometry
  void setTemperatureFactor(Real const &tempFactor);        // Sets the temperature factor of the whole
                                                            // geometry to the given value
  void setTemperatureFactor(size_t const index,             // Sets the temperature factor of the
                            Real const &tempFactor);        // index-th atom to the given value
  void setBondType(size_t const index,                      // Sets the type of the index-th bond
                   std::string const &type);
  void setBondParameters(size_t const index,                // Sets the parameters of the index-th bond
                         std::vector<Real> const &pars);

// Find bonds automatically based on interatomic distance. If the distance between two atoms is less than
// the sum of their covalent radii times the bond factor, they are considered bonded

  virtual void calculateBonds(Real const &bondFactorSquared = DEFAULT_BOND_FACTOR_SQUARED);

// Accessors

  std::string const &name() const;        // Returns the geometry label
  size_t const numAtoms() const;          // Returns the number of atoms
  size_t const numBonds() const;          // Returns the number of bonds
  Atom const &atom(size_t const i) const; // Returns the i-th atom
  Bond const &bond(size_t const i) const; // Returns the i-th bond

private:

  std::string name_;
  std::vector<Atom> geometry_;
  std::vector<Bond> bondTable_;
};

/*
** End of class Geometry
*/

// Inlines

inline void Geometry::setName(std::string const &name)
{
  name_ = name;
}

inline void Geometry::addAtom(Atom const &atom)
{
  geometry_.push_back(atom);
}

inline void Geometry::addBond(Bond const &bond)
{
  size_t numAtoms = geometry_.size();
  assert(bond.first < numAtoms && bond.second < numAtoms);
  bondTable_.push_back(bond);
}

inline void Geometry::setPosition(size_t const index, Vector3D const &newPosition)
{
  assert(index < geometry_.size());
  geometry_[index].position = newPosition;
}

inline void Geometry::setx(size_t const index, Real const &newx)
{
  assert(index < geometry_.size());
  geometry_[index].position.x = newx;
}

inline void Geometry::sety(size_t const index, Real const &newy)
{
  assert(index < geometry_.size());
  geometry_[index].position.y = newy;
}

inline void Geometry::setz(size_t const index, Real const &newz)
{
  assert(index < geometry_.size());
  geometry_[index].position.z = newz;
}

inline void Geometry::clear()
{
  name_.clear();
  geometry_.clear();
  bondTable_.clear();
}

inline void Geometry::setTemperatureFactor(Real const &tempFactor)
{
  for(size_t i = 0; i < geometry_.size(); ++i)
    geometry_[i].tempFactor = tempFactor;
}

inline void Geometry::setTemperatureFactor(size_t const index, Real const &tempFactor)
{
  assert(index < geometry_.size());
  geometry_[index].tempFactor = tempFactor;
}

inline void Geometry::setBondType(size_t const index, std::string const &type)
{
  assert(index < bondTable_.size());
  bondTable_[index].type = type;
}

inline void Geometry::setBondParameters(size_t const index, std::vector<Real> const &pars)
{
  assert(index < bondTable_.size());
  bondTable_[index].parameters = pars;
}

inline std::string const &Geometry::name() const
{
  return name_;
}

inline size_t const Geometry::numAtoms() const
{
  return geometry_.size();
}

inline size_t const Geometry::numBonds() const
{
  return bondTable_.size();
}

inline Atom const &Geometry::atom(size_t const i) const
{
  assert(i < geometry_.size());
  return geometry_[i];
}

inline Bond const &Geometry::bond(size_t const i) const
{
  assert(i < bondTable_.size());
  return bondTable_[i];
}

#endif

