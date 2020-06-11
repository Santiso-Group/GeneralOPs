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
** System class
**
** A container for multiple geometry-like objects and a unit cell,
** with a label (name).
**
** This represents a set of molecules in a box.
**
** The "Type" defined below is intended to be an object of class
** Geometry or similar (e.g. Molecule).
*/

#ifndef H_SYSTEM
#define H_SYSTEM

#include <vector>
#include <string>
#include "mymath/include/vector3D.h"
#include "mymol/include/lattice.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"

template<typename Type>
class System
{
public:

// Constructor

  System(); // Defines an empty system

// Interface

  void add(Type const &object);                 // Adds an object (e.g. molecule) to the system
  void setLattice(Lattice const &lattice);      // Sets the unit cell
  void setLatticeType(LatticeType const &type); // Sets the unit cell type only
  void setName(std::string const &name);        // Sets the name
  void clear();                                 // Clears the name and the objects in the system
  void reset();                                 // Clears the name, the objects, and the unit cell

// Access by element

  Type &operator[](size_t const i);             // Can modify
  Type const &operator[](size_t const i) const; // Cannot modify

// Accessors

  Type const &molecule(size_t const i) const; // Returns the i-th object
  Lattice const &lattice() const;             // Returns the unit cell
  std::string const &name() const;            // Returns the name
  size_t const size() const;                  // Returns the number of objects in the system

// Utilities

  void join(size_t const i); // Joins atoms from i-th object if they've been separated by wrapping
  void joinAll();            // Joins atoms from all objects if they've been separated by wrapping

private:

  std::vector<Type> system_;  // The objects in the system
  Lattice lattice_;           // The unit cell
  std::string name_;          // Name of the system
};

/*
** End of class System
*/

// Inlines

template<typename Type>
inline System<Type>::System()
:
system_(), lattice_(), name_("")
{}

template<typename Type>
inline void System<Type>::add(Type const &object)
{
  system_.push_back(object);
}

template<typename Type>
inline void System<Type>::setLattice(Lattice const &lattice)
{
  lattice_ = lattice;
}

template<typename Type>
inline void System<Type>::setLatticeType(LatticeType const &type)
{
  lattice_.setType(type);
}

template<typename Type>
inline void System<Type>::setName(std::string const &name)
{
  name_ = name;
}

template<typename Type>
inline void System<Type>::clear()
{
  name_.clear();
  system_.clear();
}

template<typename Type>
inline void System<Type>::reset()
{
  name_.clear();
  system_.clear();
  lattice_.clear();
}

template<typename Type>
inline Type &System<Type>::operator[](size_t const i)
{
  assert(i < system_.size());
  return system_[i];
}

template<typename Type>
inline Type const &System<Type>::operator[](size_t const i) const
{
  assert(i < system_.size());
  return system_[i];
}

template<typename Type>
inline Type const &System<Type>::molecule(size_t const i) const
{
  assert(i < system_.size());
  return system_[i];
}

template<typename Type>
inline Lattice const &System<Type>::lattice() const
{
  return lattice_;
}

template<typename Type>
inline std::string const &System<Type>::name() const
{
  return name_;
}

template<typename Type>
inline size_t const System<Type>::size() const
{
  return system_.size();
}

template<typename Type>
inline void System<Type>::join(size_t const i)
{
  assert(i < system_.size());
  Vector3D const pos = system_[i].atom(0).position;
  for(size_t iAtom = 0; iAtom < system_[i].numAtoms(); ++iAtom)
  {
    Vector3D diff = lattice_.difference(pos, system_[i].atom(iAtom).position);
    system_[i].setPosition(iAtom, pos + diff);
  }
}

template<typename Type>
inline void System<Type>::joinAll()
{
  for(size_t i = 0; i < system_.size(); ++i) join(i);
}

#endif

