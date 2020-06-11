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
** Simple atom class
**
** Note that all data members are public. Checking the consistency of
** the data (e.g. mass > 0, radius > 0, etc.) is the responsibility of
** the containing class/function.
**
** More format-specific properties (listed under "Other properties" below)
** can be added to handle various file formats.
*/

#ifndef H_ATOM
#define H_ATOM

#include <iostream>
#include <string>
#include "common/include/types.h"
#include "mymath/include/vector3D.h"
#include "mymol/include/datafiles.h" // Contains the name of the atom database file

class Atom  
{
public:

// Atom properties

  size_t        number;   // Atomic number
  std::string   symbol;   // Atom symbol
  std::string   name;     // Atom name
  Real          mass;     // Atom mass
  Real          radius;   // Covalent radius
  Vector3D      position; // Position

// Other properties - add more as needed

  std::string   type;         // Atom type (for force field-specific formats, e.g. mol2)
  size_t        residueID;    // Residue ID number
  std::string   residueName;  // Residue name
  std::string   roleName;     // Atom role name (used by pdb, psf formats)
  Real          tempFactor;   // Temperature factor (used by pdb)

// Constructors
//
// If only an atomic number or symbol is given, the name, mass and radius are read
// from the atom database #defined above.
//
// More constructors can be added as needed to read from various file formats or to 
// keep track of more information than the current ones.

// General:

  Atom();                           // Defines an empty atom at the origin
  Atom(std::string const &symbol,   // Defines an atom by atomic symbol 
       Vector3D const &position);   // at a given position

// For formats containing residue information

  Atom(std::string const &symbol,       // Defines an atom by symbol with
       Vector3D const &position,        // given position, residue id
       size_t const residueID,          // number and residue name
       std::string const &residueName);

  Atom(std::string const &symbol,       // Defines an atom by symbol with
       Vector3D const &position,        // given position, residue id
       size_t const residueID,          // number, residue name, role name
       std::string const &residueName,  // and temperature factor(e.g. pdb)
       std::string const &roleName,
       Real const &tempFactor);

// For force field-specific formats:

  Atom(std::string const &symbol,       // Defines an atom by symbol with
       std::string const &type,         // given atom type and position
       Vector3D const &position);

  Atom(std::string const &symbol,       // Defines an atom by symbol with
       std::string const &type,         // given atom type, position,
       Vector3D const &position,        // residue id number and residue
       size_t const residueID,          // name
       std::string const &residueName);

  Atom(std::string const &symbol,       // Defines an atom by symbol with
       std::string const &type,         // given atom type, position,
       Vector3D const &position,        // residue id number, residue
       size_t const residueID,          // name and role name
       std::string const &residueName,
       std::string const &roleName);

// Set properties using information from the atom database

  void setAtom(std::string const &symbol);  // Set by symbol

// Print atom information (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, Atom const &atom); 
};

/*
** End of class Atom
*/

#endif

