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
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "common/include/assert.h"
#include "common/include/iofile.h"
#include "mymol/include/atom.h"

// Constructors

Atom::Atom()
:
number(0), symbol(), name(), mass(0.0), radius(0.0), position(0.0, 0.0, 0.0), 
type("Any"), residueID(0), residueName(""), roleName(), tempFactor(0.0)
{}

Atom::Atom(std::string const &symbol, Vector3D const &position)
:
number(0), symbol(symbol), name(), mass(0.0), radius(0.0), position(position), 
type("Any"), residueID(0), residueName(""), roleName(), tempFactor(0.0)
{ setAtom(symbol); }

Atom::Atom(std::string const &symbol, Vector3D const &position, 
           size_t const residueID, std::string const &residueName)
:
number(0), symbol(symbol), name(), mass(0.0), radius(0.0), position(position), 
type("Any"), residueID(residueID), residueName(residueName), roleName(),
tempFactor(0.0)
{ setAtom(symbol); }

Atom::Atom(std::string const &symbol, Vector3D const &position, 
           size_t const residueID, std::string const &residueName,
           std::string const &roleName, Real const &tempFactor)
:
number(0), symbol(symbol), name(), mass(0.0), radius(0.0), position(position), 
type("Any"), residueID(residueID), residueName(residueName), roleName(roleName),
tempFactor(tempFactor)
{ setAtom(symbol); }

Atom::Atom(std::string const &symbol, std::string const &type, Vector3D const &position)
:
number(0), symbol(symbol), name(), mass(0.0), radius(0.0), position(position), 
type(type), residueID(0), residueName(""), roleName(), tempFactor(0.0)
{ setAtom(symbol); }

Atom::Atom(std::string const &symbol, std::string const &type, Vector3D const &position,
           size_t const residueID, std::string const &residueName)
:
number(0), symbol(symbol), name(), mass(0.0), radius(0.0), position(position), 
type(type), residueID(residueID), residueName(residueName), roleName(),
tempFactor(0.0)
{ setAtom(symbol); }

Atom::Atom(std::string const &symbol, std::string const &type, Vector3D const &position,
           size_t const residueID, std::string const &residueName, std::string const &roleName)
:
number(0), symbol(symbol), name(), mass(0.0), radius(0.0), position(position), 
type(type), residueID(residueID), residueName(residueName), roleName(roleName),
tempFactor(0.0)
{ setAtom(symbol); }

// Set properties using information from the atom database

void Atom::setAtom(std::string const &symbol)
{
  // Set atom properties by reading from the atom database (search by symbol)

  std::ifstream database; // Atom database
  size_t r_number;        // Atomic number read from database
  std::string r_symbol;   // Atomic symbol read from database
  std::string a_symbol;   // Symbol minus any trailing numbers (e.g. "C11" -> "C")
  std::string r_name;     // Atom name read from database
  Real r_mass;            // Atom mass read from database
  Real r_radius;          // Atom radius read from database

  size_t numIndex = symbol.find_first_of("0123456789");
  a_symbol = (numIndex != symbol.npos)?symbol.substr(0, numIndex):symbol;

  database.open(ATOM_DATABASE.c_str(), std::ios::in);  // Open database for reading

  if(!database) 
  {
    std::cerr << "Atom data file " << ATOM_DATABASE << " not found." << std::endl;
    return;
  };
  do 
  {
    database >> r_number >> r_symbol >> r_name >> r_mass >> r_radius;
    if(r_symbol == a_symbol) 
    {
      (*this).number = r_number;
      (*this).name = r_name;
      (*this).mass = r_mass;
      (*this).radius = r_radius;
      break;
    };
    database.ignore(LINE_LENGTH, '\n');
  }
  while(!database.eof());

  database.close();   // Close database

  if(database.eof())
  {
    std::cerr << "Atom " << a_symbol << " not found in database." << std::endl;
  }
}

// Output (mostly for debugging)

std::ostream &operator<<(std::ostream &outStream, Atom const &atom)
{
  size_t const width = 14;            // Width of output field
  std::string const line(width, '-'); // Line to separate output

  outStream.setf(std::ios::fixed);
  outStream.setf(std::ios::right, std::ios::adjustfield);

  return outStream << line << std::endl
                   << std::setw(width) << atom.number << std::endl
                   << std::setw(width) << atom.symbol << std::endl
                   << std::setw(width) << atom.name << std::endl
                   << std::setw(width) << atom.mass << std::endl
                   << std::setw(width) << atom.radius << std::endl
                   << line << std::endl
                   << "Atom type: " << atom.type << std::endl
                   << "Position: " << atom.position << std::endl
                   << "Residue ID number: " << atom.residueID << std::endl
                   << "Residue name: " << atom.residueName;
  
}
