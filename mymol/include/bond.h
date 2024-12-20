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
** Simple bond class
**
** Note that all data members are public. Checking the consistency of
** the data (e.g. first, second < number of atoms in molecule) is the 
** responsibility of the containing class/function.
*/

#ifndef H_BOND
#define H_BOND

#include <iostream>
#include <string>
#include <vector>

class Bond
{
public:

  size_t            first;      // First atom in the bond
  size_t            second;     // Second atom in the bond
  std::string       type;       // Bond type
  std::vector<Real> parameters; // Container for bond parameters, if used

// Constructors
  
  Bond();                                         // Defines an empty bond
  Bond(size_t const &atom1, size_t const &atom2,  // Defines a bond of given type between atoms atom1 and atom2
       std::string const &type = "un",            // with optionally defined type and parameters
       std::vector<Real> const &parameters = std::vector<Real>());

// Print bond information (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, Bond const &bond); 
};

/*
** End of class Bond
*/

// Inlines

inline Bond::Bond()
: 
first(0), second(0), type("un"), parameters()
{}

inline Bond::Bond(size_t const &atom1, size_t const &atom2, 
                  std::string const &type,
                  std::vector<Real> const &parameters)
: 
first(atom1), second(atom2), type(type), parameters(parameters)
{}

inline std::ostream &operator<<(std::ostream &outStream, Bond const &bond)
{
  outStream << bond.first << " " << bond.second << " " << bond.type;
  for(size_t i = 0; i < bond.parameters.size(); ++i)
    outStream << " " << bond.parameters[i];
  return outStream;
}

#endif

