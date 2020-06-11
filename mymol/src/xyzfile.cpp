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
** .xyz file format
*/

#include <iomanip>
#include "common/include/assert.h"
#include "mymol/include/file_formats/fileformats.h"

// Constructors

XYZFile::XYZFile()
:
IOFile()
{}

XYZFile::XYZFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

XYZFile &operator<<(XYZFile &xyzFile, Geometry const &geometry)
{
  // Write geometry to file

  assert(xyzFile.is_open() && xyzFile.mode() == OUT);
  xyzFile.setf(std::ios::fixed, std::ios::floatfield);
  xyzFile.setf(std::ios::right, std::ios::adjustfield);
  size_t const nAtoms = geometry.numAtoms();

  xyzFile << nAtoms << std::endl;
  xyzFile << geometry.name() << std::endl;

  std::streamsize precision = xyzFile.precision();
  xyzFile.precision(WRITE_PRECISION);

  for(size_t i = 0; i < nAtoms; i++)
  {
    Atom atom(geometry.atom(i));
    xyzFile << std::setw(MAX_SYMBOL_LENGTH) << atom.symbol << "  ";
    for(size_t j = 0; j < 3 ; j++)
      xyzFile << std::setw(WRITE_PRECISION + 6) << atom.position[j] << "  ";
    xyzFile << std::endl;
  } 
  xyzFile.precision(precision); // Return to original value
  return xyzFile;
}

XYZFile &operator>>(XYZFile &xyzFile, Geometry &geometry)
{  
  // Read geometry from file

  assert(xyzFile.is_open() && xyzFile.mode() == IN);
  geometry.clear();

  size_t r_nAtom;             // Number of atoms read from input file
  char r_label[LINE_LENGTH];  // Geometry label read from input file
  std::string r_symbol;       // Atomic symbol read from input file
  Vector3D r_position;        // Position read from input file

  xyzFile >> r_nAtom; // Read number of atoms
  if(!xyzFile) 
  {
    std::cerr << "Error reading file " << xyzFile.fileName() << std::endl
              << "Invalid number of atoms " << std::endl;
    return xyzFile;
  }
  xyzFile.ignore(LINE_LENGTH, '\n');      // Skip to next line
  xyzFile.getline(r_label, LINE_LENGTH);  // Read geometry label
  geometry.setName(std::string(r_label));
  while(!xyzFile.eof() && geometry.numAtoms() < r_nAtom) 
  {
    xyzFile >> r_symbol >> r_position.x >> r_position.y >> r_position.z;
    xyzFile.ignore(LINE_LENGTH, '\n');         // Skip to next line
    if(!xyzFile) 
    {
      std::cerr << "Error reading file " << xyzFile.fileName() << std::endl
                << "Atoms read: " << geometry.numAtoms() << std::endl;
      geometry.clear();
      return xyzFile;
    }
    geometry.addAtom(Atom(r_symbol, r_position));
  }
  if(geometry.numAtoms() < r_nAtom) 
  {
    std::cerr << "Error reading file " << xyzFile.fileName() << std::endl
              << "Expected " << r_nAtom << " atoms but found only " 
              << geometry.numAtoms() << std::endl;
    geometry.clear();
    return xyzFile;
  }
  geometry.calculateBonds();
  return xyzFile;
}
