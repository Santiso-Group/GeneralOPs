/*
** Copyright 2007, 2008, 2009 Erik Santiso.
** This file is part of mymol.
** mymol is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
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
** .mol2 (Sybyl) file format
**
** Notes: 
**
** The I/O routines assume that atom id numbers are consecutive and start from 1,
** and that substructures are also consecutive. They also assume that the records are 
** in this order: MOLECULE, ATOM, BOND, CRYSIN.
**
** No other fields are read/written, and substructure and charge information are ignored
** when reading to Geometry objects.
**
** CRYSIN records (if present) are only used to define the unit cell for System objects,
** the space group and setting information are discarded on reading. When writing, they're 
** both set to 1 (triclinic/P1).
**
** When reading to System<Geometry> objects, any bonds between different geometry objects
** are discarded. If using this class to represent a protein, the peptide bonds between
** residues will not be stored. If you need to keep them, read to a Geometry object
** (the residue names and IDs will be preserved anyway).
**
** Also, when reading to System<Geometry> objects, the I/O routines rely on residue ID
** information to define the different geometry objects. If no residue ID information
** is present, the whole system will contain a single geometry.
*/

// TODO for a later version: Add the automatic determination of clusters in the absence
// of residue information. This is also necessary for XYZFile.

// Note that, if no bonds are present, this may generate an invalid mol2 file (in the
// sense that it cannot be read by, e.g. VMD, Arguslab, etc...

#ifndef H_MOL2FILE
#define H_MOL2FILE

#include <sstream>
#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"

class MOL2File: public IOFile
{
public:

// Constructors

  MOL2File();                                                // Defines an empty MOL2File
  MOL2File(std::string const &fileName, IOMode const &mode); // Defines a MOL2File with a given file name and I/O mode

// Read from and write to Geometry objects

  friend MOL2File &operator<<(MOL2File &mol2File, Geometry const &geometry);  // Write geometry to file
  friend MOL2File &operator>>(MOL2File &mol2File, Geometry &geometry);        // Read geometry from file

// Read from and write to System<Geometry> objects

  friend MOL2File &operator<<(MOL2File &mol2File, System<Geometry> const &system);  // Write system to file
  friend MOL2File &operator>>(MOL2File &mol2File, System<Geometry> &system);        // Read system from file

// Read from and write to System<Molecule> objects
// (inheritance does not seem to work with the template)

  friend MOL2File &operator<<(MOL2File &mol2File, System<Molecule> const &system);  // Write system to file
  friend MOL2File &operator>>(MOL2File &mol2File, System<Molecule> &system);        // Read system from file

};

#endif

/*
** End of class MOL2File
*/

