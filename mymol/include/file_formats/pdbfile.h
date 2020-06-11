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
** .pdb file format
**
** Notes: 
**
** This class writes/reads only COMPND, ATOM/HETATM, CONECT, and 
** CRYST1 records, any other information is discarded. The class does
** not yet handle alternate coordinate systems. COMPND records (if
** any) must be first in the file, followed by CRYST1 records (if any)
** followed by ATOM/HETATM records, followed by CONECT records 
** (if any). TER records are ignored when reading, and not written.
**
** Note that reading from a pdb file does not clear the geometry/system
** object beforehand (in order to allow for the reading of connectivity
** from psf files). If that is what is needed, you need to explicitly
** call the clear() method for the geometry/system object before
** reading. Note that, if a geometry is not empty when a pdb file is
** read, CONECT information is ignored (to preserve the connectivity
** defined in the psf file).
**
** If connectivity is not defined previously, and no CONECT records
** are found in the pdb file, bonds are found automatically. This can
** be slow for large systems.
**
** Right now, this can only read properly elements with a one-letter
** symbol. For a later version, a way to parse the atom types to
** determine the chemical element should be added. Also, it would be
** good to read atom name and symbol records separately, and find
** the symbol from the name if not present.
**
** As is the case with MOL2File, when reading or writing to
** System<Geometry> or System<Molecule> objects, only intra-molecular
** bonds are read/written, not bonds between different residues. If
** you want to preserve all bonds, you need to read to a Geometry
** object.
*/

#ifndef H_PDBFILE
#define H_PDBFILE

#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"

class PDBFile: public IOFile
{
public:

// Constructors

  PDBFile();                                                // Defines an empty PDBFile
  PDBFile(std::string const &fileName, IOMode const &mode); // Defines a PDB File with a given file name and I/O mode

// Read from and write to Geometry objects

  friend PDBFile &operator<<(PDBFile &pdbFile, Geometry const &geometry); // Write geometry to file
  friend PDBFile &operator>>(PDBFile &pdbFile, Geometry &geometry);       // Read geometry from file

// Read from and write to System<Geometry> objects

  friend PDBFile &operator<<(PDBFile &pdbFile, System<Geometry> const &system); // Write system to file
  friend PDBFile &operator>>(PDBFile &pdbFile, System<Geometry> &system);       // Read system from file

// Read from and write to System<Molecule> objects (inheritance not working with template)

  friend PDBFile &operator<<(PDBFile &pdbFile, System<Molecule> const &system); // Write system to file
  friend PDBFile &operator>>(PDBFile &pdbFile, System<Molecule> &system);       // Read system from file
};

#endif

/*
** End of class PDBFile
*/

