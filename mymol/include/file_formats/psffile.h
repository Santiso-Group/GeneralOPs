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
** .psf (X-PLOR) file format
**
** Notes: 
**
** Currently it is only possible to read from psf files - writing is not
** implemented yet.
**
** This class only reads atoms and bonds, the remaining information is
** discarded. Positions are set to zero by default. Note that reading
** some file formats (e.g. xyz or mol2) after reading a psf file 
** will overwrite the information read from the psf file. Exceptions
** include pdb, dcd, and CONFIG files. Also, only one segment can be read.
**
** Right now, this can only read properly elements with a one-letter
** symbol. For a later version, a way to parse the atom types to
** determine the chemical element should be added.
**
** As is the case with MOL2File, when reading to a System<Geometry>
** only intra-molecular bonds are stored, not bonds between different
** residues. If you need to keep those bonds, you should read to a
** Geometry object.
*/

#ifndef H_PSFFILE
#define H_PSFFILE

#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"

class PSFFile: public IOFile
{
public:

// Constructors

  PSFFile();                            // Defines an empty PSFFile
  PSFFile(std::string const &fileName); // Defines a PSF File with a given file name

// Interface

  void setFile(std::string const &fileName); // Sets the file name

// Read to Geometry objects

  friend PSFFile &operator>>(PSFFile &psfFile, Geometry &geometry); // Read geometry from file

// Read to System<Geometry> objects

  friend PSFFile &operator>>(PSFFile &psfFile, System<Geometry> &system); // Read system from file

// Read to System<Molecule> objects (inheritance not working with template)

  friend PSFFile &operator>>(PSFFile &psfFile, System<Molecule> &system); // Read system from file

private:

  void setFile(std::string const &fileName, // Format is read-only - blocking
               IOMode const &mode,          // the general form from IOFile
               IOFormat const &format);
};

#endif

/*
** End of class PSFFile
*/

// Inlines

inline void PSFFile::setFile(std::string const &fileName)
{
  IOFile::setFile(fileName, IN);
}

inline void PSFFile::setFile(std::string const &fileName,
                             IOMode const &mode,
                             IOFormat const &format)
{
  std::cerr << "Error in PSFFile: Writing not implemented!" << std::endl;
}

