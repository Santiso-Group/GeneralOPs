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
** CONFIG (DL-POLY) file format
**
** Notes: 
**
** This class reads and writes only positions, not velocities or forces.
** If velocities/forces are present in an input file, the file will be
** read but those values will be discarded.
**
** The class also converts any cubic or orthorhombic periodic boundary
** conditions to parallelepiped upon reading. For systems with octahedral 
** (imcon = 4), rhombic dodecahedral (imcon = 5), and hexagonal prism
** (imcon = 7), the unit cell is read as parallelepiped, so upon writing
** in those cases the file may need to changed manually. If a system has
** only periodic boundary conditions in one direction, it is written as
** if it did not have any (imcon = 0). 
**
** Note that Geometry objects do not store periodic boundary information,
** so upon reading it is discarded, and it is written as not having any
** (imcon = 0).
**
** Note that reading from a CONFIG file does not clear the geometry/system
** object beforehand (in order to allow for the reading of connectivity
** from FIELD files). If that is what is needed, you need to explicitly
** call the clear() method for the geometry/system object before
** reading. 
**
** If connectivity is not defined previously, bonds are found automatically. 
** This can be slow for large systems.
**
** When using custom-atom types (e.g. defining dummy atoms, united atoms,
** or using a mesoscale model), the atom should be added to the atom.dat
** file to avoid errors. Masses and charges are not assigned from the 
** values in the CONFIG file.
**
** As is the case with MOL2File and PDBFile, when reading or writing to
** System<Geometry> or System<Molecule> objects, only intra-molecular
** bonds are read/written, not bonds between different residues. If
** you want to preserve all bonds, you need to read to a Geometry object.
**
** Also, when reading to a System<Geometry> or System<Molecule> object, if
** connectivity has not been defined previously, the entire system will
** contain a single geometry with all the atoms.
*/

#ifndef H_CONFIGFILE
#define H_CONFIGFILE

#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"

class CONFIGFile: public IOFile
{
public:

// Constructors

  CONFIGFile();                                                // Defines an empty CONFIGFile
  CONFIGFile(std::string const &fileName, IOMode const &mode); // Defines a CONFIG File with a given file name and I/O mode

// Read from and write to Geometry objects

  friend CONFIGFile &operator<<(CONFIGFile &configFile, Geometry const &geometry); // Write geometry to file
  friend CONFIGFile &operator>>(CONFIGFile &configFile, Geometry &geometry);       // Read geometry from file

// Read from and write to System<Geometry> objects

  friend CONFIGFile &operator<<(CONFIGFile &configFile, System<Geometry> const &system); // Write system to file
  friend CONFIGFile &operator>>(CONFIGFile &configFile, System<Geometry> &system);       // Read system from file

// Read from and write to System<Molecule> objects (inheritance not working with template)

  friend CONFIGFile &operator<<(CONFIGFile &configFile, System<Molecule> const &system); // Write system to file
  friend CONFIGFile &operator>>(CONFIGFile &configFile, System<Molecule> &system);       // Read system from file
};

#endif

/*
** End of class CONFIGFile
*/

