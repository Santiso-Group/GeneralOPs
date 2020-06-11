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
** .xyz file format
**
** Notes:
**
** This routine clears the original geometry object, and finds bonds automatically.
** This can be slow for large systems.
*/

#ifndef H_XYZFILE
#define H_XYZFILE

#include "common/include/iofile.h"
#include "mymol/include/geometry.h"

class XYZFile: public IOFile
{
public:

// Constructors

  XYZFile();                                                // Defines an empty XYZFile
  XYZFile(std::string const &fileName, IOMode const &mode); // Defines a XYZFile with a given file name and I/O mode

// Read from and write to Geometry objects

  friend XYZFile &operator<<(XYZFile &xyzFile, Geometry const &geometry); // Write geometry to file
  friend XYZFile &operator>>(XYZFile &xyzFile, Geometry &geometry);       // Read geometry from file
};

#endif

/*
** End of class XYZFile
*/

