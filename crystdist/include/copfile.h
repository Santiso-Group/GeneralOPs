/*
** Copyright 2008-2011 Erik Santiso.
** This file is part of crystdist.
** crystdist is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** crystdist is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with crystdist. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** .cop file format - used to read/write crystal order parameters to file
*/

#ifndef H_COPFILE
#define H_COPFILE

#include "common/include/iofile.h"
#include "crystdist/include/crystalops.h"

class COPFile: public IOFile
{
public:

// Constructors

  COPFile();                                                // Defines an empty COPFile
  COPFile(std::string const &fileName, IOMode const &mode); // Defines a COPFile with a given file name and I/O mode

// Read from and write to CrystalOrderParameter objects

  friend COPFile &operator<<(COPFile &copFile, CrystalOrderParameters const &parameters); // Write ops to file
  friend COPFile &operator>>(COPFile &copFile, CrystalOrderParameters &parameters);       // Read ops from file
};

/*
** End of class COPFile
*/

#endif
