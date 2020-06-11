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
** .ptm file format - used to read/write point molecule information
**
** Note that a single .ptm file may contain many frames, thus the file
** is not automatically closed after reading/writing.
*/

#ifndef H_PTMFILE
#define H_PTMFILE

#include "common/include/iofile.h"
#include "mymol/include/system.h"
#include "crystdist/include/pointmolecule.h"

class PTMFile: public IOFile
{
public:

// Constructors

  PTMFile();                                                // Defines an empty PTMFile
  PTMFile(std::string const &fileName, IOMode const &mode); // Defines a PTMFile with a given file name and I/O mode

// Read from and write to System<PointMolecule> objects

  friend PTMFile &operator<<(PTMFile &ptmFile, System<PointMolecule> const &system); // Write system to file
  friend PTMFile &operator>>(PTMFile &ptmFile, System<PointMolecule> &system);       // Read system from file
};

/*
** End of class PTMFile
*/

#endif
