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
** .mmp file format - used to read molecule map information
**
** This is only an input file class - writing is not implemented.
** Note that the reading method does *not* clear the std::vector
** of MoleculeMaps before adding information to it! This is so that
** it is possible to read molecule maps from multiple files.
*/

#ifndef H_MMPFILE
#define H_MMPFILE

#include <vector>
#include "common/include/iofile.h"
#include "crystdist/include/moleculemap.h"

class MMPFile: public IOFile
{
public:

// Constructors

  MMPFile();                            // Defines an empty MMPFile
  MMPFile(std::string const &fileName); // Defines a MMPFile with a given file name

// Interface

  void setFile(std::string const &fileName);  // Sets the file name

// Read to molecule map objects

  friend MMPFile &operator>>(MMPFile &mmpFile, std::vector<MoleculeMap> &moleculeMaps);  // Read molecule maps from file
};

/*
** End of class MMPFile
*/

// Inlines

inline void MMPFile::setFile(std::string const &fileName)
{
  IOFile::setFile(fileName, IN);
}

#endif
