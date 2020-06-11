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
** .xtp file format - used to read/write crystal distribution parameters to file
*/

#ifndef H_XTPFILE
#define H_XTPFILE

#include "common/include/iofile.h"
#include "statparameters.h"

class XTPFile: public IOFile
{
public:

// Constructors

  XTPFile();                                                // Defines an empty XTPFile
  XTPFile(std::string const &fileName, IOMode const &mode); // Defines a XTPFile with a given file name and I/O mode

// Read from and write to CrystalDistributionParameters objects

  friend XTPFile &operator<<(XTPFile &xtpFile, CrystalDistributionParameters const &parameters);
    // Write parameters to file
  friend XTPFile &operator>>(XTPFile &xtpFile, CrystalDistributionParameters &parameters);
    // Read parameters from file
};

/*
** End of class XTPFile
*/

#endif
