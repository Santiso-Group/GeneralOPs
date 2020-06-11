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
** .pkt file format - used to read/write peak group tables
*/

#ifndef H_PKTFILE
#define H_PKTFILE

#include "common/include/iofile.h"
#include "crystdist/include/peaktable.h"

class PKTFile: public IOFile
{
public:

// Constructors

  PKTFile();                                                // Defines an empty PKTFile
  PKTFile(std::string const &fileName, IOMode const &mode); // Defines a PKTFile with a given file name and I/O mode

// Read from and write to PeakTable objects

  friend PKTFile &operator<<(PKTFile &pktFile, PeakTable const &table); // Write table to file
  friend PKTFile &operator>>(PKTFile &pktFile, PeakTable &table);       // Read table from file
};

/*
** End of class PKTFile
*/

#endif
