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
** HISTORY (DL_POLY) file format
**
** Notes: 
**
** HISTORY files do not contain any bonding or topology information, if
** that is needed a FIELD file (or any other file type that contains
** bonding information) must be read before.
**
** Since a HISTORY file contains, in general, many frames, it is possible
** to read geometries from the file more than once. The member function
** numFramesRead() returns the number of frames that have already been
** read. 
**
** Currently it is only possible to read from HISTORY files. Writing is
** not implemented.
**
** As with CONFIG files, this class only reads positions, any velocity
** or force information is discarded. Charge and mass information is also
** discarded (as with CONFIG files). Finally, the supercell is read in
** the same way as with CONFIG files - see configfile.h for details.
*/

#ifndef H_HISTORYFILE
#define H_HISTORYFILE

#include <string>
#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"

class HISTORYFile: public IOFile
{
public:

// Constructors

  HISTORYFile();                            // Defines an empty HISTORYFile
  HISTORYFile(std::string const &fileName); // Defines a HISTORYFile with a given file name

// Accessors

  size_t const numFramesRead() const; // Returns the number of frames that have been already read

// Interface

   void setFile(std::string const &fileName); // Sets the file name
   void skipFrame();                          // Skips one frame
   void skipFrames(size_t const numFrames);   // Skip numFrames frames
   void clear();                              // Closes the file and clears file data

// Read to Geometry objects

  friend HISTORYFile &operator>>(HISTORYFile &historyFile, Geometry &geometry); // Read geometry from file

// Read to System<Geometry> objects

  friend HISTORYFile &operator>>(HISTORYFile &historyFile, System<Geometry> &system); // Read system from file

// Read to System<Molecule> objects (inheritance not working with template)

  friend HISTORYFile &operator>>(HISTORYFile &historyFile, System<Molecule> &system); // Read system from file

private:

  size_t numAtoms_;      // Number of atoms read from the file header (natms in DL_POLY manual)
  size_t numFramesRead_; // Number of frames that have been read from the file
  size_t trajectoryKey_; // Trajectory key (keytrj in DL_POLY manual)
  size_t supercellKey_;  // Periodic boundary key (imcon in DL_POLY manual) 
  std::string title_;    // Trajectory title

// Private functions

  void initializeHISTORYFile();              // Reads file headers and stores 
                                             // the number of atoms and 
                                             // supercell key
  void setFile(std::string const &fileName,  // Format is read-only - blocking
                IOMode const &mode,          // the general form from IOFile
                IOFormat const &format);
};

/*
** End of class HISTORYFile
*/

// Inlines

inline size_t const HISTORYFile::numFramesRead() const
{
  return numFramesRead_;
}

inline void HISTORYFile::setFile(std::string const &fileName)
{
  IOFile::setFile(fileName, IN);
  numAtoms_ = numFramesRead_ = trajectoryKey_ = supercellKey_ = 0;
  title_.clear();
  initializeHISTORYFile();
}

inline void HISTORYFile::skipFrames(size_t const numFrames)
{
  for(size_t i = 0; i < numFrames; ++i)
    skipFrame();
}

inline void HISTORYFile::clear()
{
  IOFile::clear();
  numAtoms_ = numFramesRead_ = trajectoryKey_ = supercellKey_ = 0;
  title_.clear();
}

inline void HISTORYFile::setFile(std::string const &fileName,
                                 IOMode const &mode,
                                 IOFormat const &format)
{
  std::cerr << "Error in HISTORYFile: Writing not implemented!" << std::endl;
}

#endif

