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
** .dcd (as used by NAMD/VMD) file format
**
** Notes: 
**
** Before trying to read a geometry from a .dcd file, the geometry
** must be initialized to contain the types of all atoms (since the dcd
** format does not include this information). This can be done by reading
** a file containing this information (e.g. .psf, .xyz) before reading
** the .dcd file (.psf is the most efficient). For example:
**
** XYZFile myXYZFile("file.xyz", IN)
** DCDFile myDCDFile("file.dcd", IN)
** Geometry myGeometry;
** myXYZFile >> myGeometry;
** myDCDFile >> myGeometry;
**
** Since a .dcd file contains, in general, many frames, it is possible to
** read geometries from the file more than once. The member function
** numFrames() and numFramesRead() return the total number of frames in
** the file and the number of frames that have been already read.
**
** Currently it is only possible to read from dcd files - writing is not
** implemented yet.
**
** This class only reads the number of frames and number of atoms from the
** dcd header, other variables are discarded. It also does not recognize
** files with different endianness.
*/

#ifndef H_DCDFILE
#define H_DCDFILE

#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"

class DCDFile: public IOFile
{

//using IOFile::setFile;

public:

// Constructors

  DCDFile();                            // Defines an empty DCDFile
  DCDFile(std::string const &fileName); // Defines a DCDFile with a given file name

// Accessors

  size_t const numFramesRead() const; // Returns the number of frames that have been already read
  size_t const numFrames() const;     // Returns the number of frames stored in the dcd file

// Interface

   void setFile(std::string const &fileName); // Sets the file name
   void skipFrame();                          // Skips one frame
   void skipFrames(size_t const numFrames);   // Skip numFrames frames
   void skipToLastFrame();                    // Skips to the last frame
   void clear();                              // Closes the file and clears file data

// Read to Geometry objects

  friend DCDFile &operator>>(DCDFile &dcdFile, Geometry &geometry); // Read geometry from file

// Read to System<Geometry> objects

  friend DCDFile &operator>>(DCDFile &dcdFile, System<Geometry> &system); // Read system from file

// Read to System<Molecule> objects (inheritance not working with template)

  friend DCDFile &operator>>(DCDFile &dcdFile, System<Molecule> &system); // Read system from file

private:

  size_t numAtoms_;       // Number of atoms read from the dcd file header
  size_t numFramesRead_;  // Number of frames that have been read from the dcd file
  size_t numFrames_;      // Number of frames in the dcd file

// Private functions

  void initializeDCDFile();                 // Reads DCD headers and stores 
                                            // the number of frames
  void setFile(std::string const &fileName, // Format is read-only - blocking
               IOMode const &mode,          // the general form from IOFile
               IOFormat const &format);
};

/*
** End of class DCDFile
*/

// Inlines

inline size_t const DCDFile::numFrames() const
{
  return numFrames_;
}

inline size_t const DCDFile::numFramesRead() const
{
  return numFramesRead_;
}

inline void DCDFile::setFile(std::string const &fileName)
{
  IOFile::setFile(fileName, IN, BINARY);
  numAtoms_ = numFramesRead_ = numFrames_ = 0;
  initializeDCDFile();
}

inline void DCDFile::skipFrames(size_t const numFrames)
{
  for(size_t i = 0; i < numFrames; ++i)
    skipFrame();
}

inline void DCDFile::skipToLastFrame()
{
  while(numFramesRead_ < numFrames_ - 1)
    skipFrame();
}

inline void DCDFile::clear()
{
  IOFile::clear();
  numAtoms_ = numFramesRead_ = numFrames_ = 0;
}

inline void DCDFile::setFile(std::string const &fileName,
                             IOMode const &mode,
                             IOFormat const &format)
{
  std::cerr << "Error in DCDFile: Writing not implemented!" << std::endl;
}

#endif

