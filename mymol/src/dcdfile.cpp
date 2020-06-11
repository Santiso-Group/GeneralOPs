/*
** Copyright 2007-2011 Erik Santiso.
** This file is part of mymol.
** mymol is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
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
** .dcd file format
*/

#include <cmath>
#include "mymath/include/mpconstants.h"
#include "mymol/include/file_formats/fileformats.h"

// Constructors

DCDFile::DCDFile()
:
IOFile(), numAtoms_(0), numFramesRead_(0), numFrames_(0)
{}

DCDFile::DCDFile(std::string const &fileName)
:
IOFile(fileName, IN, BINARY), numAtoms_(0), numFramesRead_(0), numFrames_(0)
{ initializeDCDFile(); }

// Interface

void DCDFile::skipFrame()
{
  assert(is_open());
  //if(numFramesRead_ >= numFrames_)
  //{
  //  std::cerr << "Error in DCDFile: Attempted to skip frames beyond the end of the file" << std::endl;
  //  return;
  //}
  size_t dummyInt;
  Real dummyReal;
  SmallReal dummySmall;
  // Length of the unit cell field
  read(reinterpret_cast<char *>(&dummyInt), sizeof(size_t));
  // Unit cell information
  read(reinterpret_cast<char *>(&dummyReal), sizeof(Real));
  read(reinterpret_cast<char *>(&dummyReal), sizeof(Real));
  read(reinterpret_cast<char *>(&dummyReal), sizeof(Real));
  read(reinterpret_cast<char *>(&dummyReal), sizeof(Real));
  read(reinterpret_cast<char *>(&dummyReal), sizeof(Real));
  read(reinterpret_cast<char *>(&dummyReal), sizeof(Real));
  // Length of the unit cell field (again)
  read(reinterpret_cast<char *>(&dummyInt), sizeof(size_t));
  // Four times the number of atoms
  for(size_t i = 0; i < 4; ++i) get();
  // x-coordinates
  for(size_t i = 0; i < numAtoms_; ++i)
    read(reinterpret_cast<char *>(&dummySmall), sizeof(SmallReal));
  // Four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) get();
  // y-coordinates
  for(size_t i = 0; i < numAtoms_; ++i)
    read(reinterpret_cast<char *>(&dummySmall), sizeof(SmallReal));
  // Four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) get();
  // z-coordinates
  for(size_t i = 0; i < numAtoms_; ++i)
    read(reinterpret_cast<char *>(&dummySmall), sizeof(SmallReal));
  // Four times the number of atoms
  for(size_t i = 0; i < 4; ++i) get();
  // Done
  ++numFramesRead_;
  return;
}

// Operators

DCDFile &operator>>(DCDFile &dcdFile, Geometry &geometry)
{
  assert(dcdFile.is_open());
  //if(dcdFile.numFramesRead_ >= dcdFile.numFrames_)
  //{
  //  std::cerr << "Error in DCDFile: Attempted to read beyond the end of the file" << std::endl;
  //  return dcdFile;
  //}
  if(geometry.numAtoms() == 0)
  {
    std::cerr << "Error in DCDFile: Attempted to read from a dcd file before getting " << std::endl
              << "atom information from another file." << std::endl;
    return dcdFile;
  }
  if(dcdFile.numAtoms_ != geometry.numAtoms())
  {
    std::cerr << "Error in DCDFile: Inconsistent number of atoms in geometry and dcd file." << std::endl;
    return dcdFile;
  }
  // Read the length of the unit cell field
  size_t length;
  dcdFile.read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 48) // Length of unit cell field should be 48
  {
    std::cerr << "Error reading .dcd file " << dcdFile.fileName() << ": Invalid file format." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Read unit cell information
  Real a, b, c, cosAlpha, cosBeta, cosGamma; // Crystallographic constants read from file
  dcdFile.read(reinterpret_cast<char *>(&a), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosGamma), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&b), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosBeta), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosAlpha), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&c), sizeof(Real));
  if(a < 0.0 || b < 0.0 || c < 0.0 || fabs(cosAlpha) > 1.0 || 
     fabs(cosBeta)  > 1.0 || fabs(cosGamma) > 1.0 || !dcdFile)
  {
    std::cerr << "Error reading dcd file " << dcdFile.fileName() 
              << ": Invalid unit cell data." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Read the length of the unit cell field (again)
  dcdFile.read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 48) // Length of unit cell field should be 48
  {
    std::cerr << "Error reading .dcd file " << dcdFile.fileName() << ": Invalid file format." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Skip four times the number of atoms
  for(size_t i = 0; i < 4; ++i) dcdFile.get();
  // Read x-coordinates
  for(size_t i = 0; i < dcdFile.numAtoms_; ++i)
  {
    SmallReal newx;
    dcdFile.read(reinterpret_cast<char *>(&newx), sizeof(SmallReal));
    geometry.setx(i, (Real)newx);
  }
  // Skip four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) dcdFile.get();
  // Read y-coordinates
  for(size_t i = 0; i < dcdFile.numAtoms_; ++i)
  {
    SmallReal newy;
    dcdFile.read(reinterpret_cast<char *>(&newy), sizeof(SmallReal));
    geometry.sety(i, (Real)newy);
  }
  // Skip four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) dcdFile.get();
  // Read z-coordinates
  for(size_t i = 0; i < dcdFile.numAtoms_; ++i)
  {
    SmallReal newz;
    dcdFile.read(reinterpret_cast<char *>(&newz), sizeof(SmallReal));
    geometry.setz(i, (Real)newz);
  }
  // Skip four times the number of atoms
  for(size_t i = 0; i < 4; ++i) dcdFile.get();
  // Done
  ++dcdFile.numFramesRead_;
  return dcdFile;
}

DCDFile &operator>>(DCDFile &dcdFile, System<Geometry> &system)
{
  assert(dcdFile.is_open());
  //if(dcdFile.numFramesRead_ >= dcdFile.numFrames_)
  //{
  //  std::cerr << "Error in DCDFile: Attempted to read beyond the end of the file" << std::endl;
  //  return dcdFile;
  //}
  if(system.size() == 0)
  {
    std::cerr << "Error in DCDFile: Attempted to read from a dcd file before getting " << std::endl
              << "atom information from another file." << std::endl;
    return dcdFile;
  }
  size_t nAtoms = 0;
  for(size_t i = 0; i < system.size(); ++i)
    nAtoms += system[i].numAtoms();
  if(dcdFile.numAtoms_ != nAtoms)
  {
    std::cerr << "Error in DCDFile: Inconsistent number of atoms in system and dcd file." << std::endl;
    return dcdFile;
  }
  // Read the length of the unit cell field
  size_t length;
  dcdFile.read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 48) // Length of unit cell field should be 48
  {
    std::cerr << "Error reading .dcd file " << dcdFile.fileName() << ": Invalid file format." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Read unit cell information
  Real a, b, c, cosAlpha, cosBeta, cosGamma; // Crystallographic constants read from file
  dcdFile.read(reinterpret_cast<char *>(&a), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosGamma), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&b), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosBeta), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosAlpha), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&c), sizeof(Real));
  if(a < 0.0 || b < 0.0 || c < 0.0 || fabs(cosAlpha) > 1.0 || 
     fabs(cosBeta)  > 1.0 || fabs(cosGamma) > 1.0 || !dcdFile)
  {
    std::cerr << "Error reading dcd file " << dcdFile.fileName() 
              << ": Invalid unit cell data." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Store unit cell information
  system.setLattice(Lattice(a, b, c,
                            RAD_TO_GRAD*acos(cosAlpha), 
                            RAD_TO_GRAD*acos(cosBeta), 
                            RAD_TO_GRAD*acos(cosGamma),
                            system.lattice().type()));  // Keep the lattice type constant - modified 01/17/09
  // Read the length of the unit cell field (again)
  dcdFile.read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 48) // Length of unit cell field should be 48
  {
    std::cerr << "Error reading .dcd file " << dcdFile.fileName() << ": Invalid file format." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Skip four times the number of atoms
  for(size_t i = 0; i < 4; ++i) dcdFile.get();
  // Read x-coordinates
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j)
  {
    SmallReal newx;
    dcdFile.read(reinterpret_cast<char *>(&newx), sizeof(SmallReal));
    system[i].setx(j, (Real)newx);
  }
  // Skip four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) dcdFile.get();
  // Read y-coordinates
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j)
  {
    SmallReal newy;
    dcdFile.read(reinterpret_cast<char *>(&newy), sizeof(SmallReal));
    system[i].sety(j, (Real)newy);
  }
  // Skip four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) dcdFile.get();
  // Read z-coordinates
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j)
  {
    SmallReal newz;
    dcdFile.read(reinterpret_cast<char *>(&newz), sizeof(SmallReal));
    system[i].setz(j, (Real)newz);
  }
  // Skip four times the number of atoms
  for(size_t i = 0; i < 4; ++i) dcdFile.get();
  // Done
  ++dcdFile.numFramesRead_;
  return dcdFile;
}

DCDFile &operator>>(DCDFile &dcdFile, System<Molecule> &system)
{
  assert(dcdFile.is_open());
  //if(dcdFile.numFramesRead_ >= dcdFile.numFrames_)
  //{
  //  std::cerr << "Error in DCDFile: Attempted to read beyond the end of the file" << std::endl;
  //  return dcdFile;
  //}
  if(system.size() == 0)
  {
    std::cerr << "Error in DCDFile: Attempted to read from a dcd file before getting " << std::endl
              << "atom information from another file." << std::endl;
    return dcdFile;
  }
  size_t nAtoms = 0;
  for(size_t i = 0; i < system.size(); ++i)
    nAtoms += system[i].numAtoms();
  if(dcdFile.numAtoms_ != nAtoms)
  {
    std::cerr << "Error in DCDFile: Inconsistent number of atoms in system and dcd file." << std::endl;
    return dcdFile;
  }
  // Read the length of the unit cell field
  size_t length;
  dcdFile.read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 48) // Length of unit cell field should be 48
  {
    std::cerr << "Error reading .dcd file " << dcdFile.fileName() << ": Invalid file format." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Read unit cell information
  Real a, b, c, cosAlpha, cosBeta, cosGamma; // Crystallographic constants read from file
  dcdFile.read(reinterpret_cast<char *>(&a), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosGamma), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&b), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosBeta), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&cosAlpha), sizeof(Real));
  dcdFile.read(reinterpret_cast<char *>(&c), sizeof(Real));
  if(a < 0.0 || b < 0.0 || c < 0.0 || fabs(cosAlpha) > 1.0 || 
     fabs(cosBeta)  > 1.0 || fabs(cosGamma) > 1.0 || !dcdFile)
  {
    std::cerr << "Error reading dcd file " << dcdFile.fileName() 
              << ": Invalid unit cell data." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Store unit cell information
  system.setLattice(Lattice(a, b, c,
                            RAD_TO_GRAD*acos(cosAlpha), 
                            RAD_TO_GRAD*acos(cosBeta), 
                            RAD_TO_GRAD*acos(cosGamma),
                            system.lattice().type()));  // Keep the lattice type constant - modified 01/17/09
  // Read the length of the unit cell field (again)
  dcdFile.read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 48) // Length of unit cell field should be 48
  {
    std::cerr << "Error reading .dcd file " << dcdFile.fileName() << ": Invalid file format." << std::endl;
    dcdFile.clear();
    return dcdFile;
  }
  // Skip four times the number of atoms
  for(size_t i = 0; i < 4; ++i) dcdFile.get();
  // Read x-coordinates
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j)
  {
    SmallReal newx;
    dcdFile.read(reinterpret_cast<char *>(&newx), sizeof(SmallReal));
    system[i].setx(j, (Real)newx);
  }
  // Skip four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) dcdFile.get();
  // Read y-coordinates
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j)
  {
    SmallReal newy;
    dcdFile.read(reinterpret_cast<char *>(&newy), sizeof(SmallReal));
    system[i].sety(j, (Real)newy);
  }
  // Skip four times the number of atoms (twice)
  for(size_t i = 0; i < 8; ++i) dcdFile.get();
  // Read z-coordinates
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j)
  {
    SmallReal newz;
    dcdFile.read(reinterpret_cast<char *>(&newz), sizeof(SmallReal));
    system[i].setz(j, (Real)newz);
  }
  // Skip four times the number of atoms
  for(size_t i = 0; i < 4; ++i) dcdFile.get();
  // Done
  ++dcdFile.numFramesRead_;
  return dcdFile;
}

// Private functions

void DCDFile::initializeDCDFile()
{
  assert(is_open());
  // Read length of dcd header
  size_t length;
  read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 84) // Length of dcd header should be 84
  {
    std::cerr << "Error reading .dcd file " << fileName() << ": Invalid file format." << std::endl;
    clear();
    return;
  }
  // Read 'CORD' field
  char cord[4];
  for(size_t i = 0; i < 4; ++i)
    get(cord[i]);
  if(cord[0] != 'C' || cord[1] != 'O' || cord[2] != 'R' || cord[3] != 'D') // Check that CORD field is OK
  {
    std::cerr << "Error reading .dcd file " << fileName() << ": Invalid file format." << std::endl;
    clear();
    return;
  }
  // Read the number of frames
  read(reinterpret_cast<char *>(&numFrames_), sizeof(size_t));
  // Skip the rest of the header
  for(size_t i = 0; i < 76; ++i) get();
  // Read another 84
  read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 84) // Length of dcd header should be 84
  {
    std::cerr << "Error reading .dcd file " << fileName() << ": Invalid file format." << std::endl;
    clear();
    return;
  }
  // Read the length of the title records
  read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length % 80 != 4) // Title record are 80-bytes long + 4 bytes to store number of records
  {
    std::cerr << "Error reading .dcd file " << fileName() << ": Invalid file format." << std::endl;
    clear();
    return;
  }
  // Skip the number of title records, and the title records
  for(size_t i = 0; i < length; ++i) get();
  // Read the length of the title records again
  read(reinterpret_cast<char *>(&length), sizeof(size_t));
  // Read the length of the number of atoms field (4)
  read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 4) // Should be 4
  {
    std::cerr << "Error reading .dcd file " << fileName() << ": Invalid file format." << std::endl;
    clear();
    return;
  }
  // Read the number of atoms
  read(reinterpret_cast<char *>(&numAtoms_), sizeof(size_t));
  // Read another 4
  read(reinterpret_cast<char *>(&length), sizeof(size_t));
  if(length != 4) // Should be 4
  {
    std::cerr << "Error reading .dcd file " << fileName() << ": Invalid file format." << std::endl;
    clear();
    return;
  }
  // Done
  if(!(*this))
  {
    std::cerr << "Error initializing .dcd file " << fileName() << std::endl;
    clear();
    return;
  }
}
