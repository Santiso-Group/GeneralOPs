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
** HISTORY file format
*/

#include "mymol/include/file_formats/fileformats.h"

// Constructors

HISTORYFile::HISTORYFile()
:
IOFile(), numAtoms_(0), numFramesRead_(0), trajectoryKey_(0), 
supercellKey_(0), title_()
{}

HISTORYFile::HISTORYFile(std::string const &fileName)
:
IOFile(fileName, IN), numAtoms_(0), numFramesRead_(0), trajectoryKey_(0),
supercellKey_(0), title_()
{ initializeHISTORYFile(); }

// Interface

void HISTORYFile::skipFrame()
{
  assert(is_open());

  // Skip header record
  ignore(LINE_LENGTH, '\n');
  if(eof()) return; // Reached end of the file, return

  // Skip unit cell information if present
  if(supercellKey_ > 0)
    for(size_t i = 0; i < 3; ++i) ignore(LINE_LENGTH, '\n');

  // Skip atom records
  for(size_t i = 0; i < numAtoms_; ++i)
  {
   ignore(LINE_LENGTH, '\n'); // Atom header record
   ignore(LINE_LENGTH, '\n'); // Position record

    // If present, skip velocity and/or force records
    for(size_t j = 0; j < trajectoryKey_; ++j)
      ignore(LINE_LENGTH, '\n');

  } // End of loop to skip atom records

 // Update number of frames read
  ++numFramesRead_;

  return; 
}

// Operators

HISTORYFile &operator>>(HISTORYFile &historyFile, Geometry &geometry)
{
  assert(historyFile.is_open());
  bool newGeometry = (geometry.numAtoms() == 0);

  if(!newGeometry && historyFile.numAtoms_ != geometry.numAtoms())
  {
    std::cerr << "Error in HISTORYFile: Inconsistent number of atoms in geometry and history file." << std::endl;
    return historyFile;
  }

  // Read configuration header record (mostly for consistency check)

  std::string line; // Line read from file

  // These variables named as in the DL_POLY manual
  std::string timestep; // The string "timestep"
  size_t nstep;         // Current time step
  size_t natms;         // Number of atoms
  size_t keytrj;        // Trajectory key
  size_t imcon;         // Periodic boundary key

  std::getline(historyFile, line);

  if(historyFile.eof()) return historyFile; // Reached end of the file, return

  std::stringstream headerStream(line);
  headerStream >> timestep;
  if(timestep.find("timestep") == timestep.npos)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": bad timestep record" << std::endl;
    return historyFile;
  }
  headerStream >> nstep >> natms >> keytrj >> imcon;
  if(headerStream.fail())
  {
    std::cerr << "Error reading timestep record from HISTORY file " 
              << historyFile.fileName() << std::endl;
    return historyFile;
  }
  if(natms != historyFile.numAtoms_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent number of atoms" << std::endl;
    return historyFile;
  }
  if(keytrj != historyFile.trajectoryKey_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent trajectory key" << std::endl;
    return historyFile;
  }
  if(imcon != historyFile.supercellKey_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent periodic boundary key" << std::endl;
    return historyFile;
  }

  // Skip unit cell information if present
  if(historyFile.supercellKey_ > 0)
    for(size_t i = 0; i < 3; ++i) historyFile.ignore(LINE_LENGTH, '\n');

  // Read atom records
  for(size_t i = 0; i < historyFile.numAtoms_; ++i)
  {
    std::string r_symbol; // Atom symbol read from file
    size_t r_index;       // Atom index read from file
    Real r_x, r_y, r_z;   // Coordinates read from file

    // Read symbol and index
    historyFile >> r_symbol >> r_index;

    if(historyFile.fail())
    {
      std::cerr << "Error reading atom record from HISTORY file " << historyFile.fileName() << std::endl
                << "Atoms read: " << i << std::endl;
      return historyFile;
    }
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    if(r_index != i + 1)
    {
      std::cerr << "Error reading HISTORY file " << historyFile.fileName()
                << ": atom indices not correlated" << std::endl;
      return historyFile;
    }

    // Read position
    historyFile >> r_x >> r_y >> r_z;

    if(historyFile.fail())
    {
      std::cerr << "Error reading atom position from file " << historyFile.fileName() << std::endl
                << "Atoms read: " << i << std::endl;
      return historyFile;
    }
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // If present, skip velocity and/or force records
    for(size_t j = 0; j < historyFile.trajectoryKey_; ++j)
      historyFile.ignore(LINE_LENGTH, '\n');

    if(newGeometry)
      geometry.addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z)));
    else
    {
      if(r_symbol.find(geometry.atom(i).symbol) == r_symbol.npos)
      {
        std::cerr << "Inconsistent data in HISTORY file " << historyFile.fileName() << std::endl;
        return historyFile;
      }
      geometry.setPosition(i, Vector3D(r_x, r_y, r_z));
    }
  } // End of loop to read atom records

 // Update number of frames read
  ++historyFile.numFramesRead_;

  return historyFile;
}

HISTORYFile &operator>>(HISTORYFile &historyFile, System<Geometry> &system)
{
  assert(historyFile.is_open());
  bool newSystem = (system.size() == 0);

  if(!newSystem) // Consistency checks
  {
    size_t nAtoms = 0;
    for(size_t i = 0; i < system.size(); ++i)
      nAtoms += system[i].numAtoms();
    if(nAtoms != historyFile.numAtoms_)
    {
      std::cerr << "Error in HISTORYFile: Inconsistent number of atoms in system and history file." << std::endl;
      return historyFile;
    }
  }

  // Read configuration header record (mostly for consistency check)

  std::string line; // Line read from file

  // These variables named as in the DL_POLY manual
  std::string timestep; // The string "timestep"
  size_t nstep;         // Current time step
  size_t natms;         // Number of atoms
  size_t keytrj;        // Trajectory key
  size_t imcon;         // Periodic boundary key

  std::getline(historyFile, line);

  if(historyFile.eof()) return historyFile; // Reached end of the file, return

  std::stringstream headerStream(line);
  headerStream >> timestep;
  if(timestep.find("timestep") == timestep.npos)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": bad timestep record" << std::endl;
    return historyFile;
  }
  headerStream >> nstep >> natms >> keytrj >> imcon;
  if(headerStream.fail())
  {
    std::cerr << "Error reading timestep record from HISTORY file " 
              << historyFile.fileName() << std::endl;
    return historyFile;
  }
  if(natms != historyFile.numAtoms_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent number of atoms" << std::endl;
    return historyFile;
  }
  if(keytrj != historyFile.trajectoryKey_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent trajectory key" << std::endl;
    return historyFile;
  }
  if(imcon != historyFile.supercellKey_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent periodic boundary key" << std::endl;
    return historyFile;
  }

  // If present, read unit cell records
  if(historyFile.supercellKey_ > 0)
  {
    Vector3D r_a1, r_a2, r_a3; // Lattice vectors read from file
    historyFile >> r_a1.x >> r_a1.y >> r_a1.z;
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    historyFile >> r_a2.x >> r_a2.y >> r_a2.z;
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    historyFile >> r_a3.x >> r_a3.y >> r_a3.z;
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    if(historyFile.fail())
    {
      std::cerr << "Error reading unit cell records from file " << historyFile.fileName() << std::endl;
      return historyFile;
    }
    if(historyFile.supercellKey_ == 6)
      system.setLattice(Lattice(r_a1, r_a2, system.lattice().type()));
    else
      system.setLattice(Lattice(r_a1, r_a2, r_a3, system.lattice().type()));
  }

  // Read atom records
  if(newSystem)
    system.add(Geometry()); // No previous info, will read to a single geometry

  size_t iMol = 0;       // Index of current molecule
  size_t iAtomInMol = 0; // Index of current atom in current molecule

  for(size_t iAtom = 0; iAtom < historyFile.numAtoms_; ++iAtom)
  {
    std::string r_symbol; // Atom symbol read from file
    size_t r_index;       // Atom index read from file
    Real r_x, r_y, r_z;   // Coordinates read from file

    // Read symbol and index
    historyFile >> r_symbol >> r_index;

    if(historyFile.fail())
    {
      std::cerr << "Error reading atom record from HISTORY file " << historyFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return historyFile;
    }
    if(r_index != iAtom + 1)
    {
      std::cerr << "Error reading HISTORY file " << historyFile.fileName()
                << ": atom indices not correlated" << std::endl;
      return historyFile;
    }
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // Read position
    historyFile >> r_x >> r_y >> r_z;

    if(historyFile.fail())
    {
      std::cerr << "Error reading atom position from file " << historyFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return historyFile;
    }
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // If present, skip velocity and/or force records
    for(size_t i = 0; i < historyFile.trajectoryKey_; ++i)
      historyFile.ignore(LINE_LENGTH, '\n');

    if(newSystem)
      system[0].addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z)));
    else
    {
      if(iMol > system.size() - 1 ||
         r_symbol.find(system[iMol].atom(iAtomInMol).symbol) == r_symbol.npos)
      {
        std::cerr << "Inconsistent data in HISTORY file " << historyFile.fileName() << std::endl;
        return historyFile;
      }
      system[iMol].setPosition(iAtomInMol, Vector3D(r_x, r_y, r_z));

      // Local index bookkeeping
      ++iAtomInMol;
      if(iAtomInMol > system[iMol].numAtoms() - 1)
      {
        ++iMol;
        iAtomInMol = 0;
      }
    }

  } // End of loop to read atom records

  // Consistency check
  if(!newSystem && iMol != system.size())
  {
    std::cerr << "Inconsistent number of atoms in HISTORY file " << historyFile.fileName() << std::endl;
    return historyFile;
  }

  // Update number of frames read
  ++historyFile.numFramesRead_;

  return historyFile;
}

HISTORYFile &operator>>(HISTORYFile &historyFile, System<Molecule> &system)
{
  assert(historyFile.is_open());
  bool newSystem = (system.size() == 0);

  if(!newSystem) // Consistency checks
  {
    size_t nAtoms = 0;
    for(size_t i = 0; i < system.size(); ++i)
      nAtoms += system[i].numAtoms();
    if(nAtoms != historyFile.numAtoms_)
    {
      std::cerr << "Error in HISTORYFile: Inconsistent number of atoms in system and history file." << std::endl;
      return historyFile;
    }
  }

  // Read configuration header record (mostly for consistency check)

  std::string line; // Line read from file

  // These variables named as in the DL_POLY manual
  std::string timestep; // The string "timestep"
  size_t nstep;         // Current time step
  size_t natms;         // Number of atoms
  size_t keytrj;        // Trajectory key
  size_t imcon;         // Periodic boundary key

  std::getline(historyFile, line);

  if(historyFile.eof()) return historyFile; // Reached end of the file, return

  std::stringstream headerStream(line);
  headerStream >> timestep;
  if(timestep.find("timestep") == timestep.npos)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": bad timestep record" << std::endl;
    return historyFile;
  }
  headerStream >> nstep >> natms >> keytrj >> imcon;
  if(headerStream.fail())
  {
    std::cerr << "Error reading timestep record from HISTORY file " 
              << historyFile.fileName() << std::endl;
    return historyFile;
  }
  if(natms != historyFile.numAtoms_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent number of atoms" << std::endl;
    return historyFile;
  }
  if(keytrj != historyFile.trajectoryKey_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent trajectory key" << std::endl;
    return historyFile;
  }
  if(imcon != historyFile.supercellKey_)
  {
    std::cerr << "Error reading HISTORY file " << historyFile.fileName()
              << ": inconsistent periodic boundary key" << std::endl;
    return historyFile;
  }

  // If present, read unit cell records
  if(historyFile.supercellKey_ > 0)
  {
    Vector3D r_a1, r_a2, r_a3; // Lattice vectors read from file
    historyFile >> r_a1.x >> r_a1.y >> r_a1.z;
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    historyFile >> r_a2.x >> r_a2.y >> r_a2.z;
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    historyFile >> r_a3.x >> r_a3.y >> r_a3.z;
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    if(historyFile.fail())
    {
      std::cerr << "Error reading unit cell records from file " << historyFile.fileName() << std::endl;
      return historyFile;
    }
    if(historyFile.supercellKey_ == 6)
      system.setLattice(Lattice(r_a1, r_a2, system.lattice().type()));
    else
      system.setLattice(Lattice(r_a1, r_a2, r_a3, system.lattice().type()));
  }

  // Read atom records
  if(newSystem)
    system.add(Molecule()); // No previous info, will read to a single molecule

  size_t iMol = 0;       // Index of current molecule
  size_t iAtomInMol = 0; // Index of current atom in current molecule

  for(size_t iAtom = 0; iAtom < historyFile.numAtoms_; ++iAtom)
  {
    std::string r_symbol; // Atom symbol read from file
    size_t r_index;       // Atom index read from file
    Real r_x, r_y, r_z;   // Coordinates read from file

    // Read symbol and index
    historyFile >> r_symbol >> r_index;

    if(historyFile.fail())
    {
      std::cerr << "Error reading atom record from HISTORY file " << historyFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return historyFile;
    }
    if(r_index != iAtom + 1)
    {
      std::cerr << "Error reading HISTORY file " << historyFile.fileName()
                << ": atom indices not correlated" << std::endl;
      return historyFile;
    }
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // Read position
    historyFile >> r_x >> r_y >> r_z;

    if(historyFile.fail())
    {
      std::cerr << "Error reading atom position from file " << historyFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return historyFile;
    }
    historyFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // If present, skip velocity and/or force records
    for(size_t i = 0; i < historyFile.trajectoryKey_; ++i)
      historyFile.ignore(LINE_LENGTH, '\n');

    if(newSystem)
      system[0].addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z)));
    else
    {
      if(iMol > system.size() - 1 ||
         r_symbol.find(system[iMol].atom(iAtomInMol).symbol) == r_symbol.npos)
      {
        std::cerr << "Inconsistent data in HISTORY file " << historyFile.fileName() << std::endl;
        return historyFile;
      }
      system[iMol].setPosition(iAtomInMol, Vector3D(r_x, r_y, r_z));

      // Local index bookkeeping
      ++iAtomInMol;
      if(iAtomInMol > system[iMol].numAtoms() - 1)
      {
        ++iMol;
        iAtomInMol = 0;
      }
    }

  } // End of loop to read atom records

  // Consistency check
  if(!newSystem && iMol != system.size())
  {
    std::cerr << "Inconsistent number of atoms in HISTORY file " << historyFile.fileName() << std::endl;
    return historyFile;
  }

  // Update number of frames read
  ++historyFile.numFramesRead_;

  return historyFile;
}

// Private functions

void HISTORYFile::initializeHISTORYFile()
{
  assert(is_open());
  
  // Read the title line
  std::getline(*this, title_);

  // Read the trajectory key, supercell key, and number of atoms

  std::string line; // Line read from file
  std::getline(*this, line);
  std::istringstream keyStream(line);
  keyStream >> trajectoryKey_;
  if(keyStream.fail() || trajectoryKey_ > 2)
  {
    std::cerr << "Error initializing HISTORY file " << fileName()
              << ": invalid keytrj record in header" << std::endl;
    clear();
    return;
  }

  keyStream >> supercellKey_;
  if(keyStream.fail() || supercellKey_ > 7)
  {
    std::cerr << "Error initializing HISTORY file " << fileName()
              << ": invalid imcon record in header" << std::endl;
    clear();
    return;
  }
  keyStream >> numAtoms_;
  if(keyStream.fail() || numAtoms_ < 1)
  {
    std::cerr << "Error initializing HISTORY file " << fileName()
              << ": invalid number of atoms in header" << std::endl;
    clear();
    return;
  }  
 
  // Done
  return;
}

