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
** CONFIG (DL-POLY) file format
*/

#include <iomanip>
#include <sstream>
#include "common/include/assert.h"
#include "mymol/include/file_formats/fileformats.h"

// Constructors

CONFIGFile::CONFIGFile()
:
IOFile()
{}

CONFIGFile::CONFIGFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

CONFIGFile &operator<<(CONFIGFile &configFile, Geometry const &geometry)
{
  assert(configFile.is_open() && configFile.mode() == OUT);

  // Set precision
  std::streamsize precision = configFile.precision();
  configFile.precision(WRITE_PRECISION);

  // Title line
  configFile.setf(std::ios::left, std::ios::adjustfield);
  configFile << std::setw(80) << geometry.name() << std::endl;

  // File key, supercell key, number of atoms
  configFile.setf(std::ios::right, std::ios::adjustfield);
  configFile << std::setw(10) << 0   // Only write positions
             << std::setw(10) << 0   // No PBCs
             << std::setw(10) << geometry.numAtoms()
             << std::endl;
  
 // Atom records
 for(size_t i = 0; i < geometry.numAtoms(); ++i)
 {
   Atom const &currentAtom = geometry.atom(i);
   configFile.setf(std::ios::left, std::ios::adjustfield);
   configFile << std::setw(8) << currentAtom.symbol;
   configFile.setf(std::ios::right, std::ios::adjustfield);
   configFile << std::setw(10) << i + 1
              << std::endl;
   configFile << std::setiosflags(std::ios::fixed)
              << std::setw(20) << currentAtom.position.x
              << std::setw(20) << currentAtom.position.y
              << std::setw(20) << currentAtom.position.z
              << std::endl;
 }

 configFile.precision(precision); // Return to original value

 return configFile;
}

CONFIGFile &operator>>(CONFIGFile &configFile, Geometry &geometry)
{
  assert(configFile.is_open() && configFile.mode() == IN);
  bool newGeometry = (geometry.numAtoms() == 0);

  // Read the title line
  std::string title;
  std::getline(configFile, title);
  geometry.setName(title);

  // Read the file key, supercell key, and number of atoms

  std::string line; // Line read from file

  // These variables follow notation from the DL_POLY manual 
  size_t levcfg; // File key
  size_t imcon;  // Periodic boundary key
  size_t natms;  // Number of atoms
  
  std::getline(configFile, line);
  std::istringstream keyStream(line);
  keyStream >> levcfg;
  if(keyStream.fail() || levcfg > 2)
  {
    std::cerr << "Error reading CONFIG file " << configFile.fileName()
              << ": invalid levcfg record" << std::endl;
    return configFile;
  }
  keyStream >> imcon;
  if(keyStream.fail() || imcon > 7)
  {
    std::cerr << "Error reading CONFIG file " << configFile.fileName()
              << ": invalid imcon record" << std::endl;
    return configFile;
  }
  keyStream >> natms;
  bool know_numAtoms = (!keyStream.fail());
  if(know_numAtoms && !newGeometry && geometry.numAtoms() != natms)
  {
    std::cerr << "Error reading CONFIG file " << configFile.fileName()
              << ": inconsistent number of atoms" << std::endl;
    return configFile;
  }
  
  // If present, skip unit cell records
  if(imcon > 0)
    for(size_t i = 0; i < 3; ++i) configFile.ignore(LINE_LENGTH, '\n');

  // Read atom records
  size_t iAtom = 0; // Index of current atom

  while(!configFile.eof())
  {
    std::string r_symbol; // Atom symbol read from file
    size_t r_index;       // Atom index read from file
    Real r_x, r_y, r_z;   // Coordinates read from file

    // Read symbol and index
    configFile >> r_symbol >> r_index;

    if(configFile.eof()) break; // Done reading
    if(configFile.fail())
    {
      std::cerr << "Error reading atom record from file " << configFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return configFile;
    }
    if(r_index != iAtom + 1)
    { 
      std::cerr << "Error reading CONFIG file " << configFile.fileName()
                << ": atom indices not correlated" << std::endl;
      return configFile;
    }
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
 
    // Read position
    configFile >> r_x >> r_y >> r_z;

    if(configFile.fail())
    {
      std::cerr << "Error reading atom position from file " << configFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return configFile;
    }
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // If present, skip velocity and/or force records
    for(size_t i = 0; i < levcfg; ++i)
      configFile.ignore(LINE_LENGTH, '\n');
    
    if(newGeometry)
      geometry.addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z))); 
    else
    {
      if(iAtom > geometry.numAtoms() - 1 ||
         r_symbol.find(geometry.atom(iAtom).symbol) == r_symbol.npos)
      {
        std::cerr << "Inconsistent data in CONFIG file " << configFile.fileName() << std::endl;
        return configFile;
      }
      geometry.setPosition(iAtom, Vector3D(r_x, r_y, r_z));   
    }
    ++iAtom;
  } // End of loop to read atom records

  if(!newGeometry && iAtom != geometry.numAtoms())
  {
    std::cerr << "Inconsistent number of atoms in CONFIG file " << configFile.fileName() << std::endl;
    return configFile;
  }

  // Find bonds automatically if this is a new geometry
  if(newGeometry) geometry.calculateBonds();

  return configFile;
}

CONFIGFile &operator<<(CONFIGFile &configFile, System<Geometry> const &system)
{
  assert(configFile.is_open() && configFile.mode() == OUT);

  // Sanity check
  if(system.size() == 0 || system[0].numAtoms() == 0)
  {
    std::cerr << "Error while writing CONFIG file " << configFile.fileName()
              << ": empty system" << std::endl;
    return configFile;
  }

  // Get number of atoms
  size_t nAtoms = 0;
  for(size_t i = 0; i < system.size(); ++i)
    nAtoms += system[i].numAtoms();

  // Get periodic boundary key (imcon)  
  size_t imcon = 0; // Following notation in DL_POLY manual
  if(system.lattice().numDimensions() == 2) imcon = 6;
  if(system.lattice().numDimensions() == 3) imcon = 3;

  // Set precision
  std::streamsize precision = configFile.precision();
  configFile.precision(WRITE_PRECISION);

  // Title line
  configFile.setf(std::ios::left, std::ios::adjustfield);
  configFile << std::setw(80) << system.name() << std::endl;

  // File key, supercell key, number of atoms
  configFile.setf(std::ios::right, std::ios::adjustfield);

  configFile << std::setw(10) << 0   // Only write positions
             << std::setw(10) << imcon
             << std::setw(10) << nAtoms
             << std::endl;

  // Supercell
  if(imcon != 0)
  {
    Vector3D const a = system.lattice().latticeVector(0);
    Vector3D const b = system.lattice().latticeVector(1);
    Vector3D const c = (imcon = 3)?system.lattice().latticeVector(2):0.0; 
   
    configFile << std::setiosflags(std::ios::fixed) 
               << std::setw(20) << a.x
               << std::setw(20) << a.y
               << std::setw(20) << a.z
               << std::endl
               << std::setw(20) << b.x
               << std::setw(20) << b.y
               << std::setw(20) << b.z
               << std::endl
               << std::setw(20) << c.x
               << std::setw(20) << c.y
               << std::setw(20) << c.z
               << std::endl;
  } 
 
  // Atom records
  size_t iAtom = 0;
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j) 
  {
    Atom const &currentAtom = system[i].atom(j);
    configFile.setf(std::ios::left, std::ios::adjustfield);
    configFile << std::setw(8) << currentAtom.symbol;
    configFile.setf(std::ios::right, std::ios::adjustfield);
    configFile << std::setw(10) << iAtom + 1
               << std::endl;
    configFile << std::setiosflags(std::ios::fixed) 
               << std::setw(20) << currentAtom.position.x
               << std::setw(20) << currentAtom.position.y
               << std::setw(20) << currentAtom.position.z
               << std::endl;
    ++iAtom;
  }

  configFile.precision(precision); // Return to original value

  return configFile;
}

CONFIGFile &operator>>(CONFIGFile &configFile, System<Geometry> &system)
{
  assert(configFile.is_open() && configFile.mode() == IN);
  bool newSystem = (system.size() == 0);

  // Read the title line
  std::string title;
  std::getline(configFile, title);
  system.setName(title);

  // Read the file key, supercell key, and number of atoms

  std::string line; // Line read from file

  // These variables follow notation from the DL_POLY manual 
  size_t levcfg; // File key
  size_t imcon;  // Periodic boundary key
  size_t natms;  // Number of atoms
  
  std::getline(configFile, line);
  std::istringstream keyStream(line);
  keyStream >> levcfg;
  if(keyStream.fail() || levcfg > 2)
  {
    std::cerr << "Error reading CONFIG file " << configFile.fileName()
              << ": invalid levcfg record" << std::endl;
    return configFile;
  }
  keyStream >> imcon;
  if(keyStream.fail() || imcon > 7)
  {
    std::cerr << "Error reading CONFIG file " << configFile.fileName()
              << ": invalid imcon record" << std::endl;
    return configFile;
  }
  keyStream >> natms;
  bool know_numAtoms = (!keyStream.fail());
  bool know_lattice = (system.lattice().numDimensions() > 0);

  // Consistency checks
  if(know_lattice && !newSystem) // In case lattice has been read before
  {
    if( (imcon == 0 && system.lattice().numDimensions() > 0) ||
        (imcon == 6 && system.lattice().numDimensions() != 2) ||
        (imcon != 0 && imcon != 6 && system.lattice().numDimensions() < 3) )
    {
      std::cerr << "Error reading CONFIG file " << configFile.fileName()
                << ": inconsistent unit cell" << std::endl;
      return configFile;
    }
  }
  if(know_numAtoms && !newSystem) 
  {
    size_t nAtoms = 0;
    for(size_t i = 0; i < system.size(); ++i)
      nAtoms += system[i].numAtoms();
    
    if(nAtoms != natms)
    {
      std::cerr << "Error reading CONFIG file " << configFile.fileName()
                << ": inconsistent number of atoms" << std::endl;
      return configFile;
    }
  }

  // If present, read unit cell records
  if(imcon > 0)
  {
    Vector3D r_a1, r_a2, r_a3; // Lattice vectors read from file  
    configFile >> r_a1.x >> r_a1.y >> r_a1.z;
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    configFile >> r_a2.x >> r_a2.y >> r_a2.z;
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    configFile >> r_a3.x >> r_a3.y >> r_a3.z;
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    if(configFile.fail())
    {
      std::cerr << "Error reading unit cell records from file " << configFile.fileName() << std::endl;
      return configFile;
    }
    if(imcon == 6)
      system.setLattice(Lattice(r_a1, r_a2, system.lattice().type()));
    else
      system.setLattice(Lattice(r_a1, r_a2, r_a3, system.lattice().type()));
  }

  // Read atom records
  size_t iAtom = 0;      // Index of current atom
  size_t iMol = 0;       // Index of current molecule
  size_t iAtomInMol = 0; // Index of current atom in current molecule

  if(newSystem)
    system.add(Geometry()); // No connectivity, will read to a single geometry

  while(!configFile.eof())
  {
    std::string r_symbol; // Atom symbol read from file
    size_t r_index;       // Atom index read from file
    Real r_x, r_y, r_z;   // Coordinates read from file

    // Read symbol and index
    configFile >> r_symbol >> r_index;

    if(configFile.eof()) break; // Done reading
    if(configFile.fail())
    {
      std::cerr << "Error reading atom record from file " << configFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return configFile;
    }
    if(r_index != iAtom + 1)
    { 
      std::cerr << "Error reading CONFIG file " << configFile.fileName()
                << ": atom indices not correlated" << std::endl;
      return configFile;
    }
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
 
    // Read position
    configFile >> r_x >> r_y >> r_z;

    if(configFile.fail())
    {
      std::cerr << "Error reading atom position from file " << configFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return configFile;
    }
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // If present, skip velocity and/or force records
    for(size_t i = 0; i < levcfg; ++i)
      configFile.ignore(LINE_LENGTH, '\n');

    if(newSystem)
      system[0].addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z)));
    else
    {
      if(iMol > system.size() - 1 || 
         r_symbol.find(system[iMol].atom(iAtomInMol).symbol) == r_symbol.npos)
      {
        std::cerr << "Inconsistent data in CONFIG file " << configFile.fileName() << std::endl;
        return configFile;
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
    ++iAtom;
  }  // End of loop to read atom positions

  // Consistency check
  if(!newSystem && iMol != system.size())
  {
    std::cerr << "Inconsistent number of atoms in CONFIG file " << configFile.fileName() << std::endl;
    return configFile;
  } 
 
  // Find bonds automatically if this is a new system (this does not separate molecules)
  if(newSystem) system[0].calculateBonds();

  return configFile;
}

CONFIGFile &operator<<(CONFIGFile &configFile, System<Molecule> const &system)
{
  assert(configFile.is_open() && configFile.mode() == OUT);

  // Sanity check
  if(system.size() == 0 || system[0].numAtoms() == 0)
  {
    std::cerr << "Error while writing CONFIG file " << configFile.fileName()
              << ": empty system" << std::endl;
    return configFile;
  }

  // Get number of atoms
  size_t nAtoms = 0;
  for(size_t i = 0; i < system.size(); ++i)
    nAtoms += system[i].numAtoms();

  // Get periodic boundary key (imcon)  
  size_t imcon = 0; // Following notation in DL_POLY manual
  if(system.lattice().numDimensions() == 2) imcon = 6;
  if(system.lattice().numDimensions() == 3) imcon = 3;

  // Set precision
  std::streamsize precision = configFile.precision();
  configFile.precision(WRITE_PRECISION);

  // Title line
  configFile.setf(std::ios::left, std::ios::adjustfield);
  configFile << std::setw(80) << system.name() << std::endl;

  // File key, supercell key, number of atoms
  configFile.setf(std::ios::right, std::ios::adjustfield);

  configFile << std::setw(10) << 0   // Only write positions
             << std::setw(10) << imcon
             << std::setw(10) << nAtoms
             << std::endl;

  // Supercell
  if(imcon != 0)
  {
    Vector3D const a = system.lattice().latticeVector(0);
    Vector3D const b = system.lattice().latticeVector(1);
    Vector3D const c = (imcon = 3)?system.lattice().latticeVector(2):0.0; 
   
    configFile << std::setiosflags(std::ios::fixed) 
               << std::setw(20) << a.x
               << std::setw(20) << a.y
               << std::setw(20) << a.z
               << std::endl
               << std::setw(20) << b.x
               << std::setw(20) << b.y
               << std::setw(20) << b.z
               << std::endl
               << std::setw(20) << c.x
               << std::setw(20) << c.y
               << std::setw(20) << c.z
               << std::endl;
  } 
 
  // Atom records
  size_t iAtom = 0;
  for(size_t i = 0; i < system.size(); ++i)
  for(size_t j = 0; j < system[i].numAtoms(); ++j) 
  {
    Atom const &currentAtom = system[i].atom(j);
    configFile.setf(std::ios::left, std::ios::adjustfield);
    configFile << std::setw(8) << currentAtom.symbol;
    configFile.setf(std::ios::right, std::ios::adjustfield);
    configFile << std::setw(10) << iAtom + 1
               << std::endl;
    configFile << std::setiosflags(std::ios::fixed) 
               << std::setw(20) << currentAtom.position.x
               << std::setw(20) << currentAtom.position.y
               << std::setw(20) << currentAtom.position.z
               << std::endl;
    ++iAtom;
  }

  configFile.precision(precision); // Return to original value

  return configFile;
}

CONFIGFile &operator>>(CONFIGFile &configFile, System<Molecule> &system)
{
  assert(configFile.is_open() && configFile.mode() == IN);
  bool newSystem = (system.size() == 0);

  // Read the title line
  std::string title;
  std::getline(configFile, title);
  system.setName(title);

  // Read the file key, supercell key, and number of atoms

  std::string line; // Line read from file

  // These variables follow notation from the DL_POLY manual 
  size_t levcfg; // File key
  size_t imcon;  // Periodic boundary key
  size_t natms;  // Number of atoms
  
  std::getline(configFile, line);
  std::istringstream keyStream(line);
  keyStream >> levcfg;
  if(keyStream.fail() || levcfg > 2)
  {
    std::cerr << "Error reading CONFIG file " << configFile.fileName()
              << ": invalid levcfg record" << std::endl;
    return configFile;
  }
  keyStream >> imcon;
  if(keyStream.fail() || imcon > 7)
  {
    std::cerr << "Error reading CONFIG file " << configFile.fileName()
              << ": invalid imcon record" << std::endl;
    return configFile;
  }
  keyStream >> natms;
  bool know_numAtoms = (!keyStream.fail());
  bool know_lattice = (system.lattice().numDimensions() > 0);

  // Consistency checks
  if(know_lattice && !newSystem) // In case lattice has been read before
  {
    if( (imcon == 0 && system.lattice().numDimensions() > 0) ||
        (imcon == 6 && system.lattice().numDimensions() != 2) ||
        (imcon != 0 && imcon != 6 && system.lattice().numDimensions() < 3) )
    {
      std::cerr << "Error reading CONFIG file " << configFile.fileName()
                << ": inconsistent unit cell" << std::endl;
      return configFile;
    }
  }
  if(know_numAtoms && !newSystem) 
  {
    size_t nAtoms = 0;
    for(size_t i = 0; i < system.size(); ++i)
      nAtoms += system[i].numAtoms();
    
    if(nAtoms != natms)
    {
      std::cerr << "Error reading CONFIG file " << configFile.fileName()
                << ": inconsistent number of atoms" << std::endl;
      return configFile;
    }
  }

  // If present, read unit cell records
  if(imcon > 0)
  {
    Vector3D r_a1, r_a2, r_a3; // Lattice vectors read from file  
    configFile >> r_a1.x >> r_a1.y >> r_a1.z;
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    configFile >> r_a2.x >> r_a2.y >> r_a2.z;
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    configFile >> r_a3.x >> r_a3.y >> r_a3.z;
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
    if(configFile.fail())
    {
      std::cerr << "Error reading unit cell records from file " << configFile.fileName() << std::endl;
      return configFile;
    }
    if(imcon == 6)
      system.setLattice(Lattice(r_a1, r_a2, system.lattice().type()));
    else
      system.setLattice(Lattice(r_a1, r_a2, r_a3, system.lattice().type()));
  }

  // Read atom records
  size_t iAtom = 0;      // Index of current atom
  size_t iMol = 0;       // Index of current molecule
  size_t iAtomInMol = 0; // Index of current atom in current molecule

  if(newSystem)
    system.add(Molecule()); // No connectivity, will read to a single molecule

  while(!configFile.eof())
  {
    std::string r_symbol; // Atom symbol read from file
    size_t r_index;       // Atom index read from file
    Real r_x, r_y, r_z;   // Coordinates read from file

    // Read symbol and index
    configFile >> r_symbol >> r_index;

    if(configFile.eof()) break; // Done reading
    if(configFile.fail())
    {
      std::cerr << "Error reading atom record from file " << configFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return configFile;
    }
    if(r_index != iAtom + 1)
    { 
      std::cerr << "Error reading CONFIG file " << configFile.fileName()
                << ": atom indices not correlated" << std::endl;
      return configFile;
    }
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line
 
    // Read position
    configFile >> r_x >> r_y >> r_z;

    if(configFile.fail())
    {
      std::cerr << "Error reading atom position from file " << configFile.fileName() << std::endl
                << "Atoms read: " << iAtom << std::endl;
      return configFile;
    }
    configFile.ignore(LINE_LENGTH, '\n'); // Skip to next line

    // If present, skip velocity and/or force records
    for(size_t i = 0; i < levcfg; ++i)
      configFile.ignore(LINE_LENGTH, '\n');

    if(newSystem)
      system[0].addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z)));
    else
    {
      if(iMol > system.size() - 1 || 
         r_symbol.find(system[iMol].atom(iAtomInMol).symbol) == r_symbol.npos)
      {
        std::cerr << "Inconsistent data in CONFIG file " << configFile.fileName() << std::endl;
        return configFile;
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
    ++iAtom;
  }  // End of loop to read atom positions

  // Consistency check
  if(!newSystem && iMol != system.size())
  {
    std::cerr << "Inconsistent number of atoms in CONFIG file " << configFile.fileName() << std::endl;
    return configFile;
  } 
 
  // Find bonds automatically if this is a new system (this does not separate molecules)
  if(newSystem) system[0].calculateBonds();

  return configFile;
}

