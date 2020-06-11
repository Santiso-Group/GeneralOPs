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
** .psf file format
*/

#include <sstream>
#include <vector>
#include "common/include/assert.h"
#include "mymol/include/file_formats/fileformats.h"

// Constructors

PSFFile::PSFFile()
:
IOFile()
{}

PSFFile::PSFFile(std::string const &fileName)
:
IOFile(fileName, IN)
{}

// Operators

PSFFile &operator>>(PSFFile &psfFile, Geometry &geometry)
{
  assert(psfFile.is_open());
  geometry.clear();

// Search for title record

  std::string line;     // Line read from file
  size_t r_numResidues; // Number of residues (?) read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NTITLE") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find title record." << std::endl;
    return psfFile;
  }

// Get the number of residues (?)

  std::istringstream titleLineStream(line);
  titleLineStream >> r_numResidues;

// Search for atom records

  size_t r_numAtoms;  // Number of atoms read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NATOM") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find atom records." << std::endl;
    return psfFile;
  }

// Get the number of atoms

  std::istringstream atomLineStream(line);
  atomLineStream >> r_numAtoms;

// Read atom records

  while(!psfFile.eof() && geometry.numAtoms() < r_numAtoms)
  {
    size_t r_atomID, r_resID; // Atom and residue IDs read from file
    std::string r_segName, r_resName; // Segment and residue name read from file
    std::string r_atomName, r_atomType;  // Atom name and type read from file
    psfFile >> r_atomID >> r_segName >> r_resID >> r_resName >> r_atomName >> r_atomType;
    psfFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    if(!psfFile)
    {
      std::cerr << "Error reading atom record from file " << psfFile.fileName() << std::endl;
      return psfFile;
    }
    std::string symbol = r_atomType.substr(0, 1); // Only for 1-letter symbols
    // Make consistent with pdb role names - note this works only for 1-letter atomic symbols
    if(r_atomName.size() == 0 ||
       r_atomName.find_first_not_of(" ") == r_atomName.npos)
        r_atomName = symbol;
    if(r_atomName.size() < 4) r_atomName = " " + r_atomName;
    while(r_atomName.size() < 4) r_atomName = r_atomName + " ";
    if(r_atomName.size() > 4) r_atomName = r_atomName.substr(0, 4);
    geometry.addAtom(Atom(symbol, r_atomType, 0.0, r_resID, r_resName, r_atomName));
  }
  if(geometry.numAtoms() < r_numAtoms)
  {
    std::cerr << "Error reading file " << psfFile.fileName() << std::endl
              << "Expected " << r_numAtoms << " atoms but found only "
              << geometry.numAtoms() << std::endl;
    geometry.clear();
    return psfFile;
  }

// Search for bond records

  size_t r_numBonds;  // Number of bonds read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NBOND") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find bond records." << std::endl;
    return psfFile;
  }

// Get the number of bonds

  std::istringstream bondLineStream(line);
  bondLineStream >> r_numBonds;

// Read bond records

  while(!psfFile.eof() && geometry.numBonds() < r_numBonds)
  {
    size_t r_first;         // First atom read from input file
    size_t r_second;        // Second atom read from input file

    psfFile >> r_first >> r_second;
    if(!psfFile)
    {
      std::cerr << "Error reading bond record from file " << psfFile.fileName() << std::endl;
      geometry.clear();
      return psfFile;
    }
    if(r_first > r_second) {  // Make first < second
        size_t temp = r_first;
        r_first = r_second;
        r_second = temp;
    }
    geometry.addBond(Bond(r_first - 1, r_second - 1));
    if(geometry.numBonds()%4 == 0) // If 4 bonds have been read, skip to next line
      psfFile.ignore(LINE_LENGTH, '\n'); 
  }
  if(geometry.numBonds() < r_numBonds)
  {
    std::cerr << "Error reading file " << psfFile.fileName() << std::endl
              << "Expected " << r_numBonds << " bonds but found only "
              << geometry.numBonds() << std::endl;
    geometry.clear();
    return psfFile;
  }
  return psfFile;
}

PSFFile &operator>>(PSFFile &psfFile, System<Geometry> &system)
{
  assert(psfFile.is_open());
  system.clear();

// Search for title record

  std::string line;     // Line read from file
  size_t r_numResidues; // Number of residues (?) read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NTITLE") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find title record." << std::endl;
    return psfFile;
  }

// Get the number of residues (?)

  std::istringstream titleLineStream(line);
  titleLineStream >> r_numResidues;

// Search for atom records

  size_t r_numAtoms;  // Number of atoms read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NATOM") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find atom records." << std::endl;
    return psfFile;
  }

// Get the number of atoms

  std::istringstream atomLineStream(line);
  atomLineStream >> r_numAtoms;

// Read atom records

  Geometry currentResidue;              // Current geometry
  size_t nAtoms = 0;                    // Total number of atoms read
  size_t previousID = 0;                // Previous residue ID read
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each geometry
                                        // (needed to convert local to global indices)

  while(!psfFile.eof() && nAtoms < r_numAtoms)
  {
    size_t r_atomID, r_resID;           // Atom and residue IDs read from file
    std::string r_segName, r_resName;   // Segment and residue name read from file
    std::string r_atomName, r_atomType; // Atom name and type read from file

    psfFile >> r_atomID >> r_segName >> r_resID >> r_resName >> r_atomName >> r_atomType;
    psfFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    if(!psfFile)
    {
      std::cerr << "Error reading atom record from file " << psfFile.fileName() << std::endl;
      return psfFile;
    }
    std::string symbol = r_atomType.substr(0, 1); // Only for 1-letter symbols
    // Make consistent with pdb role names - note this works only for 1-letter atomic symbols
    if(r_atomName.size() == 0 ||
       r_atomName.find_first_not_of(" ") == r_atomName.npos)
        r_atomName = symbol;
    if(r_atomName.size() < 4) r_atomName = " " + r_atomName;
    while(r_atomName.size() < 4) r_atomName = r_atomName + " ";
    if(r_atomName.size() > 4) r_atomName = r_atomName.substr(0, 4);
    if(!currentResidue.name().size()) currentResidue.setName(r_resName);
    if(nAtoms > 0 && r_resID != previousID)
    {
      system.add(currentResidue);
      firstAtoms.push_back(nAtoms);
      currentResidue.clear();
    }
    nAtoms++;
    currentResidue.addAtom(Atom(symbol, r_atomType, 0.0, r_resID, r_resName, r_atomName));
    previousID = r_resID;
  }
  system.add(currentResidue);
  firstAtoms.push_back(nAtoms);
  if(nAtoms < r_numAtoms)
  {
    std::cerr << "Error reading file " << psfFile.fileName() << std::endl
              << "Expected " << r_numAtoms << " atoms but found only "
              << nAtoms << std::endl;
    system.clear();
    return psfFile;
  }

// Search for bond records

  size_t r_numBonds;  // Number of bonds read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NBOND") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find bond records." << std::endl;
    return psfFile;
  }

// Get the number of bonds

  std::istringstream bondLineStream(line);
  bondLineStream >> r_numBonds;

// Read bond records

  size_t nBonds = 0;

  while(!psfFile.eof() && nBonds < r_numBonds)
  {
    size_t r_first;         // First atom read from input file
    size_t r_second;        // Second atom read from input file

    psfFile >> r_first >> r_second;
    if(!psfFile)
    {
      std::cerr << "Error reading bond record from file " << psfFile.fileName() << std::endl;
      system.clear();
      return psfFile;
    }
    if(r_first > r_second) {  // Make first < second
        size_t temp = r_first;
        r_first = r_second;
        r_second = temp;
    }
    // Find which geometry the bond belongs to
    size_t geometryIndex = 0;
    while(r_first > firstAtoms[geometryIndex + 1])
      geometryIndex++;
    if(r_second <= firstAtoms[geometryIndex + 1]) // Must be intra-geometry
      system[geometryIndex].addBond(Bond(r_first - 1 - firstAtoms[geometryIndex],
                                         r_second - 1 - firstAtoms[geometryIndex]));
    nBonds++;
    if(nBonds%4 == 0) // If 4 bonds have been read, skip to next line
      psfFile.ignore(LINE_LENGTH, '\n'); 
  }
  if(nBonds < r_numBonds)
  {
    std::cerr << "Error reading file " << psfFile.fileName() << std::endl
              << "Expected " << r_numBonds << " bonds but found only "
              << nBonds << std::endl;
    system.clear();
  }
  return psfFile;
}

PSFFile &operator>>(PSFFile &psfFile, System<Molecule> &system)
{
  assert(psfFile.is_open());
  system.clear();

// Search for title record

  std::string line;     // Line read from file
  size_t r_numResidues; // Number of residues (?) read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NTITLE") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find title record." << std::endl;
    return psfFile;
  }

// Get the number of residues

  std::istringstream titleLineStream(line);
  titleLineStream >> r_numResidues;

// Search for atom records

  size_t r_numAtoms;  // Number of atoms read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NATOM") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find atom records." << std::endl;
    return psfFile;
  }

// Get the number of atoms

  std::istringstream atomLineStream(line);
  atomLineStream >> r_numAtoms;

// Read atom records

  Molecule currentResidue;              // Current molecule
  size_t nAtoms = 0;                    // Total number of atoms read
  size_t previousID = 0;                // Previous residue ID read
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each molecule
                                        // (needed to convert local to global indices)

  while(!psfFile.eof() && nAtoms < r_numAtoms)
  {
    size_t r_atomID, r_resID;           // Atom and residue IDs read from file
    std::string r_segName, r_resName;   // Segment and residue name read from file
    std::string r_atomName, r_atomType; // Atom name and type read from file

    psfFile >> r_atomID >> r_segName >> r_resID >> r_resName >> r_atomName >> r_atomType;
    psfFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    if(!psfFile)
    {
      std::cerr << "Error reading atom record from file " << psfFile.fileName() << std::endl;
      return psfFile;
    }
    std::string symbol = r_atomType.substr(0, 1); // Only for 1-letter symbols
    // Make consistent with pdb role names - note this works only for 1-letter atomic symbols
    if(r_atomName.size() == 0 ||
       r_atomName.find_first_not_of(" ") == r_atomName.npos)
        r_atomName = symbol;
    if(r_atomName.size() < 4) r_atomName = " " + r_atomName;
    while(r_atomName.size() < 4) r_atomName = r_atomName + " ";
    if(r_atomName.size() > 4) r_atomName = r_atomName.substr(0, 4);
    if(!currentResidue.name().size()) currentResidue.setName(r_resName);
    if(nAtoms > 0 && r_resID != previousID)
    {
      system.add(currentResidue);
      firstAtoms.push_back(nAtoms);
      currentResidue.clear();
    }
    nAtoms++;
    currentResidue.addAtom(Atom(symbol, r_atomType, 0.0, r_resID, r_resName, r_atomName));
    previousID = r_resID;
  }
  system.add(currentResidue);
  firstAtoms.push_back(nAtoms);
  if(nAtoms < r_numAtoms)
  {
    std::cerr << "Error reading file " << psfFile.fileName() << std::endl
              << "Expected " << r_numAtoms << " atoms but found only "
              << nAtoms << std::endl;
    system.clear();
    return psfFile;
  }

// Search for bond records

  size_t r_numBonds;  // Number of bonds read from file

  do
  {
    std::getline(psfFile, line);
    if(line.find("!NBOND") != line.npos) break;
  }
  while(!psfFile.eof());

  if(psfFile.eof())
  {
    std::cerr << "Error reading file " << psfFile.fileName() 
              << ": could not find bond records." << std::endl;
    return psfFile;
  }

// Get the number of bonds

  std::istringstream bondLineStream(line);
  bondLineStream >> r_numBonds;

// Read bond records

  size_t nBonds = 0;

  while(!psfFile.eof() && nBonds < r_numBonds)
  {
    size_t r_first;         // First atom read from input file
    size_t r_second;        // Second atom read from input file

    psfFile >> r_first >> r_second;
    if(!psfFile)
    {
      std::cerr << "Error reading bond record from file " << psfFile.fileName() << std::endl;
      system.clear();
      return psfFile;
    }
    if(r_first > r_second) {  // Make first < second
        size_t temp = r_first;
        r_first = r_second;
        r_second = temp;
    }
    // Find which geometry the bond belongs to
    size_t geometryIndex = 0;
    while(r_first > firstAtoms[geometryIndex + 1])
      geometryIndex++;
    if(r_second <= firstAtoms[geometryIndex + 1]) // Must be intra-molecule
      system[geometryIndex].addBond(Bond(r_first - 1 - firstAtoms[geometryIndex],
                                         r_second - 1 - firstAtoms[geometryIndex]));
    nBonds++;
    if(nBonds%4 == 0) // If 4 bonds have been read, skip to next line
      psfFile.ignore(LINE_LENGTH, '\n'); 
  }
  if(nBonds < r_numBonds)
  {
    std::cerr << "Error reading file " << psfFile.fileName() << std::endl
              << "Expected " << r_numBonds << " bonds but found only "
              << nBonds << std::endl;
    system.clear();
  }
  return psfFile;
}
