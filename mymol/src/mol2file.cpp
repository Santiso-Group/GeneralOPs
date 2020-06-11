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
** .mol2 (Sybyl) file format
*/

#include <iomanip>
#include <cmath>
#include "common/include/assert.h"
#include "mymol/include/file_formats/fileformats.h"

// Constructors

MOL2File::MOL2File()
:
IOFile()
{}

MOL2File::MOL2File(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

MOL2File &operator<<(MOL2File &mol2File, Geometry const &geometry)
{
  // Write geometry to file

  assert(mol2File.is_open() && mol2File.mode() == OUT);
  mol2File.setf(std::ios::fixed, std::ios::floatfield);
  mol2File.setf(std::ios::right, std::ios::adjustfield);
  size_t const nAtoms = geometry.numAtoms();
  size_t const nBonds = geometry.numBonds();
  size_t const intWidth = (nAtoms > nBonds)?(2 + (size_t)log10((Real)nAtoms)):(2 + (size_t)log10((Real)nBonds));

  // Molecule record

  mol2File << "@<TRIPOS>MOLECULE" << std::endl;
  mol2File << geometry.name() << std::endl;
  mol2File << nAtoms << " " << nBonds << std::endl;
  mol2File << "SMALL" << std::endl;
  mol2File << "NO_CHARGES" << std::endl;

  mol2File << std::endl << std::endl;

  // Atom records

  std::streamsize precision = mol2File.precision();
  mol2File.precision(WRITE_PRECISION);

  mol2File << "@<TRIPOS>ATOM" << std::endl;

  for(size_t i = 0; i < nAtoms; i++)
  {
    Atom atom(geometry.atom(i));
    mol2File.setf(std::ios::right, std::ios::adjustfield);
    mol2File << std::setw(intWidth) << i + 1 << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield); 
    mol2File << std::setw(MAX_SYMBOL_LENGTH) 
             << atom.symbol << " ";
    mol2File.setf(std::ios::right, std::ios::adjustfield); 
    for(size_t j = 0; j < 3; j++) 
      mol2File << std::setw(WRITE_PRECISION + 6) 
               << atom.position[j] << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield);
    mol2File << std::setw(MAX_TYPE_LENGTH) 
             << atom.type << " ";
    if(atom.residueName.size() > 0)
    {
      mol2File.setf(std::ios::right, std::ios::adjustfield);
      mol2File << std::setw(intWidth) << atom.residueID << " ";
      mol2File.setf(std::ios::left, std::ios::adjustfield);
      mol2File << std::setw(MAX_RESNAME_LENGTH)
               << atom.residueName << " ";
    }
    mol2File << std::endl;
  }

  // Bond records

  mol2File << "@<TRIPOS>BOND" << std::endl;

  for(size_t i = 0; i < nBonds; i++)
  {
    Bond bond(geometry.bond(i));
    mol2File.setf(std::ios::right, std::ios::adjustfield);
    mol2File << std::setw(intWidth) 
             << i+1 << " " 
             << bond.first + 1 << " " 
             << bond.second + 1 << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield);
    mol2File << bond.type << std::endl;
  }

  mol2File << std::endl;

  mol2File.precision(precision); // Return to original value

  return mol2File;
}

MOL2File &operator>>(MOL2File &mol2File, Geometry &geometry)
{
  // Read geometry from file

  assert(mol2File.is_open() && mol2File.mode() == IN);
  geometry.clear();
  
  // Search for molecule record

  do
  {
    std::string line;
    std::getline(mol2File, line);
    if(line.find("@<TRIPOS>MOLECULE") != line.npos) break;
  }
  while(!mol2File.eof());

  if(mol2File.eof())
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": MOLECULE record not found." << std::endl;
    return mol2File;
  }

  // Read molecule record

  size_t r_nAtom;             // Number of atoms read from input file
  size_t r_nBond;             // Number of bonds read from input file
  char r_label[LINE_LENGTH];  // Geometry label read from input file

  mol2File.getline(r_label, LINE_LENGTH); // Read geometry label
  geometry.setName(std::string(r_label));

  mol2File >> r_nAtom >> r_nBond;     // Read number of atoms and number of bonds
  mol2File.ignore(LINE_LENGTH, '\n'); // Skip to next line

  if(!mol2File)
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": invalid MOLECULE record." << std::endl;
    geometry.clear();
    return mol2File;
  }

  // Search for atom record

  do
  {
    std::string line;
    std::getline(mol2File, line);
    if(line.find("@<TRIPOS>ATOM") != line.npos) break;
  }
  while(!mol2File.eof());

  if(mol2File.eof())
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": ATOM record not found." << std::endl;
    geometry.clear();
    return mol2File;
  }

  // Read atom information

  size_t r_atomID;            // Atom id number read from input file
  std::string r_symbol;       // Atomic symbol read from input file
  Vector3D r_position;        // Position read from input file
  std::string r_atomType;     // Atomic type read from input file
  size_t r_residueID;         // Residue ID number read from input file
  std::string r_residueName;  // Residue name read from input file
  std::string r_residueInfo;  // Residue info read from input file
  bool thereIsResidueInfo;    // Whether residue info is present or not

  while(!mol2File.eof() && geometry.numAtoms() < r_nAtom)
  {
    mol2File >> r_atomID >> r_symbol >> r_position.x >> r_position.y >> r_position.z >> r_atomType;
    getline(mol2File, r_residueInfo, '\n'); // Check for residue info
    thereIsResidueInfo = (r_residueInfo.find_first_not_of(" \t\n") != r_residueInfo.npos);
    if(thereIsResidueInfo)  // Read residue info
    {
      std::istringstream resInfo(r_residueInfo);
      resInfo >> r_residueID >> r_residueName;
    }
    if(!mol2File)
    {
      std::cerr << "Error reading file " << mol2File.fileName() 
                << ": invalid ATOM record." << std::endl
                << "Atoms read: " << geometry.numAtoms() << std::endl;
      geometry.clear();
      return mol2File;
    }
    if(thereIsResidueInfo)
      geometry.addAtom(Atom(r_symbol, r_atomType, r_position, r_residueID, r_residueName));
    else
      geometry.addAtom(Atom(r_symbol, r_atomType, r_position));
  }
  if(geometry.numAtoms() < r_nAtom)
  {
    std::cerr << "Error reading file " << mol2File.fileName() << std::endl
              << "Expected " << r_nAtom << " atoms but found only " 
              << geometry.numAtoms() << std::endl;
    geometry.clear();
    return mol2File;
  }

  // Search for bond and unit cell records

  bool foundBonds = false;
  bool foundCell = false;
  do
  {
    std::string line;
    std::getline(mol2File, line);
    foundBonds = (line.find("@<TRIPOS>BOND") != line.npos);
    foundCell = (line.find("@<TRIPOS>CRYSIN") != line.npos);
    if(foundBonds || foundCell) break;
  }
  while(!mol2File.eof());

  if(foundBonds)
  {
    // Read bond information

    size_t r_bondID;        // Bond id number read from input file
    size_t r_first;         // First atom read from input file
    size_t r_second;        // Second atom read from input file
    std::string r_bondType; // Bond type read from input file

    while(!mol2File.eof() && geometry.numBonds() < r_nBond)
    {
      mol2File >> r_bondID >> r_first >> r_second >> r_bondType;
      if(r_first > r_nAtom || r_second > r_nAtom)
      {
        std::cerr << "Error reading file " << mol2File.fileName() 
                  << ": invalid bond id number in BOND record" << std::endl;
        geometry.clear();
        return mol2File;
      }
      if(r_first > r_second) {  // Make first < second
        size_t temp = r_first;
        r_first = r_second;
        r_second = temp;
      }
      mol2File.ignore(LINE_LENGTH, '\n'); // Skip to next line
      if(!mol2File)
      {
        std::cerr << "Error reading file " << mol2File.fileName() 
                  << ": invalid BOND record." << std::endl
                  << "Bonds read: " << geometry.numBonds() << std::endl;
        geometry.clear();
        return mol2File;
      }
      geometry.addBond(Bond(r_first - 1, r_second - 1, r_bondType));
    }
    if(geometry.numBonds() < r_nBond)
    {
      std::cerr << "Error reading file " << mol2File.fileName() << std::endl;
      std::cerr << "Expected " << r_nBond << " bonds but found only " 
                << geometry.numBonds() << std::endl;
      geometry.clear();
    }
    // Search for unit cell record
    do
    {
      std::string line;
      std::getline(mol2File, line);
      foundCell = (line.find("@<TRIPOS>CRYSIN") != line.npos);
      if(foundCell) break;
    }
    while(!mol2File.eof());
  }

  if(foundCell)
  {
    // Read unit cell information

    Real r_a, r_b, r_c, r_alpha, r_beta, r_gamma; // Crystallographic constants
    size_t r_space, r_setting;                    // Space group and setting (not used)
    mol2File >> r_a >> r_b >> r_c
             >> r_alpha >> r_beta >> r_gamma
             >> r_space >> r_setting;
    if(!mol2File)
    {
      std::cerr << "Error reading file " << mol2File.fileName()
                << ": invalid CRYSIN record." << std::endl;
      geometry.clear();
      return mol2File;
    }
//    geometry.setLattice(Lattice(r_a, r_b, r_c, r_alpha, r_beta, r_gamma, system.lattice().type()));
  }

  if(!foundBonds)
  {
    // No bond record, find bonds automatically

    geometry.calculateBonds();
  }

  return mol2File;
}

MOL2File &operator<<(MOL2File &mol2File, System<Geometry> const &system)
{
  // Write system to file

  assert(mol2File.is_open() && mol2File.mode() == OUT);
  mol2File.setf(std::ios::fixed, std::ios::floatfield);
  mol2File.setf(std::ios::right, std::ios::adjustfield);
  size_t nAtoms = 0;
  size_t nBonds = 0;
  size_t const nMolecules = system.size();
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each geometry
  std::vector<size_t> firstBonds(1, 0); // Global index of the first bond in each geometry
                                        // (both are needed to convert local to global indices)

  // Calculate total number of atoms, bonds, and global indices
  for(size_t i = 0; i < system.size(); i++)
  {
    nAtoms += system[i].numAtoms();
    nBonds += system[i].numBonds();
    firstAtoms.push_back(nAtoms);
    firstBonds.push_back(nBonds);
  }
  size_t const intWidth = (nAtoms > nBonds)?(2 + (size_t)log10((Real)nAtoms)):(2 + (size_t)log10((Real)nBonds));

  // Molecule record

  mol2File << "@<TRIPOS>MOLECULE" << std::endl;
  mol2File << system.name() << std::endl;
  mol2File << nAtoms << " " << nBonds << std::endl;
  mol2File << "SMALL" << std::endl;
  mol2File << "NO_CHARGES" << std::endl;

  mol2File << std::endl << std::endl;

  // Atom records

  std::streamsize precision = mol2File.precision();
  mol2File.precision(WRITE_PRECISION);

  mol2File << "@<TRIPOS>ATOM" << std::endl;

  for(size_t i = 0; i < nMolecules; i++)
  for(size_t j = 0; j < system[i].numAtoms(); j++)
  {
    Atom atom(system[i].atom(j));
    mol2File.setf(std::ios::right, std::ios::adjustfield);
    mol2File << std::setw(intWidth) << firstAtoms[i] + j + 1 << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield); 
    mol2File << std::setw(MAX_SYMBOL_LENGTH) 
             << atom.symbol << " ";
    mol2File.setf(std::ios::right, std::ios::adjustfield); 
    for(size_t j = 0; j < 3; j++) 
      mol2File << std::setw(WRITE_PRECISION + 6) 
               << atom.position[j] << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield);
    mol2File << std::setw(MAX_TYPE_LENGTH) 
             << atom.type << " ";
    if(atom.residueName.size() > 0)
    {
      mol2File.setf(std::ios::right, std::ios::adjustfield);
      mol2File << std::setw(intWidth) << atom.residueID << " ";
      mol2File.setf(std::ios::left, std::ios::adjustfield);
      mol2File << std::setw(MAX_RESNAME_LENGTH)
               << atom.residueName << " ";
    }
    mol2File << std::endl;
  }

  // Bond records

  mol2File << "@<TRIPOS>BOND" << std::endl;

  for(size_t i = 0; i < nMolecules; i++)
  for(size_t j = 0; j < system[i].numBonds(); j++)
  {
    Bond bond(system[i].bond(j));
    mol2File.setf(std::ios::right, std::ios::adjustfield);
    mol2File << std::setw(intWidth) 
            << firstBonds[i] + j + 1 << " " 
            << firstAtoms[i] + bond.first + 1 << " " 
            << firstAtoms[i] + bond.second + 1 << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield);
    mol2File << bond.type << std::endl;
  }

  //mol2File << std::endl;

  // Unit cell records (only for 3D unit cells)

  if(system.lattice().numDimensions() == 3)
  {
    mol2File << "@<TRIPOS>CRYSIN" << std::endl;
    mol2File.precision(WRITE_LOW_PRECISION);
    mol2File << system.lattice().a() << " "
             << system.lattice().b() << " "
             << system.lattice().c() << " "
             << system.lattice().alpha() << " "
             << system.lattice().beta() << " "
             << system.lattice().gamma() << " "
             << "1 1" << std::endl;
  }
  mol2File.precision(precision); // Return to original value

  return mol2File;
}

MOL2File &operator>>(MOL2File &mol2File, System<Geometry> &system)
{
  // Read system from file

  assert(mol2File.is_open() && mol2File.mode() == IN);
  system.clear();
  
  // Search for molecule record

  do
  {
    std::string line;
    std::getline(mol2File, line);
    if(line.find("@<TRIPOS>MOLECULE") != line.npos) break;
  }
  while(!mol2File.eof());

  if(mol2File.eof())
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": MOLECULE record not found." << std::endl;
    return mol2File;
  }

  // Read molecule record

  size_t r_nAtom;             // Total number of atoms read from input file
  size_t r_nBond;             // Total number of bonds read from input file
  char r_label[LINE_LENGTH];  // Label read from input file

  mol2File.getline(r_label, LINE_LENGTH); // Read geometry label
  system.setName(std::string(r_label));

  mol2File >> r_nAtom >> r_nBond;     // Read number of atoms and number of bonds
  mol2File.ignore(LINE_LENGTH, '\n'); // Skip to next line

  if(!mol2File)
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": invalid MOLECULE record." << std::endl;
    system.clear();
    return mol2File;
  }

  // Search for atom record

  do
  {
    std::string line;
    std::getline(mol2File, line);
    if(line.find("@<TRIPOS>ATOM") != line.npos) break;
  }
  while(!mol2File.eof());

  if(mol2File.eof())
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": ATOM record not found." << std::endl;
    system.clear();
    return mol2File;
  }

  // Read atom information

  size_t r_atomID;                      // Atom id number read from input file
  std::string r_symbol;                 // Atomic symbol read from input file
  Vector3D r_position;                  // Position read from input file
  std::string r_atomType;               // Atomic type read from input file
  size_t r_residueID;                   // Residue ID number read from input file
  std::string r_residueName;            // Residue name read from input file
  std::string r_residueInfo;            // Residue info read from input file
  bool thereIsResidueInfo;              // Whether residue info is present or not
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each geometry
                                        // (needed to convert local to global indices)
  Geometry currentGeometry;             // Current geometry
  size_t nAtoms = 0;                    // Total number of atoms read
  size_t previousID = 0;                // Previous residue ID read

  while(!mol2File.eof() && nAtoms < r_nAtom)
  {
    mol2File >> r_atomID >> r_symbol >> r_position.x >> r_position.y >> r_position.z >> r_atomType;
    getline(mol2File, r_residueInfo, '\n'); // Check for residue info
    thereIsResidueInfo = (r_residueInfo.find_first_not_of(" \t\n") != r_residueInfo.npos);
    if(thereIsResidueInfo)  // Read residue info
    {
      std::istringstream resInfo(r_residueInfo);
      resInfo >> r_residueID >> r_residueName;
    }
    if(!mol2File)
    {
      std::cerr << "Error reading file " << mol2File.fileName() 
                << ": invalid ATOM record." << std::endl
                << "Atoms read: " << nAtoms << std::endl;
      system.clear();
      return mol2File;
    }
    if(nAtoms > 0 && r_residueID != previousID)
    {
      system.add(currentGeometry);
      firstAtoms.push_back(nAtoms);
      currentGeometry.clear();
    }
    nAtoms++;
    if(thereIsResidueInfo)
    {
      currentGeometry.addAtom(Atom(r_symbol, r_atomType, r_position, r_residueID, r_residueName));
      if(!currentGeometry.name().size()) currentGeometry.setName(r_residueName);
      previousID = r_residueID;
    }
    else
      currentGeometry.addAtom(Atom(r_symbol, r_atomType, r_position));
  }
  system.add(currentGeometry);
  firstAtoms.push_back(nAtoms);

  if(nAtoms < r_nAtom)
  {
    std::cerr << "Error reading file " << mol2File.fileName() << std::endl
              << "Expected " << r_nAtom << " atoms but found only " 
              << nAtoms << std::endl;
    system.clear();
    return mol2File;
  }

  // Search for bond and unit cell records

  bool foundBonds = false;
  bool foundCell = false;
  do
  {
    std::string line;
    std::getline(mol2File, line);
    foundBonds = (line.find("@<TRIPOS>BOND") != line.npos);
    foundCell = (line.find("@<TRIPOS>CRYSIN") != line.npos);
    if(foundBonds || foundCell) break;
  }
  while(!mol2File.eof());

  if(foundBonds)
  {
    // Read bond information

    size_t r_bondID;        // Bond id number read from input file
    size_t r_first;         // First atom read from input file
    size_t r_second;        // Second atom read from input file
    std::string r_bondType; // Bond type read from input file
    size_t nBonds = 0;      // Total number of bonds read

    while(!mol2File.eof() && nBonds < r_nBond)
    {
      mol2File >> r_bondID >> r_first >> r_second >> r_bondType;
      if(r_first > r_nAtom || r_second > r_nAtom)
      {
        std::cerr << "Error reading file " << mol2File.fileName() 
                  << ": invalid bond id number in BOND record" << std::endl;
        system.clear();
        return mol2File;
      }
      if(r_first > r_second) {  // Make first < second
        size_t temp = r_first;
        r_first = r_second;
        r_second = temp;
      }
      mol2File.ignore(LINE_LENGTH, '\n'); // Skip to next line
      if(!mol2File)
      {
        std::cerr << "Error reading file " << mol2File.fileName() 
                  << ": invalid BOND record." << std::endl
                  << "Bonds read: " << nBonds << std::endl;
        system.clear();
        return mol2File;
      }
      // Find which geometry the bond belongs to
      size_t geometryIndex = 0;
      while(r_first > firstAtoms[geometryIndex + 1])
        geometryIndex++;
      if(r_second <= firstAtoms[geometryIndex + 1]) // Must be intra-geometry
        system[geometryIndex].addBond(Bond(r_first - 1 - firstAtoms[geometryIndex],
                                           r_second - 1 - firstAtoms[geometryIndex],
                                           r_bondType));
      nBonds++;
    }
    if(nBonds < r_nBond)
    {
      std::cerr << "Error reading file " << mol2File.fileName() << std::endl;
      std::cerr << "Expected " << r_nBond << " bonds but found only " 
                << nBonds << std::endl;
      system.clear();
    }
    // Search for unit cell record
    do
    {
      std::string line;
      std::getline(mol2File, line);
      foundCell = (line.find("@<TRIPOS>CRYSIN") != line.npos);
      if(foundCell) break;
    }
    while(!mol2File.eof());
  }

  if(foundCell)
  {
    // Read unit cell information

    Real r_a, r_b, r_c, r_alpha, r_beta, r_gamma; // Crystallographic constants
    size_t r_space, r_setting;                    // Space group and setting (not used)
    mol2File >> r_a >> r_b >> r_c
             >> r_alpha >> r_beta >> r_gamma
             >> r_space >> r_setting;
    if(!mol2File)
    {
      std::cerr << "Error reading file " << mol2File.fileName()
                << ": invalid CRYSIN record." << std::endl;
      system.clear();
      return mol2File;
    }
    system.setLattice(Lattice(r_a, r_b, r_c, r_alpha, r_beta, r_gamma, system.lattice().type()));
     // Keep the lattice type constant - modified 01/17/09
  }

  if(!foundBonds)
  {
    // No bond record, find bonds automatically

    for(size_t i = 0; i < system.size(); i++)
      system[i].calculateBonds();
  }

  return mol2File;
}

MOL2File &operator<<(MOL2File &mol2File, System<Molecule> const &system)
{
  // Write system to file

  assert(mol2File.is_open() && mol2File.mode() == OUT);
  mol2File.setf(std::ios::fixed, std::ios::floatfield);
  mol2File.setf(std::ios::right, std::ios::adjustfield);
  size_t nAtoms = 0;
  size_t nBonds = 0;
  size_t const nMolecules = system.size();
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each geometry
  std::vector<size_t> firstBonds(1, 0); // Global index of the first bond in each geometry
                                        // (both are needed to convert local to global indices)

  // Calculate total number of atoms, bonds, and global indices
  for(size_t i = 0; i < system.size(); i++)
  {
    nAtoms += system[i].numAtoms();
    nBonds += system[i].numBonds();
    firstAtoms.push_back(nAtoms);
    firstBonds.push_back(nBonds);
  }
  size_t const intWidth = (nAtoms > nBonds)?(2 + (size_t)log10((Real)nAtoms)):(2 + (size_t)log10((Real)nBonds));

  // Molecule record

  mol2File << "@<TRIPOS>MOLECULE" << std::endl;
  mol2File << system.name() << std::endl;
  mol2File << nAtoms << " " << nBonds << std::endl;
  mol2File << "SMALL" << std::endl;
  mol2File << "NO_CHARGES" << std::endl;

  mol2File << std::endl << std::endl;

  // Atom records

  std::streamsize precision = mol2File.precision();
  mol2File.precision(WRITE_PRECISION);

  mol2File << "@<TRIPOS>ATOM" << std::endl;

  for(size_t i = 0; i < nMolecules; i++)
  for(size_t j = 0; j < system[i].numAtoms(); j++)
  {
    Atom atom(system[i].atom(j));
    mol2File.setf(std::ios::right, std::ios::adjustfield);
    mol2File << std::setw(intWidth) << firstAtoms[i] + j + 1 << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield); 
    mol2File << std::setw(MAX_SYMBOL_LENGTH) 
             << atom.symbol << " ";
    mol2File.setf(std::ios::right, std::ios::adjustfield); 
    for(size_t j = 0; j < 3; j++) 
      mol2File << std::setw(WRITE_PRECISION + 6) 
               << atom.position[j] << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield);
    mol2File << std::setw(MAX_TYPE_LENGTH) 
             << atom.type << " ";
    if(atom.residueName.size() > 0)
    {
      mol2File.setf(std::ios::right, std::ios::adjustfield);
      mol2File << std::setw(intWidth) << atom.residueID << " ";
      mol2File.setf(std::ios::left, std::ios::adjustfield);
      mol2File << std::setw(MAX_RESNAME_LENGTH)
               << atom.residueName << " ";
    }
    mol2File << std::endl;
  }

  // Bond records

  mol2File << "@<TRIPOS>BOND" << std::endl;

  for(size_t i = 0; i < nMolecules; i++)
  for(size_t j = 0; j < system[i].numBonds(); j++)
  {
    Bond bond(system[i].bond(j));
    mol2File.setf(std::ios::right, std::ios::adjustfield);
    mol2File << std::setw(intWidth) 
            << firstBonds[i] + j + 1 << " " 
            << firstAtoms[i] + bond.first + 1 << " " 
            << firstAtoms[i] + bond.second + 1 << " ";
    mol2File.setf(std::ios::left, std::ios::adjustfield);
    mol2File << bond.type << std::endl;
  }

  //mol2File << std::endl;

  // Unit cell records (only for 3D unit cells)

  if(system.lattice().numDimensions() == 3)
  {
    mol2File << "@<TRIPOS>CRYSIN" << std::endl;
    mol2File.precision(WRITE_LOW_PRECISION);
    mol2File << system.lattice().a() << " "
             << system.lattice().b() << " "
             << system.lattice().c() << " "
             << system.lattice().alpha() << " "
             << system.lattice().beta() << " "
             << system.lattice().gamma() << " "
             << "1 1" << std::endl;
  }
  mol2File.precision(precision); // Return to original value

  return mol2File;
}

MOL2File &operator>>(MOL2File &mol2File, System<Molecule> &system)
{
  // Read system from file

  assert(mol2File.is_open() && mol2File.mode() == IN);
  system.clear();
  
  // Search for molecule record

  do
  {
    std::string line;
    std::getline(mol2File, line);
    if(line.find("@<TRIPOS>MOLECULE") != line.npos) break;
  }
  while(!mol2File.eof());

  if(mol2File.eof())
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": MOLECULE record not found." << std::endl;
    return mol2File;
  }

  // Read molecule record

  size_t r_nAtom;             // Total number of atoms read from input file
  size_t r_nBond;             // Total number of bonds read from input file
  char r_label[LINE_LENGTH];  // Label read from input file

  mol2File.getline(r_label, LINE_LENGTH); // Read geometry label
  system.setName(std::string(r_label));

  mol2File >> r_nAtom >> r_nBond;     // Read number of atoms and number of bonds
  mol2File.ignore(LINE_LENGTH, '\n'); // Skip to next line

  if(!mol2File)
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": invalid MOLECULE record." << std::endl;
    system.clear();
    return mol2File;
  }

  // Search for atom record

  do
  {
    std::string line;
    std::getline(mol2File, line);
    if(line.find("@<TRIPOS>ATOM") != line.npos) break;
  }
  while(!mol2File.eof());

  if(mol2File.eof())
  {
    std::cerr << "Error reading file " << mol2File.fileName() << ": ATOM record not found." << std::endl;
    system.clear();
    return mol2File;
  }

  // Read atom information

  size_t r_atomID;                      // Atom id number read from input file
  std::string r_symbol;                 // Atomic symbol read from input file
  Vector3D r_position;                  // Position read from input file
  std::string r_atomType;               // Atomic type read from input file
  size_t r_residueID;                   // Residue ID number read from input file
  std::string r_residueName;            // Residue name read from input file
  std::string r_residueInfo;            // Residue info read from input file
  bool thereIsResidueInfo;              // Whether residue info is present or not
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each geometry
                                        // (needed to convert local to global indices)
  Molecule currentMolecule;             // Current molecule
  size_t nAtoms = 0;                    // Total number of atoms read
  size_t previousID = 0;                // Previous residue ID read

  while(!mol2File.eof() && nAtoms < r_nAtom)
  {
    mol2File >> r_atomID >> r_symbol >> r_position.x >> r_position.y >> r_position.z >> r_atomType;
    getline(mol2File, r_residueInfo, '\n'); // Check for residue info
    thereIsResidueInfo = (r_residueInfo.find_first_not_of(" \t\n") != r_residueInfo.npos);
    if(thereIsResidueInfo)  // Read residue info
    {
      std::istringstream resInfo(r_residueInfo);
      resInfo >> r_residueID >> r_residueName;
    }
    if(!mol2File)
    {
      std::cerr << "Error reading file " << mol2File.fileName() 
                << ": invalid ATOM record." << std::endl
                << "Atoms read: " << nAtoms << std::endl;
      system.clear();
      return mol2File;
    }
    if(nAtoms > 0 && r_residueID != previousID)
    {
      system.add(currentMolecule);
      firstAtoms.push_back(nAtoms);
      currentMolecule.clear();
    }
    nAtoms++;
    if(thereIsResidueInfo)
    {
      currentMolecule.addAtom(Atom(r_symbol, r_atomType, r_position, r_residueID, r_residueName));
      if(!currentMolecule.name().size()) currentMolecule.setName(r_residueName);
      previousID = r_residueID;
    }
    else
      currentMolecule.addAtom(Atom(r_symbol, r_atomType, r_position));
  }
  system.add(currentMolecule);
  firstAtoms.push_back(nAtoms);

  if(nAtoms < r_nAtom)
  {
    std::cerr << "Error reading file " << mol2File.fileName() << std::endl
              << "Expected " << r_nAtom << " atoms but found only " 
              << nAtoms << std::endl;
    system.clear();
    return mol2File;
  }

  // Search for bond and unit cell records

  bool foundBonds = false;
  bool foundCell = false;
  do
  {
    std::string line;
    std::getline(mol2File, line);
    foundBonds = (line.find("@<TRIPOS>BOND") != line.npos);
    foundCell = (line.find("@<TRIPOS>CRYSIN") != line.npos);
    if(foundBonds || foundCell) break;
  }
  while(!mol2File.eof());

  if(foundBonds)
  {
    // Read bond information

    size_t r_bondID;        // Bond id number read from input file
    size_t r_first;         // First atom read from input file
    size_t r_second;        // Second atom read from input file
    std::string r_bondType; // Bond type read from input file
    size_t nBonds = 0;      // Total number of bonds read

    while(!mol2File.eof() && nBonds < r_nBond)
    {
      mol2File >> r_bondID >> r_first >> r_second >> r_bondType;
      if(r_first > r_nAtom || r_second > r_nAtom)
      {
        std::cerr << "Error reading file " << mol2File.fileName() 
                  << ": invalid bond id number in BOND record" << std::endl;
        system.clear();
        return mol2File;
      }
      if(r_first > r_second) {  // Make first < second
        size_t temp = r_first;
        r_first = r_second;
        r_second = temp;
      }
      mol2File.ignore(LINE_LENGTH, '\n'); // Skip to next line
      if(!mol2File)
      {
        std::cerr << "Error reading file " << mol2File.fileName() 
                  << ": invalid BOND record." << std::endl
                  << "Bonds read: " << nBonds << std::endl;
        system.clear();
        return mol2File;
      }
      // Find which geometry the bond belongs to
      size_t geometryIndex = 0;
      while(r_first > firstAtoms[geometryIndex + 1])
        geometryIndex++;
      if(r_second <= firstAtoms[geometryIndex + 1]) // Must be intra-geometry
        system[geometryIndex].addBond(Bond(r_first - 1 - firstAtoms[geometryIndex],
                                           r_second - 1 - firstAtoms[geometryIndex],
                                           r_bondType));
      nBonds++;
    }
    if(nBonds < r_nBond)
    {
      std::cerr << "Error reading file " << mol2File.fileName() << std::endl;
      std::cerr << "Expected " << r_nBond << " bonds but found only " 
                << nBonds << std::endl;
      system.clear();
    }
    // Search for unit cell record
    do
    {
      std::string line;
      std::getline(mol2File, line);
      foundCell = (line.find("@<TRIPOS>CRYSIN") != line.npos);
      if(foundCell) break;
    }
    while(!mol2File.eof());
  }

  if(foundCell)
  {
    // Read unit cell information

    Real r_a, r_b, r_c, r_alpha, r_beta, r_gamma; // Crystallographic constants
    size_t r_space, r_setting;                    // Space group and setting (not used)
    mol2File >> r_a >> r_b >> r_c
             >> r_alpha >> r_beta >> r_gamma
             >> r_space >> r_setting;
    if(!mol2File)
    {
      std::cerr << "Error reading file " << mol2File.fileName()
                << ": invalid CRYSIN record." << std::endl;
      system.clear();
      return mol2File;
    }
    system.setLattice(Lattice(r_a, r_b, r_c, r_alpha, r_beta, r_gamma, system.lattice().type()));
      // Keep the lattice type constant - modified 01/17/09
  }

  if(!foundBonds)
  {
    // No bond record, find bonds automatically

    for(size_t i = 0; i < system.size(); i++)
      system[i].calculateBonds();
  }

  return mol2File;
}

