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
** .pdb file format
*/

#include <sstream>
#include <vector>
#include <set>
#include <iomanip>
#include "common/include/assert.h"
#include "mymol/include/file_formats/fileformats.h"

// Constructors

PDBFile::PDBFile()
:
IOFile()
{}

PDBFile::PDBFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

PDBFile &operator<<(PDBFile &pdbFile, Geometry const &geometry)
{
  assert(pdbFile.is_open() && pdbFile.mode() == OUT);

  // COMPND record
  if(geometry.name().size() > 0)
    pdbFile << "COMPND    " << geometry.name() << std::endl;

  // ATOM records
  for(size_t i = 0; i < geometry.numAtoms(); ++i)
  {
    Atom const &currentAtom = geometry.atom(i);
    std::string roleName = currentAtom.roleName;
    if(roleName.size() == 0 ||
       roleName.find_first_not_of(" ") == roleName.npos)
    {
      // No role name - create from atomic symbol
      if(currentAtom.symbol.size() == 1)
        roleName = " " + currentAtom.symbol;
      else
        roleName = currentAtom.symbol;
      while(roleName.size() < 4)
        roleName = roleName + " ";
    }
    if(roleName.size() > 4) roleName = roleName.substr(0, 4);
    std::string residueName = currentAtom.residueName;
    if(residueName.size() == 0 ||
       residueName.find_first_not_of(" ") == residueName.npos)
      residueName = "RES "; // Make up residue name
    while(residueName.size() < 4)
      residueName = residueName + " ";
    if(residueName.size() > 4) residueName = residueName.substr(0, 4);

    pdbFile << "ATOM  " 
            << std::setw(5) << i + 1
            << " "
            << std::setw(4) << roleName
            << " "
            << std::setw(4) << residueName
            << " " // Chain identifier - add if needed
            << std::setw(4) << currentAtom.residueID
            << "    " << std::setiosflags(std::ios::fixed)
            << std::setw(8) << std::setprecision(3) << currentAtom.position.x
            << std::setw(8) << std::setprecision(3) << currentAtom.position.y
            << std::setw(8) << std::setprecision(3) << currentAtom.position.z
            << "  1.00"               // Change this to use atom occupancies
            << std::setw(6) << std::setprecision(2) << currentAtom.tempFactor
            << "      "
            << std::setw(4) << residueName
            << "          "
//            << "              "       // Could add here symbol, etc.
            << std::endl; 
  }
  
  // CONECT records
  // Make list of atoms bonded to each atom
  std::vector<std::set<size_t> > conectList(geometry.numAtoms());
  for(size_t i = 0; i < geometry.numBonds(); ++i)
  {
    Bond const &currentBond = geometry.bond(i);
    conectList[currentBond.first].insert(currentBond.second);
    conectList[currentBond.second].insert(currentBond.first);
  }
  // Write CONECT records
  for(size_t i = 0; i < geometry.numAtoms(); ++i)
  {
    if(conectList[i].size() > 0)
    {
      if(conectList[i].size() < 5)
      {
        pdbFile << "CONECT"
                << std::setw(5) << i + 1;
        for(std::set<size_t>::iterator it = conectList[i].begin(); 
            it != conectList[i].end(); ++it)
          pdbFile << std::setw(5) << *it + 1;
      }
      else if(conectList[i].size() < 9)
      {
        pdbFile << "CONECT"
                << std::setw(5) << i + 1;
        std::set<size_t>::iterator it = conectList[i].begin();
        for(size_t iBond = 0; iBond < 4; ++it, ++ iBond)
          pdbFile << std::setw(5) << *it + 1;
        for(; it != conectList[i].end(); ++it)
          pdbFile << std::setw(5) << *it + 1;
      }
      else
      {
        std::cerr << "Too many bonds when writing pdb file " << pdbFile.fileName() << std::endl;
        return pdbFile;
      }
      pdbFile << std::endl;
    }
  }
  pdbFile << "END" << std::endl;

  return pdbFile;
}

PDBFile &operator>>(PDBFile &pdbFile, Geometry &geometry)
{
  assert(pdbFile.is_open() && pdbFile.mode() == IN);
  bool newGeometry = (geometry.numAtoms() == 0);

  // Search for COMPND or ATOM/HETATOM record

  std::string line;     // Line read from file
  do
  {
    std::getline(pdbFile, line);
    if(line.find("COMPND") != line.npos ||
       line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
  }
  while(!pdbFile.eof());

  if(line.find("COMPND") != line.npos) // COMPND record found: store name
    geometry.setName(line.substr(10));

  // Search for ATOM/HETATM records
  while(!pdbFile.eof())
  {
    if(line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
    std::getline(pdbFile, line);
  }

  if(pdbFile.eof())
  {
    std::cerr << "Error reading pdb file " << pdbFile.fileName()
              << ": no ATOM/HETATM records" << std::endl;
  }

  // Read ATOM/HETATM records

  size_t iAtom = 0; // Index of current atom

  while(line.find("ATOM") != line.npos ||
        line.find("HETATM") != line.npos ||
        line.find("TER") != line.npos)
  {
    if(line.find("TER") != line.npos) // Ignore TER records
    {
      std::getline(pdbFile, line);
      continue;
    }
 
    size_t r_serial;        // Atom serial number read from file
    std::string r_roleName; // Atom role name read from file
    std::string r_symbol;   // Symbol read from file
    Real r_x, r_y, r_z;     // Coordinates read from file
    std::string r_resName;  // Residue name read from file
    size_t r_resID;         // Residue ID number read from file
    Real r_tempFactor;      // Temperature factor read from file
    std::istringstream serialStream(line.substr(6, 5));
    serialStream >> r_serial;
    if(serialStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(r_serial != iAtom + 1)
    {
      std::cerr << "Error reading pdb file " << pdbFile.fileName()
                << ": atom serial numbers not correlated" << std::endl;
      return pdbFile;
    }
    r_roleName = line.substr(12, 4);
    // Find the atomic symbol
    size_t symbolIndex = r_roleName.find_first_not_of(" 1234567890");
    if(symbolIndex == r_roleName.npos)
    {
      std::cerr << "Error reading ATOM/HETATOM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    r_symbol = r_roleName.at(symbolIndex); // Note that this still works only for one-letter symbols
    /*
    // Old version below
    r_symbol = (line.substr(12, 1).find_first_of(" 1234567890") != line.substr(12, 1).npos)?
                ((line.substr(13, 1) == " ")?line.substr(14, 1):line.substr(13, 1)):
                line.substr(12, 2);
    if(r_symbol == " ")
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    */
    std::istringstream coordStream(line.substr(30, 24));
    coordStream >> r_x >> r_y >> r_z;
    if(coordStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    r_resName = line.substr(17, 4);
    std::istringstream resIDStream(line.substr(22, 5));
    resIDStream >> r_resID;
    if(resIDStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    std::istringstream tempFactorStream(line.substr(60, 6));
    tempFactorStream >> r_tempFactor;
    if(tempFactorStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(newGeometry)
      geometry.addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z), r_resID, 
                            r_resName, r_roleName, r_tempFactor));
    else
    {
      if(r_symbol.find(geometry.atom(iAtom).symbol) == r_symbol.npos ||
         r_resName.find(geometry.atom(iAtom).residueName) == r_resName.npos ||
         r_resID != geometry.atom(iAtom).residueID ||
         r_roleName != geometry.atom(iAtom).roleName)
      {
        std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
        return pdbFile;
      }
      geometry.setPosition(iAtom, Vector3D(r_x, r_y, r_z));
      geometry.setTemperatureFactor(iAtom, r_tempFactor);
    }
    ++iAtom;
    std::getline(pdbFile, line);
    if(pdbFile.eof()) break;
  } // End of ATOM/HETATM records

  if(!newGeometry) return pdbFile;  // Connectivity read only for new geometry

  // Search for CONECT records
  while(!pdbFile.eof())
  {
    if(line.find("CONECT") != line.npos) break;
    if(line.find("END") != line.npos) break; // In case there are no CONECT records
    std::getline(pdbFile, line);
  }
  if(pdbFile.eof() || (line.find("END") != line.npos))
  {
    // No connect records - find bonds automatically if this is a new geometry
    if(newGeometry) geometry.calculateBonds();
    return pdbFile;
  }

  // Read CONECT records
  while(line.find("CONECT") != line.npos)
  {
    size_t r_atom1, r_atom2; // Atoms forming a bond read from file
    std::istringstream atom1Stream(line.substr(6, 5));
    atom1Stream >> r_atom1;
    if(atom1Stream.fail())
    {
      std::cerr << "Error reading CONECT record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    for(size_t i = 11; i < 27; i += 5)
    {
      std::istringstream atom2Stream(line.substr(i, 5));
      atom2Stream >> r_atom2;
      if(atom2Stream.fail()) break;
      if(r_atom2 > r_atom1) geometry.addBond(Bond(r_atom1 - 1, r_atom2 - 1)); // CONECTs are redundant
    }
    std::getline(pdbFile, line);
    if(pdbFile.eof()) break;
  }
  return pdbFile;
}

PDBFile &operator<<(PDBFile &pdbFile, System<Geometry> const &system)
{
  assert(pdbFile.is_open() && pdbFile.mode() == OUT);

  // Sanity check
  if(system.size() == 0 || system[0].numAtoms() == 0)
  {
    std::cerr << "Error while writing PDB file " << pdbFile.fileName()
              << ": empty system" << std::endl;
    return pdbFile;
  }
  // Find out if we need to assign residue IDs
  bool needToAssignResidueIDs = false;
  for(size_t i = 1; i < system.size() && !needToAssignResidueIDs; ++i)
  {
    // Sanity check (again)
    if(system[i].numAtoms() == 0)
    {
      std::cerr << "Error while writing PDB file " << pdbFile.fileName()
                << ": empty system" << std::endl;
      return pdbFile;
    }
    if(system[i].atom(0).residueID == system[i - 1].atom(0).residueID)
      needToAssignResidueIDs = true;
  }

  // COMPND record
  if(system.name().size() > 0)
    pdbFile << "COMPND    " << system.name() << std::endl;

  // CRYST1 record
  Lattice const &lattice = system.lattice();
  if(lattice.numDimensions() > 0)
    pdbFile << "CRYST1"
            << std::setiosflags(std::ios::fixed)
            << std::setw(9) << std::setprecision(3) << lattice.a()
            << std::setw(9) << std::setprecision(3) << lattice.b()
            << std::setw(9) << std::setprecision(3) << lattice.c()
            << std::setw(7) << std::setprecision(2) << lattice.alpha()
            << std::setw(7) << std::setprecision(2) << lattice.beta()
            << std::setw(7) << std::setprecision(2) << lattice.gamma()
            << " P 1           1" << std::endl;

  // ATOM records
  size_t iAtom = 0;               // Current global atom index
  std::vector<size_t> firstAtoms; // Global atom index of first atom in each residue
  for(size_t i = 0; i < system.size(); ++i)
  {
    firstAtoms.push_back(iAtom);
    Geometry const &currentResidue = system[i];
    for(size_t j = 0; j < currentResidue.numAtoms(); ++j)
    {
      Atom const &currentAtom = currentResidue.atom(j);
      std::string roleName = currentAtom.roleName;
      if(roleName.size() == 0 ||
        roleName.find_first_not_of(" ") == roleName.npos)
      {
        // No role name - create from atomic symbol
        if(currentAtom.symbol.size() == 1)
          roleName = " " + currentAtom.symbol;
        else
          roleName = currentAtom.symbol;
        while(roleName.size() < 4)
          roleName = roleName + " ";
      }
      if(roleName.size() > 4) roleName = roleName.substr(0, 4);
      std::string residueName = currentAtom.residueName;
      if(residueName.size() == 0 ||
        residueName.find_first_not_of(" ") == residueName.npos)
      {
        if(system[i].name().size() > 0 &&
           system[i].name().find_first_not_of(" ") != residueName.npos)
          residueName = system[i].name();
        else residueName = "RES "; // Make up residue name
      }
      while(residueName.size() < 4)
        residueName = residueName + " ";
      if(residueName.size() > 4) residueName = residueName.substr(0, 4);    

      pdbFile << "ATOM  " 
              << std::setw(5) << iAtom + 1
              << " "
              << std::setw(4) << roleName
              << " "
              << std::setw(4) << residueName
              << " " // Chain identifier - add if needed
              << std::setw(4) << (needToAssignResidueIDs?(i + 1):currentAtom.residueID)
              << "    " << std::setiosflags(std::ios::fixed)
              << std::setw(8) << std::setprecision(3) << currentAtom.position.x
              << std::setw(8) << std::setprecision(3) << currentAtom.position.y
              << std::setw(8) << std::setprecision(3) << currentAtom.position.z
              << "  1.00"               // Change this to use atom occupancies
              << std::setw(6) << std::setprecision(2) << currentAtom.tempFactor
              << "      "
              << std::setw(4) << residueName
              << "          "
//            << "              "       // Could add here symbol, etc.
              << std::endl; 
      ++iAtom;
    } // End of loop over atoms
  } // End of loop over residues

  // CONECT records

  // Make list of atoms bonded to each atom
  std::vector<std::set<size_t> > conectList(iAtom);
  for(size_t i = 0; i < system.size(); ++i)
  {
    Geometry const &currentResidue = system[i];
    for(size_t j = 0; j < currentResidue.numBonds(); ++j)
    {
      Bond const &currentBond = currentResidue.bond(j);
      conectList[firstAtoms[i] + currentBond.first].insert(firstAtoms[i] + currentBond.second);
      conectList[firstAtoms[i] + currentBond.second].insert(firstAtoms[i] + currentBond.first);
    }
  } // End of loop over residues

  // Write CONECT records
  for(size_t i = 0; i < conectList.size(); ++i)
  {
    if(conectList[i].size() > 0)
    {
      if(conectList[i].size() < 5)
      {
        pdbFile << "CONECT"
                << std::setw(5) << i + 1;
        for(std::set<size_t>::iterator it = conectList[i].begin(); 
            it != conectList[i].end(); ++it)
          pdbFile << std::setw(5) << *it + 1;
      }
      else if(conectList[i].size() < 9)
      {
        pdbFile << "CONECT"
                << std::setw(5) << i + 1;
        std::set<size_t>::iterator it = conectList[i].begin();
        for(size_t iBond = 0; iBond < 4; ++it, ++ iBond)
          pdbFile << std::setw(5) << *it + 1;
        for(; it != conectList[i].end(); ++it)
          pdbFile << std::setw(5) << *it + 1;
      }
      else
      {
        std::cerr << "Too many bonds when writing pdb file " << pdbFile.fileName() << std::endl;
        return pdbFile;
      }
      pdbFile << std::endl;
    }
  } // End of loop over bonds

  pdbFile << "END" << std::endl;

  return pdbFile;
}

PDBFile &operator>>(PDBFile &pdbFile, System<Geometry> &system)
{
  assert(pdbFile.is_open() && pdbFile.mode() == IN);
  bool newSystem = (system.size() == 0);

  // Search for COMPND, CRYST1 or ATOM/HETATOM records

  std::string line; // Line read from file

  do
  {
    std::getline(pdbFile, line);
    if(line.find("COMPND") != line.npos ||
       line.find("CRYST1") != line.npos ||
       line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
  }
  while(!pdbFile.eof());

  if(line.find("COMPND") != line.npos) // COMPND record found: store name
    system.setName(line.substr(10));

  // Search for CRYST1 or ATOM/HETATM records
  while(!pdbFile.eof())
  {
    if(line.find("CRYST1") != line.npos ||
       line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
    std::getline(pdbFile, line);
  }

  if(line.find("CRYST1") != line.npos)
  {
    // CRYST1 record found: read unit cell information
    std::istringstream crystLineStream(line);
    std::string dummyString;                      // CRYST1 string
    Real r_a, r_b, r_c, r_alpha, r_beta, r_gamma; // Cell parameters read

    crystLineStream >> dummyString >> r_a >> r_b >> r_c >> r_alpha >> r_beta >> r_gamma;
    if(crystLineStream.fail())
    {  
      std::cerr << "Error reading CRYST1 record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    system.setLattice(Lattice(r_a, r_b, r_c, r_alpha, r_beta, r_gamma));
  }

  // Search for ATOM/HETATM records
  while(!pdbFile.eof())
  {
    if(line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
    std::getline(pdbFile, line);
  }
  if(pdbFile.eof())
  {
    std::cerr << "Error reading pdb file " << pdbFile.fileName()
              << ": no ATOM/HETATM records" << std::endl;
  }

  // Read ATOM/HETATM records

  Geometry currentResidue;              // Current geometry
  size_t iAtom = 0;                     // Index of current atom
  size_t iRes = 0;                      // Index of current residue
  size_t previousID = 0;                // Previous residue ID read
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each geometry
                                        // (needed to convert local to global indices)

  while(line.find("ATOM") != line.npos ||
        line.find("HETATM") != line.npos ||
        line.find("TER") != line.npos)
  {
    if(line.find("TER") != line.npos) // Ignore TER records
    {
      std::getline(pdbFile, line);
      continue;
    }
 
    size_t r_serial;        // Atom serial number read from file
    std::string r_roleName; // Atom role name read from file
    std::string r_symbol;   // Symbol read from file
    Real r_x, r_y, r_z;     // Coordinates read from file
    std::string r_resName;  // Residue name read from file
    size_t r_resID;         // Residue ID number read from file
    Real r_tempFactor;      // Temperature factor read from file
    std::istringstream serialStream(line.substr(6, 5));
    serialStream >> r_serial;
    if(serialStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(r_serial != iAtom + 1)
    {
      std::cerr << "Error reading pdb file " << pdbFile.fileName()
                << ": atom serial numbers not correlated" << std::endl;
      return pdbFile;
    }
    r_roleName = line.substr(12, 4);
    // Find the atomic symbol
    size_t symbolIndex = r_roleName.find_first_not_of(" 12345667890");
    if(symbolIndex == r_roleName.npos)
    {
      std::cerr << "Error reading ATOM/HETATOM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    r_symbol = r_roleName.at(symbolIndex); // Note that this still works only for one-letter symbols
    /*
    // Old version below
    r_symbol = (line.substr(12, 1).find_first_of(" 1234567890") != line.substr(12, 1).npos)?
                ((line.substr(13, 1) == " ")?line.substr(14, 1):line.substr(13, 1)):
                line.substr(12, 2);
    if(r_symbol == " ")
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    */
    std::istringstream coordStream(line.substr(30, 24));
    coordStream >> r_x >> r_y >> r_z;
    if(coordStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    std::istringstream tempFactorStream(line.substr(60, 6));
    tempFactorStream >> r_tempFactor;
    if(tempFactorStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    r_resName = line.substr(17, 4);
    if(!currentResidue.name().size()) currentResidue.setName(r_resName);  // Use same name for residue
    std::istringstream resIDStream(line.substr(22, 5));
    resIDStream >> r_resID;
    if(resIDStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(iAtom == 0) previousID = r_resID;  // Avoid trouble when first residue ID is not zero
    if(iAtom > 0 && r_resID != previousID)
    {
      // Finished reading current residue
      if(newSystem)
        system.add(currentResidue);
      else
      {
        if(iRes >= system.size())
        {
          std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
          return pdbFile;
        }
        if(currentResidue.numAtoms() != system[iRes].numAtoms())
        {
          std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
          return pdbFile;
        }

        for(size_t i = 0; i < currentResidue.numAtoms(); ++i)
        {
          if(currentResidue.atom(i).symbol != system[iRes].atom(i).symbol ||
             currentResidue.atom(i).residueName != system[iRes].atom(i).residueName ||
             currentResidue.atom(i).residueID != system[iRes].atom(i).residueID ||
             currentResidue.atom(i).roleName != system[iRes].atom(i).roleName)
          {
            std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
            return pdbFile;
          }
          system[iRes].setPosition(i, currentResidue.atom(i).position);
          system[iRes].setTemperatureFactor(i, currentResidue.atom(i).tempFactor);
        }
      }
      firstAtoms.push_back(iAtom);
      currentResidue.clear();
      previousID = r_resID;
      ++iRes;
    }
    currentResidue.addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z), r_resID, 
                                r_resName, r_roleName, r_tempFactor));
    ++iAtom;
    std::getline(pdbFile, line);
    if(pdbFile.eof()) break;
  } // End of ATOM/HETATM sections

  if(newSystem)
    system.add(currentResidue);
  else
  {
    if(iRes >= system.size())
    {
      std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(currentResidue.numAtoms() != system[iRes].numAtoms())
    {
      std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }

    for(size_t i = 0; i < currentResidue.numAtoms(); ++i)
    {
      if(currentResidue.atom(i).symbol != system[iRes].atom(i).symbol ||
          currentResidue.atom(i).residueName != system[iRes].atom(i).residueName ||
          currentResidue.atom(i).residueID != system[iRes].atom(i).residueID ||
          currentResidue.atom(i).roleName != system[iRes].atom(i).roleName)
      {
        std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
        return pdbFile;
      }
      system[iRes].setPosition(i, currentResidue.atom(i).position);
      system[iRes].setTemperatureFactor(i, currentResidue.atom(i).tempFactor);
    }
  }
  firstAtoms.push_back(iAtom);

  if(!newSystem) return pdbFile;  // Connectivity read only for new system

  // Search for CONECT records
  while(!pdbFile.eof())
  {
    if(line.find("CONECT") != line.npos) break;
    if(line.find("END") != line.npos) break; // In case there are no CONECT records
    std::getline(pdbFile, line);
  }
  if(pdbFile.eof() || (line.find("END") != line.npos))
  {
    // No connect records - find bonds automatically if this is a new system
    if(newSystem)
      for(size_t i = 0; i < system.size(); ++i)
        system[i].calculateBonds();
    return pdbFile;
  }

  // Read CONECT records
  while(line.find("CONECT") != line.npos)
  {
    size_t r_atom1, r_atom2; // Atoms forming a bond read from file
    std::istringstream atom1Stream(line.substr(6, 5));
    atom1Stream >> r_atom1;
    if(atom1Stream.fail())
    {
      std::cerr << "Error reading CONECT record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    for(size_t i = 11; i < 27; i += 5)
    {
      std::istringstream atom2Stream(line.substr(i, 5));
      atom2Stream >> r_atom2;
      if(atom2Stream.fail()) break;
      if(r_atom2 > r_atom1) // CONECTs are redundant
      {
        // Find which geometry the bond belongs to
        size_t geometryIndex = 0;
        while(r_atom1 > firstAtoms[geometryIndex+1])
          ++geometryIndex;
        if(r_atom2 <= firstAtoms[geometryIndex + 1]) // Must be intra-geometry
          system[geometryIndex].addBond(Bond(r_atom1 - 1 - firstAtoms[geometryIndex],
                                             r_atom2 - 1 - firstAtoms[geometryIndex]));
      }
    }
    std::getline(pdbFile, line);
    if(pdbFile.eof()) break;
  }
  return pdbFile;
}

PDBFile &operator<<(PDBFile &pdbFile, System<Molecule> const &system)
{
  assert(pdbFile.is_open() && pdbFile.mode() == OUT);

  // Sanity check
  if(system.size() == 0 || system[0].numAtoms() == 0)
  {
    std::cerr << "Error while writing PDB file " << pdbFile.fileName()
              << ": empty system" << std::endl;
    return pdbFile;
  }
  // Find out if we need to assign residue IDs
  bool needToAssignResidueIDs = false;
  for(size_t i = 1; i < system.size() && !needToAssignResidueIDs; ++i)
  {
    // Sanity check (again)
    if(system[i].numAtoms() == 0)
    {
      std::cerr << "Error while writing PDB file " << pdbFile.fileName()
                << ": empty system" << std::endl;
      return pdbFile;
    }
    if(system[i].atom(0).residueID == system[i - 1].atom(0).residueID)
      needToAssignResidueIDs = true;
  }

  // COMPND record
  if(system.name().size() > 0)
    pdbFile << "COMPND    " << system.name() << std::endl;

  // CRYST1 record
  Lattice const &lattice = system.lattice();
  if(lattice.numDimensions() > 0)
    pdbFile << "CRYST1"
            << std::setiosflags(std::ios::fixed)
            << std::setw(9) << std::setprecision(3) << lattice.a()
            << std::setw(9) << std::setprecision(3) << lattice.b()
            << std::setw(9) << std::setprecision(3) << lattice.c()
            << std::setw(7) << std::setprecision(2) << lattice.alpha()
            << std::setw(7) << std::setprecision(2) << lattice.beta()
            << std::setw(7) << std::setprecision(2) << lattice.gamma()
            << " P 1           1" << std::endl;

  // ATOM records
  size_t iAtom = 0;               // Current global atom index
  std::vector<size_t> firstAtoms; // Global atom index of first atom in each residue
  for(size_t i = 0; i < system.size(); ++i)
  {
    firstAtoms.push_back(iAtom);
    Molecule const &currentResidue = system[i];
    for(size_t j = 0; j < currentResidue.numAtoms(); ++j)
    {
      Atom const &currentAtom = currentResidue.atom(j);
      std::string roleName = currentAtom.roleName;
      if(roleName.size() == 0 ||
        roleName.find_first_not_of(" ") == roleName.npos)
      {
        // No role name - create from atomic symbol
        if(currentAtom.symbol.size() == 1)
          roleName = " " + currentAtom.symbol;
        else
          roleName = currentAtom.symbol;
        while(roleName.size() < 4)
          roleName = roleName + " ";
      }
      if(roleName.size() > 4) roleName = roleName.substr(0, 4);
      std::string residueName = currentAtom.residueName;
      if(residueName.size() == 0 ||
        residueName.find_first_not_of(" ") == residueName.npos)
      {
        if(system[i].name().size() > 0 &&
           system[i].name().find_first_not_of(" ") != residueName.npos)
          residueName = system[i].name();
        else residueName = "RES "; // Make up residue name
      }
      while(residueName.size() < 4)
        residueName = residueName + " ";
      if(residueName.size() > 4) residueName = residueName.substr(0, 4);    

      pdbFile << "ATOM  " 
              << std::setw(5) << iAtom + 1
              << " "
              << std::setw(4) << roleName
              << " "
              << std::setw(4) << residueName
              << " " // Chain identifier - add if needed
              << std::setw(4) << (needToAssignResidueIDs?(i + 1):currentAtom.residueID)
              << "    " << std::setiosflags(std::ios::fixed)
              << std::setw(8) << std::setprecision(3) << currentAtom.position.x
              << std::setw(8) << std::setprecision(3) << currentAtom.position.y
              << std::setw(8) << std::setprecision(3) << currentAtom.position.z
              << "  1.00"               // Change this to use atom occupancies
              << std::setw(6) << std::setprecision(2) << currentAtom.tempFactor
              << "      "
              << std::setw(4) << residueName
              << "          "
//            << "              "       // Could add here symbol, etc.
              << std::endl; 
      ++iAtom;
    } // End of loop over atoms
  } // End of loop over residues

  // CONECT records

  // Make list of atoms bonded to each atom
  std::vector<std::set<size_t> > conectList(iAtom);
  for(size_t i = 0; i < system.size(); ++i)
  {
    Molecule const &currentResidue = system[i];
    for(size_t j = 0; j < currentResidue.numBonds(); ++j)
    {
      Bond const &currentBond = currentResidue.bond(j);
      conectList[firstAtoms[i] + currentBond.first].insert(firstAtoms[i] + currentBond.second);
      conectList[firstAtoms[i] + currentBond.second].insert(firstAtoms[i] + currentBond.first);
    }
  } // End of loop over residues

  // Write CONECT records
  for(size_t i = 0; i < conectList.size(); ++i)
  {
    if(conectList[i].size() > 0)
    {
      if(conectList[i].size() < 5)
      {
        pdbFile << "CONECT"
                << std::setw(5) << i + 1;
        for(std::set<size_t>::iterator it = conectList[i].begin(); 
            it != conectList[i].end(); ++it)
          pdbFile << std::setw(5) << *it + 1;
      }
      else if(conectList[i].size() < 9)
      {
        pdbFile << "CONECT"
                << std::setw(5) << i + 1;
        std::set<size_t>::iterator it = conectList[i].begin();
        for(size_t iBond = 0; iBond < 4; ++it, ++ iBond)
          pdbFile << std::setw(5) << *it + 1;
        for(; it != conectList[i].end(); ++it)
          pdbFile << std::setw(5) << *it + 1;
      }
      else
      {
        std::cerr << "Too many bonds when writing pdb file " << pdbFile.fileName() << std::endl;
        return pdbFile;
      }
      pdbFile << std::endl;
    }
  } // End of loop over bonds

  pdbFile << "END" << std::endl;

  return pdbFile;
}

PDBFile &operator>>(PDBFile &pdbFile, System<Molecule> &system)
{
  assert(pdbFile.is_open() && pdbFile.mode() == IN);
  bool newSystem = (system.size() == 0);

  // Search for COMPND, CRYST1 or ATOM/HETATOM records

  std::string line; // Line read from file

  do
  {
    std::getline(pdbFile, line);
    if(line.find("COMPND") != line.npos ||
       line.find("CRYST1") != line.npos ||
       line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
  }
  while(!pdbFile.eof());

  if(line.find("COMPND") != line.npos) // COMPND record found: store name
    system.setName(line.substr(10));

  // Search for CRYST1 or ATOM/HETATM records
  while(!pdbFile.eof())
  {
    if(line.find("CRYST1") != line.npos ||
       line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
    std::getline(pdbFile, line);
  }

  if(line.find("CRYST1") != line.npos)
  {
    // CRYST1 record found: read unit cell information
    std::istringstream crystLineStream(line);
    std::string dummyString;                      // CRYST1 string
    Real r_a, r_b, r_c, r_alpha, r_beta, r_gamma; // Cell parameters read

    crystLineStream >> dummyString >> r_a >> r_b >> r_c >> r_alpha >> r_beta >> r_gamma;
    if(crystLineStream.fail())
    {  
      std::cerr << "Error reading CRYST1 record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    system.setLattice(Lattice(r_a, r_b, r_c, r_alpha, r_beta, r_gamma));
  }

  // Search for ATOM/HETATM records
  while(!pdbFile.eof())
  {
    if(line.find("ATOM") != line.npos ||
       line.find("HETATM") != line.npos) break;
    std::getline(pdbFile, line);
  }
  if(pdbFile.eof())
  {
    std::cerr << "Error reading pdb file " << pdbFile.fileName()
              << ": no ATOM/HETATM records" << std::endl;
  }

  // Read ATOM/HETATM records

  Molecule currentResidue;              // Current molecule
  size_t iAtom = 0;                     // Index of current atom
  size_t iRes = 0;                      // Index of current residue
  size_t previousID = 0;                // Previous residue ID read
  std::vector<size_t> firstAtoms(1, 0); // Global index of the first atom in each molecule
                                        // (needed to convert local to global indices)

  while(line.find("ATOM") != line.npos ||
        line.find("HETATM") != line.npos ||
        line.find("TER") != line.npos)
  {
    if(line.find("TER") != line.npos) // Ignore TER records
    {
      std::getline(pdbFile, line);
      continue;
    }
 
    size_t r_serial;        // Atom serial number read from file
    std::string r_roleName; // Atom role name read from file
    std::string r_symbol;   // Symbol read from file
    Real r_x, r_y, r_z;     // Coordinates read from file
    std::string r_resName;  // Residue name read from file
    size_t r_resID;         // Residue ID number read from file
    Real r_tempFactor;      // Temperature factor read from file
    std::istringstream serialStream(line.substr(6, 5));
    serialStream >> r_serial;
    if(serialStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(r_serial != iAtom + 1)
    {
      std::cerr << "Error reading pdb file " << pdbFile.fileName()
                << ": atom serial numbers not correlated" << std::endl;
      return pdbFile;
    }
    r_roleName = line.substr(12, 4);
    // Find the atomic symbol
    size_t symbolIndex = r_roleName.find_first_not_of(" 1234567890");
    if(symbolIndex == r_roleName.npos)
    {
      std::cerr << "Error reading ATOM/HETATOM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    r_symbol = r_roleName.at(symbolIndex); // Note that this still works only for one-letter symbols
    /*
    // Old version below
    r_symbol = (line.substr(12, 1).find_first_of(" 1234567890") != line.substr(12, 1).npos)?
                ((line.substr(13, 1) == " ")?line.substr(14, 1):line.substr(13, 1)):
                line.substr(12, 2);
    if(r_symbol == " ")
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    */
    std::istringstream coordStream(line.substr(30, 24));
    coordStream >> r_x >> r_y >> r_z;
    if(coordStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    std::istringstream tempFactorStream(line.substr(60, 6));
    tempFactorStream >> r_tempFactor;
    if(tempFactorStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    r_resName = line.substr(17, 4);
    if(!currentResidue.name().size()) currentResidue.setName(r_resName);  // Use same name for residue
    std::istringstream resIDStream(line.substr(22, 5));
    resIDStream >> r_resID;
    if(resIDStream.fail())
    {
      std::cerr << "Error reading ATOM/HETATM record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(iAtom == 0) previousID = r_resID;  // Avoid trouble when first residue ID is not zero
    if(iAtom > 0 && r_resID != previousID)
    {
      // Finished reading current residue
      if(newSystem)
        system.add(currentResidue);
      else
      {
        if(iRes >= system.size())
        {
          std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
          return pdbFile;
        }
        if(currentResidue.numAtoms() != system[iRes].numAtoms())
        {
          std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
          return pdbFile;
        }

        for(size_t i = 0; i < currentResidue.numAtoms(); ++i)
        {
          if(currentResidue.atom(i).symbol != system[iRes].atom(i).symbol ||
             currentResidue.atom(i).residueName != system[iRes].atom(i).residueName ||
             currentResidue.atom(i).residueID != system[iRes].atom(i).residueID ||
             currentResidue.atom(i).roleName != system[iRes].atom(i).roleName)
          {
            std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
            return pdbFile;
          }
          system[iRes].setPosition(i, currentResidue.atom(i).position);
          system[iRes].setTemperatureFactor(i, currentResidue.atom(i).tempFactor);
        }
      }
      firstAtoms.push_back(iAtom);
      currentResidue.clear();
      previousID = r_resID;
      ++iRes;
    }
    currentResidue.addAtom(Atom(r_symbol, Vector3D(r_x, r_y, r_z), r_resID, 
                                r_resName, r_roleName, r_tempFactor));
    ++iAtom;
    std::getline(pdbFile, line);
    if(pdbFile.eof()) break;
  } // End of ATOM/HETATM sections

  if(newSystem)
    system.add(currentResidue);
  else
  {
    if(iRes >= system.size())
    {
      std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    if(currentResidue.numAtoms() != system[iRes].numAtoms())
    {
      std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }

    for(size_t i = 0; i < currentResidue.numAtoms(); ++i)
    {
      if(currentResidue.atom(i).symbol != system[iRes].atom(i).symbol ||
          currentResidue.atom(i).residueName != system[iRes].atom(i).residueName ||
          currentResidue.atom(i).residueID != system[iRes].atom(i).residueID ||
          currentResidue.atom(i).roleName != system[iRes].atom(i).roleName)
      {
        std::cerr << "Inconsistent data in pdb file " << pdbFile.fileName() << std::endl;
        return pdbFile;
      }
      system[iRes].setPosition(i, currentResidue.atom(i).position);
      system[iRes].setTemperatureFactor(i, currentResidue.atom(i).tempFactor);
    }
  }
  firstAtoms.push_back(iAtom);

  if(!newSystem) return pdbFile;  // Connectivity read only for new system

  // Search for CONECT records
  while(!pdbFile.eof())
  {
    if(line.find("CONECT") != line.npos) break;
    if(line.find("END") != line.npos) break; // In case there are no CONECT records
    std::getline(pdbFile, line);
  }
  if(pdbFile.eof() || (line.find("END") != line.npos))
  {
    // No connect records - find bonds automatically if this is a new system
    if(newSystem)
      for(size_t i = 0; i < system.size(); ++i)
        system[i].calculateBonds();
    return pdbFile;
  }

  // Read CONECT records
  while(line.find("CONECT") != line.npos)
  {
    size_t r_atom1, r_atom2; // Atoms forming a bond read from file
    std::istringstream atom1Stream(line.substr(6, 5));
    atom1Stream >> r_atom1;
    if(atom1Stream.fail())
    {
      std::cerr << "Error reading CONECT record from file " << pdbFile.fileName() << std::endl;
      return pdbFile;
    }
    for(size_t i = 11; i < 27; i += 5)
    {
      std::istringstream atom2Stream(line.substr(i, 5));
      atom2Stream >> r_atom2;
      if(atom2Stream.fail()) break;
      if(r_atom2 > r_atom1) // CONECTs are redundant
      {
        // Find which molecule the bond belongs to
        size_t geometryIndex = 0;
        while(r_atom1 > firstAtoms[geometryIndex+1])
          ++geometryIndex;
        if(r_atom2 <= firstAtoms[geometryIndex + 1]) // Must be intra-molecule
          system[geometryIndex].addBond(Bond(r_atom1 - 1 - firstAtoms[geometryIndex],
                                             r_atom2 - 1 - firstAtoms[geometryIndex]));
      }
    }
    std::getline(pdbFile, line);
    if(pdbFile.eof()) break;
  }
  return pdbFile;
}
