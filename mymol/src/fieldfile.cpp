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
** FIELD file format
*/

#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include "common/include/assert.h"
#include "mymol/include/file_formats/fileformats.h"
#include "mymol/include/pairpotential.h"

// Constructors

FIELDFile::FIELDFile()
:
IOFile(), unitsLabel_("internal")
{}

FIELDFile::FIELDFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode), unitsLabel_("internal")
{}

// Operators

FIELDFile &operator<<(FIELDFile &fieldFile, Geometry const &geometry)
{
  assert(fieldFile.is_open() && fieldFile.mode() == OUT);
  
  if(geometry.numAtoms() == 0)
  {
    std::cerr << "Error writing to FIELD file: empty geometry" << std::endl;
    return fieldFile;
  }

  // Store precision
  std::streamsize precision = fieldFile.precision();
  fieldFile.precision(WRITE_PRECISION);

  // Header and units records
  fieldFile << geometry.name() << std::endl;
  fieldFile << "units " << fieldFile.unitsLabel_ << std::endl;
  
  // Molecule record (only one for Geometry objects)
  fieldFile << "molecular types 1" << std::endl;
  fieldFile << geometry.name() << std::endl;
  fieldFile << "nummols 1" << std::endl;

  // Atom records
  fieldFile << "atoms " << geometry.numAtoms() << std::endl;

  Atom previousAtom = geometry.atom(0); // Previous atom type
  size_t numReps = 0;                   // Number of repetitions

  for(size_t i = 0; i < geometry.numAtoms(); ++i)
  {
    Atom const &atom = geometry.atom(i);
    if(previousAtom.symbol == atom.symbol) ++numReps;
    else
    {
      fieldFile << previousAtom.symbol << " "
                << std::setiosflags(std::ios::fixed)
                << previousAtom.mass << " "
                << 0.0 << " " << numReps << std::endl;
      previousAtom = atom;
      numReps = 1;
    }
  }
  // Write the last atom type
  fieldFile << previousAtom.symbol << " "
            << std::setiosflags(std::ios::fixed)
            << previousAtom.mass << " "
            << 0.0 << " " << numReps << std::endl;

  // Separate flexible and constraint bonds
  size_t numFlexible = 0;
  size_t numConstraint = 0;
  size_t numUnknown = 0;

  static const size_t numFlexibleBondTypes = 13;
  static const std::string flexibleBondTypes[numFlexibleBondTypes] =
  {
    "harm", "-hrm", "mors", "-mrs", "12-6", "-126", "rhrm", "-rhm",
    "quar", "-qur", "buck", "-bck", "fene"
  };
  for(size_t i = 0; i < geometry.numBonds(); ++i)
  {
    std::string const &type = geometry.bond(i).type;
    bool isFlexible = false;
    for(size_t j = 0; j < numFlexibleBondTypes; ++j)
    {
      if(type.find(flexibleBondTypes[j]) != type.npos)
      {
        ++numFlexible;
        isFlexible = true;
        break;
      }
    }
    if(!isFlexible)
    {
      if(type.find("cons") != type.npos) ++numConstraint;
      else ++numUnknown;
    }
  }
  if(numUnknown > 0)
    std::cerr << "Warning from FIELDFile - some bond types being written to file "
              << fieldFile.fileName() << " are not standard DL_POLY types" << std::endl;

  // Bond records
  if(numFlexible + numUnknown > 0)
  {
    fieldFile << "bonds " << numFlexible + numUnknown << std::endl;

    for(size_t i = 0; i < geometry.numBonds(); ++i)
    {
      Bond const &bond = geometry.bond(i);
      if(bond.type.find("cons") != bond.type.npos) continue;
      fieldFile << bond.type << " "
                << bond.first + 1 << " "
                << bond.second + 1;

      for(size_t j = 0; j < bond.parameters.size(); ++j)
        fieldFile << " " << bond.parameters[j];
      fieldFile << std::endl;
    }
  }
  // Constraint records
  if(numConstraint > 0)
  {
    fieldFile << "constraints " << numConstraint << std::endl;

    for(size_t i = 0; i < geometry.numBonds(); ++i)
    {
      Bond const &bond = geometry.bond(i);
      if(bond.type.find("cons") == bond.type.npos) continue;
      Real const parameter = (bond.parameters.size() > 0)?bond.parameters[0]:0.0;
      fieldFile << bond.first + 1 << " "
                << bond.second + 1 << " "
                << parameter << std::endl;
    }
  }

  fieldFile << "finish" << std::endl;
  fieldFile.precision(precision); // Return to original value

  return fieldFile;
}

FIELDFile &operator>>(FIELDFile &fieldFile, Geometry &geometry)
{
  assert(fieldFile.is_open() && fieldFile.mode() == IN);
  geometry.clear();

  fieldFile.seekg(0, std::ios::beg); // Rewind file (in case it was read previously)

  std::string line; // Line read from file

  // Read header record
  char r_label[LINE_LENGTH]; // Header read from file
  fieldFile.getline(r_label, LINE_LENGTH);
  geometry.setName(std::string(r_label)); 

  // Read units record
  std::getline(fieldFile, line);
  if(line.find("unit") == line.npos)
  {
    std::cerr << "Error reading file " << fieldFile.fileName()
              << ": missing units record" << std::endl;
    return fieldFile;
  }
  std::stringstream unitStream(line);
  std::string keyword;
  unitStream >> keyword; // "units"
  unitStream >> keyword; // Units keyword
  fieldFile.setUnitsLabel(keyword);

  // Search for molecule types record
  do
  {
    std::getline(fieldFile, line);
    if(line.find("molecul") != line.npos) break;
  }
  while(!fieldFile.eof());
 
  if(fieldFile.eof())
  {
    std::cerr << "Error reading file " << fieldFile.fileName()
              << ": molecule types record not found." << std::endl;
    return fieldFile;
  }

  // Read molecule types
  size_t numMolTypes = 0; // Number of molecule types
  size_t numMolTypesIndex = line.find_first_of("-.1234567890");
  
  if(numMolTypesIndex == line.npos)
  {
    std::cerr << "Error reading number of molecule types from FIELD file "
              << fieldFile.fileName() << std::endl;
    return fieldFile;
  }

  std::stringstream numMolTypesStream(line.substr(numMolTypesIndex));
  numMolTypesStream >> numMolTypes;
   
  if(numMolTypesStream.fail() || numMolTypes < 1)
  {
    std::cerr << "Error reading number of molecule types from FIELD file " 
              << fieldFile.fileName() << std::endl;
    return fieldFile;
  } 

  // Read molecule fields
  for(size_t iMolType = 0; iMolType < numMolTypes; ++iMolType)
  {
    std::string keyword; // Keyword read from file

    // Get molecule name
    std::string moleculeName;
    std::getline(fieldFile, moleculeName);

    // Get the number of molecules
    size_t numMols = 0;
    std::getline(fieldFile, line);
    if(line.find("nummols") == line.npos)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName() 
                << ": nummols record not found" << std::endl;
      return fieldFile;
    }
    std::stringstream numMolsStream(line);
    numMolsStream >> keyword >> numMols;
    if(numMolsStream.fail() || numMols < 1)
    {
      std::cerr << "Error reading nummols record from FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }
    
    // Get the number of atoms in molecule
    size_t numAtomsInMol = 0;
    std::getline(fieldFile, line);
    if(line.find("atoms") == line.npos)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName()
                << ": atoms record not found" << std::endl;
      return fieldFile;
    }
    std::stringstream numAtomsStream(line);
    numAtomsStream >> keyword >> numAtomsInMol;
    if(numAtomsStream.fail() || numAtomsInMol < 1)
    {
      std::cerr << "Error reading atoms record from FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }
    
    // Read the atom information
    Geometry currentMolecule; // The reference geometry to be replicated
    currentMolecule.setName(moleculeName);
    size_t iAtom = 0;
    while(iAtom < numAtomsInMol)
    {
      // These variables follow the notation from the DL_POLY manual
      std::string sitnam; // Atomic symbol/name
      Real weight;        // Atomic mass
      Real chge;          // Atomic charge
      size_t nrept;       // Number of times the atom is repeated

      std::getline(fieldFile, line);
      std::stringstream atomStream(line);
      atomStream >> sitnam >> weight >> chge; 
      if(atomStream.fail() || weight < 0.0)
      {
        std::cerr << "Error reading atom record from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }
      Atom currentAtom(sitnam, 0.0);
      atomStream >> nrept;
      if(atomStream.fail()) nrept = 1; // Assume only one atom
      
      for(size_t i = 0; i < nrept; ++i)
        currentMolecule.addAtom(currentAtom);
      iAtom += nrept;
    }
   
    // Consistency check
    if(currentMolecule.numAtoms() != numAtomsInMol)
    {
      std::cerr << "Inconsistent molecule record in FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }

    // Search for bond records
    bool found_bonds = false;
    bool found_constraints = false;
    bool found_finish = false;
    do
    {
      std::getline(fieldFile, line);
      found_bonds = (line.find("bond") != line.npos);
      if(found_bonds) break;
      found_constraints = (line.find("constr") != line.npos);
      if(found_constraints) break;
      found_finish = (line.find("finish") != line.npos);
      if(found_finish) break;
    }
    while(!fieldFile.eof());

    // Read bond records if found
    if(found_bonds)
    {
      // Get the number of bond records
      size_t numBondRecords;
      std::stringstream bondTitleStream(line);
      bondTitleStream >> keyword >> numBondRecords;
      if(bondTitleStream.fail())
      {
        std::cerr << "Error reading number of bond records from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }
      
      // Read bond records
      std::string bondType; // Bond key
      size_t atom1, atom2;  // Atoms in bond
      Real parameter;       // Bond parameter
     
      for(size_t iBond = 0; iBond < numBondRecords; ++iBond)
      {
        std::getline(fieldFile, line);
        std::stringstream bondRecordStream(line);
        bondRecordStream >> bondType >> atom1 >> atom2;
        if(bondRecordStream.fail())
        {
          std::cerr << "Error reading bond record from FIELD file " << fieldFile.fileName() << std::endl;
          return fieldFile; 
        }
        std::vector<Real> bondParameters;
        while(!bondRecordStream.eof())
        {
          bondRecordStream >> parameter;
          if(bondRecordStream.eof() && bondRecordStream.fail()) break;
          if(bondRecordStream.fail())
          {
            std::cerr << "Error reading bond parameters from FIELD file " << fieldFile.fileName() << std::endl;
            return fieldFile;
          }
          bondParameters.push_back(parameter);
        }
        currentMolecule.addBond(Bond(atom1 - 1, atom2 - 1, bondType, bondParameters));
      }
      
      // Done reading bonds, check for constraints record
      do
      {
        std::getline(fieldFile, line);
        found_constraints = (line.find("constr") != line.npos);
        if(found_constraints) break;
        found_finish = (line.find("finish") != line.npos);
        if(found_finish) break;
      }
      while(!fieldFile.eof());
    }

    // Read constraint records if found
    if(found_constraints)
    {
      // Get the number of constraint records
      size_t numConstraintRecords;
      std::stringstream constraintTitleStream(line);
      constraintTitleStream >> keyword >> numConstraintRecords;
      if(constraintTitleStream.fail())
      {
        std::cerr << "Error reading number of contraint records from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }

      // Read constraint records
      size_t atom1, atom2; // Atoms in bond
      Real parameter;      // Bond parameter (bond length)
     
      for(size_t iBond = 0; iBond < numConstraintRecords; ++iBond)
      {
        std::getline(fieldFile, line);
        std::stringstream constraintRecordStream(line);
        constraintRecordStream >> atom1 >> atom2 >> parameter; 
        if(constraintRecordStream.fail())
        {
          std::cerr << "Error reading constraint record from FIELD file " << fieldFile.fileName() << std::endl;
          return fieldFile;
        }
        std::vector<Real> bondParameter(1, parameter);
        currentMolecule.addBond(Bond(atom1 - 1, atom2 - 1, "cons", bondParameter));
      } 
    }
  
    // Skip to finish record
    if(!found_finish)
    do
    { 
      std::getline(fieldFile, line);
      found_finish = (line.find("finish") != line.npos);
      if(found_finish) break;
    }
    while(!fieldFile.eof());
   
    if(!found_finish)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName()
                << ": finish record not found" << std::endl;
      return fieldFile;
    }

    // Replicate molecules
    size_t iGlobal = geometry.numAtoms(); // Used to convert local atom indices to global
    for(size_t iMol = 0; iMol < numMols; ++iMol)
    {
      // Atoms
      for(size_t iAtom = 0; iAtom < numAtomsInMol; ++iAtom)
        geometry.addAtom(currentMolecule.atom(iAtom));
      // Bonds
      size_t numTotalBonds = currentMolecule.numBonds();
      for(size_t iBond = 0; iBond < numTotalBonds; ++iBond)
      {
        size_t iShift = iGlobal + iMol*numAtomsInMol; // To convert local atom indices to global
        Bond const &currentBond = currentMolecule.bond(iBond);
        geometry.addBond(Bond(iShift + currentBond.first,
                              iShift + currentBond.second,
                              currentBond.type,
                              currentBond.parameters));
      }
    } 
  } // End of loop to read molecule fields

  return fieldFile;
}

FIELDFile &operator<<(FIELDFile &fieldFile, System<Geometry> const &system)
{
  assert(fieldFile.is_open() && fieldFile.mode() == OUT);
  
  if(system.size() == 0)
  {
    std::cerr << "Error writing to FIELD file: empty system" << std::endl;
    return fieldFile;
  }

  // Header and units records
  fieldFile << system.name() << std::endl;
  fieldFile << "units " << fieldFile.unitsLabel_ << std::endl; 

   // Store precision
  std::streamsize precision = fieldFile.precision();
  fieldFile.precision(WRITE_PRECISION);

  // Find distinct sets of molecules
  std::vector<size_t> distinctMoleculeIndices(1, 0);
  std::vector<size_t> distinctMoleculeCounts(1, 1);
  std::string previousName = system[0].name();
  for(size_t iMol = 1; iMol < system.size(); ++iMol)
  {
    if(system[iMol].name() != previousName)
    {
      distinctMoleculeIndices.push_back(iMol);
      distinctMoleculeCounts.push_back(1);
      previousName = system[iMol].name();
    }
    else
      ++distinctMoleculeCounts[distinctMoleculeCounts.size() - 1];
  }
  // Molecule records
  fieldFile << "molecular types " << distinctMoleculeIndices.size() << std::endl;
  for(size_t index = 0; index < distinctMoleculeIndices.size(); ++index)
  {
    size_t const iMol = distinctMoleculeIndices[index];
    Geometry const &geometry = system[iMol];
    fieldFile << geometry.name() << std::endl;
    fieldFile << "nummols " << distinctMoleculeCounts[index] << std::endl;

    // Atom records
    fieldFile << "atoms " << geometry.numAtoms() << std::endl;

    Atom previousAtom = geometry.atom(0); // Previous atom type
    size_t numReps = 0;                   // Number of repetitions

    for(size_t i = 0; i < geometry.numAtoms(); ++i)
    {
      Atom const &atom = geometry.atom(i);
      if(previousAtom.symbol == atom.symbol) ++numReps;
      else
      {
        fieldFile << previousAtom.symbol << " "
                  << std::setiosflags(std::ios::fixed)
                  << previousAtom.mass << " "
                  << 0.0 << " " << numReps << std::endl;
        previousAtom = atom;
        numReps = 1;
      }
    }
    // Write the last atom type
    fieldFile << previousAtom.symbol << " "
              << std::setiosflags(std::ios::fixed)
              << previousAtom.mass << " "
              << 0.0 << " " << numReps << std::endl;

    // Separate flexible and constraint bonds
    size_t numFlexible = 0;
    size_t numConstraint = 0;
    size_t numUnknown = 0;

    static const size_t numFlexibleBondTypes = 13;
    static const std::string flexibleBondTypes[numFlexibleBondTypes] =
    {
      "harm", "-hrm", "mors", "-mrs", "12-6", "-126", "rhrm", "-rhm",
      "quar", "-qur", "buck", "-bck", "fene"
    };
    for(size_t i = 0; i < geometry.numBonds(); ++i)
    {
      std::string const &type = geometry.bond(i).type;
      bool isFlexible = false;
      for(size_t j = 0; j < numFlexibleBondTypes; ++j)
      {
        if(type.find(flexibleBondTypes[j]) != type.npos)
        {
          ++numFlexible;
          isFlexible = true;
          break;
        }
      }
      if(!isFlexible)
      {
        if(type.find("cons") != type.npos) ++numConstraint;
        else ++numUnknown;
      }
    }
    if(numUnknown > 0)
      std::cerr << "Warning from FIELDFile - some bond types being written to file "
                << fieldFile.fileName() << " are not standard DL_POLY types" << std::endl;

    // Bond records
    if(numFlexible + numUnknown > 0)
    {
      fieldFile << "bonds " << numFlexible + numUnknown << std::endl;

      for(size_t i = 0; i < geometry.numBonds(); ++i)
      {
        Bond const &bond = geometry.bond(i);
        if(bond.type.find("cons") != bond.type.npos) continue;
        fieldFile << bond.type << " "
                  << bond.first + 1 << " "
                  << bond.second + 1;

        for(size_t j = 0; j < bond.parameters.size(); ++j)
          fieldFile << " " << bond.parameters[j];
        fieldFile << std::endl;
      }
    }
    // Constraint records
    if(numConstraint > 0)
    {
      fieldFile << "constraints " << numConstraint << std::endl;

      for(size_t i = 0; i < geometry.numBonds(); ++i)
      {
        Bond const &bond = geometry.bond(i);
        if(bond.type.find("cons") == bond.type.npos) continue;
        Real const parameter = (bond.parameters.size() > 0)?bond.parameters[0]:0.0;
        fieldFile << bond.first + 1 << " "
                  << bond.second + 1 << " "
                  << parameter << std::endl;
      }
    }
    fieldFile << "finish" << std::endl;
  }
  fieldFile.precision(precision); // Return to original value

  return fieldFile;
}

FIELDFile &operator>>(FIELDFile &fieldFile, System<Geometry> &system)
{
  assert(fieldFile.is_open() && fieldFile.mode() == IN);
  system.clear();

  fieldFile.seekg(0, std::ios::beg); // Rewind file (in case it was read previously)

  std::string line; // Line read from file

  // Read header record
  char r_label[LINE_LENGTH]; // Header read from file
  fieldFile.getline(r_label, LINE_LENGTH);
  system.setName(std::string(r_label)); 

  // Read units record
  std::getline(fieldFile, line);
  if(line.find("unit") == line.npos)
  {
    std::cerr << "Error reading file " << fieldFile.fileName()
              << ": missing units record" << std::endl;
    return fieldFile;
  }
  std::stringstream unitStream(line);
  std::string keyword;
  unitStream >> keyword; // "units"
  unitStream >> keyword; // Units keyword
  fieldFile.setUnitsLabel(keyword);

  // Search for molecule types record
  do
  {
    std::getline(fieldFile, line);
    if(line.find("molecul") != line.npos) break;
  }
  while(!fieldFile.eof());
 
  if(fieldFile.eof())
  {
    std::cerr << "Error reading file " << fieldFile.fileName()
              << ": molecule types record not found." << std::endl;
    return fieldFile;
  }

  // Read molecule types
  size_t numMolTypes = 0; // Number of molecule types
  size_t numMolTypesIndex = line.find_first_of("-.1234567890");
  
  if(numMolTypesIndex == line.npos)
  {
    std::cerr << "Error reading number of molecule types from FIELD file "
              << fieldFile.fileName() << std::endl;
    return fieldFile;
  }

  std::stringstream numMolTypesStream(line.substr(numMolTypesIndex));
  numMolTypesStream >> numMolTypes;
   
  if(numMolTypesStream.fail() || numMolTypes < 1)
  {
    std::cerr << "Error reading number of molecule types from FIELD file " 
              << fieldFile.fileName() << std::endl;
    return fieldFile;
  } 

  // Read molecule fields
  for(size_t iMolType = 0; iMolType < numMolTypes; ++iMolType)
  {
    std::string keyword; // Keyword read from file

    // Get molecule name
    std::string moleculeName;
    std::getline(fieldFile, moleculeName);

    // Get the number of molecules
    size_t numMols = 0;
    std::getline(fieldFile, line);
    if(line.find("nummols") == line.npos)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName() 
                << ": nummols record not found" << std::endl;
      return fieldFile;
    }
    std::stringstream numMolsStream(line);
    numMolsStream >> keyword >> numMols;
    if(numMolsStream.fail() || numMols < 1)
    {
      std::cerr << "Error reading nummols record from FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }
    
    // Get the number of atoms in molecule
    size_t numAtomsInMol = 0;
    std::getline(fieldFile, line);
    if(line.find("atoms") == line.npos)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName()
                << ": atoms record not found" << std::endl;
      return fieldFile;
    }
    std::stringstream numAtomsStream(line);
    numAtomsStream >> keyword >> numAtomsInMol;
    if(numAtomsStream.fail() || numAtomsInMol < 1)
    {
      std::cerr << "Error reading atoms record from FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }
    
    // Read the atom information
    Geometry currentMolecule; // The reference geometry to be replicated
    currentMolecule.setName(moleculeName);
    size_t iAtom = 0;
    while(iAtom < numAtomsInMol)
    {
      // These variables follow the notation from the DL_POLY manual
      std::string sitnam; // Atomic symbol/name
      Real weight;        // Atomic mass
      Real chge;          // Atomic charge
      size_t nrept;       // Number of times the atom is repeated

      std::getline(fieldFile, line);
      std::stringstream atomStream(line);
      atomStream >> sitnam >> weight >> chge; 
      if(atomStream.fail() || weight < 0.0)
      {
        std::cerr << "Error reading atom record from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }
      Atom currentAtom(sitnam, 0.0);
      atomStream >> nrept;
      if(atomStream.fail()) nrept = 1; // Assume only one atom
      
      for(size_t i = 0; i < nrept; ++i)
        currentMolecule.addAtom(currentAtom);
      iAtom += nrept;
    }
   
    // Consistency check
    if(currentMolecule.numAtoms() != numAtomsInMol)
    {
      std::cerr << "Inconsistent molecule record in FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }

    // Search for bond records
    bool found_bonds = false;
    bool found_constraints = false;
    bool found_finish = false;
    do
    {
      std::getline(fieldFile, line);
      found_bonds = (line.find("bond") != line.npos);
      if(found_bonds) break;
      found_constraints = (line.find("constr") != line.npos);
      if(found_constraints) break;
      found_finish = (line.find("finish") != line.npos);
      if(found_finish) break;
    }
    while(!fieldFile.eof());

    // Read bond records if found
    if(found_bonds)
    {
      // Get the number of bond records
      size_t numBondRecords;
      std::stringstream bondTitleStream(line);
      bondTitleStream >> keyword >> numBondRecords;
      if(bondTitleStream.fail())
      {
        std::cerr << "Error reading number of bond records from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }
      
      // Read bond records
      std::string bondType; // Bond key
      size_t atom1, atom2;  // Atoms in bond
      Real parameter;       // Bond parameter 
     
      for(size_t iBond = 0; iBond < numBondRecords; ++iBond)
      {
        std::getline(fieldFile, line);
        std::stringstream bondRecordStream(line);
        bondRecordStream >> bondType >> atom1 >> atom2;
        if(bondRecordStream.fail())
        {
          std::cerr << "Error reading bond record from FIELD file " << fieldFile.fileName() << std::endl;
          return fieldFile; 
        }
        std::vector<Real> bondParameters;
        while(!bondRecordStream.eof())
        {
          bondRecordStream >> parameter;
          if(bondRecordStream.eof() && bondRecordStream.fail()) break;
          if(bondRecordStream.fail())
          {
            std::cerr << "Error reading bond parameters from FIELD file " << fieldFile.fileName() << std::endl;
            return fieldFile;
          }
          bondParameters.push_back(parameter);
        }
        currentMolecule.addBond(Bond(atom1 - 1, atom2 - 1, bondType, bondParameters));
      }
      
      // Done reading bonds, check for constraints record
      do
      {
        std::getline(fieldFile, line);
        found_constraints = (line.find("constr") != line.npos);
        if(found_constraints) break;
        found_finish = (line.find("finish") != line.npos);
        if(found_finish) break;
      }
      while(!fieldFile.eof());
    }

    // Read constraint records if found
    if(found_constraints)
    {
      // Get the number of constraint records
      size_t numConstraintRecords;
      std::stringstream constraintTitleStream(line);
      constraintTitleStream >> keyword >> numConstraintRecords;
      if(constraintTitleStream.fail())
      {
        std::cerr << "Error reading number of contraint records from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }

      // Read constraint records
      size_t atom1, atom2; // Atoms in bond
      Real parameter;      // Bond parameter (bond length)
     
      for(size_t iBond = 0; iBond < numConstraintRecords; ++iBond)
      {
        std::getline(fieldFile, line);
        std::stringstream constraintRecordStream(line);
        constraintRecordStream >> atom1 >> atom2 >> parameter;
        if(constraintRecordStream.fail())
        {
          std::cerr << "Error reading constraint record from FIELD file " << fieldFile.fileName() << std::endl;
          return fieldFile;
        }
        std::vector<Real> bondParameter(1, parameter);
        currentMolecule.addBond(Bond(atom1 - 1, atom2 - 1, "cons", bondParameter));
      } 
    }
  
    // Skip to finish record
    if(!found_finish)
    do
    { 
      std::getline(fieldFile, line);
      found_finish = (line.find("finish") != line.npos);
      if(found_finish) break;
    }
    while(!fieldFile.eof());
   
    if(!found_finish)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName()
                << ": finish record not found" << std::endl;
      return fieldFile;
    }

    // Replicate molecules
    for(size_t iMol = 0; iMol < numMols; ++iMol)
      system.add(currentMolecule);

  } // End of loop to read molecule fields

  return fieldFile;
}

FIELDFile &operator<<(FIELDFile &fieldFile, System<Molecule> const &system)
{
  assert(fieldFile.is_open() && fieldFile.mode() == OUT);
  
  if(system.size() == 0)
  {
    std::cerr << "Error writing to FIELD file: empty system" << std::endl;
    return fieldFile;
  }

  // Header and units records
  fieldFile << system.name() << std::endl;
  fieldFile << "units " << fieldFile.unitsLabel_ << std::endl; 

   // Store precision
  std::streamsize precision = fieldFile.precision();
  fieldFile.precision(WRITE_PRECISION);

  // Find distinct sets of molecules
  std::vector<size_t> distinctMoleculeIndices(1, 0);
  std::vector<size_t> distinctMoleculeCounts(1, 1);
  std::string previousName = system[0].name();
  for(size_t iMol = 1; iMol < system.size(); ++iMol)
  {
    if(system[iMol].name() != previousName)
    {
      distinctMoleculeIndices.push_back(iMol);
      distinctMoleculeCounts.push_back(1);
      previousName = system[iMol].name();
    }
    else
      ++distinctMoleculeCounts[distinctMoleculeCounts.size() - 1];
  }
  // Molecule records
  fieldFile << "molecular types " << distinctMoleculeIndices.size() << std::endl;
  for(size_t index = 0; index < distinctMoleculeIndices.size(); ++index)
  {
    size_t const iMol = distinctMoleculeIndices[index];
    Molecule const &molecule = system[iMol];
    fieldFile << molecule.name() << std::endl;
    fieldFile << "nummols " << distinctMoleculeCounts[index] << std::endl;

    // Atom records
    fieldFile << "atoms " << molecule.numAtoms() << std::endl;

    Atom previousAtom = molecule.atom(0); // Previous atom type
    size_t numReps = 0;                   // Number of repetitions

    for(size_t i = 0; i < molecule.numAtoms(); ++i)
    {
      Atom const &atom = molecule.atom(i);
      if(previousAtom.symbol == atom.symbol) ++numReps;
      else
      {
        fieldFile << previousAtom.symbol << " "
                  << std::setiosflags(std::ios::fixed)
                  << previousAtom.mass << " "
                  << 0.0 << " " << numReps << std::endl;
        previousAtom = atom;
        numReps = 1;
      }
    }
    // Write the last atom type
    fieldFile << previousAtom.symbol << " "
              << std::setiosflags(std::ios::fixed)
              << previousAtom.mass << " "
              << 0.0 << " " << numReps << std::endl;

    // Separate flexible and constraint bonds
    size_t numFlexible = 0;
    size_t numConstraint = 0;
    size_t numUnknown = 0;

    static const size_t numFlexibleBondTypes = 13;
    static const std::string flexibleBondTypes[numFlexibleBondTypes] =
    {
      "harm", "-hrm", "mors", "-mrs", "12-6", "-126", "rhrm", "-rhm",
      "quar", "-qur", "buck", "-bck", "fene"
    };
    for(size_t i = 0; i < molecule.numBonds(); ++i)
    {
      std::string const &type = molecule.bond(i).type;
      bool isFlexible = false;
      for(size_t j = 0; j < numFlexibleBondTypes; ++j)
      {
        if(type.find(flexibleBondTypes[j]) != type.npos)
        {
          ++numFlexible;
          isFlexible = true;
          break;
        }
      }
      if(!isFlexible)
      {
        if(type.find("cons") != type.npos) ++numConstraint;
        else ++numUnknown;
      }
    }
    if(numUnknown > 0)
      std::cerr << "Warning from FIELDFile - some bond types being written to file "
                << fieldFile.fileName() << " are not standard DL_POLY types" << std::endl;

    // Bond records
    if(numFlexible + numUnknown > 0)
    {
      fieldFile << "bonds " << numFlexible + numUnknown << std::endl;

      for(size_t i = 0; i < molecule.numBonds(); ++i)
      {
        Bond const &bond = molecule.bond(i);
        if(bond.type.find("cons") != bond.type.npos) continue;
        fieldFile << bond.type << " "
                  << bond.first + 1 << " "
                  << bond.second + 1;

        for(size_t j = 0; j < bond.parameters.size(); ++j)
          fieldFile << " " << bond.parameters[j];
        fieldFile << std::endl;
      }
    }
    // Constraint records
    if(numConstraint > 0)
    {
      fieldFile << "constraints " << numConstraint << std::endl;

      for(size_t i = 0; i < molecule.numBonds(); ++i)
      {
        Bond const &bond = molecule.bond(i);
        if(bond.type.find("cons") == bond.type.npos) continue;
        Real const parameter = (bond.parameters.size() > 0)?bond.parameters[0]:0.0;
        fieldFile << bond.first + 1 << " "
                  << bond.second + 1 << " "
                  << parameter << std::endl;
      }
    }
    fieldFile << "finish" << std::endl;
  }
  fieldFile.precision(precision); // Return to original value

  return fieldFile;
}

FIELDFile &operator>>(FIELDFile &fieldFile, System<Molecule> &system)
{
  assert(fieldFile.is_open() && fieldFile.mode() == IN);
  system.clear();

  fieldFile.seekg(0, std::ios::beg); // Rewind file (in case it was read previously)

  std::string line; // Line read from file

  // Read header record
  char r_label[LINE_LENGTH]; // Header read from file
  fieldFile.getline(r_label, LINE_LENGTH);
  system.setName(std::string(r_label)); 

  // Read units record
  std::getline(fieldFile, line);
  if(line.find("unit") == line.npos)
  {
    std::cerr << "Error reading file " << fieldFile.fileName()
              << ": missing units record" << std::endl;
    return fieldFile;
  }
  std::stringstream unitStream(line);
  std::string keyword;
  unitStream >> keyword; // "units"
  unitStream >> keyword; // Units keyword
  fieldFile.setUnitsLabel(keyword);

  // Search for molecule types record
  do
  {
    std::getline(fieldFile, line);
    if(line.find("molecul") != line.npos) break;
  }
  while(!fieldFile.eof());
 
  if(fieldFile.eof())
  {
    std::cerr << "Error reading file " << fieldFile.fileName()
              << ": molecule types record not found." << std::endl;
    return fieldFile;
  }

  // Read molecule types
  size_t numMolTypes = 0; // Number of molecule types
  size_t numMolTypesIndex = line.find_first_of("-.1234567890");
  
  if(numMolTypesIndex == line.npos)
  {
    std::cerr << "Error reading number of molecule types from FIELD file "
              << fieldFile.fileName() << std::endl;
    return fieldFile;
  }

  std::stringstream numMolTypesStream(line.substr(numMolTypesIndex));
  numMolTypesStream >> numMolTypes;
   
  if(numMolTypesStream.fail() || numMolTypes < 1)
  {
    std::cerr << "Error reading number of molecule types from FIELD file " 
              << fieldFile.fileName() << std::endl;
    return fieldFile;
  } 

  // Read molecule fields
  for(size_t iMolType = 0; iMolType < numMolTypes; ++iMolType)
  {
    std::string keyword; // Keyword read from file

    // Get molecule name
    std::string moleculeName;
    std::getline(fieldFile, moleculeName);

    // Get the number of molecules
    size_t numMols = 0;
    std::getline(fieldFile, line);
    if(line.find("nummols") == line.npos)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName() 
                << ": nummols record not found" << std::endl;
      return fieldFile;
    }
    std::stringstream numMolsStream(line);
    numMolsStream >> keyword >> numMols;
    if(numMolsStream.fail() || numMols < 1)
    {
      std::cerr << "Error reading nummols record from FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }
    
    // Get the number of atoms in molecule
    size_t numAtomsInMol = 0;
    std::getline(fieldFile, line);
    if(line.find("atoms") == line.npos)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName()
                << ": atoms record not found" << std::endl;
      return fieldFile;
    }
    std::stringstream numAtomsStream(line);
    numAtomsStream >> keyword >> numAtomsInMol;
    if(numAtomsStream.fail() || numAtomsInMol < 1)
    {
      std::cerr << "Error reading atoms record from FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }
    
    // Read the atom information
    Molecule currentMolecule; // The reference geometry to be replicated
    currentMolecule.setName(moleculeName);
    size_t iAtom = 0;
    while(iAtom < numAtomsInMol)
    {
      // These variables follow the notation from the DL_POLY manual
      std::string sitnam; // Atomic symbol/name
      Real weight;        // Atomic mass
      Real chge;          // Atomic charge
      size_t nrept;       // Number of times the atom is repeated

      std::getline(fieldFile, line);
      std::stringstream atomStream(line);
      atomStream >> sitnam >> weight >> chge; 
      if(atomStream.fail() || weight < 0.0)
      {
        std::cerr << "Error reading atom record from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }
      Atom currentAtom(sitnam, 0.0);
      atomStream >> nrept;
      if(atomStream.fail()) nrept = 1; // Assume only one atom
      
      for(size_t i = 0; i < nrept; ++i)
        currentMolecule.addAtom(currentAtom);
      iAtom += nrept;
    }
   
    // Consistency check
    if(currentMolecule.numAtoms() != numAtomsInMol)
    {
      std::cerr << "Inconsistent molecule record in FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }

    // Search for bond records
    bool found_bonds = false;
    bool found_constraints = false;
    bool found_finish = false;
    do
    {
      std::getline(fieldFile, line);
      found_bonds = (line.find("bond") != line.npos);
      if(found_bonds) break;
      found_constraints = (line.find("constr") != line.npos);
      if(found_constraints) break;
      found_finish = (line.find("finish") != line.npos);
      if(found_finish) break;
    }
    while(!fieldFile.eof());

    // Read bond records if found
    if(found_bonds)
    {
      // Get the number of bond records
      size_t numBondRecords;
      std::stringstream bondTitleStream(line);
      bondTitleStream >> keyword >> numBondRecords;
      if(bondTitleStream.fail())
      {
        std::cerr << "Error reading number of bond records from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }
      
      // Read bond records
      std::string bondType; // Bond key
      size_t atom1, atom2;  // Atoms in bond
      Real parameter;       // Bond parameter 
     
      for(size_t iBond = 0; iBond < numBondRecords; ++iBond)
      {
        std::getline(fieldFile, line);
        std::stringstream bondRecordStream(line);
        bondRecordStream >> bondType >> atom1 >> atom2;
        if(bondRecordStream.fail())
        {
          std::cerr << "Error reading bond record from FIELD file " << fieldFile.fileName() << std::endl;
          return fieldFile; 
        }
        std::vector<Real> bondParameters;
        while(!bondRecordStream.eof())
        {
          bondRecordStream >> parameter;
          if(bondRecordStream.eof() && bondRecordStream.fail()) break;
          if(bondRecordStream.fail())
          {
            std::cerr << "Error reading bond parameters from FIELD file " << fieldFile.fileName() << std::endl;
            return fieldFile;
          }
          bondParameters.push_back(parameter);
        }
        currentMolecule.addBond(Bond(atom1 - 1, atom2 - 1, bondType, bondParameters));
      }
      
      // Done reading bonds, check for constraints record
      do
      {
        std::getline(fieldFile, line);
        found_constraints = (line.find("constr") != line.npos);
        if(found_constraints) break;
        found_finish = (line.find("finish") != line.npos);
        if(found_finish) break;
      }
      while(!fieldFile.eof());
    }

    // Read constraint records if found
    if(found_constraints)
    {
      // Get the number of constraint records
      size_t numConstraintRecords;
      std::stringstream constraintTitleStream(line);
      constraintTitleStream >> keyword >> numConstraintRecords;
      if(constraintTitleStream.fail())
      {
        std::cerr << "Error reading number of contraint records from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }

      // Read constraint records
      size_t atom1, atom2; // Atoms in bond
      Real parameter;      // Bond parameter (bond length)
     
      for(size_t iBond = 0; iBond < numConstraintRecords; ++iBond)
      {
        std::getline(fieldFile, line);
        std::stringstream constraintRecordStream(line);
        constraintRecordStream >> atom1 >> atom2 >> parameter;
        if(constraintRecordStream.fail())
        {
          std::cerr << "Error reading constraint record from FIELD file " << fieldFile.fileName() << std::endl;
          return fieldFile;
        }
        std::vector<Real> bondParameter(1, parameter);
        currentMolecule.addBond(Bond(atom1 - 1, atom2 - 1, "cons", bondParameter));
      } 
    }
  
    // Skip to finish record
    if(!found_finish)
    do
    { 
      std::getline(fieldFile, line);
      found_finish = (line.find("finish") != line.npos);
      if(found_finish) break;
    }
    while(!fieldFile.eof());
   
    if(!found_finish)
    {
      std::cerr << "Error reading FIELD file " << fieldFile.fileName()
                << ": finish record not found" << std::endl;
      return fieldFile;
    }

    // Replicate molecules
    for(size_t iMol = 0; iMol < numMols; ++iMol)
      system.add(currentMolecule);

  } // End of loop to read molecule fields

  return fieldFile;
}

FIELDFile &operator<<(FIELDFile &fieldFile, std::vector<PairPotential> const &potentials)
{
  assert(fieldFile.is_open() && fieldFile.mode() == OUT);

  if(potentials.size() < 1) 
  {
    std::cerr << "Error writing pair potentials to FIELD file " << fieldFile.fileName()
              << ": no potentials to write" << std::endl;
    return fieldFile; 
  }

  // Store precision
  std::streamsize precision = fieldFile.precision();
  fieldFile.precision(WRITE_PRECISION);

  fieldFile << "vdw " << potentials.size() << std::endl;

  // Store field widths for formatting (for readability, not really necessary)
  size_t maxLength = 0; // Length of longest atom type or potential type field

  for(size_t i = 0; i < potentials.size(); ++i)
  {
    PairPotential const &potential = potentials[i];
    if(potential.first.length() > maxLength) maxLength = potential.first.length();
    if(potential.second.length() > maxLength) maxLength = potential.second.length();
    if(maxLength < 4) maxLength = 4;
  }

  size_t nameWidth = maxLength + 2; // Width of name records 

  for(size_t i = 0; i < potentials.size(); ++i)
  {
    PairPotential const &potential = potentials[i];
    fieldFile << std::setiosflags(std::ios::right)
              << std::setw(nameWidth) << potential.first << " "
              << std::setw(nameWidth) << potential.second << " "    
              << std::setw(nameWidth) << FIELDFile::dl_polyType(potential.type);
    // Convert parameters to DL_POLY representation
    std::vector<Real> const dl_parameters =
      FIELDFile::dl_polyParameters(potential.type, potential.parameters);

    for(size_t j = 0; j < dl_parameters.size(); ++j)
    {
      fieldFile << " " << std::setiosflags(std::ios::fixed)
                << dl_parameters[j];
    }
    fieldFile << std::endl;
  }

  fieldFile.precision(precision); // Return to original value

  return fieldFile;
}

FIELDFile &operator>>(FIELDFile &fieldFile, std::vector<PairPotential> &potentials)
{
  assert(fieldFile.is_open() && fieldFile.mode() == IN);

  fieldFile.seekg(0, std::ios::beg); // Rewind file (in case it was read previously)

  std::string line; // Line read from file

  // Search for "units" record (for consistency check)
  do
  {
    std::getline(fieldFile, line);
    if(line.find("units") != line.npos) break;
  }
  while(!fieldFile.eof());

  if(fieldFile.eof())
  {
    std::cerr << "Error reading file " << fieldFile.fileName()
              << ": missing units record" << std::endl;
    return fieldFile;
  }

  // Read units record
  std::stringstream unitStream(line);
  std::string keyword;
  unitStream >> keyword; // "units"
  unitStream >> keyword; // Units keyword
  fieldFile.setUnitsLabel(keyword);

  // Search for "vdw" record
  do
  {
    std::getline(fieldFile, line);
    if(line.find("vdw") != line.npos) break;
  }
  while(!fieldFile.eof());
 
  if(fieldFile.eof())
  {
    std::cerr << "Error reading pair potentials from file " << fieldFile.fileName()
              << ": vdw record not found." << std::endl;
    return fieldFile;
  }

  // Store the number of pair potential records
  std::stringstream headerStream(line);
  size_t numPotentials; // Number of pair potential records

  headerStream >> keyword; // "vdw"
  headerStream >> numPotentials;

  if(headerStream.fail())
  {
    std::cerr << "Error reading number of pair potential records from file " << fieldFile.fileName() << std::endl;
    return fieldFile;
  }

  // Read pair potential records
  for(size_t i = 0; i < numPotentials; ++i)
  {
    std::string atom1, atom2; // Atom types 
    std::string typeString;   // Potential type read from file
    PotentialType type;       // Pair potential type
    Real parameter;           // Pair potential parameter

    std::getline(fieldFile, line);
    std::stringstream potentialRecordStream(line);
    potentialRecordStream >> atom1 >> atom2 >> typeString;
    if(potentialRecordStream.fail())
    {
      std::cerr << "Error reading pair potential record from FIELD file " << fieldFile.fileName() << std::endl;
      return fieldFile;
    }
    std::vector<Real> r_parameters; // Parameters read from file
    while(!potentialRecordStream.eof())
    {
      potentialRecordStream >> parameter;
      if(potentialRecordStream.eof() && potentialRecordStream.fail()) break;
      if(potentialRecordStream.fail())
      {
        std::cerr << "Error reading pair potential parameters from FIELD file " << fieldFile.fileName() << std::endl;
        return fieldFile;
      }
      r_parameters.push_back(parameter);
    }
    type = FIELDFile::potentialType(typeString);
    std::vector<Real> const parameters =
      FIELDFile::potentialParameters(typeString, r_parameters);
    potentials.push_back(PairPotential(atom1, atom2, type, parameters));
  }
  return fieldFile;
}

// Private functions

std::string const FIELDFile::dl_polyType(PotentialType const &type)
{
  switch(type)
  {
    case LENNARD_JONES:
      return "lj";
    case MIE:
      return "nm";
    case BUCKINGHAM:
      return "buck";
    case BORN_HUGGINS_MEYER:
      return "bhm";
    case MORSE:
      return "mors";
    case WCA:
      return "wca";
    case UNKNOWN:
    default:
      std::cerr << "Error in FIELDFile::dl_polyType - unknown potential" 
                << std::endl;
      return std::string();
  }
}

PotentialType const FIELDFile::potentialType(std::string const &dl_polyType)
{
  if((dl_polyType.find("lj") != dl_polyType.npos) || 
     (dl_polyType.find("12-6") != dl_polyType.npos))
    return LENNARD_JONES;
  else if(dl_polyType.find("nm") != dl_polyType.npos)
    return MIE;
  else if(dl_polyType.find("buck") != dl_polyType.npos)
    return BUCKINGHAM;
  else if(dl_polyType.find("bhm") != dl_polyType.npos)
    return BORN_HUGGINS_MEYER;
  else if(dl_polyType.find("mors") != dl_polyType.npos)
    return MORSE;
  else if(dl_polyType.find("wca") != dl_polyType.npos)
    return WCA;
  else return UNKNOWN;
}

std::vector<Real> const 
  FIELDFile::dl_polyParameters(PotentialType const &type, 
                               std::vector<Real> const &parameters)
{
  std::vector<Real> dl_parameters(parameters.size(), 0.0);
  switch(type)
  {
    case LENNARD_JONES:
      // Same format, no change
      assert(parameters.size() == 2);
      dl_parameters[0] = parameters[0];
      dl_parameters[1] = parameters[1];
      break;
    case MIE:
      // Need to convert sigma to r_0 and reorder
      assert(parameters.size() == 4);
      assert(parameters[2] != parameters[3]);
      assert((parameters[2] != 0) && parameters[3] != 0);
      dl_parameters[0] = parameters[0];
      dl_parameters[1] = parameters[2];
      dl_parameters[2] = parameters[3];
      dl_parameters[3] = parameters[1]*
                           pow(parameters[2]/parameters[3], 
                               1.0/(parameters[2] - parameters[3]));
      break;
    case BUCKINGHAM:
      // Need to convert B to rho
      assert(parameters.size() == 3);
      assert(parameters[1] != 0);
      dl_parameters[0] = parameters[0];
      dl_parameters[1] = 1.0/parameters[1];
      dl_parameters[2] = parameters[2];
      break;
    case BORN_HUGGINS_MEYER:
      // Need to reorder
      assert(parameters.size() == 5);
      dl_parameters[0] = parameters[0];
      dl_parameters[1] = parameters[1];
      dl_parameters[2] = parameters[4];
      dl_parameters[3] = parameters[2];
      dl_parameters[4] = parameters[3];
      break;
    case MORSE:
      // Need to reorder
      assert(parameters.size() == 3);
      dl_parameters[0] = parameters[0];
      dl_parameters[1] = parameters[2];
      dl_parameters[2] = parameters[1];
      break;
    case WCA:
      // Same format, no change
      assert(parameters.size() == 2);
      dl_parameters[0] = parameters[0];
      dl_parameters[1] = parameters[1];
      break;
    case UNKNOWN:
    default:
      // Potential is unknown, just copy parameters
      dl_parameters = parameters;
  }
  return dl_parameters;
}
 
std::vector<Real> const 
  FIELDFile::potentialParameters(std::string const &dl_polyType,
                                 std::vector<Real> const &dl_polyParameters)
{
  std::vector<Real> parameters(dl_polyParameters.size(), 0.0);
  if(dl_polyType.find("lj") != dl_polyType.npos)
  {
    // Parameters are in the same format
    if(dl_polyParameters.size() != 2)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": bad number of parameters" << std::endl;
      return parameters;
    }
    parameters = dl_polyParameters;
  }
  else if(dl_polyType.find("12-6") != dl_polyType.npos)
  {
    // Parameters are given for A/r^12 - B/r^6
    if(dl_polyParameters.size() != 2)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": bad number of parameters" << std::endl;
      return parameters;
    }
    if(dl_polyParameters[1] == 0)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": invalid potential parameters" << std::endl;
      return parameters;
    }
    parameters[0] = (dl_polyParameters[1]*dl_polyParameters[1]/(4.0*dl_polyParameters[0]));
    parameters[1] = pow(dl_polyParameters[0]/dl_polyParameters[1], 1.0/6.0);
  }
  else if(dl_polyType.find("nm") != dl_polyType.npos)
  {
    // Parameters are given for E_0/(n-m)*[m*(r0/r)^n - n*(r0/r)^m]
    if(dl_polyParameters.size() != 4)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": bad number of parameters" << std::endl;
      return parameters;
    }
    if((dl_polyParameters[1] == dl_polyParameters[2]) ||
       (dl_polyParameters[1] == 0) ||
       (dl_polyParameters[2] == 0))
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": invalid potential parameters" << std::endl;
      return parameters;
    }
    parameters[0] = dl_polyParameters[0];
    parameters[1] = dl_polyParameters[3]*
                     pow(dl_polyParameters[2]/dl_polyParameters[1],
                         1.0/(dl_polyParameters[1]-dl_polyParameters[2]));
    parameters[2] = dl_polyParameters[1];
    parameters[3] = dl_polyParameters[2];
  }
  else if(dl_polyType.find("buck") != dl_polyType.npos)
  {
    // Parameters are given for A*exp(-r/rho) - C/r^6
    if(dl_polyParameters.size() != 3)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": bad number of parameters" << std::endl;
      return parameters;
    }
    if(dl_polyParameters[1] == 0)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": invalid potential parameters" << std::endl;
      return parameters;
    }
    parameters[0] = dl_polyParameters[0];
    parameters[1] = 1.0/dl_polyParameters[1];
    parameters[2] = dl_polyParameters[2];
  }
  else if(dl_polyType.find("bhm") != dl_polyType.npos)
  {
    // Parameters are shifted (A, B, sigma, C, D instead of A, B, C, D, sigma)
    if(dl_polyParameters.size() != 5)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": bad number of parameters" << std::endl;
      return parameters;
    }
    parameters[0] = dl_polyParameters[0];
    parameters[1] = dl_polyParameters[1];
    parameters[2] = dl_polyParameters[3];
    parameters[3] = dl_polyParameters[4];
    parameters[4] = dl_polyParameters[2];
  }
  else if(dl_polyType.find("mors") != dl_polyType.npos)
  {
    // Parameters are shifted (D, r0, a instead of D, a, r0)
    if(dl_polyParameters.size() != 3)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": bad number of parameters" << std::endl;
      return parameters;
    }
    parameters[0] = dl_polyParameters[0];
    parameters[1] = dl_polyParameters[2];
    parameters[2] = dl_polyParameters[1];
  }
  else if(dl_polyType.find("wca") != dl_polyType.npos)
  {
    // Parameters are in the same format
    if(dl_polyParameters.size() != 2)
    {
      std::cerr << "Error in FIELDFile::potentialParameters"
                << ": bad number of parameters" << std::endl;
      return parameters;
    }
  }
  else // No checks done here, potential is unknown
    parameters = dl_polyParameters;

  return parameters;
}

