/*
** Copyright 2009-2011 Erik Santiso.
** This file is part of pkt2pdb.
** pkt2pdb is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** pkt2pdb is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a pkty of the GNU Lesser General Public License
** along with pkt2pdb. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** pkt2pdb version 0.1
**
** Given a geometry file in pdb format and a pkt file containing order
** parameters, pkt2pdb generates several pdb files containing the
** various order parameters stored in the "temperature factor" fields.
** This is useful for visualization.
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include "common/include/types.h"
#include "mymol/include/file_formats/pdbfile.h"
#include "mymol/include/system.h"
#include "mymol/include/molecule.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/relativeconfiguration.h"
#include "crystdist/include/peaktable.h"
#include "crystdist/include/pktfile.h"

// Some tolerances for comparing internal DOFs
#define DISTANCE_TOL 0.001
#define ANGLE_TOL 0.01

void printHelpMessage()
{
  std::cout << std::endl
            << "pkt2pdb - Visualization of peak table entries." << std::endl << std::endl;
  std::cout << "Usage: pkt2pdb [-help] " << std::endl
            << "               [-pdb pdb_filename] " << std::endl
            << "               [-mmp mmp_filename] " << std::endl
            << "               [-pkt pkt_filename] " << std::endl
            << "               [-out out_prefix] " << std::endl << std::endl
            << "Description: pkt2pdb can be used to generate pdb files containing" << std::endl
            << "             configurations consistent with the entries in the" << std::endl
            << "             the given pkt file. This can be used for visualization." << std::endl
            << "             in the \"temperature factor\" (aka \"beta\") field of" << std::endl
            << "             the ATOM records. This can be used for visualization." << std::endl << std::endl
            << "Options: " << std::endl 
            << "         -help                     - Print this help message" << std::endl
            << "         -pdb pdb_filename         - Input geometry in pdb format" << std::endl
            << "         -mmp mmp_filename         - Molecule map - defines the" << std::endl
            << "                                     relevant degrees of freedom" << std::endl
            << "                                     of the molecule" << std::endl
            << "         -pkt pkt_filename         - Input file containing the peak" << std::endl
            << "                                     table to visualize" << std::endl
            << "         -out out_prefix           - Prefix for the output files" << std::endl 
            << "         The -pdb, -mmp and -pkt options are mandatory." << std::endl
            << "         The input pdb file must contain the geometry of a" << std::endl
            << "         single molecule of the crystal." << std::endl
            << "         If the mmp file contains more than one molecule, only" << std::endl
            << "         first one is used." << std::endl
            << "         If the -out option is not given, the output prefix will" << std::endl
            << "         be the same as the input pdb file name (without the " << std::endl
            << "         .pdb extension)," << std::endl << std::endl;
}

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    // No input options, print help message
    printHelpMessage();
    return 0;
  }
  // Copy argument list to a vector of strings
  std::vector<std::string> arguments;
  for(int i = 0; i < argc; ++i)
    arguments.push_back(std::string(argv[i]));
  // Search for command line options
  bool found_pdb = false;
  bool found_mmp = false;
  bool found_pkt = false;
  bool found_out = false;
  
  std::string pdbFileName, mmpFileName, pktFileName, outPrefix;

  for(size_t i = 0; i < arguments.size(); ++i)
  {
    if(arguments[i].find("-") != arguments[i].npos)
    {
      if(arguments[i] == "-help")
      {
        printHelpMessage();
        return 0;
      } // End of "-help" option
      else if(arguments[i] == "-pdb")
      {
        found_pdb = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of pdb file must be given after -pdb option." << std::endl
                    << "Use pkt2pdb -help for help." << std::endl;
          return -1;
        }
        else pdbFileName = arguments[i + 1];
      } // End of "-pdb" option
      else if(arguments[i] == "-mmp")
      {
        found_mmp = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of mmp file must be given after -mmp option." << std::endl
                    << "Use pkt2pdb -help for help." << std::endl;
          return -1;
        }
        else mmpFileName = arguments[i + 1];
      } // End of "-mmp" option
      else if(arguments[i] == "-pkt")
      {
        found_pkt = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of pkt file must be given after -pkt option." << std::endl
                    << "Use pkt2pdb -help for help." << std::endl;
          return -1;
        }
        else pktFileName = arguments[i + 1];
      } // End of "-pkt" option
      else if(arguments[i] == "-out")
      {
        found_out = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Output prefix must be given after -out option." << std::endl
                    << "Use pkt2pdb -help for help." << std::endl;
          return -1;
        }
        else outPrefix = arguments[i + 1];
      } // End of "-out" option
      else
      {
        std::cerr << "Unknown option: " << arguments[i] << std::endl
                  << "Use pkt2pdb -help for help." << std::endl;
        return -1;
      }
    }
  } // End of loop over arguments
  
  if(!found_pdb || !found_mmp || !found_pkt)
  {
    std::cerr << "Error: Need to specify pdb, mmp and pkt files." << std::endl
              << "Use pkt2pdb -help for help." << std::endl;
    return -1;
  }
  
  if (!found_out) 
  {
    size_t pdbNameSize = pdbFileName.size();
    outPrefix = pdbFileName.substr(0, pdbNameSize - 4);
  }
  
  // Read input files
  System<Molecule> mySystem;
  PDBFile myPDBFile(pdbFileName, IN);
  myPDBFile >> mySystem;

  if((!myPDBFile && !myPDBFile.eof()) ||
     mySystem.size() < 1 ||
     (mySystem.size() > 0 && mySystem[0].numAtoms() == 0))
  {
    std::cerr << "Error reading geometry file " << pdbFileName << std::endl
              << "Use pkt2pdb -help for help." << std::endl;
    return -1;
  }

  if(mySystem.size() > 1)
  {
    std::cerr << "pdb file must contain only one molecule" << std::endl
              << "Use pkt2pdb -help for help." << std::endl;
  }

  std::vector<MoleculeMap> myMoleculeMaps;
  MMPFile myMMPFile(mmpFileName);
  myMMPFile >> myMoleculeMaps;
  if((!myMMPFile && !myMMPFile.eof()) ||
     myMoleculeMaps.size() < 1 ||
     (myMoleculeMaps.size() > 0 && myMoleculeMaps[0].frameAtoms.size() < 1))
  {
    std::cerr << "Error reading molecule map file " << mmpFileName << std::endl
              << "Use getpeaks -help for help." << std::endl;
    return -1;
  }

  PeakTable myPeakTable;
  PKTFile myPKTFile(pktFileName, IN);
  myPKTFile >> myPeakTable;
  if((!myPKTFile && !myPKTFile.eof()) ||
     myPeakTable.numGroups() < 1)
  {
    std::cerr << "Error reading peak table file " << pktFileName << std::endl
              << "Use getpeaks -help for help." << std::endl;
    return -1;
  }

  MoleculeMap const &molMap = myMoleculeMaps[0];
  // Loop over peak groups
  for(size_t iGroup = 0; iGroup < myPeakTable.numGroups(); ++iGroup)
  {
    PeakGroup const &currentGroup = myPeakTable.peaks(iGroup);
    // Loop over peaks
    for(size_t iPeak = 0; iPeak < currentGroup.size(); ++iPeak)
    {
      RelativeConfiguration const &conf = currentGroup[iPeak];
      // Build reference structures
      Molecule &firstMolecule = mySystem[0];
      std::vector<InternalDOF> const &firstInternalDOFs = conf.internalDOFs.first;
      if(firstInternalDOFs.size() > 0)
      {
        // Make internal DOFs consistent
        for(size_t iIntDOF = 0; iIntDOF < firstInternalDOFs.size(); ++iIntDOF)
        {
          InternalDOF const &intDOF = firstInternalDOFs[iIntDOF];
          InternalDOFMap const &intMap = molMap.internalDOFs[iIntDOF];

          if(intDOF.type != intMap.type)
          {
            std::cerr << "Inconsistent data in pkt file and mmp file " << std::endl
                      << "Exiting..." << std::endl;
            return -1;
          }

          switch(intDOF.type)
          {
          case DISTANCE:
          {
            size_t const atom1 = intMap.atoms[0];
            size_t const atom2 = intMap.atoms[1];
            Real currentDistance = firstMolecule.distance(atom1, atom2);
            if(fabs(currentDistance - intDOF.value) > DISTANCE_TOL)
            {
              // Adjust distance
              
            }
          }
            break;
          case ANGLE:
            // TODO: Add code here
            break;
          case DIHEDRAL:
          default:
            // TODO: Add code here
            break;
          } // End of switch for internal DOF type
        } // End of loop over internal DOFs of the first molecule
      } // End of internal degrees of freedom for first molecule

    }
  }
/*
  CrystalOrderParameters myOrderParameters;
  COPFile myCOPFile(pktFileName, IN);
  myCOPFile >> myOrderParameters;

  // Internal degree of freedom test
  if(found_i && intIndx >= myOrderParameters.internal[0].size())
  {
    std::cerr << "Internal degree of freedom index too large: file contains only "
              << myOrderParameters.internal[0].size() << " internal DOFs" << std::endl
              << "Use pkt2pdb -help for help." << std::endl;
    return -1;
  }
  
  std::vector<Real> *opsPtr;  // Pointer to the OPs to write
  if(found_d) opsPtr = &myOrderParameters.distance;
  else if(found_bo) opsPtr = &myOrderParameters.bondOrientation;
  else if(found_ro) opsPtr = &myOrderParameters.relativeOrientation;
  else if(found_i) opsPtr = &myOrderParameters.internal[intIndx];
  
  // If needed, add dummy atoms
  if(found_add)
  {
    // Find maximum residue ID number
    size_t maxResID = 0;
    for(size_t i = 0; i < mySystem.size(); ++i)
    {
      size_t currentID = mySystem[i].atom(0).residueID;
      if(currentID > maxResID) maxResID = currentID;
    }
    
    // Dummy atom properties
    std::string const dummySymbol("H");
    Vector3D const dummyPosition(LARGE_COORDINATE);
    size_t const dummyResID = maxResID + 1;
    std::string const dummyResName("DUM");
    std::string const dummyRoleName(" H0 ");
    
    Molecule dummyResidue; // Dummy residue
    for(size_t i = 0; i < opsToAdd.size(); ++i)
    {
      dummyResidue.addAtom(Atom(dummySymbol, dummyPosition, dummyResID,
                                dummyResName, dummyRoleName, opsToAdd[i]));
      (*opsPtr).push_back(opsToAdd[i]);
    }
    mySystem.add(dummyResidue);
  }

  // If needed, normalize OPs
  if(found_norm)
  {
    Real maxOP = 0.0;
    for(size_t i = 0; i < (*opsPtr).size(); ++i)
      if((*opsPtr)[i] > maxOP) maxOP = (*opsPtr)[i];
    if(maxOP == 0.0)
      std::cerr << "Warning: maximum value of order parameters is zero - will not normalize" << std::endl;
    else
    {
      Real normFactor = maxValue/maxOP;
      for(size_t i = 0; i < (*opsPtr).size(); ++i)
        (*opsPtr)[i] *= normFactor;
      if(found_add)
      {
        Molecule &dummyResidue = mySystem[mySystem.size() - 1];
        for(size_t i = 0; i < dummyResidue.numAtoms(); ++i)
          dummyResidue.setTemperatureFactor(i, normFactor*dummyResidue.atom(i).tempFactor);
      }      
    }
  }
  
  // Write order parameters to pdb file
  size_t const numCells = myOrderParameters.grid.numCells();
  size_t const numResidues = (found_add)?(mySystem.size() - 1):mySystem.size();
  
  for(size_t i = 0; i < numResidues; ++i)
  {
    size_t opIndex;
    Molecule &currentMolecule = mySystem[i];
    if(numCells > 0) // Use grid
      opIndex = myOrderParameters.grid.getIndex(currentMolecule.centerOfMass(), mySystem.lattice());
    else opIndex = i; // Use molecule indices
    if(opIndex >= (*opsPtr).size())
    {
      std::cerr << "Order parameter index out of bounds - " << std::endl
                << "pkt file may not correspond to input pdb file." << std::endl;
      return -1;
    }
    currentMolecule.setTemperatureFactor((*opsPtr)[opIndex]);
  }
  PDBFile outputFile(outFileName, OUT);
  outputFile << mySystem;
*/  
  return 0;
}
