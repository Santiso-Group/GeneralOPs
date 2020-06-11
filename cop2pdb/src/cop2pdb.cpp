/*
** Copyright 2009-2011 Erik Santiso.
** This file is part of cop2pdb.
** cop2pdb is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** cop2pdb is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with cop2pdb. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** cop2pdb version 0.1
**
** Given a geometry file in pdb format and a cop file containing order
** parameters, cop2pdb generates several pdb files containing the
** various order parameters stored in the "temperature factor" fields.
** This is useful for visualization.
**
** Notes:
** 
** - This assumes that the "solute" residues (i.e. the ones for which
**   the OPs have been written in the cop file) are at the beginning
**   of the pdb file.
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
#include "crystdist/include/crystalops.h"
#include "crystdist/include/copfile.h"
#define LARGE_COORDINATE 9999.0

void printHelpMessage()
{
  std::cout << std::endl
            << "cop2pdb - Write crystal order parameters to pdb files." << std::endl << std::endl;
  std::cout << "Usage: cop2pdb [-help] " << std::endl
            << "               [-pdb pdb_filename] " << std::endl
            << "               [-cop cop_filename] " << std::endl
            << "               [-out out_filename] " << std::endl 
            << "               [-n {max_value}] " << std::endl << std::endl
            << "Description: cop2pdb can be used to generate pdb files containing" << std::endl
            << "             the crystal order parameters in the given cop file" << std::endl
            << "             in the \"temperature factor\" (aka \"beta\") field of" << std::endl
            << "             the ATOM records. This can be used for visualization." << std::endl << std::endl
            << "Options: " << std::endl
            << "         -help                     - Print this help message" << std::endl
            << "         -d, -bo, -ro, -i {index}, - Type of order parameter to write" << std::endl
            << "         -t, -ld                     (distance, bond orientation, relative" << std::endl
            << "                                     orientation, index-th internal)" << std::endl
            << "         -pdb pdb_filename         - Input geometry in pdb format" << std::endl
            << "         -cop cop_filename         - Input file containing the values of" << std::endl
            << "                                     the crystal order parameters" << std::endl
            << "         -out out_filename         - Name of the output pdb file" << std::endl
            << "         -n {max_value}            - Normalize the order parameters to the" << std::endl
            << "                                     (optionally) given maximum value" << std::endl
            << "         -add op1 {op2 ...}        - Adds one or more dummy atoms with the" << std::endl
            << "                                     op values given before normalization" << std::endl << std::endl
            << "         The -pdb and -cop options are mandatory." << std::endl
            << "         Exactly one of -d, -bo, -ro, -i, -t, -ld must be given." << std::endl
            << "         If no index is given after the -i option, the first internal" << std::endl
            << "         degree of freedom is written." << std::endl
            << "         If the -out option is not given, the output file name will" << std::endl
            << "         be the same as the input pdb file name, with a suffix" << std::endl 
            << "         indicating the type of order parameter written." << std::endl
            << "         If the -n option is not given, the raw order parameter values" << std::endl
            << "         will be written. Note that they may overflow the temperature" << std::endl
            << "         factor record of the pdb file if they are too large." << std::endl
            << "         If the -n option is given with no arguments, order parameters" << std::endl
            << "         are normalized to be between 0 and 1. Otherwise they are" << std::endl
            << "         normalized to be between 0 and max_value." << std::endl
            << "         The -add option adds one or more atoms with the op values" << std::endl
            << "         given. The dummy atoms have residue name DUM and are placed" << std::endl
            << "         at 9999.0, 9999.0, 9999.0. This is useful to define a scale" << std::endl
            << "         for the order parameters. The normalization is done after"  << std::endl
            << "         adding the dummy atoms." << std::endl << std::endl;
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
  bool found_cop = false;
  bool found_out = false;
  bool found_norm = false;
  bool found_add = false;
  bool found_d = false;
  bool found_bo = false;
  bool found_ro = false;
  bool found_i = false;
  bool found_t = false;
  bool found_ld = false;
  
  size_t intIndx = 0;
  Real maxValue = 1.0;
  std::vector<Real> opsToAdd;
  std::string pdbFileName, copFileName, outFileName;

  for(size_t i = 0; i < arguments.size(); ++i)
  {
    if(arguments[i].find("-") != arguments[i].npos)
    {
      if(arguments[i] == "-help")
      {
        printHelpMessage();
        return 0;
      }
      else if(arguments[i] == "-d") found_d = true;
      else if(arguments[i] == "-bo") found_bo = true;
      else if(arguments[i] == "-ro") found_ro = true;
      else if(arguments[i] == "-i")
      {
        found_i = true;
        if(i + 1 < arguments.size() &&
           arguments[i + 1].find("-") == arguments[i + 1].npos)
        {
          std::stringstream intIndxStream(arguments[i + 1]);
          intIndxStream >> intIndx;
          if(!intIndxStream && !intIndxStream.eof())
          {
            std::cerr << "Index of internal degree of freedom must be a nonnegative integer." << std::endl;
            return -1;
          }
        }
      } // End of "-i" option
      else if(arguments[i] == "-t") found_t = true;
      else if(arguments[i] == "-ld") found_ld = true;
      else if(arguments[i] == "-pdb")
      {
        found_pdb = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of pdb file must be given after -pdb option." << std::endl
                    << "Use cop2pdb -help for help." << std::endl;
          return -1;
        }
        else pdbFileName = arguments[i + 1];
      } // End of "-pdb" option
      else if(arguments[i] == "-cop")
      {
        found_cop = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of cop file must be given after -cop option." << std::endl
                    << "Use cop2pdb -help for help." << std::endl;
          return -1;
        }
        else copFileName = arguments[i + 1];
      } // End of "-cop" option
      else if(arguments[i] == "-out")
      {
        found_out = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Output file name must be given after -out option." << std::endl
                    << "Use cop2pdb -help for help." << std::endl;
          return -1;
        }
        else outFileName = arguments[i + 1];
      } // End of "-out" option
      else if(arguments[i] == "-n")
      {
        found_norm = true;
        if(i + 1 < arguments.size() &&
           arguments[i + 1].find("-") == arguments[i + 1].npos)
        {
          // Read maximum value
          std::stringstream maxValueStream(arguments[i + 1]);
          maxValueStream >> maxValue;
          if((!maxValueStream && !maxValueStream.eof()) || maxValue <= 0.0)
          {
            std::cerr << "Maximum value for normalization must be a positive real number." << std::endl
                      << "Use cop2pdb -help for help." << std::endl;
            return -1;
          }
        }
      } // End of "-n" option
      else if(arguments[i] == "-add")
      {
        found_add = true;
        if(i + 1 >= arguments.size() ||
           (arguments[i + 1].find("-") != arguments[i + 1].npos) &&
           (arguments[i + 1].find_first_not_of("-.0123456789") != arguments[i + 1].npos))
            // (To allow for negative numbers)
        {
          std::cerr << "Order parameter values must follow the -add option." << std::endl
                    << "Use cop2pdb -help for help." << std::endl;
          return -1;
        }
        for(size_t indx = i + 1;
            (indx < arguments.size()) && 
            (arguments[indx].find_first_not_of("-.0123456789") == arguments[indx].npos);
            ++indx)
        {
          std::stringstream opStream(arguments[indx]);
          Real opValue;
          opStream >> opValue;
          if(!opStream && !opStream.eof()) // Removed the condition that they should be >0 (logs)
          {
            std::cerr << "Order parameters must be real numbers." << std::endl
                      << "Use cop2pdb -help for help." << std::endl;
            return -1;
          }
          opsToAdd.push_back(opValue);
        }
        for(size_t indx = 0; indx < opsToAdd.size(); ++indx) ++i;
      } // End of "-add" option
      else
      {
        std::cerr << "Unknown option: " << arguments[i] << std::endl
                  << "Use cop2pdb -help for help." << std::endl;
        return -1;
      }
    }
  } // End of loop over arguments
  
  if(!found_pdb || !found_cop)
  {
    std::cerr << "Error: Need to specify pdb and cop files." << std::endl
              << "Use cop2pdb -help for help." << std::endl;
    return -1;
  }
  
  if((!found_d && !found_bo && !found_ro && !found_i && !found_t && !found_ld) ||
     (found_d && found_bo) || (found_d && found_ro) || (found_d && found_i) ||
     (found_d && found_t) || (found_d && found_ld) || (found_bo && found_ro) ||
     (found_bo && found_i) || (found_bo && found_t) || (found_bo && found_ld) ||
     (found_ro && found_i) || (found_ro && found_t) || (found_ro && found_ld) ||
     (found_i && found_t) || (found_i && found_ld) || (found_t && found_ld)) // I know, crazy
  {
    std::cerr << "Error: Exactly one of the -d, -bo, -ro, or -i options must be given." << std::endl
              << "Use cop2pdb -help for help." << std::endl;
    return -1;
  }
  
  if (!found_out) 
  {
    size_t pdbNameSize = pdbFileName.size();
    outFileName = pdbFileName.substr(0, pdbNameSize - 4);
    if(found_d) outFileName = outFileName + "_d.pdb";
    else if(found_bo) outFileName = outFileName + "_bo.pdb";
    else if(found_ro) outFileName = outFileName + "_ro.pdb";
    else if(found_i)
    {
      std::stringstream outFileStream;
      outFileStream << outFileName << "_i" << std::setfill('0') << std::setw(2) << intIndx << ".pdb";
      outFileName = outFileStream.str();
    }
    else if(found_t) outFileName = outFileName + "_t.pdb";
    else if(found_ld) outFileName = outFileName + "_ld.pdb";
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
              << "Use cop2pdb -help for help." << std::endl;
    return -1;
  }

  CrystalOrderParameters myOrderParameters;
  COPFile myCOPFile(copFileName, IN);
  myCOPFile >> myOrderParameters;

  // Internal degree of freedom test
  if(found_i && intIndx >= myOrderParameters.internal[0].size())
  {
    std::cerr << "Internal degree of freedom index too large: file contains only "
              << myOrderParameters.internal[0].size() << " internal DOFs" << std::endl
              << "Use cop2pdb -help for help." << std::endl;
    return -1;
  }

  // Total OPs test
  if(found_t && myOrderParameters.total.size() < 1)
  {
    std::cerr << "Error: Order parameter file does not contain total order parameters" << std::endl
              << "Use cop2pdb -help for help." << std::endl;
    return -1;
  }

  // Local density test
  if(found_ld && myOrderParameters.localDensity.size() < 1)
  {
    std::cerr << "Error: Order parameter file does not contain local densities" << std::endl
              << "Use cop2pdb -help for help." << std::endl;
    return -1;
  }

  // Get number of solute residues
  size_t numSoluteResidues = mySystem.size();
  for(size_t i = 1; i < mySystem.size(); ++i)
    if(mySystem[i].name() != mySystem[i-1].name())
    {
      numSoluteResidues = i;
      break;
    }
 
  std::vector<Real> internalOPs;
  std::vector<Real> *opsPtr;  // Pointer to the OPs to write
  if(found_d) opsPtr = &myOrderParameters.distance;
  else if(found_bo) opsPtr = &myOrderParameters.bondOrientation;
  else if(found_ro) opsPtr = &myOrderParameters.relativeOrientation;
  else if(found_i) // Copy values - not very efficient
  {
    internalOPs.clear();
    for(size_t i = 0; i < myOrderParameters.internal.size(); ++i)
      internalOPs.push_back(myOrderParameters.internal[i][intIndx]);
    opsPtr = &internalOPs;
  }
  else if(found_t) opsPtr = &myOrderParameters.total;
  else if(found_ld) opsPtr = &myOrderParameters.localDensity;
 
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
        for(size_t i = 0; i < opsToAdd.size(); ++i)
        {
          Molecule &dummyResidue = mySystem[mySystem.size() - 1 - i];
          for(size_t j = 0; j < dummyResidue.numAtoms(); ++j)
            dummyResidue.setTemperatureFactor(j, normFactor*dummyResidue.atom(j).tempFactor);
        }
      }      
    }
  }
  
  // Write order parameters to pdb file
  std::string soluteResName = mySystem[0].name();
  size_t const numCells = myOrderParameters.grid.numCells();
  for(size_t i = 0; i < numSoluteResidues; ++i)
  {
    size_t opIndex;
    Molecule &currentMolecule = mySystem[i];
    if(numCells > 0) // Use grid
      opIndex = myOrderParameters.grid.getIndex(currentMolecule.centerOfMass(), mySystem.lattice());
    else opIndex = i; // Use molecule indices
    if(opIndex >= (*opsPtr).size())
    {
      std::cerr << "Order parameter index out of bounds - " << std::endl
                << "cop file may not correspond to input pdb file." << std::endl;
      return -1;
    }
    currentMolecule.setTemperatureFactor((*opsPtr)[opIndex]);
  }
  PDBFile outputFile(outFileName, OUT);
  outputFile << mySystem;
  
  return 0;
}

