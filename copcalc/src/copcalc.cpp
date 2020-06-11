/*
** Copyright 2009-2011 Erik Santiso.
** This file is part of copcalc.
** copcalc is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** copcalc is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with copcalc. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** copcalc version 0.1
**
** Calculate crystal order parameters for a given geometry - really just
** a command-line wrapper for class opcalculator in the crystdist library.
*/

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "common/include/types.h"
#include "mymath/include/vector3D.h"
#include "mymol/include/file_formats/pdbfile.h"
#include "mymol/include/system.h"
#include "mymol/include/molecule.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/xtpfile.h"
#include "crystdist/include/statparameters.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/crystalops.h"
#include "crystdist/include/opcalculator.h"
#include "crystdist/include/copfile.h"

void printHelpMessage()
{
  std::cout << std::endl
            << "copcalc - Calculate crystal order parameters for a given geometry." << std::endl << std::endl;
  std::cout << "Usage: copcalc [-help] " << std::endl
            << "               [-pdb pdb_filename] " << std::endl
            << "               [-mmp mmp_filename] " << std::endl
            << "               [-xtp xtp_filename] " << std::endl
            << "               [-numcells nx {ny nz}] " << std::endl
            << "               [-origin ox oy oz] " << std::endl
            << "               [-cutoff cutoff] " << std::endl
            << "               [-sfw width]" << std::endl
            << "               [-cop cop_filename] " << std::endl << std::endl
            << "Description: copcalc can be used to compute crystal order parameters" << std::endl
            << "             for a given geometry using a molecule map, a file" << std::endl
            << "             containing the crystal distribution parameters," << std::endl 
            << "             and the definition of the order parameter grid." << std::endl << std::endl
            << "Options: " << std::endl
            << "         -help                - Print this help message" << std::endl
            << "         -pdb pdb_filename    - Input geometry in pdb format" << std::endl
            << "         -mmp mmp_filename    - Molecule map - defines the species for" << std::endl
            << "                                which order parameters will be calculated," << std::endl
            << "                                and how. If it contains more than one" << std::endl
            << "                                species, only the first one is read" << std::endl
            << "         -xtp xtp_filename    - Crystal distribution parameters" << std::endl
            << "         -numcells nx {ny nz} - Number of cells in the directions of each" << std::endl
            << "                                of the unit cell axes" << std::endl
            << "         -origin ox oy oz     - Origin of the cell" << std::endl
            << "         -sfw width           - Switching function width" << std::endl
            << "         -cutoff cutoff       - Distance cutoff" << std::endl
            << "         -cop cop_filename    - Output file containing target values of" << std::endl
            << "                                the crystal order parameters" << std::endl 
            << "         -log                 - If given, the program writes the natural" << std::endl
            << "                                logarithm of the OPs" << std::endl
            << "         -pnorm p_value       - If given, the program uses a p-norm to" << std::endl
            << "                                calculate averages" << std::endl << std::endl
            << "         The -pdb, -mmp and -xtp options are mandatory." << std::endl
            << "         The mmp and xtp file formats are defined in the crystdist library." << std::endl
            << "         If the -numcells option is not given, a single order parameter" << std::endl
            << "         is calculated for the whole cell." << std::endl
            << "         If only nx is given in the -numcells option, ny and nz are" << std::endl
            << "         assumed to be equal to nx." << std::endl
            << "         If the origin is not given, it is assumed to be (0, 0, 0)." << std::endl
            << "         If the switching function width is not given, it is assumed" << std::endl
            << "         to be zero (i.e. no switching function is used)." << std::endl
            << "         If the distance cutoff is not given, it is assumed to be" << std::endl
            << "         infinity (i.e. no cutoff is used)." << std::endl
            << "         If the -cop option is not given, the output file name will be" << std::endl
            << "         the same as the input file name, with extension .cop" << std::endl << std::endl;
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
  bool found_xtp = false;
  bool found_numcells = false;
  bool found_origin = false;
  bool found_sfw = false;
  bool found_cutoff = false;
  bool found_cop = false;
  bool found_log = false;
  bool found_pnorm = false;
  std::string pdbFileName, mmpFileName, xtpFileName, copFileName;
  size_t nx, ny, nz; nx = ny = nz = 1;
  Real ox, oy, oz; ox = oy = oz = 0.0;
  Real switchWidth = 0.0;
  Real cutoff = 0.0;
  Real pValue = 0.0;
  for(size_t i = 0; i < arguments.size(); ++i)
  {
    if(arguments[i].find("-") != arguments[i].npos)
    {
      if(arguments[i] == "-help")
      {
        printHelpMessage();
        return 0;
      }
      else if(arguments[i] == "-pdb")
      {
        found_pdb = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of pdb file must be given after -pdb option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
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
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        else mmpFileName = arguments[i + 1];
      } // End of "-mmp" option
      else if(arguments[i] == "-xtp")
      {
        found_xtp = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of xtp file must be given after -xtp option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        else xtpFileName = arguments[i + 1];
      } // End of "-xtp" option
      else if(arguments[i] == "-numcells")
      {
        found_numcells = true;
        if(i + 1 == arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Number(s) of cells must follow the -numcells option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
        }
        else if(i + 2 >= arguments.size() ||
                arguments[i + 2].find("-") != arguments[i + 2].npos)
        {
          // Only one number given
          std::stringstream numCellStream(arguments[i + 1]);
          numCellStream >> nx;
          if(!numCellStream && !numCellStream.eof())
          {
            std::cerr << "Number of cells must be a positive integer." << std::endl
                      << "Use copcalc -help for help." << std::endl;
            return -1;
          }
          ny = nz = nx;
        }
        else
        {
          if(i + 3 >= arguments.size() ||
             arguments[i + 3].find("-") != arguments[i + 3].npos)
          {
            std::cerr << "If ny is given after the -numcells option, nz must be given as well." << std::endl
                      << "Use copcalc -help for help." << std::endl;
            return -1;
          }
          // Three numbers given
          std::stringstream numCellStream(arguments[i + 1] + " " + arguments[i + 2] + " " + arguments[i + 3]);
          numCellStream >> nx >> ny >> nz;
          if(!numCellStream && !numCellStream.eof())
          {
            std::cerr << "Numbers of cells must be positive integers." << std::endl
                      << "Use copcalc -help for help." << std::endl;
            return -1;
          }
        }
      } // End of "-numcells" option
      else if(arguments[i] == "-origin")
      {
        found_origin = true;
        if(i + 3 >= arguments.size())
        {
          std::cerr << "Unit cell origin must follow the -origin option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        std::stringstream originStream(arguments[i + 1] + " " + arguments[i + 2] + " " + arguments[i + 3]);
        originStream >> ox >> oy >> oz;
        if(!originStream && !originStream.eof())
        {
          std::cerr << "Unit cell origin must follow the -origin option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
      } // End of "-origin" option
      else if(arguments[i] == "-sfw")
      {
        found_sfw = true;
        if(i + 1 >= arguments.size())
        {
          std::cerr << "Switching function width must follow the -sfw option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        std::stringstream switchWidthStream(arguments[i + 1]);
        switchWidthStream >> switchWidth;
        if(!switchWidthStream && !switchWidthStream.eof())
        {
          std::cerr << "Switching function width must follow the -sfw option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
      } // End of "-sfw" option
      else if(arguments[i] == "-cutoff")
      {
        found_cutoff = true;
        if(i + 1 >= arguments.size())
        {
          std::cerr << "Distance cutoff must follow the -cutoff option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        std::stringstream cutoffStream(arguments[i + 1]);
        cutoffStream >> cutoff;
        if(!cutoffStream && !cutoffStream.eof())
        {
          std::cerr << "Distance cutoff must follow the -cutoff option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
      } // End of "-cutoff" option
      else if(arguments[i] == "-cop")
      {
        found_cop = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of cop file must be given after -cop option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        else copFileName = arguments[i + 1];
      } // End of "-cop" option
      else if(arguments[i] == "-log")
        found_log = true;
      else if(arguments[i] == "-pnorm")
      {
        found_pnorm = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Value of p must be given after the -pnorm option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        std::stringstream pnormStream(arguments[i + 1]);
        pnormStream >> pValue;
        if(!pnormStream && !pnormStream.eof())
        {
          std::cerr << "Value of p must be given after the -pnorm option." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
        if(pValue < 1.0)
        {
          std::cerr << "Value of p must be at least 1." << std::endl
                    << "Use copcalc -help for help." << std::endl;
          return -1;
        }
      }
      else
      {
        std::cerr << "Unknown option: " << arguments[i] << std::endl
                  << "Use copcalc -help for help." << std::endl;
        return -1;
      }
    }
  } // End of loop over arguments

  if(!found_pdb || !found_mmp || !found_xtp)
  {
    std::cerr << "Error: Need to specify pdb, mmp and xtp files." << std::endl
              << "Use copcalc -help for help." << std::endl;
    return -1;
  }
  if(!found_cop)
  {
    size_t pdbNameSize = pdbFileName.size();
    copFileName = pdbFileName.substr(0, pdbNameSize - 4) + ".cop";
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
              << "Use copcalc -help for help." << std::endl;
    return -1;
  }

  std::vector<MoleculeMap> myMoleculeMaps;
  MMPFile myMMPFile(mmpFileName);
  myMMPFile >> myMoleculeMaps;
  if((!myMMPFile && !myMMPFile.eof()) ||
     myMoleculeMaps.size() < 1 ||
     (myMoleculeMaps.size() > 0 && myMoleculeMaps[0].framePoints.size() < 1))
  {
    std::cerr << "Error reading molecule map file " << mmpFileName << std::endl
              << "Use copcalc -help for help." << std::endl;
    return -1;
  }

  CrystalDistributionParameters myCOPParameters;
  XTPFile myXTPFile(xtpFileName, IN);
  myXTPFile >> myCOPParameters;
  if((!myXTPFile && !myXTPFile.eof()) ||
     myCOPParameters.concentrationParameters.size() < 1 ||
     myCOPParameters.means.size() < 1 ||
     myCOPParameters.normalizationFactors.size() < 1)
  {
    std::cerr << "Error reading crystal distribution parameters file " << xtpFileName << std::endl
              << "Use copcalc -help for help." << std::endl;
    return -1;
  }

  // Calculate order parameters
  mySystem.setLatticeType(CRYSTAL); // To get PBC's right
  System<PointMolecule> myPointMolecules = getPointMolecules(mySystem, myMoleculeMaps[0]);
  OrderParameterGrid myCOPGrid(nx, ny, nz, Vector3D(ox, oy, oz));
  OrderParameterCalculator myCOPCalculator(myPointMolecules, myCOPParameters, myCOPGrid, switchWidth, cutoff, found_log, pValue);

  // Write order parameters to file
  COPFile myCOPFile(copFileName, OUT);
  myCOPFile << myCOPCalculator.orderParameters();

  return 0;
}
