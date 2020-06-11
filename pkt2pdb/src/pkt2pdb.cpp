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
**
** Notes:
**
** - The linear cases haven't been debugged yet.
** - Internal degrees of freedom have not been implemented yet. The problem
**   is, it is not trivial to implement them for the case of cyclic or
**   long molecules. Probably the best solution is to use a pseudo-force
**   field to generate a reasonable structure, as it is done in form2geom.
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "common/include/types.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/quaternion.h"
#include "mymol/include/file_formats/pdbfile.h"
#include "mymol/include/system.h"
#include "mymol/include/molecule.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/relativeconfiguration.h"
#include "crystdist/include/peaktable.h"
#include "crystdist/include/pktfile.h"

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
            << "         .pdb extension)." << std::endl  
            << "         One output file is generated per group and peak in the" << std::endl
            << "         peak table. Each contains two molecules in the" << std::endl
            << "         corresponding relative configuration." << std::endl << std::endl;
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
     (myMoleculeMaps.size() > 0 && myMoleculeMaps[0].framePoints.size() < 1))
  {
    std::cerr << "Error reading molecule map file " << mmpFileName << std::endl
              << "Use getpeaks -help for help." << std::endl;
    return -1;
  }

  PeakTable myPeakTable;
  PKTFile myPKTFile(pktFileName, IN);
  myPKTFile >> myPeakTable;
  if((!myPKTFile && !myPKTFile.eof()) ||
     myPeakTable.numGroups() < 1 ||
     myPeakTable.peaks(0).size() < 1)
  {
    std::cerr << "Error reading peak table file " << pktFileName << std::endl
              << "Use getpeaks -help for help." << std::endl;
    return -1;
  }

  MoleculeMap const &map = myMoleculeMaps[0];
  if(map.internalDOFs.size() > 0)
  {
    std::cerr << "Error: Internal degrees of freedom not implemented yet" << std::endl
              << "Exiting..." << std::endl;
    return -1;
  }

  // Prepare reference configuration
  Molecule refMolecule = mySystem[0];
  PointMolecule pointMolecule(refMolecule, map);
  if(map.type == GENERAL)
    refMolecule.rotate(~pointMolecule.orientation, pointMolecule.position);
  else
  {
    Vector3D const currentAxis = vector(pointMolecule.orientation);
    Vector3D const zAxis(0.0, 0.0, 1.0);
    Vector3D const rotAxis = cross(currentAxis, zAxis);
    Real cosAngle = currentAxis.z/norm(currentAxis);
    if(cosAngle > 1.0) cosAngle = 1.0;    // Prevent acos...
    if(cosAngle < -1.0) cosAngle = -1.0;  // ... from puking
    Real const angle = acos(cosAngle);
    Quaternion rotation(rotAxis, angle);
    refMolecule.rotate(rotation, pointMolecule.position);
  }
  refMolecule.translate(-pointMolecule.position);

  // Width of number fields for file names
  size_t groupWidth = 1 + (size_t)floor(log10((Real)myPeakTable.numGroups()));
  size_t peakWidth = 1 + (size_t)floor(log10((Real)myPeakTable.peaks(0).size()));

  // Loop over peak groups
  for(size_t iGroup = 0; iGroup < myPeakTable.numGroups(); ++iGroup)
  {
    PeakGroup const &currentGroup = myPeakTable.peaks(iGroup);
    // Loop over peaks
    for(size_t iPeak = 0; iPeak < currentGroup.size(); ++iPeak)
    {
      RelativeConfiguration const &conf = currentGroup[iPeak];
      // Build second molecule
      Molecule neighbor = refMolecule;
      switch(map.type)
      {
        case GENERAL:
        {
          Vector3D displacement = conf.distance*vector(conf.bondOrientation);
          Quaternion rotation = conf.relativeOrientation;
          neighbor.rotate(conf.relativeOrientation, 0.0);
          neighbor.translate(displacement);
        }
        break;
        case LINEAR_SYMMETRIC:
        case LINEAR_ASYMMETRIC:
        case PLANAR_SYMMETRIC:
        default: // Should all be the same?
        {
          Real bondAngle = scalar(conf.bondOrientation);
          Real relAngle = scalar(conf.relativeOrientation);
          Vector3D displacement = conf.distance*Vector3D(0.0, sin(bondAngle), cos(bondAngle));
          Quaternion rotation(Vector3D(1.0, 0.0, 0.0), relAngle);
          neighbor.rotate(rotation, 0.0);
          neighbor.translate(displacement);
        }
        break;
      }
      // Write neighbor to file
      std::stringstream fileNameStream;
      fileNameStream << outPrefix;
      fileNameStream << "_" << std::setw(groupWidth) << std::setfill('0') << iGroup;
      fileNameStream << "_" << std::setw(peakWidth) << std::setfill('0') << iPeak;
      PDBFile confPDBFile(fileNameStream.str()+".pdb", OUT);
      System<Molecule> outSystem;
      outSystem.add(refMolecule);
      outSystem.add(neighbor);
      confPDBFile << outSystem;
      confPDBFile.close();
    } // End of loop over peaks
  } // End of loop over groups

  return 0;
}
