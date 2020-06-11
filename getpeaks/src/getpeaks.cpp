/*
** Copyright 2009-2011 Erik Santiso.
** This file is part of getpeaks.
** getpeaks is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** getpeaks is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with getpeaks. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** getpeaks version 0.1
**
** For a given geometry (or set of geometries), finds peaks and peak groups 
** in the pair distribution function of distances, bond orientations, 
** relative orientations and internal degrees of freedom. This is a wrapper
** for classes GlobalHistogram and PeakTable in the crystdist library.
**
** Notes:
**
** - It would be good (for later) to include a utility to read, trim, and
**   then write global histograms and peak tables.
*/

#include <iostream>
#include <string>
#include <vector>
#include "common/include/types.h"
#include "common/include/iofile.h"
#include "mymol/include/file_formats/pdbfile.h"
#include "mymol/include/system.h"
#include "mymol/include/molecule.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/ptmfile.h"
#include "crystdist/include/globalhistogram.h"
#include "crystdist/include/peaktable.h"
#include "crystdist/include/pktfile.h"
#include "getpeaks/include/getpeaksettings.h"

void printHelpMessage()
{
  std::cout << std::endl
            << "getpeaks - Find peaks and peak groups in the pair distribution function of" << std::endl
            << "           distances, bond orientations, relative orientations and internal" << std::endl
            << "           degrees of freedom of an ordered system." << std::endl << std::endl;
  std::cout << "Usage: getpeaks [-help] " << std::endl
            << "                [-pdb pdb_filename] " << std::endl
            << "                [-mmp mmp_filename] " << std::endl
            << "                [-hsf hsf_filename] " << std::endl
            << "                [-out out_prefix] " << std::endl << std::endl
            << "Description: getpeaks can be used to find peaks in the distribution of" << std::endl 
            << "             distances, bond orientations, relative orientations and" << std::endl
            << "             internal degrees of freedom for a particular species given" << std::endl
            << "             a geometry (or more) and a molecule map." << std::endl << std::endl
            << "Options: " << std::endl
            << "         -help                - Print this help message" << std::endl
            << "         -pdb pdb_filename    - Input geometry in pdb format. If the file" << std::endl
            << "                                contains multiple frames, all are used." << std::endl
            << "         -mmp mmp_filename    - Molecule map - defines the species for" << std::endl
            << "                                which order parameters will be calculated," << std::endl
            << "                                and how. If it contains more than one" << std::endl
            << "                                species, only the first one is read." << std::endl
            << "         -hsf hsf_filename    - Histogram settings file - contains parameters" << std::endl
            << "                                used to construct the global histogram of" << std::endl
            << "                                distances, bond orientations, relative" << std::endl
            << "                                orientations, and internal degrees of" << std::endl
            << "                                freedom, and to construct peak groups." << std::endl
            << "         -out out_prefix      - Prefix for output file names." << std::endl << std::endl
            << "         The -pdb and -mmp options are mandatory." << std::endl
            << "         The pdb file should contain unit cell parameters in a CRYST1" << std::endl
            << "         record - otherwise default cell parameters will be used." << std::endl
            << "         The mmp file format is defined in the crystdist library." << std::endl
            << "         If the histogram settings file name is not given, getpeaks uses" << std::endl
            << "         default values for all the settings. See getpeaksettings.h and" << std::endl
            << "         globalhistogram.h in the crystdist library for details." << std::endl
            << "         If the output file prefix is not given, a default prefix is" << std::endl
            << "         constructed from the pdb file name." << std::endl << std::endl;
}

int main(int argc, char* argv[])
{
  typedef IOFile HSFFile;
  typedef IOFile HistogramFile;

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
  bool found_hsf = false;
  bool found_out = false;
  std::string pdbFileName, mmpFileName, hsfFileName, outPrefix; 
  
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
                    << "Use getpeaks -help for help." << std::endl;
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
                    << "Use getpeaks -help for help." << std::endl;
          return -1;
        }
        else mmpFileName = arguments[i + 1];
      } // End of "-mmp" option
      else if(arguments[i] == "-hsf")
      {
        found_hsf = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of hsf file must be given after -hsf option." << std::endl
                    << "Use getpeaks -help for help." << std::endl;
          return -1;
        }
        else hsfFileName = arguments[i + 1];
      } // End of "-hsf" option
      else if(arguments[i] == "-out")
      {
        found_out = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Output prefix must be given after -out option." << std::endl
                    << "Use getpeaks -help for help." << std::endl;
          return -1;
        }
        else outPrefix = arguments[i + 1];
      } // End of "-out" option
      else
      {
        std::cerr << "Unknown option: " << arguments[i] << std::endl
                  << "Use getpeaks -help for help." << std::endl;
        return -1;
      }
    }
  } // End of loop over arguments

  if(!found_pdb || !found_mmp) 
  {
    std::cerr << "Error: Need to specify pdb and mmp files." << std::endl
              << "Use getpeaks -help for help." << std::endl;
    return -1;
  }
  if(!found_out)
  {
    size_t pdbNameSize = pdbFileName.size();
    outPrefix = pdbFileName.substr(0, pdbNameSize - 4); 
  }

  // Read molecule map and histogram settings files
  
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

  GetpeakSettings mySettings;
  HSFFile myHSFFile(hsfFileName, IN);
  myHSFFile >> mySettings;
  if(!myHSFFile && !myHSFFile.eof())
  {
    std::cerr << "Error reading histogram settings file " << hsfFileName << std::endl
              << "Use getpeaks -help for help." << std::endl;
    return -1;
  } 

  System<Molecule> mySystem;
  PDBFile myPDBFile(pdbFileName, IN);
  std::string ptmFileName = outPrefix + ".ptm"; 
  PTMFile myPTMFile(ptmFileName, OUT);
  
  GlobalHistogram myHistogram;
  myHistogram.changeSettings(mySettings.histogramSettings);

  // Loop over geometry frames
  do
  {
    mySystem.clear();

    // Avoid error message from PDBfile if we are at the end
    if(myPDBFile.peek() == EOF) break; 

    myPDBFile >> mySystem;

    if((!myPDBFile && !myPDBFile.eof()) ||
       mySystem.size() < 1 ||
       (mySystem.size() > 0 && mySystem[0].numAtoms() == 0))
    {
      std::cerr << "Error reading geometry file " << pdbFileName << std::endl
                << "Use getpeaks -help for help." << std::endl;
      return -1;
    }

    if(myPDBFile.eof() && mySystem.size() == 0) break; // Done reading frames
  
    // Convert to System<PointMolecule>
    mySystem.setLatticeType(CRYSTAL); // To get PBC's right
    System<PointMolecule> myPointMolecules = getPointMolecules(mySystem, myMoleculeMaps[0]);

    // Add to ptm file to use later
    myPTMFile << myPointMolecules;

    // Add to global histogram
    myHistogram.addSystem(myPointMolecules, myMoleculeMaps[0].name);
  }
  while(!myPDBFile.eof());
 
  // Done adding frames, trim and write global histogram to file

  myHistogram.reduce(mySettings.normCutoff, mySettings.histFreqCutoff);
  std::string histogramFileName = outPrefix + ".hist"; 
  HistogramFile myHistogramFile(histogramFileName, OUT);
  myHistogramFile << myHistogram;

  // Build peak table

  myPTMFile.close();
  PTMFile myPTMFileAgain(ptmFileName, IN);
  PeakTable myPeakTable;
  do
  {
    System<PointMolecule> myPointMolecules;
   
    // Avoid error message from PTMFile if we are at the end
    if(myPTMFileAgain.peek() == EOF) break;

    myPTMFileAgain >> myPointMolecules; 

    if(!myPTMFileAgain && !myPTMFileAgain.eof())
    {
      std::cerr << "Internal error reading point molecule file " << ptmFileName << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }

    if(myPTMFileAgain.eof() && myPointMolecules.size() == 0) break; // Done reading frames

    // Add to global peak table
    myPeakTable.addSystem(myPointMolecules, myHistogram);
  }
  while(!myPTMFileAgain.eof()); 

  // Done adding frames, trim, sort and write peak table to file

  myPeakTable.reduce(mySettings.groupFreqCutoff);
  myPeakTable.sort();
  std::string pktFileName = outPrefix + ".pkt"; 
  PKTFile myPKTFile(pktFileName, OUT);
  myPKTFile << myPeakTable;

  return 0;
}
