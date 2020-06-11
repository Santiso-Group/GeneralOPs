/*
** Copyright 2009-2011 Erik Santiso.
** This file is part of peakdist.
** peakdist is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** peakdist is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with peakdist. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** peakdist version 0.1
**
** Calculate the parameters appearing in the crystal order parameters by
** using a clustering algorithm on the relative configurations. This
** works better for partially or loosely ordered systems, but it does not
** separate peaks by groups.
*/

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "common/include/types.h"
#include "mymol/include/file_formats/pdbfile.h"
#include "mymol/include/file_formats/psffile.h"
#include "mymol/include/file_formats/dcdfile.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/peakclustering.h"
#include "crystdist/include/xtpfile.h"

void printHelpMessage()
{
  std::cout << std::endl
            << "peakdist - Calculate means and concentration parameters for the pair" << std::endl
            << "           distribution function of distances, bond orientations," << std::endl
            << "           relative orientations and internal degrees of freedom of" << std::endl
            << "           a partially or loosely ordered system." << std::endl << std::endl;
  std::cout << "Usage: peakdist [-help] " << std::endl
            << "                [-psf psf_filename] " << std::endl
            << "                [-dcd dcd_filename] " << std::endl
            << "                [-pdb pdb_filename] " << std::endl
            << "                [-mmp mmp_filename] " << std::endl
            << "                [-xtp xtp_filename] " << std::endl << std::endl
            << "Description: peakdist can be used to compute means and concentration" << std::endl
            << "             parameters in the pair distribution function of distances," << std::endl
            << "             bond orientations, relative orientations and internal" << std::endl
            << "             degrees of freedom for a particular species given a" << std::endl
            << "             structure file in psf or pdb format, a dcd trajectory file," << std::endl
            << "             and a molecule map. Unlike cryststat, peakdist does not" << std::endl
            << "             separate peaks into groups." << std::endl << std::endl
            << "Options: " << std::endl
            << "         -help                - Print this help message" << std::endl
            << "         -psf psf_filename    - Structure file in X-PLOR psf format" << std::endl 
            << "         -dcd dcd_filename    - Trajectory file in dcd format" << std::endl
            << "         -pdb pdb_filename    - Coordinates in pdb format (can be used" << std::endl
            << "                                instead of the psf file above" << std::endl
            << "         -mmp mmp_filename    - Molecule map - defines the species for" << std::endl
            << "                                which order parameters will be calculated," << std::endl
            << "                                and how. If it contains more than one" << std::endl
            << "                                species, only the first one is read" << std::endl
            << "         -cutoff cutoff       - Specify the distance cutoff for center-of-" << std::endl
            << "                                mass distances" << std::endl
            << "         -numpeaks num_peaks  - Specify the number of peaks to group" << std::endl
            << "                                the data into" << std::endl
            << "         -numtries num_tries  - Specify the number of clustering attempts" << std::endl
            << "         -verbose             - Print information to cout on the progress" << std::endl
            << "                                of the clustering algorithm" << std::endl
            << "         -xtp xtp_filename    - Output file containing the crystal" << std::endl
            << "                                distribution parameters." << std::endl << std::endl
            << "         The -dcd, -mmp, -numpeaks and one of the -psf or -pdb options are mandatory." << std::endl
            << "         If the -cutoff option is not given, all pairs of molecules are" << std::endl
            << "         used. This may cause the program to crash if the system is large." << std::endl
            << "         If the -numtries option is not given, 50 trials will be done." << std::endl
            << "         The mmp and xtp formats are defined in the crystdist library." << std::endl
            << "         If the -xtp option is not given, the output file name will be" << std::endl
            << "         the same as the psf/pdb file name, with extension .xtp" << std::endl << std::endl;
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
  bool found_psf = false;
  bool found_pdb = false;
  bool found_dcd = false;
  bool found_mmp = false;
  bool found_cutoff = false;
  bool found_numpeaks = false;
  bool found_numtries = false;
  bool found_verbose = false;
  bool found_xtp = false;
  size_t numPeaks, numTries = 50;
  Real cutoff;
  std::string psfFileName, pdbFileName, dcdFileName, mmpFileName, xtpFileName;
  for(size_t i = 0; i < arguments.size(); ++i)
  {
    if(arguments[i].find("-") != arguments[i].npos)
    {
      if(arguments[i] == "-help")
      {
        printHelpMessage();
        return 0;
      }
      else if(arguments[i] == "-psf")
      {
        found_psf = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of psf file must be given after -psf option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else psfFileName = arguments[i + 1];
      } // End of "-psf" option
      else if(arguments[i] == "-dcd")
      {
        found_dcd = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of dcd file must be given after -dcd option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else dcdFileName = arguments[i + 1];
      } // End of "-dcd" option
      else if(arguments[i] == "-pdb")
      {
        found_pdb = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of pdb file must be given after -pdb option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else pdbFileName = arguments[i + 1];
      }
      else if(arguments[i] == "-mmp")
      {
        found_mmp = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of mmp file must be given after -mmp option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else mmpFileName = arguments[i + 1];
      } // End of "-mmp" option
      else if (arguments[i] == "-cutoff")
      {
        found_cutoff = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Distance cutoff must be given after -num option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else
        {
          std::stringstream cutoffStream(arguments[i + 1]);
          cutoffStream >> cutoff;
          if((!cutoffStream && !cutoffStream.eof()) ||
             (cutoff <= 0.0))
          {
            std::cerr << "Distance cutoff must be a positive real number." << std::endl
                      << "Use peakdist -help for help." << std::endl;
            return -1;
          }
        }
      } // End of "-cutoff" option
      else if (arguments[i] == "-numpeaks")
      {
        found_numpeaks = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Number of peaks must be given after the -num option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else
        {
          std::stringstream numPeakStream(arguments[i + 1]);
          numPeakStream >> numPeaks;
          if((!numPeakStream && !numPeakStream.eof())
             || (numPeaks == 0))
          {
            std::cerr << "Number of peaks must be a positive integer." << std::endl
                      << "Use peakdist -help for help." << std::endl;
            return -1;
          }
        }
      } // End of "-numpeaks" option
      else if(arguments[i] == "-numtries")
      {
        found_numtries = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Number of tries must be given after the -num option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else
        {
          std::stringstream numTriesStream(arguments[i + 1]);
          numTriesStream >> numTries;
          if((!numTriesStream && !numTriesStream.eof())
             || (numTries == 0))
          {
            std::cerr << "Number of tries must be a positive integer." << std::endl
                      << "Use peakdist -help for help." << std::endl;
            return -1;
          }
        }
      } // End of "-numtries" option
      else if(arguments[i] == "-verbose")
        found_verbose = true;
      else if(arguments[i] == "-xtp")
      {
        found_xtp = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of xtp output file must be given after -xtp option." << std::endl
                    << "Use peakdist -help for help." << std::endl;
          return -1;
        }
        else xtpFileName = arguments[i + 1];
      } // End of "-xtp" option
      else
      {
        std::cerr << "Unknown option: " << arguments[i] << std::endl
                  << "Use peakdist -help for help." << std::endl;
        return -1;
      }
    }
  } // End of loop over arguments

  if((!found_psf && !found_pdb) || !found_dcd || !found_mmp || !found_numpeaks)
  {
    std::cerr << "Error: Need to specify number of clusters, dcd, mmp and either pdb or psf files." << std::endl
              << "Use peakdist -help for help." << std::endl;
    return -1;
  }
  if(found_psf && found_pdb)
  {
    std::cerr << "Error: Either a psf or a pdb file must be specified, not both." << std::endl
              << "Use peakdist -help for help." << std::endl;
    return -1;
  }
  if(!found_xtp)
  {
    if(found_psf)
    {
      size_t psfNameSize = psfFileName.size();
      xtpFileName = psfFileName.substr(0, psfNameSize - 4) + ".xtp";
    }
    else
    {
      size_t pdbNameSize = pdbFileName.size();
      xtpFileName = pdbFileName.substr(0, pdbNameSize - 4) + ".xtp";
    }
  }

  // Read molecule map file
  
  std::vector<MoleculeMap> myMoleculeMaps;
  MMPFile myMMPFile(mmpFileName);
  myMMPFile >> myMoleculeMaps;
  if((!myMMPFile && !myMMPFile.eof()) ||
     myMoleculeMaps.size() < 1 ||
     (myMoleculeMaps.size() > 0 && myMoleculeMaps[0].framePoints.size() < 1))
  {
    std::cerr << "Error reading molecule map file " << mmpFileName << std::endl
              << "Use peakdist -help for help." << std::endl;
    return -1;
  }
  
  // Read structure
  
  System<Molecule> mySystem;
  if(found_psf)
  {
    PSFFile myPSFFile(psfFileName);
    myPSFFile >> mySystem;
     
    if((!myPSFFile && !myPSFFile.eof()) ||
       mySystem.size() < 1)
    {
      std::cerr << "Error reading psf file " << psfFileName << std::endl
                << "Use peakdist -help for help." << std::endl;
      return -1;
    }
  }
  else
  {
    PDBFile myPDBFile(pdbFileName, IN);
    myPDBFile >> mySystem;
    if((!myPDBFile && !myPDBFile.eof()) ||
       mySystem.size() < 1)
    {
      std::cerr << "Error reading pdb file " << pdbFileName << std::endl
                << "Use peakdist -help for help." << std::endl;
      return -1;
    }
  }

  // Initialize dcd, point molecule objects, peak clustering

  DCDFile myDCDFile(dcdFileName);

  if(!myDCDFile || myDCDFile.numFrames() < 1)
  {
    std::cerr << "Error opening dcd file " << dcdFileName << std::endl
              << "Use peakdist -help for help." << std::endl;
    return -1;
  }

  System<PointMolecule> myPointMolecules;
  PeakClustering myClustering;
  if(found_cutoff) myClustering.setCutoff(cutoff);
  myClustering.setNumClusters(numPeaks);
  myClustering.setNumTrials(numTries);
  myClustering.verboseOutput(found_verbose);
  mySystem.setLatticeType(CRYSTAL); // To get PBC's right  

  // Loop over frames to get data points
  while(myDCDFile.numFramesRead() < myDCDFile.numFrames()) 
  {
    // Read frame
    myDCDFile >> mySystem;

    if(!myDCDFile && !myDCDFile.eof())
    {
      std::cerr << "Error reading dcd file " << dcdFileName << std::endl
                << "Use peakdist -help for help." << std::endl;
      return -1;
    }

    // Convert to point molecules
    myPointMolecules.clear();
    myPointMolecules = getPointMolecules(mySystem, myMoleculeMaps[0]);
   
    // Accumulate statistics
    myClustering.addFrame(myPointMolecules); 
  } // End of loop over dcd file frames

  // Do the clustering and write results to output file
  XTPFile myXTPFile(xtpFileName, OUT);
  myXTPFile << myClustering.parameters();
  
  std::cout << "Davies-Bouldin index for best clustering: " << myClustering.dbIndex() << std::endl;

  return 0;
}