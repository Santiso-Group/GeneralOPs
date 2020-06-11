/*
** Copyright 2008-2011 Erik Santiso.
** This file is part of crystdist.
** crystdist is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** crystdist is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with crystdist. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** .mmp file format - used to read molecule map information
*/

#include <sstream>
#include "crystdist/include/mmpfile.h"

// Constructors

MMPFile::MMPFile()
:
IOFile()
{}

MMPFile::MMPFile(std::string const &fileName)
:
IOFile(fileName, IN)
{}

MMPFile &operator>>(MMPFile &mmpFile, std::vector<MoleculeMap> &moleculeMaps)
{
  assert(mmpFile.is_open());

  size_t r_numEntries;                    // Number of molecule maps in file

  mmpFile.ignore(LINE_LENGTH, '\n');  // Skip title line
  mmpFile >> r_numEntries;            // Get the number of entries
  mmpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line

  for(size_t i = 0; i < r_numEntries; i++)  // Loop over entries
  {
    MoleculeMap currentMap;             // Current molecule map
    std::string r_name;                 // Name read from file
    std::string r_type;                 // Type read from file
    std::string r_frameLine;            // Line containing the frame points (for parsing)
    FramePoint r_framePoint;            // Frame point read from file
    size_t r_frameAtom = 0;             // Frame atom read from file
    Real r_frameWeight = 0.0;           // Frame weight read from file
    size_t r_numIntDOFs;                // Number of internal degrees of freedom read from file
    std::string r_intDOFType;           // Internal DOF type read from file
    std::vector<size_t> r_intDOFAtoms;  // Atom indices for internal degree of freedom read from file
    size_t r_symmetryNumber;            // Symmetry number read from file

    mmpFile.ignore(LINE_LENGTH, '\n');  // Skip title line
    mmpFile >> r_name;                  // Get the molecule name
    currentMap.name = r_name;
    mmpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    mmpFile >> r_type;                  // Get the molecule type
    if(r_type.find("LINEAR_SYMMETRIC") != r_type.npos)
      currentMap.type = LINEAR_SYMMETRIC;
    else if(r_type.find("PLANAR_SYMMETRIC") != r_type.npos)
      currentMap.type = PLANAR_SYMMETRIC;
    else if(r_type.find("LINEAR_ASYMMETRIC") != r_type.npos)
      currentMap.type = LINEAR_ASYMMETRIC;
    else if(r_type.find("GENERAL") != r_type.npos)
      currentMap.type = GENERAL;
    else
    {
      std::cerr << "Error reading file " << mmpFile.fileName() << std::endl
                << "Invalid molecule type " << r_type << ": Assuming GENERAL" << std::endl;
      currentMap.type = GENERAL;
    }
    mmpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    // Parse the line containing the frame points
    std::getline(mmpFile, r_frameLine);
    size_t numPoints;
    if((currentMap.type == LINEAR_SYMMETRIC) || (currentMap.type == LINEAR_ASYMMETRIC))
      numPoints = 2;
    else
      numPoints = 3;
    if(r_frameLine.find("O") != r_frameLine.npos)
    {
      // New version of mmp file - use linear combination of atomic positions
      std::string pointLabel;
      for(size_t i = 0; i < numPoints; ++i)
      {
        // Allow for comments
        size_t commentIndex = r_frameLine.find_first_of("!#/");
        if(commentIndex != r_frameLine.npos)
          r_frameLine = r_frameLine.substr(0, commentIndex);
        std::istringstream pointStream(r_frameLine);
        pointStream >> pointLabel; // Read the label (first one must contain "O" or "o")
        r_framePoint.clear();
        while(!pointStream.eof())
        {
          pointStream >> r_frameWeight >> r_frameAtom;
          if(!pointStream && !pointStream.eof())
          { 
            std::cerr << "Wrong format for frame points in file "
                      << mmpFile.fileName() << std::endl;
            return mmpFile;
          }
          if(pointStream.eof()) break;
          r_framePoint.atoms.push_back(r_frameAtom);
          r_framePoint.weights.push_back(r_frameWeight);
        }
        // Normalize the weights
        Real totalWeight = 0.0;
        for(size_t j = 0; j < r_framePoint.weights.size(); ++j)
          totalWeight += r_framePoint.weights[j];
        if(!totalWeight)
        {
          std::cerr << "Error reading map file " << mmpFile.fileName()
                    << ": Null total weight" << std::endl;
          return mmpFile;
        }
        for(size_t j = 0; j < r_framePoint.weights.size(); ++j)
          r_framePoint.weights[j] /= totalWeight;
        if(i < numPoints - 1) std::getline(mmpFile, r_frameLine);
        currentMap.framePoints.push_back(r_framePoint);
      }
    }
    else
    {
      // For compatibility with old version - read only frame atom indices
      std::istringstream frameStream(r_frameLine);
      for(size_t i = 0; i < numPoints; ++i)
      {
        r_framePoint.clear();
        frameStream >> r_frameAtom;
        r_framePoint.atoms.push_back(r_frameAtom);
        r_framePoint.weights.push_back(1.0);
        if((!mmpFile && !mmpFile.eof()) ||
           (mmpFile.eof() && i < numPoints - 1))
        {
          std::cerr << "Error reading frame point data from file " 
                    << mmpFile.fileName() << std::endl;
          return mmpFile;
        }
        currentMap.framePoints.push_back(r_framePoint);
      }
    }
//    mmpFile >> r_frameAtoms[0] >> r_frameAtoms[1] >> r_frameAtoms[2]; // Get the frame atoms
//    currentMap.frameAtoms = r_frameAtoms;
    mmpFile >> r_numIntDOFs;            // Read the number of internal degrees of freedom
    mmpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    for(size_t j = 0; j < r_numIntDOFs; j++)
    {
      InternalDOFMap currentInternalDOF;
      mmpFile >> r_intDOFType;            // Read the type
      mmpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
      if(r_intDOFType.find("DISTANCE") != r_intDOFType.npos)
      {
        currentInternalDOF.type = DISTANCE;
        r_intDOFAtoms.resize(2);
        mmpFile >> r_intDOFAtoms[0] >> r_intDOFAtoms[1];  // Read atom indices
        mmpFile.ignore(LINE_LENGTH, '\n');                // Skip to next line
        currentInternalDOF.atoms = r_intDOFAtoms;
      }
      else if(r_intDOFType.find("ANGLE") != r_intDOFType.npos)
      {
        currentInternalDOF.type = ANGLE;
        r_intDOFAtoms.resize(3);
        mmpFile >> r_intDOFAtoms[0] >> r_intDOFAtoms[1] >> r_intDOFAtoms[2];  // Read atom indices
        mmpFile.ignore(LINE_LENGTH, '\n');                                    // Skip to next line
        currentInternalDOF.atoms = r_intDOFAtoms;
      }
      else if(r_intDOFType.find("DIHEDRAL") != r_intDOFType.npos)
      {
        currentInternalDOF.type = DIHEDRAL;
        r_intDOFAtoms.resize(4);
        mmpFile >> r_intDOFAtoms[0] >> r_intDOFAtoms[1] 
                >> r_intDOFAtoms[2] >> r_intDOFAtoms[3];  // Read atom indices
        mmpFile.ignore(LINE_LENGTH, '\n');                // Skip to next line
        currentInternalDOF.atoms = r_intDOFAtoms;
        mmpFile >> r_symmetryNumber;                      // Read symmetry number
        mmpFile.ignore(LINE_LENGTH, '\n');                // Skip to next line
        currentInternalDOF.symmetryNumber = r_symmetryNumber;
      }
      else
      {
        std::cerr << "Error reading molecule map file " << mmpFile.fileName() << ": " << std::endl
                  << "Unknown internal degree of freedom type: " << r_intDOFType << std::endl;
        return mmpFile;
      }
      currentMap.internalDOFs.push_back(currentInternalDOF);
    }
    if(!mmpFile)
    {
      std::cerr << "Error reading file " << mmpFile.fileName() << ": " << std::endl
                << "Incomplete or invalid molecule map information." << std::endl;
      return mmpFile;
    }
    moleculeMaps.push_back(currentMap);
  }
  return mmpFile;
}
