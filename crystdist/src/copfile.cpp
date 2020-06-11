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
** .cop file format
*/

#include <cmath>
#include "common/include/assert.h"
#include "crystdist/include/copfile.h"

// Constructors

COPFile::COPFile()
:
IOFile()
{}

COPFile::COPFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

COPFile &operator<<(COPFile &copFile, CrystalOrderParameters const &parameters)
{
  // Write order parameters to file
  assert(copFile.is_open() && copFile.mode() == OUT);
  copFile.setf(std::ios::right, std::ios::adjustfield);
  std::streamsize precision = copFile.precision();
  copFile.precision(WRITE_PRECISION);
  copFile << "!SWITCH_WIDTH: " << std::endl;
  copFile << parameters.switchWidth << std::endl;
  copFile << "!CUTOFF: " << std::endl;
  copFile << sqrt(parameters.cutoffSq) << std::endl;
  copFile << "!GRID: " << std::endl;
  copFile << parameters.grid << std::endl;
  copFile << "!NCELLS: " << std::endl;
  copFile << parameters.distance.size() << std::endl;
  copFile << "!DISTANCE_OPS: " << std::endl;
  for(size_t i = 0; i < parameters.distance.size(); ++i)
    copFile << parameters.distance[i] << std::endl;
  copFile << "!BOND_ORIENTATION_OPS: " << std::endl;
  for(size_t i = 0; i < parameters.bondOrientation.size(); ++i)
    copFile << parameters.bondOrientation[i] << std::endl;
  copFile << "!RELATIVE_ORIENTATION_OPS: " << std::endl;
  for(size_t i = 0; i < parameters.relativeOrientation.size(); ++i)
    copFile << parameters.relativeOrientation[i] << std::endl;
  copFile << "!INTERNAL_OPS: " << std::endl;
  for(size_t i = 0; i < parameters.internal.size(); ++i)
  {
    copFile << parameters.internal[i].size() << " ";
    for(size_t j = 0; j < parameters.internal[i].size(); ++j)
      copFile << parameters.internal[i][j] << " ";
    copFile << std::endl;
  }
  copFile << "!LOCAL_DENSITIES: " << std::endl;
  for(size_t i = 0; i < parameters.localDensity.size(); ++i)
    copFile << parameters.localDensity[i] << std::endl;
  copFile << "!TOTAL_OPS: " << std::endl;
  for(size_t i = 0; i < parameters.total.size(); ++i)
    copFile << parameters.total[i] << std::endl;
  copFile.precision(precision); // Return to original value
  return copFile;
}

COPFile &operator>>(COPFile &copFile, CrystalOrderParameters &parameters)
{
  // Read order parameters from file
  parameters.clear();

  std::string r_label;  // Label read from file

  // Search for switching function width label
  do
  {
    std::getline(copFile, r_label);
    if(r_label.find("!SWITCH_WIDTH:") != r_label.npos) break;
  }
  while(!copFile.eof());
  if(copFile.eof())
  {
    std::cerr << "Error reading file " << copFile.fileName() << ": SWITCH_WIDTH record not found." << std::endl;
    return copFile;
  }

  // Get the switching function width

  Real r_switchWidth; // Switching function width read from file
  copFile >> r_switchWidth;
  if(!copFile)
  {
    std::cerr << "Error reading switching function width from file " << copFile.fileName() << std::endl;
    return copFile;
  }
  if(r_switchWidth < 0.0 || r_switchWidth > 0.5)
  {
    std::cerr << "Invalid switching function width: " << r_switchWidth << " in file " << copFile.fileName() << std::endl;
    return copFile;
  }
  parameters.switchWidth = r_switchWidth;
  copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line

  // Search for distance cutoff label
  do
  {
    std::getline(copFile, r_label);
    if(r_label.find("!CUTOFF:") != r_label.npos) break;
  }
  while(!copFile.eof());
  if(copFile.eof())
  {
    std::cerr << "Error reading file " << copFile.fileName() << ": CUTOFF record not found." << std::endl;
    return copFile;
  }
  
  // Get the distance cutoff
  
  Real r_cutoff; // Distance cutoff read from file
  copFile >> r_cutoff;
  if(!copFile || r_cutoff < 0.0)
  {
    std::cerr << "Invalid distance cutoff: " << r_cutoff << " in file " << copFile.fileName() << std::endl;
    return copFile;
  }
  parameters.cutoffSq = r_cutoff*r_cutoff;
  copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  
  // Search for grid label
  do
  {
    std::getline(copFile, r_label);
    if(r_label.find("!GRID:") != r_label.npos) break;
  }
  while(!copFile.eof());
  if(copFile.eof())
  {
    std::cerr << "Error reading file " << copFile.fileName() << ": GRID record not found." << std::endl;
    return copFile;
  }
  
  // Get grid

  OrderParameterGrid r_grid;  // Grid read from file
  copFile >> r_grid;
  if(!copFile)
  {
    std::cerr << "Error reading grid from file " << copFile.fileName() << std::endl;
    return copFile;
  }
  parameters.grid = r_grid;
  copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line

  // Search for number of cells label
  do
  {
    std::getline(copFile, r_label);
    if(r_label.find("!NCELLS:") != r_label.npos) break;
  }
  while(!copFile.eof());
  if(copFile.eof())
  {
    std::cerr << "Error reading file " << copFile.fileName() << ": NCELLS record not found." << std::endl;
    return copFile;
  }

  // Get number of cells

  size_t r_nCells;  //  Number of cells read from file

  copFile >> r_nCells;
  if(!copFile)
  {
    std::cerr << "Error reading number of cells from file " << copFile.fileName() << std::endl;
    return copFile;
  }
  copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line

  // Read distance order parameters
  std::getline(copFile, r_label);
  if(r_label.find("!DISTANCE_OPS:") == r_label.npos)
  {
    std::cerr << "Error reading file " << copFile.fileName() 
              << ": DISTANCE_OPS record not found" << std::endl;
    parameters.clear();
    return copFile;
  }
  parameters.distance.reserve(r_nCells);
  for(size_t i = 0; i < r_nCells; ++i)
  {
    Real r_distanceOP;
    copFile >> r_distanceOP;
    if(!copFile)
    {
      std::cerr << "Error reading distance order parameters from file "
                << copFile.fileName() << std::endl;
      parameters.clear();
      return copFile;
    }
    parameters.distance.push_back(r_distanceOP);
    copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  }
  // Read bond orientation order parameters
  std::getline(copFile, r_label);
  if(r_label.find("!BOND_ORIENTATION_OPS:") == r_label.npos)
  {
    std::cerr << "Error reading file " << copFile.fileName() 
              << ": BOND_ORIENTATION_OPS record not found" << std::endl;
    parameters.clear();
    return copFile;
  }
  parameters.bondOrientation.reserve(r_nCells);
  for(size_t i = 0; i < r_nCells; ++i)
  {
    Real r_bondOrientationOP;
    copFile >> r_bondOrientationOP;
    if(!copFile)
    {
      std::cerr << "Error reading bond orientation order parameters from file "
                << copFile.fileName() << std::endl;
      parameters.clear();
      return copFile;
    }
    parameters.bondOrientation.push_back(r_bondOrientationOP);
    copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  }

  // Read relative orientation order parameters
  std::getline(copFile, r_label);
  if(r_label.find("!RELATIVE_ORIENTATION_OPS:") == r_label.npos)
  {
    std::cerr << "Error reading file " << copFile.fileName() 
              << ": RELATIVE_ORIENTATION_OPS record not found" << std::endl;
    parameters.clear();
    return copFile;
  }
  parameters.relativeOrientation.reserve(r_nCells);
  for(size_t i = 0; i < r_nCells; ++i)
  {
    Real r_relativeOrientationOP;
    copFile >> r_relativeOrientationOP;
    if(!copFile)
    {
      std::cerr << "Error reading relative orientation order parameters from file "
                << copFile.fileName() << std::endl;
      parameters.clear();
      return copFile;
    }
    parameters.relativeOrientation.push_back(r_relativeOrientationOP);
    copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  }

  // Read internal order parameters
  std::getline(copFile, r_label);
  if(r_label.find("!INTERNAL_OPS:") == r_label.npos)
  {
    std::cerr << "Error reading file " << copFile.fileName() 
              << ": INTERNAL_OPS record not found" << std::endl;
    parameters.clear();
    return copFile;
  }
  parameters.internal.reserve(r_nCells);
  for(size_t i = 0; i < r_nCells; ++i)
  {
    std::vector<Real> r_internalOPs;
    size_t r_numInternalOPs;
    copFile >> r_numInternalOPs;
    if(!copFile)
    {
      std::cerr << "Error reading internal order parameters from file "
                << copFile.fileName() << std::endl;
      parameters.clear();
      return copFile;
    }
    r_internalOPs.reserve(r_numInternalOPs);
    for(size_t j = 0; j < r_numInternalOPs; ++j)
    {
      Real r_internalOP;
      copFile >> r_internalOP;
      if(!copFile)
      {
        std::cerr << "Error reading internal order parameters from file "
                  << copFile.fileName() << std::endl;
        parameters.clear();
        return copFile;
      }
      r_internalOPs.push_back(r_internalOP);
    }
    parameters.internal.push_back(r_internalOPs);
    copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  }

  // If present, read local densities
  std::getline(copFile, r_label);
  if(r_label.find("!LOCAL_DENSITIES:") == r_label.npos)
    return copFile; // Old version of cop file - no local densities
  parameters.localDensity.reserve(r_nCells);
  for(size_t i = 0; i < r_nCells; ++i)
  {
    Real r_localDensity;
    copFile >> r_localDensity;
    if(i == 0 && copFile.eof()) return copFile; // Old version of the OPs - no local densities
    if(!copFile)
    {
      std::cerr << "Error reading local densities from file "
                << copFile.fileName() << std::endl;
      parameters.clear();
      return copFile;
    }
    parameters.localDensity.push_back(r_localDensity);
    copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  }

  // If present, read total order parameters
  std::getline(copFile, r_label);
  if(r_label.find("!TOTAL_OPs:") == r_label.npos)
    return copFile; // Old version of cop file - no total OPs
  parameters.total.reserve(r_nCells);
  for(size_t i = 0; i < r_nCells; ++i)
  {
    Real r_totalOP;
    copFile >> r_totalOP;
    if(i == 0 && copFile.eof()) return copFile; // Old version of the OPs - no total OPs
    if(!copFile)
    {
      std::cerr << "Error reading total order parameter from file "
                << copFile.fileName() << std::endl;
      parameters.clear();
      return copFile;
    }
    parameters.total.push_back(r_totalOP);
    copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  }
  return copFile;
}

