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
** .xtp file format
*/

#include "crystdist/include/xtpfile.h"

// Constructors

XTPFile::XTPFile()
:
IOFile()
{}

XTPFile::XTPFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

XTPFile &operator<<(XTPFile &xtpFile, CrystalDistributionParameters const &parameters)
{
  // Write parameters to file
  assert(xtpFile.is_open() && xtpFile.mode() == OUT);
  // Check consistency
  size_t numGroups = parameters.means.size();
  if((parameters.concentrationParameters.size() != numGroups) ||
     (parameters.normalizationFactors.size() != numGroups))
  {
    std::cerr << "Internal error writing file " << xtpFile.fileName()
              << ": Inconsistent parameter set." << std::endl;
    return xtpFile;
  }
  size_t numPeaksPerGroup = parameters.means[0].size();
  if((parameters.concentrationParameters[0].size() != numPeaksPerGroup) ||
     (parameters.normalizationFactors[0].size() != numPeaksPerGroup))
  {
    std::cerr << "Internal error writing file " << xtpFile.fileName()
              << ": Inconsistent parameter set." << std::endl;
    return xtpFile;
  }
//  xtpFile.setf(std::ios::scientific, std::ios::floatfield);
  xtpFile.setf(std::ios::right, std::ios::adjustfield);
  std::streamsize precision = xtpFile.precision();
  xtpFile.precision(WRITE_PRECISION);
  xtpFile << "!NGROUPS: " << std::endl;
  xtpFile << numGroups << std::endl;
  for(size_t i = 0; i < numGroups; i++)
  {
    xtpFile << "!GROUP: " << i << std::endl;
    size_t numPeaks = parameters.means[i].size();
    xtpFile << "!PEAKS: " << std::endl;
    xtpFile << numPeaks << std::endl;
    xtpFile << "!MEANS: " << std::endl;
    for(size_t j = 0; j < numPeaks; j++)
      xtpFile << parameters.means[i][j] << std::endl;
    if(parameters.concentrationParameters[i].size() != numPeaks)
      std::cerr << "Warning: Inconsistent data when writing to file "
                << xtpFile.fileName() << std::endl;
    xtpFile << "!CONCENTRATION_PARAMETERS: " << std::endl;
    for(size_t j = 0; j < numPeaks; j++)
      xtpFile << parameters.concentrationParameters[i][j] << std::endl;
    if(parameters.normalizationFactors[i].size() != numPeaks)
      std::cerr << "Warning: Inconsistent data when writing to file "
                << xtpFile.fileName() << std::endl;
    xtpFile << "!NORMALIZATION_FACTORS: " << std::endl;
    for(size_t j = 0; j < numPeaks; j++)
      xtpFile << parameters.normalizationFactors[i][j] << std::endl;
  }
  xtpFile.precision(precision); // Return to original value
  return xtpFile;
}

XTPFile &operator>>(XTPFile &xtpFile, CrystalDistributionParameters &parameters)
{
  // Read parameters from file
  parameters.clear();
  
  std::string r_label;  // Label read from file

  // Search for number of groups label
  do
  {
    std::getline(xtpFile, r_label);
    if(r_label.find("!NGROUPS:") != r_label.npos) break;
  }
  while(!xtpFile.eof());
  if(xtpFile.eof())
  {
    std::cerr << "Error reading file " << xtpFile.fileName() << ": NGROUPS record not found." << std::endl;
    return xtpFile;
  }

  // Read number of groups
  size_t r_nGroups;
  xtpFile >> r_nGroups;
  if(!xtpFile)
  {
    std::cerr << "Error reading number of groups from file " << xtpFile.fileName() << std::endl;
    return xtpFile;
  }
  xtpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line

  // Read parameters
  for(size_t i = 0; i < r_nGroups; i++)
  {
    xtpFile.ignore(LINE_LENGTH, '\n');  // Skip group label
    std::getline(xtpFile, r_label);
    if(r_label.find("!PEAKS") == r_label.npos)
    {
      std::cerr << "Error reading file " << xtpFile.fileName() << ": PEAKS record not found." << std::endl;
      parameters.clear();
      return xtpFile;
    }
    size_t r_nPeaks;  // Number of peaks per group
    xtpFile >> r_nPeaks;
    if(!xtpFile)
    {
      std::cerr << "Error reading number of peaks from file " << xtpFile.fileName() << std::endl;
      parameters.clear();
      return xtpFile;
    }
    xtpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    // Read averages
    std::getline(xtpFile, r_label);
    if(r_label.find("!MEANS") == r_label.npos)
    {
      std::cerr << "Error reading file " << xtpFile.fileName() << ": MEANS record not found." << std::endl;
      parameters.clear();
      return xtpFile;
    }
    PeakGroup meanGroup;
    for(size_t j = 0; j < r_nPeaks; j++)
    {
      RelativeConfiguration r_mean;
      xtpFile >> r_mean;
      xtpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
      if(!xtpFile)
      {
        std::cerr << "Error reading mean peaks from file " << xtpFile.fileName() << std::endl;
        parameters.clear();
        return xtpFile;
      }
      meanGroup.push_back(r_mean);
    }
    parameters.means.push_back(meanGroup);
    // Read concentration parameters
    std::getline(xtpFile, r_label);
    if(r_label.find("!CONCENTRATION_PARAMETERS") == r_label.npos)
    {
      std::cerr << "Error reading file " << xtpFile.fileName() 
                << ": CONCENTRATION_PARAMETERS record not found." << std::endl;
      parameters.clear();
      return xtpFile;
    }
    GroupParameters concentrationParameters;
    for(size_t j = 0; j < r_nPeaks; j++)
    {
      DistributionParameters r_concentrationParameter;
      xtpFile >> r_concentrationParameter;
      xtpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
      if(!xtpFile)
      {
        std::cerr << "Error reading concentration parameters from file " 
                  << xtpFile.fileName() << std::endl;
        parameters.clear();
        return xtpFile;
      }
      concentrationParameters.push_back(r_concentrationParameter);
    }
    parameters.concentrationParameters.push_back(concentrationParameters);
    // Read normalization factors
    std::getline(xtpFile, r_label);
    if(r_label.find("!NORMALIZATION_FACTORS") == r_label.npos)
    {
      std::cerr << "Error reading file " << xtpFile.fileName()
                << ": NORMALIZATION_FACTORS record not found." << std::endl;
      parameters.clear();
      return xtpFile;
    }
    GroupParameters normalizationFactors;
    for(size_t j = 0; j < r_nPeaks; j++)
    {
      DistributionParameters r_normalizationFactor;
      xtpFile >> r_normalizationFactor;
      xtpFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
      if(!xtpFile)
      {
        std::cerr << "Error reading normalization factors from file " 
                  << xtpFile.fileName() << std::endl;
        parameters.clear();
        return xtpFile;
      }
      normalizationFactors.push_back(r_normalizationFactor);
    }
    parameters.normalizationFactors.push_back(normalizationFactors);
  }
  return xtpFile;
}
