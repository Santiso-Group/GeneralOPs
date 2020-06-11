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
** Crystal statistical parameters containers
**
** These are simple containers for the peak means, concentration
** parameters, and normalization factors for a particular crystal
*/

#ifndef H_STAT_PARAMETERS
#define H_STAT_PARAMETERS

#include "crystdist/include/relativeconfiguration.h"

typedef std::vector<RelativeConfiguration> PeakGroup;

struct DistributionParameters // To store variances, normalization factors
{
  Real distance;
  Real bondOrientation;
  Real relativeOrientation;
  std::pair<std::vector<Real>, std::vector<Real> > internalDOFs;

  DistributionParameters()
  :
  distance(0.0), bondOrientation(0.0), relativeOrientation(0.0), internalDOFs()
  {}

  void clear()
  {
    distance = bondOrientation = relativeOrientation = 0.0;
    internalDOFs.first.clear();
    internalDOFs.second.clear();
  }

  // Distribution parameters output
  friend std::ostream &operator<<(std::ostream &outStream, 
                                  DistributionParameters const &parameters)
  {
    outStream << parameters.distance << " "
              << parameters.bondOrientation << " "
              << parameters.relativeOrientation << " ";
    outStream << parameters.internalDOFs.first.size() << " ";
    for(size_t i = 0; i < parameters.internalDOFs.first.size(); i++)
      outStream << parameters.internalDOFs.first[i] << " ";
    outStream << parameters.internalDOFs.second.size() << " ";
    for(size_t i = 0; i < parameters.internalDOFs.second.size(); i++)
      outStream << parameters.internalDOFs.second[i] << " ";
    return outStream;
  }

  // Distribution parameters input
  friend std::istream &operator>>(std::istream &inStream,
                                  DistributionParameters &parameters)
  {
    parameters.clear();
    inStream >> parameters.distance
             >> parameters.bondOrientation
             >> parameters.relativeOrientation;
    size_t numInternalDOFs;
    inStream >> numInternalDOFs;  // First
    for(size_t i = 0; i < numInternalDOFs; i++)
    {
      Real r_parameter; // Internal DOF parameter read from file
      inStream >> r_parameter;
      parameters.internalDOFs.first.push_back(r_parameter);
    }
    inStream >> numInternalDOFs;  // Second
    for(size_t i = 0; i < numInternalDOFs; i++)
    {
      Real r_parameter; // Internal DOF parameter read from file
      inStream >> r_parameter;
      parameters.internalDOFs.second.push_back(r_parameter);
    }
    return inStream;
  }
};

typedef std::vector<DistributionParameters> GroupParameters;

struct CrystalDistributionParameters
{
  std::vector<PeakGroup> means;                         // Average peaks
  std::vector<GroupParameters> concentrationParameters; // Concentration parameters
  std::vector<GroupParameters> normalizationFactors;    // Normalization factors

  void clear()
  {
    means.clear();
    concentrationParameters.clear();
    normalizationFactors.clear();
  }

  // Statistical parameters output
  friend std::ostream &operator<<(std::ostream &outStream,
                                  CrystalDistributionParameters const &parameters)
  {
    outStream << "Averages: " << std::endl;
    for(size_t i = 0; i < parameters.means.size(); i++)
    {
      outStream << "Group: " << i << std::endl;
      for(size_t j = 0; j < parameters.means[i].size(); j++)
        outStream << parameters.means[i][j] << std::endl;
    }
    outStream << "Concentration parameters: " << std::endl;
    for(size_t i = 0; i < parameters.concentrationParameters.size(); i++)
    {
      outStream << "Group: " << i << std::endl;
      for(size_t j = 0; j < parameters.concentrationParameters[i].size(); j++)
        outStream << parameters.concentrationParameters[i][j] << std::endl;
    }
    outStream << "Normalization factors: " << std::endl;
    for(size_t i = 0; i < parameters.normalizationFactors.size(); i++)
    {
      outStream << "Group: " << i << std::endl;
      for(size_t j = 0; j < parameters.normalizationFactors[i].size(); j++)
        outStream << parameters.normalizationFactors[i][j] << std::endl;
    }
    return outStream;
  }
};

#endif

