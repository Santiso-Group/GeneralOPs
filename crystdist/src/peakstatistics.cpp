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
** Peak statistics class
*/

#include "mymath/include/matrixND.h"
#include "mymath/include/sorting.h"
#include "mymath/include/hungarian.h"
#include "crystdist/include/peakstatistics.h"

// Constructors

PeakStatistics::PeakStatistics()
:
numPoints_(0), peakAccumulators_()
{}

// Interface

void PeakStatistics::addPoint(System<PointMolecule> const &system, PeakTable const &table)
{
  size_t const numGroups = table.numGroups();
  size_t const numPeaksPerGroup = table.peaks(0).size();

  if(numPoints_ == 0)
  {
    // Copy peak table information
    peakAccumulators_.resize(numGroups);
    for(size_t i = 0; i < numGroups; i++)
    {
      if(table.peaks(i).size() != numPeaksPerGroup)
      {
        std::cerr << "Error in PeakStatistics::addPoint: Variable peaks per group not implemented";
        clear();
        return;
      }
      peakAccumulators_[i].resize(numPeaksPerGroup);
    }
  }
  else  // Make sure data is consistent
  {
    if(table.numGroups() != peakAccumulators_.size())
    {
      std::cerr << "Error in PeakStatistics::addPoint: Inconsistent number of groups in peak table";
      return;
    }
    for(size_t i = 0; i < table.numGroups(); i++)
      if(table.peaks(i).size() != peakAccumulators_[i].size())
      {
        std::cerr << "Error in PeakStatistics::addPoint: Inconsistent group data in peak table";
        return;
      }
  }

  // Main loop over molecules
  for(size_t i = 0; i < system.size(); i++)
  {
    // Find neighbors
    std::vector<Real> distancesSq;
    for(size_t j = 0; j < system.size(); j++)
    {
      if(i == j) distancesSq.push_back(REAL_VERYBIG);
      else
        distancesSq.push_back(norm2(system.lattice().difference(system[i].position, system[j].position)));
    }
    // Sort by ascending distance
    std::vector<size_t> neighbors = introSort<std::vector<Real>, Real>(distancesSq);
    // Assign neighbors to peaks
    std::vector<MatrixND> costMatrices(numGroups, MatrixND(numPeaksPerGroup, 0.0)); // One cost matrix per group
    std::vector<RelativeConfiguration> relativeConfigurations;
    // Calculate cost matrices
    for(size_t peak = 0; peak < numPeaksPerGroup; peak++)
    {
      size_t j = neighbors[peak];
      relativeConfigurations.push_back(RelativeConfiguration(system[i], system[j], system.lattice()));
      for(size_t group = 0; group < numGroups; group++)
        for(size_t tablePeak = 0; tablePeak < numPeaksPerGroup; tablePeak++)
          costMatrices[group](peak, tablePeak) = squareDeviation(relativeConfigurations[peak],
                                                                 table.peaks(group)[tablePeak]);
    }
    // Find the right group
    std::vector<std::vector<size_t> > matches;
    std::vector<Real> totalCosts;
    for(size_t group = 0; group < numGroups; group++)
    {
      matches.push_back(hungarianMatch(costMatrices[group]));
      Real totalCost = 0.0;
      for(size_t peak = 0; peak < numPeaksPerGroup; peak++)
        totalCost += costMatrices[group](peak, matches[group][peak]);
      totalCosts.push_back(totalCost);
    }
    size_t rightGroup = 0;
    Real minCost = totalCosts[0];
    for(size_t group = 1; group < numGroups; group++)
      if(totalCosts[group] < minCost)
      {
        rightGroup = group;
        minCost = totalCosts[group];
      }
    // Update accumulators
    for(size_t peak = 0; peak < numPeaksPerGroup; peak++)
      peakAccumulators_[rightGroup][matches[rightGroup][peak]].addValue(relativeConfigurations[peak]);
  } // End of main loop over molecules
  ++numPoints_;
}

// Accessors

CrystalDistributionParameters const PeakStatistics::parameters()
{
  CrystalDistributionParameters parameters;
  for(size_t i = 0; i < peakAccumulators_.size(); i++)
  {
    PeakGroup meanGroup;
    for(size_t j = 0; j < peakAccumulators_[i].size(); j++)
      meanGroup.push_back(peakAccumulators_[i][j].mean());
    parameters.means.push_back(meanGroup);

    GroupParameters currentGroupParameters;
    for(size_t j = 0; j < peakAccumulators_[i].size(); j++)
      currentGroupParameters.push_back(peakAccumulators_[i][j].concentrationParameters());
    parameters.concentrationParameters.push_back(currentGroupParameters);

    currentGroupParameters.clear();
    for(size_t j = 0; j < peakAccumulators_[i].size(); j++)
      currentGroupParameters.push_back(peakAccumulators_[i][j].normalizationFactors());
    parameters.normalizationFactors.push_back(currentGroupParameters);
  }
  return parameters;
}

std::ostream &operator<<(std::ostream &outStream, PeakStatistics &peakStats)
{
  outStream << "Number of points: " << peakStats.numPoints() << std::endl;
  outStream << peakStats.parameters() << std::endl;
  return outStream;
}
