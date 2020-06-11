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
** Peak clustering class
*/

#include "crystdist/include/peakclustering.h"
#include <time.h> // For random seed

// Constructors

PeakClustering::PeakClustering()
:
knowClustering_(false), verbose_(false), numClusters_(DEFAULT_NUMCLUSTERS),
numFrames_(0), numTrials_(DEFAULT_NUMTRIALS), cutoffSq_(0.0),
dbIndex_(REAL_VERYBIG), data_(), peakAccumulators_(), ranGen_(0)
{
  for(size_t i = 0; i < numClusters_; ++i)
    peakAccumulators_.push_back(PeakAccumulator());
  ranGen_.RandomInit((int)time(0));
}

// Interface

void PeakClustering::addFrame(System<PointMolecule> const &system)
{
  // Main loop over molecules
  for(size_t i = 0; i < system.size(); i++)
  {
    // Find neighbors
    for(size_t j = 0; j < system.size(); j++)
    {
      if(i == j) continue;
      if(cutoffSq_ &&
         norm2(system.lattice().difference(system[i].position, system[j].position)) > cutoffSq_)
        continue;
      data_.push_back(RelativeConfiguration(system[i], system[j], system.lattice()));
    }
  } // End of main loop over molecules
  knowClustering_ = false;
  ++numFrames_;
}

// Accessors

CrystalDistributionParameters const PeakClustering::parameters()
{
  if(!knowClustering_) findClusters();
  CrystalDistributionParameters parameters;
  PeakGroup meanGroup;
  for(size_t i = 0; i < peakAccumulators_.size(); i++)
    meanGroup.push_back(peakAccumulators_[i].mean());
  parameters.means.push_back(meanGroup);

  GroupParameters currentGroupParameters;
  for(size_t i = 0; i < peakAccumulators_.size(); i++)
    currentGroupParameters.push_back(peakAccumulators_[i].concentrationParameters());
  parameters.concentrationParameters.push_back(currentGroupParameters);

  currentGroupParameters.clear();
  for(size_t i = 0; i < peakAccumulators_.size(); i++)
    currentGroupParameters.push_back(peakAccumulators_[i].normalizationFactors());
  parameters.normalizationFactors.push_back(currentGroupParameters);
  return parameters;
}

std::ostream &operator<<(std::ostream &outStream, PeakClustering &peakStats)
{
  if(!peakStats.knowClustering_) peakStats.findClusters();
  outStream << "Number of points: " << peakStats.numPoints() << std::endl;
  outStream << peakStats.parameters() << std::endl;
  return outStream;
}

// Private functions

void PeakClustering::findClusters()
{
  dbIndex_ = REAL_VERYBIG;
  for(size_t i = 0; i < numClusters_; ++i)
    peakAccumulators_[i].clear();
  size_t numDataPoints = data_.size();
  if(!numFrames_ || !numDataPoints)
  {
    std::cerr << "Error in PeakClustering::findClusters() - no data to cluster!" << std::endl;
    return;
  }
  if(!numTrials_)
  {
    std::cerr << "Error in PeakClustering::findClusters() - number of trials is zero!" << std::endl;
    return;
  }
  if(verbose_)
    std::cout << "Clustering " << numDataPoints << " data points from "
              << numFrames_ << " frames." << std::endl;

  ClusterAccumulator currentClusters;
  for(size_t i = 0; i < numClusters_; ++i) currentClusters.push_back(PeakAccumulator());
  // Main loop to find best clustering
  for(size_t trial = 0; trial < numTrials_; ++trial)
  {
    if(verbose_)
      std::cout << "K-Means trial number " << trial << std::endl;
    for(size_t i = 0; i < numClusters_; ++i)
      currentClusters[i].clear();
    // K-Means++ - choose initial cluster centers
    currentClusters[0].addValue(data_[ranGen_.IRandomX(0, numDataPoints - 1)]);
    for(size_t i = 1; i < numClusters_; ++i)
    {
      // Find distance between each data point and the nearest cluster center
      std::vector<Real> distSqs;
      for (size_t j = 0; j < numDataPoints; ++j)
      {
        Real minDistSq = REAL_VERYBIG;
        for(size_t k = 0; k < i; ++k)
        {
          Real currDistSq = squareDeviation(data_[j], currentClusters[k].mean());
          if(currDistSq < minDistSq) minDistSq = currDistSq;
        }
        distSqs.push_back(minDistSq);
      }
      std::vector<Real> cumulativeDists(numDataPoints);
      cumulativeDists[0] = distSqs[0];
      for(size_t j = 1; j < numDataPoints; ++j)
        cumulativeDists[j] = cumulativeDists[j - 1] + distSqs[j];
      Real rand = cumulativeDists[numDataPoints - 1]*ranGen_.Random();
      size_t chosenOne = 0;
      if(rand > cumulativeDists[0])
        for(size_t j = 1; j < numDataPoints; ++j)
        {
          if((rand > cumulativeDists[j - 1]) &&
             (rand <= cumulativeDists[j]))
          {
            chosenOne = j;
            break;
          }
        }
      currentClusters[i].addValue(data_[chosenOne]);
    } // End of loop over clusters for K-Means++ initialization
    
    if(verbose_) 
    {
      std::cout << "Initial clusters: " << std::endl;
      for(size_t i = 0; i < numClusters_; ++i)
        std::cout << currentClusters[i].mean() << std::endl;
    }

    // K-Means loop
    size_t numIters = 0;
    Real error = REAL_VERYBIG;
    while(error > K_MEANS_TOL)
    {
      ClusterAccumulator previousClusters = currentClusters;
      for(size_t i = 0; i < numClusters_; ++i) currentClusters[i].clear();
      // Update clusters
      for(size_t i = 0; i < numDataPoints; ++i)
      {
        // Find nearest center
        Real minDistSq = REAL_VERYBIG;
        size_t nearest;
        for(size_t j = 0; j < numClusters_; ++j)
        {
          Real distSq = squareDeviation(data_[i], previousClusters[j].mean());
          if(distSq < minDistSq)
          {
            minDistSq = distSq;
            nearest = j;
          }
        }
        // Add point to nearest center
        currentClusters[nearest].addValue(data_[i]);
      } // End of cluster update loop
      // Compute the error
      error = 0.0;
      for(size_t i = 0; i < numClusters_; ++i)
        error += squareDeviation(previousClusters[i].mean(), currentClusters[i].mean());
      error = sqrt(error);
      ++numIters;
      if(verbose_)
      {
        std::cout << "New clusters: " << std::endl;
        for(size_t i = 0; i < numClusters_; ++i)
        {
          std::cout << currentClusters[i].numValues() << " : "
                    << currentClusters[i].mean() << std::endl;
        }
        std::cout << "K-Means number of iterations: " << numIters << " error = " << error << std::endl;
      }
      if(numIters > K_MEANS_MAXITER)
      {
        std::cerr << "Warning: Reached maximum number of iterations in K-Means" << std::endl;
        break;
      }
    } // End of K-Means loop

    // Calculate Davies-Bouldin index
    std::vector<Real> scatters;
    std::vector<std::vector<Real> > separations;
    for(size_t i = 0; i < numClusters_; ++i)
    {
      // Scatters - use sum of reciprocal concentration parameters
      if(currentClusters[i].numValues() < 2)
      {
        scatters.push_back(0.0);
        std::cerr << "Warning: A cluster has less than 2 members" << std::endl;
      }
      else
      {
        DistributionParameters const &concPars = currentClusters[i].concentrationParameters();
        Real scatter = 0.0;
        scatter += 1.0/concPars.distance;
        scatter += 1.0/concPars.bondOrientation;
        scatter += 1.0/concPars.relativeOrientation;
        for(size_t j = 0; j < concPars.internalDOFs.first.size(); ++j)
          scatter += 1.0/concPars.internalDOFs.first[j];
        for(size_t j = 0; j < concPars.internalDOFs.second.size(); ++j)
          scatter += 1.0/concPars.internalDOFs.second[j];
        scatters.push_back(scatter);
      }
      // Separations - use the squareDeviation norm from RelativeConfiguration
      std::vector<Real> distances;
      for(size_t j = 0; j < numClusters_; ++j) // I know, twice the work, but numClusters_ should be small
      {
        if(j == i) distances.push_back(0.0);
        else distances.push_back(squareDeviation(currentClusters[i].mean(), 
                                                 currentClusters[j].mean()));
      }
      separations.push_back(distances);
    }
    Real currentdbIndex = 0.0;
    for(size_t i = 0; i < numClusters_; ++i)
    {
      Real maxOverlap = 0.0;
      for(size_t j = 0; j < numClusters_; ++j)
      {
        if(j == i) continue;
        Real overlap = (scatters[i] + scatters[j])/separations[i][j];
        if(overlap > maxOverlap) maxOverlap = overlap;
      }
      currentdbIndex += maxOverlap;
    }
    currentdbIndex /= numClusters_;
    if(verbose_)
      std::cout << "Davies-Bouldin index of current clustering: " << currentdbIndex << std::endl;

    if(currentdbIndex < dbIndex_)
    {
      dbIndex_ = currentdbIndex;
      peakAccumulators_ = currentClusters;
    }
  } // End of main loop to find best clustering

  if(verbose_)
  {
    std::cout << "Best clustering: " << std::endl;
    for(size_t i = 0; i < numClusters_; ++i)
      std::cout << peakAccumulators_[i].mean() << std::endl;
    std::cout << "Davies-Bouldin index for best clustering: " << dbIndex_ << std::endl;
  }
  knowClustering_ = true;
}


