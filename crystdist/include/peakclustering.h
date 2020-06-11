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
**
** This class estimates parameters in an approximate distribution
** function of relative configurations given the results of a
** finite-temperature simulation. Unlike PeakStatistics, this
** class uses a clustering algorithm (K-Means++) to build the
** peak table, and does not separate the tables into groups.
**
** To use, keep adding frames from the finite-temperature simulation
** using the addPoint() method, and then use the accessors to
** obtain averages, concentration parameters, and normalization 
** factors for each peak.
**
** In general, it is a good idea to use a cutoff to keep the amount
** of data points manageable. If it is not set (or set to 0), all
** molecules are included in the calculation.
*/

/*
** The same notes in class PeakClustering apply here.
** Also, for a future version it may be a good idea to implement
** K-Means++ as a template in mymath and modify this to call it.
**
** This uses the randomc library by Agner Fog
*/

#ifndef H_PEAK_CLUSTERING
#define H_PEAK_CLUSTERING
#define DEFAULT_NUMTRIALS 50
#define DEFAULT_NUMCLUSTERS 3
#define K_MEANS_TOL 1.E-3
#define K_MEANS_MAXITER 1000

#include "common/include/assert.h"
#include "common/include/types.h"
#include "mymol/include/system.h"
#include "mymath/include/hungarian.h"
#include "crystdist/include/statparameters.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/peakaccumulator.h"
#include "randomc/include/randomc.h"

// The notion of "Cluster" in this class is the same as "Group" in PeakStatistics

typedef std::vector<PeakAccumulator> ClusterAccumulator;
typedef std::vector<RelativeConfiguration> RelativeConfigurationData;

class PeakClustering
{
public:

// Constructors

  PeakClustering(); // Defines an empty peak clustering object

// Interface

  void addFrame(System<PointMolecule> const &system); // Adds a frame to the clustering
  void setNumClusters(size_t const numClusters);      // Set the number of clusters
  void setNumTrials(size_t const numTrials);          // Sets the number of K-Means trials
  void setCutoff(Real const &cutoff);                 // Set the cutoff
  void clear();                                       // Clears the clustering (all except the
                                                      // cutoff and numbers of clusters and trials)
  void verboseOutput(bool const shouldITalk);         // Turns output to std::cout while clustering
                                                      // on or off.

// Accessors

  size_t const numClusters() const;                 // Returns the number of clusters
  size_t const numPoints() const;                   // Returns the total number of data points
  size_t const numFrames() const;                   // Returns the number of frames added so far
  size_t const numTrials() const;                   // Returns the number of K-Means trials
  Real const cutoff() const;                        // Returns the cutoff used to calculate data
  Real const &dbIndex();                            // Returns the Davies-Bouldin index
  CrystalDistributionParameters const parameters(); // Returns the statistical parameters from the
                                                    // data accumulated so far

  // Need accessors for other properties, e.g. DB index

// Print peak statistics (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, PeakClustering &peakClustering);

private:

  bool knowClustering_;                 // Whether the clustering is known
  bool verbose_;                        // Produce verbose output to std::cout
  size_t numClusters_;                  // Number of clusters for K-Means
  size_t numFrames_;                    // Number of frames added so far
  size_t numTrials_;                    // Number of K-Means trials to use
  Real cutoffSq_;                       // Center-of-mass distance cutoff, squared
  Real dbIndex_;                        // Davies-Bouldin index for the clustering
  RelativeConfigurationData data_;      // Relative configurations added so far 
  ClusterAccumulator peakAccumulators_; // Accumulators for all the peaks
  CRandomMersenne ranGen_;              // Random number generator

// Private functions
  void findClusters(); // Do the clustering
};

/*
** End of class PeakClustering
*/

// Inlines

inline void PeakClustering::setNumClusters(size_t const numClusters)
{
  if(!numClusters)
  {
    std::cerr << "Error in PeakClustering::setNumClusters - number of clusters must be positive"
              << std::endl;
    return;
  }
  if(numClusters_ != numClusters) 
  {
    knowClustering_ = false;
    peakAccumulators_.clear();
    for(size_t i = 0; i < numClusters; ++i)
      peakAccumulators_.push_back(PeakAccumulator());
  }
  numClusters_ = numClusters;
}

inline void PeakClustering::setNumTrials(size_t const numTrials)
{
  if(numTrials_ != numTrials) knowClustering_ = false;
  numTrials_ = numTrials;
}

inline void PeakClustering::setCutoff(Real const &cutoff)
{
  if(cutoff*cutoff != cutoffSq_) knowClustering_ = false;
  cutoffSq_ = cutoff*cutoff;
}

inline void PeakClustering::clear()
{
  knowClustering_ = false;
  numFrames_ = 0;
  data_.clear();
  dbIndex_ = REAL_VERYBIG;
  for(size_t i = 0; i < numClusters_; ++i)
    peakAccumulators_[i].clear();
}

inline void PeakClustering::verboseOutput(bool const shouldITalk)
{
  verbose_ = shouldITalk;
}

inline size_t const PeakClustering::numClusters() const
{
  return numClusters_;
}

inline size_t const PeakClustering::numPoints() const
{
  return data_.size();
}

inline size_t const PeakClustering::numFrames() const
{
  return numFrames_;
}

inline size_t const PeakClustering::numTrials() const
{
  return numTrials_;
}

inline Real const PeakClustering::cutoff() const
{
  return sqrt(cutoffSq_);
}

inline Real const &PeakClustering::dbIndex()
{
  if(!knowClustering_) findClusters();
  return dbIndex_;
}

#endif

