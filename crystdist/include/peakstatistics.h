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
**
** This class estimates parameters in an approximate distribution
** function of relative configurations given a peak table and the
** results of a finite-temperature simulation.
**
** To use, keep adding frames from the finite-temperature simulation
** using the addPoint() method, and then use the accessors to
** obtain averages, concentration parameters, and normalization 
** factors for each peak.
*/

/*
** Some notes:
** 
** - This class, like CrystalOrderParameters, was written for
**   the case of a homogeneous system. For hybrid relative
**   configurations this should be changed (see note in 
**   crystdist.cpp).
** - This class assigns peaks in the finite-temperature data
**   to the peaks in the zero-temperature peak table by
**   using Hungarian matching on the N nearest neighbors of
**   each molecule, where N is the number of peaks per group.
**   Note that this can leave out one or more of the right
**   neighbors depending on the choice of cutoff used to
**   generate the group table. Ideally, the cutoff should be
**   chosen so that there is a gap as wide as possible between
**   the last peak counted and the next one. This problem
**   could be avoided by implementing a rectangular-matrix
**   version of the Hungarian algorithm in mymath.
** - Technically, it would make more sense to implement the ANGLE
**   type for DistributionVariable in mymath and create a class to
**   keep track of accumulators for different DistributionVariables.
**   This would take quite a bit of re-coding, however, so it'd be
**   better to leave this for a future version (the current 
**   implementation of DistributionVariable, FrequencyTable, and
**   GlobalHistogram is quite cumbersome anyway).
*/

#ifndef H_PEAK_STATISTICS
#define H_PEAK_STATISTICS

#include "common/include/assert.h"
#include "mymol/include/system.h"
#include "mymath/include/hungarian.h"
#include "crystdist/include/statparameters.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/peakaccumulator.h"
#include "crystdist/include/peaktable.h"

typedef std::vector<PeakAccumulator> GroupAccumulator;
typedef std::vector<GroupAccumulator> TableAccumulator;

class PeakStatistics
{
public:

// Constructors

  PeakStatistics(); // Defines an empty peak statistics object

// Interface

  void addPoint(System<PointMolecule> const &system,  // Adds a point to the peak statistics
                PeakTable const &table);
  void clear();                                       // Clears the statistics

// Accessors

  size_t const numPoints() const;                   // Returns the number of points added so far
  CrystalDistributionParameters const parameters(); // Returns the statistical parameters from the
                                                    // data accumulated so far

// Print peak statistics (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, PeakStatistics &peakStats);

private:

  size_t numPoints_;                  // Number of values added so far
  TableAccumulator peakAccumulators_; // Accumulators for all the peaks
};

/*
** End of class PeakStatistics
*/

// Inlines

inline void PeakStatistics::clear()
{
  numPoints_ = 0;
  peakAccumulators_.clear();
}

inline size_t const PeakStatistics::numPoints() const
{
  return numPoints_;
}

#endif
