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
** Peak accumulator class
**
** This class estimates parameters in an approximate distribution
** function of relative configurations given a peak table and the
** results of a finite-temperature simulation.
*/

/*
** Some notes:
** 
** - Technically, it would make more sense to implement the ANGLE
**   type for DistributionVariable in mymath and create a class to
**   keep track of accumulators for different DistributionVariables.
**   This would take quite a bit of re-coding, however, so it'd be
**   better to leave this for a future version (the current 
**   implementation of DistributionVariable, FrequencyTable, and
**   GlobalHistogram is quite cumbersome anyway).
** - Note that addValue() only checks that the name of the relative
**   configuration being added is consistent with the data already
**   in the accumulator. If you try to add a different type of data
**   with the same name string, results are unpredictable
*/

#ifndef H_PEAK_ACCUMULATOR
#define H_PEAK_ACCUMULATOR

#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include "common/include/assert.h"
#include "mymol/include/system.h"
#include "mymath/include/matrix4D.h"
#include "crystdist/include/relativeconfiguration.h"
#include "crystdist/include/statparameters.h"

struct InternalDOFAccumulator // Accumulator for internal DOFs
{
  InternalDOFType type;   // Type of internal DOF (see moleculemap.h)
  size_t symmetryNumber;  // Symmetry number (only meaningful for dihedrals)
  Real sumC, sumS;        // Accumulators

  // The accumulators can have two interpretations:
  // - If type = DISTANCE, sumC is the sum of the distances and sumS is
  //   the sum of the distances squared. The underlying distribution is Gaussian.
  // - If type = ANGLE or DIHEDRAL, sumC is the sum of the cosines and sumS
  //   the sum of the sines (with the angle multiplied by the appropriate
  //   factor). The underlying distribution is von Mises.

  InternalDOFAccumulator()
  :
  type(DISTANCE), symmetryNumber(1), sumC(0.0), sumS(0.0)
  {}

  InternalDOFAccumulator(InternalDOFType const &type, size_t const symmetryNumber)
  :
  type(type), symmetryNumber(symmetryNumber), sumC(0.0), sumS(0.0)
  {}

  void clear()
  {
    type = DISTANCE;
    symmetryNumber = 1;
    sumC = sumS = 0.0;
  }
};

typedef std::pair<std::vector<InternalDOFAccumulator>, 
                  std::vector<InternalDOFAccumulator> > 
                  InternalDOFDataAccumulator;

class PeakAccumulator  // Accumulates relative configuration values
{
public:

// Constructors

  PeakAccumulator();  // Defines an empty peak accumulator

// Interface

  void addValue(RelativeConfiguration const &conf); // Adds a relative configuration
  void clear();                                     // Clears the accumulator

// Accessors

  size_t const numValues() const;                           // Returns the number of values accumulated
  RelativeConfiguration const &mean();                      // Returns the mean relative configuration
  DistributionParameters const &concentrationParameters();  // Returns the concentration parameters
  DistributionParameters const &normalizationFactors();     // Returns the normalization factors

// Target functions for maximum likelihood estimation

  Real const vonMisesTarget(Real const &concentrationParameter);
  Real const vonMisesFisherTarget(Real const &concentrationParameter);
  Real const bipolarWatsonTarget(Real const &concentrationParameter);

// Print current statistics (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, PeakAccumulator &accumulator);

private:

  bool knowParameters_;                             // Whether the distribution parameters are known
  size_t numValues_;                                // Number of values added to the accumulator
  std::string name_;                                // Name of the data in the accumulator
  TypeData types_;                                  // Type data in the accumulator (see relativeconfiguration.h)
  Real sumDistances_;                               // Distance accumulator
  Real sumDistancesSquared_;                        // Distance squared accumulator
  Quaternion sumBondOrientations_;                  // Bond orientation accumulator
  Matrix4D sumRelativeOrientations_;                // Relative orientation accumulator
  InternalDOFDataAccumulator sumInternalDOFData_;   // Accumulator for Internal DOFs
  RelativeConfiguration mean_;                      // Mean relative configuration
  DistributionParameters concentrationParameters_;  // Concentration parameters
  DistributionParameters normalizationFactors_;     // Normalization factors
  Real lengthOfMeanVector_;                         // Length of mean vector (used for likelihood maximization)

// Private functions

  void calculateParameters(); // Calculates the distribution parameters
};

/*
** End of class PeakAccumulator
*/

// Inlines

inline size_t const PeakAccumulator::numValues() const
{
  return numValues_;
}

inline RelativeConfiguration const &PeakAccumulator::mean()
{
  if(!knowParameters_) calculateParameters();
  return mean_;
}

inline DistributionParameters const &PeakAccumulator::concentrationParameters()
{
  if(!knowParameters_) calculateParameters();
  return concentrationParameters_;
}

inline DistributionParameters const &PeakAccumulator::normalizationFactors()
{
  if(!knowParameters_) calculateParameters();
  return normalizationFactors_;
}

#endif
