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
** Global histogram class
**
** This implements a frequency table for the distribution of relative
** configurations in a system of point molecules.
*/

/*
** Some notes:
**
** - Currently it is only possible to define one cutoff for distance
**   internal DOFs, and the class uses the same number of bins for 
**   all internal DOFs of the same kind. If needed, it would be better
**   to add the possibility to change these values independently for
**   each internal DOF.
**
** - For a future version, it might be better to define ANGLE and AXIS
**   distribution variable types (see mymath/distributionvariable.h),
**   and collect statistics as 2D vectors/axes.
*/

#ifndef H_GLOBAL_HISTOGRAM
#define H_GLOBAL_HISTOGRAM

#include <utility>
#include "common/include/iofile.h"
#include "mymath/include/distributionvariable.h"
#include "mymath/include/frequencytable.h"
#include "mymol/include/system.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/relativeconfiguration.h"

size_t const DEFAULT_NUM_DISTANCE_BINS = 100;
size_t const DEFAULT_NUM_VECTOR_BINS = 21;
size_t const DEFAULT_NUM_QUATERNION_BINS = 11;
size_t const DEFAULT_NUM_ANGLE_BINS = 361;
size_t const DEFAULT_NUM_AXIS_BINS = 181;
Real const DEFAULT_DISTANCE_CUTOFF = 10.0;

struct InternalDOFDefinition  // Definition of individual internal DOFs
{
  InternalDOFType type;       // Type of internal DOF
  size_t symmetryNumber;      // Symmetry number (only meaningful for dihedrals)

  InternalDOFDefinition()
  :
  type(DISTANCE), symmetryNumber(1)
  {}
};

typedef std::vector<InternalDOFDefinition> InternalDOFInfo; // Internal DOF info for a single molecule

struct HistogramSettings
{
  size_t numDistanceBins;         // Number of bins for intermolecular distances
  size_t numAngleBins;            // Number of bins for angles (used for bond and relative 
                                  // orientations in the planar/linear cases)
  size_t numVectorBins;           // Number of bins in each direction for vectors
                                  // (used for bond orientations in the general case)
  size_t numQuaternionBins;       // Number of bins in each direction for quaternions
                                  // (used for relative orientations in the general case)
  size_t numInternalDistanceBins; // Number of bins for distance internal DOFs
  size_t numInternalAngleBins;    // Number of bins for angle internal DOFs
  size_t numInternalDihedralBins; // Number of bins for dihedral internal DOFs
  Real cutoff;                    // Cutoff for intermolecular distances
  Real internalCutoff;            // Cutoff for internal distance DOFs

  HistogramSettings() // Assign default values
  :
  numDistanceBins(DEFAULT_NUM_DISTANCE_BINS),
  numAngleBins(DEFAULT_NUM_AXIS_BINS),
  numVectorBins(DEFAULT_NUM_VECTOR_BINS),
  numQuaternionBins(DEFAULT_NUM_QUATERNION_BINS),
  numInternalDistanceBins(DEFAULT_NUM_DISTANCE_BINS),
  numInternalAngleBins(DEFAULT_NUM_AXIS_BINS),
  numInternalDihedralBins(DEFAULT_NUM_ANGLE_BINS),
  cutoff(DEFAULT_DISTANCE_CUTOFF),
  internalCutoff(DEFAULT_DISTANCE_CUTOFF)
  {}

  void reset()  // Assigns default values
  {
    numDistanceBins = DEFAULT_NUM_DISTANCE_BINS;
    numAngleBins = DEFAULT_NUM_AXIS_BINS;
    numVectorBins = DEFAULT_NUM_VECTOR_BINS;
    numQuaternionBins = DEFAULT_NUM_QUATERNION_BINS;
    numInternalDistanceBins = DEFAULT_NUM_DISTANCE_BINS;
    numInternalAngleBins = DEFAULT_NUM_AXIS_BINS;
    numInternalDihedralBins = DEFAULT_NUM_ANGLE_BINS;
    cutoff = DEFAULT_DISTANCE_CUTOFF;
    internalCutoff = DEFAULT_DISTANCE_CUTOFF;
  }

  // Input
  friend std::istream &operator>>(std::istream &inStream, HistogramSettings &settings)
  {
    inStream >> settings.numDistanceBins; inStream.ignore(LINE_LENGTH, '\n');
    inStream >> settings.numAngleBins; inStream.ignore(LINE_LENGTH, '\n');
    inStream >> settings.numVectorBins; inStream.ignore(LINE_LENGTH, '\n'); 
    inStream >> settings.numQuaternionBins; inStream.ignore(LINE_LENGTH, '\n'); 
    inStream >> settings.numInternalDistanceBins; inStream.ignore(LINE_LENGTH, '\n'); 
    inStream >> settings.numInternalAngleBins; inStream.ignore(LINE_LENGTH, '\n'); 
    inStream >> settings.numInternalDihedralBins; inStream.ignore(LINE_LENGTH, '\n'); 
    inStream >> settings.cutoff; inStream.ignore(LINE_LENGTH, '\n'); 
    inStream >> settings.internalCutoff; inStream.ignore(LINE_LENGTH, '\n'); 
    if( (!inStream && !inStream.eof()) ||
        (settings.cutoff < 0.0) ||
        (settings.internalCutoff < 0.0) )
    {
      std::cerr << "Error reading histogram settings - reverting to default values" << std::endl;
      settings.reset();
    }
    return inStream;
  }

  // Output
  friend std::ostream &operator<<(std::ostream &outStream, HistogramSettings const &settings)
  {
    outStream << settings.numDistanceBins << " ! Number of distance bins" << std::endl;
    outStream << settings.numAngleBins << " ! Number of angle bins" << std::endl;
    outStream << settings.numVectorBins << " ! Number of vector bins" << std::endl;;
    outStream << settings.numQuaternionBins << " ! Number of quaternion bins" << std::endl;
    outStream << settings.numInternalDistanceBins << " ! Number of bins for distance DOFs" << std::endl;
    outStream << settings.numInternalAngleBins << " ! Number of bins for angle DOFs" << std::endl;
    outStream << settings.numInternalDihedralBins << " ! Number of bins for dihedral DOFs" << std::endl;
    outStream << settings.cutoff << " ! Distance cutoff" << std::endl;
    outStream << settings.internalCutoff << " ! Internal DOF cutoff" << std::endl;
    return outStream;
  }
};


class GlobalHistogram
{
public:

// Constructors

  GlobalHistogram();  // Define an empty global histogram

// Interface

  // Note that changing histogram settings clears all data from the histogram - set
  // them before adding data points.
  void changeSettings(HistogramSettings const &settings); // Changes the histogram settings
  void addPoint(RelativeConfiguration const &conf);       // Adds data from a single pair
  void addSystem(System<PointMolecule> const &system,     // Adds data for all the molecules in the given
                 std::string const &name);                // system that have the given name
  void clearData();                                       // Clears all data and variable types
  void clear();                                           // Clears all data and variable types and resets 
                                                          // cutoff and numbers of bins to default values

// Utilities

  void reduce(Real const &normCutoff,       // Merges all peaks for which the total difference
              size_t const freqCutoff = 0); // between entries is less than the given cutoff and,
                                            // optionally, removes entries with frequency less than
                                            // the given frequency cutoff afterwards.

// Accessors

  TypeData const &entryTypes() const;                       // Returns the molecule types stored in the histogram
  std::string const &entryLabel() const;                    // Returns the entry label
  HistogramSettings const &settings() const;                // Returns the histogram settings
  size_t const numEntries() const;                          // Returns the number of entries in the histogram
  RelativeConfiguration const entry(size_t const i) const;  // Returns a relative configuration corresponding to
                                                            // the i-th entry in the histogram
  std::pair<size_t, size_t> const numIntDOFs() const;       // Returns the numbers of internal degrees of freedom

// Print histogram data

  friend std::ostream &operator<<(std::ostream &outStream, GlobalHistogram const &hist); 

private:

  TypeData entryTypes_;                                         // Molecule types stored in the histogram
  std::string entryLabel_;                                      // Histogram entry label (residue IDs)
  HistogramSettings settings_;                                  // Histogram settings (see struct above)
  FrequencyTable frequencyTable_;                               // The histogram
  std::pair<InternalDOFInfo, InternalDOFInfo> internalDOFInfo_; // Definitions of the internal DOFs
};

/*
** End of class GlobalHistogram
*/

// Inlines

inline void GlobalHistogram::changeSettings(HistogramSettings const &settings)
{
  settings_ = settings;
  clearData();
}

inline void GlobalHistogram::clearData()
{
  entryTypes_ = TypeData(GENERAL, GENERAL);
  entryLabel_.clear();
  frequencyTable_.clear();
  internalDOFInfo_ = std::pair<InternalDOFInfo, InternalDOFInfo>();
}

inline void GlobalHistogram::clear()
{
  entryTypes_ = TypeData(GENERAL, GENERAL);
  entryLabel_.clear();
  settings_.reset();
  frequencyTable_.clear();
  internalDOFInfo_ = std::pair<InternalDOFInfo, InternalDOFInfo>();
}

inline TypeData const &GlobalHistogram::entryTypes() const
{
  return entryTypes_;
}

inline std::string const &GlobalHistogram::entryLabel() const
{
  return entryLabel_;
}

inline HistogramSettings const &GlobalHistogram::settings() const
{
  return settings_;
}

inline size_t const GlobalHistogram::numEntries() const
{
  return frequencyTable_.numEntries();
}

inline std::pair<size_t, size_t> const GlobalHistogram::numIntDOFs() const
{
  return std::pair<size_t, size_t>(internalDOFInfo_.first.size(), internalDOFInfo_.second.size());
}

#endif

