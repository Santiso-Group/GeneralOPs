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
** Peak table class
**
** A peak table contains the locations of all peaks in the distribution
** function of positions, bond orientations, and relative orientations
** for a crystal or liquid crystal (usually at zero temperature). The
** peaks are separated in groups, which usually correspond to the
** various symmetry-unique molecules in the primitive cell.
**
** This class implements methods to build a peak table given a global
** histogram and a system of point molecules corresponding to that
** histogram.
*/

/*
** Some notes:
**
** - Like PeakStatistics and CrystalOrderParameters, this class was 
**   written for the case of a homogeneous system. The case of
**   hybrid relative configurations is not implemented yet.
*/

#ifndef H_PEAK_TABLE
#define H_PEAK_TABLE

#include <algorithm>
#include "common/include/assert.h"
#include "mymol/include/system.h"
#include "crystdist/include/statparameters.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/relativeconfiguration.h"
#include "crystdist/include/globalhistogram.h"

class PeakTable
{
public:

// Constructors

  PeakTable();                                    // Defines an empty group table
  PeakTable(System<PointMolecule> const &system,  // Defines a group table for the given system
            GlobalHistogram const &hist);         // and the given global histogram

// Interface

  void addGroup(PeakGroup const &group, size_t const &freq);  // Adds the given group with the given
                                                              // frequency (useful for file input)
  void findGroups(System<PointMolecule> const &system,        // Builds the group table for the given system
                  GlobalHistogram const &hist);               // and the given global histogram
  void addSystem(System<PointMolecule> const &system,         // Adds a system snapshot to an existing group
                 GlobalHistogram const &hist);                // table - not that this does *not* check that
                                                              // the histogram is consistent with the table.
  void reduce(size_t const freqCutoff);                       // Removes all entries from the group table
                                                              // with a frequency smaller than the given cutoff
  void sort();                                                // Sorts the peaks in the group table
  void clear();                                               // Clears the group table

// Accessors

  PeakGroup const &peaks(size_t const index) const; // Returns the peaks in the index-th group
  size_t const frequency(size_t const index) const; // Returns the frequency of the index-th group
  size_t const numGroups() const;                   // Returns the number of peak groups

// Print peak table (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, PeakTable const &peakTable); 

private:

  std::vector<PeakGroup> peakTable_;  // Peak groups in the system
  std::vector<size_t> frequencies_;   // Frequency of each group
};

/*
** End of class PeakTable
*/

// Inlines

inline void PeakTable::addGroup(PeakGroup const &group, size_t const &freq)
{
  peakTable_.push_back(group);
  frequencies_.push_back(freq);
}

inline void PeakTable::sort()
{
  for(size_t i = 0; i < peakTable_.size(); i++)
    std::stable_sort(peakTable_[i].begin(), peakTable_[i].end());
}

inline void PeakTable::clear()
{
  peakTable_.clear();
  frequencies_.clear();
}

inline PeakGroup const &PeakTable::peaks(size_t const index) const
{
  assert(index < peakTable_.size());
  return peakTable_[index];
}

inline size_t const PeakTable::frequency(size_t const index) const
{
  assert(index < frequencies_.size());
  return frequencies_[index];
}

inline size_t const PeakTable::numGroups() const
{
  return peakTable_.size();
}

#endif
