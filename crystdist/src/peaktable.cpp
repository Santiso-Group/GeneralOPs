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
** Peak Table class
**
** Notes:
**
** - Storing the peak table using relative configurations is convenient
**   for i/o, but inefficient - this is especially obvious in the addSystem
**   method. For a future version it would be better to store the table
**   using a std::map<std::multiset<size_t>, size_t> (i.e. as the
**   groupTable variable in findGroups and addSystem), and add a method
**   to return the set of relative configurations.
*/

#include <set>
#include <map>
#include "crystdist/include/peaktable.h"

// Constructors

PeakTable::PeakTable()
:
peakTable_(), frequencies_()
{}

PeakTable::PeakTable(System<PointMolecule> const &system, GlobalHistogram const &hist)
:
peakTable_(), frequencies_()
{ findGroups(system, hist); }

// Interface

void PeakTable::findGroups(System<PointMolecule> const &system, GlobalHistogram const &hist)
{
  if(system.size() < 2)
  {
    std::cerr << "Error in findGroups: System too small" << std::endl;
    clear();
    return;
  }
  if(hist.numEntries() == 0)
  {
    std::cerr << "Error in findGroups: Empty histogram" << std::endl;
    clear();
    return;
  }

  Lattice const &lattice = system.lattice();
  TypeData const &entryTypes = hist.entryTypes();
  std::string const &entryLabel = hist.entryLabel();

  // Check consistency
  RelativeConfiguration conf_s(system[0], system[1], lattice);
  RelativeConfiguration conf_h = hist.entry(0);

  if( (conf_s.types != entryTypes) ||
      (conf_s.name != entryLabel) ||
      (conf_s.internalDOFs.first.size() != conf_h.internalDOFs.first.size()) ||
      (conf_s.internalDOFs.second.size() != conf_h.internalDOFs.second.size()) )
  {
    std::cerr << "Error in findGroups: system and histogram are inconsistent" << std::endl;
    clear();
    return;
  }

  std::map<std::multiset<size_t>, size_t> groupTable;

  // Loop over molecules to find groups
  for(size_t i = 0; i < system.size(); ++i)
  {
    std::multiset<size_t> peakIndices;
    PointMolecule const &currentMolecule = system[i];
    Real cutoff = hist.settings().cutoff;
    Real cutoff2 = cutoff*cutoff;

    for(size_t j = 0; j < system.size(); ++j) // Loop over neighbors
    {
      if(j == i) continue;
      PointMolecule const &currentNeighbor = system[j];
      if(norm2(lattice.difference(currentMolecule.position, currentNeighbor.position)) > cutoff2) continue;
      RelativeConfiguration conf(currentMolecule, currentNeighbor, lattice);
      if(conf.distance > hist.settings().cutoff) continue;
      // Search for closest entry in the histogram
      size_t iMin = 0;  Real minNorm2 = 1E50;
      for(size_t iEntry = 0; iEntry < hist.numEntries(); ++iEntry)
      {
        Real const norm2 = squareDeviation(conf, hist.entry(iEntry));
        if(norm2 < minNorm2)
        {
          iMin = iEntry;
          minNorm2 = norm2;
        }
      }
      peakIndices.insert(iMin);
    } // End of loop over neighbors
    std::map<std::multiset<size_t>, size_t>::iterator it = groupTable.find(peakIndices);
    if(it == groupTable.end())
      groupTable.insert(std::pair<std::multiset<size_t>, size_t>(peakIndices, 1));
    else
      it->second++;
  } // End of loop over molecules

  // Store peak table and frequencies
  clear();
  for(std::map<std::multiset<size_t>, size_t>::iterator it = groupTable.begin();
      it != groupTable.end(); ++it)
  {
    PeakGroup currentGroup;
    std::multiset<size_t> peakIndices = it->first;
    for(std::multiset<size_t>::iterator jt = peakIndices.begin();
        jt != peakIndices.end(); jt++)
      currentGroup.push_back(hist.entry(*jt));
    peakTable_.push_back(currentGroup);
    frequencies_.push_back(it->second);
  }
}

void PeakTable::addSystem(System<PointMolecule> const &system, GlobalHistogram const &hist)
{
  if(system.size() < 2)
  {
    std::cerr << "Error in addSystem: System too small" << std::endl;
    return;
  }
  if(hist.numEntries() == 0)
  {
    std::cerr << "Error in addSystem: Empty histogram" << std::endl;
    return;
  }
  
  Lattice const &lattice = system.lattice();
  TypeData const &entryTypes = hist.entryTypes();
  std::string const &entryLabel = hist.entryLabel();

  // Check consistency
  RelativeConfiguration conf_s(system[0], system[1], lattice);
  RelativeConfiguration conf_h = hist.entry(0);

  if( (conf_s.types != entryTypes) ||
      (conf_s.name != entryLabel) ||
      (conf_s.internalDOFs.first.size() != conf_h.internalDOFs.first.size()) ||
      (conf_s.internalDOFs.second.size() != conf_h.internalDOFs.second.size()) )
  {
    std::cerr << "Error in addSystem: system and histogram are inconsistent" << std::endl;
    clear();
    return;
  }

  // Convert existing peak table to a map
  // Note: this is wildly inefficient, see note at the top
  // It can also crash or do funny stuff if the histogram is
  // not the same as the original one
  std::map<std::multiset<size_t>, size_t> groupTable;
  
  for(size_t i = 0; i < peakTable_.size(); ++i)
  {
    PeakGroup const &currentGroup = peakTable_[i];
    std::multiset<size_t> peakIndices;
    // Loop over group entries
    for(size_t j = 0; j < currentGroup.size(); ++j)
    {
      RelativeConfiguration const &conf = currentGroup[j];
      
      // Find closest entry in the histogram
      size_t iMin = 0; Real minNorm2 = 1E50;
      for(size_t iEntry = 0; iEntry < hist.numEntries(); ++iEntry)
      {
        Real const norm2 = squareDeviation(conf, hist.entry(iEntry));
        if(norm2 < minNorm2)
        {
          iMin = iEntry;
          minNorm2 = norm2;
        }
      }
      peakIndices.insert(iMin);
    } // End of loop over group entries
    std::map<std::multiset<size_t>, size_t>::iterator it = groupTable.find(peakIndices);
    if(it == groupTable.end())
      groupTable.insert(std::pair<std::multiset<size_t>, size_t>(peakIndices, frequencies_[i]));
    else
    {
      std::cerr << "Warning: Found a repeated group in addSystem - peak table may be corrupted" << std::endl;
      it->second += frequencies_[i];
    }
  } // End of loop over peak table entries
  
  // Loop over molecules to find groups
  for(size_t i = 0; i < system.size(); ++i)
  {
    std::multiset<size_t> peakIndices;
    PointMolecule const &currentMolecule = system[i];
    Real cutoff = hist.settings().cutoff;
    Real cutoff2 = cutoff*cutoff;

    for(size_t j = 0; j < system.size(); ++j) // Loop over neighbors
    {
      if(j == i) continue;
      PointMolecule const &currentNeighbor = system[j];
      if(norm2(lattice.difference(currentMolecule.position, currentNeighbor.position)) > cutoff2) continue;
      RelativeConfiguration conf(currentMolecule, currentNeighbor, lattice);
      if(conf.distance > hist.settings().cutoff) continue;
      // Search for closest entry in the histogram
      size_t iMin = 0;  Real minNorm2 = 1E50;
      for(size_t iEntry = 0; iEntry < hist.numEntries(); ++iEntry)
      {
        Real const norm2 = squareDeviation(conf, hist.entry(iEntry));
        if(norm2 < minNorm2)
        {
          iMin = iEntry;
          minNorm2 = norm2;
        }
      }
      peakIndices.insert(iMin);
    } // End of loop over neighbors
    std::map<std::multiset<size_t>, size_t>::iterator it = groupTable.find(peakIndices);
    if(it == groupTable.end())
      groupTable.insert(std::pair<std::multiset<size_t>, size_t>(peakIndices, 1));
    else
      it->second++;
  } // End of loop over molecules

  // Store peak table and frequencies
  clear();
  for(std::map<std::multiset<size_t>, size_t>::iterator it = groupTable.begin();
      it != groupTable.end(); ++it)
  {
    PeakGroup currentGroup;
    std::multiset<size_t> peakIndices = it->first;
    for(std::multiset<size_t>::iterator jt = peakIndices.begin();
        jt != peakIndices.end(); jt++)
      currentGroup.push_back(hist.entry(*jt));
    peakTable_.push_back(currentGroup);
    frequencies_.push_back(it->second);
  }
}

void PeakTable::reduce(size_t const freqCutoff)
{
  for(size_t i = 0; i < peakTable_.size(); ++i)
    if(frequencies_[i] < freqCutoff)
    {
      std::vector<PeakGroup>::iterator it = peakTable_.begin();
      std::vector<size_t>::iterator jt = frequencies_.begin();
      for(size_t j = 0; j < i; j++, it++, jt++) {}
      peakTable_.erase(it);
      frequencies_.erase(jt);
      i--; 
    }
}

// Print peak table (mostly for debugging)

std::ostream &operator<<(std::ostream &outStream, PeakTable const &peakTable)
{
  outStream << "Number of groups: " << peakTable.numGroups() << std::endl;
  for(size_t i = 0; i < peakTable.numGroups(); ++i)
  {
    outStream << "Frequency: " << peakTable.frequency(i) << std::endl;
    outStream << "Peaks: " << std::endl;
    for(size_t j = 0; j < peakTable.peaks(i).size(); ++j)
      outStream << peakTable.peaks(i)[j] << std::endl;
  }
  return outStream;
}
