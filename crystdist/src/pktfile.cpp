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
** .pkt file format
*/

#include "common/include/assert.h"
#include "crystdist/include/pktfile.h"

// Constructors

PKTFile::PKTFile()
:
IOFile()
{}

PKTFile::PKTFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

PKTFile &operator<<(PKTFile &pktFile, PeakTable const &table)
{
  // Write table to file
  assert(pktFile.is_open() && pktFile.mode() == OUT);

  pktFile.setf(std::ios::fixed, std::ios::floatfield);
  pktFile.setf(std::ios::right, std::ios::adjustfield);
  std::streamsize precision = pktFile.precision();
  pktFile.precision(WRITE_PRECISION);

  // Number of groups
  pktFile << "!NGROUPS: " << std::endl;
  pktFile << table.numGroups() << std::endl;
  
  // Group information
  for(size_t i = 0; i < table.numGroups(); i++)
  {
    pktFile << "!FREQUENCY: " << std::endl;
    pktFile << table.frequency(i) << std::endl;
    pktFile << "!PEAKS: " << std::endl;
    pktFile << table.peaks(i).size() << std::endl;
    for(size_t j = 0; j < table.peaks(i).size(); j++)
      pktFile << table.peaks(i)[j] << std::endl;
  }

  pktFile.precision(precision); // Return to original value
  return pktFile;
}

PKTFile &operator>>(PKTFile &pktFile, PeakTable &table)
{
  // Read table from file
  assert(pktFile.is_open() && pktFile.mode() == IN);
  table.clear();

  std::string r_label;  // Label read from file

  // Search for number of groups label
  do
  {
    std::getline(pktFile, r_label);
    if(r_label.find("!NGROUPS:") != r_label.npos) break;
  }
  while(!pktFile.eof());
  if(pktFile.eof())
  {
    std::cerr << "Error reading file " << pktFile.fileName() << ": NGROUPS record not found." << std::endl;
    return pktFile;
  }

  // Read number of groups
  size_t r_nGroups;
  pktFile >> r_nGroups;
  if(!pktFile)
  {
    std::cerr << "Error reading number of groups from file " << pktFile.fileName() << std::endl;
    return pktFile;
  }

  // Read group information
  for(size_t i = 0; i < r_nGroups; i++)
  {
    size_t r_freq;                // Frequency read from file
    size_t r_nPeaks;              // Number of peaks read from file
    RelativeConfiguration r_conf; // Relative configuration read from file
    PeakGroup r_group;            // Peak group as read from file
    pktFile >> r_label;
    pktFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    pktFile >> r_freq;
    pktFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    if((r_label.find("!FREQUENCY:") == r_label.npos) || !pktFile)
    {
      std::cerr << "Error reading frequency record from file " << pktFile.fileName() << std::endl;
      table.clear();
      return pktFile;
    }
    pktFile >> r_label;
    pktFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    pktFile >> r_nPeaks;
    pktFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    if((r_label.find("!PEAKS:") == r_label.npos) || !pktFile)
    {
      std::cerr << "Error reading peak record from file " << pktFile.fileName() << std::endl;
      table.clear();
      return pktFile;
    }
    r_group.clear();
    for(size_t j = 0; j < r_nPeaks; j++)
    {
      pktFile >> r_conf;
      pktFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
      if(!pktFile)
      {
        std::cerr << "Error reading peak record from file " << pktFile.fileName() << std::endl;
        table.clear();
        return pktFile;
      }
      r_group.push_back(r_conf);
    }
    table.addGroup(r_group, r_freq);
  }
  return pktFile;
}
