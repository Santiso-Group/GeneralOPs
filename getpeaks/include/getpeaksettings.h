/*
** Copyright 2008-2011 Erik Santiso.
** This file is part of getpeaks.
** getpeaks is free software: you can redistribute it and/or modify it
** under the terms of the GNU Lesser General Public License as published
** version 2.1 as published by the Free Software Foundation.
** 
**
** getpeaks is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with getpeaks. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** Getpeak Settings class
**
** This implements a container for settings used to construct the
** global histogram of distances, bond orientations, relative
** orientations and internal degrees of freedom in an ordered system.
** Most of the quantities in this class are defined in class
** HistogramSettings in the crystdist library, with the exception of:
**
** - The norm cutoff for merging peaks in the global histogram, and
**   the frequency cutoff for trimming the global histogram. These
**   are defined in the reduce() method of class GlobalHistogram.
** - The frequncy cutoff for trimming the final group table. This
**   is defined in the reduce() method of class PeakTable.
**
** See the files globalhistogram.h and peaktable.h in the include/
** directory of class crystdist for more details. 
*/

#ifndef H_GET_PEAK_SETTINGS
#define H_GET_PEAK_SETTINGS

#include "common/include/types.h"
#include "common/include/iofile.h"
#include "crystdist/include/globalhistogram.h"

Real const DEFAULT_NORM_CUTOFF = 0.1;
size_t const DEFAULT_FREQ_CUTOFF_GLOBAL_HISTOGRAM = 2;
size_t const DEFAULT_FREQ_CUTOFF_PEAK_TABLE = 2;

struct GetpeakSettings
{
  HistogramSettings histogramSettings;
  Real normCutoff;
  size_t histFreqCutoff;
  size_t groupFreqCutoff;

  GetpeakSettings() // Assign default values
  :
  histogramSettings(), normCutoff(DEFAULT_NORM_CUTOFF),
  histFreqCutoff(DEFAULT_FREQ_CUTOFF_GLOBAL_HISTOGRAM),
  groupFreqCutoff(DEFAULT_FREQ_CUTOFF_PEAK_TABLE)
  {}
 
  void reset() // Assigns default values
  {
    histogramSettings.reset();
    normCutoff = DEFAULT_NORM_CUTOFF;
    histFreqCutoff = DEFAULT_FREQ_CUTOFF_GLOBAL_HISTOGRAM;
    groupFreqCutoff = DEFAULT_FREQ_CUTOFF_PEAK_TABLE;
  }

  // Input
  friend std::istream &operator>>(std::istream &inStream, GetpeakSettings &settings)
  {
    // Read histogram settings
    inStream >> settings.histogramSettings;
    // Read the cutoffs
    inStream >> settings.normCutoff; inStream.ignore(LINE_LENGTH, '\n');
    inStream >> settings.histFreqCutoff; inStream.ignore(LINE_LENGTH, '\n');
    inStream >> settings.groupFreqCutoff; inStream.ignore(LINE_LENGTH, '\n');
    if( (!inStream && !inStream.eof()) ||
        (settings.normCutoff < 0.0) )
    {
      std::cerr << "Error reading settings - reverting to default values" << std::endl;
      settings.reset();
    }
    return inStream;
  }
  
  // Output
  friend std::ostream &operator<<(std::ostream &outStream, GetpeakSettings const &settings)
  {
    outStream << settings.histogramSettings;
    outStream << settings.normCutoff << " ! Norm cutoff for global histogram" << std::endl;
    outStream << settings.histFreqCutoff << " ! Frequency cutoff for global histogram" << std::endl;
    outStream << settings.groupFreqCutoff << " ! Frequency cutoff for group table" << std::endl;
    return outStream;
  }
};

/*
** End of class GetpeakSettings
*/

#endif

