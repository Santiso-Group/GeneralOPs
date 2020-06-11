/*
** Copyright 2012 Erik Santiso.
** This file is part of mymath.
** mymath is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** mymath is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with mymath. If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * bootstrap - a command-line utility to generate bootstrap distributions
 * 
 */

// At this stage this is not production-level code, just for testing purposes

#include <vector>
#include <string>
#include <sstream>
#include "common/include/types.h"
#include "common/include/iofile.h"
#include "mymath/include/statistics.h"

#define MAX_NUMSAMPLES 10000

int main(int argc, char* argv[])
{
  if(argc < 3)
  {
    std::cerr << "Usage: bootstrap <file> <num_samples>" << std::endl;
    return -1;
  }

  size_t numSamples;
  std::string const numString(argv[2]);
  std::stringstream numSamplesStream(numString);
  numSamplesStream >> numSamples;

  if(numSamplesStream.fail())
  {
    std::cerr << "Usage: bootstrap <file> <num_samples>" << std::endl;
    return -1;
  }

  if(numSamples > MAX_NUMSAMPLES)
  {
    std::cerr << "Number of samples must be > 0 and <= "
              << MAX_NUMSAMPLES << std::endl;
    return -1;
  }

  std::string const fileName(argv[1]);
  IOFile dataFile(fileName, IN);
  if(dataFile.fail())
  {
    std::cerr << "File " << fileName << " not found!" << std::endl;
    return -1;
  }
  
  std::vector<Real> data;

  while(!dataFile.eof())
  {
    Real point;
    dataFile >> point;
    if(dataFile.eof()) break;
    if(dataFile.fail()) 
    {
      std::cerr << "Error reading file " << fileName << std::endl;
      return -1;
    }
    data.push_back(point);
  }

  std::vector<Real> meanSample = 
    bootstrap<std::vector<Real>, Real>(data, numSamples);

  for(size_t i = 0; i < numSamples; ++i)
    std::cout << meanSample[i] << std::endl;
}

