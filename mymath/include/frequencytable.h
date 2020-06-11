/*
** Copyright 2007-2011 Erik Santiso.
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

/*
** Frequency table class
**
** This implements a frequency table for an arbitrary number of 
** variables. When the number of variables is large, the available
** memory limits this to sparse distributions.
**
** Note that the addValue method without bin indices is variadic -
** The first argument is the number of variables (which is redundant,
** since it should be variables_.size(), but <cstdarg> requires a
** first argument). The remaining arguments are the values to
** add, whose types should match those defined in variables_.
** Unfortunately, this is not very safe - if you do not supply the
** right number of arguments, or they are not the right type,
** the result is undefined.
**
** In order to recover the values of the variables that correspond
** to a set of bin indices, you need to invoke the appropriate 
** method from DistributionVariable. For example, if the second
** variable is a unit 3D vector, the value corresponding to the
** i-th frequency table entry for this variable is:
** variable(2).unitVector3DValue(entry(i).indices)
*/

/*
** Note: It would be good to implement this in a better way, e.g.
** using a variadic template. I'm leaving this for a future version.
*/

#ifndef H_FREQUENCY_TABLE
#define H_FREQUENCY_TABLE

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <cstdarg>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/distributionvariable.h"

struct FrequencyEntry // Entry in the frequency table
{
  std::vector<BinIndices> indices;
  size_t frequency;

  FrequencyEntry()
  :
  indices(), frequency(0)
  {}

  FrequencyEntry(std::vector<BinIndices> const &indices, size_t const frequency)
  :
  indices(indices), frequency(frequency)
  {}
};

class FrequencyTable
{
public:

// Constructors

  FrequencyTable(); // Define an empty frequency table

// Interface

  void addVariable(DistributionVariable const &var);      // Add an independent variable - note that
                                                          // this clears all the frequency entries
  void deleteVariable(size_t const index);                // Deletes the index-th independent variable
                                                          // this clears all the frequency entries
  void addValue(std::vector<BinIndices> const &indices);  // Add a value as a set of bin indices
/*
  // The variadic method below makes g++ puke - can be used with Microsoft's C++ compiler though
  void addValue(size_t const num, ...);                   // Add a value as a set of values of the
                                                          // independent variables
*/
  void setVariable(size_t const index,                    // Sets the index-th variable to the given
                   DistributionVariable const &var);      // one - this clears the frequency entries
  void mergeEntries(size_t const i, size_t const j);      // Merges the i-th and j-th entries, i.e. deletes
                                                          // the one with the smallest frequency and sets
                                                          // the other's frequency to their total.
  void deleteEntry(size_t const i);                       // Deletes the i-th entry
  void clearEntries();                                    // Clears all the table entries
  void clear();                                           // Clears the table

// Note: In mergeEntries, the resulting entry is stored in the i-th position, and 
// the j-th position is erased.

// Accessors

  size_t const numEntries() const;                                // Returns the number of entries
  size_t const numVariables() const;                              // Returns the number of variables
  FrequencyEntry const &entry(size_t const index) const;          // Returns the index-th table entry
  DistributionVariable const &variable(size_t const index) const; // Returns the range and number of bins
                                                                  // of the index-th variable

// Output

  friend std::ostream& operator<<(std::ostream &str, FrequencyTable const &table);

private:

  std::vector<DistributionVariable> variables_; // The types, ranges and bin numbers of each variable
  std::vector<FrequencyEntry> entries_;         // Frequency entries
};

/*
** End of class FrequencyTable
*/

// Inlines

inline FrequencyTable::FrequencyTable()
:
variables_(), entries_()
{}

inline void FrequencyTable::addVariable(DistributionVariable const &var)
{
  entries_.clear();
  variables_.push_back(var);
}

inline void FrequencyTable::deleteVariable(size_t const index)
{
  assert(index < variables_.size());
  entries_.clear();
  std::vector<DistributionVariable>::iterator it = variables_.begin();
  for(size_t i = 0; i < index; i++, it++);
  variables_.erase(it);
}
 
inline void FrequencyTable::addValue(std::vector<BinIndices> const &binIndices)
{
  if(entries_.size() > 0 && binIndices.size() != entries_[0].indices.size())
  {
    std::cerr << "Error in FrequencyTable::addValue: " << std::endl
              << "New entry is inconsistent" << std::endl;
    return;
  }
  // Search table for entry
  bool foundEntry = false; size_t i = 0;
  while(!foundEntry && i < entries_.size())
  {
    if(entries_[i].indices == binIndices) foundEntry = true;
    if(!foundEntry) i++;
  }
  if(foundEntry) entries_[i].frequency++;
  else entries_.push_back(FrequencyEntry(binIndices, 1));
}

/*
// See note above about this method
inline void FrequencyTable::addValue(size_t const num, ...)
{
  if(num != variables_.size())
  {
    std::cerr << "Error in FrequencyTable::addValue: " << std::endl
              << "Wrong number of parameters" << std::endl;
    return;
  }

  std::vector<BinIndices> indices;  // Indices to add to the table

  va_list argList;        // List of function arguments
  va_start(argList, num); // Initialize argument list

  for(size_t i = 0; i < num; i++)
  {
    switch(variables_[i].type)
    {
      case(LINEAR):
        {
          Real linearValue = va_arg(argList, Real);
          indices.push_back(variables_[i].linearIndex(linearValue));
        }
        break;
      case(UNIT_VECTOR_3D):
        {
          Vector3D vectorValue = va_arg(argList, Vector3D);
          indices.push_back(variables_[i].unitVector3DIndices(vectorValue));
        }
        break;
      case(UNIT_QUATERNION):
        {
          Quaternion quaternionValue = va_arg(argList, Quaternion);
          indices.push_back(variables_[i].unitQuaternionIndices(quaternionValue));
        }
        break;
      default:
        std::cerr << "Internal error in FrequencyTable::addValue!" << std::endl;
        return;
    }
  }

  va_end(argList);  // Close argument list
  addValue(indices);
}
*/

inline void FrequencyTable::setVariable(size_t const index, DistributionVariable const &var)
{
  assert(index < variables_.size());
  variables_[index] = var;
  entries_.clear();
}

inline void FrequencyTable::mergeEntries(size_t const i, size_t const j)
{
  assert(i < entries_.size() && j < entries_.size());
  size_t jFreq = entries_[j].frequency;
  if(jFreq > entries_[i].frequency) entries_[i].indices = entries_[j].indices;
  entries_[i].frequency += jFreq;
  std::vector<FrequencyEntry>::iterator it = entries_.begin();
  for(size_t k = 0; k < j; k++) it++;
  entries_.erase(it);
}

inline void FrequencyTable::deleteEntry(size_t const i)
{
  assert(i < entries_.size());
  std::vector<FrequencyEntry>::iterator it = entries_.begin();
  for(size_t j = 0; j < i; j++) it++;
  entries_.erase(it);
}

inline void FrequencyTable::clearEntries()
{
  entries_.clear();
}

inline void FrequencyTable::clear()
{
  variables_.clear();
  entries_.clear();
}

inline size_t const FrequencyTable::numEntries() const
{
  return entries_.size();
}

inline size_t const FrequencyTable::numVariables() const
{
  return variables_.size();
}

inline FrequencyEntry const &FrequencyTable::entry(size_t const index) const
{
  assert(index < entries_.size());
  return entries_[index];
}

inline DistributionVariable const &FrequencyTable::variable(size_t const index) const
{
  assert(index < variables_.size());
  return variables_[index];
}

inline std::ostream& operator<<(std::ostream &str, FrequencyTable const &table)
{
  str << std::endl << "Variables: ";
  if(table.numVariables() == 0) str << "None.";
  for(size_t i = 0; i < table.numVariables(); i++)
  {
    str << std::endl << "Index: " << i << table.variable(i);
  }
  str << std::endl << "Entries: ";
  if(table.numEntries() == 0) str << "None.";
  for(size_t i = 0; i < table.numEntries(); i++)
  {
    str << std::endl;
    //for(size_t j = 0; j < table.entry(i).indices.size(); j++)
    //{
    //  std::cout << table.entry(i).indices[j] << " ";
    //}
    for(size_t j = 0; j < table.numVariables(); j++)
    {
      switch(table.variable(j).type)
      {
        case(LINEAR):
          str << table.variable(j).linearValue(table.entry(i).indices[j]) << " ";
          break;
        case(UNIT_VECTOR_3D):
          str << table.variable(j).unitVector3DValue(table.entry(i).indices[j]) << " ";
          break;
        case(UNIT_QUATERNION):
          str << table.variable(j).unitQuaternionValue(table.entry(i).indices[j]) << " ";
          break;
        default:
          std::cerr << "Internal error in operator<<(std::ostream &str, FrequencyTable const &table!"
                    << std::endl;
          return str;
      }
    }
    str << table.entry(i).frequency;
  }

  return str;
}

#endif
