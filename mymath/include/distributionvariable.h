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
** Distribution variable
**
** This is used to do statistics on continuous variables. A distribution
** variable can be linear, a vector, a quaternion, etc. This class
** implements the methods to discretize such variables and convert back
** and forth between bin indices and the continuous representation.
**
** Note that all the data members are public. Unit vectors and unit 
** quaternions are patched onto a cube/tesseract. More types can be added
** by defining the appropriate specialized ...Value() and ...Indices()
** methods.
**
** Notes: 
**
** - Is there a more elegant/better way to do this? Maybe using a 
**   variadic template?
** - For future versions, it would be better to implement a vector,
**   quaternion, etc. grid class and use that to convert between indices
**   and values.
*/

#ifndef H_DISTRIBUTION_VARIABLE
#define H_DISTRIBUTION_VARIABLE

#include <iostream>
#include <cmath>
#include <vector>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/quaternion.h"

enum VariableType { LINEAR, UNIT_VECTOR_3D, UNIT_QUATERNION };  // Add more as needed

// I/O for variable types - add more as needed

inline std::ostream &operator<<(std::ostream &str, VariableType const &var)
{
  switch(var)
  {
    case(LINEAR):
      str << "LINEAR";
      break;
    case(UNIT_VECTOR_3D):
      str << "UNIT_VECTOR_3D";
      break;
    case(UNIT_QUATERNION):
      str << "UNIT_QUATERNION";
      break;
  }
  return str;
}

typedef std::vector<size_t> BinIndices;

// I/O for bin indices

inline std::ostream &operator<<(std::ostream &str, BinIndices const &indices)
{
  for(size_t i = 0; i < indices.size(); i++)
    str << indices[i] << " ";
  return str;
}

class DistributionVariable
{
public:

  VariableType type;  // Type of variable
  Real minValue;      // Minimum value (not meaningful for all types)
  Real maxValue;      // Maximum value (not meaningful for all types)
  size_t numBins;     // Number of bins (in each dimension)

// Constructors

  DistributionVariable();                                   // Define an "empty" variable
  DistributionVariable(Real const &minValue,                // Define by giving the minimum
                       Real const &maxValue,                // and maximum values, the number of
                       size_t const numBins,                // bins and, optionally, the type
                       VariableType const &type = LINEAR);
  DistributionVariable(VariableType const &type,            // Define by givind the variable 
                       size_t const numBins,                // type,the number of bins, and,
                       Real const &minValue = -1.0,         // optionally, the minimum and 
                       Real const &maxValue = 1.0);         // maximum values

// Interface - add more as needed
  
  Real const linearValue(BinIndices const &indices) const;                // Returns the real value associated
                                                                          // with the given index
  Vector3D const unitVector3DValue(BinIndices const &indices) const;      // Returns the 3D unit vector associated
                                                                          // with the given indices
  Quaternion const unitQuaternionValue(BinIndices const &indices) const;  // Returns the unit quaternion associated
                                                                          // with the given indices
  BinIndices const linearIndex(Real const &value) const;                  // Returns the index associated
                                                                          // with the given real value
  BinIndices const unitVector3DIndices(Vector3D const &v) const;          // Returns the indices associated
                                                                          // with the given 3D unit vector
  BinIndices const unitQuaternionIndices(Quaternion const &q) const;      // Returns the indices associated
                                                                          // with the given unit quaternion

// Output

  friend std::ostream &operator<<(std::ostream &str, DistributionVariable const &var);
};

/*
** End of class DistributionVariable
*/

// Inlines

inline DistributionVariable::DistributionVariable()
:
minValue(0.0), maxValue(0.0), numBins(0)
{}

inline DistributionVariable::DistributionVariable(Real const &minValue, 
                                                  Real const &maxValue, 
                                                  size_t const numBins,
                                                  VariableType const &type)
:
type(type), minValue(minValue), maxValue(maxValue), numBins(numBins)
{}

inline DistributionVariable::DistributionVariable(VariableType const &type,
                                                  size_t const numBins,
                                                  Real const &minValue, 
                                                  Real const &maxValue)
:
type(type), minValue(minValue), maxValue(maxValue), numBins(numBins)
{}

// Interface

inline Real const DistributionVariable::linearValue(BinIndices const &indices) const
{
  assert(type == LINEAR);
  assert(indices.size() == 1);
  if(indices[0] >= numBins) return maxValue;
  if(numBins <= 1) return (minValue + maxValue)/2.0;
  return minValue + ((Real)indices[0] + 0.5)*(maxValue - minValue)/((Real)numBins);
}

inline BinIndices const DistributionVariable::linearIndex(Real const &value) const
{
  assert(type == LINEAR);
  if(value > maxValue) return BinIndices(1, numBins - 1);
  if(value < minValue) return BinIndices(1, 0);
  if(numBins <= 1) return BinIndices(1, 0);
  if(minValue == maxValue) return BinIndices(1, 0);
  size_t index = (size_t)floor((Real)numBins*(value - minValue)/(maxValue - minValue));
  if(index > numBins - 1) index = numBins - 1;
  return BinIndices(1, index);
}

inline Vector3D const DistributionVariable::unitVector3DValue(BinIndices const &indices) const
{
  assert(type == UNIT_VECTOR_3D);
  assert(indices.size() == 3);
  if(numBins <= 1) return 0.0;
  size_t faceIndex = indices[0];  // Get the cube face
  Vector3D v = 0.0;
  Real deltav = 2.0/(Real)numBins;
  switch(faceIndex)
  {
    case(0):
      v.x = 1.0;
      v.y = ((Real)indices[1] + 0.5)*deltav - 1.0;
      v.z = ((Real)indices[2] + 0.5)*deltav - 1.0;
      break;
    case(1):
      v.x = ((Real)indices[1] + 0.5)*deltav - 1.0;
      v.y = 1.0;
      v.z = ((Real)indices[2] + 0.5)*deltav - 1.0;
      break;
    case(2):
      v.x = ((Real)indices[1] + 0.5)*deltav - 1.0;
      v.y = ((Real)indices[2] + 0.5)*deltav - 1.0;
      v.z = 1.0;
      break;
    case(3):
      v.x = -1.0;
      v.y = ((Real)indices[1] + 0.5)*deltav - 1.0;
      v.z = ((Real)indices[2] + 0.5)*deltav - 1.0;
      break;
    case(4):
      v.x = ((Real)indices[1] + 0.5)*deltav - 1.0;
      v.y = -1.0;
      v.z = ((Real)indices[2] + 0.5)*deltav - 1.0;
      break;
    case(5):
      v.x = ((Real)indices[1] + 0.5)*deltav - 1.0;
      v.y = ((Real)indices[2] + 0.5)*deltav - 1.0;
      v.z = -1.0;
      break;
    default:
      std::cerr << "Internal error in DistributionVariable::unitVector3DValue!" << std::endl;
      return 0.0;
  }
  v.normalize();
  return v;
}

inline BinIndices const DistributionVariable::unitVector3DIndices(Vector3D const &v) const
{
  assert(type == UNIT_VECTOR_3D);
  BinIndices indices(3, 0);
  size_t faceIndex = (v.x > 0.0)?0:3; // Defines the cube face
  Real vMax = fabs(v.x);
  if(fabs(v.y) > vMax)
  {
    vMax = fabs(v.y);
    faceIndex = (v.y > 0.0)?1:4;
  }
  if(fabs(v.z) > vMax)
  {
    vMax = fabs(v.z);
    faceIndex = (v.z > 0.0)?2:5;
  }
  Real deltav = 2.0/(Real)numBins;
  
  size_t v1Index, v2Index;  // Indices within face of cubic patch

  switch(faceIndex)
  {
    case(0):
    case(3):
      v1Index = (size_t)floor((v.y/vMax + 1.0)/deltav);
      v2Index = (size_t)floor((v.z/vMax + 1.0)/deltav);
      break;
    case(1):
    case(4):
      v1Index = (size_t)floor((v.x/vMax + 1.0)/deltav);
      v2Index = (size_t)floor((v.z/vMax + 1.0)/deltav);
      break;
    case(2):
    case(5):
      v1Index = (size_t)floor((v.x/vMax + 1.0)/deltav);
      v2Index = (size_t)floor((v.y/vMax + 1.0)/deltav);
      break;
    default:
      std::cerr << "Internal error in DistributionVariable::unitVector3DIndices!" << std::endl;
      return BinIndices(3, 0);
  }
  indices[0] = faceIndex; 
  indices[1] = (v1Index < numBins - 1)?v1Index:(numBins - 1);
  indices[2] = (v2Index < numBins - 1)?v2Index:(numBins - 1);
  return indices;
} 

inline Quaternion const DistributionVariable::unitQuaternionValue(BinIndices const &indices) const
{
  assert(type == UNIT_QUATERNION);
  assert(indices.size() == 4);
  if(numBins <= 1) return 0.0;
  size_t faceIndex = indices[0];  // Get the tesseract face
  Quaternion q = 0.0;
  Real deltaq = 2.0/(Real)numBins;
  switch(faceIndex)
  {
    case(0):
      q.w = 1.0;
      q.x = ((Real)indices[1] + 0.5)*deltaq - 1.0;
      q.y = ((Real)indices[2] + 0.5)*deltaq - 1.0;
      q.z = ((Real)indices[3] + 0.5)*deltaq - 1.0;
      break;
    case(1):
      q.w = ((Real)indices[1] + 0.5)*deltaq - 1.0;
      q.x = 1.0;
      q.y = ((Real)indices[2] + 0.5)*deltaq - 1.0;
      q.z = ((Real)indices[3] + 0.5)*deltaq - 1.0;
      break;
    case(2):
      q.w = ((Real)indices[1] + 0.5)*deltaq - 1.0;
      q.x = ((Real)indices[2] + 0.5)*deltaq - 1.0;
      q.y = 1.0;
      q.z = ((Real)indices[3] + 0.5)*deltaq - 1.0;
      break;
    case(3):
      q.w = ((Real)indices[1] + 0.5)*deltaq - 1.0;
      q.x = ((Real)indices[2] + 0.5)*deltaq - 1.0;
      q.y = ((Real)indices[3] + 0.5)*deltaq - 1.0;
      q.z = 1.0;
      break;
    default:
      std::cerr << "Internal error in DistributionVariable::unitQuaternionValue!" << std::endl;
      return 0.0;
  }
  q.normalize();
  return q;
}

inline BinIndices const DistributionVariable::unitQuaternionIndices(Quaternion const &q) const
{
  assert(type == UNIT_QUATERNION);
  BinIndices indices(4, 0);
  size_t faceIndex = 0; // Defines the tesseract face
  Real qMax = q.w;
  if(fabs(q.x) > fabs(qMax))
  {
    qMax = q.x;
    faceIndex = 1;
  }
  if(fabs(q.y) > fabs(qMax))
  {
    qMax = q.y;
    faceIndex = 2;
  }
  if(fabs(q.z) > fabs(qMax))
  {
    qMax = q.z;
    faceIndex = 3;
  }
  Real deltaq = 2.0/(Real)numBins;
  
  size_t q1Index, q2Index, q3Index; // Indices within face of tesseract

  switch(faceIndex)
  {
    case(0):
      q1Index = (size_t)floor((q.x/qMax + 1.0)/deltaq);
      q2Index = (size_t)floor((q.y/qMax + 1.0)/deltaq);
      q3Index = (size_t)floor((q.z/qMax + 1.0)/deltaq);
      break;
    case(1):
      q1Index = (size_t)floor((q.w/qMax + 1.0)/deltaq);
      q2Index = (size_t)floor((q.y/qMax + 1.0)/deltaq);
      q3Index = (size_t)floor((q.z/qMax + 1.0)/deltaq);
      break;
    case(2):
      q1Index = (size_t)floor((q.w/qMax + 1.0)/deltaq);
      q2Index = (size_t)floor((q.x/qMax + 1.0)/deltaq);
      q3Index = (size_t)floor((q.z/qMax + 1.0)/deltaq);
      break;
    case(3):
      q1Index = (size_t)floor((q.w/qMax + 1.0)/deltaq);
      q2Index = (size_t)floor((q.x/qMax + 1.0)/deltaq);
      q3Index = (size_t)floor((q.y/qMax + 1.0)/deltaq);
      break;
    default:
      std::cerr << "Internal error in DistributionVariable::unitQuaternionIndices!" << std::endl;
      return BinIndices(4, 0);
  }

  indices[0] = faceIndex; 
  indices[1] = (q1Index < numBins - 1)?q1Index:(numBins - 1);
  indices[2] = (q2Index < numBins - 1)?q2Index:(numBins - 1);
  indices[3] = (q3Index < numBins - 1)?q3Index:(numBins - 1);
  return indices;
}

inline std::ostream &operator<<(std::ostream &str, DistributionVariable const &var)
{
  str << std::endl << "Variable type: " << var.type;
  if(var.type == LINEAR)
  {
    str << std::endl << "Minimum value: " << var.minValue
        << std::endl << "Maximum value: " << var.maxValue;
  }
  str << std::endl << "Number of bins: " << var.numBins;
  return str;
}

#endif
