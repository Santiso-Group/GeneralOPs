/*
** Copyright 2008, 2009 Erik Santiso.
** This file is part of crystdist.
** crystdist is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
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
*/

// Uncomment to implement the water temporary fix
#define WATER_FIX

#include "mymath/include/mpconstants.h"
#include "crystdist/include/globalhistogram.h"

// Constructors

GlobalHistogram::GlobalHistogram()
:
entryTypes_(GENERAL, GENERAL), entryLabel_(), settings_(), frequencyTable_(), internalDOFInfo_()
{}

// Interface

void GlobalHistogram::addPoint(RelativeConfiguration const &conf)
{
  if(conf.distance > settings_.cutoff) return;
  if(frequencyTable_.numEntries() == 0)
  {
    // Store entry types and label
    entryTypes_ = conf.types;
    entryLabel_ = conf.name;

    // Set variable types 
    if(entryTypes_.first == GENERAL && entryTypes_.second == GENERAL)
    {
      // General case: Bond orientation = Vector3D, Relative orientation = Quaternion
      frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numDistanceBins, 0.0, settings_.cutoff)); 
      frequencyTable_.addVariable(DistributionVariable(UNIT_VECTOR_3D, settings_.numVectorBins)); 
      frequencyTable_.addVariable(DistributionVariable(UNIT_QUATERNION, settings_.numQuaternionBins));
    }
    else if(entryTypes_.first == LINEAR_ASYMMETRIC && entryTypes_.second == LINEAR_ASYMMETRIC)
    {
      // Axisymmetric case: Bond orientation and relative orientation are both angles
      frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numDistanceBins, 0.0, settings_.cutoff)); 
      frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numAngleBins, 0.0, PI)); 
      frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numAngleBins, 0.0, PI));
    }
    else if((entryTypes_.first == LINEAR_SYMMETRIC && entryTypes_.second == LINEAR_SYMMETRIC) ||
            (entryTypes_.first == PLANAR_SYMMETRIC && entryTypes_.second == PLANAR_SYMMETRIC))
    {
      // Axisymmetric case with symmetry plane: Bond orientation and relative orientation are both angles
      frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numDistanceBins, 0.0, settings_.cutoff)); 
      frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numAngleBins, 0.0, PI/2.0)); 
      frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numAngleBins, 0.0, PI/2.0));
    }
    else
    {
      std::cerr << "Error in GlobalHistogram::addPoint: Combination of molecule types not implemented" << std::endl;
      clear();
      return;
    }

    // Store internal DOF info
    for(size_t i = 0; i < conf.internalDOFs.first.size(); i++)
    {
      InternalDOFDefinition currentDefinition;
      currentDefinition.type = conf.internalDOFs.first[i].type;
      currentDefinition.symmetryNumber = conf.internalDOFs.first[i].symmetryNumber;
      internalDOFInfo_.first.push_back(currentDefinition);
      switch(conf.internalDOFs.first[i].type)
      {
        case(DISTANCE):
          frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numInternalDistanceBins, 
                                                           0.0, settings_.internalCutoff));
          break;
        case(ANGLE):
          frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numInternalAngleBins, 0.0, PI));
          break;
        case(DIHEDRAL):
          {
          Real maxValue = TWO_PI/(Real)conf.internalDOFs.first[i].symmetryNumber;
          frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numInternalDihedralBins, 
                                                           0.0, maxValue));
          }
          break;
      }
    }
    for(size_t i = 0; i < conf.internalDOFs.second.size(); i++)
    {
      InternalDOFDefinition currentDefinition;
      currentDefinition.type = conf.internalDOFs.second[i].type;
      currentDefinition.symmetryNumber = conf.internalDOFs.second[i].symmetryNumber;
      internalDOFInfo_.second.push_back(currentDefinition);
      switch(conf.internalDOFs.second[i].type)
      {
        case(DISTANCE):
          frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numInternalDistanceBins, 
                                                           0.0, settings_.internalCutoff));
          break;
        case(ANGLE):
          frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numInternalAngleBins, 0.0, PI));
          break;
        case(DIHEDRAL):
          {
          Real maxValue = TWO_PI/(Real)conf.internalDOFs.second[i].symmetryNumber;
          frequencyTable_.addVariable(DistributionVariable(LINEAR, settings_.numInternalDihedralBins, 
                                                           0.0, maxValue));
          }
          break;
      }
    }
  }
  else
  {
    bool isConsistent = true;
    isConsistent &= (conf.types == entryTypes_);
    isConsistent &= (conf.name == entryLabel_);
    isConsistent &= (conf.internalDOFs.first.size() == internalDOFInfo_.first.size());
    if(isConsistent)
      for(size_t i = 0; i < internalDOFInfo_.first.size(); i++)
      {
        isConsistent &= (conf.internalDOFs.first[i].type == internalDOFInfo_.first[i].type);
        isConsistent &= (conf.internalDOFs.first[i].symmetryNumber == internalDOFInfo_.first[i].symmetryNumber);
      }
    isConsistent &= (conf.internalDOFs.second.size() == internalDOFInfo_.second.size());
    if(isConsistent)
      for(size_t i = 0; i < internalDOFInfo_.second.size(); i++)
      {
        isConsistent &= (conf.internalDOFs.second[i].type == internalDOFInfo_.second[i].type);
        isConsistent &= (conf.internalDOFs.second[i].symmetryNumber == internalDOFInfo_.second[i].symmetryNumber);
      }
    if(!isConsistent)
    {
      std::cerr << "Error in GlobalHistogram::addPoint: Configuration added is inconsistent "
                << std::endl << "with data already on the histogram." << std::endl;
      return;
    }
  }

  std::vector<BinIndices> indices;  // Indices for frequency entry

  // Find bin indices for distance, orientation and relative orientation
  if(entryTypes_.first == GENERAL && entryTypes_.second == GENERAL)
  {
    indices.push_back(frequencyTable_.variable(0).linearIndex(conf.distance));
    indices.push_back(frequencyTable_.variable(1).unitVector3DIndices(vector(conf.bondOrientation)));
    indices.push_back(frequencyTable_.variable(2).unitQuaternionIndices(conf.relativeOrientation));
  }
  else if((entryTypes_.first == LINEAR_ASYMMETRIC && entryTypes_.second == LINEAR_ASYMMETRIC) ||
          (entryTypes_.first == LINEAR_SYMMETRIC && entryTypes_.second == LINEAR_SYMMETRIC) ||
          (entryTypes_.first == PLANAR_SYMMETRIC && entryTypes_.second == PLANAR_SYMMETRIC))
  {
    indices.push_back(frequencyTable_.variable(0).linearIndex(conf.distance));
    indices.push_back(frequencyTable_.variable(1).linearIndex(scalar(conf.bondOrientation)));
    indices.push_back(frequencyTable_.variable(2).linearIndex(scalar(conf.relativeOrientation)));
  }
  else  // Should never get here
  {
    std::cerr << "Error in GlobalHistogram::addPoint: Combination of molecule types not implemented" << std::endl;
    return;
  }
  // Find bin indices for internal DOFs
  size_t nextIndex = 3;
  for(size_t i = 0; i < conf.internalDOFs.first.size(); i++)
  {
    Real value = conf.internalDOFs.first[i].value;
    if(conf.internalDOFs.first[i].type == ANGLE)
    {
      value -= TWO_PI*floor(value/TWO_PI);
      if(value > PI) value = TWO_PI - value;  // Impose proper symmetry
    }
    else if(conf.internalDOFs.first[i].type == DIHEDRAL)
    {
      Real sigma = (Real)conf.internalDOFs.first[i].symmetryNumber;
      value -= TWO_PI/sigma*floor(value*sigma/TWO_PI);  // Impose proper symmetry
    }
    indices.push_back(frequencyTable_.variable(i + nextIndex).linearIndex(value));
  }
  nextIndex = indices.size();
  for(size_t i = 0; i < conf.internalDOFs.second.size(); i++)
  {
    Real value = conf.internalDOFs.second[i].value;
    if(conf.internalDOFs.second[i].type == ANGLE)
    {
      value -= TWO_PI*floor(value/TWO_PI);
      if(value > PI) value = TWO_PI - value;  // Impose proper symmetry
    }
    else if(conf.internalDOFs.second[i].type == DIHEDRAL)
    {
      Real sigma = (Real)conf.internalDOFs.second[i].symmetryNumber;
      value -= TWO_PI/sigma*floor(value*sigma/TWO_PI);  // Impose proper symmetry
    }
    indices.push_back(frequencyTable_.variable(i + nextIndex).linearIndex(value));
  }
  frequencyTable_.addValue(indices);
}

void GlobalHistogram::addSystem(System<PointMolecule> const &system, std::string const &name)
{
  Real cutoff2 = settings_.cutoff; cutoff2 *= cutoff2;
  for(size_t i = 0; i < system.size(); i++) // Loop over molecules
  {
    if(system[i].name.find(name) == system[i].name.npos) 
      continue;  // Skip if not the named molecule

    for(size_t j = i + 1; j < system.size(); j++) // Loop over neighbors
    {
      if(system[j].name.find(name) == system[j].name.npos)
        continue; // Skip if not the named molecule
      if(norm2(system.lattice().difference(system[i].position, system[j].position)) > cutoff2) 
        continue; // Skip if distance > cutoff
      addPoint(RelativeConfiguration(system[i], system[j], system.lattice()));
      addPoint(RelativeConfiguration(system[j], system[i], system.lattice()));
    }
  }
}

void GlobalHistogram::reduce(Real const &normCutoff, size_t const freqCutoff)
{
  // Merge similar entries
  Real const cutoffSq = normCutoff*normCutoff;
  size_t i = 0;
  while(i < frequencyTable_.numEntries())
  {
    size_t j = i + 1;
    while(j < frequencyTable_.numEntries())
    {
      Real norm2 = 0.0;
      for(size_t k = 0; k < frequencyTable_.numVariables(); k++)
      {
        switch(frequencyTable_.variable(k).type)
        {
          case(LINEAR):
          {
            Real value1 = frequencyTable_.variable(k).linearValue(frequencyTable_.entry(i).indices[k]);
            Real value2 = frequencyTable_.variable(k).linearValue(frequencyTable_.entry(j).indices[k]);
            Real diff = value2 - value1;
            norm2 += diff*diff; 
          }
          break;
          case(UNIT_VECTOR_3D):
          {
            Vector3D value1 = frequencyTable_.variable(k).unitVector3DValue(frequencyTable_.entry(i).indices[k]);
            Vector3D value2 = frequencyTable_.variable(k).unitVector3DValue(frequencyTable_.entry(j).indices[k]);
#ifdef WATER_FIX            
            // Temporary fix for water 01/25/10
            Real dotProduct = value1.x*value2.x + fabs(value1.y*value2.y) + fabs(value1.z*value2.z);
#else
            Real dotProduct = value1*value2;
#endif
            // Erase the two lines above to remove the fix
            norm2 += 1.0 - dotProduct; // 2*sin^2(theta/2)
          }
          break;
          case(UNIT_QUATERNION):
          {
            Quaternion value1 = frequencyTable_.variable(k).unitQuaternionValue(frequencyTable_.entry(i).indices[k]);
            Quaternion value2 = frequencyTable_.variable(k).unitQuaternionValue(frequencyTable_.entry(j).indices[k]);
            Real dotProduct = dot(value1, value2);
#ifdef WATER_FIX
            // Temporary fix for water 01/25/10
            Quaternion q_i(0.0, 1.0, 0.0, 0.0);
            Real dotProduct2 = dot(-q_i*value1, value2);
            Real dotProduct3 = dot(value1*q_i, value2);
            if(fabs(dotProduct2) > fabs(dotProduct)) dotProduct = dotProduct2;
            if(fabs(dotProduct3) > fabs(dotProduct)) dotProduct = dotProduct3;
#endif
            // Erase the 5 lines above to remove the fix
            norm2 += 1.0 - fabs(dotProduct);  // 1 - |cos(theta/2)|
          }
          break;
        }
      }
      if(norm2 < cutoffSq) frequencyTable_.mergeEntries(i, j);
      else j++;
    }
    i++;
  }
  // Remove low-frequency entries
  if(freqCutoff > 0)
  {
    size_t i = 0;
    while(i < frequencyTable_.numEntries())
    {
      if(frequencyTable_.entry(i).frequency < freqCutoff) frequencyTable_.deleteEntry(i);
      else i++;
    }
  }
}

RelativeConfiguration const GlobalHistogram::entry(size_t const i) const
{
  RelativeConfiguration conf;
  // Molecule types and labels
  conf.types = entryTypes_;
  conf.name = entryLabel_;
  // Distance
  conf.distance = frequencyTable_.variable(0).linearValue(frequencyTable_.entry(i).indices[0]);
  // Bond orientation
  if(frequencyTable_.variable(1).type == LINEAR)
    conf.bondOrientation = frequencyTable_.variable(1).linearValue(frequencyTable_.entry(i).indices[1]);
  else if(frequencyTable_.variable(1).type == UNIT_VECTOR_3D)
    conf.bondOrientation = frequencyTable_.variable(1).unitVector3DValue(frequencyTable_.entry(i).indices[1]);
  else
  {
    std::cerr << "Internal error in GlobalHistogram::entry: Bad bond orientation data." << std::endl;
    conf.clear();
    return conf;
  }
  // Relative orientation
  if(frequencyTable_.variable(2).type == LINEAR)
    conf.relativeOrientation = frequencyTable_.variable(2).linearValue(frequencyTable_.entry(i).indices[2]);
  else if(frequencyTable_.variable(2).type == UNIT_QUATERNION)
    conf.relativeOrientation = frequencyTable_.variable(2).unitQuaternionValue(frequencyTable_.entry(i).indices[2]);
  else
  {
    std::cerr << "Internal error in GlobalHistogram::entry: Bad relative orientation data." << std::endl;
    conf.clear();
    return conf;
  }
  size_t currentIndex = 3;
  // Internal DOFs
  for(size_t j = 0; j < internalDOFInfo_.first.size(); j++)
  {
    InternalDOF currentDOF;
    currentDOF.type = internalDOFInfo_.first[j].type;
    currentDOF.symmetryNumber = internalDOFInfo_.first[j].symmetryNumber;
    currentDOF.value = frequencyTable_.variable(currentIndex + j).linearValue(frequencyTable_.entry(i).indices[currentIndex + j]);
    conf.internalDOFs.first.push_back(currentDOF);
  }
  currentIndex = 3 + internalDOFInfo_.first.size();
  for(size_t j = 0; j < internalDOFInfo_.second.size(); j++)
  {
    InternalDOF currentDOF;
    currentDOF.type = internalDOFInfo_.second[j].type;
    currentDOF.symmetryNumber = internalDOFInfo_.second[j].symmetryNumber;
    currentDOF.value = frequencyTable_.variable(currentIndex + j).linearValue(frequencyTable_.entry(i).indices[currentIndex + j]);
    conf.internalDOFs.second.push_back(currentDOF);
  }
  return conf;
}

// Print histogram data

std::ostream &operator<<(std::ostream &outStream, GlobalHistogram const &hist)
{
  // Write table header
  outStream << "Distance ";
  if(hist.frequencyTable_.variable(1).type == UNIT_VECTOR_3D)
    outStream << "v.x v.y v.z ";
  else
    outStream << "Bond_Angle ";
  if(hist.frequencyTable_.variable(2).type == UNIT_QUATERNION)
    outStream << "q.w q.x q.y q.z ";
  else
    outStream << "Relative_Angle ";
  if(hist.frequencyTable_.numVariables() > 3)
    outStream << "Internals(1) Internals(2) ";
  outStream << "Frequency";
  // Write table
  if(hist.frequencyTable_.numEntries() == 0) outStream << "No entries.";
  for(size_t i = 0; i < hist.frequencyTable_.numEntries(); i++)
  {
    outStream << std::endl;

    for(size_t j = 0; j < hist.frequencyTable_.numVariables(); j++)
    {
      switch(hist.frequencyTable_.variable(j).type)
      {
        case(LINEAR):
          outStream << hist.frequencyTable_.variable(j).linearValue(hist.frequencyTable_.entry(i).indices[j])
                    << " ";
          break;
        case(UNIT_VECTOR_3D):
          outStream << hist.frequencyTable_.variable(j).unitVector3DValue(hist.frequencyTable_.entry(i).indices[j]) 
                    << " ";
          break;
        case(UNIT_QUATERNION):
          {
          Quaternion q = hist.frequencyTable_.variable(j).unitQuaternionValue(hist.frequencyTable_.entry(i).indices[j]);
          outStream << q.w << " " << q.x << " " << q.y << " " << q.z << " ";
          }
          break;
      }
    }
    outStream << hist.frequencyTable_.entry(i).frequency;
  }
  return outStream;
}

#undef WATER_FIX

