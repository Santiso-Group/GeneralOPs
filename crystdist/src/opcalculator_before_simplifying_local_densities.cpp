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
** Crystal Order Parameters class
**
** Notes:
**
** - The implementation of the "per-molecule" case is not very elegant,
**   and it wastes some memory. It might be a good idea to change it for
**   a later version.
** - Note that p-norms are not applied to local densities
*/

#include "crystdist/include/opcalculator.h"

// Constructors

OrderParameterCalculator::OrderParameterCalculator()
:
ops_()
{}

OrderParameterCalculator::OrderParameterCalculator(System<PointMolecule> const &system, 
                                                   CrystalDistributionParameters const &parameters, 
                                                   OrderParameterGrid const &grid,
                                                   Real const &switchWidth,
                                                   Real const &cutoff, 
                                                   bool const useLogs,
                                                   Real const &pValue)
:
ops_()
{ 
  ops_.switchWidth = switchWidth;
  ops_.cutoffSq = cutoff*cutoff;
  ops_.grid = grid;
  calculate(system, parameters, useLogs, pValue); 
}

// Interface

void OrderParameterCalculator::calculate(System<PointMolecule> const &system, 
                                         CrystalDistributionParameters const &parameters,
                                         bool const useLogs,
                                         Real const &pValue)
{
  // Note that this assumes (without verifying) that the system is homogeneous (e.g. internal DOFs are
  // the same for all molecules)
  if(system.size() < 1)
  {
    std::cerr << "Error in CrystalOrderParameters::calculate(): Empty system" << std::endl;
    return;
  }

  // Store number of groups
  size_t const numGroups = parameters.means.size();
  if(numGroups < 1)
  {
    std::cerr << "Error in CrystalOrderParameters::calculate: Empty parameter set" << std::endl;
    clearOPs();
    return;
  }

  // Sanity check for p-norms
  if(pValue > 0.0 && pValue < 1.0)
  {
    std::cerr << "Error in CrystalOrderParameters::calculate: p value must be at least 1" << std::endl;
    clearOPs();
    return;
  }

  // Initialize variables
  clearOPs();

  size_t const numCells = ops_.grid.numCells();
  bool const useGrid = (numCells > 0);
  std::vector<Real> totalWeights((numCells > 0)?numCells:system.size(), 0);
  size_t const numIntDOFs = system[0].numInternalDOFs();
  if(useGrid) ops_.resize(numCells, numIntDOFs);
  else ops_.resize(system.size(), numIntDOFs);  // Do OPs per molecule
  
  // Main loop over molecules
  for(size_t i = 0; i < system.size(); ++i)
  {
    // Find cell index and update number of molecules
    CellData iData;
    if(useGrid) iData = cellData(system[i].position, system.lattice());
    else  // Per-molecule case
    {
      iData.clear();
      iData.cellIndices.push_back(i);
      iData.cellWeights.push_back(1.0);
    }
    size_t const iDataSize = iData.size();
    
    for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
      totalWeights[iData.cellIndices[iIndex]] += iData.cellWeights[iIndex];
    
    for(size_t j = i + 1; j < system.size(); ++j) // Loop over neighbors
    {
      // Check that neighbor is within cutoff
      if(ops_.cutoffSq && 
         norm2(system.lattice().difference(system[i].position, system[j].position)) > ops_.cutoffSq)
        continue;
      // Find cell index and relative configurations
      CellData jData;
      if(useGrid) jData = cellData(system[j].position, system.lattice());
      else
      {
        jData.clear();
        jData.cellIndices.push_back(j);
        jData.cellWeights.push_back(1.0);
      }
      size_t const jDataSize = jData.size();
      RelativeConfiguration const confij(system[i], system[j], system.lattice());
      RelativeConfiguration const confji(system[j], system[i], system.lattice());
      // Loop over groups and peaks
      for(size_t iGroup = 0; iGroup < numGroups; ++iGroup)
      {
        size_t numPeaks = parameters.means[iGroup].size();
        for(size_t iPeak = 0; iPeak < numPeaks; ++iPeak)
        {
          // Distance OPs
          Real deltaDistance = confij.distance - parameters.means[iGroup][iPeak].distance;
          Real distanceOP = parameters.normalizationFactors[iGroup][iPeak].distance*
                            exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].distance*
                                deltaDistance*deltaDistance);
          if(pValue > 0.0)
          {
            for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
              ops_.distance[iData.cellIndices[iIndex]] += 
                pow(iData.cellWeights[iIndex]*distanceOP, pValue);
            for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
              ops_.distance[jData.cellIndices[jIndex]] += 
                pow(jData.cellWeights[jIndex]*distanceOP, pValue);
          }
          else
          {
            for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
              ops_.distance[iData.cellIndices[iIndex]] += iData.cellWeights[iIndex]*distanceOP;
            for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
              ops_.distance[jData.cellIndices[jIndex]] += jData.cellWeights[jIndex]*distanceOP;
          }

          // Bond orientation and relative orientation OPs
          if(confij.types.first == GENERAL && confij.types.second == GENERAL)
          {
            // General case: Bond orientation = Vector3D, Relative orientation = Quaternion
            Real bondDotij = vector(confij.bondOrientation)*
                             vector(parameters.means[iGroup][iPeak].bondOrientation);
            Real bondDotji = vector(confji.bondOrientation)*
                             vector(parameters.means[iGroup][iPeak].bondOrientation);
            Real bondOPij = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                bondDotij);
            Real bondOPji = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                bondDotji);
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*bondOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*bondOPji;
            }
            Real relDotij = dot(parameters.means[iGroup][iPeak].relativeOrientation,
                                confij.relativeOrientation);
            Real relDotji = dot(parameters.means[iGroup][iPeak].relativeOrientation,
                                confji.relativeOrientation);
            Real relOPij = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               relDotij*relDotij);
            Real relOPji = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               relDotji*relDotji);
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*relOPji;
            }
            // Total OPs added 07-14-10
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji;
            }
            // Local densities added 06-09-11
            for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
              ops_.localDensity[iData.cellIndices[iIndex]] += 
                iData.cellWeights[iIndex];
            for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
              ops_.localDensity[jData.cellIndices[jIndex]] +=
                jData.cellWeights[jIndex];
          }
          else if(confij.types.first == LINEAR_ASYMMETRIC && confij.types.second == LINEAR_ASYMMETRIC)
          {
            // Axisymmetric case: Bond orientation and relative orientation are both angles
            Real deltaBondij = confij.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real deltaBondji = confji.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real bondOPij = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(deltaBondij));
            Real bondOPji = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(deltaBondji));
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*bondOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*bondOPji;
            }
            Real deltaRelij = confij.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real deltaRelji = confji.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real relOPij = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(deltaRelij));
            Real relOPji = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(deltaRelji));
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*relOPji;
            }
            // Total OPs added 07-14-10
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji;
            }
            // Local densities added 06-09-11
            for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
              ops_.localDensity[iData.cellIndices[iIndex]] +=
                iData.cellWeights[iIndex];
            for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
              ops_.localDensity[jData.cellIndices[jIndex]] +=
                jData.cellWeights[jIndex];
          }
          else if((confij.types.first == LINEAR_SYMMETRIC && confij.types.second == LINEAR_SYMMETRIC) ||
                  (confij.types.first == PLANAR_SYMMETRIC && confij.types.second == PLANAR_SYMMETRIC))
          {
            // Axisymmetric case with symmetry plane: Bond orientation and relative orientation are both angles
            Real deltaBondij = confij.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real deltaBondji = confji.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real bondOPij = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(2.0*deltaBondij));  // 2 due to symmetry
            Real bondOPji = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(2.0*deltaBondji));  // 2 due to symmetry
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*bondOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*bondOPji;
            }

            Real deltaRelij = confij.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real deltaRelji = confji.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real relOPij = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(2.0*deltaRelij));  // 2 due to symmetry
            Real relOPji = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(2.0*deltaRelji));  // 2 due to symmetry
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*relOPji;
            }
            // Total OPs added 07-14-10
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji;
            }
            // Local densities added 06-09-11
            for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
              ops_.localDensity[iData.cellIndices[iIndex]] +=
                iData.cellWeights[iIndex];
            for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
              ops_.localDensity[jData.cellIndices[jIndex]] +=
                jData.cellWeights[jIndex];
          }
          else
          {
            std::cerr << "Error in CrystalOrderParameters::calculate(): "
                      << "Combination of molecule types not implemented" << std::endl;
            return;
          }
          // Internal DOF order parameters
          // Note that this assumes that both molecules have the same types and number
          // of internal DOFs

          for(size_t iDOF = 0; iDOF < numIntDOFs; ++iDOF)
          {
            switch(confij.internalDOFs.first[iDOF].type)
            {
              case DISTANCE:
              {
                // Gaussian
                Real deltaIntDOFij1 = confij.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFij2 = confij.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPij1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                         deltaIntDOFij1*deltaIntDOFij1);
                Real intOPij2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                         deltaIntDOFij2*deltaIntDOFij2);
                if(pValue > 0.0)
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      pow(iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2, pValue);
                }
                else
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2;
                }
                
                Real deltaIntDOFji1 = confji.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFji2 = confji.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPji1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                         deltaIntDOFji1*deltaIntDOFji1);
                Real intOPji2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                         deltaIntDOFji2*deltaIntDOFji2);
                if(pValue > 0.0)
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      pow(jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2, pValue);
                }
                else
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2;
                }
              }
                break;
              case ANGLE:
              {
                // von Mises
                Real deltaIntDOFij1 = confij.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFij2 = confij.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPij1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFij1));
                Real intOPij2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFij2));
                if(pValue > 0.0)
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      pow(iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2, pValue);
                }
                else
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2;
                }

                Real deltaIntDOFji1 = confji.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFji2 = confji.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPji1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFji1));
                Real intOPji2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFji2));
                if(pValue > 0.0)
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      pow(jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2, pValue);
                }
                else
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2;
                }
              }
                break;
              case DIHEDRAL:
              {
                // von Mises with symmetry factor
                Real deltaIntDOFij1 = confij.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                deltaIntDOFij1 *= (Real)confij.internalDOFs.first[iDOF].symmetryNumber;
                Real deltaIntDOFij2 = confij.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                deltaIntDOFij2 *= (Real)confij.internalDOFs.second[iDOF].symmetryNumber;
                Real intOPij1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFij1));
                Real intOPij2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFij2));
                if(pValue > 0.0)
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      pow(iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2, pValue);
                }
                else
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2;
                }

                Real deltaIntDOFji1 = confji.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                deltaIntDOFji1 *= (Real)confji.internalDOFs.first[iDOF].symmetryNumber;
                Real deltaIntDOFji2 = confji.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                deltaIntDOFji2 *= (Real)confji.internalDOFs.second[iDOF].symmetryNumber;
                Real intOPji1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFji1));
                Real intOPji2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFji2));
                if(pValue > 0.0)
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      pow(jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2, pValue);
                }
                else
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2;
                }
              }
                break;
            }
          }
        } // End of loop over peaks
      } // End of loop over groups
    } // End of loop over neighbors
  } // End of main loop over molecules
  
  // Normalize order parameters
  Real cellVolume = useGrid?(system.lattice().volume()/ops_.grid.numCells()):1.0;
  for(size_t i = 0; i < numCells; i++) // Won't run if numCells = 0 (per molecule case)
  {
    Real totalWeight = totalWeights[i];
    if(totalWeight == 0) continue;
    if(pValue > 0.0)
    {
      ops_.distance[i] = pow(ops_.distance[i], 1.0/pValue);
      ops_.bondOrientation[i] = pow(ops_.bondOrientation[i], 1.0/pValue);
      ops_.relativeOrientation[i] = pow(ops_.relativeOrientation[i], 1.0/pValue);
      for(size_t j = 0; j < ops_.internal[i].size(); ++j)
        ops_.internal[i][j] = pow(ops_.internal[i][j], 1.0/pValue);
      ops_.total[i] = pow(ops_.total[i], 1.0/pValue);
    }
    else
    {
      ops_.distance[i] /= totalWeight;
      ops_.bondOrientation[i] /= totalWeight;
      ops_.relativeOrientation[i] /= totalWeight;
      ops_.total[i] /= totalWeight;
      for(size_t j = 0; j < ops_.internal[i].size(); ++j)
        ops_.internal[i][j] /= totalWeight;
    }
    size_t totalPeaks = 0;
    for(size_t iGroup = 0; iGroup < numGroups; ++iGroup)
      totalPeaks += parameters.means[iGroup].size();
    std::cout << "Cell " << i << " : num = " << ops_.localDensity[i] << " ; totalPeaks = " << totalPeaks << " ; cellVol = " << cellVolume << std::endl; // TEST
    ops_.localDensity[i] /= (Real)totalPeaks*cellVolume;
  }

  // Switch to logs if needed
  if(useLogs)
  {
    for(size_t i = 0; i < ops_.distance.size(); ++i)
      ops_.distance[i] = (ops_.distance[i] > 0.0)?log(ops_.distance[i]):-99999;
    for(size_t i = 0; i < ops_.bondOrientation.size(); ++i)
      ops_.bondOrientation[i] = (ops_.bondOrientation[i] > 0.0)?log(ops_.bondOrientation[i]):-99999;
    for(size_t i = 0; i < ops_.relativeOrientation.size(); ++i)
      ops_.relativeOrientation[i] = (ops_.relativeOrientation[i] > 0.0)?log(ops_.relativeOrientation[i]):-99999;
    for(size_t i = 0; i < ops_.total.size(); ++i)
      ops_.total[i] = (ops_.total[i] > 0.0)?log(ops_.total[i]):-99999;
    for(size_t i = 0; i < ops_.localDensity.size(); ++i)
      ops_.localDensity[i] = (ops_.localDensity[i] > 0.0)?log(ops_.localDensity[i]):-99999;
    for(size_t i = 0; i < ops_.internal.size(); ++i)
      for(size_t j = 0; j < ops_.internal[i].size(); ++j)
        ops_.internal[i][j] = (ops_.internal[i][j] > 0.0)?log(ops_.internal[i][j]):-99999;
  }
}

// Private functions

CellData const OrderParameterCalculator::cellData(Vector3D const &position, Lattice const &lattice)
{
  CellData cellData;

  if(lattice.numDimensions() != 3)
  {
    std::cerr << "Error in OrderParameterCalculator::cellData - invalid number of dimensions " << std::endl;
    return cellData;
  }

  Vector3D const pos = position - ops_.grid.origin; // For compatibility with NAMD
  Vector3D scaled(lattice.reciprocalVector(0)*pos,
                  lattice.reciprocalVector(1)*pos,
                  lattice.reciprocalVector(2)*pos);

  // In case molecule has drifted out of simulation box
  scaled.x -= 1 + floor(scaled.x - 0.5);
  scaled.y -= 1 + floor(scaled.y - 0.5);
  scaled.z -= 1 + floor(scaled.z - 0.5);

  // Numbers of cells
  size_t const nx = ops_.grid.x;
  size_t const ny = ops_.grid.y;
  size_t const nz = ops_.grid.z;

  // Indices that contribute to average
  std::vector<size_t> indx, indy, indz;

  // Indices of nearest cell center
  size_t ixc = (size_t)floor(nx*(scaled.x + 0.5 - floor(scaled.x + 0.5)));
  size_t iyc = (size_t)floor(ny*(scaled.y + 0.5 - floor(scaled.y + 0.5)));
  size_t izc = (size_t)floor(nz*(scaled.z + 0.5 - floor(scaled.z + 0.5)));
  if(ixc > nx - 1) ixc = nx - 1;
  if(iyc > ny - 1) iyc = ny - 1;
  if(izc > nz - 1) izc = nz - 1;

  indx.push_back(ixc);
  indy.push_back(iyc);
  indz.push_back(izc);

  // Nearest cell center in scaled coordinates and difference
  Vector3D const center(((Real)ixc+0.5)/(Real)nx - 0.5,
                        ((Real)iyc+0.5)/(Real)ny - 0.5, 
                        ((Real)izc+0.5)/(Real)nz - 0.5);

  Vector3D delta = scaled - center;

  // Partial contributions to weights
  std::vector<Real> fx, fy, fz;

  Real scaledDiff = (ops_.switchWidth > 0.0)?(fabs(nx*delta.x) - 0.5)/ops_.switchWidth:-2.0;

  if(scaledDiff > -1.0)
  {
    if(delta.x > 0.0)
    {
      if(ixc == nx - 1) indx.push_back(0);
      else indx.push_back(ixc + 1);
    }
    else
    {
      if(ixc == 0) indx.push_back(nx - 1);
      else indx.push_back(ixc - 1);
    }
    fx.push_back(0.5 + 0.25*scaledDiff*(scaledDiff*scaledDiff - 3));
    fx.push_back(1.0 - fx[0]);
  }
  else
    fx.push_back(1.0);

  scaledDiff = (ops_.switchWidth > 0.0)?(fabs(ny*delta.y) - 0.5)/ops_.switchWidth:-2.0;

  if(scaledDiff > -1.0)
  {
    if(delta.y > 0.0)
    {
      if(iyc == ny - 1) indy.push_back(0);
      else indy.push_back(iyc + 1);
    }
    else
    {
      if(iyc == 0) indy.push_back(ny - 1);
      else indy.push_back(iyc - 1);
    }
    fy.push_back(0.5 + 0.25*scaledDiff*(scaledDiff*scaledDiff - 3));
    fy.push_back(1.0 - fy[0]);
  }
  else
    fy.push_back(1.0);

  scaledDiff = (ops_.switchWidth > 0.0)?(fabs(nz*delta.z) - 0.5)/ops_.switchWidth:-2.0;

  if(scaledDiff > -1.0)
  {
    if(delta.z > 0.0)
    {
      if(izc == nz - 1) indz.push_back(0);
      else indz.push_back(izc + 1);
    }
    else
    {
      if(izc == 0) indz.push_back(nz - 1);
      else indz.push_back(izc - 1);
    }
    fz.push_back(0.5 + 0.25*scaledDiff*(scaledDiff*scaledDiff - 3));
    fz.push_back(1.0 - fz[0]);
  }
  else
    fz.push_back(1.0);

  for(size_t inx = 0; inx < indx.size(); ++inx)
  for(size_t iny = 0; iny < indy.size(); ++iny)
  for(size_t inz = 0; inz < indz.size(); ++inz)
  {
    cellData.cellIndices.push_back(indx[inx] + nx*(indy[iny] + ny*indz[inz]));
    cellData.cellWeights.push_back(fx[inx]*fy[iny]*fz[inz]);
  }

  return cellData;
}
