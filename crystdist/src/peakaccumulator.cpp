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
** Peak Accumulator class
*/

#include "mymath/include/mpconstants.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/vector4D.h"
#include "mymath/include/matrix4D.h"
#include "mymath/include/quaternion.h"
#include "mymath/include/rootfinding.h"
#include "mymath/include/specialfunctions.h"
#include "mymath/include/symmetriceigensystem.h"
#include "crystdist/include/peakaccumulator.h"

// Constructors

PeakAccumulator::PeakAccumulator()
:
knowParameters_(false), numValues_(0), name_(), types_(), sumDistances_(0.0), 
sumDistancesSquared_(0.0), sumBondOrientations_(0.0), sumRelativeOrientations_(0.0), 
sumInternalDOFData_(), mean_(), concentrationParameters_(), normalizationFactors_(),
lengthOfMeanVector_(0.0)
{}

// Interface

void PeakAccumulator::addValue(RelativeConfiguration const &conf)
{
  if(numValues_ == 0)  // Store name, types and symmetry information
  {
    name_ = conf.name;
    types_ = conf.types;
    for(size_t i = 0; i < conf.internalDOFs.first.size(); i++)
      sumInternalDOFData_.first.push_back(InternalDOFAccumulator(conf.internalDOFs.first[i].type,
                                                                 conf.internalDOFs.first[i].symmetryNumber));
    for(size_t i = 0; i < conf.internalDOFs.second.size(); i++)
      sumInternalDOFData_.second.push_back(InternalDOFAccumulator(conf.internalDOFs.second[i].type,
                                                                  conf.internalDOFs.second[i].symmetryNumber));
  }
  else if(conf.name != name_) // Check consistency
  {
    std::cerr << "Error in PeakAccumulator::addValue: Data added is incompatible with data "
              << "already accumulated" << std::endl;
    return;
  }
  // Distance
  sumDistances_ += conf.distance;
  sumDistancesSquared_ += conf.distance*conf.distance;

  // Bond orientation and relative orientation
  if(conf.types.first == GENERAL && conf.types.second == GENERAL)
  {
    // General case: Bond orientation = Vector3D, Relative orientation = Quaternion
    sumBondOrientations_ += conf.bondOrientation;
    sumRelativeOrientations_ += outer(conf.relativeOrientation, conf.relativeOrientation);
  }
  else if(conf.types.first == LINEAR_ASYMMETRIC && conf.types.second == LINEAR_ASYMMETRIC)
  {
    // Axisymmetric case: Bond orientation and relative orientation are both angles
    Real bondAngle = conf.bondOrientation.w;
    sumBondOrientations_.w += cos(bondAngle);  // Not very elegant: w and x
    sumBondOrientations_.x += sin(bondAngle);  // store sine and cosine
    Real relativeAngle = conf.relativeOrientation.w;
    sumRelativeOrientations_(0, 0) += cos(relativeAngle); // Ditto: Here they are stored
    sumRelativeOrientations_(1, 1) += sin(relativeAngle); // in the (0, 0) and (1, 1) elements

  }
  else if((conf.types.first == LINEAR_SYMMETRIC && conf.types.second == LINEAR_SYMMETRIC) ||
          (conf.types.first == PLANAR_SYMMETRIC && conf.types.second == PLANAR_SYMMETRIC))
  {
    // Axisymmetric case with symmetry plane: Bond orientation and relative orientation are both angles
    Real bondAngle = 2.0*conf.bondOrientation.w;  // Due to the axial symmetry
    sumBondOrientations_.w += cos(bondAngle);  // Not very elegant: w and x
    sumBondOrientations_.x += sin(bondAngle);  // store sine and cosine
    Real relativeAngle = 2.0*conf.relativeOrientation.w;  // Due to the axial symmetry
    sumRelativeOrientations_(0, 0) += cos(relativeAngle); // Ditto: Here they are stored
    sumRelativeOrientations_(1, 1) += sin(relativeAngle); // in the (0, 0) and (1, 1) elements
  }
  else
  {
    std::cerr << "Error in PeakAccumulator::addValue(): Combination of molecule types not implemented" << std::endl;
    return;
  }

  // Internal degrees of freedom
  for(size_t i = 0; i < conf.internalDOFs.first.size(); i++)
  {
    switch(conf.internalDOFs.first[i].type)
    {
      case DISTANCE:
      {
        Real distance = conf.internalDOFs.first[i].value;
        sumInternalDOFData_.first[i].sumC += distance;
        sumInternalDOFData_.first[i].sumS += distance*distance;
      }
        break;
      case ANGLE:
      {
        Real angle = conf.internalDOFs.first[i].value;
        sumInternalDOFData_.first[i].sumC += cos(angle);
        sumInternalDOFData_.first[i].sumS += sin(angle);
      }
        break;
      case DIHEDRAL:
      {
        Real angle = conf.internalDOFs.first[i].value*
                     (Real)conf.internalDOFs.first[i].symmetryNumber;
        sumInternalDOFData_.first[i].sumC += cos(angle);
        sumInternalDOFData_.first[i].sumS += sin(angle);
      }
        break;
    }
  }

  for(size_t i = 0; i < conf.internalDOFs.second.size(); i++)
  {
    switch(conf.internalDOFs.second[i].type)
    {
      case DISTANCE:
      {
        Real distance = conf.internalDOFs.second[i].value;
        sumInternalDOFData_.second[i].sumC += distance;
        sumInternalDOFData_.second[i].sumS += distance*distance;
      }
        break;
      case ANGLE:
      {
        Real angle = conf.internalDOFs.second[i].value;
        sumInternalDOFData_.second[i].sumC += cos(angle);
        sumInternalDOFData_.second[i].sumS += sin(angle);
      }
        break;
      case DIHEDRAL:
      {
        Real angle = conf.internalDOFs.second[i].value*
                     (Real)conf.internalDOFs.second[i].symmetryNumber;
        sumInternalDOFData_.second[i].sumC += cos(angle);
        sumInternalDOFData_.second[i].sumS += sin(angle);
      }
        break;
    }
  }
  ++numValues_;
  knowParameters_ = false;
}

void PeakAccumulator::clear()
{
  knowParameters_ = false;
  numValues_ = 0;
  name_ = "";
  sumDistances_ = 0.0;
  sumDistancesSquared_ = 0.0;
  sumBondOrientations_ = 0.0;
  sumRelativeOrientations_ = 0.0;
  sumInternalDOFData_.first.clear();
  sumInternalDOFData_.second.clear();
  lengthOfMeanVector_ = 0.0;
}

// Target functions for maximum likelihood estimation

Real const PeakAccumulator::vonMisesTarget(Real const &concentrationParameter)
{
  if(concentrationParameter <= 0.0) return -lengthOfMeanVector_;
  else return bessel_I1_over_I0(concentrationParameter) - lengthOfMeanVector_;
}

Real const PeakAccumulator::vonMisesFisherTarget(Real const &concentrationParameter)
{
  if(concentrationParameter <= 0.0) return -lengthOfMeanVector_;
  else 
  {
    Real const expmtwoXi = exp(-concentrationParameter - concentrationParameter);
    return ((1.0 + expmtwoXi)/(1.0 - expmtwoXi) - 1.0/concentrationParameter
            - lengthOfMeanVector_);
  }
}

Real const PeakAccumulator::bipolarWatsonTarget(Real const &concentrationParameter)
{
  if(concentrationParameter <= 0.0) return -lengthOfMeanVector_;
  else return d_ln_confluent(0.5, 2.0, concentrationParameter) -
              lengthOfMeanVector_;
}

// Print current statistics (mostly for debugging)

std::ostream &operator<<(std::ostream &outStream, PeakAccumulator &accumulator)
{
  outStream << std::endl << "Number of data points: " << accumulator.numValues();
  outStream << std::endl << "Mean relative configuration: " << accumulator.mean();
  outStream << std::endl << "Concentration parameters: " << accumulator.concentrationParameters();
  outStream << std::endl << "Normalization factors: " << accumulator.normalizationFactors();
  return outStream;
}

// Private functions

void PeakAccumulator::calculateParameters()
{
  // Check that there's actually data in the accumulator
  if(numValues_ < 1)
  {
    std::cerr << "Error in PeakAccumulator::calculateParameters(): Not enough data added" << std::endl;
    return;
  }

  // Store types and name into mean
  mean_.types = types_;
  mean_.name = name_;

  // Distance
  Real const meanDistance = sumDistances_/(Real)numValues_;
  Real const meanDistanceSq = sumDistancesSquared_/(Real)numValues_;
  Real const distanceVariance = meanDistanceSq - meanDistance*meanDistance;
  /* // Removed for compatibility with PeakClustering 
  if(distanceVariance == 0.0)
  {
    std::cerr << "Warning from PeakAccumulator::calculateParameters(): Zero distance variance" << std::endl;
    return;
  }
  */
  mean_.distance = meanDistance;
  if(distanceVariance != 0.0)
  {
    concentrationParameters_.distance = 1.0/distanceVariance;
    normalizationFactors_.distance = 1.0/sqrt(TWO_PI*distanceVariance);
  }
  else
  {
    concentrationParameters_.distance = REAL_VERYBIG;
    normalizationFactors_.distance = REAL_VERYBIG;
  }

  // Bond orientation and relative orientation
  if(types_.first == GENERAL && types_.second == GENERAL)
  {
    // General case: Bond orientation = Vector3D, Relative orientation = Quaternion
    Vector3D const meanBondOrientation = vector(sumBondOrientations_)/(Real)numValues_;
    lengthOfMeanVector_ = norm(meanBondOrientation);
    mean_.bondOrientation = unit(meanBondOrientation);
    // Likelihood maximization for von Mises-Fisher distribution
    // see A. Tanabe et al., Comput. Stat. 22, 145 (2007)
    if(lengthOfMeanVector_ == 0.0) 
    {
      concentrationParameters_.bondOrientation = 0.0;
      normalizationFactors_.bondOrientation = 1.0/(4.0*PI);
    }
    else if((lengthOfMeanVector_ == 1.0) ||
            (numValues_ == 1)) 
    {
      concentrationParameters_.bondOrientation = REAL_VERYBIG;
      normalizationFactors_.bondOrientation = 0.0;
    }
    else 
    {
      Real const lowerBound = lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
      Real const upperBound = 3.0*lowerBound;
      concentrationParameters_.bondOrientation = 
        bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesFisherTarget, 
                                   lowerBound, upperBound, 1.E-6);
      normalizationFactors_.bondOrientation = concentrationParameters_.bondOrientation/
                                              (4.0*PI*sinh(concentrationParameters_.bondOrientation));
    }
    // Likelihood maximization for bipolar Watson distribution
    // see A. Figueiredo and P. Gomes, Stat. Prob. Lett. 76, 142 (2006)
    // also Mardia and Jupp, "Directional Statistics", Wiley (2000)
    SymmetricEigensystem<Matrix4D, Vector4D> eig(sumRelativeOrientations_);
    Vector4D const &eigenvalues = eig.eigenvalues();
    Matrix4D const &eigenvectors = eig.eigenvectors();
    Real maxEigenvalue = eigenvalues[0]; size_t maxEigenvalueIndex = 0;
    for(size_t i = 1; i < 4; i++)
      if(eigenvalues[i] > maxEigenvalue)
      {
        maxEigenvalue = eigenvalues[i];
        maxEigenvalueIndex = i;
      }
    mean_.relativeOrientation = eigenvectors.column(maxEigenvalueIndex);
    lengthOfMeanVector_ = maxEigenvalue/(Real)numValues_;
    if(lengthOfMeanVector_ < 0.25)
    {
      std::cerr << "Bad relative orientation data in PeakAccumulator::calculateParameters" << std::endl;
      concentrationParameters_.relativeOrientation = 0.0;
      normalizationFactors_.relativeOrientation = 0.0;
    }
    else if(numValues_ == 1)
    {
      concentrationParameters_.relativeOrientation = REAL_VERYBIG;
      normalizationFactors_.relativeOrientation = 0.0;
    }
    else  // Bounds are not rigorous, just found through quick numerical experimentation
    {
      Real const lowerBound = (lengthOfMeanVector_ < 0.8)?
                              (12*(lengthOfMeanVector_ - 0.25)):
                              1.0/(1.0-lengthOfMeanVector_);
      Real const upperBound = (lengthOfMeanVector_ < 0.8)?
                              (16*(lengthOfMeanVector_ - 0.25)):
                              2.0/(1.0-lengthOfMeanVector_);
      concentrationParameters_.relativeOrientation =
        bisection<PeakAccumulator>(this, &PeakAccumulator::bipolarWatsonTarget,
                                   lowerBound, upperBound, 1.E-6);
      normalizationFactors_.relativeOrientation = 
        1.0/confluent(0.5, 2.0, concentrationParameters_.relativeOrientation);
    }
  }
  else if(types_.first == LINEAR_ASYMMETRIC && types_.second == LINEAR_ASYMMETRIC)
  {
    // Axisymmetric case: Bond orientation and relative orientation are both angles
    Real meanCos = sumBondOrientations_.w/(Real)numValues_; // As above, cos and sin
    Real meanSin = sumBondOrientations_.x/(Real)numValues_; // stored in w and x
    mean_.bondOrientation = atan2(meanSin, meanCos);
    // Likelihood maximization for von Mises distribution
    // See Singh et al., Comm. Stat. Simulat. Comput. 34, 21 (2005)
    // also Mardia and Jupp, "Directional Statistics", Wiley (2000)
    lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
    if(lengthOfMeanVector_ == 0.0) 
    {
      concentrationParameters_.bondOrientation = 0.0;
      normalizationFactors_.bondOrientation = 1/TWO_PI;
    }
    else if((lengthOfMeanVector_ == 1.0) ||
            (numValues_ == 1))
    {
      concentrationParameters_.bondOrientation = REAL_VERYBIG;
      normalizationFactors_.bondOrientation = 0.0;
    }
    else 
    {
      Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
      concentrationParameters_.bondOrientation = 
        bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
      normalizationFactors_.bondOrientation = 1.0/(TWO_PI*bessel_I0(concentrationParameters_.bondOrientation));
    } 
    // Relative orientation
    meanCos = sumRelativeOrientations_(0, 0)/(Real)numValues_;  // As above, cos and sin are
    meanSin = sumRelativeOrientations_(1, 1)/(Real)numValues_;  // stored in (0, 0) and (1, 1)
    mean_.relativeOrientation = atan2(meanSin, meanCos);
    lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
    if(lengthOfMeanVector_ == 0.0) 
    {
      concentrationParameters_.relativeOrientation = 0.0;
      normalizationFactors_.relativeOrientation = 1.0/TWO_PI;
    }
    else if((lengthOfMeanVector_ == 1.0) ||
            (numValues_ == 1))
    {
      concentrationParameters_.relativeOrientation = REAL_VERYBIG;
      normalizationFactors_.relativeOrientation = 0.0;
    }
    else
    {
      Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
      concentrationParameters_.relativeOrientation = 
        bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
      normalizationFactors_.relativeOrientation = 
        1.0/(TWO_PI*bessel_I0(concentrationParameters_.relativeOrientation));
    }
  }
  else if((types_.first == LINEAR_SYMMETRIC && types_.second == LINEAR_SYMMETRIC) ||
          (types_.first == PLANAR_SYMMETRIC && types_.second == PLANAR_SYMMETRIC))
  {
    // Axisymmetric case with symmetry plane: Bond orientation and relative orientation are both angles
    Real meanCos = sumBondOrientations_.w/(Real)numValues_; // As above, cos and sin
    Real meanSin = sumBondOrientations_.x/(Real)numValues_; // stored in w and x
    mean_.bondOrientation = atan2(meanSin, meanCos)/2.0;    // Factor of 2 from symmetry
    // Likelihood maximization for von Mises distribution
    // (same as axisymmetric case, but with double angle)
    lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
    if(lengthOfMeanVector_ == 0.0) 
    {
      concentrationParameters_.bondOrientation = 0.0;
      normalizationFactors_.bondOrientation = 1/TWO_PI;
    }
    else if((lengthOfMeanVector_ == 1.0) ||
            (numValues_ == 1))
    {
      concentrationParameters_.bondOrientation = REAL_VERYBIG;
      normalizationFactors_.bondOrientation = 0.0;
    }
    else 
    {
      Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
      concentrationParameters_.bondOrientation = 
        bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
      normalizationFactors_.bondOrientation = 1.0/(TWO_PI*bessel_I0(concentrationParameters_.bondOrientation));
    }
    
    // Relative orientation
    meanCos = sumRelativeOrientations_(0, 0)/(Real)numValues_;  // As above, cos and sin are
    meanSin = sumRelativeOrientations_(1, 1)/(Real)numValues_;  // stored in (0, 0) and (1, 1)
    mean_.relativeOrientation = atan2(meanSin, meanCos)/2.0;    // Factor of 2 from symmetry
    lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
    if(lengthOfMeanVector_ == 0.0) 
    {
      concentrationParameters_.relativeOrientation = 0.0;
      normalizationFactors_.relativeOrientation = 1.0/TWO_PI;
    }
    else if((lengthOfMeanVector_ == 1.0) ||
            (numValues_ == 1))
    {
      concentrationParameters_.relativeOrientation = REAL_VERYBIG;
      normalizationFactors_.relativeOrientation = 0.0;
    }
    else
    {
      Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
      concentrationParameters_.relativeOrientation = 
        bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
      normalizationFactors_.relativeOrientation = 
        1.0/(TWO_PI*bessel_I0(concentrationParameters_.relativeOrientation));
    }
  }
  else  // Should never get here
  {
    std::cerr << "Error in PeakAccumulator::calculateParameters: Combination of molecule types not implemented" << std::endl;
    return;
  }

  // Internal degrees of freedom
  mean_.internalDOFs.first.clear();
  concentrationParameters_.internalDOFs.first.clear();
  normalizationFactors_.internalDOFs.first.clear();
  for(size_t i = 0; i < sumInternalDOFData_.first.size(); i++)
  {
    switch(sumInternalDOFData_.first[i].type)
    {
      case DISTANCE:
      {
        // Gaussian
        Real const meanDistance = sumInternalDOFData_.first[i].sumC/(Real)numValues_;
        Real const meanDistanceSq = sumInternalDOFData_.first[i].sumS/(Real)numValues_;
        Real const distanceVariance = meanDistanceSq - meanDistance*meanDistance;
        /* // Removed for compatibility with PeakClustering
        if(distanceVariance == 0.0)
        {
          std::cerr << "Warning from PeakAccumulator::calculateParameters(): "
                    << "Zero internal distance variance" << std::endl;
          return;
        }
        */
        mean_.internalDOFs.first.push_back(InternalDOF(DISTANCE ,meanDistance, 1));
        if(distanceVariance != 0.0)
        {
          concentrationParameters_.internalDOFs.first.push_back(1.0/distanceVariance);
          normalizationFactors_.internalDOFs.first.push_back(1.0/sqrt(TWO_PI*distanceVariance));
        }
        else
        {
          concentrationParameters_.internalDOFs.first.push_back(REAL_VERYBIG);
          normalizationFactors_.internalDOFs.first.push_back(REAL_VERYBIG);
        }
      }
        break;
      case ANGLE:
      {
        // von Mises
        Real const meanCos = sumInternalDOFData_.first[i].sumC/(Real)numValues_;
        Real const meanSin = sumInternalDOFData_.first[i].sumS/(Real)numValues_;
        mean_.internalDOFs.first.push_back(InternalDOF(ANGLE, atan2(meanSin, meanCos), 1));
        Real concentrationParameter, normalizationFactor;
        lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
        if(lengthOfMeanVector_ == 0.0) 
        {
          concentrationParameter = 0.0;
          normalizationFactor = 1.0/TWO_PI;
        }
        else if((lengthOfMeanVector_ == 1.0) ||
                (numValues_ == 1))
        {
          concentrationParameter = REAL_VERYBIG;
          normalizationFactor = 0.0;
        }
        else
        {
          Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
          concentrationParameter = 
            bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
          normalizationFactor = 1.0/(TWO_PI*bessel_I0(concentrationParameter));
        }
        
        concentrationParameters_.internalDOFs.first.push_back(concentrationParameter);
        normalizationFactors_.internalDOFs.first.push_back(normalizationFactor);
      }
        break;
      case DIHEDRAL:
      {
        // von Mises with symmetry factor
        Real const meanCos = sumInternalDOFData_.first[i].sumC/(Real)numValues_;
        Real const meanSin = sumInternalDOFData_.first[i].sumS/(Real)numValues_;
        Real const meanAngle = atan2(meanSin, meanCos)/(Real)sumInternalDOFData_.first[i].symmetryNumber;
        mean_.internalDOFs.first.push_back(
          InternalDOF(DIHEDRAL, meanAngle, sumInternalDOFData_.first[i].symmetryNumber));
        Real concentrationParameter, normalizationFactor;
        lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
        if(lengthOfMeanVector_ == 0.0) 
        {
          concentrationParameter = 0.0;
          normalizationFactor = 1.0/TWO_PI;
        }
        else if((lengthOfMeanVector_ == 1.0) ||
                (numValues_ == 1))
        {
          concentrationParameter = REAL_VERYBIG;
          normalizationFactor = 0.0;
        }
        else
        {
          Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
          concentrationParameter = 
            bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
          normalizationFactor = 1.0/(TWO_PI*bessel_I0(concentrationParameter));
        }
        concentrationParameters_.internalDOFs.first.push_back(concentrationParameter);
        normalizationFactors_.internalDOFs.first.push_back(normalizationFactor);
      }
        break;
    }
  }

  for(size_t i = 0; i < sumInternalDOFData_.second.size(); i++)
  {
    switch(sumInternalDOFData_.second[i].type)
    {
      case DISTANCE:
      {
        // Gaussian
        Real const meanDistance = sumInternalDOFData_.second[i].sumC/(Real)numValues_;
        Real const meanDistanceSq = sumInternalDOFData_.second[i].sumS/(Real)numValues_;
        Real const distanceVariance = meanDistanceSq - meanDistance*meanDistance;
        if(distanceVariance == 0.0)
        {
          std::cerr << "Warning from PeakAccumulator::calculateParameters(): "
                    << "Zero internal distance variance" << std::endl;
          return;
        }
        mean_.internalDOFs.second.push_back(InternalDOF(DISTANCE ,meanDistance, 1));
        if(distanceVariance != 0.0)
        {
          concentrationParameters_.internalDOFs.second.push_back(1.0/distanceVariance);
          normalizationFactors_.internalDOFs.second.push_back(1.0/sqrt(TWO_PI*distanceVariance));
        }
        else
        {
          concentrationParameters_.internalDOFs.second.push_back(REAL_VERYBIG);
          normalizationFactors_.internalDOFs.second.push_back(REAL_VERYBIG);
        }
      }
        break;
      case ANGLE:
      {
        // von Mises
        Real const meanCos = sumInternalDOFData_.second[i].sumC/(Real)numValues_;
        Real const meanSin = sumInternalDOFData_.second[i].sumS/(Real)numValues_;
        mean_.internalDOFs.second.push_back(InternalDOF(ANGLE, atan2(meanSin, meanCos), 1));
        Real concentrationParameter, normalizationFactor;
        lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
        if(lengthOfMeanVector_ == 0.0) 
        {
          concentrationParameter = 0.0;
          normalizationFactor = 1.0/TWO_PI;
        }
        else if((lengthOfMeanVector_ == 1.0) ||
                (numValues_ == 1))
        {
          concentrationParameter = REAL_VERYBIG;
          normalizationFactor = 0.0;
        }
        else
        {
          Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
          concentrationParameter = 
            bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
          normalizationFactor = 1.0/(TWO_PI*bessel_I0(concentrationParameter));
        }
        concentrationParameters_.internalDOFs.second.push_back(concentrationParameter);
        normalizationFactors_.internalDOFs.second.push_back(normalizationFactor);
      }
        break;
      case DIHEDRAL:
      {
        // von Mises with symmetry factor
        Real const meanCos = sumInternalDOFData_.second[i].sumC/(Real)numValues_;
        Real const meanSin = sumInternalDOFData_.second[i].sumS/(Real)numValues_;
        Real const meanAngle = atan2(meanSin, meanCos)/(Real)sumInternalDOFData_.second[i].symmetryNumber;
        mean_.internalDOFs.second.push_back(
          InternalDOF(DIHEDRAL, meanAngle, sumInternalDOFData_.second[i].symmetryNumber));
        Real concentrationParameter, normalizationFactor;
        lengthOfMeanVector_ = sqrt(meanCos*meanCos + meanSin*meanSin);
        if(lengthOfMeanVector_ == 0.0) 
        {
          concentrationParameter = 0.0;
          normalizationFactor = 1.0/TWO_PI;
        }
        else if((lengthOfMeanVector_ == 1.0) ||
                (numValues_ == 1))
        {
          concentrationParameter = REAL_VERYBIG;
          normalizationFactor = 0.0;
        }
        else
        {
          Real const upperBound = 2.0*lengthOfMeanVector_/(1.0 - lengthOfMeanVector_*lengthOfMeanVector_);
          concentrationParameter = 
            bisection<PeakAccumulator>(this, &PeakAccumulator::vonMisesTarget, 0.0, upperBound, 1.E-6);
          normalizationFactor = 1.0/(TWO_PI*bessel_I0(concentrationParameter));
        }
        concentrationParameters_.internalDOFs.second.push_back(concentrationParameter);
        normalizationFactors_.internalDOFs.second.push_back(normalizationFactor);
      }
        break;
    }
  }
  knowParameters_ = true;
}

