/*
** Copyright 2007-2011 Erik Santiso.
** This file is part of mymol.
** mymol is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** mymol is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with mymol. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** PairPotential class
*/

#include <cmath>
#include "mymath/include/vector3D.h"
#include "mymath/include/mpconstants.h"
#include "mymol/include/pairpotential.h"

Real const PairPotential::energy(Real const &r) const
{
  if(r == 0)
  {
    std::cerr << "Error in PairPotential: Overlapping atoms" << std::endl;
    return 0.0;
  }
  Real energy = 0.0;

  switch(type)
  {
    case LENNARD_JONES:
      {
        Real const sigma_r = parameters[1]/r; // sigma/r
        Real const sigma_r2 = sigma_r*sigma_r; // (sigma/r)^2
        Real const sigma_r6 = sigma_r2*sigma_r2*sigma_r2; // (sigma/r)^6
        energy = prefactor*sigma_r6*(sigma_r6 - 1.0);
      }
      break;
    case MIE:
      {
        Real const sigma_r = parameters[1]/r; // sigma/r
        energy = prefactor*(pow(sigma_r, parameters[2]) -
                            pow(sigma_r, parameters[3]));
      }
      break;
    case BUCKINGHAM:
      {
        Real const one_r = 1.0/r;
        Real const one_r2 = one_r*one_r;
        Real const one_r6 = one_r2*one_r2*one_r2;
        energy = parameters[0]*exp(-parameters[1]*r) - parameters[2]*one_r6;
      }
      break;
    case BORN_HUGGINS_MEYER:
      {
        Real const one_r = 1.0/r;
        Real const one_r2 = one_r*one_r;
        Real const one_r4 = one_r2*one_r2;
        Real const one_r6 = one_r2*one_r4;
        Real const one_r8 = one_r4*one_r4;
        energy = parameters[0]*exp(parameters[1]*(parameters[4] - r)) -
                 parameters[2]*one_r6 - parameters[3]*one_r8;
      }
      break;
    case MORSE:
      {
        Real const expfactor = exp(-parameters[1]*(r - parameters[2]));
        energy = prefactor*expfactor*(expfactor - 2.0);
      }
      break;
    case WCA:
      {
        Real const sigma_r = parameters[1]/r; // sigma/r
        if(1.0/sigma_r < TWO_TO_ONE_SIXTH)
        {
          Real const sigma_r2 = sigma_r*sigma_r; // (sigma/r)^2
          Real const sigma_r6 = sigma_r2*sigma_r2*sigma_r2; // (sigma/r)^6
          energy = prefactor*sigma_r6*(sigma_r6 - 1.0);
        }
        else energy = 0.0;
      }
      break;
    case UNKNOWN:
    default:
      std::cerr << "Error in PairPotential::energy: "
                << "unknown pair potential." << std::endl;
  }

  return energy;
}

Real const PairPotential::energy2(Real const &r2) const
{
  if(r2 == 0)
  {
    std::cerr << "Error in PairPotential: Overlapping atoms" << std::endl;
    return 0.0;
  }
  Real energy = 0.0;

  switch(type)
  {
    case LENNARD_JONES:
      {
        Real const sigma_r2 = parameters[1]*parameters[1]/r2; // (sigma/r)^2
        Real const sigma_r6 = sigma_r2*sigma_r2*sigma_r2; // (sigma/r)^6
        energy = prefactor*sigma_r6*(sigma_r6 - 1.0);
      }
      break;
    case MIE:
      {
        Real const sigma_r2 = parameters[1]*parameters[1]/r2; // (sigma/r)^2
        energy = prefactor*(pow(sigma_r2, 0.5*parameters[2]) -
                            pow(sigma_r2, 0.5*parameters[3]));
      }
      break;
    case BUCKINGHAM:
      {
        Real const r = sqrt(r2);
        Real const one_r2 = 1.0/r2;
        Real const one_r6 = one_r2*one_r2*one_r2;
        energy = parameters[0]*exp(-parameters[1]*r) - parameters[2]*one_r6;
      }
      break;
    case BORN_HUGGINS_MEYER:
      {
        Real const r = sqrt(r2);
        Real const one_r2 = 1.0/r2;
        Real const one_r4 = one_r2*one_r2;
        Real const one_r6 = one_r2*one_r4;
        Real const one_r8 = one_r4*one_r4;
        energy = parameters[0]*exp(parameters[1]*(parameters[4] - r)) -
                 parameters[2]*one_r6 - parameters[3]*one_r8;
      }
      break;
    case MORSE:
      {
        Real const r = sqrt(r2);
        Real const expfactor = exp(-parameters[1]*(r - parameters[2]));
        energy = prefactor*expfactor*(expfactor - 2.0);
      }
      break;
    case WCA:
      {
        Real const sigma_r2 = parameters[1]*parameters[1]/r2; // (sigma/r)^2
        if(1.0/sigma_r2 < TWO_TO_ONE_THIRD)
        {
          Real const sigma_r6 = sigma_r2*sigma_r2*sigma_r2; // (sigma/r)^6
          energy = prefactor*sigma_r6*(sigma_r6 - 1.0);
        }
        else energy = 0.0;
      }
      break;
    case UNKNOWN:
    default:
      std::cerr << "Error in PairPotential::energy: "
                << "unknown pair potential." << std::endl;
  }

  return energy;
}

Vector3D const PairPotential::force(Vector3D const &r12) const
{
  Real const r2 = norm2(r12);
  if(r2 == 0)
  {
    std::cerr << "Error in PairPotential: Overlapping atoms" << std::endl;
    return Vector3D(0.0);
  }
  Vector3D force = 0.0;

  switch(type)
  {
    case LENNARD_JONES:
      {
        Real const sigma_r2 = parameters[1]*parameters[1]/r2; // (sigma/r)^2
        Real const sigma_r6 = sigma_r2*sigma_r2*sigma_r2; // (sigma/r)^6
        force = 6.0*prefactor*sigma_r6*(2.0*sigma_r6 - 1.0)*r12/r2;
      }
      break;
    case MIE:
      {
        Real const sigma_r = parameters[1]/sqrt(r2); // sigma/r
        force = prefactor*(parameters[2]*pow(sigma_r, parameters[2]) -
                           parameters[3]*pow(sigma_r, parameters[3]))*r12/r2;
      }
      break;
    case BUCKINGHAM:
      {
        Real const r = sqrt(r2);
        Real const one_r2 = 1.0/r2;
        Real const one_r6 = one_r2*one_r2*one_r2;
        force = (parameters[0]*parameters[1]*r*exp(-parameters[1]*r) - 
                 6.0*parameters[2]*one_r6)*r12/r2;
      }
      break;
    case BORN_HUGGINS_MEYER:
      {
        Real const r = sqrt(r2);
        Real const one_r2 = 1.0/r2;
        Real const one_r4 = one_r2*one_r2;
        Real const one_r6 = one_r2*one_r4;
        Real const one_r8 = one_r4*one_r4;
        force = (parameters[0]*parameters[1]*r*
                    exp(parameters[1]*(parameters[4] - r)) -
                  6.0*parameters[2]*one_r6 - 
                  8.0*parameters[3]*one_r8)*r12/r2;
      }
      break;
    case MORSE:
      {
        Real const r = sqrt(r2);
        Real const expfactor = exp(-parameters[1]*(r - parameters[2]));
        force = 2.0*parameters[1]*prefactor*expfactor*(expfactor - 1.0)*r12/r;
      }
      break;
    case WCA:
      {
        Real const sigma_r2 = parameters[1]*parameters[1]/r2; // (sigma/r)^2
        if(1.0/sigma_r2 < TWO_TO_ONE_THIRD)
        {
          Real const sigma_r6 = sigma_r2*sigma_r2*sigma_r2; // (sigma/r)^6
          force = 6.0*prefactor*sigma_r6*(2.0*sigma_r6 - 1.0)*r12/r2;
        }
        else force = 0.0;
      }
      break;
    case UNKNOWN:
    default:
      std::cerr << "Error in PairPotential::force: "
                << "unknown pair potential." << std::endl;
  }

  return force;
}

Real const PairPotential::tailCorrection(Real const &r_c) const
{
  if(r_c == 0)
  {
    std::cerr << "Error in PairPotential: Zero cutoff" << std::endl;
    return 0.0;
  }
  Real integral = 0.0;

  switch(type)
  {
    case LENNARD_JONES:
      {
        Real const sigma_r = parameters[1]/r_c; // sigma/r_c
        Real const sigma_r2 = sigma_r*sigma_r; // (sigma/r_c)^2
        Real const sigma_r6 = sigma_r2*sigma_r2*sigma_r2; // (sigma/r_c)^6
        Real const r_c3 = r_c*r_c*r_c;
        integral = prefactor*r_c3*sigma_r6*(sigma_r6 - 3.0)/9.0;
      }
      break;
    case MIE:
      {
        if(parameters[2] <= 3.0 || parameters[3] <= 3.0)
        {
          std::cerr << "Error in PairPotential: invalid Mie exponent"
                    << std::endl;
          return 0.0;
        }
        Real const sigma_r = parameters[1]/r_c; // sigma/r_c
        Real const r_c3 = r_c*r_c*r_c;
/* This is ignoring the repulsive contribution, as DL_POLY does "sometimes" :/
        integral = prefactor*r_c3*
                        (-pow(sigma_r, parameters[3])/(parameters[3] - 3.0));
*/
        integral = prefactor*r_c3*
                        (pow(sigma_r, parameters[2])/(parameters[2] - 3.0) -
                         pow(sigma_r, parameters[3])/(parameters[3] - 3.0));
      }
      break;
    case BUCKINGHAM:
      {
        Real const one_r = 1.0/r_c;
        Real const one_r3 = one_r*one_r*one_r;
        Real const b_r = parameters[1]*r_c;
        Real const b3 = parameters[1]*parameters[1]*parameters[1];
        integral = parameters[0]/b3*exp(-b_r)*(2.0 + b_r*(2.0 + b_r)) - 
                   parameters[2]*one_r3/3.0;
      }
      break;
    case BORN_HUGGINS_MEYER:
      {
        Real const one_r = 1.0/r_c;
        Real const one_r2 = one_r*one_r;
        Real const one_r3 = one_r2*one_r;
        Real const one_r7 = one_r2*one_r2*one_r3;
        Real const b_r = parameters[1]*r_c;
        Real const b3 = parameters[1]*parameters[1]*parameters[1];
        integral = parameters[0]/b3*exp(parameters[1]*(parameters[4] - r_c))*
                   (2.0 + b_r*(2.0 + b_r)) -
                   parameters[2]*one_r3/3.0 - 
                   parameters[3]*one_r7/7.0;
      }
      break;
    case MORSE:
      {
        Real const expfactor = exp(-parameters[1]*(r_c - parameters[2]));
        Real const a3 = parameters[1]*parameters[1]*parameters[1];
        Real const a_r = parameters[1]*r_c;
        integral = prefactor/(4.0*a3)*expfactor*
          (expfactor*(1.0 + a_r*(2.0 + 2.0*a_r) - 2.0) -
           8.0*(2.0 + a_r*(2.0 + a_r)));
      }
      break;
    case WCA:
      {
        Real const sigma_r = parameters[1]/r_c; // sigma/r
        if(1.0/sigma_r < TWO_TO_ONE_SIXTH)
        {
          std::cerr << "Error in PairPotential:tailCorrection: "
                    << "cutoff is too small" << std::endl;
          return 0.0;
        }
        else integral = 0.0;
      }
      break;
    case UNKNOWN:
    default:
      std::cerr << "Error in PairPotential::tailCorrection: "
                << "unknown pair potential." << std::endl;
  }
  return TWO_PI*integral;
}

// Private methods

void PairPotential::calculatePrefactor()
{
  switch(type)
  {
    case LENNARD_JONES:
      prefactor = 4.0*parameters[0];
      break;
    case MIE:
      {
        Real const &n = parameters[2];
        Real const &m = parameters[3];
        prefactor = parameters[0]*(n/(n - m))*
                    pow(n/m, m/(n - m));
      }
      break;
    case BUCKINGHAM:
      prefactor = 1.0; // Not used
      break;
    case BORN_HUGGINS_MEYER:
      prefactor = 1.0; // Not used
      break;
    case MORSE:
      prefactor = parameters[0];
      break;
    case WCA:
      prefactor = 4.0*parameters[0];
      break;
    case UNKNOWN:
    default:
      prefactor = 1.0;
      std::cerr << "Error in PairPotential::prefactor: "
                << "unknown pair potential." << std::endl;
  }
}
