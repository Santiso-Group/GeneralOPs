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
** Simple pair potential class
**
** Note that all data members are public. Checking the consistency of
** the data is the responsibility of the containing class/function.
**
** There is no explicit cutoff in the energy and force methods.
** Any cutoff needs to be handled externally - in-class handling may be a good
** feature for a later version.
**
** Note for a later version that it may be good to have a separate prefactor 
** for forces too.
*/

#ifndef H_PAIRPOTENTIAL
#define H_PAIRPOTENTIAL

#include <iostream>
#include <string>
#include <vector>
#include "mymath/include/vector3D.h"

// Types of pair potentials - Add as needed

enum PotentialType
{
  LENNARD_JONES,       // U(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
                       // Parameters: epsilon, sigma
  MIE,                 // U(r) = epsilon*[n/(n-m)]*(n/m)^[m/(n-m)]*
                       //        [(sigma/r)^n - (sigma/r)^m]
                       // Parameters: epsilon, sigma, n, m
  BUCKINGHAM,          // U(r) = A*exp(-B*r) - C/r^6
                       // Parameters: A, B, C
  BORN_HUGGINS_MEYER,  // U(r) = A*exp[B*(sigma-r)] - C/r^6 - D/r^8
                       // Parameters: A, B, C, D, sigma
  MORSE,               // U(r) = D*{exp[-2*a*(r-r0)] - 2*exp[-a*(r-r0)]}
                       // Parameters: D, a, r0
  WCA,                 // U(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6] +
                       // epsilon; for r < sigma*2^(1/6) 
                       // U(r) = 0 for r >= sigma*2^(1/6)
                       // Parameters: epsilon, sigma
  UNKNOWN              // Denotes an unknown potential function
};

// Output for pair potential types - mostly for debugging, add as needed

inline std::ostream &operator<<(std::ostream &outStream, 
                                PotentialType const &type)
{
  switch(type)
  {
    case LENNARD_JONES:
      outStream << "Lennard-Jones";
      break;
    case MIE:
      outStream << "Mie";
      break;
    case BUCKINGHAM:
      outStream << "Buckingham";
      break;
    case BORN_HUGGINS_MEYER:
      outStream << "Born-Huggins-Meyer";
      break;
    case MORSE:
      outStream << "Morse";
      break;
    case WCA:
      outStream << "Weeks-Chandler-Andersen (WCA)";
      break;
    case UNKNOWN:
    default:
      outStream << "Unrecognized potential";
  }
  return outStream;
}

// Pair potential class

class PairPotential
{

public:

  std::string       first;      // First atom type
  std::string       second;     // Second atom type
  PotentialType     type;       // Potential type
  std::vector<Real> parameters; // Container for potential parameters, if used
  Real              prefactor;  // Prefactor for potential energy

// Constructors

  PairPotential();                                   // Defines an empty 
                                                     // pair potential
  PairPotential(std::string const &first,            // Defines a pair potential
                std::string const &second,           // between two given atom
                PotentialType const &type = UNKNOWN, // types and, optionally,
                std::vector<Real> const              // gives the type and the
                &pars = std::vector<Real>());        // parameters

// Calculate potential energy and force

  Real const energy(Real const &r) const; // U(r)
  Real const energy2(Real const &r2) const; // U(r^2) -> could be faster
  Vector3D const force(Vector3D const &r12) const; // -grad_r2(U(r2 - r1))
  Real const tailCorrection(Real const &r_c) const; 
    // 2*pi*integral(r_c, infinity, r^2*U(r)*dr)
                                                    

// Print pair potential information (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, 
                                  PairPotential const &potential); 

// Calculate the prefactor from the parameters

  void calculatePrefactor();
};

/*
** End of class PairPotential
*/

// Inlines

inline PairPotential::PairPotential()
: 
first(""), second(""), type(UNKNOWN), parameters(), prefactor(1)
{}

inline PairPotential::PairPotential(std::string const &first, 
                                    std::string const &second,
                                    PotentialType const &type, 
                                    std::vector<Real> const &parameters)
:
first(first), second(second), type(type), parameters(parameters)
{ calculatePrefactor(); }

inline std::ostream &operator<<(std::ostream &outStream, 
                                PairPotential const &potential)
{
  outStream << potential.first << " " 
            << potential.second << " " 
            << potential.type;
  for(size_t i = 0; i < potential.parameters.size(); ++i)
    outStream << " " << potential.parameters[i];
  return outStream;
}

#endif

