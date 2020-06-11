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
** Simple lattice class
*/

#include <iostream>
#include "mymath/include/mpconstants.h"
#include "common/include/assert.h"
#include "mymol/include/lattice.h"

// Constructors

Lattice::Lattice(LatticeType const &latticeType)
:
latticeType_(latticeType)
{}

Lattice::Lattice(Vector3D const &a, LatticeType const &latticeType)
:
latticeType_(latticeType)
{
  assert(norm2(a));
  latticeVectors_.push_back(a);
  updateLattice();
}

Lattice::Lattice(Vector3D const &a1, Vector3D const &a2, LatticeType const &latticeType)
:
latticeType_(latticeType)
{
  assert(norm2(a1) + norm2(a2));
  latticeVectors_.push_back(a1);
  latticeVectors_.push_back(a2);
  updateLattice();
}

Lattice::Lattice(Vector3D const &a1, Vector3D const &a2, Vector3D const &a3, 
                 LatticeType const &latticeType)
:
latticeType_(latticeType)
{
  assert(norm2(a1) + norm2(a2) + norm2(a3));
  latticeVectors_.push_back(a1);
  latticeVectors_.push_back(a2);
  latticeVectors_.push_back(a3);
  updateLattice();
}

Lattice::Lattice(Real const &a, Real const &b, Real const &c,
                 Real const &alpha, Real const &beta, Real const &gamma, 
                 LatticeType const &latticeType)
:
latticeType_(latticeType)
{
  calculateLatticeVectors(a, b, c, alpha, beta, gamma);
  updateLattice();
}

// Interface 

void Lattice::setType(LatticeType const &latticeType)
{
  latticeType_ = latticeType;
  updateLattice();
}

void Lattice::setLattice(Vector3D const &a, LatticeType const &latticeType)
{
  latticeType_ = latticeType;
  assert(norm2(a));
  latticeVectors_.clear();
  latticeVectors_.push_back(a);
  updateLattice();
}

void Lattice::setLattice(Vector3D const &a1, Vector3D const &a2, LatticeType const &latticeType)
{
  latticeType_ = latticeType;
  assert(norm2(a1) + norm2(a2));
  latticeVectors_.clear();
  latticeVectors_.push_back(a1);
  latticeVectors_.push_back(a2);
  updateLattice();
}

void Lattice::setLattice(Vector3D const &a1, Vector3D const &a2, Vector3D const &a3, 
                         LatticeType const &latticeType)
{
  latticeType_ = latticeType;
  assert(norm2(a1) + norm2(a2) + norm2(a3));
  latticeVectors_.clear();
  latticeVectors_.push_back(a1);
  latticeVectors_.push_back(a2);
  latticeVectors_.push_back(a3);
  updateLattice();
}

void Lattice::setLattice(Real const &a, Real const &b, Real const &c,
                         Real const &alpha, Real const &beta, Real const &gamma,
                         LatticeType const &latticeType)
{
  latticeType_ = latticeType;
  calculateLatticeVectors(a, b, c, alpha, beta, gamma);
  updateLattice();
}

void Lattice::clear(LatticeType const &latticeType)
{
  latticeType_ = latticeType;
  latticeVectors_.clear();
  reciprocalVectors_.clear();
  imagePositions_.clear();
}

// Accessors

Vector3D const &Lattice::latticeVector(size_t const i) const
{
  assert(i < latticeVectors_.size());
  return latticeVectors_[i];
}

Vector3D const &Lattice::reciprocalVector(size_t const i) const
{
  assert(i < reciprocalVectors_.size());
  return reciprocalVectors_[i];
}

LatticeType const &Lattice::type() const
{
  return latticeType_;
}

// Other functions

Vector3D const Lattice::difference(Vector3D const &v1, Vector3D const &v2) const
{
  if(!isPeriodic()) return v2 - v1;
  Vector3D origDiff, diff;
  origDiff = diff = v2 - v1;
  for(size_t i = 0; i < latticeVectors_.size(); i++)
    diff -= latticeVectors_[i]*floor(reciprocalVectors_[i]*origDiff + 0.5);
  if(latticeType_ == CRYSTAL)
  {
    Real const numImages = imagePositions_.size();
    Real minDistance2 = norm2(diff);
    bool foundCloserImage = false;
    do
    {
      foundCloserImage = false;
      for(size_t i = 0; i < numImages; i++)
      {
        Vector3D newDiff = diff + imagePositions_[i];
        Real newDistance2 = norm2(newDiff);
        if(newDistance2 < minDistance2)
        {
          foundCloserImage = true;
          minDistance2 = newDistance2;
          diff = newDiff;
        }
      }
    }
    while(foundCloserImage);
  }
  return diff;
}

size_t const Lattice::numDimensions() const
{
  return latticeVectors_.size();
}

Real const Lattice::a() const
{
  return (latticeVectors_.size() > 0)?norm(latticeVectors_[0]):0.0;
}

Real const Lattice::b() const
{
  return (latticeVectors_.size() > 1)?norm(latticeVectors_[1]):0.0;
}

Real const Lattice::c() const
{
  return (latticeVectors_.size() > 2)?norm(latticeVectors_[2]):0.0;
}

Real const Lattice::alpha() const
{
  return (latticeVectors_.size() > 2)?
          RAD_TO_GRAD*acos((latticeVectors_[1]*latticeVectors_[2])/
                           (norm(latticeVectors_[1])*norm(latticeVectors_[2])))
          :0.0;
}

Real const Lattice::beta() const
{
  return (latticeVectors_.size() > 2)?
          RAD_TO_GRAD*acos((latticeVectors_[0]*latticeVectors_[2])/
                           (norm(latticeVectors_[0])*norm(latticeVectors_[2])))
          :0.0;
}
Real const Lattice::gamma() const
{
  return (latticeVectors_.size() > 1)?
          RAD_TO_GRAD*acos((latticeVectors_[0]*latticeVectors_[1])/
                           (norm(latticeVectors_[0])*norm(latticeVectors_[1])))
          :0.0;
}

Real const Lattice::length() const
{
  return (latticeVectors_.size() == 1)?norm(latticeVectors_[0]):0.0;
}

Real const Lattice::area() const
{
  return (latticeVectors_.size() == 2)?
          norm(cross(latticeVectors_[0], latticeVectors_[1]))
          :0.0;
}

Real const Lattice::volume() const
{
  return (latticeVectors_.size() == 3)?
          cross(latticeVectors_[0], latticeVectors_[1])*latticeVectors_[2]:
          0.0;
}

bool const Lattice::isPeriodic() const
{
  return (latticeVectors_.size() > 0)?true:false;
}

// Print lattice information (mostly for debugging)

std::ostream &operator<<(std::ostream &outStream, Lattice const &lattice)
{
  if(!lattice.isPeriodic()) return outStream << "System does not have periodic boundary conditions";
  outStream << "Lattice type: " << ((lattice.type() == SUPERCELL)?"SUPERCELL":"CRYSTAL") << std::endl;
  size_t nDim = lattice.numDimensions();
  outStream << "Number of dimensions: " << nDim << std::endl;
  outStream << "Lattice vectors: " << std::endl;
  for(size_t i = 0; i < nDim; i++)
    outStream << lattice.latticeVectors_[i] << std::endl;
  outStream << "Reciprocal lattice vectors: " << std::endl;
  for(size_t i = 0; i < nDim; i++)
    outStream << lattice.reciprocalVectors_[i] << std::endl;
  if(lattice.type() == CRYSTAL)
  {
    outStream << "Positions of central and neighboring cells: " << std::endl;
    for(size_t i = 0; i < lattice.imagePositions_.size(); i++)
      outStream << lattice.imagePositions_[i] << std::endl; 
  }
  outStream << "Crystallographic constants:" << std::endl;
  outStream << "a = " << lattice.a() << std::endl;
  outStream << "b = " << lattice.b() << std::endl;
  outStream << "c = " << lattice.c() << std::endl;
  outStream << "alpha = " << lattice.alpha() << std::endl;
  outStream << "beta = " << lattice.beta() << std::endl;
  outStream << "gamma = " << lattice.gamma() << std::endl;
  switch(nDim)
  {
  case 1:
    outStream << "Length of the unit cell: " << lattice.length();
    break;
  case 2:
    outStream << "Area of the unit cell: " << lattice.area();
    break;
  case 3:
    outStream << "Volume of the unit cell: " << lattice.volume();
    break;
  default:
    outStream << "Error: Wrong number of dimensions.";
  }
  return outStream;
}

// Private methods

void Lattice::updateLattice()
{
  reciprocalVectors_.clear();
  imagePositions_.clear();

  switch(latticeVectors_.size())
  {
  case 0:
    break;
  case 1:
    {
      Real den = norm2(latticeVectors_[0]);
      assert(den);
      reciprocalVectors_.push_back(latticeVectors_[0]/den);
      if(latticeType_ == CRYSTAL)
        for(int i = -1; i <= 1; i++)
          imagePositions_.push_back(i*latticeVectors_[0]);
    }
    break;
  case 2:
    {
      Vector3D normal = cross(latticeVectors_[0], latticeVectors_[1]);
      Vector3D reciprocal = cross(latticeVectors_[1], normal);
      Real den = latticeVectors_[0]*reciprocal;
      assert(den);
      reciprocalVectors_.push_back(reciprocal/den);
      reciprocal = cross(normal, latticeVectors_[0]);
      den = latticeVectors_[1]*reciprocal;
      assert(den);
      reciprocalVectors_.push_back(reciprocal/den);
      if(latticeType_ == CRYSTAL)
        for(int i = -1; i <= 1; i++)
          for(int j = -1; j <= 1; j++)
            imagePositions_.push_back(i*latticeVectors_[0] + j*latticeVectors_[1]);
    }
    break;
  case 3:
    {
      Vector3D reciprocal = cross(latticeVectors_[1], latticeVectors_[2]);
      Real den = latticeVectors_[0]*reciprocal;
      assert(den);
      reciprocalVectors_.push_back(reciprocal/den);
      reciprocal = cross(latticeVectors_[2], latticeVectors_[0]);
      den = latticeVectors_[1]*reciprocal;
      assert(den);
      reciprocalVectors_.push_back(reciprocal/den);
      reciprocal = cross(latticeVectors_[0], latticeVectors_[1]);
      den = latticeVectors_[2]*reciprocal;
      assert(den);
      reciprocalVectors_.push_back(reciprocal/den);
      if(latticeType_ == CRYSTAL)
        for(int i = -1; i <= 1; i++)
          for(int j = -1; j <= 1; j++)
            for(int k = -1; k <= 1; k++)
              imagePositions_.push_back(i*latticeVectors_[0] + j*latticeVectors_[1] + k*latticeVectors_[2]);
    }
    break;
  default:
    latticeVectors_.clear();
    std::cerr << "Error in Lattice::updateLattice - wrong number of lattice vectors" << std::endl;
  }
}

void Lattice::calculateLatticeVectors(Real const &a, Real const &b, Real const &c,
                                      Real const &alpha, Real const &beta, Real const &gamma)
{
  assert(a >= 0.0 && b >= 0.0 && c >= 0.0 && 
         alpha >= 0.0  && alpha < 180.0 &&
          beta >= 0.0  &&  beta < 180.0 &&
         gamma >= 0.0  && gamma < 180.0);
  latticeVectors_.clear();
  // Not periodic
  if(a == 0.0 && b == 0.0 && c == 0.0)
    return;
  // Linear cases
  else if(a != 0.0 && b == 0.0 && c == 0.0) 
  {
    latticeVectors_.push_back(Vector3D(a, 0.0, 0.0));
  }
  else if(a == 0.0 && b != 0.0 && c == 0.0) 
  {
    latticeVectors_.push_back(Vector3D(0.0, b, 0.0));
  }
  else if(a == 0.0 && b == 0.0 && c != 0.0) 
  {
    latticeVectors_.push_back(Vector3D(0.0, 0.0, c));
  }
  // Surface cases
  else if(a != 0.0 && b != 0.0 && c == 0.0)
  {
    assert(gamma != 0.0);
    latticeVectors_.push_back(Vector3D(a, 0.0, 0.0));
    latticeVectors_.push_back(Vector3D(b*cos(GRAD_TO_RAD*gamma), b*sin(GRAD_TO_RAD*gamma), 0.0));
  }
  else if(a != 0.0 && b == 0.0 && c != 0.0)
  {
    assert(beta != 0.0);
    latticeVectors_.push_back(Vector3D(a, 0.0, 0.0));
    latticeVectors_.push_back(Vector3D(c*cos(GRAD_TO_RAD*beta), 0.0, c*sin(GRAD_TO_RAD*beta)));
  }
  else if(a == 0.0 && b != 0.0 && c != 0.0) 
  {
    assert(alpha != 0.0);
    latticeVectors_.push_back(Vector3D(0.0, b, 0.0));
    latticeVectors_.push_back(Vector3D(0.0, c*cos(GRAD_TO_RAD*alpha), c*sin(GRAD_TO_RAD*alpha)));
  }
  // 3D case
  else
  {
    assert (alpha != 0.0 && beta != 0.0 && gamma != 0.0);
    Real const cosAlpha = cos(GRAD_TO_RAD*alpha); Real const cosBeta = cos(GRAD_TO_RAD*beta);
    Real const cosGamma = cos(GRAD_TO_RAD*gamma); Real const sinGamma = sin(GRAD_TO_RAD*gamma);
    Real const volume = a*b*c*
        sqrt(1 - cosAlpha*cosAlpha - cosBeta*cosBeta - cosGamma*cosGamma + 2.0*cosAlpha*cosBeta*cosGamma);
    latticeVectors_.push_back(Vector3D(a, 0.0, 0.0));
    latticeVectors_.push_back(Vector3D(b*cosGamma, b*sinGamma, 0.0));
    latticeVectors_.push_back(Vector3D(c*cosBeta, 
                                       c*(cosAlpha - cosBeta*cosGamma)/sinGamma, 
                                       volume/(a*b*sinGamma)));
  }
}
