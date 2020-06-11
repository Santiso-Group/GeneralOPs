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
** Molecule class
*/

#include "mymath/include/symmetriceigensystem.h"
#include "mymol/include/molecule.h"

// Constructors

Molecule::Molecule()
:
Geometry(), mass_(0), center_(0.0), centerOfMass_(0.0), momentsOfInertia_(0.0),
principalAxes_(0.0), knowMass_(false), knowCenter_(false), knowCenterOfMass_(false),
knowInertia_(false)
{}

// Destructor

Molecule::~Molecule()
{}

// Interface

void Molecule::translate(Vector3D const &translation)
{
  for(size_t i = 0; i < numAtoms(); ++i)
    setPosition(i, atom(i).position + translation);
  if(knowCenter_) center_ += translation;
  if(knowCenterOfMass_) centerOfMass_ += translation;
}

void Molecule::rotate(Quaternion const &quaternion, Vector3D const &origin)
{
  for(size_t i = 0; i < numAtoms(); ++i)
    setPosition(i, origin + rotatedVector(quaternion, atom(i).position - origin));
  if(knowCenter_) center_ = origin + rotatedVector(quaternion, center_ - origin);
  if(knowCenterOfMass_) centerOfMass_ = origin + rotatedVector(quaternion, centerOfMass_ - origin);
  if(knowInertia_)
  {
    Vector3D v1 = rotatedVector(quaternion, principalAxes_.column(0));
    Vector3D v2 = rotatedVector(quaternion, principalAxes_.column(1));
    Vector3D v3 = rotatedVector(quaternion, principalAxes_.column(2));
    principalAxes_ = Matrix3D(v1, v2, v3);
  }
  if(knowGyration_)
  {
    Vector3D v1 = rotatedVector(quaternion, gyrationAxes_.column(0));
    Vector3D v2 = rotatedVector(quaternion, gyrationAxes_.column(1));
    Vector3D v3 = rotatedVector(quaternion, gyrationAxes_.column(2));
    gyrationAxes_ = Matrix3D(v1, v2, v3);
  }
}

Real const Molecule::momentOfInertia(Vector3D const &axis)
{
  if(!knowCenterOfMass_) calculateCenterOfMass();
  Vector3D unitAxis = unit(axis);
  Real momentOfInertia = 0.0;
  for(size_t i = 0; i < numAtoms(); ++i)
  {
    Vector3D position = atom(i).position - centerOfMass_;
    momentOfInertia += atom(i).mass*position*(position - unitAxis);
  }
  return momentOfInertia;
}

Real const Molecule::radiusOfGyration(Vector3D const &axis)
{
  return sqrt(momentOfInertia(axis)/mass_); // Mass will be calculated in momentOfInertia()
}

Real const Molecule::angle(size_t const atom1, size_t const atom2, size_t const atom3) const
{
  Vector3D v12 = atom(atom1).position - atom(atom2).position;
  Vector3D v32 = atom(atom3).position - atom(atom2).position;
  Real arg = v12*v32/(norm(v12)*norm(v32));
  if(arg < -1.0) arg = -1.0;  // Avoid roundoff problems
  if(arg > 1.0) arg = 1.0;
  return acos(arg);
}

Real const Molecule::dihedral(size_t const atom1, size_t const atom2, size_t const atom3, size_t const atom4) const
{
  Vector3D v12 = atom(atom1).position - atom(atom2).position;
  Vector3D v23 = atom(atom2).position - atom(atom3).position;
  Vector3D v34 = atom(atom3).position - atom(atom4).position;
  Vector3D n1 = cross(v12, v23);
  Vector3D n2 = cross(v23, v34);
  Vector3D n3 = cross(v23, n1);
  Real normn2 = norm(n2);
  return -atan2(n3*n2/(normn2*norm(n3)), n1*n2/(normn2*norm(n1)));
}

GradientVector const Molecule::distanceGradients(size_t const atom1, size_t const atom2) const
{
  GradientVector distanceGradients;
  Vector3D grad2 = unit(atom(atom2).position - atom(atom1).position);
  distanceGradients.push_back(-grad2);
  distanceGradients.push_back(grad2);
  return distanceGradients;
}

GradientVector const Molecule::angleGradients(size_t const atom1, size_t const atom2, size_t const atom3) const
{
  GradientVector angleGradients;
  Vector3D v12 = atom(atom1).position - atom(atom2).position;
  Vector3D v32 = atom(atom3).position - atom(atom2).position;
  Real invd12 = 1.0/norm(v12); Real invd32 = 1.0/norm(v32);
  Real cosTheta = (invd12*invd32)*(v12*v32);
  if(cosTheta < -1.0) cosTheta = -1.0;  // Avoid roundoff problems
  if(cosTheta > 1.0) cosTheta = 1.0;
  Real invSinTheta = 1.0/sqrt(1.0 - cosTheta*cosTheta);
  Vector3D grad1 = invd12*invSinTheta*((invd12*cosTheta)*v12 - invd32*v32);
  Vector3D grad3 = invd32*invSinTheta*((invd32*cosTheta)*v32 - invd12*v12);
  angleGradients.push_back(grad1);
  angleGradients.push_back(-grad1 - grad3);
  angleGradients.push_back(grad3);
  return angleGradients;
}

GradientVector const Molecule::dihedralGradients(size_t const atom1, size_t const atom2, 
                                                 size_t const atom3, size_t const atom4) const
{
  GradientVector dihedralGradients;
  Vector3D v12 = atom(atom1).position - atom(atom2).position;
  Vector3D v23 = atom(atom2).position - atom(atom3).position;
  Vector3D v34 = atom(atom3).position - atom(atom4).position;
  Vector3D n1 = cross(v12, v23);
  Vector3D n2 = cross(v23, v34);
  Vector3D n3 = cross(v23, n1);
  Real invdn1 = 1.0/norm(n1); Real invdn2 = 1.0/norm(n2); Real invdn3 = 1.0/norm(n3);
  Real cosTheta = (invdn1*invdn2)*(n1*n2); Real sinTheta = (invdn3*invdn2)*(n3*n2);
  Vector3D grad1, grad12, grad3;
  if(fabs(sinTheta) > 0.1)
  {
    n1 *= invdn1; n2 *= invdn2;
    Vector3D dcosThetadn1 = invdn1*(cosTheta*n1 - n2);
    Vector3D dcosThetadn2 = invdn2*(cosTheta*n2 - n1);
    Real invSin = 1.0/sinTheta;
    grad1 = -invSin*cross(v23, dcosThetadn1);
    grad3 = -invSin*cross(v23, dcosThetadn2);
    grad12 = invSin*(cross(v12, dcosThetadn1) - cross(v34, dcosThetadn2));
  }
  else
  {
    n2 *= invdn2; n3 *= invdn3;
    Vector3D dsinThetadn2 = n2*(sinTheta*n2 - n3);
    Vector3D dsinThetadn3 = n3*(sinTheta*n3 - n2);
    Real invCos = 1.0/cosTheta;
    Vector3D v23crossdsinThetadn3 = cross(v23, dsinThetadn3);
    grad1 = invCos*cross(v23crossdsinThetadn3, v23);
    grad3 = -invCos*cross(dsinThetadn2, v23);
    grad12 = invCos*(cross(v12, v23crossdsinThetadn3) + cross(n1, dsinThetadn3) + cross(v34, dsinThetadn2));
  }
  dihedralGradients.push_back(grad1);
  dihedralGradients.push_back(grad12 - grad1);
  dihedralGradients.push_back(-grad3 - grad12);
  dihedralGradients.push_back(grad3);
  return dihedralGradients;
}

// Private functions to calculate various properties

void Molecule::calculateMass()
{
  mass_ = 0.0;
  for(size_t i = 0; i < numAtoms(); ++i)
    mass_ += atom(i).mass;
  knowMass_ = true;
}

void Molecule::calculateCenter()
{
  center_ = 0.0;
  size_t nAtom = numAtoms();
  for(size_t i = 0; i < nAtom; ++i)
    center_ += atom(i).position;
  center_ /= (nAtom > 0)?(Real)nAtom:1.0;
  knowCenter_ = true;
}

void Molecule::calculateCenterOfMass()
{
  if(!knowMass_) calculateMass();
  centerOfMass_ = 0.0;
  for(size_t i = 0; i < numAtoms(); ++i)
    centerOfMass_ += atom(i).mass*atom(i).position;
  centerOfMass_ /= (mass_ > 0.0)?mass_:1.0;
  knowCenterOfMass_ = true;
}

void Molecule::calculateInertia()
{
  if(!knowCenterOfMass_) calculateCenterOfMass();
  Matrix3D inertiaTensor = 0.0;
  for(size_t i = 0; i < numAtoms(); ++i)
  {
    Real mass = atom(i).mass;
    Vector3D position = atom(i).position - centerOfMass_;
    Real x2 = position.x*position.x;
    Real y2 = position.y*position.y;
    Real z2 = position.z*position.z;
    inertiaTensor(0, 0) += mass*(y2 + z2);
    inertiaTensor(1, 1) += mass*(x2 + z2);
    inertiaTensor(2, 2) += mass*(x2 + y2);
    inertiaTensor(0, 1) -= mass*(position.x*position.y);
    inertiaTensor(0, 2) -= mass*(position.x*position.z);
    inertiaTensor(1, 2) -= mass*(position.y*position.z);
  }
  inertiaTensor(1, 0) = inertiaTensor(0, 1);
  inertiaTensor(2, 0) = inertiaTensor(0, 2);
  inertiaTensor(2, 1) = inertiaTensor(1, 2);
  SymmetricEigensystem<Matrix3D, Vector3D> inertia(inertiaTensor);
  momentsOfInertia_ = inertia.eigenvalues();
  principalAxes_ = inertia.eigenvectors();
  knowInertia_ = true;
}

void Molecule::calculateGyration()
{
  if(!knowCenter_) calculateCenter();
  Matrix3D gyrationTensor = 0.0;
  for(size_t i = 0; i < numAtoms(); ++i)
  {
    Vector3D position = atom(i).position - center_;
    gyrationTensor(0, 0) += position.x*position.x;
    gyrationTensor(0, 1) += position.x*position.y;
    gyrationTensor(0, 2) += position.x*position.z;
    gyrationTensor(1, 1) += position.y*position.y;
    gyrationTensor(1, 2) += position.y*position.z;
    gyrationTensor(2, 2) += position.z*position.z;
  }
  gyrationTensor(1, 0) = gyrationTensor(0, 1);
  gyrationTensor(2, 0) = gyrationTensor(0, 2);
  gyrationTensor(2, 1) = gyrationTensor(1, 2);
  SymmetricEigensystem<Matrix3D, Vector3D> gyration(gyrationTensor);
  gyrationMoments_ = gyration.eigenvalues();
  gyrationAxes_ = gyration.eigenvectors();
  Real const radiusOfGyration2 = gyrationMoments_[0] + gyrationMoments_[1] + gyrationMoments_[2];
  radiusOfGyration_ = sqrt(radiusOfGyration2);
  asphericity_ = gyrationMoments_[2] - 0.5*(gyrationMoments_[0] + gyrationMoments_[1]);
  acylindricity_ = gyrationMoments_[1] - gyrationMoments_[0];
  anisotropy_ = (asphericity_*asphericity_ + 0.75*acylindricity_*acylindricity_)/
                (radiusOfGyration2*radiusOfGyration2);
  knowGyration_ = true;
}

