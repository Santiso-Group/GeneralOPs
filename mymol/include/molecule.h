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
**
** A specialized Geometry that includes methods to calculate some
** properties, and to manipulate the geometry (see below).
**
** Note that this class, being derived from Geometry, does not include 
** any unit cell information, so it assumes that molecules have been 
** wrapped around the unit cell are joined when periodic boundary
** conditions are used (see the join and joinAll methods in system.h).
*/

// This class is under construction.
// Things to add:
// - Symmetry

#ifndef H_MOLECULE
#define H_MOLECULE

#include "common/include/types.h"
#include "mymath/include/vector3D.h"
#include "mymath/include/quaternion.h"
#include "mymol/include/geometry.h"

typedef std::vector<Vector3D> GradientVector;

class Molecule: public Geometry
{
public:

// Constructors

  Molecule();

// Destructor

  ~Molecule();

// Modified methods from base class

  void addAtom(Atom const &atom);                   // Adds an atom to the molecule
  void setPosition(size_t const index,              // Set the position of the index-th atom
                   Vector3D const &newPosition);    // As in class Geometry, it does not recalculate bonds
  void setx(size_t const index, Real const &newx);  // Set the x-coordinate of the index-th atom to newx
                                                    // As in class Geometry, it does not recalculate bonds
  void sety(size_t const index, Real const &newy);  // Set the y-coordinate of the index-th atom to newy
                                                    // As in class Geometry, it does not recalculate bonds
  void setz(size_t const index, Real const &newz);  // Set the z-coordinate of the index-th atom to newz
                                                    // As in class Geometry, it does not recalculate bonds
  void clear();                                     // Clears the molecule

// Interface

  void translate(Vector3D const &translation);  // Translate the molecule by the given vector
  void rotate(Quaternion const &quaternion,     // Rotate the molecule around the given
              Vector3D const &origin);          // origin by the given quaternion

// Accessors

  Real const &mass();                               // Returns the total mass of the molecule
  Real const &radiusOfGyration();                   // Returns the radius of gyration
  Real const &asphericity();                        // Returns the asphericity
  Real const &acylindricity();                      // Returns the acylindricity
  Real const &anisotropy();                         // Returns the relative shape anisotropy
  Vector3D const &center();                         // Returns the geometric center
  Vector3D const &centerOfMass();                   // Returns the center of mass
  Vector3D const &momentsOfInertia();               // Returns the principal moments of inertia
  Vector3D const &gyrationMoments();                // Returns the principal moments of gyration
  Matrix3D const &principalAxes();                  // Returns the principal axes as columns of a matrix
  Matrix3D const &gyrationAxes();                   // Returns the axes of gyration as columns of a matrix

// Calculate other properties

  Real const momentOfInertia(Vector3D const &axis);  // Calculates and returns the moment of inertia 
                                                     // around a specific axis
  Real const radiusOfGyration(Vector3D const &axis); // Calculates and returns the radius of gyration
                                                     // around a specific axis
  Real const distance(size_t const atom1,            // Calculates and returns the distance between
                      size_t const atom2) const;     // atoms atom1 and atom2
  Real const angle(size_t const atom1,               // Calculates and returns the angle (in radians)
                   size_t const atom2,               // defined by atoms atom1, atom2, and atom3.
                   size_t const atom3) const;        // The value returned is between 0 and pi.
  Real const dihedral(size_t const atom1,            // Calculates and returns the dihedral angle
                      size_t const atom2,            // atom1-atom2-atom3-atom4 (in radians)
                      size_t const atom3,            // The value returned is between -pi and pi.
                      size_t const atom4) const;

  // The following methods return the gradient of distances, angles, and dihedrals with respect to
  // the coordinates of each atom as a std::vector of Vector3Ds. The first component of the std::vector
  // is the gradient with respect to the first atom, the second one the gradient with respect to the
  // second atom, and so on.

  GradientVector const distanceGradients(size_t const atom1,        // Calculates and returns the gradients of the
                                         size_t const atom2) const; // distance between atom1 and atom2
  GradientVector const angleGradients(size_t const atom1,           // Calculates and returns the gradients of the
                                      size_t const atom2,           // angle defined by atom1, atom2 and atom3
                                      size_t const atom3) const;    // (in radians)
  GradientVector const dihedralGradients(size_t const atom1,        // Calculates and returns the gradients of the
                                         size_t const atom2,        // dihedral angle atom1-atom2-atom3-atom4
                                         size_t const atom3,        // (in radians)
                                         size_t const atom4) const;

private:

  Real mass_;                 // Mass of the molecule
  Real radiusOfGyration_;     // Radius of gyration
  Real asphericity_;          // Asphericity of the molecule
  Real acylindricity_;        // Acylindricity of the molecule
  Real anisotropy_;           // Relative shape anisotropy 
  Vector3D center_;           // Geometric center of the molecule
  Vector3D centerOfMass_;     // Center of mass
  Vector3D momentsOfInertia_; // Principal moments of inertia
  Vector3D gyrationMoments_;  // Principal moments of gyration
  Matrix3D principalAxes_;    // Principal axes
  Matrix3D gyrationAxes_;     // Axes of gyration
  bool knowMass_;             // Has the total mass been calculated?
  bool knowCenter_;           // Has the geometric center been calculated?
  bool knowCenterOfMass_;     // Has the center of mass been calculated?
  bool knowInertia_;          // Have the moments of inertia and principal axes been calculated?
  bool knowGyration_;         // Have the axes of gyration and related properties been calculated?

// Private functions to calculate various properties

  void calculateMass();         // Calculate the total mass
  void calculateCenter();       // Calculate the geometric center
  void calculateCenterOfMass(); // Calculate the center of mass
  void calculateInertia();      // Calculates the moments of inertia and principal axes
  void calculateGyration();     // Calculates the axes of gyration and related properties
};

/*
** End of class Molecule
*/ 

// Inlines

inline void Molecule::addAtom(Atom const &atom)
{
  Geometry::addAtom(atom);
  knowMass_ = knowCenter_ = knowCenterOfMass_ = knowInertia_ = false;
}

inline void Molecule::setPosition(size_t const index, Vector3D const &newPosition)
{
  Geometry::setPosition(index, newPosition);
  knowCenter_ = knowCenterOfMass_ = knowInertia_ = false;
}
 
inline void Molecule::setx(size_t const index, Real const &newx)
{
  Geometry::setx(index, newx);
  knowCenter_ = knowCenterOfMass_ = knowInertia_ = false;
}

inline void Molecule::sety(size_t const index, Real const &newy)
{
  Geometry::sety(index, newy);
  knowCenter_ = knowCenterOfMass_ = knowInertia_ = false;
}

inline void Molecule::setz(size_t const index, Real const &newz)
{
  Geometry::setz(index, newz);
  knowCenter_ = knowCenterOfMass_ = knowInertia_ = false;
}

inline void Molecule::clear()
{
  Geometry::clear();
  knowMass_ = knowCenter_ = knowCenterOfMass_ = knowInertia_ = false;
}

inline Real const &Molecule::mass()
{
  if(!knowMass_) calculateMass();
  return mass_;
}

inline Real const &Molecule::radiusOfGyration()
{
  if(!knowGyration_) calculateGyration();
  return radiusOfGyration_; 
}

inline Real const &Molecule::asphericity()
{
  if(!knowGyration_) calculateGyration();
  return asphericity_;
}

inline Real const &Molecule::acylindricity()
{
  if(!knowGyration_) calculateGyration();
  return acylindricity_;
}

inline Real const &Molecule::anisotropy()
{
  if(!knowGyration_) calculateGyration();
  return anisotropy_;
}

inline Vector3D const &Molecule::center()
{
  if(!knowCenter_) calculateCenter();
  return center_;
}

inline Vector3D const &Molecule::centerOfMass()
{
  if(!knowCenterOfMass_) calculateCenterOfMass();
  return centerOfMass_;
}

inline Vector3D const &Molecule::momentsOfInertia()
{
  if(!knowInertia_) calculateInertia();
  return momentsOfInertia_;
}

inline Vector3D const &Molecule::gyrationMoments()
{
  if(!knowGyration_) calculateGyration();
  return gyrationMoments_;
}

inline Matrix3D const &Molecule::principalAxes()
{
  if(!knowInertia_) calculateInertia();
  return principalAxes_;
}

inline Matrix3D const &Molecule::gyrationAxes()
{
  if(!knowGyration_) calculateGyration();
  return gyrationAxes_;
}

inline Real const Molecule::distance(size_t const atom1, size_t const atom2) const
{
  return norm(atom(atom2).position - atom(atom1).position);
}

#endif

