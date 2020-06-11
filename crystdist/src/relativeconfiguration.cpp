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
** Relative configuration class
*/

#include <cmath>
#include "crystdist/include/relativeconfiguration.h"

// Constructors

RelativeConfiguration::RelativeConfiguration()
:
types(GENERAL, GENERAL), name(), internalDOFs(), distance(0.0), 
bondOrientation(1.0), relativeOrientation(1.0)
{}

RelativeConfiguration::RelativeConfiguration(PointMolecule const &first,
                                             PointMolecule const &second,
                                             Lattice const &lattice)
:
types(first.type, second.type), internalDOFs(), distance(0.0), 
bondOrientation(1.0), relativeOrientation(1.0)
{ calculate(first, second, lattice); }

void RelativeConfiguration::calculate(PointMolecule const &first, 
                                      PointMolecule const &second,
                                      Lattice const &lattice)
{
  if((first.type != second.type) || (first.name != second.name))
  {
    std::cerr << "Error in RelativeConfiguration::calculate: "
              << "Hybrid types not implemented yet." << std::endl;
    clear();
    return;
  }

  types.first = first.type;
  types.second = second.type;
  name = first.name + "-" + second.name;
  internalDOFs.first = first.internalDOFs;
  internalDOFs.second = second.internalDOFs;

  if(first.type == GENERAL && second.type == GENERAL)
  {
    // General case (quaternion vs. quaternion):
    //
    // Bond vector = q_i_*r_ij*q_i, where r_ij is the bond in the lab frame
    // Relative orientation = q_i_*q_j (on the i-th molecule's frame)
    // -> equivalent to q_i*q_j_ on the lab frame.
    
    Vector3D bond = lattice.difference(first.position, second.position);
    rotateVector(~first.orientation, bond); // Project onto first molecule frame
    distance = norm(bond);
    bondOrientation = bond/distance;
    relativeOrientation = ~first.orientation*second.orientation;
  }
  else if(first.type == LINEAR_ASYMMETRIC && second.type == LINEAR_ASYMMETRIC)
  {
    // Vector vs. vector:
    //
    // Bond orientation = Angle between r_ij and v_i
    // Relative orientation = Angle between v_i and v_j (for a different way
    // see http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/index.htm
    // - that way is not invariant w.r.t. a change of coordinate frame, though).

    Vector3D bond = lattice.difference(first.position, second.position);
    Vector3D v_i = vector(first.orientation);
    Vector3D v_j = vector(second.orientation);
    distance = norm(bond);
    Real dotProduct = bond*v_i/distance;
    if(dotProduct > 1.0) dotProduct = 1.0;   // Avoid NaN due to
    if(dotProduct < -1.0) dotProduct = -1.0; // roundoff
    bondOrientation = acos(dotProduct);
    dotProduct = v_i*v_j;
    if(dotProduct > 1.0) dotProduct = 1.0;   // Avoid NaN due to
    if(dotProduct < -1.0) dotProduct = -1.0; // roundoff
    relativeOrientation = acos(dotProduct);
    //bondOrientation = (v_i*bond/distance)*v_i;
    //relativeOrientation = 1.0 + v_i*v_j + cross(v_i, v_j);
    //relativeOrientation.normalize();
  }
  else if((first.type == LINEAR_SYMMETRIC && second.type == LINEAR_SYMMETRIC) ||
          (first.type == PLANAR_SYMMETRIC && second.type == PLANAR_SYMMETRIC))
  {
    // Axis vs. axis:
    //
    // Bond orientation and relative orientation are defined like in the 
    // vector case, only in this case they are between 0 and 90 degrees

    Vector3D bond = lattice.difference(first.position, second.position);
    Vector3D v_i = vector(first.orientation);
    Vector3D v_j = vector(second.orientation);
    distance = norm(bond);
    Real dotProduct = fabs(bond*v_i/distance);
    if(dotProduct > 1.0) dotProduct = 1.0;   // Avoid NaN due to roundoff
    bondOrientation = acos(dotProduct);
    dotProduct = fabs(v_i*v_j);
    if(dotProduct > 1.0) dotProduct = 1.0;   // Avoid NaN due to roundoff
    relativeOrientation = acos(dotProduct);
    //bondOrientation = (v_i*bond/distance)*v_i;
    //relativeOrientation = v_i*v_j + cross(v_i, v_j);
    //relativeOrientation.normalize();
  }
  else  // Should never get here
  {
    std::cerr << "Error in RelativeConfiguration::calculate: "
              << "Unknown combination of point molecule types." << std::endl;
    clear();
    return;
  }
}

Real const squareDeviation(RelativeConfiguration const &conf1, RelativeConfiguration const &conf2)
{
  assert( (conf1.types == conf2.types) ||
          (conf1.name == conf2.name) ||
          (conf1.internalDOFs.first.size() == conf2.internalDOFs.first.size()) ||
          (conf1.internalDOFs.second.size() == conf2.internalDOFs.second.size()) );
  Real norm2 = 0.0;
  Real diffDist = conf1.distance - conf2.distance;
  norm2 += diffDist*diffDist;
  if(conf1.types.first == GENERAL && conf1.types.second == GENERAL)
  {
    Real dotProduct = vector(conf1.bondOrientation)*vector(conf2.bondOrientation);
    norm2 += 1 - dotProduct;
    dotProduct = dot(conf1.relativeOrientation, conf2.relativeOrientation);
    norm2 += 1 - fabs(dotProduct);
  }
  else if(conf1.types.first == LINEAR_ASYMMETRIC && conf1.types.second == LINEAR_ASYMMETRIC)
  {
    Real angleDiff = scalar(conf1.bondOrientation) - scalar(conf2.bondOrientation);
    Real diff = 1.0 - cos(angleDiff);
    norm2 += diff*diff;
    angleDiff = scalar(conf1.relativeOrientation) - scalar(conf2.relativeOrientation);
    diff = 1.0 - cos(angleDiff);
    norm2 += diff*diff;
  }
  else if((conf1.types.first == LINEAR_SYMMETRIC && conf1.types.second == LINEAR_SYMMETRIC) ||
          (conf1.types.first == PLANAR_SYMMETRIC && conf1.types.second == PLANAR_SYMMETRIC))
  {
    Real angleDiff = 2.0*(scalar(conf1.bondOrientation) - scalar(conf2.bondOrientation)); // 2 due to symmetry
    Real diff = 1.0 - cos(angleDiff);
    norm2 += diff*diff;
    angleDiff = 2.0*(scalar(conf1.relativeOrientation) - scalar(conf2.relativeOrientation)); // Ditto above
    diff = 1.0 - cos(angleDiff);
    norm2 += diff*diff;
  }
  else  // Should never get here
  {
    std::cerr << "Internal error in squareDeviation: combination of molecule types not implemented" << std::endl;
    return 0.0;
  }
  for(size_t i = 0; i < conf1.internalDOFs.first.size(); i++)
  {
    if(conf1.internalDOFs.first[i].type == DISTANCE)
    {
      Real diff = conf1.internalDOFs.first[i].value - conf2.internalDOFs.first[i].value;
      norm2 += diff*diff;
    }
    else if(conf1.internalDOFs.first[i].type == ANGLE)
    {
      Real angleDiff = conf1.internalDOFs.first[i].value - conf2.internalDOFs.first[i].value;
      Real diff = 1.0 - cos(angleDiff);
      norm2 += diff*diff;
    }
    else if(conf1.internalDOFs.first[i].type == DIHEDRAL)
    {
      Real angleDiff = (conf1.internalDOFs.first[i].value - conf2.internalDOFs.first[i].value)*
                       (Real)conf1.internalDOFs.first[i].symmetryNumber;
      Real diff = 1.0 - cos(angleDiff);
      norm2 += diff*diff;
    }
    else  // Should never get here
    {
      std::cerr << "Internal error in squareDeviation: internal degree of freedom not implemented" << std::endl;
      return 0.0;
    }
  }
  for(size_t i = 0; i < conf1.internalDOFs.second.size(); i++)
  {
    if(conf1.internalDOFs.second[i].type == DISTANCE)
    {
      Real diff = conf1.internalDOFs.second[i].value - conf2.internalDOFs.second[i].value;
      norm2 += diff*diff;
    }
    else if(conf1.internalDOFs.second[i].type == ANGLE)
    {
      Real angleDiff = conf1.internalDOFs.second[i].value - conf2.internalDOFs.second[i].value;
      Real diff = 1.0 - cos(angleDiff);
      norm2 += diff*diff;
    }
    else if(conf1.internalDOFs.second[i].type == DIHEDRAL)
    {
      Real angleDiff = (conf1.internalDOFs.second[i].value - conf2.internalDOFs.second[i].value)*
                       (Real)conf1.internalDOFs.second[i].symmetryNumber;
      Real diff = 1.0 - cos(angleDiff);
      norm2 += diff*diff;
    }
    else  // Should never get here
    {
      std::cerr << "Internal error in squareDeviation: internal degree of freedom not implemented" << std::endl;
      return 0.0;
    }
  }
  return norm2;
}

bool operator<(RelativeConfiguration const &conf1, RelativeConfiguration const &conf2)
{
  // Arbitrary ordering used for sorting peak tables
  if(conf1.distance < conf2.distance) return true;
  if(conf1.distance == conf2.distance)
  {
    Real norm1 = norm(conf1.bondOrientation);
    Real norm2 = norm(conf2.bondOrientation);
    if(norm1 < norm2) return true;
    if(norm1 == norm2)
    {
      norm1 = norm(conf1.relativeOrientation);
      norm2 = norm(conf2.relativeOrientation);
      if(norm1 < norm2) return true;
      if(norm1 == norm2)
      {
        if(conf1.internalDOFs.first.size() > 0 &&
           conf2.internalDOFs.first.size() > 0)
        {
          norm1 = conf1.internalDOFs.first[0].value;
          norm2 = conf2.internalDOFs.first[0].value;
          if(norm1 < norm2) return true;
          if(norm1 == norm2)
          {
            if(conf1.internalDOFs.second.size() > 0 &&
               conf2.internalDOFs.second.size() > 0)
            {
              norm1 = conf1.internalDOFs.second[0].value;
              norm2 = conf2.internalDOFs.second[1].value;
              if(norm1 < norm2) return true;
            }
          }
        }
      }
    }
  }
  return false;
}

std::istream &operator>>(std::istream &inStream, RelativeConfiguration &conf)
{
  conf.clear();
  inStream >> conf.types.first >> conf.types.second;
  inStream >> conf.name;
  inStream >> conf.distance;
  inStream >> conf.bondOrientation.w >> conf.bondOrientation.x
           >> conf.bondOrientation.y >> conf.bondOrientation.z;
  inStream >> conf.relativeOrientation.w >> conf.relativeOrientation.x
           >> conf.relativeOrientation.y >> conf.relativeOrientation.z;
  size_t numIntDOFs;
  inStream >> numIntDOFs;
  conf.internalDOFs.first.clear();
  for(size_t i = 0; i < numIntDOFs; i++)
  {
    InternalDOF currentDOF;
    inStream >> currentDOF.type
             >> currentDOF.value
             >> currentDOF.symmetryNumber;
    if(!inStream)
    {
      std::cerr << "Error reading internal degree of freedom" << std::endl;
      return inStream;
    }
    conf.internalDOFs.first.push_back(currentDOF);
  }
  inStream >> numIntDOFs;
  conf.internalDOFs.second.clear();
  for(size_t i = 0; i < numIntDOFs; i++)
  {
    InternalDOF currentDOF;
    inStream >> currentDOF.type
             >> currentDOF.value
             >> currentDOF.symmetryNumber;
    if(!inStream)
    {
      std::cerr << "Error reading internal degree of freedom" << std::endl;
      return inStream;
    }
    conf.internalDOFs.second.push_back(currentDOF);
  }
  return inStream;
}

std::ostream &operator<<(std::ostream &outStream, RelativeConfiguration const &conf)
{
  // Short version
  outStream << conf.types.first << " " << conf.types.second << " "
            << ((conf.name.size() != 0)?conf.name:".") << " " << conf.distance << " "
            << conf.bondOrientation.w << " " << conf.bondOrientation.x << " "
            << conf.bondOrientation.y << " " << conf.bondOrientation.z << " "
            << conf.relativeOrientation.w << " " << conf.relativeOrientation.x << " "
            << conf.relativeOrientation.y << " " << conf.relativeOrientation.z << " ";
  outStream << conf.internalDOFs.first.size() << " ";
  for(size_t i = 0; i < conf.internalDOFs.first.size(); i++)
    outStream << conf.internalDOFs.first[i].type << " "
              << conf.internalDOFs.first[i].value << " "
              << conf.internalDOFs.first[i].symmetryNumber << " ";
  outStream << conf.internalDOFs.second.size() << " ";
  for(size_t i = 0; i < conf.internalDOFs.second.size(); i++)
    outStream << conf.internalDOFs.second[i].type << " "
              << conf.internalDOFs.second[i].value << " "
              << conf.internalDOFs.second[i].symmetryNumber << " ";

  // The long version below was used for debugging, it is no longer used
  /*
  outStream << std::endl << "Name: " << conf.name
            << std::endl << "Distance: " << conf.distance
            << std::endl << "Bond orientation: " << conf.bondOrientation
            << std::endl << "Relative orientation: " << conf.relativeOrientation;
  outStream << std::endl << "Internal degrees of freedom - first:";
  for(size_t i = 0; i < conf.internalDOFs.first.size(); i++)
  {
    outStream << std::endl << i << " - Type: ";
    switch(conf.internalDOFs.first[i].type)
    {
      case DISTANCE:
        outStream << "Distance.";
        break;
      case ANGLE:
        outStream << "Angle.";
        break;
      case DIHEDRAL:
        outStream << "Dihedral - symmetry number: "
                  << conf.internalDOFs.first[i].symmetryNumber;
        break;
      default:
        outStream << "Unknown! Internal error";
    }
    outStream << std::endl << "Value: " << conf.internalDOFs.first[i].value;
  }
  outStream << std::endl << "Internal degrees of freedom - second:";
  for(size_t i = 0; i < conf.internalDOFs.second.size(); i++)
  {
    outStream << std::endl << i << " - Type: ";
    switch(conf.internalDOFs.second[i].type)
    {
      case DISTANCE:
        outStream << "Distance.";
        break;
      case ANGLE:
        outStream << "Angle.";
        break;
      case DIHEDRAL:
        outStream << "Dihedral - symmetry number: "
                  << conf.internalDOFs.second[i].symmetryNumber;
        break;
      default:
        outStream << "Unknown! Internal error";
    }
    outStream << std::endl << "Value: " << conf.internalDOFs.second[i].value;
  }
  */
  return outStream;
}
