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
** Point molecule class
*/

#include "crystdist/include/pointmolecule.h"

// Constructors

PointMolecule::PointMolecule(MoleculeType const &type)
: 
name(), position(0.0), orientation(1.0), type(type), internalDOFs()
{}

PointMolecule::PointMolecule(Molecule &molecule, MoleculeMap const &map)
:
name(), position(0.0), orientation(1.0), type(GENERAL), internalDOFs()
{ set(molecule, map); }

// Interface

void PointMolecule::set(Molecule &molecule, MoleculeMap const &map)
{
  if(molecule.name().find(map.name) == molecule.name().npos)
    std::cerr << "Warning from PointMolecule::setPointMolecule: " << std::endl
              << "Molecule and map residue name are not the same." << std::endl;

  // Name
  name = molecule.name();

  // Remove trailing space (e.g. from pdb residue name)
  size_t indx = name.find_last_not_of(" \n") + 1;
  size_t len = name.size() - indx;
  if(len > 0) name.erase(indx, len);  

  // Type
  type = map.type;

  // Position
  position = molecule.centerOfMass();

  // Orientation
  Vector3D origin = 0.0, xpoint = 0.0, ypoint = 0.0;
  switch(map.type)
  {
  case LINEAR_SYMMETRIC:
  case LINEAR_ASYMMETRIC:
    assert(map.framePoints.size() == 2);
    for(size_t i = 0; i < map.framePoints[0].numAtoms(); ++i)
      origin += map.framePoints[0].weights[i]*
                molecule.atom(map.framePoints[0].atoms[i]).position;
    for(size_t i = 0; i < map.framePoints[1].numAtoms(); ++i)
      xpoint += map.framePoints[1].weights[i]*
                molecule.atom(map.framePoints[1].atoms[i]).position;
    orientation = xpoint - origin;
    orientation.normalize();
    break;
  case PLANAR_SYMMETRIC:
    {
      assert(map.framePoints.size() == 3);
      for(size_t i = 0; i < map.framePoints[0].numAtoms(); ++i)
        origin += map.framePoints[0].weights[i]*
                  molecule.atom(map.framePoints[0].atoms[i]).position;
      for(size_t i = 0; i < map.framePoints[1].numAtoms(); ++i)
        xpoint += map.framePoints[1].weights[i]*
                  molecule.atom(map.framePoints[1].atoms[i]).position;
      for(size_t i = 0; i < map.framePoints[2].numAtoms(); ++i)
        ypoint += map.framePoints[2].weights[i]*
                  molecule.atom(map.framePoints[2].atoms[i]).position;
      Vector3D vx = xpoint - origin;
      vx.normalize();
      Vector3D vy = ypoint - origin;
      vy -= (vx*vy)*vx;
      vy.normalize();
      orientation = cross(vx, vy);
    }
    break;
  case GENERAL:
    {
      assert(map.framePoints.size() == 3);
      for(size_t i = 0; i < map.framePoints[0].numAtoms(); ++i)
        origin += map.framePoints[0].weights[i]*
                  molecule.atom(map.framePoints[0].atoms[i]).position;
      for(size_t i = 0; i < map.framePoints[1].numAtoms(); ++i)
        xpoint += map.framePoints[1].weights[i]*
                  molecule.atom(map.framePoints[1].atoms[i]).position;
      for(size_t i = 0; i < map.framePoints[2].numAtoms(); ++i)
        ypoint += map.framePoints[2].weights[i]*
                  molecule.atom(map.framePoints[2].atoms[i]).position;
      Vector3D vx;
      vx = xpoint - origin;
      vx.normalize();
      Vector3D vy = ypoint - origin;
      vy -= (vx*vy)*vx;
      vy.normalize();
      Vector3D vz = cross(vx, vy);
      orientation = Quaternion(Matrix3D(vx, vy, vz));
    }
    break;
  }
  // Internal degrees of freedom
  internalDOFs.clear();
  for(size_t i = 0; i < map.internalDOFs.size(); i++)
  {
    InternalDOF currentInternalDOF;
    currentInternalDOF.type = map.internalDOFs[i].type;
    switch(currentInternalDOF.type)
    {
    case DISTANCE:
      assert(map.internalDOFs[i].atoms.size() == 2);
      currentInternalDOF.value = molecule.distance(map.internalDOFs[i].atoms[0], 
                                                   map.internalDOFs[i].atoms[1]);
      break;
    case ANGLE:
      assert(map.internalDOFs[i].atoms.size() == 3);
      currentInternalDOF.value = molecule.angle(map.internalDOFs[i].atoms[0],
                                                map.internalDOFs[i].atoms[1],
                                                map.internalDOFs[i].atoms[2]);
      break;
    case DIHEDRAL:
      assert(map.internalDOFs[i].atoms.size() == 4);
      currentInternalDOF.value = molecule.dihedral(map.internalDOFs[i].atoms[0],
                                                   map.internalDOFs[i].atoms[1],
                                                   map.internalDOFs[i].atoms[2],
                                                   map.internalDOFs[i].atoms[3]);
      currentInternalDOF.symmetryNumber = map.internalDOFs[i].symmetryNumber;
      break;
    default:  // Should never get here
      std::cerr << "Internal error!" << std::endl
                << "Unknown internal degree of freedom in PointMolecule::setPointMolecule" << std::endl;
    }
    internalDOFs.push_back(currentInternalDOF);
  }
}

// Build a System<PointMolecule> from a System<Molecule> and a molecule map

System<PointMolecule> const getPointMolecules(System<Molecule> &system, MoleculeMap const &map)
{
  System<PointMolecule> mySystem;
  mySystem.setName(map.name + ":" + system.name());
  mySystem.setLattice(system.lattice());
  for(size_t i = 0; i < system.size(); i++)
    if(system[i].name().find(map.name) != system[i].name().npos)
      mySystem.add(PointMolecule(system[i], map));
  return mySystem;
}

// Print point molecule information (mostly for debugging)

std::ostream &operator<<(std::ostream &outStream, PointMolecule const &pointMolecule)
{
  outStream << std::endl << "Point molecule: " << pointMolecule.name;
  outStream << std::endl << "Type: ";

  switch(pointMolecule.type)
  {
    case LINEAR_SYMMETRIC:
      outStream << "Linear and symmetric.";
      break;
    case PLANAR_SYMMETRIC:
      outStream << "Planar and symmetric.";
      break;
    case LINEAR_ASYMMETRIC:
      outStream << "Linear and asymmetric.";
      break;
    default:
      outStream << "General.";
  }

  outStream << std::endl << "Position: " << pointMolecule.position
            << std::endl << "Orientation: " << pointMolecule.orientation;
  outStream << std::endl << "Internal degrees of freedom: ";
  
  for(size_t i = 0; i < pointMolecule.numInternalDOFs(); i++)
  {
    outStream << std::endl << i << " - Type: ";
    switch(pointMolecule.internalDOF(i).type)
    {
      case DISTANCE:
        outStream << "Distance.";
        break;
      case ANGLE:
        outStream << "Angle.";
        break;
      case DIHEDRAL:
        outStream << "Dihedral - symmetry number: "
                  << pointMolecule.internalDOF(i).symmetryNumber;
        break;
      default:
        outStream << "Unknown! Internal error";
    }
    outStream << std::endl << "Value: " << pointMolecule.internalDOF(i).value;
  }
  return outStream;
}
