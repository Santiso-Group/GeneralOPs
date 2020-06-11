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
** Molecule map class.
**
** A molecule map is used to convert molecules to point molecule
** objects. The map contains the name of the molecule, the
** local indices of the atoms defining the molecule-centered
** frame, and information on the internal degrees of freedom -
** their type, and the local indices of the atoms that define 
** each one of them.
**
** Note that all the data members are public.
*/

#ifndef H_MOLECULE_MAP
#define H_MOLECULE_MAP

#include <string>
#include <vector>
#include "common/include/assert.h"

// Internal degree of freedom types and molecule types

enum InternalDOFType { DISTANCE, ANGLE, DIHEDRAL };
enum MoleculeType { LINEAR_SYMMETRIC, PLANAR_SYMMETRIC, LINEAR_ASYMMETRIC, GENERAL };

// I/O for DOF and Molecule types

inline std::ostream &operator<<(std::ostream &outStream, InternalDOFType const &type)
{
  switch(type)
  {
    case DISTANCE:
      outStream << "DISTANCE";
      break;
    case ANGLE:
      outStream << "ANGLE";
      break;
    case DIHEDRAL:
      outStream << "DIHEDRAL";
      break;
  }
  return outStream;
}

inline std::istream &operator>>(std::istream &inStream, InternalDOFType &type)
{
  std::string r_type; // Type read from stream
  inStream >> r_type;
  if(r_type.find("DISTANCE") != r_type.npos)
    type = DISTANCE;
  else if(r_type.find("ANGLE") != r_type.npos)
    type = ANGLE;
  else if(r_type.find("DIHEDRAL") != r_type.npos)
    type = DIHEDRAL;
  else
  {
    std::cerr << "I/O error: Invalid internal degree of freedom type: " 
              << r_type << std::endl;
    inStream.setstate(std::ios::failbit);
  }
  return inStream;
}

inline std::ostream &operator<<(std::ostream &outStream, MoleculeType const &type)
{
  switch(type)
  {
    case LINEAR_SYMMETRIC:
      outStream << "LINEAR_SYMMETRIC";
      break;
    case PLANAR_SYMMETRIC:
      outStream << "PLANAR_SYMMETRIC";
      break;
    case LINEAR_ASYMMETRIC:
      outStream << "LINEAR_ASYMMETRIC";
      break;
    default:
      outStream << "GENERAL";
  }
  return outStream;
}

inline std::istream &operator>>(std::istream &inStream, MoleculeType &type)
{
  std::string r_type; // Type read from stream
  inStream >> r_type;
  if(r_type.find("LINEAR_SYMMETRIC") != r_type.npos)
    type = LINEAR_SYMMETRIC;
  else if(r_type.find("PLANAR_SYMMETRIC") != r_type.npos)
    type = PLANAR_SYMMETRIC;
  else if(r_type.find("LINEAR_ASYMMETRIC") != r_type.npos)
    type = LINEAR_ASYMMETRIC;
  else if(r_type.find("GENERAL") != r_type.npos)
    type = GENERAL;
  else
  {
    std::cerr << "I/O error: Invalid molecule type: " 
              << r_type << std::endl;
    inStream.setstate(std::ios::failbit);
  }
  return inStream;
}

// Mapping of internal degrees of freedom

struct InternalDOFMap
{
  InternalDOFType type;       // Type of internal DOF
  std::vector<size_t> atoms;  // Atoms defining the internal DOF
  size_t symmetryNumber;      // Symmetry number (only meaningful for dihedrals)

  InternalDOFMap()
  :
  type(DISTANCE), atoms(2, 0), symmetryNumber(1)
  {}
};

// To allow for linear combinations of atomic coordinates
// in defining the molecule map

struct FramePoint
{
  std::vector<size_t> atoms;  // The atoms whose coordinates make up the point
  std::vector<Real> weights;  // The weights on those coordinates

  size_t const numAtoms() const
  {
    return atoms.size();
  }

  void clear()
  {
    atoms.clear(); weights.clear();
  }
};

// Molecule map

class MoleculeMap
{
public:

  std::string name;                         // Name/residue name
  std::vector<FramePoint> framePoints;      // Points defining the molecule-centered frame
  MoleculeType type;                        // Molecule type
  std::vector<InternalDOFMap> internalDOFs; // Definitions of the internal degrees of freedom

// Constructors

  MoleculeMap(); // Defines an empty molecule map

// Print molecule map information (mostly for debugging)

  friend std::ostream &operator<<(std::ostream &outStream, MoleculeMap const &moleculeMap); 
};

/*
** End of class MoleculeMap
*/

// Inlines

inline MoleculeMap::MoleculeMap()
:
name(""), framePoints(), type(GENERAL), internalDOFs()
{}

inline std::ostream &operator<<(std::ostream &outStream, MoleculeMap const &moleculeMap)
{
  outStream << std::endl << "Name: " << moleculeMap.name;
  outStream << std::endl << "Type: " << moleculeMap.type;
  outStream << std::endl << "Frame points: ";
  for(size_t i = 0; i < moleculeMap.framePoints.size(); i++)
  {
    outStream << std::endl;
    for(size_t j = 0; j < moleculeMap.framePoints[i].numAtoms(); ++j)
      outStream << moleculeMap.framePoints[i].weights[j] << "*" 
                << moleculeMap.framePoints[i].atoms[j] << " ";
  }
  
  outStream << std::endl << "Internal degrees of freedom: " << moleculeMap.internalDOFs.size();
  for(size_t i = 0; i < moleculeMap.internalDOFs.size(); i++)
  {
    outStream << std::endl << "Type: " << moleculeMap.internalDOFs[i].type;
    outStream << std::endl << "Atoms: ";
    for(size_t j = 0; j < moleculeMap.internalDOFs[i].atoms.size(); j++)
      outStream << moleculeMap.internalDOFs[i].atoms[j] << " ";
    if(moleculeMap.internalDOFs[i].type == DIHEDRAL)
      outStream << std::endl << "Symmetry number: " << moleculeMap.internalDOFs[i].symmetryNumber;
  }
  return outStream;
}

#endif
