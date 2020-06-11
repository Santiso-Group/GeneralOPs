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
** .ptm file format
*/

#include "common/include/assert.h"
#include "crystdist/include/ptmfile.h"

// Constructors

PTMFile::PTMFile()
:
IOFile()
{}

PTMFile::PTMFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}

// Operators

PTMFile &operator<<(PTMFile &ptmFile, System<PointMolecule> const &system)
{
  // Write system to file
  assert(ptmFile.is_open() && ptmFile.mode() == OUT);
  ptmFile.setf(std::ios::fixed, std::ios::floatfield);
  ptmFile.setf(std::ios::right, std::ios::adjustfield);
  std::streamsize precision = ptmFile.precision();
  ptmFile.precision(WRITE_PRECISION);

  // System label and name
  ptmFile << "!SYSTEM: " << std::endl;
  ptmFile << system.name() << std::endl;
  // Write unit lattice information
  ptmFile << "!LATTICE: " << std::endl;
  ptmFile << system.lattice().a() << " " 
          << system.lattice().b() << " " 
          << system.lattice().c() << std::endl;
  ptmFile << system.lattice().alpha() << " " 
          << system.lattice().beta() << " " 
          << system.lattice().gamma() << std::endl;
  ptmFile << ((system.lattice().type() == SUPERCELL)?"SUPERCELL":"CRYSTAL") << std::endl;
  ptmFile << "!MOLECULES: " << std::endl;
  // Write point molecule information
  ptmFile << system.size() << std::endl;
  for(size_t i = 0; i < system.size(); i++)
  {
    ptmFile << system[i].name << std::endl;
    ptmFile << system[i].type << std::endl;
    ptmFile << system[i].position << std::endl;
    ptmFile << system[i].orientation.w << " "
            << system[i].orientation.x << " "
            << system[i].orientation.y << " "
            << system[i].orientation.z << " " << std::endl;
    ptmFile << system[i].numInternalDOFs() << std::endl;
    for(size_t j = 0; j < system[i].numInternalDOFs(); j++)
    {
      ptmFile << system[i].internalDOF(j).type << " "
              << system[i].internalDOF(j).value << " ";
      if(system[i].internalDOF(j).type == DIHEDRAL)
        ptmFile << system[i].internalDOF(j).symmetryNumber << " ";
      ptmFile << std::endl;
    }
  }
  ptmFile.precision(precision); // Return to original value
  return ptmFile;
}

PTMFile &operator>>(PTMFile &ptmFile, System<PointMolecule> &system)
{
  // Read system from file
  assert(ptmFile.is_open() && ptmFile.mode() == IN);
  system.clear();

  std::string r_label;  // Label read from file

  // Search for system label
  do
  {
    std::getline(ptmFile, r_label);
    if(r_label.find("!SYSTEM:") != r_label.npos) break;
  }
  while(!ptmFile.eof());
  if(ptmFile.eof())
  {
    std::cerr << "Error reading file " << ptmFile.fileName() << ": SYSTEM record not found." << std::endl;
    return ptmFile;
  }

  // Read system name;
  std::string r_name;   // Name read from file
  std::getline(ptmFile, r_name);
  //ptmFile >> r_name;
  //ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  if(!ptmFile)
  {
    std::cerr << "Error reading file " << ptmFile.fileName()
              << ": Invalid system name." << std::endl;
    system.clear();
    return ptmFile;
  }
  if(r_name.find("!LATTICE:") != r_name.npos)  // In case system name is empty
    system.setName("");
  else 
    system.setName(r_name);

  // Read unit lattice information
  Real r_a, r_b, r_c, r_alpha, r_beta, r_gamma; // Unit cell parameters read from file
  std::string r_latticeType;                    // Lattice type read from file

  if(r_name.find("!LATTICE:") != r_name.npos) // In case system name is empty
    r_label = r_name;
  else
  {
    ptmFile >> r_label;
    ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  }
  if(r_label.find("!LATTICE:") == r_label.npos)
  {
    std::cerr << "Error reading file " << ptmFile.fileName()
              << ": Invalid file format." << std::endl;
    system.clear();
    return ptmFile;
  }
  ptmFile >> r_a >> r_b >> r_c;
  ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  ptmFile >> r_alpha >> r_beta >> r_gamma;
  ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  ptmFile >> r_latticeType;
  if(!ptmFile)
  {
    std::cerr << "Error reading file " << ptmFile.fileName()
              << ": Invalid unit cell data." << std::endl;
    system.clear();
    return ptmFile;
  }
  if(r_latticeType.find("CRYSTAL") != r_latticeType.npos)
    system.setLattice(Lattice(r_a, r_b, r_c, r_alpha, r_beta, r_gamma, CRYSTAL));
  else if (r_latticeType.find("SUPERCELL") != r_latticeType.npos)
    system.setLattice(Lattice(r_a, r_b, r_c, r_alpha, r_beta, r_gamma, SUPERCELL));
  else
  {
    std::cerr << "Error reading file " << ptmFile.fileName()
              << ": Invalid unit cell type." << std::endl;
    system.clear();
    return ptmFile;
  }
  
  // Read molecule information
  size_t r_nMol;        // Number of molecules read from file
  ptmFile >> r_label;
  ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  if(r_label.find("!MOLECULES:") == r_label.npos)
  {
    std::cerr << "Error reading file " << ptmFile.fileName()
              << ": Invalid file format." << std::endl;
    system.clear();
    return ptmFile;
  }
  ptmFile >> r_nMol;
  ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
  if(!ptmFile)
  {
    std::cerr << "Error reading file " << ptmFile.fileName()
              << ": Invalid number of molecules." << std::endl;
    system.clear();
    return ptmFile;
  }
  for(size_t i = 0; i < r_nMol; i++)
  {
    PointMolecule currentPointMolecule; // Current point molecule
    MoleculeType r_type;                // Molecule type read from file
    Vector3D r_pos;                     // Position read from file
    Quaternion r_quat;                  // Orientation read from file
    size_t r_numIntDOFs;                // Number of internal DOFs read from file

    currentPointMolecule.clear();
    //std::getline(ptmFile, r_name); // Avoid issues with whitespace in residue names
    ptmFile >> r_name;
    ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    currentPointMolecule.name = r_name;
    ptmFile >> r_type;
    ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    currentPointMolecule.type = r_type;
    ptmFile >> r_pos.x >> r_pos.y >> r_pos.z;
    ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    currentPointMolecule.position = r_pos;
    ptmFile >> r_quat.w >> r_quat.x >> r_quat.y >> r_quat.z;
    r_quat.normalize();
    ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    currentPointMolecule.orientation = r_quat;
    ptmFile >> r_numIntDOFs;
    ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    if(!ptmFile)
    {
      std::cerr << "Error reading file " << ptmFile.fileName()
                << ": Invalid molecule data" << std::endl;
      system.clear();
      return ptmFile;
    }
    for(size_t j = 0; j < r_numIntDOFs; j++)
    {
      InternalDOF currentDOF;
      ptmFile >> currentDOF.type >> currentDOF.value;
      if(currentDOF.type == DIHEDRAL)
        ptmFile >> currentDOF.symmetryNumber;
      if(!ptmFile)
      {
        std::cerr << "Error reading file " << ptmFile.fileName() << std::endl
                  << ": Invalid internal degree of freedom record." << std::endl;
        system.clear();
        return ptmFile;
      }
      currentPointMolecule.internalDOFs.push_back(currentDOF);
      ptmFile.ignore(LINE_LENGTH, '\n');  // Skip to next line
    }
    system.add(currentPointMolecule);
  }
  return ptmFile;
}
