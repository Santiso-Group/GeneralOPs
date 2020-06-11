/*
** Copyright 2008-2011 Erik Santiso.
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

/********************************************************************/
/*                                                                  */
/*  mollib                                                          */
/*                                                                  */
/*  Simple molecular modeling library                               */
/*                                                                  */
/*                                                                  */
/********************************************************************/

// LOOK AT THE TO-DO'S BELOW AND UPDATE TODO LIST...

// Some notes:
//
// All indices are internally defined to count from zero.
//
// Things that may be good:
//
// - ****** Add a method to reserve() the atom and bond arrays in geometry
//   and system classes - makes reading from some file formats faster *******
//
// - Add a method to geometry to update bonds after moving a single atom?
//   (more efficient if moving only one or a few atoms)
//
// - In xyzfile - add a variable to toggle automatic determination of
//   bonds. Otherwise it can take very long to load large files
//
// - Return to const and non-const versions of the atom method in 
//   geometry? Maybe more elegant than setx, sety, setz?
//
// ******************************
// ******************************
// ADD A FUNCTION TO WRAP A MOLECULE AROUND A PERIODIC BOX!!!!
// ******************************
// ******************************
// Things that would be good to add:
//
// - Wrap everything in a namespace
//
// - Eventually, replace most asserts and outputs to std::cerr by proper exception handling
//
// - Read pdb and dcd files - DCD and PSF are done, BUT PSF ONLY READS 1-LETTER SYMBOLS
// - PDB ALSO DONE, BUT ALSO ONLY READS 1-LETTER SYMBOLS
//
// - Add atomgraph/system class (with reverse tables etc.)? - check comment at
// the top of geometry.h
// ***** REGARDING THE ABOVE - MAYBE WHAT IS JUST NEEDED IS AN UTILITY TO CONVERT GEOMETRY
// TO SYSTEM<GEOMETRY> OBJECTS, I.E. TO FIND CONNECTIVITY. THIS WOULD ALSO BE USEFUL
// WHEN READING XYZ FILES, ETC... ALSO UTILITY FOR VICE-VERSA (GEOMETRY FROM SYSTEM)
// -> THAT ONE IS TRIVIAL
// ***** ALSO UTILITY TO FIND ANGLES/DIHEDRALS/ETC. - THIS WOULD BE GOOD TO WRITE PSF'S
// ***** THE UTILITY SHOULD BE ABLE TO FIND THE CONNECTIVITY EITHER FROM DISTANCES
// OR FROM RESIDUE IDS... <==== NOT REALLY NECESSARY FROM RESID's (USE GEOMETRY or SYSTEM
// DEPENDING ON NEEDS)
// ***** THIS WOULD ALLOW READING FROM XYZ FILES TO SYSTEM<*> OBJECTS, ETC...
// ***** ADD A PROPERTY TO TURN ON/OFF THE AUTOMATIC DETERMINATION OF CONNECTIVITY
// - FOR LARGE SYSTEMS IT IS CRAZY SLOW...
//
//************
    // NOTES MAYBE FOR NEXT VERSION - READ THE MASS FROM THE PSF.
    // ALSO, MAKE THE ATOM LIBRARY A STATIC THAT IS READ ONLY ONCE, AS IT IS DONE
    // IN THE CHROMOSOME CLASS (FORM2GEOM) FOR THE FORBIDDEN LIST. MUCH FASTER...
// NEED SOME RELIABLE WAY TO TELL WHETHER A DCD FILE HAS REACHED THE END
// - RIGHT NOW THE CHECKS FOR FRAMESREAD<NUMFRAMES ARE TURNED OFF
  // ******************* NEED TO MODIFY ALL WRITING ROUTINES TO GIVE AN ERROR WHEN WRITING AN EMPTY GEOMETRY/SYSTEM!!!!!

//************
// - Non-member functions to calculate things for geometry objects?
// (actually maybe it's better to leave these for a more specific class,
// like molecule/atomgraph...).
// (1) Mass, center of mass, moment of inertia, symmetry, etc.
// (2) Connectivity (atoms bonded to, clusters)
// (3) Rigid rotation and translation (requires setPosition method) -
// ideally this would allow also to align bonds with axes, angles with planes, etc.
// (4) Change distance, angle, dihedral (maybe)
//
// - Determination of bond orders (see list of mol2 bond types below). This
// requires adding a "valence" property to atoms.
//
// - OpenGL functions to draw geometry objects.
//
// - For the crystal project: - BOTH DONE
// (1) A class "simple molecule" (or similar) containing a center of mass and
// and absolute orientation (vector or quaternion). This should allow as input a
// molecule, the atoms that define the molecule-centered frame, and the information 
// (symmetry) necessary to decide how to represent absolute orientations.
// (2) Then a container for "simple molecules" that can calculate all
// distances, bond orientations, and relative orientations, and (non-member?) 
// functions to do statistics on them, and calculate order parameters and derivatives.
//
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
// TO DO:
// CHANGING THE LATTICE TYPE IS A BIT AWKWARD RIGHT NOW - NEED TO CREATE A
// COPY OF THE LATTICE, CHANGE ITS TYPE, AND USE SETLATTICE(). PLUS, WHEN
// READING FROM FILES CONTAINING UNIT CELL INFO (E.G. MOL2, DCD), TYPE
// IS RESET TO SUPERCELL. IT'D BE GOOD TO:
// 1. ADD A METHOD TO JUST SET THE TYPE OF THE LATTICE IN CLASS SYSTEM<> - DONE
// 2. MAKE DCDFILE, MOL2FILE, ETC. NOT CHANGE THE LATTICE TYPE WHEN
//    READING... - DONE
// BOTH OF THESE CHANGES ARE SMALL AND EASY.
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
/* 

Mol2 Bond types are:

1 = single
2 = double
3 = triple
am = amide (peptide bond)
ar = aromatic
du = dummy
un = unknown
nc = not connected

*/

/* 
** ADD THE GNU COPYRIGHT INFO ON TOP OF EVERYTHING
** (CHECK LINK ON DESKTOP)
**/

#include "stdlib.h"
#include "time.h"
#include "../../common/include/iofile.h"
#include "../include/discretefield.h"
#include "../include/grid.h"
#include "../include/file_formats/fileformats.h"
#include "../include/molecule.h"
#include "../include/pairpotential.h"

using namespace std;


int main(int argc, char* argv[])
{

// TEST
/*
  PDBFile myOtherPDBFile("my.pdbCOMB1.101", IN);
  System<Molecule> myOtherSystem;
  myOtherPDBFile >> myOtherSystem;
  PDBFile myOther1File("pdb_1.pdb", OUT);
  myOther1File << myOtherSystem;
  myOtherSystem.reset();
  myOtherPDBFile >> myOtherSystem;
  PDBFile myOther2File("pdb_2.pdb", OUT);
  myOther2File << myOtherSystem;
  return 0;
*/

  // TPA dcd test
/*
  //PSFFile myTPAPSF("f1_512.xplor.psf");
  //DCDFile myTPADCD1("terephthalic_cry.dcd");
  PSFFile myTPAPSF("glycine.psf");
  DCDFile myTPADCD1("glycine_trimmed.dcd");
  //DCDFile myTPADCD2("form1_trimmed.dcd");
  System<Molecule> myTPASystem1;
  //System<Molecule> myTPASystem2;
  myTPAPSF >> myTPASystem1;
  //myTPAPSF >> myTPASystem2;
  myTPADCD1 >> myTPASystem1;
  //myTPADCD2 >> myTPASystem2;
  //std::cout << myTPASystem1.lattice() << std::endl;
  std::cout << myTPASystem1.lattice() << std::endl;
  PDBFile myPDBFileTPA("gly1.pdb", OUT);
  myPDBFileTPA << myTPASystem1 << std::endl;
  return 0;


  // Lattice test
  Real a = 78.825;
  Real b = 50.210;
  Real c = 30.809;
  Real alpha = 94.28;
  Real beta = 107.47;
  Real gamma = 51.13;

  Lattice myLattice1(a, b, c, alpha, beta, gamma);
  std::cout << myLattice1 << std::endl;
*/
/*
  std::cout << "a = " << myLattice1.a() << std::endl;
  std::cout << "b = " << myLattice1.b() << std::endl;
  std::cout << "c = " << myLattice1.c() << std::endl;
  std::cout << "alpha = " << myLattice1.alpha() << std::endl;
  std::cout << "beta = " << myLattice1.beta() << std::endl;
  std::cout << "gamma = " << myLattice1.gamma() << std::endl;
*/  
//  return 0;

/*
  // TEST - rotate TPA
  PDBFile TPAPDB("form1.pdb", IN);
  System<Molecule> TPASystem;
  TPAPDB >> TPASystem;

  Vector3D currentv(7.062, 6.408, 0.0); // a lattice vector
  Vector3D finalv(1.0, 0.0, 0.0); // x-axis
  Vector3D rotaxis = cross(currentv, finalv);
  Real cosangle = currentv*finalv/(norm(currentv)*norm(finalv));
  if(cosangle > 1.0) cosangle = 1.0;
  if(cosangle < -1.0) cosangle = -1.0;
  Real angle = acos(cosangle);
  Quaternion myQuat(rotaxis, angle);
  Vector3D origin(20.7472, 18.1435, 12.3483);
  for(size_t i = 0; i < TPASystem.size(); ++i)
    TPASystem[i].rotate(myQuat, origin);
  Vector3D currenty(3.811, 5.194, 0.0); // rotated lattice vector
  Vector3D currentxy = cross(finalv, currenty);
  Vector3D finaln(0.0, 0.0, 1.0); // z-axis
  rotaxis = cross(currentxy, finaln);
  cosangle = currentxy*finaln/norm(currentxy);
  if(cosangle > 1.0) cosangle = 1.0;
  if(cosangle < -1.0) cosangle = -1.0;
  cosangle =0.0 ;
  angle = acos(cosangle);
  rotaxis = finalv;
  Quaternion myQuat2(rotaxis, angle);
  //for(size_t i = 0; i < TPASystem.size(); ++i)
    //TPASystem[i].rotate(myQuat2, origin);
  PDBFile TPArPDB("form1_r.pdb", OUT);
  TPArPDB << TPASystem;

  
  return 0;
*/
  // TEST - Ignoring TER records
  //System<Molecule> thisIsMyGeometry;
  //PDBFile thisIsMyPDBFile("mypdb40wat.350.pdb", IN);
  //PDBFile thisIsAnother("test350.pdb", OUT);
  //thisIsMyPDBFile >> thisIsMyGeometry;
  //thisIsAnother << thisIsMyGeometry;
  //return 0;

  // TEST - PDB Temperature factors
  // Do tests below (see where it says "Done")
  // See if temp factors work, then write
  // a code to write the cop's to the
  // temp. factors...

  // TO DO - RUN THESE TESTS!!! -> DONE

  // Tests to do:
  //
  // xyz to everything
  // mol2 to everything
  // pdb alone to everything
  // psf + pdb to everything
  // psf + dcd to everything
  // Including checking that wrong psf gives error!
  //
  // where everything = xyz, mol2, pdb
  // for both Geometry, System<Geometry> and System<Molecule> objects
  // Also check writing to *same* type! => DONE
  // And reading/writing files with different residue types...
  // Also CONFIG to everything and back, and CONFIG+FIELD to everything
  // Here it goes:

  std::string base_path("/Users/erik/Projects/mymol/debug_files/");

  Geometry myGeometry;
  System<Geometry> myGeometries;
  System<Molecule> myMolecules;

  // Works for all:
  //XYZFile myXYZInputFile(base_path+"benzene_crystal.xyz", IN);
  //myXYZInputFile >> myGeometry;

  //PSFFile myPSFFile(base_path+"benzene_cry.psf");
  //PSFFile myPSFFile(base_path+"glycine.psf");
  //myPSFFile >> myGeometry;
  //myPSFFile >> myMolecules;
  //CONFIGFile myCONFIGInputFile(base_path+"CONFIG", IN);
  //CONFIGFile myCONFIGInputFile(base_path+"REVCON", IN);
  //myCONFIGInputFile >> myGeometries;
  //myCONFIGInputFile >> myMolecules;

  //FIELDFile myFIELDFile(base_path+"TEST"); 
  //FIELDFile myFIELDFile(base_path+"FIELD");
  //FIELDFile myFIELDFile(base_path+"FIELD_1W");
  //myFIELDFile >> myGeometry;
  //myFIELDFile >> myGeometries;
  //myFIELDFile >> myMolecules;

  // TEST
  //for(size_t i = 0; i < myMolecules[21000].numBonds(); ++i)
  //  std::cout << myMolecules[21000].bond(i) << std::endl;
  //return 0;

  // TEST - writing FIELD files
  //XYZFile myXYZFile(base_path+"surfactant.xyz", IN);
  /*
  XYZFile myXYZFile(base_path+"one_PEG.xyz", IN);
  XYZFile myOtherXYZFile(base_path+"one_water.xyz", IN);
  FIELDFile myFIELDFile(base_path+"FIELD_TEST", OUT);
  myFIELDFile.setUnitLabel("kelvin");
  myGeometries.add(Geometry());
  myGeometries.add(Geometry());
  myGeometries.add(Geometry());
  myXYZFile >> myGeometry;
  myOtherXYZFile >> myGeometries[2];
  */
  // Fix bonds
/*
  for(size_t i = 0; i < myGeometry.numBonds(); ++i)
  {
    Bond const &bond = myGeometry.bond(i);
    myGeometry.setBondType(i, "harm");
    std::vector<Real> parameters(1, 937.5);
    Real length = myGeometry.atom(bond.first).radius +
                  myGeometry.atom(bond.second).radius;
    parameters.push_back(length);
    myGeometry.setBondParameters(i, parameters);
  }
  */
  /*
  // TEST - bonds and constraints
  for(size_t i = 0; i < 2; ++i)
  {
    myGeometry.setBondType(i, "12-6");
    std::vector<Real> parameters;
    parameters.push_back(100.0);
    parameters.push_back(50.0);
    myGeometry.setBondParameters(i, parameters);
  }
  for(size_t i = 2; i < 4; ++i)
  {
    myGeometry.setBondType(i, "cons");
    std::vector<Real> parameters(1, 1.0);
    myGeometry.setBondParameters(i, parameters);
  }
  for(size_t i = 4; i < myGeometry.numBonds(); ++i)
  {
    myGeometry.setBondType(i, "harm");
    std::vector<Real> parameters;
    parameters.push_back(100.0);
    parameters.push_back(50.0);
    myGeometry.setBondParameters(i, parameters);
  }
  */
  /*
  myGeometries[0] = myGeometry;
  myGeometries[1] = myGeometry;
  myFIELDFile << myGeometries;
  */

  // Pair potential tests:

  /*
  PairPotential myPotential("CM", "EO", "lj");
  std::vector<Real> parameters;
  parameters.push_back(200.0);
  parameters.push_back(1.5);
  myPotential.parameters = parameters;
  std::cout << myPotential << std::endl;
  return(0);
  */

  /*
  std::vector<PairPotential> myPotentials;
  std::vector<Real> parameters(4, 0.0);
  parameters[0] = 461.11;
  parameters[1] = 19.00;
  parameters[2] = 6.00;
  parameters[3] = 4.0697;
  myPotentials.push_back(PairPotential("OA", "OA", "nm", parameters));
  parameters[0] = 318.67;
  parameters[1] = 16.86;
  parameters[2] = 6.00;
  parameters[3] = 4.4476;
  myPotentials.push_back(PairPotential("OA", "CM", "nm", parameters));
  parameters[0] = 344.42;
  parameters[1] = 15.00;
  parameters[2] = 6.00;
  parameters[3] = 4.8312;
  myPotentials.push_back(PairPotential("CM", "CM", "nm", parameters));

  FIELDFile myFIELDFile(base_path+"FIELD_PP", OUT);
  myFIELDFile << myPotentials;
  myFIELDFile.close();
  return 0;
  */

  /*
  std::vector<PairPotential> myPotentials;
  FIELDFile myFIELDInputFile(base_path+"FIELD", IN);
  FIELDFile myFIELDOutputFile(base_path+"FIELD_PP", OUT);
  myFIELDInputFile >> myPotentials;
  myFIELDOutputFile << myPotentials;
  myFIELDOutputFile.close();
  return 0;
  */

  // TEST
  //DCDFile myDCDFile;
  //myDCDFile.setFile("pupu", OUT);

  // Reading/writing entire field file

  FIELDFile myFIELDInputFile(base_path+"FIELD", IN);
  FIELDFile myFIELDOutputFile(base_path+"FIELD_2", OUT);
  myFIELDInputFile >> myGeometries;
  myFIELDOutputFile.setUnitsLabel(myFIELDInputFile.unitsLabel());
  std::vector<PairPotential> myPotentials;
  myFIELDInputFile >> myPotentials;
  myFIELDOutputFile << myGeometries;
  myFIELDOutputFile << myPotentials;
  myFIELDOutputFile.close();

  // TEST
  std::cout << "Potentials: " << std::endl;
  for(size_t i = 0; i < myPotentials.size(); ++i)
    std::cout << myPotentials[i] << std::endl;

  // Energy and force test (limited!)
  PairPotential const &potential = myPotentials[myPotentials.size() - 1];
  Vector3D const r1 = 0.0;
  Vector3D const r2 = 3.0;
  Vector3D const r12 = r2 - r1;
  Real const r = norm(r12);
  std::cout << "Energy: " << potential.energy(r) << std::endl;
  std::cout << "Force: " << potential.force(r12) << std::endl;

  Real const tiny = 0.0001;
  Vector3D const direction(1.0, -1.0, 2.0);
  Real const r12p = norm(r12 + unit(direction)*tiny);
  Real const r12m = norm(r12 - unit(direction)*tiny);
  Real const forceNum = -(potential.energy(r12p) - potential.energy(r12m))/
                      (tiny + tiny);
  std::cout << "Compare: " << std::endl;
  std::cout << potential.force(r12)*unit(direction) << std::endl;
  std::cout << forceNum << std::endl;

/*
  IOFile myPlot("plot", OUT);
  for(size_t i = 0; i < 1001; ++i)
  {
    Real const r_i = 1.0 + 10.0*(Real)i/1000.0;
    myPlot << r_i  << " " << potential.energy(r_i) << std::endl;
  }
  myPlot.close();
*/
  return 0;

  /*
  FIELDFile myFIELDInputFile(base_path+"FIELD", IN);
  FIELDFile myFIELDOutputFile(base_path+"FIELD_TEST", OUT);
  myFIELDOutputFile.setUnitsLabel("kelvin");
  myFIELDInputFile >> myMolecules;
  std::cout << "Units are: " << myFIELDInputFile.unitsLabel() << std::endl;
  myFIELDOutputFile << myMolecules;
  
  return 0;
  */
  //CONFIGFile myCONFIGFile(base_path+"CONFIG", IN);
  //myCONFIGFile >> myGeometry;
  //myCONFIGFile >> myGeometries;

  // Add residue info

  //HISTORYFile myHISTORYFile(base_path+"HISTORY");
  //myHISTORYFile >> myGeometry;
  //myHISTORYFile.skipFrames(3);
  //myHISTORYFile >> myGeometry;
  //myHISTORYFile >> myGeometries;
  //myHISTORYFile >> myGeometries;
  //myHISTORYFile >> myMolecules;

  // TEST - wrap all except first residue
  /*
  for(size_t i = 1; i < myGeometries.size(); ++i)
  {
    Geometry &current = myGeometries[i];
    Lattice const &lattice = myGeometries.lattice();
    Vector3D const pos = current.atom(0).position;
    for(size_t iAtom = 0; iAtom < current.numAtoms(); ++iAtom)
    {
      Vector3D diff = lattice.difference(pos, current.atom(iAtom).position);
      current.setPosition(iAtom, pos + diff);
    }
  }
  */
  //for(size_t i = 1; i < myGeometries.size(); ++i)
   // myGeometries.join(i);
  //myGeometries.joinAll();
  //std::cout << myGeometries[27436].atom(2) << std::endl;

  //PDBFile myPDBFile(base_path+"test_config.pdb", OUT);
  //myPDBFile << myMolecules;

  // Works for all:
  //MOL2File myMOL2InputFile("benzene_crystal.mol2", IN);

  // Works for all:
  //PDBFile myPDBInputFile(base_path+"frame_00.pdb", IN);
  //PDBFile myPDBInputFile(base_path+"test_pdb.pdb", IN);
  //myPDBInputFile >> myGeometry; 
  //std::cout << myGeometry.atom(0).tempFactor << std::endl;
  
  //PSFFile myPSFInputFile(base_path+"benzene_cry.psf");
  //myPSFInputFile >> myGeometries;
  //myPDBInputFile >> myGeometries;

  //CONFIGFile myCONFIGFile(base_path+"CONFIG_out", OUT);
  //myCONFIGFile << myGeometries;
  //myCONFIGFile << myGeometry;
 
  //DCDFile myDCDInputFile(base_path+"benzene_cry.dcd");
  //PDBFile myPDBInputFile("solv_x_min.pdb", IN);
  //PDBFile myPDBInputFile("test_xyz.pdb", IN);

  // Works for all:
  //myMOL2InputFile >> myGeometry;

  // Works for all:
  //myMOL2InputFile >> myGeometries;

  //myPSFInputFile >> myGeometries;
  //myDCDInputFile >> myGeometry;

  // Works:
  //XYZFile myXYZFile(base_path+"test_xyz.xyz", OUT);
  //myXYZFile << myGeometry;

  // Works:
  //MOL2File myMOL2File(base_path+"test_xyz.mol2", OUT);
  //myMOL2File << myGeometry;
  //myMOL2File << myGeometries;
  //myMOL2File << myMolecules;

  // Works:
  
  //for(size_t i = 0; i < myGeometry.numAtoms(); ++i)
  //  myGeometry.setTemperatureFactor(i, i/10.0);
  //PDBFile myPDBFile(base_path+"test_pdb.pdb", OUT);
  //PDBFile myPDBFile(base_path+"test_pdb2.pdb", OUT);
  //myPDBFile << myGeometry;
  //myPDBFile << myGeometries;
  //myPDBFile << myMolecules;
  
  // OK, THIS WORKS...
  
  // End of file type tests

  // Tests for packmol
  /*
  CONFIGFile myCONFIGFile(base_path+"CONFIG_one_PEG", IN);
  FIELDFile myFIELDFile(base_path+"FIELD_one_PEG");
  PDBFile myPDBFile(base_path+"one_PEG.pdb", OUT);
  myFIELDFile >> myGeometries;
  myCONFIGFile >> myGeometries;
  myPDBFile << myGeometries;
  */
  // TEST
  //std::cout << myGeometries[0].atom(1) <<std::endl; 
  /*
  Geometry one_water;
  one_water.addAtom(Atom("W", 0.0));
  PDBFile myPDBFile(base_path+"one_W.pdb", OUT);
  myPDBFile << one_water;
  */
 
  /* 
  FIELDFile myFIELDFile(base_path+"FIELD_system_box");
  //PDBFile myPDBFile(base_path+"system_box.pdb", IN);
  XYZFile myXYZFile(base_path+"system_box.xyz", IN);
  CONFIGFile myCONFIGFile(base_path+"CONFIG_system", OUT);
  myFIELDFile >> myGeometry;
  myXYZFile >> myGeometry;
  System<Geometry> oneGeometry;
  oneGeometry.add(myGeometry);
  Lattice pbcLattice(85.015, 85.015, 341.0);
  oneGeometry.setLattice(pbcLattice);
  myCONFIGFile << oneGeometry;
  */
  /*
  FIELDFile myFIELDFile(base_path+"FIELD_water_box");
  //PDBFile myPDBFile(base_path+"system_box.pdb", IN);
  XYZFile myXYZFile(base_path+"water_box.xyz", IN);
  CONFIGFile myCONFIGFile(base_path+"CONFIG_water", OUT);
  myFIELDFile >> myGeometry;
  myXYZFile >> myGeometry;
  System<Geometry> oneGeometry;
  oneGeometry.add(myGeometry);
  Lattice pbcLattice(85.015, 85.015, 341.0);
  oneGeometry.setLattice(pbcLattice);
  myCONFIGFile << oneGeometry;

  exit(0);
  */

  // Grid/DiscreteField tests


  // Big field test

  //srand(time(NULL));
  srand(18593727);

  Lattice myLattice(1.0, 1.0, 1.0, 90, 90, 90);
  Vector3D origin(0.5, 0.5, 0.5);
  Grid myGrid(5, myLattice, origin);


  // TEST

  size_t const indx = 16;
  std::vector<size_t> const adj = myGrid.adjacentCells(indx);
  std::cout << "Adjacent indices: ";
  for(size_t i = 0; i < adj.size(); ++i)
    std::cout << adj[i] << " ";
  std::cout << std::endl;
  exit(0);
  // END TEST

  DiscreteField<Real> myField(myGrid, 0.0, "The data");

  size_t numPoints = 10000;
  for(size_t i = 0; i < numPoints; ++i)
  {
    Real x = (Real)(rand()%1001)/1000.0;
    Real y = (Real)(rand()%1001)/1000.0;
    Real z = (Real)(rand()%1001)/1000.0;
    Real val = x*x + y*y + z*z;
    myField.addValue(Vector3D(x, y, z), val); 
  }

  std::cout << myField << std::endl;

  VTKFile myVTKFile("field.vtk");
  myVTKFile << myField;
  exit(0);

  // end of big field test
/*
  Lattice myLattice(2.0, 1.0, 1.0, 90, 90, 45);
  Vector3D origin = 0.0;
//  Vector3D origin(1.0, 0.5, 0.5);
//  Grid myGrid(5, myLattice, origin);
//  Grid myGrid(10, 4, 2, myLattice);
  Grid myGrid(4, myLattice, origin);

  std::cout << myGrid << std::endl;

  VTKFile myVTKFile("grid.vtk");
  myVTKFile << myGrid << std::endl;

  Vector3D testV(0.1, 0.99, 0.6);
  std::cout << "Vector: " << testV << std::endl;
  std::cout << "Cell index: " << myGrid.cellIndex(testV) << std::endl;

  DiscreteField<Real> myField(myGrid, 0.05);
  myField.addValue(testV, 0.0);

  exit(0);

  size_t indx = 79;
  Vector3D center = myGrid.cellCenter(indx);

  std::cout << "Index: " << indx << std::endl;
  std::cout << "Cell center: " << center << std::endl;
  std::cout << "Index of cell: " << myGrid.cellIndex(center) << std::endl;

  Lattice myEmptyLattice;
  myGrid.setLattice(myEmptyLattice);

  std::cout << myGrid << std::endl;
  std::cout << "Vector: " << testV << std::endl;
  std::cout << "Cell index: " << myGrid.cellIndex(testV) << std::endl;

  myGrid.setLattice(myLattice);
  myGrid.setNumDivisions(4, 2, 1);

  std::cout << myGrid << std::endl;
  std::cout << "Vector: " << testV << std::endl;
  std::cout << "Cell index: " << myGrid.cellIndex(testV) << std::endl;

  Lattice myFlatLattice(1.0, 2.0, 0.0, 90, 90, 45);
  myGrid.setLattice(myFlatLattice);

  std::cout << myGrid << std::endl;
  std::cout << "Vector: " << testV << std::endl;
  std::cout << "Cell index: " << myGrid.cellIndex(testV) << std::endl;

  Lattice myLineLattice(0.0, 0.0, 3.0);
  myGrid.setLattice(myLineLattice);

  std::cout << myGrid << std::endl;
  std::cout << "Vector: " << testV << std::endl;
  std::cout << "Cell index: " << myGrid.cellIndex(testV) << std::endl;

  //DiscreteField<Real> myField(myGrid, 0.1);
  //myField.addValue(testV, 0.0);

  exit(0);



  myGrid.clear();

  std::cout << myGrid << std::endl;
  std::cout << "Vector: " << testV << std::endl;
  std::cout << "Cell index: " << myGrid.cellIndex(testV) << std::endl;

  exit(0);
*/
// TEST - SKIP TO END OF DCD
//
//  PSFFile myPSFFile("benzene_cry.psf");
//  DCDFile myDCDFile("benzene_cry.dcd");
//  std::cout << myDCDFile.numFrames() << std::endl;
//  System<Molecule> mySystem;
//  myPSFFile >> mySystem;
//  myDCDFile.skipToLastFrame();
//  std::cout << myDCDFile.numFramesRead() << std::endl;
////  myDCDFile.skipFrames(300);
////  std::cout << myDCDFile.numFramesRead() << std::endl;
//  myDCDFile >> mySystem;
//  PDBFile myOutPDBFile("benzene_cry.pdb", OUT);
//  myOutPDBFile << mySystem;
//  return 0;

  // TEST (PDB)
  //PSFFile myPSFFile("benzene_cry.psf");
  //PDBFile myPDBFile("test.pdb");
  //PDBFile myPDBFile("p-didehydrobenzene_noCONECT.pdb");
  //Geometry myGeometry;
  //myPSFFile >> myGeometry;
  //myPDBFile >> myGeometry;
  //myGeometry.calculateBonds();
  //MOL2File myMOL2File("test3.mol2", OUT);
  //myMOL2File << myGeometry;
  //PSFFile myPSFFile("benzene_cry.psf");
  //PDBFile myPDBFile("test.pdb");
  //System<Geometry> mySystem;
  //MOL2File myMOL2FileTest("test5.mol2",IN);
  //myMOL2FileTest >> mySystem;
  //PDBFile myPDBFile("p-didehydrobenzene_noCONECT.pdb", IN);
  //myPSFFile >> mySystem;
  //myPDBFile >> mySystem;
  //MOL2File myMOL2File("test8.mol2", OUT);
  //myMOL2File << mySystem;
  //Geometry myGeometry;
  //PSFFile myPSFFile("benzene_cry.psf");
  //myPSFFile >> myGeometry;
  //PDBFile myOtherPDBFile("test.pdb", IN);
  //myOtherPDBFile >> myGeometry;
  //return 0;
  //MOL2File myMOL2File("test8.mol2", IN);
  //myMOL2File >> myGeometry;
  //PDBFile myPDBFile("testpdb2.pdb", OUT);
  //myPDBFile << myGeometry;
  //std::cout << myGeometry.atom(0).type << std::endl;
/*
  // Tests to do:
  //
  // xyz to everything        => DONE
  // mol2 to everything       => DONE
  // pdb alone to everything  => DONE
  // psf + pdb to everything  => DONE
  // psf + dcd to everything  => DONE
  // Including checking that wrong psf gives error! => DONE
  //
  // where everything = xyz, mol2, pdb
  // for both Geometry, System<Geometry> and System<Molecule> objects
  // Also check writing to *same* type! => DONE
  // And reading/writing files with different residue types...
  // Here it goes:
  
  Geometry myGeometry;
  System<Geometry> myGeometries;
  System<Molecule> myMolecules;

  // Works for all:
  //XYZFile myXYZInputFile("benzene_crystal.xyz", IN);
  //myXYZInputFile >> myGeometry;

  // Works for all:
  //MOL2File myMOL2InputFile("benzene_crystal.mol2", IN);

  // Works for all:
  //PDBFile myPDBInputFile("frame_00.pdb", IN);

  PSFFile myPSFInputFile("benzene_cry.psf");

  DCDFile myDCDInputFile("benzene_cry.dcd");
  //PDBFile myPDBInputFile("solv_x_min.pdb", IN);
  //PDBFile myPDBInputFile("test_xyz.pdb", IN);

  // Works for all:
  //myMOL2InputFile >> myGeometry;

  // Works for all:
  //myMOL2InputFile >> myGeometries;

  myPSFInputFile >> myGeometries;
  myDCDInputFile >> myGeometries;

  // Works:
  //XYZFile myXYZFile("test_xyz.xyz", OUT);
  //myXYZFile << myGeometry;

  // Works:
  MOL2File myMOL2File("test_xyz.mol2", OUT);
  myMOL2File << myGeometries;

  // Works:
  PDBFile myPDBFile("test_pdb.pdb", OUT);
  myPDBFile << myGeometries;


  // End of file type tests

  return 0;
*/
// Test

  //Atom myAtom("P", Vector3D(1.0, 2.0, 1.0));
  //cout << myAtom << endl;

  //XYZFile myFile("GEOMETRY.xyz", IN);
  //Geometry myGeometry;
  //myFile >> myGeometry;

  //myGeometry.setLattice(Lattice(1.0, 2.0, 1.0, 85.0, 120.0, 90.0));
  //cout << myGeometry.lattice() << endl;
  //MOL2File moFile("pupu.mol2", OUT);
  //moFile << myGeometry;

  //MOL2File pu1File("GEOMETRY.mol2", IN);
  //pu1File >> myGeometry;
  //myGeometry.setLattice(Lattice(1.0, 2.0, 1.0, 85.0, 120.0, 90.0));
  //MOL2File pu2File("pupu2.mol2", OUT);
  //pu2File << myGeometry;
  //for(size_t i = 0; i < myGeometry.numBonds(); i++) { cout << myGeometry.bond(i) << endl; }
  //MOL2File moFile("pupu.mol2", OUT);
  //moFile << myGeometry;
  //XYZFile pu1File("pupu1.xyz", OUT);
  //XYZFile pu2File("pupu2.xyz", OUT);
  //myGeometry.addAtom(myAtom);
  //pu1File << myGeometry;
  //myGeometry.deleteAtom(2);
  //pu2File << myGeometry;
  //MOL2File mo2File("GEOMETRY.mol2", IN);
  //mo2File >> myGeometry;
  //MOL2File mo3File("pu.mol2", OUT);
  //mo3File << myGeometry;

  //Vector3D a(1.0, 0.0, 0.0);
  //Vector3D b(15.0, 1.0, 0.0);
  //Vector3D v1(8.0, 0.5, 0.0);
  //Vector3D v2(14.9, 0.9, 0.0);

  //Lattice lat(a, b);
  //cout << lat.difference(v1, v2) << endl;

  // TEST
/*
  XYZFile myFile("debug_files/formaldehyde.xyz", IN);
  Molecule myMolecule;
  myFile >> myMolecule;
*/
/*
  cout << "Mass: " << myMolecule.mass() << endl;
  cout << "Center: " << myMolecule.center() << endl;
  cout << "Center of Mass: " << myMolecule.centerOfMass() << endl;
  cout << "Moments of Inertia: " << myMolecule.momentsOfInertia() << endl;
  cout << "Principal Axes: " << myMolecule.principalAxes() << endl;
 
  // New: gyration tensor et al.

  cout << "Radius of Gyration: " << myMolecule.radiusOfGyration() << endl;
  cout << "Asphericity: " <<  myMolecule.asphericity() << endl;
  cout << "Acylindricity: " << myMolecule.acylindricity() << endl;
  cout << "Anisotropy: " << myMolecule.anisotropy() << endl;
  cout << "Gyration Moments: " << myMolecule.gyrationMoments() << endl;
  cout << "Gyration Axes: " <<  myMolecule.gyrationAxes() << endl;
*/ // FOR TEST
/*
  Vector3D axis(1.0, 2.0, -1.0);
  Quaternion q(axis, 60.0*3.1415926535897932384/180.0);
  myMolecule.rotate(q, 0.0);
  cout << "Rotated Principal Axes: " << myMolecule.principalAxes() << endl;
  cout << "Rotated Axes of Gyration: " << myMolecule.gyrationAxes() << endl;

  cout << "Moment of inertia around axis: " << myMolecule.momentOfInertia(axis) << endl;
  cout << "Radius of Gyration around axis: " << myMolecule.radiusOfGyration(axis) << endl;
*/
  //MOL2File my2File("test_files/hydrogen.mol2", IN);
  //Molecule my2Molecule;
  //my2File >> my2Molecule;
  //cout << my2Molecule.name() << endl;
  //cout << my2Molecule.name().length() << endl;

  // Test (distance, angle, dihedrals)

  //XYZFile gaucheFile("test_files/GEOMETRY_gauche.xyz", IN);
  //XYZFile cisFile("test_files/GEOMETRY_cis.xyz", IN);

  //Molecule gaucheMolecule, cisMolecule;

  //gaucheFile >> gaucheMolecule;
  //cisFile >> cisMolecule;

  //cout << "Gauche: distance(0, 1) = " << gaucheMolecule.distance(0, 1) << endl;
  //cout << "Gauche: angle(0, 1, 2) = " << gaucheMolecule.angle(0, 1, 2)*180.0/3.1415926535897932384 << endl;
  //cout << "Gauche: dihedral(0, 1, 2, 3) = " << gaucheMolecule.dihedral(0, 1, 2, 3)*180.0/3.1415926535897932384 << endl;

  //GradientVector myGradients;

  //myGradients = gaucheMolecule.distanceGradients(0, 1);

  //cout << "Gauche: distanceGradients(0, 1) = " << endl;
  //cout << myGradients[0] << endl;
  //cout << myGradients[1] << endl;

  //cout << "Gauche: Numerical distanceGradients(0, 1) = " << endl;

  //Vector3D grad0, grad1;
  //Real dplus, dminus;
  //Real delta = 0.001;
  //Vector3D origPos;

  //origPos = gaucheMolecule.atom(0).position;
  //gaucheMolecule.setPosition(0, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.distance(0, 1);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.distance(0, 1);
  //grad0.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.distance(0, 1);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.distance(0, 1);
  //grad0.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.distance(0, 1);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.distance(0, 1);
  //grad0.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos);

  //origPos = gaucheMolecule.atom(1).position;
  //gaucheMolecule.setPosition(1, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.distance(0, 1);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.distance(0, 1);
  //grad1.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.distance(0, 1);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.distance(0, 1);
  //grad1.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.distance(0, 1);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.distance(0, 1);
  //grad1.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos);

  //cout << grad0 << endl;
  //cout << grad1 << endl;

  //myGradients.clear();
  //myGradients = gaucheMolecule.angleGradients(0, 1, 2);

  //cout << "Gauche: angleGradients(0, 1, 2) = " << endl;
  //cout << myGradients[0] << endl;
  //cout << myGradients[1] << endl;
  //cout << myGradients[2] << endl;

  //cout << "Gauche: Numerical angleGradients(0, 1, 2) = " << endl;

  //Vector3D grad2;

  //origPos = gaucheMolecule.atom(0).position;
  //gaucheMolecule.setPosition(0, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad0.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad0.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad0.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos);

  //origPos = gaucheMolecule.atom(1).position;
  //gaucheMolecule.setPosition(1, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad1.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad1.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad1.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos);

  //origPos = gaucheMolecule.atom(2).position;
  //gaucheMolecule.setPosition(2, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(2, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad2.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(2, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(2, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad2.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(2, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.angle(0, 1, 2);
  //gaucheMolecule.setPosition(2, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.angle(0, 1, 2);
  //grad2.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(2, origPos);

  //cout << grad0 << endl;
  //cout << grad1 << endl;
  //cout << grad2 << endl;

  //myGradients.clear();
  //myGradients = gaucheMolecule.dihedralGradients(0, 1, 2, 3);

  //cout << "Gauche: dihedralGradients(0, 1, 2, 3) = " << endl;
  //cout << myGradients[0] << endl;
  //cout << myGradients[1] << endl;
  //cout << myGradients[2] << endl;
  //cout << myGradients[3] << endl;

  //cout << "Gauche: Numerical dihedralGradients(0, 1, 2, 3) = " << endl;

  //Vector3D grad3;

  //origPos = gaucheMolecule.atom(0).position;
  //gaucheMolecule.setPosition(0, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad0.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad0.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(0, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad0.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(0, origPos);

  //origPos = gaucheMolecule.atom(1).position;
  //gaucheMolecule.setPosition(1, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad1.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad1.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(1, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad1.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(1, origPos);

  //origPos = gaucheMolecule.atom(2).position;
  //gaucheMolecule.setPosition(2, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(2, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad2.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(2, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(2, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad2.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(2, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(2, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad2.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(2, origPos);

  //origPos = gaucheMolecule.atom(3).position;
  //gaucheMolecule.setPosition(3, origPos + Vector3D(delta, 0.0, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(3, origPos - Vector3D(delta, 0.0, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad3.x = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(3, origPos + Vector3D(0.0, delta, 0.0));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(3, origPos - Vector3D(0.0, delta, 0.0));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad3.y = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(3, origPos + Vector3D(0.0, 0.0, delta));
  //dplus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //gaucheMolecule.setPosition(3, origPos - Vector3D(0.0, 0.0, delta));
  //dminus = gaucheMolecule.dihedral(0, 1, 2, 3);
  //grad3.z = (dplus - dminus)/(delta + delta);
  //gaucheMolecule.setPosition(3, origPos);

  //cout << grad0 << endl;
  //cout << grad1 << endl;
  //cout << grad2 << endl;
  //cout << grad3 << endl;

  //MOL2File myMOL2File("benzene_crystal.mol2", IN);
  //DCDFile myDCDFile("benzene_cry_50K.dcd");
  //DCDFile* test = &myDCDFile;

  //Geometry myGeometry;
  //myMOL2File >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //myDCDFile >> myGeometry;

  // PSF test

  //PSFFile myPSFFile("benzene_cry.psf");
  //DCDFile myDCDFile("benzene_cry.dcd");
  //Geometry myGeometry;
  //myPSFFile >> myGeometry;
  //myDCDFile >> myGeometry;
  //MOL2File myMol2File("first_frame.mol2", OUT);
  //myMol2File << myGeometry;
  
  // New MOL2File test

  //MOL2File myFile("acetylenyl.mol2", IN);
  //Geometry myGeometry;
  //myFile >> myGeometry;
  //MOL2File myOtherFile("out.mol2", OUT);
  //myOtherFile << myGeometry;

  // New PSF and DCD test

  //PSFFile myFile("benzene_cry.psf");
  //Geometry myGeometry;
  //myFile >> myGeometry;
  //std::cout << myGeometry.numResidues();
  //DCDFile dcdFile("benzene_cry.dcd");
  //dcdFile >> myGeometry;
  //MOL2File outFile("test.mol2", OUT);
  //outFile << myGeometry;

  // MOL2File test with System<Geometry>

  //System<Molecule> mySystem;
  //MOL2File myFile("benzene_crystal.mol2", IN);
  //myFile >> mySystem;
  //MOL2File myOtherFile("out.mol2", OUT);
  //myOtherFile << mySystem;

  // New PSF and DCD test with System<Geometry>

  //System<Geometry> mySystem;
  //PSFFile myFile("benzene_cry.psf");
  //DCDFile myDCDFile("benzene_cry.dcd");
  //myFile >> mySystem;
  //myDCDFile.skipFrame();
  //myDCDFile >> mySystem;
  //MOL2File myOtherFile("out.mol2", OUT);
  //myOtherFile << mySystem;

  // TEST
  //for(size_t i = 0; i < mySystem.size(); i++)
  //  std::cout << mySystem.molecule(i).name() << std::endl;

  return 0;
}

