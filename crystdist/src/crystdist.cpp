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
** crystdist v. 0.1
**
** crystdist is an utility to analyze the structure of molecular crystals
** and liquid crystals.
**
** TODO:
**
** - Add description of features.
**
** - Finally, a trajectory analysis class to take over what crystal_trajectory does now
**   - maybe this should be a separate program, including the calculation of OPs etc.
** - Also, need somehow the equivalent of FTS_averages and FTS_step, but this probably is for
**   another code done especially for FTSM (with a "String" class, etc.)
**
*/

/*
** For future versions:
**
** - Use a more general treatment of molecule symmetry.
** - Need a more general way to define the molecule map (e.g. use linear combinations
**   of atomic coordinates to define the axes).
** - Right now the classes GlobalHistogram, PeakStatistics, and CrystalOrderParameters
**   are implemented for homogeneous relative configurations (both molecules are
**   of the same kind). They should be extended for the more general case of mixed
**   relative configurations.
** - It would be useful to have a function (in class RelativeConfiguration) to extract
**   all relative configurations of a given kind from a System<PointMolecule>. For
**   example, return a std::vector<RelativeConfiguration> for the pairs "GLY:GLY"
** - It'd may be better to implement histograms, peak tables, etc. using multiset/map.
**   The current way wastes time copying information to std::vectors.
** - It's also probably easier to implement operator << and >> directly for classes
**   such as PointMolecule, PeakTable, etc. so that the file I/O classes (PTMFile,
**   PKTFile, etc.) are not necessary.
** - It'd also be better to implement the internal DOF data in RelativeConfiguration
**   using something other than std::pair - it's a pain to have to repeat the same
**   code all the time for the .first and .second members. This, of course, affects
**   GlobalHistogram, PeakAccumulator, and PeakTable, at least.
** - The linear/planar case would also run much faster if, instead of storing
**   angles, the program stored their cosines or tangents. This would remove a lot
**   of calls to acos, atan2, cos and sin (the latter in PeakAccumulator). This
**   applies to angular internal DOFs as well.
** - The method calculateParameters() in PeakAccumulator is a bit of a monster.
**   It'd probably be better to move the various maximum likelihood estimation
**   routines to separate functions (in mymath?) to make it cleaner.
** - It'd make more sense, for angular quantities, to use an ANGLE type - this
**   should be implemented in mymath (DistributionVariable).
** - Regarding the two points above: DistributionVariable et al. in mymath could
**   use a lot of improvement - e.g. define an accumulator for DistributionVariable
**   that does the job of calculating means, concentration parameters, etc.
**   and adding the ANGLE type, etc. This would greatly simplify PeakAccumulator
**   and other things here.
** - Use rectangular-matrix version of the Hungarian algorithm to match peaks
**   at finite temperature to the zero temperature peaks - see note in
**   peakstatistics.h.
** - Just like calculateParameters(), the method calculate() in OrderParameterCalculator
**   is huge. This would also benefit from a better implementation of class
**   DistributionVariable in mymath, as well as the various distributions. See the
**   above comment about implementing accumulators in mymath as well. For
**   the order parameters it would be useful to have a way to calculate likelihoods
**   for a DistributionVariable.
** - There are also too many typedefs - this needs to be simplified.
**
*/

// - Think about the normalization - in glycine the relative orientation
// OPs are huge
// 1. Make a script to visualize .ptm files in vmd
// 2. Is there a faster way to skip frames in dcdfile?
// 3. Add gradients of OPs(?)


#include <fstream>

#include "mymol/include/file_formats/fileformats.h"
#include "mymol/include/system.h"
#include "crystdist/include/pointmolecule.h"
#include "crystdist/include/ptmfile.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/relativeconfiguration.h"
#include "crystdist/include/globalhistogram.h"
#include "crystdist/include/peaktable.h"
#include "crystdist/include/pktfile.h"
#include "crystdist/include/peakstatistics.h"
#include "crystdist/include/xtpfile.h"
#include "crystdist/include/opcalculator.h"
#include "crystdist/include/copfile.h"
#include "crystdist/include/peakclustering.h"

#include <ctime>

void printHelpMessage()
{
  // Prints information on how to use the program to standard output
  // OPTIONS (TILL NOW) ARE:
  // 0. CONVERT A DCD FILE TO PTM FILE.
  // 1. WHAT IS DONE IN SAVE_FILES BELOW, I.E. TAKE A MAP FILE PLUS A ZERO-TEMPERATURE DCD FILE
  // AND BUILD A GLOBALHISTOGRAM THEN A PEAK TABLE, AND WRITE THIS ONE TO FILE.
  // 2. SAME AS ABOVE BUT WITH A PTM FILE AS INPUT
  // 3. DO FINITE-TEMPERATURE STATISTICS, I.E. TAKE A FINITE-TEMPERATURE DCD FILE PLUS 
  // PEAK TABLE FILE AND WRITE THE MEANS, CONCENTRATION PARAMETERS, NORMALIZATION FACTOR(S) TO FILE
  // 4. SAME AS ABOVE BUT WITH PTM FILE INSTEAD OF DCD FILE (?)
  // 5. (FOR LATER) TAKE A DCD FILE PLUS A FILE WITH STAT PARAMETERS AND CALCULATE COP'S.
  // NOTE: TO MAKE THINGS FASTER, LET'S REMOVE OPTIONS 2 AND 4, OR AT LEAST 4 (DOESN'T MAKE A LOT
  // OF SENSE ANYWAY)
  std::cout << "crystdist - A utility to analyze molecular clusters." << std::endl << std::endl;
  std::cout << "Usage: crystdist [-h] " << std::endl
            << "                 [-d psf_filename dcd_filename] " << std::endl
            << "                 [-]" << std::endl  // VOY POR AQUI 
            << "                 [-l bead_indices]" << std::endl 
            << "                 [-o geometry_file] " << std::endl << std::endl
            << "Description: form2geom generates a geometry file (in .mol2 format)" << std::endl
            << "             corresponding to a given list of bead indices (as " << std::endl
            << "             defined in the bead_index.dat file in the beadlib" << std::endl
            << "             directory). Alternatively, it can be used to generate" << std::endl
            << "             a random molecule of given length." << std::endl << std::endl
            << "Options: " << std::endl
            << "         -h                 - Print this help message" << std::endl
            << "         -r chain_length    - Generate a random chain with the" << std::endl
            << "                              given length" << std::endl
            << "         -i bead_index_file - Generate a chain using the bead" << std::endl
            << "                              indices in the given file" << std::endl
            << "         -l bead_indices    - Generate a chain with the given" << std::endl
            << "                              bead indices" << std::endl
            << "         -o geometry_file   - Write the geometry (in .mol2" << std::endl
            << "                              format) to the given file."   << std::endl << std::endl
            << "         Exactly one of the -h, -i, -r, or -l options must be given." << std::endl
            << "         If the -o option is omitted, chain information is written" << std::endl
            << "         to standard output." << std::endl << std::endl
            << "         If a list or input file is given, the initial and final" << std::endl
            << "         beads must be terminal beads. See the files bead_index.dat" << std::endl
            << "         and readme.txt in the beadlib directory for details." << std::endl << std::endl;
}

//#define TEST_COP_IO
//#define MAKE_PTM_FILE_FOR_DEBUGGING
//#define ORDER_PARAMETERS
#define FINITE_TEMPERATURE
//#define SAVE_FILES
//#define WATER
//#define GLYCINE
#define BENZENE
//#define TPA1
//#define TPA2
#if (defined(TPA1) || defined (TPA2))
#define TPA
#endif

int main(int argc, char* argv[])
{

// Test for new molecule map

  MMPFile myMapFile("mapsnew.mmp");
  std::vector<MoleculeMap> myFMaps;
  myMapFile >> myFMaps;
  std::cout << myFMaps[0] << std::endl;
  std::cout << myFMaps[1] << std::endl;
  exit(0);

// K-Means test

  PTMFile myPTMFile("water.ptm", IN);
  System<PointMolecule> waterSystem;
  myPTMFile >> waterSystem;
  std::cout << "Done reading!" << std::endl;
  PeakClustering myClustering;
  myClustering.setCutoff(5.0);
  myClustering.setNumTrials(50);
  myClustering.addFrame(waterSystem);
  std::cout << "Done adding frames!" << std::endl;
  std::cout << myClustering.parameters() << std::endl;
  std::cout << "DBIndex: " << myClustering.dbIndex() << std::endl;
  exit(0); 

// End of K-Means test
/*
  // PTM I/O TEST 2

  PTMFile benzeneTestFile("benzene_BENZ.ptm", IN);
  System<PointMolecule> benzeneTestSystem;
  benzeneTestFile >> benzeneTestSystem;
  PTMFile copyTest("test000.ptm", OUT);
  copyTest << benzeneTestSystem ;
  return 0;
*/
/*// QUICK TEST
  MMPFile myMPFile("maps.mmp");
  std::vector<MoleculeMap> myMps;
  myMPFile >> myMps;
  MOL2File myTestFile("frame_00.mol2", IN);
  System<Molecule> myTestSystem;
  Lattice lattice(Vector3D(45.1791, 0.0, 0.0),
                  Vector3D(0.0, 48.3266, 0.0),
                  Vector3D(0.0, 0.0, 42.0214), CRYSTAL);
  myTestSystem.setLattice(lattice);
  myTestFile >> myTestSystem;
  System<PointMolecule> myPMTestSystem;
  myPMTestSystem = getPointMolecules(myTestSystem, myMps[0]);
  PTMFile myPMTestFile("benzene_T.ptm", OUT);
  myPMTestFile << myPMTestSystem;
  myPMTestFile << myPMTestSystem;
  myPMTestFile.close();
  PTMFile myPMTestFile2("benzene_T.ptm", IN);
  myPMTestSystem.clear();
  myPMTestFile2 >> myPMTestSystem;
  System<PointMolecule> myPMTestSystem2;
  myPMTestFile2 >> myPMTestSystem2;
  PTMFile myPMTestFileOut("benzene_T_1.ptm", OUT);
  PTMFile myPMTestFileOut2("benzene_T_2.ptm", OUT);
  myPMTestFileOut << myPMTestSystem;
  myPMTestFileOut2 << myPMTestSystem2;

  return 0;
*/
  // Non-planar TEST
  //myMaps[0].type = GENERAL;

  // Finite-differences test
  /*
  size_t iMol = 530;
  size_t iIndex = 2;

  Real delta = 1.E-8;
  for(size_t iCoord = 0; iCoord < 3; ++iCoord)
  {
    size_t iAtom = myMaps[0].frameAtoms[iIndex];
    Vector3D pos = myTestSystem[iMol].atom(iIndex).position;
    Vector3D incr = 0.0; incr[iCoord] += delta;
    pos += incr;
    myTestSystem[iMol].setPosition(iAtom, pos);

    System<PointMolecule> myPMTestSystem;
    myPMTestSystem = getPointMolecules(myTestSystem, myMaps[0]);
    Quaternion q_plus = myPMTestSystem[iMol].orientation;
    myPMTestSystem.clear();

    pos -= 2.0*incr;
    myTestSystem[iMol].setPosition(iAtom, pos);
    myPMTestSystem = getPointMolecules(myTestSystem, myMaps[0]);
    Quaternion q_minus = myPMTestSystem[iMol].orientation;
    myPMTestSystem.clear();

    std::cout << (q_plus - q_minus)/(2.0*delta) << std::endl;

    pos += incr;
    myTestSystem[iMol].setPosition(iAtom, pos);
  }
  */
/*  
  System<PointMolecule> myPMTestSystem;
  myPMTestSystem = getPointMolecules(myTestSystem, myMaps[0]);
  std::ofstream myTestOutput("comTest.out", std::ios::out);
  for(size_t i = 0; i < myPMTestSystem.size(); ++i)
    myTestOutput << i << " : " << myPMTestSystem[i].position << std::endl;
  myTestOutput.close();
  std::ofstream myOtherTestOutput("quatTest.out", std::ios::out);
  for(size_t i = 0; i < myPMTestSystem.size(); ++i)
    myOtherTestOutput << i << " : " << myPMTestSystem[i].orientation << std::endl;
  myOtherTestOutput.close();
  
  return 0;
*/
// END OF QUICK TEST

  enum CommandLineOption { HELP, READ_DCD, ZERO_K, FINITE_T };

#ifdef SAVE_FILES

  // For time testing
  time_t prevTime, currTime;  // Previous + current time

  std::cout << "Reading map file..." << std::endl;
  time(&prevTime);
  MMPFile myMMPFile("maps.mmp");
  std::vector<MoleculeMap> myMaps;
  myMMPFile >> myMaps;
// TEST:
  //std::cout << myMaps[2]; 

  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
  
  System<PointMolecule> system;

#ifdef GLYCINE
  PTMFile myPTMFile("glycine.ptm", OUT);
  PSFFile myPSFFile("alpha_crystal.psf");
  DCDFile myDCDFile("glycine_cry_frozen_new.dcd");
#elif defined BENZENE
  PTMFile myPTMFile("benzene.ptm", OUT);
  PSFFile myPSFFile("benzene_cry.psf");
  DCDFile myDCDFile("benzene_anneal.dcd");
#elif defined WATER
  PTMFile myPTMFile("water.ptm", OUT);
  PDBFile myPDBFile("water.pdb.10", IN);
#elif defined TPA1
  PTMFile myPTMFile("form1.ptm", OUT);
  PSFFile myPSFFile("f1_512.xplor.psf");
  DCDFile myDCDFile("f1_anneal.dcd");
#elif defined TPA2
  PTMFile myPTMFile("form2.ptm", OUT);
  PSFFile myPSFFile("f1_512.xplor.psf");
  DCDFile myDCDFile("f2_anneal.dcd");
#endif

  System<Molecule> mySystem;
  mySystem.setLatticeType(CRYSTAL); // To get PBC's right
  std::cout << "Reading psf..." << std::endl;
  time(&prevTime);
#ifndef WATER
  myPSFFile >> mySystem;
#else
  myPDBFile >> mySystem;
#endif
  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
  std::cout << "Skipping dcd frames..." << std::endl;
  time(&prevTime);
/*
#ifdef GLYCINE
  myDCDFile.skipFrames(659);  // For glycine
#elif defined BENZENE
  myDCDFile.skipFrames(1104);  // For benzene
#elif defined TPA
  myDCDFile.skipToLastFrame(); // For TPA (could use this for the others)
#endif
*/
  myDCDFile.skipToLastFrame();
  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
  std::cout << "Reading dcd..." << std::endl;
  time(&prevTime);

#ifndef WATER
  myDCDFile >> mySystem;
#endif

// TEST 2
/*
  PDBFile myYAPDBFile("benzene_last.pdb", IN);
  mySystem.clear();
  myYAPDBFile >> mySystem;
*/
// TEST!
/*
  PDBFile myYAPDBFile("g_last.pdb", OUT);
  myYAPDBFile << mySystem;
*/
  // TEST3
  PDBFile myYAPDBFile("g_last.pdb", IN);
  mySystem.clear();
  myYAPDBFile >> mySystem;
  // TEST
  //std::cout << "Lattice: " << mySystem.lattice() << std::endl;
  // OJO - WHY ARE ALPHA AND GAMMA SWAPPED?? CRAZY!
  // TEMPORARY FIX:
  /*
  Lattice latt(mySystem.lattice().a(), mySystem.lattice().b(),
               mySystem.lattice().c(), mySystem.lattice().gamma(),
               mySystem.lattice().beta(), mySystem.lattice().alpha());
  */
  //mySystem.setLattice(latt);
  // END OF TEMPORARY FIX
  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
  system.clear();
  std::cout << "Converting to System<PointMolecule>..." << std::endl;
  time(&prevTime);
#ifdef GLYCINE
  system = getPointMolecules(mySystem, myMaps[1]);
#elif defined BENZENE
  system = getPointMolecules(mySystem, myMaps[0]);
#elif defined WATER
  system = getPointMolecules(mySystem, myMaps[2]);
#elif defined TPA
  system = getPointMolecules(mySystem, myMaps[3]);
#endif
  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
  std::cout << "Writing to ptm file " << myPTMFile.fileName() << " ..." << std::endl;
  time(&prevTime);
  myPTMFile << system;
  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
  GlobalHistogram myHistogram;
  HistogramSettings mySettings = myHistogram.settings();
#ifdef GLYCINE
  mySettings.cutoff = 6.5;  // For glycine
#elif defined BENZENE
  mySettings.cutoff = 7.5;  // For benzene
#elif defined WATER
  mySettings.cutoff = 4.0; // For water
#elif defined TPA
  mySettings.cutoff = 7.5; // For TPA
#endif
  // TEST
  //std::cout << mySystem.lattice() << std::endl;
  myHistogram.changeSettings(mySettings);

// TEST
/*
  std::ifstream histSet("histSet.txt", std::ios::in);
  histSet >> mySettings;
  std::cout << mySettings; 
  return 0;
*/
  std::cout << "Building histogram..." << std::endl;
  time(&prevTime);
  // TEST
  //std::cout << system.size() << std::endl;
  //std::cout << system[0] << std::endl;
  //return 0;
#ifdef GLYCINE
  myHistogram.addSystem(system, "GLY");
#elif defined BENZENE
  myHistogram.addSystem(system, "BENZ");
#elif defined WATER
  myHistogram.addSystem(system, "WAT");
#elif defined TPA
  myHistogram.addSystem(system, "TPA");
#endif
  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
  std::cout << "Trimming histogram..." << std::endl;
  time(&prevTime);
  // TEST1
  std::ofstream myFile1("histogram1", std::ios::out);
  myFile1 << myHistogram;
#ifdef GLYCINE
  myHistogram.reduce(0.5, 2); // For glycine
#elif defined BENZENE
  myHistogram.reduce(0.05); // For benzene
#elif defined WATER
  myHistogram.reduce(0.5, 10); // For water
#elif defined TPA
  myHistogram.reduce (0.1); // For TPA
#endif
  // TEST2
  std::ofstream myFile2("histogram2", std::ios::out);
  myFile2 << myHistogram;
  time(&currTime);
  std::cout << "Time: " << difftime(currTime, prevTime) << std::endl;
// THE BELOW IS HOW TO CREATE, TRIM, AND PRINT GROUP TABLES...
  std::cout << "Building group table..." << std::endl;
  time(&prevTime);
// TEST: addSystem method...

  PeakTable myPeakTable(system, myHistogram);
  //myPeakTable.addSystem(system, myHistogram);
  //PeakTable myPeakTable;
  //myPeakTable.addSystem(system, myHistogram);

// END TEST
  //PeakTable myPeakTable(system, myHistogram);
  myPeakTable.sort();
#ifdef GLYCINE
  PKTFile myPeakFileBefore("glycine_groups_before.out", OUT);
#elif defined BENZENE
  PKTFile myPeakFileBefore("benzene_groups_before.out", OUT);
#elif defined WATER
  PKTFile myPeakFileBefore("water_groups_before.out", OUT);
#elif defined TPA1
  PKTFile myPeakFileBefore("TPA1_groups_before.out", OUT);
#elif defined TPA2
  PKTFile myPeakFileBefore("TPA2_groups_before.out", OUT);
#endif
  myPeakFileBefore << myPeakTable;
  myPeakTable.reduce(20); // <-- This minimum frequency should be defined somewhere...
  myPeakTable.sort();
#ifdef GLYCINE
  PKTFile myPeakFile("glycine_groups.out", OUT);
#elif defined BENZENE
  PKTFile myPeakFile("benzene_groups.out", OUT);
#elif defined WATER
  PKTFile myPeakFile("water_groups.out", OUT);
#elif defined TPA1
  PKTFile myPeakFile("TPA1_groups.out", OUT);
#elif defined TPA2
  PKTFile myPeakFile("TPA2_groups.out", OUT);
#endif
  myPeakFile << myPeakTable;
  myPeakFile.close();
// END OF THE GROUP STUFF

#elif defined FINITE_TEMPERATURE

  PeakTable myPeakTable;
#ifdef GLYCINE
  PKTFile myPKTFile("glycine_groups.out", IN);
#elif defined BENZENE
  PKTFile myPKTFile("benzene_groups.out", IN);
#endif

  myPKTFile >> myPeakTable;

  std::cout << "Reading map file..." << std::endl;
  MMPFile myMMPFile("maps.mmp");
  std::vector<MoleculeMap> myMaps;
  myMMPFile >> myMaps;

#ifdef GLYCINE
  PSFFile myPSFFile("alpha_crystal.psf");
  DCDFile myDCDFile("glycine_cry_300K_new.dcd");
#elif defined BENZENE
  PSFFile myPSFFile("benzene_cry.psf");
  //DCDFile myDCDFile("benzene_cry_230K.dcd");
  DCDFile myDCDFile("benzene_trimmed.dcd");
#endif

  System<PointMolecule> system;
  PeakStatistics myPeakStatistics;
  System<Molecule> mySystem;
  mySystem.setLatticeType(CRYSTAL); // To get PBC's right
  std::cout << "Reading psf..." << std::endl;
  myPSFFile >> mySystem;
  std::cout << "Skipping dcd frames..." << std::endl;
  
//#ifdef GLYCINE
//  myDCDFile.skipFrames(100);  // For glycine
//#elif defined BENZENE
//  myDCDFile.skipFrames(99);  // For benzene
//#endif
  while(myDCDFile.numFramesRead() < myDCDFile.numFrames())
//  while(myDCDFile.numFramesRead() < myDCDFile.numFrames() &&
//        myDCDFile.numFramesRead() < 600)  // To keep things quick
  {
    std::cout << "Reading dcd..." << std::endl;
    myDCDFile >> mySystem;
    system.clear();
    std::cout << "Converting to System<PointMolecule>..." << std::endl;
#ifdef GLYCINE
    system = getPointMolecules(mySystem, myMaps[1]);
#elif defined BENZENE
    system = getPointMolecules(mySystem, myMaps[0]);
#endif
    std::cout << "Adding point to statistics..." << std::endl;
    myPeakStatistics.addPoint(system, myPeakTable);
    std::cout << "Points added: " << myPeakStatistics.numPoints() << std::endl;
    /*
    std::cout << "Skipping dcd frames..." << std::endl;
    size_t skipped = 0;
    while(skipped < 9)
    {
      myDCDFile.skipFrame();
      ++skipped;
      if(myDCDFile.numFramesRead() >= myDCDFile.numFrames()) break;
    }
    */ // Commented for debugging wrapper
  }
  std::cout << "Writing parameters to file..." << std::endl;
#ifdef GLYCINE
  XTPFile myXTPFile("par_glycine.xtp", OUT);
#elif defined BENZENE
  XTPFile myXTPFile("par_benzene.xtp", OUT);
#endif
  myXTPFile << myPeakStatistics.parameters();

//  PeakStatistics myPeakStatistics;
//  myPeakStatistics.addPoint(mySystem, myPeakTable); // For debugging
//  XTPFile myXTPFile("parameters.xtp", OUT);
//  myXTPFile << myPeakStatistics.parameters();

#elif defined ORDER_PARAMETERS

  CrystalDistributionParameters myParameters;
#ifdef GLYCINE
  XTPFile myXTPFile("par_glycine.xtp", IN);
#elif defined BENZENE
  XTPFile myXTPFile("par_benzene.xtp", IN);
#endif
  std::cout << "Reading parameters..." << std::endl;
  myXTPFile >> myParameters;

  std::cout << "Reading map file..." << std::endl;
  MMPFile myMMPFile("maps.mmp");
  std::vector<MoleculeMap> myMaps;
  myMMPFile >> myMaps;
#ifdef GLYCINE
  PSFFile myPSFFile("alpha_crystal.psf");
  DCDFile myDCDFile("glycine_cry_300K_new.dcd");
#elif defined BENZENE
  PSFFile myPSFFile("benzene_cry.psf");
  //DCDFile myDCDFile("benzene_cry_230K.dcd");
  //DCDFile myDCDFile("benzene_cry.dcd");
  //DCDFile myDCDFile("frame_00.dcd");
  DCDFile myDCDFile("benzene_mess_cry.dcd");
#endif

  System<PointMolecule> system;
  PeakStatistics myPeakStatistics;
  System<Molecule> mySystem;
  mySystem.setLatticeType(CRYSTAL); // To get PBC's right
  std::cout << "Reading psf..." << std::endl;
  myPSFFile >> mySystem;
  std::cout << "Skipping dcd frames..." << std::endl;
#ifdef GLYCINE
  myDCDFile.skipFrames(100);
#elif defined BENZENE
  //myDCDFile.skipFrames(1);
//  myDCDFile.skipFrames(100);
  myDCDFile.skipFrames(150);
#endif
  std::cout << "Reading dcd..." << std::endl;
  myDCDFile >> mySystem;
  system.clear();
  std::cout << "Converting to System<PointMolecule>..." << std::endl;
#ifdef GLYCINE
  system = getPointMolecules(mySystem, myMaps[1]);
#elif defined BENZENE
  system = getPointMolecules(mySystem, myMaps[0]);
#endif

/*
  // Below is code for debugging with PTM file
  System<PointMolecule> system;
#ifdef GLYCINE
  PTMFile myPTMFile("glycine_T.ptm", IN);
#elif defined BENZENE
  PTMFile myPTMFile("benzene_T.ptm", IN);
#endif
  myPTMFile >> system;
  // End of debugging with PTM file
*/

  // TEST
  //std::cout << "TEST! " << std::endl;
  //std::cout.precision(7);
  //for(size_t i = 0; i < system.size(); ++i)
  //  std::cout << i << " " << system[i].orientation << std::endl;
  //std::cout << mySystem[0].atom(0).position << std::endl;

  // TEST
/*
  std::ofstream myCrappyFile("whatever.out", std::ios::out);
  for(size_t ii = 0; ii < system.size(); ++ii)
    myCrappyFile << system[ii].position << std::endl;
  myCrappyFile.close();

  // TEST
  PDBFile myOtherPDBFile("test_frame.pdb", OUT);
  myOtherPDBFile << mySystem;
  exit(0);
*/
  std::cout << "Calculating order parameters..." << std::endl;
  Vector3D origin(20.46, 21.4875, 19.03);
  //OrderParameterCalculator myCalculator(system, myParameters, 0.15, OrderParameterGrid(5, 5, 5, origin));
  // TEST: Per-molecule OPs
  OrderParameterCalculator myCalculator(system, myParameters, OrderParameterGrid(0, 0, 0, origin), 0.15, 0.0);
  std::cout << "Writing order parameters to file..." << std::endl;
#ifdef GLYCINE
  COPFile myCOPFile("glycine_ops.cop", OUT);
#elif defined BENZENE
//  COPFile myCOPFile("benzene_ops.cop", OUT);
  COPFile myCOPFile("benzene_0.15.cop", OUT);
#endif
  myCOPFile << myCalculator.orderParameters();

#elif defined MAKE_PTM_FILE_FOR_DEBUGGING

  std::cout << "Reading map file..." << std::endl;
  MMPFile myMMPFile("maps.mmp");
  std::vector<MoleculeMap> myMaps;
  myMMPFile >> myMaps;
#ifdef GLYCINE
  PSFFile myPSFFile("alpha_crystal.psf");
  DCDFile myDCDFile("glycine_cry_300K_new.dcd");
#elif defined BENZENE
  PSFFile myPSFFile("benzene_cry.psf");
  DCDFile myDCDFile("benzene_cry_230K.dcd");
#endif
  System<PointMolecule> system;
  PeakStatistics myPeakStatistics;
  System<Molecule> mySystem;
  mySystem.setLatticeType(CRYSTAL); // To get PBC's right
  std::cout << "Reading psf..." << std::endl;
  myPSFFile >> mySystem;
  std::cout << "Skipping dcd frames..." << std::endl;
#ifdef GLYCINE
  myDCDFile.skipFrames(100);
#elif defined BENZENE
  myDCDFile.skipFrames(100);
#endif
  std::cout << "Reading dcd..." << std::endl;
  myDCDFile >> mySystem;
  system.clear();
  std::cout << "Converting to System<PointMolecule>..." << std::endl;
#ifdef GLYCINE
  system = getPointMolecules(mySystem, myMaps[1]);
#elif defined BENZENE
  system = getPointMolecules(mySystem, myMaps[0]);
#endif
  std::cout << "Writing to .ptm file..." << std::endl;
#ifdef GLYCINE
  PTMFile myPTMFile("glycine_T.ptm", OUT);
#elif defined BENZENE
  PTMFile myPTMFile("benzene_T.ptm", OUT);
#endif
  myPTMFile << system;

#elif defined TEST_COP_IO

#ifdef GLYCINE
  COPFile myCOPInFile("glycine_ops.cop", IN);
  COPFile myCOPOutFile("glycine_ops_copy.cop", OUT);
#elif defined BENZENE
//  COPFile myCOPInFile("benzene_ops.cop", IN);
//  COPFile myCOPOutFile("benzene_ops_copy.cop", OUT);
  COPFile myCOPInFile("benzene_0.15.cop", IN);
  COPFile myCOPOutFile("benzene_0.15_copy.cop", OUT);
#endif

  CrystalOrderParameters myParameters;
  myCOPInFile >> myParameters;
  myCOPOutFile << myParameters;

  // Test operator []
  std::cout << "ops[" << 0 << "] = " << myParameters[0] << std::endl;
  std::cout << "ops[" << 17 << "] = " << myParameters[17] << std::endl;
  std::cout << "ops[" << 125 << "] = " << myParameters[125] << std::endl;
  std::cout << "ops[" << 200 << "] = " << myParameters[200] << std::endl;
  std::cout << "ops[" << 250 << "] = " << myParameters[250] << std::endl;
  std::cout << "ops[" << 300 << "] = " << myParameters[300] << std::endl;
  std::cout << "ops[" << 374 << "] = " << myParameters[374] << std::endl;
  std::cout << "ops[" << 375 << "] = " << myParameters[375] << std::endl;
  std::cout << "ops[" << 400 << "] = " << myParameters[400] << std::endl;
  std::cout << "ops[" << 499 << "] = " << myParameters[499] << std::endl;
  std::cout << "ops[" << 500 << "] = " << myParameters[500] << std::endl;

#endif

  return 0;
}
