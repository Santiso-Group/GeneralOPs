/*
** Copyright 2007-2011 Erik Santiso.
** This file is part of mymath.
** mymath is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** mymath is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with mymath. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** NOTE: This file is mostly for debugging/notes to self purposes.
** It is not really part of the library
*/

/********************************************************************/
/*                                                                  */
/*  mathlib                                                         */
/*                                                                  */
/*  Simple math library                                             */
/*                                                                  */
/********************************************************************/
/*                                                                  */
/*  Some notes:                                                     */
/*                                                                  */
/*  - Matrix3D and Matrix4D could be improved by removing the calls */
/*    to operators [] and () (e.g. in the definition of matrix      */
/*    products, outer product, etc.). The way they are now is OK    */
/*    for applications that don't perform large numbers of matrix   */
/*    operations.                                                   */
/*                                                                  */
/*  - In the future, it'd probably be a good idea to replace this   */
/*    by, e.g.  blitz++.                                            */
/*                                                                  */
/********************************************************************/

#include "mymath/include/vector3D.h"
#include "mymath/include/matrix3D.h"
#include "mymath/include/vector4D.h"
#include "mymath/include/matrix4D.h"
#include "mymath/include/quaternion.h"
#include "mymath/include/vectorND.h"
#include "mymath/include/matrixND.h"
#include "mymath/include/optimization.h"
#include "mymath/include/distributionvariable.h"
#include "mymath/include/frequencytable.h"
#include "mymath/include/sorting.h"
#include "mymath/include/statistics.h"
#include "mymath/include/hungarian.h"
#include "mymath/include/rootfinding.h"
#include "mymath/include/specialfunctions.h"
#include "mymath/include/integration.h"
#include "mymath/include/interpolation.h"
#include "mymath/include/metrics.h"
#include "mymath/include/tree.h"
//#include "mymath/include/clustering.h"
#include "mymath/include/linearalgebra.h"
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

//#define VECMATIO

//#define DEBUG_CLUSTERING
//#define DEBUG_EIGENSYSTEM
//#define DEBUG_STATISTICS

#define DEBUG_TREE

#ifdef DEBUG_CLUSTERING

#include "mymath/include/clustering_example.h"

#endif

#ifdef OPTIMIZATION_EXAMPLE

#include "mymath/optimization_example.h"

#endif

#ifdef ROOT_FINDING_EXAMPLE

#include "mymath/rootfinding_examples.h"

#endif

// About sorting:
//
// Need to code a real introsort (see sorting.h)

// Things to do (maybe):
//
// 1. For next version - separate the concepts of point and vector (i.e. have separate Point3D and Vector3D)?
// This would require changes in mymol (class Geometry, Molecule - rotation of a Point3D
// requires an origin, go through the class codes).
// 2. Add utility functions for eigenvalues, etc. (Matrix3D class - advanced functions?) - done
// 2.6 Make diagonal function return the diagonal of the matrix as a vector? (both 3D and 4D) - done
// 3. Add lattice class - done (in mymol)
// 4. Add quaternion class - done
// 5. Eventually, add VectorND and MatrixND - done
// 5.1 For later, clean up eigenvalue routines - try jacobi, direct diagonalization (for 3D)...
// ******************************************************* //
// ADD MATRIX NORMS, p-NORMS OF VECTORS, ETC.
// *******************************************************//
// 5.2. Create a version that only calculates eigenvalues (useful, e.g. when only moments of inertia
// are needed, etc.)
// 6. Add utilities (as needed) to find roots of functions, optimization (BFGS, etc.) - check the
// fortran modules/parfums.

// To finish debugging quaternion class - check that matrix<->quaternion works same as fortran code - done
// One test: Check which of rotate and rotatedvector is faster... - done, fixed

// SO: TODO:
//
// - Add symmetric eigensystem to Matrix4D (this can be done later, as it is not necessary now) - done
// - Start optimization module - this is more urgent.
// - Eventually build linear algebra module - add things as needed.


// Check this:
// OJO:
//
// How to organize the code?
//
// Probably best is to keep Matrix3D and Matrix4D as they are, add SymmetricEigensystem to Matrix4D
// Then add the complicated things (like linSolve, etc.) to MatrixND.
// The way Matrix3D and Matrix4D are defined, they are really tensors, so they're conceptually different
// from a MatrixND. To make them compatible, although not pretty, there could be a constructor in
// MatrixND that constructs from Matrix3D, Matrix4D (- done).
//
// linSolve, LU, and stuff like that should not be here - put them in a separate header file
// (linear_algebra.h or something like that). Same applies to Matrix3D/4D, keep only
// containers and very simple operations that require direct access to components
// (row, columns, maybe transpose/diagonal) and move inverses etc. to the linear algebra file.
// Keep this class as close as possible to Matrix3D/4D.
// Things to move to linear_algebra:
//
// - Linsolve, LU, eigensystem, etc...

/*
// Other operations

  friend VectorND const linSolve(MatrixND const &a, VectorND const &v); // Returns the solution of the linear system
                                                                        // a*x = v by Gaussian elimination
*/

int main(int argc, char* argv[])
{
  //clock_t start, finish;
  //double time;

  //Vector4D v1(1.0, 2.0, 3.0, -5.0);
  //Vector4D v2(-2.0, 1.0, -3.0, 3.0);
  //Vector4D v3(4.0, -1.0, -2.0, 7.0);
  //Vector4D v4(2.0, 7.0, -1.0, 4.0);
  //Matrix4D a4(v1, v2, v3, v4);

  //start = clock();

  //for(int i = 0; i < 1000000; i++) { inverse(a4); }
  //
  //finish = clock();

  //time = (double)(finish-start)/CLOCKS_PER_SEC;
  //cout << time << endl;

// Mostly for debugging
    
#ifdef DEBUG_VECTOR3D

  Real arr1[3];
  SmallReal arr2[3];

  arr1[0] = -2.0; arr1[1] = 0.0; arr1[2] = 3.0;
  arr2[0] = -3.0; arr2[1] = -2.0; arr2[2] = -4.0;

  Vector3D v1;
  Vector3D v2(1.0);
  Vector3D v3(1.0, -1.0, 2.0);
  Vector3D v4(arr1);
  Vector3D v5(arr2);
  Vector3D v6(v3);

  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << "v3 = " << v3 << endl;
  cout << "v4 = " << v4 << endl;
  cout << "v5 = " << v5 << endl;
  cout << "v6 = " << v6 << endl;
  cout << " --- " << endl;
  cout << "v2 + v3 = " << v2 + v3 << endl;
  cout << "-v2 = " << -v2 << endl;
  cout << "v3 - v4 = " << v3 - v4 << endl;
  cout << "2.0*v3 = " << 2.0*v3 << endl;
  cout << "v3*3.0 = " << v3*3.0 << endl;
  cout << "v5/2.0 = " << v5/2.0 << endl;
  cout << "v3*v5 = " << v3*v5 << endl;
  cout << "dot(v3, v5) = " << dot(v3, v5) << endl;
  cout << "cross(v3, v5) = " << cross(v3, v5) << endl;
  cout << " v1 == v2 " << (v1 == v2) << endl;
  cout << " v2 != v3 " << (v2 != v3) << endl;
  cout << " v1 == 0.0 " << (v1 == 0.0) << endl;
  cout << " Setting v1 = v2 " << endl;
  v1 = v2;
  cout << " v1 == v2 " << (v1 == v2) << endl;
  cout << " v1 != v2 " << (v1 != v2) << endl;
  cout << " v2 != 1.0 " << (v2 != 1.0) << endl;
  cout << " v1 == 0.0 " << (v1 == 0.0) << endl;
  cout << " --- " << endl;
  cout << " Setting v1 = 0.0 " << endl;
  v1 = 0.0;
  cout << " v1 = " << v1 << endl;
  cout << " v2 += v3 " << endl;
  v2 += v3;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " v2 -= v3 " << endl;
  v2 -= v3;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " v2 *= 3.0 " << endl;
  v2 *= 3.0;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " v2 /= 3.0 " << endl;
  v2 /= 3.0;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " Setting v2[1] = 2.0 " << endl;
  v2[1] = 2.0;
  cout << " v2[0] = " << v2[0] << endl;
  cout << " v2[1] = " << v2[1] << endl;
  cout << " v2[2] = " << v2[2] << endl;
  cout << " Setting v2.x = 3.0 " << endl;
  v2.x = 3.0;
  cout << "v2.x = " << v2.x << endl;
  cout << " Setting v2.y = 3.0 " << endl;
  v2.y = 3.0;
  cout << "v2.y = " << v2.y << endl;
  cout << " Setting v2.z = 3.0 " << endl;
  v2.z = 3.0;
  cout << "v2.z = " << v2.z << endl;
  cout << " --- " << endl;
  cout << "Norm of v3 " << norm(v3) << endl;
  cout << "Norm2 of v3 " << norm2(v3) << endl;
  cout << "Unit vector in the direction of v3 " << unit(v3) << endl;
  cout << "Normalizing v3" << endl;
  v3.normalize();
  cout << " v3 = " << v3 << endl;
  cout << "Copying v3 into arr1" << endl;
  v3.getArray(arr1);
  cout << " arr1 = " << arr1[0] << " " << arr1[1] << " " << arr1[2] << " " << endl;
  cout << "Copying v3 into arr2" << endl;
  v3.getArray(arr2);
  cout << " arr2 = " << arr2[0] << " " << arr2[1] << " " << arr2[2] << " " << endl;

#endif

#ifdef DEBUG_MATRIX3D

  Matrix3D a1;
  Matrix3D a2(2);
  Vector3D v1(1.0, 2.0, 3.0);
  Vector3D v2(-2.0, 1.0, -3.0);
  Vector3D v3(4.0, -1.0, -2.0);
  Matrix3D a3(v1);
  Matrix3D a4(v1, v2, v3);

  Real arr1[9];
  SmallReal arr2[9];
  Real arr3[3][3];
  SmallReal arr4[3][3];

  for(int i = 0; i < 9; i++) {
    arr1[i] = (Real)i;
    arr2[i] = (SmallReal)(i*i);
  }

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      arr3[i][j] = (Real)(i + j);
      arr4[i][j] = (SmallReal)(i*i + j*j);
    }
  }
  Matrix3D a5(arr1);
  Matrix3D a6(arr2);
  Matrix3D a7(arr3);
  Matrix3D a8(arr4);
  Matrix3D a9(a4);

  cout << "a1 = " << a1 << endl;
  cout << "a2 = " << a2 << endl;
  cout << "a3 = " << a3 << endl;
  cout << "a4 = " << a4 << endl;
  cout << "a5 = " << a5 << endl;
  cout << "a6 = " << a6 << endl;
  cout << "a7 = " << a7 << endl;
  cout << "a8 = " << a8 << endl;
  cout << "a9 = " << a9 << endl;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << " --- " << endl;
  cout << "a4 + a5 = " << a4 + a5 << endl;
  cout << "-a4 = " << -a4 << endl;
  cout << "a4 - a5 = " << a4 - a5 << endl;
  cout << "2*a4 = " << 2*a4 << endl;
  cout << "a4*3 = " << a4*3 << endl;
  cout << "a4/4 = " << a4/4 << endl;
  cout << "a4*a5 = " << a4*a5 << endl;
  cout << "a4*a4 = " << a4*a4 << endl;
  cout << "a4*v1 = " << a4*v1 << endl;
  cout << "v1*a5 = " << v1*a5 << endl;
  cout << "a5^2 = " << (a5^2) << endl;
  cout << "a5*a5 = " << a5*a5 << endl;
  cout << "a4^3 = " << (a4^3) << endl;
  cout << "a4*a4*a4 = " << a4*a4*a4 << endl;
  cout << "a4^0 = " << (a4^0) << endl;
  cout << "a4^-2 = " << (a4^-2) << endl;
  cout << "inverse(a4)*inverse(a4) = " << inverse(a4)*inverse(a4) << endl;
  cout << "outer(v1, v2) = " << outer(v1, v2) << endl;
  cout << "a2 == a3 " << (a2 == a3) << endl;
  cout << "a2 != a3 " << (a2 != a3) << endl;
  cout << "a1 == 0.0 " << (a1 == 0.0) << endl;
  cout << "a2 == 2.0 " << (a2 == 2.0) << endl;
  cout << "a2 != 1.0 " << (a2 != 1.0) << endl;
  cout << "Setting a1 = a2" << endl;
  a1 = a2;
  cout << "a1 == a2 " << (a1 == a2) << endl;
  cout << "a1 != a2 " << (a1 != a2) << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 += a5" << endl;
  a6 += a5;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 -= a5" << endl;
  a6 -= a5;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 *= 3.0" << endl;
  a6 *= 3.0;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 /= 3.0" << endl;
  a6 /= 3.0;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6[5] = 1.5" << endl;
  a6[5] = 1.5;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6(0, 2) = 2.5" << endl;
  a6(0, 2) = 2.5;
  cout << "a6 = " << a6 << endl;
  cout << "a6[5] = " << a6[5] << endl;
  cout << "a6(0, 2) = " << a6(0, 2) << endl;
  cout << "Row 2 of a4 = " << a4.row(2) << endl;
  cout << "Column 1 of a4 = " << a4.column(1) << endl;
  cout << "Setting row 1 of a6 to v2 " << endl;
  a6.setRow(1, v2);
  cout << "Row 1 of a6 = " << a6.row(1) << endl;
  cout << "Setting column 2 of a6 to v1 " << endl;
  a6.setColumn(2, v1);
  cout << "Column 2 of a6 = " << a6.column(2) << endl;
  cout << " --- " << endl;
  cout << "trace of a4 = " << trace(a4) << endl;
  cout << "determinant of a4 = " << determinant(a4) << endl;
  cout << "a4 transposed = " << transpose(a4) << endl;
  cout << "a4 inverse = " << inverse(a4) << endl;
  cout << "trace of a5 = " << trace(a5) << endl;
  cout << " --- " << endl;
  cout << "Diagonal (v2) = " << diagonal(v2) << endl;
  cout << "Setting v2 = diagonal(a4)" << endl;
  v2 = diagonal(a4);
  cout << "v2 = " << v2 << endl;
  cout << "Matrix3D::Identity() = " << Matrix3D::identity() << endl;
  cout << "Matrix3D::Identity(4.0) = " << Matrix3D::identity(4.0) << endl;
  cout << " --- " << endl;
  cout << "Copying a4 into arr1" << endl;
  a4.getArray(arr1);
  cout << "arr1 = " << endl;
  cout << arr1[0] << " " << arr1[3] << " " << arr1[6] << endl;
  cout << arr1[1] << " " << arr1[4] << " " << arr1[7] << endl;
  cout << arr1[2] << " " << arr1[5] << " " << arr1[8] << endl;
  cout << "Copying a4 into arr2" << endl;
  a4.getArray(arr2);
  cout << "arr2 = " << endl;
  cout << arr2[0] << " " << arr2[3] << " " << arr2[6] << endl;
  cout << arr2[1] << " " << arr2[4] << " " << arr2[7] << endl;
  cout << arr2[2] << " " << arr2[5] << " " << arr2[8] << endl;
  cout << "Copying a4 into arr3" << endl;
  a4.getArray(arr3);
  cout << "arr3 = " << endl;
  cout << arr3[0][0] << " " << arr3[0][1] << " " << arr3[0][2] << endl;
  cout << arr3[1][0] << " " << arr3[1][1] << " " << arr3[1][2] << endl;
  cout << arr3[2][0] << " " << arr3[2][1] << " " << arr3[2][2] << endl;
  cout << "Copying a4 into arr3" << endl;
  a4.getArray(arr4);
  cout << "arr4 = " << endl;
  cout << arr4[0][0] << " " << arr4[0][1] << " " << arr4[0][2] << endl;
  cout << arr4[1][0] << " " << arr4[1][1] << " " << arr4[1][2] << endl;
  cout << arr4[2][0] << " " << arr4[2][1] << " " << arr4[2][2] << endl;

#endif

#ifdef DEBUG_QUATERNION

  Vector3D v(1.0, 2.0, -1.0);
  Vector4D v4(-1.0, 3.0, 1.0, 5.0);
  Quaternion q1;
  Quaternion q2(1.0);
  Quaternion q3(v);
  Quaternion q4(3.0, v);
  Quaternion q5(v, 3.1415926535897932384/4.0);
  Quaternion q6(-4.0, 5.0, 7.0, 3.0);
  Matrix3D a(rotationMatrix(q6));
  Quaternion q7(a);
  Quaternion q8(v4);

  cout << "v = " << v << endl;
  cout << "a = " << a << endl;
  cout << "q1 = " << q1 << endl;
  cout << "q2 = " << q2 << endl;
  cout << "q3 = " << q3 << endl;
  cout << "q4 = " << q4 << endl;
  cout << "q5 = " << q5 << endl;
  cout << "q6 = " << q6 << endl;
  cout << "q7 = " << q7 << endl;
  cout << "q8 = " << q8 << endl;
  cout << "rotationMatrix(q7) = " << rotationMatrix(q7) << endl;
  cout << "Setting q7 = 3.0" << endl;
  q7 = 3.0;
  cout << "q7 = " << q7 << endl;
  cout << "Setting q7 = v" << endl;
  q7 = v;
  cout << "q7 = " << q7 << endl;
  cout << "q4 += q6" << endl;
  q4 += q6;
  cout << "q4 = " << q4 << endl;
  cout << "q4 += 2.5" << endl;
  q4 += 2.5;
  cout << "q4 = " << q4 << endl;
  cout << "q4 += v" << endl;
  q4 += v;
  cout << "q4 = " << q4 << endl;
  cout << "q4 -= q6" << endl;
  q4 -= q6;
  cout << "q4 = " << q4 << endl;
  cout << "q4 -= 2.5" << endl;
  q4 -= 2.5;
  cout << "q4 = " << q4 << endl;
  cout << "q4 -= v" << endl;
  q4 -= v;
  cout << "q4 = " << q4 << endl;
  cout << "q4 *= 2.0" << endl;
  q4 *= 2.0;
  cout << "q4 = " << q4 << endl;
  cout << "q4 /= 2.0" << endl;
  q4 /= 2.0;
  cout << "q4 = " << q4 << endl;
  cout << "q4.normalize()" << endl;
  q4.normalize();
  cout << "q4 = " << q4 << endl;
  cout << "q4 = 3.0 + v" << endl;
  q4 = 3.0 + v;
  cout << "q4 = " << q4 << endl;
  cout << "q6 = " << q6 << endl;
  cout << "v = " << v << endl;
  cout << "q1 = q4 +  q6" << endl;
  q1 = q4 + q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4 +  3.5" << endl;
  q1 = q4 + 3.5;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = 2.5 +  q6" << endl;
  q1 = 2.5 + q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4 +  v" << endl;
  q1 = q4 + v;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = v +  q6" << endl;
  q1 = v + q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = 2.0 +  v" << endl;
  q1 = 2.0 + v;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = v +  1.0" << endl;
  q1 = v + 1.0;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = -q4" << endl;
  q1 = -q4;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = -2.0" << endl;
  q1 = -2.0;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = -v" << endl;
  q1 = -v;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4 - q6" << endl;
  q1 = q4 - q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4 - 1.0" << endl;
  q1 = q4 - 1.0;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = 7.0 - q6" << endl;
  q1 = 7.0 - q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4 - v" << endl;
  q1 = q4 - v;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = v - q6" << endl;
  q1 = v - q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = 4.0 - v" << endl;
  q1 = 4.0 - v;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = v - 3.0" << endl;
  q1 = v - 3.0;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = 5.0*q4" << endl;
  q1 = 5.0*q4;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4*3.0" << endl;
  q1 = q4*3.0;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4*q6" << endl;
  q1 = q4*q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4*v" << endl;
  q1 = q4*v;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = v*q6" << endl;
  q1 = v*q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4/4.0" << endl;
  q1 = q4/4.0;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = 2.0/q4" << endl;
  q1 = 2.0/q4;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = ~q4" << endl;
  q1 = ~q4;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4%q6" << endl;
  q1 = q4%q6;
  cout << "q1 = " << q1 << endl;
  cout << "q1 = q4*~q6" << endl;
  q1 = q4*~q6;
  cout << "q1 = " << q1 << endl;
  cout << " --- " << endl;
  cout << "q1 = q4" << endl;
  q1 = q4;
  cout << "q1 == q4 " << (q1 == q4) << endl;
  cout << "q1 != q4 " << (q1 != q4) << endl;
  cout << "q1 = q6" << endl;
  q1 = q6;
  cout << "q1 == q4 " << (q1 == q4) << endl;
  cout << "q1 != q4 " << (q1 != q4) << endl;
  cout << "q1 == 2.0 " << (q1 == 2.0) << endl;
  cout << "q1 != 2.0 " << (q1 != 2.0) << endl;
  cout << "q1 == v " << (q1 == v) << endl;
  cout << "q1 != v " << (q1 != v) << endl;
  cout << "q1 = 5.0" << endl;
  q1 = 5.0;
  cout << "q1 == 5.0 " << (q1 == 5.0) << endl;
  cout << "q1 != 5.0 " << (q1 != 5.0) << endl;
  cout << "q1 = v" << endl;
  q1 = v;
  cout << "q1 == v " << (q1 == v) << endl;
  cout << "q1 != v " << (q1 != v) << endl;
  cout << " --- " << endl;
  cout << "q4 = " << q4 << endl;
  cout << "q6 = " << q6 << endl;
  cout << "v = " << v << endl;
  cout << "dot(q4, q6) = " << dot(q4, q6) << endl;
  cout << "norm(q4) = " << norm(q4) << endl;
  cout << "norm2(q4) = " << norm2(q4) << endl;
  cout << "conjugate(q4) = " << conjugate(q4) << endl;
  cout << "scalar(q4) = " << scalar(q4) << endl;
  cout << "vector(q4) = " << vector(q4) << endl;
  cout << "axisAngle(v, pi/3) = " << axisAngle(v, 3.141592653589793/3.0) << endl;
  cout << "a = rotationMatrix(q4)" << endl;
  a = rotationMatrix(q4);
  cout << "a = " << a << endl;
  cout << "quaternion(a) = " << quaternion(a) << endl;
  cout << "rotatedVector(q4, v) = " << rotatedVector(q4, v) << endl;
  cout << "a*v = " << a*v << endl;
  cout << "rotateVector(q4, v)" << endl;
  cout << "v = " << v << endl;
  cout << "q4*v*~q4/norm2(q4) = " << q4*v*~q4/norm2(q4) << endl;

#endif

#ifdef DEBUG_VECTOR4D

  Real arr1[4];
  SmallReal arr2[4];

  arr1[0] = -2.0; arr1[1] = 0.0; arr1[2] = 3.0; arr1[3] = 5.0;
  arr2[0] = -3.0; arr2[1] = -2.0; arr2[2] = -4.0; arr2[3] = 7.0;

  Vector4D v1;
  Vector4D v2(1.0);
  Vector4D v3(1.0, -1.0, 2.0, -3.0);
  Vector4D v4(arr1);
  Vector4D v5(arr2);
  Vector4D v6(v3);

  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << "v3 = " << v3 << endl;
  cout << "v4 = " << v4 << endl;
  cout << "v5 = " << v5 << endl;
  cout << "v6 = " << v6 << endl;
  cout << " --- " << endl;
  cout << "v2 + v3 = " << v2 + v3 << endl;
  cout << "-v2 = " << -v2 << endl;
  cout << "v3 - v4 = " << v3 - v4 << endl;
  cout << "2.0*v3 = " << 2.0*v3 << endl;
  cout << "v3*3.0 = " << v3*3.0 << endl;
  cout << "v5/2.0 = " << v5/2.0 << endl;
  cout << "v3*v5 = " << v3*v5 << endl;
  cout << "dot(v3, v5) = " << dot(v3, v5) << endl;
  cout << " v1 == v2 " << (v1 == v2) << endl;
  cout << " v2 != v3 " << (v2 != v3) << endl;
  cout << " v1 == 0.0 " << (v1 == 0.0) << endl;
  cout << " Setting v1 = v2 " << endl;
  v1 = v2;
  cout << " v1 == v2 " << (v1 == v2) << endl;
  cout << " v1 != v2 " << (v1 != v2) << endl;
  cout << " v2 != 1.0 " << (v2 != 1.0) << endl;
  cout << " v1 == 0.0 " << (v1 == 0.0) << endl;
  cout << " --- " << endl;
  cout << " Setting v1 = 0.0 " << endl;
  v1 = 0.0;
  cout << " v1 = " << v1 << endl;
  cout << " v2 += v3 " << endl;
  v2 += v3;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " v2 -= v3 " << endl;
  v2 -= v3;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " v2 *= 3.0 " << endl;
  v2 *= 3.0;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " v2 /= 3.0 " << endl;
  v2 /= 3.0;
  cout << " v2 = " << v2 << endl;
  cout << " Setting v2 = 1.0 " << endl;
  v2 = 1.0;
  cout << " Setting v2[1] = 2.0 " << endl;
  v2[1] = 2.0;
  cout << " v2[0] = " << v2[0] << endl;
  cout << " v2[1] = " << v2[1] << endl;
  cout << " v2[2] = " << v2[2] << endl;
  cout << " v2[3] = " << v2[3] << endl;
  cout << " Setting v2.w = -4.0 " << endl;
  v2.w = -4.0;
  cout << "v2.w = " << v2.w << endl;
  cout << " Setting v2.x = 3.0 " << endl;
  v2.x = 3.0;
  cout << "v2.x = " << v2.x << endl;
  cout << " Setting v2.y = 3.0 " << endl;
  v2.y = 3.0;
  cout << "v2.y = " << v2.y << endl;
  cout << " Setting v2.z = 3.0 " << endl;
  v2.z = 3.0;
  cout << "v2.z = " << v2.z << endl;
  cout << " --- " << endl;
  cout << "Norm of v3 " << norm(v3) << endl;
  cout << "Norm2 of v3 " << norm2(v3) << endl;
  cout << "Unit vector in the direction of v3 " << unit(v3) << endl;
  cout << "Normalizing v3" << endl;
  v3.normalize();
  cout << " v3 = " << v3 << endl;
  cout << "Copying v3 into arr1" << endl;
  v3.getArray(arr1);
  cout << " arr1 = " << arr1[0] << " " << arr1[1] << " " << arr1[2] << " " << arr1[3] << endl;
  cout << "Copying v3 into arr2" << endl;
  v3.getArray(arr2);
  cout << " arr2 = " << arr2[0] << " " << arr2[1] << " " << arr2[2] << " " << arr2[3] << endl;

#endif

#ifdef DEBUG_MATRIX4D

  Matrix4D a1;
  Matrix4D a2(2);
  Vector4D v1(1.0, 2.0, 3.0, -5.0);
  Vector4D v2(-2.0, 1.0, -3.0, 3.0);
  Vector4D v3(4.0, -1.0, -2.0, 7.0);
  Vector4D v4(2.0, 7.0, -1.0, 4.0);
  Matrix4D a3(v1);
  Matrix4D a4(v1, v2, v3, v4);

  Real arr1[16];
  SmallReal arr2[16];
  Real arr3[4][4];
  SmallReal arr4[4][4];

  for(int i = 0; i < 16; i++) {
    arr1[i] = (Real)i;
    arr2[i] = (SmallReal)(i*i);
  }

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      arr3[i][j] = (Real)(i + j);
      arr4[i][j] = (SmallReal)(i*i + j*j);
    }
  }
  Matrix4D a5(arr1);
  Matrix4D a6(arr2);
  Matrix4D a7(arr3);
  Matrix4D a8(arr4);
  Matrix4D a9(a4);

  cout << "a1 = " << a1 << endl;
  cout << "a2 = " << a2 << endl;
  cout << "a3 = " << a3 << endl;
  cout << "a4 = " << a4 << endl;
  cout << "a5 = " << a5 << endl;
  cout << "a6 = " << a6 << endl;
  cout << "a7 = " << a7 << endl;
  cout << "a8 = " << a8 << endl;
  cout << "a9 = " << a9 << endl;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << " --- " << endl;
  cout << "a4 + a5 = " << a4 + a5 << endl;
  cout << "-a4 = " << -a4 << endl;
  cout << "a4 - a5 = " << a4 - a5 << endl;
  cout << "2*a4 = " << 2*a4 << endl;
  cout << "a4*3 = " << a4*3 << endl;
  cout << "a4/4 = " << a4/4 << endl;
  cout << "a4*a5 = " << a4*a5 << endl;
  cout << "a4*a4 = " << a4*a4 << endl;
  cout << "a4*v1 = " << a4*v1 << endl;
  cout << "v1*a5 = " << v1*a5 << endl;
  cout << "a5^2 = " << (a5^2) << endl;
  cout << "a5*a5 = " << a5*a5 << endl;
  cout << "a4^3 = " << (a4^3) << endl;
  cout << "a4*a4*a4 = " << a4*a4*a4 << endl;
  cout << "a4^0 = " << (a4^0) << endl;
  cout << "a4^-2 = " << (a4^-2) << endl;
  cout << "inverse(a4)*inverse(a4) = " << inverse(a4)*inverse(a4) << endl;
  cout << "outer(v1, v2) = " << outer(v1, v2) << endl;
  Quaternion q1(5.0, 1.0, -6.0, 9.0);
  Quaternion q2(1.0, -2.0, 7.0, 3.0);
  cout << "q1 = " << q1 << endl;
  cout << "q2 = " << q2 << endl;
  cout << "outer(q1, q2) = " << outer(q1, q2) << endl;
  cout << "a2 == a3 " << (a2 == a3) << endl;
  cout << "a2 != a3 " << (a2 != a3) << endl;
  cout << "a1 == 0.0 " << (a1 == 0.0) << endl;
  cout << "a2 == 2.0 " << (a2 == 2.0) << endl;
  cout << "a2 != 1.0 " << (a2 != 1.0) << endl;
  cout << "Setting a1 = a2" << endl;
  a1 = a2;
  cout << "a1 == a2 " << (a1 == a2) << endl;
  cout << "a1 != a2 " << (a1 != a2) << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 += a5" << endl;
  a6 += a5;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 -= a5" << endl;
  a6 -= a5;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 *= 3.0" << endl;
  a6 *= 3.0;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6 = a4" << endl;
  a6 = a4;
  cout << "a6 /= 3.0" << endl;
  a6 /= 3.0;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6[5] = 1.5" << endl;
  a6[5] = 1.5;
  cout << "a6 = " << a6 << endl;
  cout << "Setting a6(0, 2) = 2.5" << endl;
  a6(0, 2) = 2.5;
  cout << "a6 = " << a6 << endl;
  cout << "a6[5] = " << a6[5] << endl;
  cout << "a6(0, 2) = " << a6(0, 2) << endl;
  cout << "Row 2 of a4 = " << a4.row(2) << endl;
  cout << "Column 1 of a4 = " << a4.column(1) << endl;
  cout << "Setting row 1 of a6 to v2 " << endl;
  a6.setRow(1, v2);
  cout << "Row 1 of a6 = " << a6.row(1) << endl;
  cout << "Setting column 2 of a6 to v1 " << endl;
  a6.setColumn(2, v1);
  cout << "Column 2 of a6 = " << a6.column(2) << endl;
  cout << " --- " << endl;
  cout << "trace of a4 = " << trace(a4) << endl;
  cout << "determinant of a4 = " << determinant(a4) << endl;
  cout << "a4 transposed = " << transpose(a4) << endl;
  cout << "a4 inverse = " << inverse(a4) << endl;
  cout << "trace of a5 = " << trace(a5) << endl;
  cout << " --- " << endl;
  cout << "Diagonal (v2) = " << diagonal(v2) << endl;
  cout << "Setting v2 = diagonal(a4)" << endl;
  v2 = diagonal(a4);
  cout << "v2 = " << v2 << endl;
  cout << "Matrix4D::Identity() = " << Matrix4D::identity() << endl;
  cout << "Matrix4D::Identity(4.0) = " << Matrix4D::identity(4.0) << endl;
  cout << " --- " << endl;
  cout << "Copying a4 into arr1" << endl;
  a4.getArray(arr1);
  cout << "arr1 = " << endl;
  cout << arr1[0] << " " << arr1[4] << " " << arr1[ 8] << " " << arr1[12] << endl;
  cout << arr1[1] << " " << arr1[5] << " " << arr1[ 9] << " " << arr1[13] << endl;
  cout << arr1[2] << " " << arr1[6] << " " << arr1[10] << " " << arr1[14] << endl;
  cout << arr1[3] << " " << arr1[7] << " " << arr1[11] << " " << arr1[15] << endl;
  cout << "Copying a4 into arr2" << endl;
  a4.getArray(arr2);
  cout << "arr2 = " << endl;
  cout << arr2[0] << " " << arr2[4] << " " << arr2[ 8] << " " << arr2[12] << endl;
  cout << arr2[1] << " " << arr2[5] << " " << arr2[ 9] << " " << arr2[13] << endl;
  cout << arr2[2] << " " << arr2[6] << " " << arr2[10] << " " << arr2[14] << endl;
  cout << arr2[3] << " " << arr2[7] << " " << arr2[11] << " " << arr2[15] << endl;
  cout << "Copying a4 into arr3" << endl;
  a4.getArray(arr3);
  cout << "arr3 = " << endl;
  cout << arr3[0][0] << " " << arr3[0][1] << " " << arr3[0][2] << " " << arr3[0][3] << endl;
  cout << arr3[1][0] << " " << arr3[1][1] << " " << arr3[1][2] << " " << arr3[1][3] << endl;
  cout << arr3[2][0] << " " << arr3[2][1] << " " << arr3[2][2] << " " << arr3[2][3] << endl;
  cout << arr3[3][0] << " " << arr3[3][1] << " " << arr3[3][2] << " " << arr3[3][3] << endl;
  cout << "Copying a4 into arr3" << endl;
  a4.getArray(arr4);
  cout << "arr4 = " << endl;
  cout << arr4[0][0] << " " << arr4[0][1] << " " << arr4[0][2] << " " << arr4[0][3] << endl;
  cout << arr4[1][0] << " " << arr4[1][1] << " " << arr4[1][2] << " " << arr4[1][3] << endl;
  cout << arr4[2][0] << " " << arr4[2][1] << " " << arr4[2][2] << " " << arr4[2][3] << endl;
  cout << arr4[3][0] << " " << arr4[3][1] << " " << arr4[3][2] << " " << arr4[3][3] << endl;

#endif

#ifdef DEBUG_EIGENSYSTEM3D

  Matrix3D sym1;
  sym1(0, 0) = 1.0;        sym1(0, 1) = -3.0;       sym1(0, 2) = 2.0;
  sym1(1, 0) = sym1(0, 1); sym1(1, 1) = -2.0;       sym1(1, 2) = 4.0;
  sym1(2, 0) = sym1(0, 2); sym1(2, 1) = sym1(1, 2); sym1(2, 2) = 3.0;
  Matrix3D eigenvectors; Vector3D eigenvalues;
  symmetricEigensystem(sym1, eigenvalues, eigenvectors);

  cout << "Matrix:" << endl;
  cout << sym1 << endl;
  cout << "Eigenvalues:" << endl;
  cout << eigenvalues << endl;
  cout << "Eigenvectors:"  << endl;
  cout << eigenvectors << endl;

#endif

#ifdef DEBUG_EIGENSYSTEM4D

  Matrix4D sym1;
  sym1(0, 0) = 1.0;        sym1(0, 1) = -3.0;       sym1(0, 2) = 2.0;        sym1(0, 3) = 7.0;
  sym1(1, 0) = sym1(0, 1); sym1(1, 1) = -2.0;       sym1(1, 2) = 4.0;        sym1(1, 3) = 5.0;
  sym1(2, 0) = sym1(0, 2); sym1(2, 1) = sym1(1, 2); sym1(2, 2) = 3.0;        sym1(2, 3) = 3.0;
  sym1(3, 0) = sym1(0, 3); sym1(3, 1) = sym1(1, 3); sym1(3, 2) = sym1(2, 3); sym1(3, 3) = 1.0;
  Matrix4D eigenvectors; Vector4D eigenvalues;
  symmetricEigensystem(sym1, eigenvalues, eigenvectors);

  cout << "Matrix:" << endl;
  cout << sym1 << endl;
  cout << "Eigenvalues:" << endl;
  cout << eigenvalues << endl;
  cout << "Eigenvectors:"  << endl;
  cout << eigenvectors << endl;

#endif

#ifdef DEBUG_VECTORND

// Construct from std::vector<Real> - added 12/02/08

  std::vector<Real> vs1;
  vs1.push_back(1.0);
  vs1.push_back(-2.0);
  vs1.push_back(-3.0);
  vs1.push_back(1.5);

  VectorND vs2 = vs1;

  cout << vs2 << endl;

// Construct from Vector3D, Vector4D - added 09/10/08

  Vector3D v31(-1.0, 2.0, 1.0);
  Vector4D v41(1.0, -3.0, 2.0, 7.0);

  VectorND vt1(v31);
  VectorND vt2(v41);

  cout << "vt1 = " << vt1 << endl;
  cout << "vt2 = " << vt2 << endl;

  VectorND v1;
  VectorND v2(4);
  VectorND v3(5, 1.1);

  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "v2: " << endl;
  cout << v2 << endl;
  cout << "v3: " << endl;
  cout << v3 << endl;
  cout << "Resizing v2 (5, 2.0)" << endl;
  v2.resize(5, 2.0);
  cout << "v2: " << endl;
  cout << v2 << endl;
  cout << "Resizing v1 (5, 3.2)" << endl;
  v1.resize(5, 3.2);
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "Resizing v1(3)" << endl;
  v1.resize(3);
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "Setting v1 = -1.0" << endl;
  v1 = -1.0;
  cout << "v1: " << endl;
  cout << v1 << endl;
//  cout << "v2 += v1" << endl;
//  v2 += v1;  // Should crash
//  cout << "v2 -= v1" << endl;
//  v2 -= v1;  // Should crash
  cout << "v2 += v3" << endl;
  v2 += v3;
  cout << "v2: " << endl;
  cout << v2 << endl;
  cout << "Changing v2 and v3" << endl;
  for(size_t i = 0; i < 5; i++)
  {
    v2[i] = (Real)i;
    v3[i] = -(Real)i*i;
  }
  cout << "v2: " << endl;
  cout << v2 << endl;
  cout << "v3: " << endl;
  cout << v3 << endl;
  cout << "Setting v1 = v2" << endl;
  v1 = v2;
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "v1 += v3" << endl;
  v1 += v3;
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "v1 -= v2" << endl;
  v1 -= v2;
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "v1 *= 3.0" << endl;
  v1 *= 3.0;
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "v1 /= 2.0" << endl;
  v1 /= 2.0;
  cout << "v1: " << endl;
  cout << v1 << endl;
//  cout << "v1 /= 0.0" << endl;
//  v1 /= 0.0;  // Should crash
  cout << "v1[3] = " << v1[3] << endl;
  cout << "v2[0] = " << v2[0] << endl;
//  cout << "v2[9] = " << v2[9] << endl;  // Should crash
  cout << "v3[4] = " << v3[4] << endl;
  cout << "v1(2) = " << v1(2) << endl;
  cout << "v2(3) = " << v2(3) << endl;
  cout << "Setting v1[3] = 9.0" << endl;
  v1[3] = 9.0;
  cout << "v1[3] = " << v1[3] << endl;
  cout << "Setting v2(0) = 1.5" << endl;
  v2(0) = 1.5;
  cout << "v2: " << endl;
  cout << v2 << endl;
  cout << "Clearing v1" << endl;
  v1.clear();
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "Adding 1.0 and -2.0 to v1" << endl;
  v1.push_back(1.0);
  v1.push_back(-2.0);
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "v1.size() = " << v1.size() << endl;
  cout << "v2.size() = " << v2.size() << endl;
  cout << "Normalizing v1" << endl;
  v1.normalize();
  cout << "v1: " << endl;
  cout << v1 << endl;
  cout << "Normalizing v2" << endl;
  v2.normalize();
  cout << "v2: " << endl;
  cout << v2 << endl;
  cout << "Changing v2" << endl;
  v2[0] = 1.0; v2[1] = -2.0; v2[2] = 1.5; v2[3] = 4.0; v2[4] = -5.0;
  cout << "v2: " << endl;
  cout << v2 << endl;
  cout << "v3: " << endl;
  cout << v3 << endl;
//  cout << "v1 + v2: " << endl;
//  cout << v1 + v2 << endl;  // Should crash
  cout << "v2 + v3: " << endl;
  cout << v2 + v3 << endl;
  cout << "-v2: " << endl;
  cout << -v2 << endl;
//  cout << "v2 - v1" << endl;
//  cout << v2 - v1 << endl;  // Should crash
  cout << "v3 - v2: " << endl;
  cout << v3 - v2 << endl;
  cout << "3.5*v3: " << endl;
  cout << 3.5*v3 << endl;
  cout << "v2*1.5: " << endl;
  cout << v2*1.5 << endl;
//  cout << "v1*v2 = " << v1*v2 << endl;  // Should crash
  cout << "v2*v3 = " << v2*v3 << endl;
  cout << "2.5*v2 + 1.5*v3: " << endl;
  cout << 2.5*v2 + 1.5*v3 << endl;
  cout << "v1 == v2 =" << (v1 == v2) << endl;
  cout << "v1 != v2 =" << (v1 != v2) << endl;
  cout << "v2 == v3 =" << (v2 == v3) << endl;
  cout << "v2 != v3 =" << (v2 != v3) << endl;
  cout << "Setting v1 = v2" << endl;
  v1 = v2;
  cout << "v1 == v2 =" << (v1 == v2) << endl;
  cout << "v1 != v2 =" << (v1 != v2) << endl;
  cout << "v1 == 1.0 =" << (v1 == 1.0) << endl;
  cout << "v1 != 0.0 =" << (v1 != 0.0) << endl;
  cout << "Setting v1 = -1.0" << endl;
  v1 = -1.0;
  cout << "v1 == -1.0 =" << (v1 == -1.0) << endl;
  cout << "v1 != -1.0 =" << (v1 != -1.0) << endl;
  cout << "v1 == 0.0 =" << (v1 == 0.0) << endl;
  cout << "v1 != 0.0 =" << (v1 != 0.0) << endl;
//  v1.clear();
//  cout << "dot(v1, v2) = " << dot(v1, v2) << endl;  // Should crash
  cout << "dot(v2, v3) = " << dot(v2, v3) << endl;
  cout << "norm(v2) = " << norm(v2) << endl;
  cout << "norm2(v2) = " << norm2(v2) << endl;
  cout << "unit(v3): " << endl;
  cout << unit(v3) << endl;

#endif

#ifdef DEBUG_MATRIXND

// Construct from Matrix3D, Matrix4D - added 09/10/08

  Vector3D v31(1.0, 2.0, 3.0);
  Vector3D v32(-2.0, 1.0, -3.0);
  Vector3D v33(4.0, -1.0, -2.0);
  Matrix3D a31(v31, v32, v33);

  Vector4D v41(1.0, 2.0, 3.0, -5.0);
  Vector4D v42(-2.0, 1.0, -3.0, 3.0);
  Vector4D v43(4.0, -1.0, -2.0, 7.0);
  Vector4D v44(2.0, 7.0, -1.0, 4.0);
  Matrix4D a41(v41, v42, v43, v44);

  MatrixND at1(a31);
  MatrixND at2(a41);

  cout << "at1 = " << at1 << endl;
  cout << "at2 = " << at2 << endl;

  MatrixND a1;
  MatrixND a2(5);
  VectorND v1(3);
  v1[0] = 1.0; v1[1] = 2.0; v1[2] = -3.0;
  MatrixND a3(v1);
  MatrixND a4(5, 3);

  cout << "a1 = " << a1 << endl;
  cout << "a2 = " << a2 << endl;
  cout << "a3 = " << a3 << endl;
  cout << "a4 = " << a4 << endl;
  cout << "v1 = " << v1 << endl;
  cout << "Setting a1 = a3" << endl;
  a1 = a3;
  cout << "a1 = " << a1 << endl;
  cout << "Setting a2 = -1.0" << endl;
  a2 = -1.0;
  cout << "a2 = " << a2 << endl;
  cout << "Setting a4 = v1" << endl;
  a4 = v1;
  cout << "a4 = " << a4 << endl;
  cout << "Resizing a1 to 4x4" << endl;
  a1.resize(4);
  cout << "a1 = " << a1 << endl;
  cout << "Changing matrices..." << endl;
  a1(0, 0) = -1.0; a1(0, 1) =  2.0; a1(0, 2) = -1.0; a1(0, 3) =  3.0;
  a1(1, 0) =  1.0; a1(1, 1) = -3.0; a1(1, 2) =  2.0; a1(1, 3) =  5.0;
  a1(2, 0) = -2.0; a1(2, 1) =  2.0; a1(2, 2) =  7.0; a1(2, 3) = -1.0;
  a1(3, 0) =  3.0; a1(3, 1) = -5.0; a1(3, 2) =  4.0; a1(3, 3) =  4.0;
  a2.resize(4);
  a2(0, 0) =  1.0; a2(0, 1) =  5.0; a2(0, 2) = -3.0; a2(0, 3) =  2.0;
  a2(1, 0) = -3.0; a2(1, 1) =  3.0; a2(1, 2) = -7.0; a2(1, 3) =  3.0;
  a2(2, 0) =  2.0; a2(2, 1) = -2.0; a2(2, 2) =  3.0; a2(2, 3) = -4.0;
  a2(3, 0) = -5.0; a2(3, 1) =  1.0; a2(3, 2) =  1.0; a2(3, 3) =  5.0;
  a3 = MatrixND::identity(4, -3.0);
  v1.resize(4);
  v1[0] = -3.0; v1[1] = 4.0; v1[2] = -1.0; v1[3] = 2.0;
  MatrixND a5(a1);
  cout << "a1 = " << a1 << endl;
  cout << "a2 = " << a2 << endl;
  cout << "a3 = " << a3 << endl;
  cout << "a4 = " << a4 << endl;
  cout << "a5 = " << a5 << endl;
  cout << "v1 = " << v1 << endl;
  cout << "a1 += a2" << endl;
  a1 += a2;
  cout << "a1 = " << a1 << endl;
  cout << "a1 -= a2" << endl;
  a1 -= a2;
  cout << "a1 = " << a1 << endl;
  cout << "a1 *= 3.0" << endl;
  a1 *= 3.0;
  cout << "a1 = " << a1 << endl;
  cout << "a1 /= 3.0" << endl;
  a1 /= 3.0;
  cout << "a1 = " << a1 << endl;
//  cout << "a1 /= 0.0" << endl;
//  a1 /= 0.0;  // Should crash
  cout << "a1[1] = " << a1[1] << endl;
  cout << "a1(2, 3) = " << a1(2, 3) << endl;
  cout << "Setting a1(2, 3) = -9.0" << endl;
  a1(2, 3) = -9.0;
  cout << "a1 = " << a1 << endl;
  cout << "a1.row(3) = " << a1.row(3) << endl;
  cout << "a1.column(0) = " << a1.column(0) << endl;
  cout << "a1.setRow(3, v1)" << endl;
  a1.setRow(3, v1);
  cout << "a1.row(3) = " << a1.row(3) << endl;
  cout << "a1.setColumn(2, v1)" << endl;
  a1.setColumn(2, v1);
  cout << "a1.column(2) = " << a1.column(2) << endl;
  // Some crash tests: all should crash from here...
  //a1[5] = v1;
  //cout << a1.row(5);
  //cout << a1.column(4);
  //a1.setRow(6, v1);
  //a1.setRow(1, VectorND(3, 1.0));
  //a1.setColumn(8, v1);
  //a1.setColumn(1, VectorND(9, 1.0));
  // ... to here.
  cout << "a1.clear()" << endl;
  a1.clear();
  cout << "a1 = " << a1 << endl;
  cout << "Setting a1 = a5" << endl;
  a1 = a5;
  cout << "a1 = " << a1 << endl;
  cout << "a1.size() = " << a1.size() << endl;
  cout << "a1 + a2 = " << a1 + a2 << endl;
  //cout << "a1 + a4 = " << a1 + a4 << endl;  // Should crash
  cout << "-a1 = " << -a1 << endl;
  cout << "a1 - a2 = " << a1 - a2 << endl;
  //cout << "a1 - a4 = " << a1 - a4 << endl;  // Should crash
  cout << "2.5*a1 = " << 2.5*a1 << endl;
  cout << "a2*-3.5 = " << a2*-3.5 << endl;
  cout << "a1/-1.5 = " << a1/-1.5 << endl;
  //cout << "a1/0.0 = " << a1/0.0 << endl;  // Should crash
  cout << "a1*a2 = " << a1*a2 << endl;
  //cout << "a1*a4 = " << a1*a4 << endl;  // Should crash
  cout << "a1*v1 = " << a1*v1 << endl;
  cout << "v1*a2 = " << v1*a2 << endl;
  // More crash tests from here...
  //v1.resize(3);
  //cout << "v1 = " << v1 << endl;
  //cout << "a1*v1 = " << a1*v1 << endl;
  //cout << "v1*a1 = " << v1*a1 << endl;
  // to here
  cout << "a1 == a2 = " << (a1 == a2) << endl;
  cout << "a1 != a2 = " << (a1 != a2) << endl;
  cout << "a1 == a5 = " << (a1 == a5) << endl;
  cout << "a1 != a5 = " << (a1 != a5) << endl;
  cout << "a5 == 1.0 = " << (a5 == 1.0) << endl;
  cout << "a5 != 1.0 = " << (a5 != 1.0) << endl;
  cout << "Setting a5 = 2.0" << endl;
  a5 = 2.0;
  cout << "a5 == 2.0 = " << (a5 == 2.0) << endl;
  cout << "a5 != 2.0 = " << (a5 != 2.0) << endl;
  VectorND v2(4);
  v2[0] = 2.0; v2[1] = -7.0; v2[2] = 3.0; v2[3] = -1.5;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << "outer(v1, v2) = " << outer(v1, v2) << endl;
  cout << "trace(a1) = " << trace(a1) << endl;
  cout << "transpose(a1) = " << transpose(a1) << endl;
  cout << "diagonal(a1) = " << diagonal(a1) << endl;
  cout << "diagonal(v2) = " << diagonal(v2) << endl;
  cout << "identity(5, 3.0) = " << MatrixND::identity(5, 3.0) << endl;

#endif

#ifdef OPTIMIZATION_EXAMPLE

  OptimizationExample myExample;
  VectorND v_init(2, 1);  // Change this to find other minima
  VectorND v_opt = myExample.minimize(v_init);
  cout << v_opt << endl;

#endif

#ifdef DEBUG_FREQUENCY_TABLE

// These are for the old version:
/*
// Testing "Variable" struct

  std::cout << "\"Variable\" test: " << std::endl;

  size_t numBins = 10;

  Variable myVariable(0.0, 10.0, numBins);
  for(size_t i = 0; i < numBins; i++)
    std::cout << myVariable.value(i) << std::endl;

  std::cout << myVariable.minValue << " " 
            << myVariable.maxValue << " " 
            << myVariable.numBins << std::endl;

  for(Real i = 0.0; i < 10.0; i += 0.25)
    std::cout << i << " " << myVariable.index(i) << std::endl;

// Testing the FrequencyTable class

  std::cout << std::endl << "\"FrequencyTable\" test: " << std::endl;

  FrequencyTable myTable;
  myTable.addVariable(Variable(0.0, 1.0, 10));
  myTable.addVariable(Variable(0.0, 10.0, 10));

  std::vector<size_t> binIndices;
  binIndices.push_back(3);
  binIndices.push_back(7);
  myTable.addValue(binIndices);
  myTable.addValue(binIndices);
  binIndices[0] = 5;
  binIndices[1] = 0;
  myTable.addValue(binIndices);
  std::cout << myTable << std::endl;

  myTable.clearEntries();
  std::vector<Real> data;
  data.push_back(0.61);
  data.push_back(5.3);
  myTable.addValue(data);
  data.clear();
  data.push_back(0.63);
  data.push_back(5.7);
  myTable.addValue(data);
  data.clear();
  data.push_back(0.12);
  data.push_back(9.8);
  myTable.addValue(data);
  std::cout << myTable << std::endl;
*/

// Distribution variable test

  //DistributionVariable myVariable(UNIT_VECTOR_3D, 11);
  //BinIndices myIndices;
  //Vector3D myVector(-1.0, -1.0, -1.0);
  //myVector.normalize();
  //std::cout << myVector << std::endl;
  //myIndices = myVariable.unitVector3DIndices(myVector);
  //for(size_t i = 0; i < myIndices.size(); i++)
  //  std::cout << myIndices[i] << " ";
  //std::cout << std::endl;
  //std::cout << myVariable.unitVector3DValue(myIndices) << std::endl;
  //std::cout << myVariable << std::endl;
  //Quaternion myQuaternion(-1.0, -1.0, 2.0, 1.0);
  //myQuaternion.normalize();
  //std::cout << myQuaternion << std::endl;
  //myVariable.type = UNIT_QUATERNION;
  //myIndices = myVariable.unitQuaternionIndices(myQuaternion);
  //for(size_t i = 0; i < myIndices.size(); i++)
  //  std::cout << myIndices[i] << " ";
  //std::cout << std::endl;
  //std::cout << myVariable.unitQuaternionValue(myIndices) << std::endl;
  //myQuaternion = -myQuaternion;
  //std::cout << myQuaternion << std::endl;
  //myIndices = myVariable.unitQuaternionIndices(myQuaternion);
  //for(size_t i = 0; i < myIndices.size(); i++)
  //  std::cout << myIndices[i] << " ";
  //std::cout << std::endl;
  //std::cout << myVariable.unitQuaternionValue(myIndices) << std::endl;

  //std::cout << myVariable << std::endl;

// Testing the new FrequencyTable class

  FrequencyTable myTable;
  myTable.addVariable(DistributionVariable(LINEAR, 10, 0.0, 10.0));
  myTable.addVariable(DistributionVariable(UNIT_VECTOR_3D, 10));
  myTable.addVariable(DistributionVariable(UNIT_QUATERNION, 10));

  Vector3D v(1.0, -1.0, 4.0);
  v.normalize();
  Quaternion q(1.0, 2.0, -1.0, 4.0);
  q.normalize();
  std::cout << "v = " << v << std::endl;
  std::cout << "q = " << q << std::endl;
  myTable.addValue(3, 5.4, v, q);
  myTable.addValue(3, 5.3, v, q);
  v.x = -4.0; v.y = 1.0; v.z = 0.0;
  v.normalize();
  q.w = -3.0; q.x = 1.0; q.y = 2.0; q.z = 0.0;
  q.normalize();
  std::cout << "v = " << v << std::endl;
  std::cout << "q = " << q << std::endl;
  myTable.addValue(3, 10.0, v, q);
  std::cout << myTable << std::endl;
  // TEST ADDED 12/08/08
  v.x = 2.0; v.y = 2.0; v.z = 1.0;
  v.normalize();
  myTable.addValue(3, 1.3, v, q);
  myTable.deleteEntry(0);
  std::cout << ":::::" << myTable << std::endl;
  myTable.mergeEntries(0, 1);
  std::cout << myTable << std::endl;
  myTable.setVariable(1, DistributionVariable(LINEAR, 5, 0.0, 10.0));
  std::cout << myTable << std::endl;
  std::cout << "Deleting variable..." << std::endl;
  myTable.deleteVariable(1);
  std::cout << myTable << std::endl;

// A NOTE ABOUT FREQUENCYTABLE:
//
// IT WOULD BE MUCH BETTER TO GENERALIZE THE "VARIABLE" TYPE
// USING TEMPLATES, SO THAT THE TABLE COULD CONTAIN HISTOGRAMS
// WITH DIFFERENT KINDS OF OBJECTS (VECTORS, QUATERNIONS, ETC.)
// TRY THIS FOR LATER... LOOK AT distributionvariable.h FOR AN
// EXAMPLE THAT DIDN'T GO TOO WELL...

#endif

#ifdef DEBUG_SORTING

// Test insertion sort

  srand(time(0));
  std::vector<size_t> data1;
  for(size_t i = 0; i < 15; i++)
    data1.push_back(rand()%120);
  IndexTable index1 = insertionSort<std::vector<size_t>, size_t>(data1);
  std::cout << "Data: " << std::endl;
  for(size_t i = 0; i < data1.size(); i++)
    std::cout << data1[i] << " ";
  std::cout << std::endl << "Index table: " << std::endl;
  for(size_t i = 0; i < data1.size(); i++)
    std::cout << index1[i] << " ";
  std::cout << std::endl << "Sorted data: " << std::endl;
  for(size_t i = 0; i < data1.size(); i++)
    std::cout << data1[index1[i]] << " ";
  std::cout << std::endl << "Ranks: " << std::endl;
  std::vector<size_t> rank1;
  rank1 = rankTable<std::vector<size_t>, size_t>(data1, index1);
  for(size_t i = 0; i < data1.size(); i++)
    std::cout << rank1[i] << " ";
  std::cout << std::endl;

// Test introsort

  index1 = introSort<std::vector<size_t>, size_t>(data1);
  std::cout << "Index table (introsort): " << std::endl;
  for(size_t i = 0; i < index1.size(); i++)
    std::cout << index1[i] << " ";
  std::cout << std::endl;

#endif

#ifdef DEBUG_STATISTICS

  srand(time(0));
  std::vector<Real> data1;
  for(size_t i = 0; i < 15; i++)
    data1.push_back((Real)(rand()%20));

  std::cout << "Data: " << std::endl;
  for(size_t i = 0; i < 15; i++)
    std::cout << data1[i] << " ";
  std::cout << std::endl;
  std::cout << "Average: " << average<std::vector<Real>, Real>(data1) << std::endl;
  std::cout << "Variance: " << variance<std::vector<Real>, Real>(data1) << std::endl;
  std::vector<Real> meanSample = bootstrap<std::vector<Real>, Real>(data1, 200);
  std::cout << "Bootstrap distribution of the mean: " << std::endl;
  for(size_t i = 0; i < meanSample.size(); ++i)
    std::cout << meanSample[i] << std::endl;

#endif

#ifdef DEBUG_HUNGARIAN

  size_t test_case = 2; // Change this for different tests
  switch(test_case)
  {
    case(0):
    {
      MatrixND cost(3);
      for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 3; j++)
          cost(i, j) = (i+1)*(j+1);
      std::cout << hungarianMatch(cost) << std::endl;
        // Should be 2 1 0
    }
      break;
    case(1):
    {
      MatrixND cost(4);
      cost(0, 0) = 90.0;  cost(0, 1) = 75.0;  cost(0, 2) = 75.0; cost(0, 3) = 80.0;
      cost(1, 0) = 35.0;  cost(1, 1) = 85.0;  cost(1, 2) = 55.0; cost(1, 3) = 65.0;
      cost(2, 0) = 125.0; cost(2, 1) = 95.0;  cost(2, 2) = 90.0; cost(2, 3) = 105.0;
      cost(3, 0) = 45.0;  cost(3, 1) = 110.0; cost(3, 2) = 95.0; cost(3, 3) = 115.0;
      std::cout << hungarianMatch(cost) << std::endl;
        // Should be 3 2 1 0 or 1 3 2 0
    }
      break;
    case(2):
    {
      MatrixND c(6);
      c(0, 0) = 41; c(0, 1) = 72; c(0, 2) = 39; c(0, 3) = 52; c(0, 4) = 25; c(0, 5) = 51;
      c(1, 0) = 22; c(1, 1) = 29; c(1, 2) = 49; c(1, 3) = 65; c(1, 4) = 81; c(1, 5) = 50;
      c(2, 0) = 27; c(2, 1) = 39; c(2, 2) = 60; c(2, 3) = 51; c(2, 4) = 32; c(2, 5) = 32;
      c(3, 0) = 45; c(3, 1) = 50; c(3, 2) = 48; c(3, 3) = 52; c(3, 4) = 37; c(3, 5) = 43;
      c(4, 0) = 29; c(4, 1) = 40; c(4, 2) = 39; c(4, 3) = 26; c(4, 4) = 30; c(4, 5) = 33;
      c(5, 0) = 82; c(5, 1) = 40; c(5, 2) = 40; c(5, 3) = 60; c(5, 4) = 51; c(5, 5) = 30;
      std::cout << hungarianMatch(c) << std::endl;
        // Should be 4 1 0 2 3 5
    }
    default:
      break;
  }   

#endif

#ifdef ROOT_FINDING_EXAMPLE

  RootFinding1DExample myExample;
  cout << myExample.findRootNewton(0.5) << endl;
  cout << myExample.findRootSecant(0.5, 1.0) << endl;
  cout << myExample.findRootBisection(0.5, 1.0) << endl;

#endif

#ifdef DEBUG_SPECIAL_FUNCTIONS
  
  std::cout.precision(15);
  std::cout << "Testing I0: " << std::endl;
  std::cout << "I0(0.01) = " << bessel_I0(0.01) << std::endl;
  std::cout << "I0(0.1) = " << bessel_I0(0.1) << std::endl;
  std::cout << "I0(1.0) = " << bessel_I0(1.0) << std::endl;
  std::cout << "I0(10.0) = " << bessel_I0(10.0) << std::endl;
  std::cout << "I0(100.0) = " << bessel_I0(100.0) << std::endl;
  std::cout << "Testing I1: " << std::endl;
  std::cout << "I1(0.01) = " << bessel_I1(0.01) << std::endl;
  std::cout << "I1(0.1) = " << bessel_I1(0.1) << std::endl;
  std::cout << "I1(1.0) = " << bessel_I1(1.0) << std::endl;
  std::cout << "I1(10.0) = " << bessel_I1(10.0) << std::endl;
  std::cout << "I1(100.0) = " << bessel_I1(100.0) << std::endl;
  std::cout << "Testing I1/I0: " << std::endl;
  std::cout << "I1(0.01)/I0(0.01) = " << bessel_I1_over_I0(0.01) << std::endl;
  std::cout << "I1(0.1)/I0(0.1) = " << bessel_I1_over_I0(0.1) << std::endl;
  std::cout << "I1(1.0)/I0(1.0) = " << bessel_I1_over_I0(1.0) << std::endl;
  std::cout << "I1(10.0)/I0(10.0) = " << bessel_I1_over_I0(10.0) << std::endl;
  std::cout << "I1(100.0)/I0(100.0) = " << bessel_I1_over_I0(100.0) << std::endl;
  std::cout << "I1(1000.0)/I0(1000.0) = " << bessel_I1_over_I0(1000.0) << std::endl;
  std::cout << "Testing 1F1: " << std::endl;
  std::cout << "1F1(0.5, 2.0, -70.0) = " << confluent(0.5, 2.0, -70.0) << std::endl;
  std::cout << "1F1(0.5, 2.0, -10.0) = " << confluent(0.5, 2.0, -10.0) << std::endl;
  std::cout << "1F1(0.5, 2.0, 1E-6) = " << confluent(0.5, 2.0, 1E-6) << std::endl;
  std::cout << "1F1(0.5, 2.0, 0.1) = " << confluent(0.5, 2.0, 0.1) << std::endl;
  std::cout << "1F1(0.5, 2.0, 1.0) = " << confluent(0.5, 2.0, 1.0) << std::endl;
  std::cout << "1F1(0.5, 2.0, 5.0) = " << confluent(0.5, 2.0, 5.0) << std::endl;
  std::cout << "1F1(0.5, 2.0, 15.0) = " << confluent(0.5, 2.0, 15.0) << std::endl;
  std::cout << "1F1(0.5, 2.0, 20.01) = " << confluent(0.5, 2.0, 20.01) << std::endl;
  std::cout << "1F1(0.5, 2.0, 25.0) = " << confluent(0.5, 2.0, 25.0) << std::endl;
  std::cout << "1F1(0.5, 2.0, 50.0) = " << confluent(0.5, 2.0, 50.0) << std::endl;
  std::cout << "1F1(0.5, 2.0, 200.0) = " << confluent(0.5, 2.0, 200.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, -70.0) = " << confluent(1.3, 2.8, -70.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, -10.0) = " << confluent(1.3, 2.8, -10.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, 1E-6) = " << confluent(1.3, 2.8, 1E-6) << std::endl;
  std::cout << "1F1(1.3, 2.8, 0.1) = " << confluent(1.3, 2.8, 0.1) << std::endl;
  std::cout << "1F1(1.3, 2.8, 1.0) = " << confluent(1.3, 2.8, 1.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, 5.0) = " << confluent(1.3, 2.8, 5.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, 15.0) = " << confluent(1.3, 2.8, 15.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, 20.01) = " << confluent(1.3, 2.8, 20.01) << std::endl;
  std::cout << "1F1(1.3, 2.8, 25.0) = " << confluent(1.3, 2.8, 25.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, 50.0) = " << confluent(1.3, 2.8, 50.0) << std::endl;
  std::cout << "1F1(1.3, 2.8, 200.0) = " << confluent(1.3, 2.8, 200.0) << std::endl;
  std::cout << "Testing dln(1F1)/dx: " << std::endl;
  //std::cout << "dln(1F1)/dx(0.5, 2.0, -1200.0) = " << d_ln_confluent(0.5, 2.0, -1200.0) << std::endl;
  // ^^^ Should give warning...
  std::cout << "dln(1F1)/dx(0.5, 2.0, -500.0) = " << d_ln_confluent(0.5, 2.0, -500.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, -70.0) = " << d_ln_confluent(0.5, 2.0, -70.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, -10.0) = " << d_ln_confluent(0.5, 2.0, -10.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 1E-6) = " << d_ln_confluent(0.5, 2.0, 1E-6) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 0.1) = " << d_ln_confluent(0.5, 2.0, 0.1) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 1.0) = " << d_ln_confluent(0.5, 2.0, 1.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 5.0) = " << d_ln_confluent(0.5, 2.0, 5.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 15.0) = " << d_ln_confluent(0.5, 2.0, 15.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 20.01) = " << d_ln_confluent(0.5, 2.0, 20.01) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 25.0) = " << d_ln_confluent(0.5, 2.0, 25.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 50.0) = " << d_ln_confluent(0.5, 2.0, 50.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 200.0) = " << d_ln_confluent(0.5, 2.0, 200.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 900.0) = " << d_ln_confluent(0.5, 2.0, 900.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 1500.0) = " << d_ln_confluent(0.5, 2.0, 1500.0) << std::endl;
  std::cout << "dln(1F1)/dx(0.5, 2.0, 10000.0) = " << d_ln_confluent(0.5, 2.0, 10000.0) << std::endl;

  std::cout << "Testing Gamma: " << std::endl;
  std::cout << "Gamma(-21.5) = " << gamma(-21.5) << std::endl;
  std::cout << "Gamma(-11.3) = " << gamma(-11.3) << std::endl;
  std::cout << "Gamma(-1.7) = " << gamma(-1.7) << std::endl;
  std::cout << "Gamma(-1.0) = " << gamma(-1.0) << std::endl;
  std::cout << "Gamma(-0.3) = " << gamma(-0.3) << std::endl;
  std::cout << "Gamma(0.001) = " << gamma(0.001) << std::endl;
  std::cout << "Gamma(0.6) = " << gamma(0.6) << std::endl;
  std::cout << "Gamma(1.4) = " << gamma(1.4) << std::endl;
  std::cout << "Gamma(7) = " << gamma(7) << std::endl;
  std::cout << "Gamma(19.6) = " << gamma(19.6) << std::endl;
  std::cout << "Gamma(23.5) = " << gamma(23.5) << std::endl;
  std::cout << "Gamma(60) = " << gamma(60) << std::endl;
  std::cout << "Gamma(76.3) = " << gamma(76.3) << std::endl;

#endif

#ifdef DEBUG_INTEGRATION

  std::vector<std::vector<Real> > myPoints;
  std::vector<std::vector<Real> > myField;
  
  size_t numPoints = 10;
  for(size_t i = 0; i < numPoints; ++i)
  {
    Real alpha = (Real)i/(Real)(numPoints - 1);
    std::vector<Real> point(2, 0), field(2, 0);
    Real x = point[0] = cos(2*PI*alpha);
    Real y = point[1] = sin(2*PI*alpha);
    field[0] = exp(-x);
    field[1] = y;
    myPoints.push_back(point);
    myField.push_back(field);
  }
 
  typedef std::vector<std::vector<Real> > Trajectory;
 
  Real integral = trapezoidLineIntegral<Trajectory, Trajectory, Real>(myPoints, myField); 
  std::cout << "Integral = " << setprecision(15) << integral << std::endl;
  Real length = trapezoidArcLength<Trajectory, Real>(myPoints);
  std::cout << "Arc length = " << setprecision(15) << length << std::endl;

  std::vector<Real> cumulativeIntegral = 
    trapezoidCumulativeIntegral<Trajectory, Trajectory, Real>(myPoints, myField);
  std::vector<Real> cumulativeArcLength =
    trapezoidCumulativeArcLength<Trajectory, Real>(myPoints);
  std::cout << "Potential along line: " << std::endl;
  for(size_t i = 0; i < numPoints; ++i) 
    std::cout << setprecision(15)
              << cumulativeArcLength[i] << " "
              << cumulativeIntegral[i] << std::endl;

#endif

#ifdef DEBUG_INTERPOLATION

// This example is from the examples by John Burkardt.
// See http://people.sc.fsu.edu/~burkardt/cpp_src/spline/spline.html

  std::vector<Real> tdata;
  std::vector<VectorND> ydata;
  size_t numPoints = 11; 

  std::cout << "Data points: " << std::endl; 
  for(size_t i = 0; i < numPoints; ++i) 
  {
    tdata.push_back((Real)i);
    ydata.push_back(VectorND(1,
                    sin(2.0*PI*(Real)i/(Real)(numPoints - 1))));
    std::cout << setprecision(15)
              << tdata[i] << " " << ydata[i][0] << std::endl;
  }

  // Test findBracket

  std::cout << "findBracket(-0.5) = " << findBracket(tdata, -0.5) << std::endl;
  std::cout << "findBracket(0.7) = " << findBracket(tdata, 0.7) << std::endl;
  std::cout << "findBracket(3.5) = " << findBracket(tdata, 3.5) << std::endl;
  std::cout << "findBracket(8.0) = " << findBracket(tdata, 8.0) << std::endl;
  std::cout << "findBracket(10.3) = " << findBracket(tdata, 10.3) << std::endl;
  std::cout << "findBracket(15.7) = " << findBracket(tdata, 15.7) << std::endl; 

  for(size_t i = 0; i < 11; ++i)
  {
    std::cout << "findbracket(" << i << ") = "
              << findBracket(tdata, (Real)i) << std::endl;
  }

  // Test evaluateBSpline
 
  std::cout << std::endl;
  std::cout << "Interpolated values (B-spline):" << std::endl;

  Real tinit = -0.5;
  Real tfinal = 10.5;
  Real t = tinit;
  while(t <= tfinal)
  {
    std::cout << t << " " << evaluateBSpline<VectorND>(tdata, ydata, t)[0] << std::endl;
    t += 0.125;
  }

  // Test evaluateLinearSpline
   
  std::cout << std::endl;
  std::cout << "Interpolated values (linear spline):" << std::endl;

  tinit = -0.5;
  tfinal = 10.5;
  t = tinit;
  while(t <= tfinal)
  {
    std::cout << t << " " << evaluateLinearSpline<VectorND>(tdata, ydata, t)[0] << std::endl;
    t += 0.125;
  }


#endif

#ifdef DEBUG_METRICS

  // Frechet distance test
  // Example from A. Mascret et al., "Coastline Matching Process Based on the
  // Discrete Frechet Distance", in "Progress in Spatial Data Handling",
  // Part 7, Springer (2006), ISBN 978-3-540-35588-5 (Print) 978-3-540-35589-2 (Online)
  // DOI 10.1007/3-540-35589-8_25, Pages 383-400
  
  typedef std::vector<Real> Point;
  typedef std::vector<Point> Line;
  Line line1(8, std::vector<Real>(2, 0));
  Line line2(7, std::vector<Real>(2, 0));
   
  line1[0][0] = 0.2; line1[0][1] = 2.0;
  line1[1][0] = 1.5; line1[1][1] = 2.8;
  line1[2][0] = 2.3; line1[2][1] = 1.6;
  line1[3][0] = 2.9; line1[3][1] = 1.8;
  line1[4][0] = 4.1; line1[4][1] = 3.1;
  line1[5][0] = 5.6; line1[5][1] = 2.9;
  line1[6][0] = 7.2; line1[6][1] = 1.3;
  line1[7][0] = 8.2; line1[7][1] = 1.1;

  line2[0][0] = 0.3; line2[0][1] = 1.6;
  line2[1][0] = 3.2; line2[1][1] = 3.4;
  line2[2][0] = 3.8; line2[2][1] = 1.8;
  line2[3][0] = 5.2; line2[3][1] = 3.1;
  line2[4][0] = 6.5; line2[4][1] = 2.8;
  line2[5][0] = 7.0; line2[5][1] = 0.8;
  line2[6][0] = 8.9; line2[6][1] = 0.6;
/*
  //TEST - distances
  std::cout << "Distance matrix: " << std::endl;
  for(size_t i = 0; i < 8; ++i)
  for(size_t j = 0; j < 7; ++j)
    std::cout << "d[" << i << "][" << j << "] = " << distance<Point, Real>(line1[i], line2[j]) << std::endl;
  // END OF TEST - distances
*/
  std::cout << "Frechet distance = "
            << frechetDistance<Line, Point, Real>(line1, line2) << std::endl;

#endif

#ifdef DEBUG_CLUSTERING

// TODO: Add code here
  ClusteringExample myExample;
  // TEST
  for(size_t i = 0; i < myExample.numDataPoints(); ++i)
  {
    std::cout << myExample.dataPoint(i) << std::endl;
  }

#endif

#ifdef DEBUG_LU
  Vector4D v1(1.0, 2.0, 3.0, -1.0);
  Vector4D v2(1.0, 1.0, -1.0, 2.0);
  Vector4D v3(0.0, -1.0, -1.0, 3.0);
  Vector4D v4(3.0, 1.0, 2.0, -1.0);
  Matrix4D A(v1, v2, v3, v4);
  Vector4D b(4, 1, -3, 4);
  LUDecomposition<Matrix4D,Vector4D> LU(A);
  std::cout << "A = ";
  std::cout << A << std::endl;
  std::cout << "L = ";
  std::cout << LU.getL() << std::endl;
  std::cout << "U = ";
  std::cout << LU.getU() << std::endl;
  std::cout << "Pivot: " << std::endl;
  std::vector<size_t> pivot = LU.pivot();
  for(size_t i = 0; i < A.size(); ++i)
    std::cout << pivot[i] << " ";
  std::cout << std::endl;
  std::cout << "Solution of Ax = b: " << LU.solve(b) << std::endl;
  Matrix4D I = 1.0;
  Matrix4D Ainv = LU.solve(I);
  std::cout << "A^-1: " << Ainv << std::endl;
  std::cout << "A*A^-1: " << A*Ainv << std::endl;
  Vector4D a(1, -1, 3, 1);
  std::vector<Vector4D> B;
  B.push_back(a); B.push_back(b);
  std::vector<Vector4D> X = LU.solve(B);
  std::cout << "Solutions: " << std::endl;
  for(size_t i = 0; i < X.size(); ++i)
    std::cout << X[i] << std::endl;
  std::cout << "Tests (should be zero): " << std::endl;
  for(size_t i = 0; i < X.size(); ++i)
    std::cout << A*X[i] - B[i] << std::endl;
  // Testing with MatrixND
  MatrixND AA(3);
  AA(0, 0) = 1.0; AA(0, 1) = 0.0; AA(0, 2) = 2.0;
  AA(1, 0) = -1.0; AA(1, 1) = 1.0; AA(1, 2) = 2.0;
  AA(2, 0) = 2.0; AA(2, 1) = 3.0; AA(2, 2) = -2.0;
  VectorND vv(3);
  vv(0) = 1.0; vv(1) = 2.0; vv(2) = -1.0;
  LUDecomposition<MatrixND, VectorND> LLUU(AA);
  std::cout << "U = " << LLUU.getU() << std::endl;
  std::cout << "L = " << LLUU.getL() << std::endl;
  std::cout << "Solution: " << LLUU.solve(vv) << std::endl;
  std::cout << "Test (should be zero): " << AA*LLUU.solve(vv) - vv << std::endl;
#endif

#ifdef VECMATIO
  // Test input for vectorND, matrixND
  std::ifstream testV("testV.in", std::ios::in);
  VectorND vec;
  testV >> vec;
  std::cout << "V = " << std::endl << vec;
  std::ifstream testA("testA.in", std::ios::in);
  MatrixND A;
  testA >> A;
  std::cout << "A = " << std::endl << A;
#endif

#ifdef DEBUG_EIGENSYSTEM
  Vector4D v1(1.0, 2.0, 3.0, -1.0);
  Vector4D v2(2.0, 2.0, -1.0, 2.0);
  Vector4D v3(3.0, -1.0, -1.0, 3.0);
  Vector4D v4(-1.0, 2.0, 3.0, -1.0);
  Matrix4D A(v1, v2, v3, v4);
  SymmetricEigensystem<Matrix4D, Vector4D> eig(A);
  std::cout << "A = ";
  std::cout << A << std::endl;
  std::cout << "eigenvalues = ";
  std::cout << eig.eigenvalues() << std::endl;
  std::cout << "eigenvectors = ";
  std::cout << eig.eigenvectors() << std::endl;
  std::cout << "eigenvalue(2) = " << eig.eigenvalue(2) << std::endl;
  std::cout << "eigenvector(2) = " << eig.eigenvector(2) << std::endl;
  std::cout << "determinant = " << eig.determinant() << std::endl;
  std::cout << "trace = " << eig.trace() << std::endl;
  std::cout << "is it singular? " << eig.isSingular() << std::endl;

  // Testing with MatrixND
  MatrixND AA(3);
  AA(0, 0) = 1.0; AA(0, 1) = 0.0; AA(0, 2) = 2.0;
  AA(1, 0) = 0.0; AA(1, 1) = 1.0; AA(1, 2) = 3.0;
  AA(2, 0) = 2.0; AA(2, 1) = 3.0; AA(2, 2) = -2.0;
  SymmetricEigensystem<MatrixND, VectorND> eig2(AA);
  std::cout << "AA = " << std::endl << AA << std::endl;
  std::cout << "eigenvalues: " << std::endl << eig2.eigenvalues() << std::endl;
  std::cout << "eigenvectors: " << std::endl << eig2.eigenvectors() << std::endl;

#endif

#ifdef DEBUG_TREE

  // TEST - vector iterators
  /*
  std::vector<size_t> myVector(4, 1);
  myVector[1] = 2; myVector[2] = 3; myVector[3] = 4;
  std::vector<size_t>::iterator myIt;
  for(std::vector<size_t>::iterator it = myVector.begin(); 
      it != myVector.end(); ++it)
    if(*it == 2)
    {
      myIt = it;
      break;
    }
  std::cout << "Vector is: " << std::endl;
  for(size_t i = 0; i < myVector.size(); ++i)
    std::cout << myVector[i] << std::endl;
  std::cout << "Iterator points to: " << std::endl;
  std::cout << *myIt << std::endl;
  myVector.erase(myVector.begin() + 1);
  std::cout << "Vector is now: " << std::endl;
  for(size_t i = 0; i < myVector.size(); ++i)
    std::cout << myVector[i] << std::endl;
  std::cout << "Iterator points to: " << std::endl;
  std::cout << *myIt << std::endl;
  */

  
  Tree<size_t> myTree(1000);
  myTree.addChild(0, 2000);
  myTree.addChild(1, 3000);
  myTree.addChild(1, 3001);
  myTree.addChild(1, 3002);
  myTree.addChild(2, 4001);
  myTree.addChild(3, 4002);
  myTree.addChild(3, 4003);
  myTree.addChild(3, 4004);
  myTree.addChild(5, 5000);
  myTree.addChild(5, 5001);
  myTree.addChild(0, 2001);
  myTree.addChild(4, 4005);

  /*
  std::vector<size_t> leaves = myTree.getLeaves();
  std::cout << "Leaves: " << std::endl;
  for(size_t i = 0; i < leaves.size(); ++i)
    std::cout << leaves[i] << " ";
  std::cout << std::endl;
  std::cout << "Depth of the tree: " << myTree.getDepth() << std::endl;
  std::vector<size_t> nodes = myTree.getNodes();
  std::cout << "Nodes: " << std::endl;
  for(size_t i = 0; i < nodes.size(); ++i)
    std::cout << nodes[i] << " ";
  std::cout << std::endl;
  nodes = myTree.getNodesBreadthFirst();
  std::cout << "Nodes (breadth first): " << std::endl;
  for(size_t i = 0; i < nodes.size(); ++i)
    std::cout << nodes[i] << " ";
  std::cout << std::endl;
  nodes = myTree.getNodesDepthFirst();
  std::cout << "Nodes (depth first): " << std::endl;
  for(size_t i = 0; i < nodes.size(); ++i)
    std::cout << nodes[i] << " ";
  std::cout << std::endl;
  size_t const iLvl = 3;
  std::vector<size_t> level = myTree.getLevel(iLvl);
  std::cout << "Nodes at level " << iLvl << ": " << std::endl;
  for(size_t i = 0; i < level.size(); ++i)
    std::cout << level[i] << " ";
  std::cout << std::endl;
  nodes = myTree.getNodes();
  std::cout << "Parents of each node: " << std::endl;
  for(size_t i = 1; i < nodes.size(); ++i)
    std::cout << nodes[i] << ": " << myTree.getParent(nodes[i]) << std::endl;
  std::cout << "Children of each node: " << std::endl;
  for(size_t i = 0; i < nodes.size(); ++i)
  {
    std::cout << nodes[i] << ": ";
    std::vector<size_t> const children = myTree.getChildren(nodes[i]);
    for(size_t i = 0; i < children.size(); ++i)
      std::cout << children[i] << " ";
    std::cout << std::endl;
  } 
  std::cout << "Siblings of each node: " << std::endl;
  for(size_t i = 1; i < nodes.size(); ++i)
  {
    std::cout << nodes[i] << ": ";
    std::vector<size_t> const siblings = myTree.getSiblings(nodes[i]);
    for(size_t i = 0; i < siblings.size(); ++i)
      std::cout << siblings[i] << " ";
    std::cout << std::endl;
  } 
  std::cout << "Ancestors of each node: " << std::endl;
  for(size_t i = 1; i < nodes.size(); ++i)
  {
    std::cout << nodes[i] << ": ";
    std::vector<size_t> const ancestors = myTree.getAncestors(nodes[i]);
    for(size_t i = 0; i < ancestors.size(); ++i)
      std::cout << ancestors[i] << " ";
    std::cout << std::endl;
  } 
  std::cout << "Descendants of each node: " << std::endl;
  for(size_t i = 0; i < nodes.size(); ++i)
  {
    std::cout << nodes[i] << ": ";
    std::vector<size_t> const descendants = myTree.getDescendants(nodes[i]);
    for(size_t i = 0; i < descendants.size(); ++i)
      std::cout << descendants[i] << " ";
    std::cout << std::endl;
  } 
  */
  /*
  std::cout << "Tree: " << std::endl;
  std::cout << myTree;
  std::cout << "Size: " << myTree.getNumNodes() << std::endl;

  Tree<size_t> mySubtree = myTree.getSubtree(2);
  std::cout << "Subtree: " << std::endl;
  std::cout << mySubtree;
  std::cout << "Size: " << mySubtree.getNumNodes() << std::endl;

  myTree.addSubtree(4, mySubtree);
  std::cout << "After adding: " << std::endl;
  std::cout << myTree;

  myTree.prune(2);

  std::cout << "After pruning: " << std::endl;
  std::cout << myTree;
  std::cout << "Size: " << myTree.getNumNodes() << std::endl;

  myTree.addSubtree(1, mySubtree);

  std::cout << "After adding: " << std::endl;
  std::cout << myTree;
  std::cout << "Size: " << myTree.getNumNodes() << std::endl;
  */

  Tree<size_t> myOtherTree(101);
  myOtherTree.addChild(0, 201);
  myOtherTree.addChild(0, 202);
  myOtherTree.addChild(0, 203);
  myOtherTree.addChild(1, 301);
  myOtherTree.addChild(1, 302);
  myOtherTree.addChild(3, 303);
  myOtherTree.addChild(3, 304);
  myOtherTree.addChild(7, 401);

  std::cout << "Tree 1: " << std::endl;
  std::cout << myTree;
  std::cout << "Tree 2: " << std::endl;
  std::cout << myOtherTree;
  /*
  std::pair<Tree<size_t>, Tree<size_t> > newTrees = 
    Tree<size_t>::crossover(myTree, 2, myOtherTree, 3);

  std::cout << "After crossover: " << std::endl;
  
  std::cout << "Tree 1: " << std::endl;
  std::cout << newTrees.first;
  std::cout << "Tree 2: " << std::endl;
  std::cout << newTrees.second;
  */

#endif

}
