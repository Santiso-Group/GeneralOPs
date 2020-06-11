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
** Linear algebra routines - add as needed
**
** This is just a general header that automatically includes all the 
** linear algebra-related classes.
**
** Most of these routines are implemented as templates to be compatible with
** the various matrix/vector types. See the individual files for details.
**
** Many of these are adapted from the public domain JAMA library implementation 
** of the EISPACK routines available on the netlib repository.
** JAMA can be obtained at http://math.nist.gov/javanumerics/jama/, the netlib
** repository is at http://www.netlib.org.
*/

#ifndef H_LINEARALGEBRA
#define H_LINEARALGEBRA

#include "mymath/include/vector3D.h"
#include "mymath/include/vector4D.h"
#include "mymath/include/vectorND.h"
#include "mymath/include/matrix3D.h"
#include "mymath/include/matrix4D.h"
#include "mymath/include/matrixND.h"
#include "mymath/include/ludecomposition.h"
#include "mymath/include/symmetriceigensystem.h"

#endif

