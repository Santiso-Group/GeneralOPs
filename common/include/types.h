/*
** Copyright 2007-2011 Erik Santiso.
** This is a common file included in mymath, mymol, etc...
** This is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
**
** This is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with this. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** Type definitions
*/

#ifndef H_TYPES
#define H_TYPES

#include <float.h>

// typedef long double BigReal;
typedef double Real;
typedef float  SmallReal;

Real const REAL_EPSILON = DBL_EPSILON;
SmallReal const SMALLREAL_EPSILON = FLT_EPSILON;
Real const REAL_VERYBIG = DBL_MAX;
SmallReal const SMALLREAL_VERYBIG = FLT_MAX;

#endif

