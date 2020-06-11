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
** Some mathematical + physical constants and conversion factors
*/

#ifndef H_MPCONSTANTS
#define H_MPCONSTANTS

#include "common/include/types.h"

// Math

Real const PI = 3.1415926535897932385;
Real const TWO_PI = 6.2831853071795864770;
Real const LOG_TWO_PI = 1.8378770664093453;
Real const LOG_SQRT_TWO_PI = 0.9189385332046727;
Real const TWO_TO_ONE_THIRD = 1.25992104989487316476;
Real const TWO_TO_ONE_SIXTH = 1.12246204830937298142;

// Physics

// Conversions

Real const RAD_TO_GRAD = 57.295779513082320877;
Real const GRAD_TO_RAD = 1.7453292519943295769E-2;

#endif
