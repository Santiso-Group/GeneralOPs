/*
** Copyright 2007-2011 Erik Santiso.
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

/*
** Geometry class
*/

#include "common/include/assert.h"
#include "mymol/include/geometry.h"

typedef std::vector<Atom>::iterator AtomIterator;
typedef std::vector<Bond>::iterator BondIterator;

//Constructors

Geometry::Geometry()
{}

// Destructor

Geometry::~Geometry()
{}

// Find bonds automatically based on interatomic distance

void Geometry::calculateBonds(Real const &bondFactorSquared)
{
  bondTable_.clear();
  size_t nAtom = geometry_.size();

  for(size_t i = 0; i < nAtom; i++)
  for(size_t j = i+1; j < nAtom; j++)
  {
    Real sumRadius = geometry_[i].radius + geometry_[j].radius;
    Real distanceSquared = norm2(geometry_[j].position - geometry_[i].position);
    if(distanceSquared < sumRadius*sumRadius*bondFactorSquared) bondTable_.push_back(Bond(i, j));
  }
}
