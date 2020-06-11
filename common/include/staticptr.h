/*
 * Copyright (C) 2012 Erik Santiso.
 * This is a common file included in mymath, mymol, etc...
 *
 * This is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef H_STATICPTR
#define H_STATICPTR

#include <tr1/memory>

/**
 * Null deleter
 *
 * Used to prevent shared_ptr from attempting to delete a statically allocated
 * object, as explained in the Boost library documentation. See notes in the 
 * README file.
 */

struct NullDeleter
{
  void operator()(void const *) const {}
};


/**
 * Static pointer
 *
 * Implements a shared_ptr to a statically allocated object, as explained in
 * the Boost library documentation.
 */

template<typename Type>
std::tr1::shared_ptr<Type> StaticPtr(Type *object)
{
  return std::tr1::shared_ptr<Type>(object, NullDeleter());
}

#endif

