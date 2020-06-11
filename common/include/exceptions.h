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

#ifndef H_EXCEPTIONS
#define H_EXCEPTIONS

#include <string>
#include <sstream>
#include <stdexcept>

/**
 * SimpleException class
 *
 * A simple exception class containing just an error message.
 * Can derive from this one to make more complex exceptions (e.g. see
 * DebugException below).
 */

class SimpleException: public std::exception
{
public:

// Constructors

  /// Direct constructor
  /**
   * Defines an exception with a given error message.
   */
  explicit SimpleException(std::string const &message = "An error happened!")
  :
  message_(message)
  {}

  /// Destructor
  virtual ~SimpleException() throw() {}

// Accessors

  /// Get the error message
  virtual std::string const message() const
  { return message_; }

  /// Get the error message using what()
  virtual const char* what() const throw()
  { return message_.c_str(); }

protected:

// Members

  std::string message_; // The error message
};

/*
 * End of class SimpleException
 */

/**
 * DebugException class
 *
 * An exception class that contains some extra information for tracing, such
 * as the file and the line where the error happened.
 */

class DebugException: public SimpleException
{
public:

// Constructors

  /// Direct constructor
  /**
   * Defines an exception with a file name and line number where the error
   * happened, plus a message.
   */
  DebugException(std::string const &file, 
                 int const line,
                 std::string const &message)
  :
  SimpleException(message), file_(file), line_(line)
  {}

 /// Destructor
  ~DebugException() throw() {}

// Accessors

  /// Get the full error message
  std::string const message() const
  { 
    std::stringstream msg;
    msg << message_ << " in file " << file_ << ", line " << line_;
    return msg.str();
  }

  /// Get the full error message using what()
  const char* what() const throw()
  { 
    std::stringstream msg;
    msg << message_ << " in file " << file_ << ", line " << line_;
    return msg.str().c_str();
  }

private:

  std::string file_;
  int line_;
};

/*
 * End of class DebugException
 */

#endif

