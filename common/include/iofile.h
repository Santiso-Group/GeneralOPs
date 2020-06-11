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
** Base class for file I/O
*/

#ifndef H_IOFILE
#define H_IOFILE

#define MAX_SYMBOL_LENGTH 4   // Maximum length of atomic symbol written to file
#define MAX_TYPE_LENGTH 5     // Maximum length of atomic type written to file
#define MAX_RESNAME_LENGTH 5  // Maximum length of residue name written to file (for free-format)
#define LINE_LENGTH 10000       // Maximum length of a line for reading
#define WRITE_PRECISION 12    // Precision to write doubles to file
#define WRITE_LOW_PRECISION 6 // Low precision to write doubles to file

#include <iostream>
#include <string>
#include <fstream>
#include "types.h"

enum IOMode { IN, OUT };
enum IOFormat { FORMATTED, BINARY };

class IOFile: public std::fstream
{
public:

// Constructors

  IOFile();                                                 // Defines an empty IOFile
  IOFile(std::string const &fileName, IOMode const &mode,   // Defines a file with a given file name,
         IOFormat const &format = FORMATTED);               // I/O mode and, optionally, I/O format

// Destructor

  virtual ~IOFile();  // For virtual derivation

// Interface

  virtual void setFile(std::string const &fileName, IOMode const &mode, // Sets the file name, I/O mode
                       IOFormat const &format = FORMATTED);             // and, optionally, I/O format
  virtual void clear();                                                 // Closes the file and clears file data

// Accessors

  std::string const &fileName() const;  // Returns the file name
  IOMode const &mode() const;           // Returns the I/O mode
  IOFormat const &format() const;       // Returns the I/O format

private:

  std::string fileName_;
  IOMode mode_;
  IOFormat format_;

// Prevent copying (not meaningful for fstream objects)

  IOFile(IOFile const &iofile);
  IOFile &operator=(IOFile const &iofile);
};

/*
** End of class IOFile
*/

// Inlines

inline IOFile::IOFile()
:
mode_(IN), format_(FORMATTED)
{}

inline IOFile::IOFile(std::string const &fileName, IOMode const &mode, IOFormat const &format)
{
  fileName_ = fileName;
  mode_ = mode;
  format_ = format;
  if(mode_ == IN)
  {
    if(format_ == FORMATTED)
      open(fileName.c_str(), std::ios::in);
    else
      open(fileName_.c_str(), std::ios::in|std::ios::binary);
  }
  else
  {
    if(format_ == FORMATTED)
      open(fileName.c_str(), std::ios::out);
    else
      open(fileName_.c_str(), std::ios::out|std::ios::binary);
  }
  if(fail()) 
  {
    std::cerr << "Error in IOFile: Cannot open file " << fileName << std::endl;
  }
}

inline IOFile::~IOFile()
{}

inline void IOFile::setFile(std::string const &fileName, IOMode const &mode, IOFormat const &format)
{
  if(is_open()) close();

  fileName_ = fileName;
  mode_ = mode;
  format_ = format;
  if(mode_ == IN)
  {
    if(format_ == FORMATTED)
      open(fileName.c_str(), std::ios::in);
    else
      open(fileName_.c_str(), std::ios::in|std::ios::binary);
  }
  else
  {
    if(format_ == FORMATTED)
      open(fileName.c_str(), std::ios::out);
    else
      open(fileName_.c_str(), std::ios::out|std::ios::binary);
  }
  if(fail()) 
  {
    std::cerr << "Error in IOFile: Cannot open file " << fileName << std::endl;
  }
}

inline void IOFile::clear()
{
  // Closes the file and clears file data

  if(is_open()) close();
  fileName_.clear();
  mode_ = OUT;
  format_ = FORMATTED;
}

inline std::string const &IOFile::fileName() const
{
  return fileName_;
}

inline IOMode const &IOFile::mode() const
{
  return mode_;
}

inline IOFormat const &IOFile::format() const
{
  return format_;
}

#endif
