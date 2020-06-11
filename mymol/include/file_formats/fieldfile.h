/*
** Copyright 2007, 2008, 2009 Erik Santiso.
** This file is part of mymol.
** mymol is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
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
** FIELD (DL_POLY) file format
**
** Notes: 
**
** When reading to Geometry/Molecule/System type objects, this class 
** only reads/writes atoms and bonds, the remaining information is
** discarded. Positions are set to zero by default. Note that reading
** some file formats (e.g. xyz or mol2) after reading a FIELD file
** will overwrite the information read from the field file. Exceptions
** include pdb, dcd, and CONFIG files. 
**
** Information on atom charges, frozen atoms or neutral/charge group numbers
** is discarded when reading atom records. The weight is also discarded,
** and silently replaced with the corresponding weight on the atom database.
**
** Constraint bonds are read as regular bonds, with type "cons".
**
** When writing from Geometry/Molecule/System type objects, this class 
** writes only atom and bond information. Consecutive molecules with the 
** same name are assumed to have the same geometry and topology even if 
** the name is an empty string. In all other cases repeated molecules 
** may be written to the FIELD file. Writing from a single geometry
** makes sense really only if the geometry contains a single molecule.
**
** Note that, even though this class reads and writes bond parameters, it
** does not check that the parameters make sense - that is the responsibility
** of the calling routine. If a constraint bond is defined with no bond
** length, the length is written as zero.
**
** The "units" label is used to allow for setting it when writing, it has no
** other role. It is set to "internal" by default. Note that this is a file
** attribute, not a geometry/molecule attribute (i.e. if reading to a geometry
** and then writing back to another file, it will be set to the default value).
**
** When reading or writing to/from PairPotential objects, this class
** ignores all the information except that in the vdw section of the FIELD
** file. Information is appended when writing, so that a FIELD file containing
** containing both the geometry and the pair potentials can be created by
** writing first from a System type object and then writing from a
** std::vector<PairPotential>. The calling routine would have to check that
** the two are consistent, however.
**
** When reading to a std::vector<PairPotential>, data is appended to the
** vector, allowing to read parameters from multiple FIELD files.
**
** Note that, as it is the case with bond information, this class does not
** check that the pair potential parameters make sense, except to the extent
** needed to convert to the internal representation.
**
** When writing The "close" keyword at the end of the file is written
** by calling the close() method - otherwise it is not. This allows to
** write multiple types of information to the same FIELD file.
**
** Given that multiple types of objects can be read from the same file,
** with overlapping information, the reading operators rewind the file
** each time they are called.
*/

#ifndef H_FIELDFILE
#define H_FIELDFILE

#include <string>
#include "common/include/iofile.h"
#include "mymol/include/geometry.h"
#include "mymol/include/molecule.h"
#include "mymol/include/system.h"
#include "mymol/include/pairpotential.h"

class FIELDFile: public IOFile
{
public:

// Constructors

  FIELDFile();                                                // Defines an empty FIELDFile
  FIELDFile(std::string const &fileName, IOMode const &mode); // Defines a FIELD File with a given file name and I/O mode

// Interface

  void setUnitsLabel(std::string const &label); // Set the "units" label (for writing)
  virtual void close();                         // Close the file

// Accessors

  std::string const &unitsLabel() const;        // Returns the "units" label


// Read to and write from Geometry objects

  friend FIELDFile &operator<<(FIELDFile &fieldFile, Geometry const &geometry); // Write geometry to file
  friend FIELDFile &operator>>(FIELDFile &fieldFile, Geometry &geometry);       // Read geometry from file

// Read to and write from System<Geometry> objects
 
  friend FIELDFile &operator<<(FIELDFile &fieldFile, System<Geometry> const &system); // Write system to file
  friend FIELDFile &operator>>(FIELDFile &fieldFile, System<Geometry> &system);       // Read system from file

// Read to and write from System<Molecule> objects 

  friend FIELDFile &operator<<(FIELDFile &fieldFile, System<Molecule> const &system); // Write system to file
  friend FIELDFile &operator>>(FIELDFile &fieldFile, System<Molecule> &system);       // Read system from file

// Read to and write from sets of PairPotentials

  friend FIELDFile &operator<<(FIELDFile &fieldFile, std::vector<PairPotential> const &potentials); // Write pair potentials to file
  friend FIELDFile &operator>>(FIELDFile &fieldFile, std::vector<PairPotential> &potentials);       // Read pair potentials from file

private:

  std::string unitsLabel_; // Label for energy units

  // Utilities to convert to/from DL_POLY names and potential types
  static std::string const dl_polyType(PotentialType const &type);
  static PotentialType const potentialType(std::string const &dl_polyType);
  // Utilities to convert potential parameters between mymol and DL_POLY formats
  static std::vector<Real> const 
    dl_polyParameters(PotentialType const &type, 
                      std::vector<Real> const &parameters);
  static std::vector<Real> const
    potentialParameters(std::string const &dl_polyType,
                        std::vector<Real> const &dl_polyParameters);
};

#endif

/*
** End of class FIELDFile
*/

// Inlines

inline void FIELDFile::setUnitsLabel(std::string const &label)
{
  unitsLabel_ = label;
}

inline std::string const &FIELDFile::unitsLabel() const
{
  return unitsLabel_;
}

inline void FIELDFile::close()
{
  if(this->mode() == OUT) *this << "close" << std::endl;
  IOFile::close();
}

