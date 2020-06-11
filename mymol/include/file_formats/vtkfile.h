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
** .vtk (Visualization Toolkit)
**
** Notes: 
**
** This class allows the visualization of grid data with vtk-enabled
** software (e.g. ParaView, VisIt). It uses only a subset of the possible
** features of the vtk format, and is intended only for writing data,
** not for reading vtk files. Data is only written in "legacy" vtk, ASCII
** format.
**
** Note that DiscreteField is only implemented so far for Real and Vector3D
** objects, additional types are possible, add as needed.
*/

#ifndef H_VTKFILE
#define H_VTKFILE

#include "common/include/iofile.h"
#include "mymol/include/grid.h"
#include "mymol/include/discretefield.h"

class VTKFile: public IOFile
{

public:

// Constructors

  VTKFile();                            // Defines an empty VTKFile
  VTKFile(std::string const &fileName); // Defines a VTKFile with a given file name

// Write data from Grid objects

  friend VTKFile &operator<<(VTKFile &vtkFile, Grid const &grid);  // Write the grid only to file

// Write scalar field data from DiscreteField<Real> objects

  friend VTKFile &operator<<(VTKFile &vtkFile, DiscreteField<Real> const &field); // Write scalar field to file

// Write vector field data from DiscreteField<Vector3D> objects

  friend VTKFile &operator<<(VTKFile &vtkFile, DiscreteField<Vector3D> const &field); // Write vector field to file

private:

  void setFile(std::string const &fileName, // Format is write-only - blocking
               IOMode const &mode,          // the general form from IOFile
               IOFormat const &format);
};

#endif

/*
** End of class VTKFile
*/

// Inlines 

inline void VTKFile::setFile(std::string const &fileName,
                             IOMode const &mode,
                             IOFormat const &format)
{
  std::cerr << "Error in VTKFile: Reading not implemented!" << std::endl;
}


