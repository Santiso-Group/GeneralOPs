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
** .vtk (Visualization Toolkit) format
*/

#include <iomanip>
#include "common/include/assert.h"
#include "mymol/include/file_formats/fileformats.h"

// Constructors

VTKFile::VTKFile()
:
IOFile()
{}

VTKFile::VTKFile(std::string const &fileName)
:
IOFile(fileName, OUT)
{}

// Operators

VTKFile &operator<<(VTKFile &vtkFile, Grid const &grid)
{
  assert(vtkFile.is_open());

  // Write the header, title, data type, and dataset type
  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << "Grid data" << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET STRUCTURED_GRID" << std::endl;

  // Write the grid data

  vtkFile << "DIMENSIONS ";
  for(size_t i = 0; i < 3; ++i)
    vtkFile << grid.numPoints(i) << " ";
  vtkFile << std::endl;  

  size_t numPoints = grid.numPoints();
  vtkFile << "POINTS " << numPoints << " double" << std::endl;

  std::streamsize precision = vtkFile.precision(); // Store original precision
  vtkFile.precision(WRITE_PRECISION);

  for(size_t i = 0; i < numPoints; ++i)
  {
    Vector3D const &point = grid.point(i);
    vtkFile << point.x << " " << point.y << " " << point.z << std::endl;
  }

  vtkFile.precision(precision); // Return to original value

  return vtkFile;
}

VTKFile &operator<<(VTKFile &vtkFile, DiscreteField<Real> const &field)
{
  assert(vtkFile.is_open());

  // Write the header, title, data type, and dataset type
  vtkFile << "# vtk DataFile Version 2.0" << std::endl;
  vtkFile << field.name().substr(0, 255) << std::endl;
  vtkFile << "ASCII" << std::endl;
  vtkFile << "DATASET STRUCTURED_GRID" << std::endl;

  // Write the grid data
  Grid const &grid = field.grid();

  vtkFile << "DIMENSIONS ";
  for(size_t i = 0; i < 3; ++i)
    vtkFile << grid.numPoints(i) << " ";
  vtkFile << std::endl;  

  size_t numPoints = grid.numPoints();
  vtkFile << "POINTS " << numPoints << " double" << std::endl;

  std::streamsize precision = vtkFile.precision(); // Store original precision
  vtkFile.precision(WRITE_PRECISION);

  for(size_t i = 0; i < numPoints; ++i)
  {
    Vector3D const &point = grid.point(i);
    vtkFile << point.x << " " << point.y << " " << point.z << std::endl;
  }

  // Write the field data

  Real const numCells = grid.numCells();

  vtkFile << "CELL_DATA " << numCells << std::endl;
  vtkFile << "SCALARS Values double 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  for(size_t i = 0; i < numCells; ++i)
    vtkFile << field.value(i) << std::endl;

  vtkFile << "SCALARS Totals double 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  for(size_t i = 0; i < numCells; ++i)
    vtkFile << field.totalValue(i) << std::endl;

  vtkFile << "SCALARS Weights double 1" << std::endl;
  vtkFile << "LOOKUP_TABLE default" << std::endl;
  for(size_t i = 0; i < numCells; ++i)
    vtkFile << field.totalWeight(i) << std::endl;

  vtkFile.precision(precision); // Return to original value

  return vtkFile;
}

VTKFile &operator<<(VTKFile &vtkFile, DiscreteField<Vector3D> const &field)
{
  assert(vtkFile.is_open());

  // TODO: Add code here

  return vtkFile;
}

