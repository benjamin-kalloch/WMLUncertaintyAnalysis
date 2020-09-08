/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkInterpolateToSurfaceFilter - interpolate a cell-data field onto a surface as point-data

// .SECTION Description
// A Class to interpolate the cell data of a volume to the vertices of a
// surface mesh that positioned inside the volume mesh.
// Inputs: 1) vtkPolyData = surface, 2) vtk 

// .SECTION Thanks
// noone

#ifndef vtkInterpolateToSurfaceFilter_h
#define vtkInterpolateToSurfaceFilter_h

#include "vtkDataSetAlgorithm.h"

#include "vtkSmartPointer.h"
#include <set>
#include <vector>

class vtkInterpolateToSurfaceFilter : public vtkDataSetAlgorithm
{
public:
  vtkGetMacro(selected_field_data, int);
  vtkSetMacro(selected_field_data, int);

  vtkTypeMacro(vtkInterpolateToSurfaceFilter, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkInterpolateToSurfaceFilter *New();

protected:
  vtkInterpolateToSurfaceFilter();
  ~vtkInterpolateToSurfaceFilter();

  /* implementation of algorithm */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

protected:
  double selected_field_data;

private:
  vtkInterpolateToSurfaceFilter(const vtkInterpolateToSurfaceFilter&);  // Not implemented.
  void operator=(const vtkInterpolateToSurfaceFilter&);  // Not implemented.
};

#endif
