/*=========================================================================

 A filter to compute the intermediate surface between two surfaces
 ( one nested into another). However, none of the approaches I have
 implemented so far is working well enough...

=========================================================================*/

#include "vtkInterpolateToSurfaceFilter.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#include "vtkCellDataToPointData.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkIdList.h"
//#include "vtkCellArray.h"
//#include "vtkAppendPolyData.h"
//#include "vtkMath.h"

#include "vtkModifiedBSPTree.h"
#include "vtkCellIterator.h"
#include "vtkGenericCell.h"

#include <algorithm>
#include <numeric>
#include <functional>

vtkStandardNewMacro(vtkInterpolateToSurfaceFilter);

vtkInterpolateToSurfaceFilter::vtkInterpolateToSurfaceFilter()
{
    this->SetNumberOfInputPorts(2);
}

vtkInterpolateToSurfaceFilter::~vtkInterpolateToSurfaceFilter()
{
}

int vtkInterpolateToSurfaceFilter::FillInputPortInformation
(
  int               port,
  vtkInformation*   info
)
{
    switch( port )
    {
        case 0:
            info->Remove( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() );
            info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData" );
        break;

        case 1:
            info->Append( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
        break;

        default:
            return 0;
    }

  return 1;
}

int vtkInterpolateToSurfaceFilter::RequestData
(
  vtkInformation        *vtkNotUsed(request),
  vtkInformationVector  **inputVector,
  vtkInformationVector  *outputVector
  )
{
    auto debugging_ctr = [](std::string extra_tag=""){static size_t ctr = 0; std::cout<<"[DEBUG] " << extra_tag << ctr << std::endl; ctr++; };

    int numInputs = this->GetNumberOfInputPorts();

    if( numInputs < 2 )
    {
        vtkErrorMacro(<<"Too few inputs provided! Need a surface to interpolate to and a cell-data field!");
        return 0;
    }

    // get the info objects
    vtkSmartPointer<vtkInformation> inInfoSurface( inputVector[0]->GetInformationObject(0) ),
                                    inInfoField( inputVector[1]->GetInformationObject(0) ),
                                    outInfoSurfaceWithData( outputVector->GetInformationObject(0) );

    // get the input and output
    vtkSmartPointer<vtkPolyData> surface( vtkPolyData::SafeDownCast( inInfoSurface->Get( vtkDataObject::DATA_OBJECT() ) ) ),
                                 surfaceWithData( vtkPolyData::SafeDownCast( outInfoSurfaceWithData->Get( vtkDataObject::DATA_OBJECT() ) ) );
    vtkSmartPointer<vtkDataSet> field( vtkDataSet::SafeDownCast( inInfoField->Get( vtkDataObject::DATA_OBJECT() ) ) );
        
    vtkIdType numPts_surf = surface->GetNumberOfPoints();
    
    std::string selected_field_name( this->GetInputArrayToProcess( 0, field )->GetName() );

    // create output surface (with point data) as copy from input surface 
    surfaceWithData->CopyStructure( surface );
    vtkSmartPointer<vtkPointData> surfaceWithData_PD( surfaceWithData->GetPointData() );
    vtkSmartPointer<vtkDataArray> surfacePointData( vtkDataArray::CreateDataArray(VTK_DOUBLE)  );
    surfacePointData->SetNumberOfComponents( 1 );
    surfacePointData->SetNumberOfTuples( numPts_surf );
    surfacePointData->SetName( (selected_field_name + "_interpolated").c_str() );
    surfaceWithData_PD->SetScalars( surfacePointData );
    surfacePointData->Fill( -1. );

    // Interpolate cell data of the volume mesh to its mesh points using VTK's built-in filter (vtkCellDataToPointData)
    vtkSmartPointer<vtkCellDataToPointData> cellDToPointDConverter = vtkSmartPointer<vtkCellDataToPointData>::New();
    cellDToPointDConverter->SetInputData( field );
    cellDToPointDConverter->Update();
    vtkSmartPointer<vtkDataSet> field_copy( cellDToPointDConverter->GetOutput() );
    vtkSmartPointer<vtkPointData> field_copy_PD( field_copy->GetPointData() ); 
    // use the array (denoted by its name) specified in the GUI (StringProperty, SetInputDataToProcess)
    int dump = 0;
    vtkSmartPointer<vtkDataArray> fieldPointData( field_copy_PD->GetArray(selected_field_name.c_str(),  dump));

    // (dummy) variables to populate during iteration// getCell(0) - lets hope there will always be at least one cell...
    vtkIdType numPointsPerCell = field_copy->GetCell(0)->GetNumberOfPoints(), dummyCellID = 0, containingCellID = 0; int subID = 0;
    double currentPoint[3], pCoords [3], weights[ numPointsPerCell ],pointVal[1] = {0}, interp_val=0;// pCoords=parametric coordinates, weights=interpolation weights
    vtkSmartPointer<vtkCell> dummyCell(nullptr);
    vtkSmartPointer<vtkCell> containingCell(nullptr);
    vtkSmartPointer<vtkIdList> pointIDsOfContainingCell(nullptr);

    // for every point of the surface, interpolate from the field's data
    for( vtkIdType pointSeqNo = 0; pointSeqNo < numPts_surf; pointSeqNo++ )
    {
        surface->GetPoint( pointSeqNo, currentPoint );

        //std::cout << "pt(" << currentPoint[0] << "," << currentPoint[1] << "," << currentPoint[2] << ")" << std::endl;
        containingCellID = field_copy->FindCell( currentPoint, dummyCell, dummyCellID, 1e-05, subID, pCoords, weights );
        if( containingCellID >= 0 )
        {
            containingCell = field_copy->GetCell( containingCellID );
            if( containingCell )
            {
                pointIDsOfContainingCell = containingCell->GetPointIds();
                
                // loop over points of cell
                for( vtkIdType cellPointID = 0; cellPointID < pointIDsOfContainingCell->GetNumberOfIds(); cellPointID++ )
                {
                    // obtain point-data value from global point ID
                    fieldPointData->GetTuple( pointIDsOfContainingCell->GetId( cellPointID ), pointVal );

                    interp_val += pointVal[0] * weights[ cellPointID ];
                }
                surfacePointData->SetTuple1( pointSeqNo, interp_val );
                
                // Try to make it abit more efficient. From the docs:
                // "If cell and cellId is non-nullptr, then search starts from this cell and looks at immediate neighbors."
                dummyCell = containingCell;
                dummyCellID = containingCellID;
                interp_val = 0.;
            }
            else
            {
          //      std::cout << "containingCell = nullptr" << std::endl;
            }
        }
        else
        {
            //std::cout << "Cell-ID < 0" << std::endl;
        }
    }

    return 1; // strangely ParaView/VTK interprets 1=EXIT_FAILURE as success, and 0=EXIT_SUCCESS as failure...
}

void vtkInterpolateToSurfaceFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
