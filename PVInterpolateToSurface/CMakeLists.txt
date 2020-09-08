cmake_minimum_required(VERSION 3.10)

FIND_PACKAGE(ParaView REQUIRED)

INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN(InterpolateToSurface "1.0" 
   SERVER_MANAGER_XML vtkInterpolateToSurfaceFilter_SM.xml
   SERVER_MANAGER_SOURCES vtkInterpolateToSurfaceFilter.cxx
   REQUIRED_ON_SERVER
)

# let cmake decide which c++-version is needed to satisfy all the features used in the source-code
target_compile_features(InterpolateToSurface PRIVATE cxx_range_for)
