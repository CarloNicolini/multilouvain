cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(multilouvain)

set(LOUVAIN_SRCS
CPMVertexPartition.cpp
GraphHelper.cpp
LinearResolutionParameterVertexPartition.cpp
ModularityVertexPartition.cpp
MutableVertexPartition.cpp
Optimiser.cpp
RBConfigurationVertexPartition.cpp
RBERVertexPartition.cpp
SignificanceVertexPartition.cpp
SurpriseVertexPartition.cpp
DCSurpriseVertexPartition.cpp
igraph_utils.cpp
)

set(LOUVAIN_HDRS
CPMVertexPartition.h
GraphHelper.h
ModularityVertexPartition.h
LinearResolutionParameterVertexPartition.h
MutableVertexPartition.h
Optimiser.h 
RBConfigurationVertexPartition.h
RBERVertexPartition.h
SignificanceVertexPartition.h
SurpriseVertexPartition.h
DCSurpriseVertexPartition.h
igraph_utils.h
)

include_directories(eigen) # for linear algebra operations
include_directories(iglfind)

add_library(LOUVAIN SHARED ${LOUVAIN_SRCS} ${LOUVAIN_HDRS})
target_link_libraries(LOUVAIN ${IGRAPH_LIBRARIES})
if (MATLAB_SUPPORT)
    target_link_libraries(LOUVAIN ${MATLAB_LIBRARIES} ${MATLAB_LIBRARIES_EXTRA})
endif()

add_executable(test_multilouvain test_multilouvain.cpp)
target_link_libraries(test_multilouvain LOUVAIN )

add_executable(test_multilouvain2 test_multilouvain2.cpp)
target_link_libraries(test_multilouvain2 LOUVAIN )

add_executable(run_multilouvain run_multilouvain.cpp)
target_link_libraries(run_multilouvain LOUVAIN )

if(MATLAB_SUPPORT)
    add_mex(multilouvain multilouvain.cpp)
    target_link_libraries(multilouvain LOUVAIN ${MATLAB_LIBRARIES})
endif(MATLAB_SUPPORT)

if(OCTAVE_SUPPORT)
    octave_add_oct(multilouvain SOURCES multilouvain.cpp LINK_LIBRARIES LOUVAIN EXTENSION "mex")
endif(OCTAVE_SUPPORT)
