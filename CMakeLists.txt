# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(inclusive_diffraction CXX C Fortran)
set (CMAKE_CXX_STANDARD 11)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)
git_local_changes(GIT_LOCAL_CHANGES)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Directory where you have build AmplitudeLib libraries 
#link_directories("../amplitudelib_v2/build/lib/")
link_directories("../subnucleondiffraction/build/bin/")
link_directories("../subnucleondiffraction/build/src/")

link_directories("../subnucleondiffraction/build/libColorDipole/src/")
link_directories("../subnucleondiffraction/build/lib/")


# AmplitudeLib headers
#include_directories("../amplitudelib_v2/")
include_directories("../subnucleondiffraction/src/")

add_subdirectory(src)

set_target_properties(inclusive_diffraction PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
