configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	inclusivediffraction
	STATIC
	gitsha1.h
	inclusive_diffraction.cpp
	qqg.cpp
	dummy.f
)

add_executable(
	inclusive_diffraction	
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
	main.cpp
	dummy.f
)
# Fortran file is here as a hack, so that CMake includes -lgfortran


include_directories(.)

find_package(GSL REQUIRED)    #
target_link_libraries(
	inclusive_diffraction
	PRIVATE
	inclusivediffraction
	amplitude
	ddis
	ColorDipole
	mzipsat
	GSL::gsl
	GSL::gslcblas
)


