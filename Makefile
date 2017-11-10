CXXFLAGS = `gsl-config --cflags` -I ../amplitudelib_v2/ -I ../subnucleondiffraction/src/ -O2
LDFLAGS = `gsl-config --libs` ../amplitudelib_v2/libamplitude.a ../subnucleondiffraction/libColorDipole/libraries/libColorDipole.a -lgfortran

include filelist.m

CXX = g++-7

all: diffraction

diffraction: $(OBJECTS)
	$(CXX) $(CXXFLAGS)  $(OBJECTS) $(LDFLAGS) -o inclusive_diffraction

.cpp.o: src/subnucleon_config.hpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f tools/*.o
	rm -f inclusive_diffraction
