#
#
#

include ../Makefile.params

PLUGINS = test.so xyz.so skewstart.so crystalfcc.so tempscaling.so berendsen.so cellscaling.so testvisualizer.so nullintegrator.so euler.so verlet.so velocityverlet.so beeman.so nullpotential.so nullpairpotential.so lennardjones.so fastlj.so morse.so harmonic.so constantforce.so pressure.so gdr.so angdist.so cordnum.so cordnumfunc.so vacf.so buckingham.so lpmd.so ewald.so ewald-2.so povray.so povray2.so

all: $(PLUGINS)

%.so: %.cc %.h
	$(CXX) $(CXXFLAGS) $(PLUGINFLAGS) -I../lib -shared $< -o $@

clean:
	-rm -fr $(PLUGINS) *~ 

