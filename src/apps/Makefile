#Makefile for e4nu apps
# Julia Tena Vidal, at jtenavidal \at tauex.ac.il

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

include $(E4NUANALYSIS)/src/make/Make.include

E4NU_LIBS = -L$(E4NUANALYSIS)/lib -lE4NuConf -lE4NuUtils -lE4NuPhysics -lE4NuAnalysisTools
LIBRARIES := $(ROOTLIBS) $(E4NU_LIBS)
LDFLAGS = -g -Wl,--no-as-needed -Wl,--no-undefined
all: e4nuanalysis

e4nuanalysis: FORCE	
	@echo "** Building e4nuanalysis ... "
	g++ -W -Wall -fPIC -Wshadow $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)  -c e4nuanalysis.cxx $(ROOT_INCLUDES) $(E4NUANALYSIS_INCLUDES)
	g++ $(LDFLAGS) e4nuanalysis.o $(LIBRARIES) -o $(E4NUANALYSIS)/e4nuanalysis

purge: FORCE
	rm -f *.o *~ *.d core

clean: FORCE
	rm -f *.o *~ *.d core
	rm -f $(E4NUANALYSIS)/e4nuanalysis

FORCE:

# DO NOT DELETE
