ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

INCLUDES    := -I../include

CXX       := g++
CXXFLAGS  += -std=c++11 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

OBJECTS1   := run_e2a_ep_neutrino6_united4_radphot.o Fiducial.o e2a_ep_neutrino6_united4_radphot.o
OBJECTS2   := run_genie_analysis.o Fiducial.o genie_analysis.o


all: run_e2a_ep run_genie_analysis

genie_analysis: $(OBJECTS2)
		$(CXX) -o genie_analysis $(OBJECTS2) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

run_e2a_ep: $(OBJECTS1)
	$(CXX) -o e2a_ep $(OBJECTS1) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o run_e2a_ep run_genie_analysis *~

%.o: %.C
	$(CXX) -c $< -O2 $(CXXFLAGS) $(INCLUDES)
