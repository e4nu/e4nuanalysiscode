ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

INCLUDES    := -I../include

CXX       := g++
CXXFLAGS  += -std=c++11 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

OBJECTS1   := FilterData.o Subtraction.o run_e2a_ep_neutrino6_united4_radphot.o Fiducial.o e2a_ep_neutrino6_united4_radphot.o GetCharge_FilterData.o
OBJECTS2   := Subtraction.o run_genie_analysis.o Fiducial.o genie_analysis.o
OBJECTS3   := Subtraction.o run_systematics.o Fiducial.o systematics.o


all: e2a_ep genie_analysis systematics

genie_analysis: $(OBJECTS2)
		$(CXX) -o genie_analysis $(OBJECTS2) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

e2a_ep: $(OBJECTS1)
	$(CXX) -o e2a_ep $(OBJECTS1) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

systematics: $(OBJECTS3)
		$(CXX) -o systematics $(OBJECTS3) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o run_e2a_ep run_genie_analysis run_systematics *~

%.o: %.C
	$(CXX) -c $< -O2 $(CXXFLAGS) $(INCLUDES)
