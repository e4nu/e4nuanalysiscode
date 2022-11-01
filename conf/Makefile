#
# Makefile for configuration files
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

INCLUDES    := -I../include

CXX       := g++
CXXFLAGS  += -std=c++11 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

PACKAGE           = E4nuConf
DICTIONARY        = _ROOT_DICT_E4nuConf
LIBNAME           = libE4nuConf
EXTRA_EXT_LIBS    =

all     : E4nuConf

OBJ := ParticleI.o TargetI.o AnalysisCuts.o FiducialCuts.o AcceptanceMaps.o
E4nuConf: $(OBJ)
	$(CXX) -o E4nuConf $(OBJ) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)
FORCE:

%.o: %.h	
	$(CXX) -c $< -O2 $(CXXFLAGS) $(INCLUDES)
