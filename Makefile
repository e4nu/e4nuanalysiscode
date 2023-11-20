# Directories
SRCDIR := src
OBJDIR := obj
BINDIR := bin

# ROOT                                                                                                                                                                                               
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

# Compiler options
CXX := g++
CXXFLAGS := -Wall -O2 -I$(SRCDIR) $(ROOTCFLAGS)

# Linker options
LD := g++
LDFLAGS := $(ROOTLIBS) -lstdc++ -Wno-undef

# Source files
UTILS_SRCS := $(wildcard $(SRCDIR)/utils/*.cxx)
PLOTTING_SRCS := $(wildcard $(SRCDIR)/plotting/*.cxx)
CONF_SRCS := $(wildcard $(SRCDIR)/conf/*.cxx)
PHYSICS_SRCS := $(wildcard $(SRCDIR)/physics/*.cxx)
ANALYSIS_SRCS := $(wildcard $(SRCDIR)/analysis/*.cxx)
APP_SRCS := $(wildcard $(SRCDIR)/app/*.cxx)

# Object files
UTILS_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(UTILS_SRCS))
PLOTTING_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(PLOTTING_SRCS))
CONF_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(CONF_SRCS))
PHYSICS_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(PHYSICS_SRCS))
ANALYSIS_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(ANALYSIS_SRCS))
APP_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(APP_SRCS)) 

all: e4nuanalysis plot_e4nuanalysis plot_radiative radiate_flux plot_bkg
 
e4nuanalysis: $(SRCDIR)/apps/e4nuanalysis.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_e4nuanalysis: $(SRCDIR)/apps/plote4nuanalysis.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_radiative: $(SRCDIR)/apps/plot_radiative.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

radiate_flux: $(SRCDIR)/apps/radiate_flux.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_bkg: $(SRCDIR)/apps/plot_BkgCheck.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	@mkdir -p $(@D)	
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/* e4nuanalysis
	rm -rf $(OBJDIR)/* plot_e4nuanalysis
	rm -rf $(OBJDIR)/* plot_radiative
	rm -rf $(OBJDIR)/* plot_bkg

.PHONY: test

test: e4nuanalysis
	./e4nuanalysis
