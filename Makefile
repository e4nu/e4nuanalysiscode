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

all: e4nuanalysis plot_e4nuanalysis plot_bkg plot_histograms plot_data plot_syst plot_analised_data

e4nuanalysis: $(SRCDIR)/apps/e4nuanalysis.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_e4nuanalysis: $(SRCDIR)/plotting_apps/plot_e4nuanalysis.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_bkg: $(SRCDIR)/plotting_apps/plot_BkgCheck.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_histograms: $(SRCDIR)/plotting_apps/plot_histograms.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_data: $(SRCDIR)/plotting_apps/plot_data.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_syst: $(SRCDIR)/plotting_apps/plot_syst.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@

plot_analised_data: $(SRCDIR)/plotting_apps/plot_analised_data.cxx $(UTILS_OBJS) $(PLOTTING_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS)
		@mkdir -p $(@D)
		$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(OBJDIR)/utils/*.o  $(OBJDIR)/plotting/*.o $(OBJDIR)/physics/*.o $(OBJDIR)/conf/*.o $(OBJDIR)/analysis/*.o $< -o $@


$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/* e4nuanalysis
	rm -rf $(OBJDIR)/* plot_e4nuanalysis
	
.PHONY: test

test: e4nuanalysis
	./e4nuanalysis
