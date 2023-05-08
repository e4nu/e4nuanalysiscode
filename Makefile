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
CONF_SRCS := $(wildcard $(SRCDIR)/conf/*.cxx)
PHYSICS_SRCS := $(wildcard $(SRCDIR)/physics/*.cxx)
ANALYSIS_SRCS := $(wildcard $(SRCDIR)/analysis/*.cxx)
APP_SRCS := $(wildcard $(SRCDIR)/app/*.cxx)

# Object files
UTILS_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(UTILS_SRCS))
CONF_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(CONF_SRCS))
PHYSICS_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(PHYSICS_SRCS))
ANALYSIS_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(ANALYSIS_SRCS))
APP_OBJS := $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(APP_SRCS)) $(OBJDIR)/apps/e4nuanalysis.o

all: e4nuanalysis
 
e4nuanalysis: $(UTILS_OBJS) $(CONF_OBJS) $(PHYSICS_OBJS) $(ANALYSIS_OBJS) $(APP_OBJS) $(APP_OBJS)
	@mkdir -p $(@D)
	$(LD) $(LDFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	@mkdir -p $(@D)	
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/apps/e4nuanalysis.o: $(SRCDIR)/apps/e4nuanalysis.cxx
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o e4nuanalysis

.PHONY: test

test: e4nuanalysis
	./e4nuanalysis
