#
# Makefile for Utilities
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

include $(E4NUANALYSIS)/src/make/Make.include

PACKAGE           = E4NuPhysics
DICTIONARY        = _ROOT_DICT_E4NuPhysics
LIBNAME           = libE4NuPhysics
EXTRA_EXT_LIBS    =

all     : rootcint lib
install : install-lib

# Include standard package makefile targets
#
include $(E4NUANALYSIS)/src/make/Make.targets

FORCE:
