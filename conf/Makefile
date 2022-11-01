#
# Makefile for configuration files
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

include $(E4NUANALYSIS)/src/make/Make.include

PACKAGE           = E4nuConf
DICTIONARY        = _ROOT_DICT_E4NuConf
LIBNAME           = libE4nuConf
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(E4NUANALYSIS)/src/make/Make.targets

FORCE:
