#
# This is the main e4nu project Makefile.
# Author: Julia Tena Vidal, jtenavidal \at tauex.tau.ac.il
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

include $(E4NUANALYSIS)/src/make/Make.include

BUILD_TARGETS = configuration \
		utilities \
		physics \
		analysis \
		apps \

# ...
# here add main targets to build

INSTALL_TARGETS =  make-install-dirs \
		   copy-install-files

all:     $(BUILD_TARGETS)
install: $(INSTALL_TARGETS)

configuration: FORCE
	@echo " "
	@echo "** Building Configuration..."
	cd ${E4NUANALYSIS}/src/conf && \
	make && \
	cd ${E4NUANALYSIS}
	@echo " Done."

utilities: FORCE
	@echo " "
	@echo "** Building Utilities..."
	cd ${E4NUANALYSIS}/src/utils && \
	make && \
	cd ${E4NUANALYSIS}
physics : FORCE
	@echo " "
	@echo "** Building Physics..."
	cd ${E4NUANALYSIS}/src/physics && \
	make && \
	cd ${E4NUANALYSIS}

analysis : FORCE
	@echo " "
	@echo "** Building Analysis..."
	cd ${E4NUANALYSIS}/src/analysis && \
	make && \
	cd ${E4NUANALYSIS}

apps : FORCE
	@echo " "
	@echo "** Building Apps..."
	cd ${E4NUANALYSIS}/src/apps && \
	make && \
	cd ${E4NUANALYSIS}

clean: FORCE
	@echo " "
	@echo "** Cleaning ..."
	cd ${E4NUANALYSIS}/src/conf &&  rm -f *.o *~ _*  && \
	cd ${E4NUANALYSIS}/src/utils && rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS}/src/apps && rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS}/src/physics &&  rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS}/src/analysis &&  rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS} && rm -f e4nuanalysis && \
	cd ${E4NUANALYSIS}/lib && rm -f *.so &&\
	cd ${E4NUANALYSIS}
FORCE: 	
