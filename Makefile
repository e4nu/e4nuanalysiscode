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
	cd ${E4NUANALYSIS}/conf && \
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
	cd ${E4NUAnalysis}

apps : FORCE
	@echo " "
	@echo "** Building Apps..."
	cd ${E4NUANALYSIS}/src/apps && \
	make && \
	cd ${E4NUANALYSIS}

make-install-dirs: FORCE
	@echo " "
	@echo "** Creating directory structure for E4NUANALYSIS installation..."
	[ -d ${E4NUANALYSIS_INSTALLATION_PATH} ] || mkdir ${E4NUANALYSIS_INSTALLATION_PATH}
	cd ${E4NUANALYSIS_INSTALLATION_PATH}
	[ -d ${E4NUANALYSIS_BIN_INSTALLATION_PATH} ] || mkdir ${E4NUANALYSIS_BIN_INSTALLATION_PATH}
	[ -d ${E4NUANALYSIS_LIB_INSTALLATION_PATH} ] || mkdir ${E4NUANALYSIS_LIB_INSTALLATION_PATH}
	[ -d ${E4NUANALYSIS_INC_INSTALLATION_PATH} ] || mkdir ${E4NUANALYSIS_INC_INSTALLATION_PATH}
	mkdir ${E4NUANALYSIS_INC_INSTALLATION_PATH}/conf
	mkdir ${E4NUANALYSIS_INC_INSTALLATION_PATH}/src
	mkdir ${E4NUANALYSIS_INC_INSTALLATION_PATH}/src/utils
	mkdir ${E4NUANALYSIS_INC_INSTALLATION_PATH}/src/physics
	mkdir ${E4NUANALYSIS_INC_INSTALLATION_PATH}/src/analysis
	mkdir ${E4NUANALYSIS_INC_INSTALLATION_PATH}/src/apps

copy-install-files: FORCE	
	@echo " "
	@echo "** Copying libraries/binaries/headers to installation location..."
	cp ${E4NUANALYSIS_BIN_PATH}/* ${E4NUANALYSIS_BIN_INSTALLATION_PATH} && \
	cd ${E4NUANALYSIS}/conf && make install && cd .. \
	cd ${E4NUANALYSIS}/src/utils && make install && cd ../.. \
	cd ${E4NUANALYSIS}/src/analysis && make install && cd ../.. \
	cd ${E4NUANALYSIS}/src/physics && make install && cd ../.. \
	cd ${E4NUANALYSIS}/src/apps && make install && cd ../.. \
	cd ${E4NUANALYSIS}

clean: FORCE
	@echo " "
	@echo "** Cleaning ..."
	cd ${E4NUANALYSIS}/conf &&  rm -f *.o *~ _*  && \
	cd ${E4NUANALYSIS}/src/utils && rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS}/src/apps && rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS}/src/analysis && rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS}/src/physics &&  rm -f *.o *~ _* && \
	cd ${E4NUANALYSIS} && rm -f e4nutest && \
	cd ${E4NUANALYSIS}/lib && rm -f *.so &&\
	cd ${E4NUANALYSIS}
FORCE: 	
