#!/bin/bash

# Set up the UPS products needed to build and use GENIE
source /grid/fermiapp/products/larsoft/setups
setup root v6_12_06a -q e17:debug
setup lhapdf v5_9_1k -q e17:debug
setup log4cpp v1_1_3a -q e17:debug
setup pdfsets v5_9_1b

# Set up much more recent versions of gdb and git
setup gdb v8_1
setup git v2_15_1

echo "Setting E4NU environment variables..."

# Finds the directory where this script is located. This method isn't
# foolproof. See https://stackoverflow.com/a/246128/4081973 if you need
# something more robust for edge cases (e.g., you're calling the script using
# symlinks).
THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export E4NUANALYSIS=$THIS_DIRECTORY

export LD_LIBRARY_PATH="/cvmfs/fermilab.opensciencegrid.org/products/common/db/../prd/curl/v7_64_1/Linux64bit-3-10/lib:/genie/app/users/jtenavid/genie-v3/Generator/lib:/genie/app/users/jtenavid/genie-v3/Reweight/lib:/grid/fermiapp/products/larsoft/git/v2_15_1/Linux64bit+3.10-2.17/lib64:/grid/fermiapp/products/larsoft/gcc/v7_3_0/Linux64bit+3.10-2.17/lib64:/grid/fermiapp/products/larsoft/gcc/v7_3_0/Linux64bit+3.10-2.17/lib:/grid/fermiapp/products/larsoft/log4cpp/v1_1_3a/Linux64bit+3.10-2.17-e17-debug/lib:/grid/fermiapp/products/larsoft/lhapdf/v5_9_1k/Linux64bit+3.10-2.17-e17-debug/lib:/grid/fermiapp/products/larsoft/libxml2/v2_9_5/Linux64bit+3.10-2.17-debug/lib:/grid/fermiapp/products/larsoft/tbb/v2018_2a/Linux64bit+3.10-2.17-e17-debug/lib:/grid/fermiapp/products/larsoft/xrootd/v4_8_0b/Linux64bit+3.10-2.17-e17-debug/lib64:/grid/fermiapp/products/larsoft/mysql_client/v5_5_58a/Linux64bit+3.10-2.17-e17/lib:/grid/fermiapp/products/larsoft/sqlite/v3_20_01_00/Linux64bit+3.10-2.17/lib:/grid/fermiapp/products/larsoft/python/v2_7_14b/Linux64bit+3.10-2.17/lib:/grid/fermiapp/products/larsoft/postgresql/v9_6_6a/Linux64bit+3.10-2.17-p2714b/lib:/grid/fermiapp/products/larsoft/pythia/v6_4_28k/Linux64bit+3.10-2.17-gcc730-debug/lib:/grid/fermiapp/products/larsoft/gsl/v2_4/Linux64bit+3.10-2.17-debug/lib:/grid/fermiapp/products/larsoft/fftw/v3_3_6_pl2/Linux64bit+3.10-2.17-debug/lib:/grid/fermiapp/products/larsoft/clhep/v2_3_4_6/Linux64bit+3.10-2.17-e17-debug/lib:/grid/fermiapp/products/larsoft/root/v6_12_06a/Linux64bit+3.10-2.17-e17-debug/lib:/usr/local/lib:/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vmaster/lib"

export ROOT_INCLUDE_PATH="/grid/fermiapp/products/larsoft/tbb/v2018_2a/Linux64bit+3.10-2.17-e17-debug/include:/grid/fermiapp/products/larsoft/xrootd/v4_8_0b/Linux64bit+3.10-2.17-e17-debug/include:/grid/fermiapp/products/larsoft/mysql_client/v5_5_58a/Linux64bit+3.10-2.17-e17/include:/grid/fermiapp/products/larsoft/postgresql/v9_6_6a/Linux64bit+3.10-2.17-p2714b/include:/grid/fermiapp/products/larsoft/pythia/v6_4_28k/Linux64bit+3.10-2.17-gcc730-debug/include:/grid/fermiapp/products/larsoft/gsl/v2_4/Linux64bit+3.10-2.17-debug/include:/grid/fermiapp/products/larsoft/fftw/v3_3_6_pl2/Linux64bit+3.10-2.17-debug/include:/grid/fermiapp/products/larsoft/clhep/v2_3_4_6/Linux64bit+3.10-2.17-e17-debug/include"
