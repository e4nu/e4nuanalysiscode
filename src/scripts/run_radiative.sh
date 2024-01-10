#!/bin/bash 
source /cvmfs/fermilab.opensciencegrid.org/products/genie/bootstrap_genie_ups.sh 
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups 
setup ifdhc v2_6_6
export IFDH_CP_MAXRETRIES=0

setup pdfsets v5_9_1b 
setup gdb v8_1 

cd $CONDOR_DIR_INPUT ;
ifdh cp -D /pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/RadCorr/master-routine_validation_01-eScattering/G18_10a_H_4325MeV_RadFlux.gst.root  $CONDOR_DIR_INPUT ;

git clone https://github.com/e4nu/e4nuanalysiscode.git -b develop/Systematics ;
cd e4nuanalysiscode ; 
source e4nu_gpvm_env.sh ;
make;

./process_radweights --input-gst-file $CONDOR_DIR_INPUT/G18_10a_H_4325MeV_RadFlux.gst.root --output-gst-file $CONDOR_DIR_INPUT/e_on_1000010010_4325MeV_H_G18_10a_radiated.gst.root --true-EBeam 4.325 --nevents 1000 --target 1000010010 --rad-model simple
ifdh cp -D $CONDOR_DIR_INPUT/e_on_1000010010_4325MeV_H_G18_10a_radiated.gst.root /pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/G18_10a/MonoFlux/master-routine_validation_01-eScattering/ ;
