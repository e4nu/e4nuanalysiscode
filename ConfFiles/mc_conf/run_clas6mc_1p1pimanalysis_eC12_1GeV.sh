#!/bin/bash

declare -a InputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/G18_10a_Dipole_LFG_Q2_01_1GeV_eCarbon.gst.root"
		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/G18_10b_Dipole_LFG_Q2_01_1GeV_eCarbon.gst.root"
		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/G18_10c_Dipole_LFG_Q2_01_1GeV_eCarbon.gst.root"
		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/G18_10d_Dipole_LFG_Q2_01_1GeV_eCarbon.gst.root")
#"/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/GEM21_11a_Dipole_LFG_NoFSI_Q2_01_1GeV_eCarbon.gst.root")
#                       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/G18_10a_Rarita_LFG_Q2_01_1GeV_eCarbon.gst.root"

declare -a OutputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/e4nuanalysis_1p1pimanalysis_G18_10a_Dipole_LFG_Q2_01_1GeV_eCarbon_NoRad"
    "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/e4nuanalysis_1p1pimanalysis_G18_10b_Dipole_LFG_Q2_01_1GeV_eCarbon_NoRad"
    "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/e4nuanalysis_1p1pimanalysis_G18_10c_Dipole_LFG_Q2_01_1GeV_eCarbon_NoRad"
    "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/e4nuanalysis_1p1pimanalysis_G18_10d_Dipole_LFG_Q2_01_1GeV_eCarbon_NoRad")
 #   "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/e4nuanalysis_1p1pimanalysis_GEM21_11a_Dipole_LFG_NoFSI_Q2_01_1GeV_eCarbon_NoRad"
#     "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/e4nuanalysis_1p1pimanalysis_G18_10a_Rarita_LFG_Q2_01_1GeV_eCarbon_NoRad"

declare -a XSecFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/G18_10a_Dipole_Q2_01_eCarbon.root"
    "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/G18_10a_Dipole_Q2_01_eCarbon.root"
    "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/G18_10a_Dipole_Q2_01_eCarbon.root"
    "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/G18_10a_Dipole_Q2_01_eCarbon.root"
    "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/G18_10a_Dipole_Q2_01_eCarbon.root")

cd $E4NUANALYSIS
source e4nu_gpvm_env.sh

number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
  ./e4nuanalysis --conf-file ConfFiles/mc_conf/clas6mc_1p1pimanalysis_eC12_1GeV.txt --root-file ${InputFiles[$i]} --output-file ${OutputFiles[$i]} --analysis-type ComputeTrueAccCorr --xsec-file ${XSecFiles[$i]}
  ./e4nuanalysis --conf-file ConfFiles/mc_conf/clas6mc_1p1pimanalysis_eC12_1GeV.txt --root-file ${InputFiles[$i]} --output-file ${OutputFiles[$i]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${XSecFiles[$i]}
done
