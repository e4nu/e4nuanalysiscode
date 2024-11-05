#!/bin/bash
path_mc_files="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/"
path_output="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/"
path_xsec="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/"

# Unradiated input 
declare -a InputFiles=("G18_10a_Dipole_LFG_Q2_04_2GeV_eCarbon.gst.root"
    "GEM21_11a_Dipole_LFG_Q2_04_2GeV_eCarbon.gst.root"
)

declare -a OutputFiles=("e4nuanalysis_1p1pip_G18_10a_Dipole_LFG_Q2_04_2GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pip_GEM21_11a_Dipole_LFG_Q2_04_2GeV_eCarbon_NoRad"
)

declare -a XSecFiles=("G18_10a_Dipole_Q2_04_eCarbon.root"
    "GEM21_11a_Dipole_Q2_04_eCarbon.root"
)

## Radiated input
declare -a RadiatedFiles=("G18_10a_Dipole_LFG_Q2_04_2GeV_eCarbon_radiated.gst.root"
    "GEM21_11a_Dipole_LFG_Q2_04_2GeV_eCarbon_radiated.gst.root")

declare -a OutputFilesRadiated=("e4nuanalysis_1p1pip_G18_10a_Dipole_LFG_Q2_04_2GeV_eCarbon_Rad"
    "e4nuanalysis_1p1pip_GEM21_11a_Dipole_LFG_Q2_04_2GeV_eCarbon_Rad"
)

conf_file="ConfFiles/mc_conf/clas6mc_1p1pipanalysis_eC12_2GeV.txt"

cd $E4NUANALYSIS

source e4nu_gpvm_env.sh

number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
  ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[$i]} --output-file ${path_output}${OutputFiles[$i]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[$i]}
  ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[$i]} --output-file ${path_output}${OutputFiles[$i]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[$i]}
done

# Radiated
number_inputs_rad=${#RadiatedFiles[@]}
for (( i=0; i<${number_inputs_rad}; i++ ));
do
    ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${RadiatedFiles[0]} --output-file ${path_output}${OutputFilesRadiated[$i]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --rad-corr true
done

# For systematics

#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --phi-shift 3
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --phi-shift 3
