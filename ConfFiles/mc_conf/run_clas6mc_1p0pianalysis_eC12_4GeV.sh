#!/bin/bash
path_mc_files="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/FinalSPSPiAnalysis/"
path_output="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/"
path_xsec="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/FinalSPSPiAnalysis/"

# Unradiated input 
declare -a InputFiles=(    
    "Rarita/Carbon/4GeV/"
    "Rarita/Carbon/NoFsi/4GeV/"
    "Dipole/Carbon/4GeV/"
)

declare -a OutputFiles=(
    "e4nuanalysis_1p0pi_GEM21_11a_Rarita_LFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p0pi_GEM21_11a_Rarita_LFG_Q2_08_4GeV_eCarbon_NoRad_NoFSI"
    "e4nuanalysis_1p0pi_GEM21_11a_Dipole_LFG_Q2_08_4GeV_eCarbon_NoRad"
)

declare -a XSecFiles=(
    "GEM21_11a_Rarita_Q2_08_C.gst.root"
    "GEM21_11a_Rarita_Q2_08_C.gst.root"
    "GEM21_11a_Dipole_Q2_08_C.gst.root"
)

## Radiated input
declare -a RadiatedFiles=("Rarita/Carbon/Radiated/4GeV/")

declare -a OutputFilesRadiated=("e4nuanalysis_1p0pi_GEM21_11a_Rarita_LFG_Q2_08_4GeV_eCarbon_Rad")
                                 
conf_file="ConfFiles/mc_conf/clas6mc_1p0pianalysis_eC12_4GeV.txt"

cd $E4NUANALYSIS

source e4nu_gpvm_env.sh

number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
  ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[$i]} --output-file ${path_output}${OutputFiles[$i]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[$i]}
  ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[$i]} --output-file ${path_output}${OutputFiles[$i]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[$i]}
done

# Closure Test 

  ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ClosureTest --xsec-file ${path_xsec}${XSecFiles[0]}

# Radiated
./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${RadiatedFiles[0]} --output-file ${path_output}${OutputFilesRadiated[0]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --rad-corr true

# For systematics

./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --phi-shift 3
./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --phi-shift 3


