#!/bin/bash
path_mc_files="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/FinalSPSPiAnalysis/"
path_output="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/"
path_xsec="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/FinalSPSPiAnalysis/"

# Unradiated input 
declare -a InputFiles=(    
    "Rarita/Carbon/1GeV/"
    "Rarita/Carbon/2GeV/"
    "Rarita/Carbon/4GeV/"
)

declare -a OutputFiles=(
    "e4nuanalysis_inclusive_GEM21_11a_Rarita_LFG_Q2_01_1GeV_eCarbon_NoRad"
    "e4nuanalysis_inclusive_GEM21_11a_Rarita_LFG_Q2_04_2GeV_eCarbon_NoRad"
    "e4nuanalysis_inclusive_GEM21_11a_Rarita_LFG_Q2_08_4GeV_eCarbon_NoRad"
)

declare -a XSecFiles=(
    "GEM21_11a_Rarita_Q2_01_C.gst.root"
    "GEM21_11a_Rarita_Q2_04_C.gst.root"
    "GEM21_11a_Rarita_Q2_08_C.gst.root"
)

## Radiated input
declare -a RadiatedFiles=(
    "Rarita/Carbon/Radiated/1GeV/"
    "Rarita/Carbon/Radiated/2GeV/"
    "Rarita/Carbon/Radiated/4GeV/"
)

declare -a OutputFilesRadiated=(
    "e4nuanalysis_inclusive_GEM21_11a_Rarita_LFG_Q2_01_1GeV_eCarbon_Rad"
    "e4nuanalysis_inclusive_GEM21_11a_Rarita_LFG_Q2_04_2GeV_eCarbon_Rad"
    "e4nuanalysis_inclusive_GEM21_11a_Rarita_LFG_Q2_08_4GeV_eCarbon_Rad"
)

cd $E4NUANALYSIS

source e4nu_gpvm_env.sh

conf_file="ConfFiles/mc_conf/clas6mc_Inclusive_eC12_1GeV.txt"
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[0]}
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[0]}
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${RadiatedFiles[0]} --output-file ${path_output}${OutputFilesRadiated[0]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --rad-corr true

conf_file="ConfFiles/mc_conf/clas6mc_Inclusive_eC12_2GeV.txt"
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[1]} --output-file ${path_output}${OutputFiles[1]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[1]}
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[1]} --output-file ${path_output}${OutputFiles[1]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[1]}
./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${RadiatedFiles[1]} --output-file ${path_output}${OutputFilesRadiated[1]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[1]} --rad-corr true

conf_file="ConfFiles/mc_conf/clas6mc_Inclusive_eC12_4GeV.txt"
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[2]} --output-file ${path_output}${OutputFiles[2]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[2]}
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[2]} --output-file ${path_output}${OutputFiles[2]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[2]}
#./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${RadiatedFiles[2]} --output-file ${path_output}${OutputFilesRadiated[2]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[2]} --rad-corr true

