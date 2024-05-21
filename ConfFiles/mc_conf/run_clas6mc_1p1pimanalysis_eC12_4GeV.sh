#!/bin/bash
path_mc_files="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/"
path_output="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024AnalisedFiles/"
path_xsec="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Splines/"

declare -a InputFiles=("G18_10a_Dipole_LFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11a_Dipole_CFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11a_Dipole_LFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11a_Dipole_RFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11b_Dipole_LFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11c_Dipole_CFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11c_Dipole_LFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11d_Dipole_CFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11d_Dipole_LFG_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11a_Dipole_LFG_NoFSI_Q2_08_4GeV_eCarbon.gst.root"
    "GEM21_11a_Rarita_LFG_Q2_08_4GeV_eCarbon.gst.root")

declare -a OutputFiles=("e4nuanalysis_1p1pim_G18_10a_Dipole_LFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11a_Dipole_CFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11a_Dipole_LFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11a_Dipole_RFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11b_Dipole_LFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11c_Dipole_CFG_Q2_08_4GeV_eCarbon.gst.root"
    "e4nuanalysis_1p1pim_GEM21_11c_Dipole_LFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11d_Dipole_CFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11d_Dipole_LFG_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11a_Dipole_LFG_NoFSI_Q2_08_4GeV_eCarbon_NoRad"
    "e4nuanalysis_1p1pim_GEM21_11a_Rarita_LFG_Q2_08_4GeV_eCarbon_NoRad"
)

declare -a XSecFiles=("G18_10a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Dipole_Q2_08_eCarbon.root"
		      "GEM21_11a_Rarita_Q2_08_eCarbon.root"
)

conf_file="ConfFiles/mc_conf/clas6mc_1p1pimanalysis_eC12_4GeV.txt"

cd $E4NUANALYSIS

source e4nu_gpvm_env.sh

number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
  ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[$i]} --output-file ${path_output}${OutputFiles[$i]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[$i]}
  ./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[$i]} --output-file ${path_output}${OutputFiles[$i]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[$i]}
done

./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --phi-shift 3
./e4nuanalysis --conf-file ${conf_file} --root-file ${path_mc_files}${InputFiles[0]} --output-file ${path_output}${OutputFiles[0]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${path_xsec}${XSecFiles[0]} --phi-shift 3
