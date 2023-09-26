#!/bin/bash

declare -a InputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_04_e_on_1000060120_2261MeV_NoRad.gst.root"
                       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10b_00_000/G18_10b_Q2_04_e_on_1000060120_2261MeV_NoRad.gst.root"
		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_04_NoFSI_e_on_1000060120_2261MeV_NoRad.gst.root")
#		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/GEM21_11a/GEM21_11a_Q2_04_e_on_1000060120_2216MeV_NoRad.gst.root"
#                       "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/apapadop_SuSav2_C12_2261GeV_master.root"
#                       "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/apapadop_SuSav2_C12_2261GeV_master.root")

declare -a OutputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10a_Q2_04_e_on_1000060120_2261MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10b_Q2_04_e_on_1000060120_2261MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10a_NoFSI_Q2_04_e_on_1000060120_2261MeV_NoRad")
#                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_GEM21_11a_Q2_04_e_on_1000060120_2261MeV_NoRad"
#                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_GEM21_11a_Q2_04_Afros_e_on_1000060120_2261MeV_NoRad"
#                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10a_Q2_04_Afros_e_on_1000060120_2261MeV_NoRad")

declare -a XSecFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_04_total_xsec.root"
                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_04_total_xsec.root"
                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_04_total_xsec.root")
#                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/GEM21_11a_00_000/GEM21_11a_00_000_eC12_Q2min_04_total_xsec.root"
#                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/GEM21_11a_00_000/GEM21_11a_00_000_eC12_Q2min_04_total_xsec.root"
#                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_04_total_xsec.root")

cd $E4NUANALYSIS
source e4nu_gpvm_env.sh
number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
  ./e4nuanalysis --conf-file ConfFiles/mc_conf/clas6mc_1p1pimanalysis_eC12_2GeV.txt --root-file ${InputFiles[$i]} --output-file ${OutputFiles[$i]} --analysis-type ComputeTrueAccCorr --xsec-file ${XSecFiles[$i]}
  ./e4nuanalysis --conf-file ConfFiles/mc_conf/clas6mc_1p1pimanalysis_eC12_2GeV.txt --root-file ${InputFiles[$i]} --output-file ${OutputFiles[$i]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${XSecFiles[$i]}
done
