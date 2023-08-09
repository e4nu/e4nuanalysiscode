#!/bin/bash

declare -a InputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad.gst.root"
                       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10b_00_000/G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad.gst.root"
		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_01_NoFSI_e_on_1000060120_1161MeV_NoRad.gst.root"
                       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/GEM21_11a/GEM21_11a_Q2_01_e_on_1000060120_1161MeV_NoRad.gst.root"
                       "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_1161GeV/apapadop_SuSav2_C12_1161GeV_master.root"
                       "/pnfs/genie/persistent/users/apapadop/e4v_G2018/Exclusive/electrons/C12_1161GeV/apapadop_G2018_C12_1161GeV_master.root")

declare -a OutputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_G18_10a_NoFSI_Q2_01_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_GEM21_11a_Q2_01_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_G18_10a_Q2_01_Afros_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_GEM21_11a_Q2_01_Afros_e_on_1000060120_1161MeV_NoRad")

declare -a XSecFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_01_total_xsec.root"
                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_01_total_xsec.root"
                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_01_total_xsec.root"
		      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/GEM21_11a_00_000/GEM21_11a_00_000_eC12_Q2min_01_total_xsec.root"
                      "/pnfs/genie/persistent/users/apapadop/mySplines/master_Q2_0_1/v3_0_6_WithSuSav2/RootFiles/xsec_carbon12_spline_CCMEC_GTEST19_10b_00_000_Q2_0_1.root"
                      "/pnfs/genie/persistent/users/apapadop/mySplines/master_Q2_0_1/v3_0_6_WithSuSav2/RootFiles/xsec_carbon12_spline_CCMEC_G18_10a_02_11a_Q2_0_1.root")

cd $E4NUANALYSIS
source e4nu_gpvm_env.sh

number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
  ./e4nuanalysis ConfFiles/mc_conf/clas6mc_1pimanalysis_eC12_1GeV.txt ${InputFiles[$i]} ${OutputFiles[$i]} ComputeTrueAccCorr ${XSecFiles[$i]}
  ./e4nuanalysis ConfFiles/mc_conf/clas6mc_1pimanalysis_eC12_1GeV.txt ${InputFiles[$i]} ${OutputFiles[$i]} ComputeTrueRecoAccCorr ${XSecFiles[$i]}
done
