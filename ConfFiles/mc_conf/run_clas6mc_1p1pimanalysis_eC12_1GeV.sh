#!/bin/bash

declare -a InputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad.gst.root"
                       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10b_00_000/G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad.gst.root"
		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_04_NoFSI_e_on_1000060120_2261MeV_NoRad.gst.root"
                       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/GEM21_11a/GEM21_11a_Q2_01_e_on_1000060120_1161MeV_NoRad.gst.root")

declare -a OutputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10a_NoFSI_Q2_01_e_on_1000060120_1161MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_GEM21_11a_Q2_01_e_on_1000060120_1161MeV_NoRad")

declare -a XSecFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad"
                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad"
		      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_G18_10a_NoFSI_Q2_01_e_on_1000060120_1161MeV_NoRad"
                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/e4nuanalysis_1p1pimanalysis_GEM21_11a_Q2_01_e_on_1000060120_1161MeV_NoRad")

cd $E4NUANALYSIS
number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
  ./e4nuanalysis ConfFiles/mc_conf/clas6mc_1p1pimanalysis_eC12_1GeV.txt ${InputFiles[$i]} ${OutputFiles[$i]} ComputeTrueAccCorr ${XSecFiles[$i]}
  ./e4nuanalysis ConfFiles/mc_conf/clas6mc_1p1pimanalysis_eC12_1GeV.txt ${InputFiles[$i]} ${OutputFiles[$i]} ComputeTrueRecoAccCorr ${XSecFiles[$i]}
done
