#!/bin/bash

declare -a InputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_08_e_1000060120_4461MeV_NoRad.gst.root"
                       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10b_00_000/G18_10b_Q2_08_e_on_1000060120_4461MeV_NoRad.gst.root"
		       "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/EventGeneration/G18_10a_00_000/G18_10a_Q2_08_NoFSI_e_1000060120_4461MeV_NoRad.gst.root")

declare -a OutputFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_G18_10a_Q2_08_e_on_1000060120_4461MeV_NoRad"
                        "/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_G18_10b_Q2_08_e_on_1000060120_4461MeV_NoRad"
			"/pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1pim/e4nuanalysis_1pimanalysis_G18_10a_NoFSI_Q2_08_e_on_1000060120_4461MeV_NoRad")

declare -a XSecFiles=("/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_08_total_xsec.root"
                      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_08_total_xsec.root"
		      "/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/Splines/G18_10a_00_000/G18_10a_00_000_eC12_Q2min_08_total_xsec.root")

cd $E4NUANALYSIS
source e4nu_gpvm_env.sh

number_inputs=${#InputFiles[@]}
for (( i=0; i<${number_inputs}; i++ ));
do
    ./e4nuanalysis --conf-file ConfFiles/mc_conf/clas6mc_1pimanalysis_eC12_4GeV.txt --root-file ${InputFiles[$i]} --output-file ${OutputFiles[$i]} --analysis-type ComputeTrueAccCorr --xsec-file ${XSecFiles[$i]}
    ./e4nuanalysis --conf-file ConfFiles/mc_conf/clas6mc_1pimanalysis_eC12_4GeV.txt --root-file ${InputFiles[$i]} --output-file ${OutputFiles[$i]} --analysis-type ComputeTrueRecoAccCorr --xsec-file ${XSecFiles[$i]}
done
