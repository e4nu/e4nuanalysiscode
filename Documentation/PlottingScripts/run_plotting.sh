#!/bin/bash

cd $E4NUANALYSIS
source e4nu_gpvm_env.sh

./plote4nuanalysis --mc_location /pnfs/genie/persistent/users/jtenavid/e4nu_files/AnalisedFiles/1p1pim/ --output_location /genie/app/users/jtenavid/Software/e4v/E4NuAnalysis/Source/e4nuanalysiscode/ --output_name test --input_mc_files e4nuanalysis_1p1pimanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad,e4nuanalysis_1p1pimanalysis_G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad --model_names G18_10a,G18_10b --nofsi_file e4nuanalysis_1p1pimanalysis_G18_10a_NoFSI_Q2_01_e_on_1000060120_1161MeV_NoRad --title test --observable_list HadAlphaT