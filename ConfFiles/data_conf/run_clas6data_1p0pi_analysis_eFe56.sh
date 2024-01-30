#!/bin/bash

cd $E4NUANALYSIS
source e4nu_ifarm_env.csh

# 1p0pi analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p0pianalysis_eFe56_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_56Fe_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p0pianalysis_e_on_1000260560_2261MeV" --analysis-type "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1p0pianalysis_eFe56_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_56Fe_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p0pianalysis_e_on_1000260560_4461MeV" --analysis-type "IsData"

