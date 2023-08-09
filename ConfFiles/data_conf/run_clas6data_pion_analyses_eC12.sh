#!/bin/bash

cd $E4NUANALYSIS
source e4nu_ifarm_env.csh

# 1p1pim analysis
./e4nuanalysis ConfFiles/data_conf/clas6data_1p1pimanalysis_eC12_1GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_1161_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1p1pimanalysis_e_on_1000060120_1161MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1p1pimanalysis_eC12_2GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_2261_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1p1pimanalysis_e_on_1000060120_2261MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1p1pimanalysis_eC12_4GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_4461_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1p1pimanalysis_e_on_1000060120_4461MeV" "IsData"

# 1p1pip analysis
./e4nuanalysis ConfFiles/data_conf/clas6data_1p1pipanalysis_eC12_1GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_1161_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1p1pipanalysis_e_on_1000060120_1161MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1p1pipanalysis_eC12_2GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_2261_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1p1pipanalysis_e_on_1000060120_2261MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1p1pipanalysis_eC12_4GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_4461_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1p1pipanalysis_e_on_1000060120_4461MeV" "IsData"

# 1pim analysis
./e4nuanalysis ConfFiles/data_conf/clas6data_1pimanalysis_eC12_1GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_1161_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1pimanalysis_e_on_1000060120_1161MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1pimanalysis_eC12_2GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_2261_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1pimanalysis_e_on_1000060120_2261MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1pimanalysis_eC12_4GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_4461_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1pimanalysis_e_on_1000060120_4461MeV" "IsData"

# 1pip analysis
./e4nuanalysis ConfFiles/data_conf/clas6data_1pipanalysis_eC12_1GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_1161_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1pipanalysis_e_on_1000060120_1161MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1pipanalysis_eC12_2GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_2261_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1pipanalysis_e_on_1000060120_2261MeV" "IsData"
./e4nuanalysis ConfFiles/data_conf/clas6data_1pipanalysis_eC12_4GeV.txt "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_4461_neutrino6_united4_radphot_test_100M.root" "/u/home/jtena/Software/e4nuanalysiscode/e4nuanalysis_1pipanalysis_e_on_1000060120_4461MeV" "IsData"
