
cd $E4NUANALYSIS
source e4nu_ifarm_env.csh

# Inclusive analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_Inclusive_eC12_1GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_1161_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_Inclusive_analysis_e_on_1000060120_1161MeV" --analysis-type "IsData" 
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_Inclusive_eC12_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_Inclusive_analysis_e_on_1000060120_2261MeV" --analysis-type "IsData" 
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_Inclusive_eC12_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_C12_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_Inclusive_analysis_e_on_1000060120_4461MeV" --analysis-type "IsData" 

