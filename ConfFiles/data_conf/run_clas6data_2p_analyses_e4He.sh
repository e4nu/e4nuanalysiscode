
cd $E4NUANALYSIS
source e4nu_ifarm_env.csh

# 2p analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_2panalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_2panalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 5
