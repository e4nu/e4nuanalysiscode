#!/bin/bash

cd $E4NUANALYSIS
source e4nu_ifarm_env.csh

# 1p1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 4
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 4

# 1p1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 4
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 4

# 1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 4
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 4

# 1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 4
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 4


# Systematics - background multiplicity changes
# Multiplicity 5
# 1p1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 5
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 5

# 1p1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 5
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 5

# 1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 5
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 5

# 1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 5
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 5

# Systematics - background multiplicity changes
# Multiplicity 3
# 1p1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 3
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 3

# 1p1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 3
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 3

# 1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 3
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 3

# 1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 3
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 3


# Multiplicity 0 - No subtraction
# 1p1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 0
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 0

# 1p1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 0
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1p1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1p1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 0

# 1pim analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 0
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pimanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pimanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 0

# 1pip analysis
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_2GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_2261_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_2261MeV" --analysis-type "IsData" --bkg-mult 0
./e4nuanalysis --conf-file ConfFiles/data_conf/clas6data_1pipanalysis_eHe4_4GeV.txt --root-file "/w/hallb-scshelf2102/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_4He_4461_neutrino6_united4_radphot_test_100M.root" --output-file "/work/clas12/jtena/e4nuanalysis_1pipanalysis_e_on_1000020040_4461MeV" --analysis-type "IsData" --bkg-mult 0
