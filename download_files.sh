#!/bin/bash

xB=NoxBCut
#xB=xBCut

LocalPath=/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles

#JLabAccount=apapadop@ftp.jlab.org
#JLabPath=/u/home/apapadop/e4nu
DataExtension=Data_Final

GENIEAccount=apapadop@geniegpvm02.fnal.gov
GENIEPath=/genie/app/users/apapadop/e4nu

#GENIEOnlineExtension=G18_10a_02_11a
#GENIELocalExtension=hA2018_Final_RadCorr_LFGM

#GENIEOnlineExtension=SuSav2
#GENIELocalExtension=SuSav2_RadCorr_LFGM

#GENIEOnlineExtension=SuSav2
#GENIELocalExtension=SuSav2_RadCorr_LFGM_NoAccMaps

##### Genie Samples

echo "$LocalPath/1_161/${GENIELocalExtension}/$xB/12C_1_161_${GENIELocalExtension}_Plots_FSI_em.root"

scp $GENIEAccount:$GENIEPath/genie_e2a_ep_C12_1161_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/1_161/${GENIELocalExtension}/$xB/12C_1_161_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/genie_e2a_ep_4He_2261_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/2_261/${GENIELocalExtension}/$xB/4He_2_261_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/genie_e2a_ep_C12_2261_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/2_261/${GENIELocalExtension}/$xB/12C_2_261_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/genie_e2a_ep_56Fe_2261_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/2_261/${GENIELocalExtension}/$xB/56Fe_2_261_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/genie_e2a_ep_4He_4461_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/4_461/${GENIELocalExtension}/$xB/4He_4_461_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/genie_e2a_ep_C12_4461_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/4_461/${GENIELocalExtension}/$xB/12C_4_461_${GENIELocalExtension}_Plots_FSI_em.root
scp $GENIEAccount:$GENIEPath/genie_e2a_ep_56Fe_4461_neutrino6_united4_radphot_test_${GENIEOnlineExtension}.root $LocalPath/4_461/${GENIELocalExtension}/$xB/56Fe_4_461_${GENIELocalExtension}_Plots_FSI_em.root

##### Data Samples

#scp $JLabAccount:$JLabPath/data_e2a_ep_C12_1161_neutrino6_united4_radphot_test.root $LocalPath/1_161/$DataExtension/$xB/12C_1_161_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_4He_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/4He_2_261_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_C12_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/12C_2_261_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_56Fe_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/56Fe_2_261_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_4He_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/4He_4_461_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_C12_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/12C_4_461_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_56Fe_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/56Fe_4_461_${DataExtension}_Plots_FSI_em.root
