#!/bin/bash

xB=NoxBCut
#xB=xBCut

LocalPath=/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles

#JLabAccount=apapadop@ftp.jlab.org
#JLabPath=/u/home/apapadop/e4nu

JLabAccount=apapadop@uboonegpvm01.fnal.gov
JLabPath=/uboone/app/users/apapadop/e4nu

DataExtension=Data_Final
GENIEExtension=hA2018_Final_RadCorr_LFGM

##### Genie Samples

echo "$LocalPath/1_161/${GENIEExtension}/$xB/12C_1_161_${GENIEExtension}_Plots_FSI_em.root"

scp $JLabAccount:$JLabPath/genie_e2a_ep_C12_1161_neutrino6_united4_radphot_test.root $LocalPath/1_161/${GENIEExtension}/$xB/12C_1_161_${GENIEExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/genie_e2a_ep_4He_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/${GENIEExtension}/$xB/4He_2_261_${GENIEExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/genie_e2a_ep_C12_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/${GENIEExtension}/$xB/12C_2_261_${GENIEExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/genie_e2a_ep_56Fe_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/${GENIEExtension}/$xB/56Fe_2_261_${GENIEExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/genie_e2a_ep_4He_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/${GENIEExtension}/$xB/4He_4_461_${GENIEExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/genie_e2a_ep_C12_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/${GENIEExtension}/$xB/12C_4_461_${GENIEExtension}_Plots_FSI_em.root
scp $JLabAccount:$JLabPath/genie_e2a_ep_56Fe_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/${GENIEExtension}/$xB/56Fe_4_461_${GENIEExtension}_Plots_FSI_em.root

##### Data Samples

#scp $JLabAccount:$JLabPath/data_e2a_ep_C12_1161_neutrino6_united4_radphot_test.root $LocalPath/1_161/$DataExtension/$xB/12C_1_161_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_4He_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/4He_2_261_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_C12_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/12C_2_261_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_56Fe_2261_neutrino6_united4_radphot_test.root $LocalPath/2_261/$DataExtension/$xB/56Fe_2_261_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_4He_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/4He_4_461_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_C12_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/12C_4_461_${DataExtension}_Plots_FSI_em.root
#scp $JLabAccount:$JLabPath/data_e2a_ep_56Fe_4461_neutrino6_united4_radphot_test.root $LocalPath/4_461/$DataExtension/$xB/56Fe_4_461_${DataExtension}_Plots_FSI_em.root
