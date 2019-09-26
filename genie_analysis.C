
#define GENIE_ANALYSIS_C

#include "genie_analysis.h"
#include "Constants.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <exception>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH3.h>
#include <TGraph.h>

using namespace std;

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en

void genie_analysis::Loop()
{

  std::map<std::string,double>bind_en;
  std::map<std::string,double>target_mass;
  std::map<std::string,double>residual_target_mass;
  std::map<std::string, double> Ecal_offset; //that might not be necessary for simulation data


  target_name = ftarget;   //std string for target name
  en_beam["1161"]=1.161;
  en_beam["2261"]=2.261;
  en_beam["4461"]=4.461;

  en_beam_Ecal["1161"]=1.161;
  en_beam_Ecal["2261"]=2.261;
  en_beam_Ecal["4461"]=4.461;

  en_beam_Eqe["1161"]=1.161;
  en_beam_Eqe["2261"]=2.261;
  en_beam_Eqe["4461"]=4.461;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //nentries =8000000;

  double reso_p = 0.01; // smearing for the proton
  double reso_e = 0.005; // smearing for the electrons
//  double reso_pipl = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
//  double reso_pimi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
  double reso_pi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)

  double N_prot1 = 0, N_prot2 = 0,N_prot_both = 0;
  const int n_slice=3;
  const double pperp_min[n_slice]={0.,0.2,0.4};
  const double pperp_max[n_slice]={0.2,0.4,10.};
  TVector3 V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;
  TVector3 V3_phot_angles;

  TString E_acc_file;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV  Configuration parameters and cuts
  {
      E_acc_file="1_161";
  }


  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV  Configuration parameters and cuts
  {
      E_acc_file="2_261";
  }

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5)  //4.4 GeV  Configuration parameters and cuts
  {
      E_acc_file="4_461";
  }
  //Further constants for binding energies and target masses
  Ecal_offset["3He"]=0.004;
  Ecal_offset["4He"]=0.005;
  Ecal_offset["C12"]=0.005;
  Ecal_offset["56Fe"]=0.011;

  bind_en["3He"] = He3_bind_en-D2_bind_en + Ecal_offset["3He"]; //the offset is used to shift the peak to be at 0
  bind_en["4He"] = He4_bind_en-H3_bind_en + Ecal_offset["4He"];
  bind_en["C12"] = C12_bind_en-B_bind_en  + Ecal_offset["C12"];
  bind_en["56Fe"]= Fe_bind_en-Mn_bind_en  + Ecal_offset["56Fe"];
  bind_en["CH2"] = C12_bind_en-B_bind_en;

  target_mass["3He"] = 2*m_prot+m_neut-He3_bind_en;
  target_mass["4He"] = 2*m_prot+2*m_neut-He4_bind_en;
  target_mass["C12"] = 6*m_prot+6*m_neut-C12_bind_en;
  target_mass["56Fe"]= 26*m_prot+30*m_neut-Fe_bind_en;
  target_mass["CH2"] = 6*m_prot+6*m_neut-C12_bind_en;

  residual_target_mass["3He"] = m_prot+m_neut-D2_bind_en;
  residual_target_mass["4He"] = m_prot+2*m_neut-H3_bind_en;
  residual_target_mass["C12"] = 5*m_prot+6*m_neut-B_bind_en;
  residual_target_mass["56Fe"]= 25*m_prot+30*m_neut-Mn_bind_en;
  residual_target_mass["CH2"] = 25*m_prot+30*m_neut-Mn_bind_en;

  //Definition of Histograms
  TH1F *h1_Etot_p_bkgd_slice[n_slice], *h1_Erec_p_bkgd_slice[n_slice];
  TH1F *h1_Etot_Npi0[n_slice],*h1_Erec_Npi0[n_slice];
  TH1F *h1_Etot_Npi1[n_slice],*h1_Erec_Npi1[n_slice];

  TH1F *h1_Erec_bkgd_pipl_pimi_new_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact[n_slice];
  TH1F *h1_Etot_bkgd_pipl_pimi_fact_pipl[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pimi[n_slice];

  TH1F *h1_Etot_bkgd_1p2pi[n_slice],      *h1_Erec_bkgd_1p2pi[n_slice];
  TH1F *h1_Etot_bkgd_1p2pi_1p0pi[n_slice],*h1_Erec_bkgd_1p2pi_1p0pi[n_slice];
  TH1F *h1_Etot_bkgd_1p3pi[n_slice],      *h1_Erec_bkgd_1p3pi[n_slice];

  TH1F *h1_Etot_p_bkgd_slice_2p2pi[n_slice],        *h1_Erec_p_bkgd_slice_2p2pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[n_slice];

  TH1F *h1_Etot_piplpimi_subtruct_fact[n_slice],*h1_Erec_piplpimi_subtruct_fact[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub[n_slice],*h1_Erec_p_bkgd_slice_sub[n_slice];
  TH1F *h1_Etot_3pto1p_slice[n_slice],*h1_Erec_3pto1p_slice[n_slice];
  TH1F *h1_Etot_3pto2p_slice[n_slice],*h1_Erec_3pto2p_slice[n_slice];
  TH1F *h1_Etot_3p1pi_slice[n_slice], *h1_Erec_3p1pi_slice[n_slice];
  TH1F *h1_Etot_4pto3p_slice[n_slice],*h1_Erec_4pto3p_slice[n_slice];
  TH1F *h1_Etot_4pto1p_slice[n_slice],*h1_Erec_4pto1p_slice[n_slice];
  TH1F *h1_Etot_4pto2p_slice[n_slice],*h1_Erec_4pto2p_slice[n_slice];
  TH1F *h1_Etot_43pto1p_slice[n_slice],*h1_Erec_43pto1p_slice[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub1p2pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p2pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_1p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_2p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_2p[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub1p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p2pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub1p3pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p3pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p2pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub3p1pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub3p1pi_0pi[n_slice];

  TH1F *h1_Etot_p_bkgd_slice_sub32[n_slice],*h1_Erec_p_bkgd_slice_sub32[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub31[n_slice],*h1_Erec_p_bkgd_slice_sub31[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub43[n_slice],*h1_Erec_p_bkgd_slice_sub43[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub41[n_slice],*h1_Erec_p_bkgd_slice_sub41[n_slice];
  TH1F *h1_Erec_p_bkgd_slice_sub42[n_slice],*h1_Etot_p_bkgd_slice_sub42[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub431[n_slice],*h1_Erec_p_bkgd_slice_sub431[n_slice];
  TH2F *h2_N_pi_phot[20];


  gRandom = new TRandom3();
  gRandom->SetSeed(10);

  TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
  TLorentzVector V4_target(0,0,0,target_mass[ftarget]);
  //Acceptance Maps
  TString WhichMap = "e2a_maps";
  TFile* file_acceptance = TFile::Open(WhichMap+"/"+WhichMap+"_12C_E_"+E_acc_file+".root");
  TFile* file_acceptance_p = TFile::Open(WhichMap+"/"+WhichMap+"_12C_E_"+E_acc_file+"_p.root");
  TFile* file_acceptance_pip = TFile::Open(WhichMap+"/"+WhichMap+"_12C_E_"+E_acc_file+"_pip.root");


  //Output file definition
  TFile *file_out = new TFile(Form("e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");

  //initialize Fiducial functions for EC limits
  fiducialcut->InitEClimits();
  std::cout << " Test InitEClimits GenieLoop " << fiducialcut->up_lim1_ec->Eval(60) << std::endl;




  //Definition and initialization of Histograms
  TH1F *h1_el_Mott_crosssec = new TH1F("h1_el_Mott_crosssec","",200,0.,0.01);
  TH1F *h1_Wvar = new TH1F("h1_Wvar","",400,0,3);
  TH1F *h1_Wvar_weight = new TH1F("h1_Wvar_weight","",400,0,3);
  TH1F *h1_xbjk = new TH1F("h1_xbjk","",400,0,3);
  TH1F *h1_xbjk_weight = new TH1F("h1_xbjk_weight","",400,0,3);
  TH1F *h1_Q2 = new TH1F("h1_Q2","",400,0,6);
  TH1F *h1_Q2_weight = new TH1F("h1_Q2_weight","",400,0,6);
  TH1F *h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
  TH1F *h1_Nprot=new TH1F("h1_Nprot","",10,0,5);
  TH1F *h1_Nphot=new TH1F("h1_Nphot","",10,0,5);
  TH1F *h1_Npiphot=new TH1F("h1_Npiphot","",10,0,5);
  TH1F *h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,0,5);
  TH1F *h1_Npi=new TH1F("h1_Npi","",10,0,5);
  TH1F *h1_Npipl=new TH1F("h1_Npipl","",10,0,5);
  TH1F *h1_Npimi=new TH1F("h1_Npimi","",10,0,5);
  TH1F *h1_el_mom = new TH1F("h1_el_mom","",100,1.2,6);
  TH1F *h1_el_mom_corr = new TH1F("h1_el_mom_corr","",100,1.2,6);
  TH1F *h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
  TH1F *h1_prot_mom = new TH1F("h1_prot_mom","",300,0,3);
  TH1F *h1_prot_mom_ratio = new TH1F("h1_prot_mom_ratio","",50,0.97,1.2);

  //Binning for fractional feed down energy histograms
  int N_qe;
  double *x_qe;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
    N_qe=109;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=64;i++)  x_qe[i]=-1+i*0.015;
    for (int i=0;i<=44;i++)  x_qe[i+65]=-0.04+(i+1)*0.01;
  }

  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
    N_qe=109;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=64;i++)  x_qe[i]=-1+i*0.015;
    for (int i=0;i<=44;i++)  x_qe[i+65]=-0.04+(i+1)*0.01;
  }

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
    N_qe=82;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=29;i++)  x_qe[i]=-1+i*0.03;
    for (int i=0;i<=52;i++)  x_qe[i+30]=-0.13+(i+1)*0.01;
  }

 //Definitions of further Histograms
  TH1F *h1_E_rec_1pi_weight_frac_feed=new TH1F("h1_E_rec_1pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_2pi_weight_frac_feed=new TH1F("h1_E_rec_2pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_3pi_weight_frac_feed=new TH1F("h1_E_rec_3pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pi_weight_frac_feed=new TH1F("h1_E_rec_4pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_0pi_frac_feed=new TH1F("h1_E_rec_0pi_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_tot_cut2_fracfeed = new TH1F("h1_E_tot_cut2_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_cut2_new_fracfeed = new TH1F("h1_E_rec_cut2_new_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_p_bkgd_fracfeed = new TH1F("h1_E_tot_p_bkgd_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_p_bkgd_fracfeed = new TH1F("h1_E_rec_p_bkgd_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_3pto2p_fracfeed = new TH1F("h1_E_tot_3pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_3pto2p_fracfeed = new TH1F("h1_E_rec_3pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_3pto1p_fracfeed = new TH1F("h1_E_tot_3pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_3pto1p_fracfeed = new TH1F("h1_E_rec_3pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_4pto3p_fracfeed = new TH1F("h1_E_tot_4pto3p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pto3p_fracfeed = new TH1F("h1_E_rec_4pto3p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_43pto1p_fracfeed = new TH1F("h1_E_tot_43pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_43pto1p_fracfeed = new TH1F("h1_E_rec_43pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_4pto2p_fracfeed = new TH1F("h1_E_tot_4pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pto2p_fracfeed = new TH1F("h1_E_rec_4pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_4pto1p_fracfeed = new TH1F("h1_E_tot_4pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pto1p_fracfeed = new TH1F("h1_E_rec_4pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_1p2pi_fracfeed = new TH1F("h1_E_tot_1p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_1p2pi_fracfeed = new TH1F("h1_E_rec_1p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_1p3pi_fracfeed = new TH1F("h1_E_tot_1p3pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_1p3pi_fracfeed = new TH1F("h1_E_rec_1p3pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p2pi_fracfeed = new TH1F("h1_E_tot_2p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p2pi_fracfeed = new TH1F("h1_E_rec_2p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_3p1pi_fracfeed = new TH1F("h1_E_tot_3p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_3p1pi_fracfeed = new TH1F("h1_E_rec_3p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_tot_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_rec_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_undetfactor_fracfeed = new TH1F("h1_E_rec_undetfactor_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_undetfactor_fracfeed = new TH1F("h1_E_tot_undetfactor_fracfeed","",N_qe,x_qe);

  TH1F *h1_theta0=new TH1F("h1_theta0","",300,0,180);
  TH2F *h2_Ecal_Eqe=new TH2F("h2_Ecal_Eqe","",600,0,5,600,0,5);
  TH2F *h2_EqeEcalratio_Eqe=new TH2F("h2_EqeEcalratio_Eqe","",600,0,5,300,0,2);
  TH2F *h2_EqeEcaldiff_Eqe=new TH2F("h2_EqeEcaldiff_Eqe","",600,0,5,300,-3,3);
  TH2F *h2_N_prot_pi=new TH2F("h2_N_prot_pi","",10,0,5,10,0,5);
  TH2F *h2_N_prot_pi_phot=new TH2F("h2_N_prot_pi_phot","",10,0,5,10,0,5);
  TH2F *h2_N_prot_pi_phot_nonrad=new TH2F("h2_N_prot_pi_phot_nonrad","",10,0,5,10,0,5);
//  TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,0,180);
  TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,10,60); // apapadop
  TH2F *h2_el_mom_diff = new TH2F("h2_el_mom_diff","",500,0.,1.,500,-0.1,0.1);
  TH2F *h2_Q2_nu = new TH2F("h2_Q2_nu","",200,0,3.5,200,0,5);
  TH2F *h2_Q2_nu_weight = new TH2F("h2_Q2_nu_weight","",200,0,3.5,200,0,5);
  TH2F *h2_Q2_xbjk_weight = new TH2F("h2_Q2_xbjk_weight","",200,0,3,200,0,5);
  TH2F *h2_Q2_W=new TH2F("h2_Q2_W","",200,0,3,200,0,5);
  TH2F *h2_xB_W=new TH2F("h2_xB_W","",200,0,3,200,0,3);
  TH2F *h2_Q2_W_weight=new TH2F("h2_Q2_W_weight","",200,0,3,200,0,5);
  TH2F *h2_el_pcorr_puncorr = new TH2F("h2_el_pcorr_puncorr","",100,0,1,100,0,3);
  TH2F *h2_Erec_pperp = new TH2F("h2_Erec_pperp","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_newcut2 = new TH2F("h2_Erec_pperp_newcut2","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_cut3 = new TH2F("h2_Erec_pperp_cut3","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p = new TH2F("h2_Erec_pperp_2p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_321p = new TH2F("h2_Erec_pperp_321p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_31p = new TH2F("h2_Erec_pperp_31p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_4321p = new TH2F("h2_Erec_pperp_4321p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_431p = new TH2F("h2_Erec_pperp_431p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_421p = new TH2F("h2_Erec_pperp_421p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_41p = new TH2F("h2_Erec_pperp_41p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p1pi = new TH2F("h2_Erec_pperp_1p1pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p2pi_1p0pi = new TH2F("h2_Erec_pperp_1p2pi_1p0pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p2pi_1p1pi = new TH2F("h2_Erec_pperp_1p2pi_1p1pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p1pi_2p0pi = new TH2F("h2_Erec_pperp_2p1pi_2p0pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p1pi_1p1pi = new TH2F("h2_Erec_pperp_2p1pi_1p1pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p1pi_1p0pi = new TH2F("h2_Erec_pperp_2p1pi_1p0pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p3pi = new TH2F("h2_Erec_pperp_1p3pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p2pi = new TH2F("h2_Erec_pperp_2p2pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_3p1pi = new TH2F("h2_Erec_pperp_3p1pi","",400,0,1,400,0,6.);
  TH2F *h2_pperp_W=new TH2F("h2_pperp_W","",200,0,3,200,0,2);

  TH2F *h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);

  //Binning for energy reconstruction histograms
  int n_bins;
  double *x_values;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
    n_bins=38;
    x_values=new double[n_bins+1];
    for (int i=0;i<=17;i++) x_values[i]=0.4+i*0.04;
    for (int i=0;i<=20;i++) x_values[i+18]=1.08+(i+1)*0.02;
  }

  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
    n_bins=54;
    x_values=new double[n_bins+1];
    for (int i=0;i<=23;i++) x_values[i]=i*0.09;
    for (int i=0;i<=30;i++) x_values[i+24]=2.07+(i+1)*0.03;
  }

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
    n_bins=38;
    x_values=new double[n_bins+1];
    for (int i=0;i<=21;i++)  x_values[i]=i*0.2;
    for (int i=0;i<=16;i++)  x_values[i+22]=4.2+(i+1)*0.05;
  }

 //Definitions of further Histograms
  TH1F *h1_E_rec_2p_det = new TH1F("h1_E_rec_2p_det","",n_bins,x_values);
  TH1F *h1_E_tot_2p_det = new TH1F("h1_E_tot_2p_det","",n_bins,x_values);
  TH1F *h1_E_tot_p_bkgd = new TH1F("h1_E_tot_p_bkgd","",n_bins,x_values);
  TH1F *h1_E_rec_p_bkgd = new TH1F("h1_E_rec_p_bkgd","",n_bins,x_values);
  TH1F *h1_E_tot_3pto1p = new TH1F("h1_E_tot_3pto1p","",n_bins,x_values);
  TH1F *h1_E_rec_3pto1p = new TH1F("h1_E_rec_3pto1p","",n_bins,x_values);
  TH1F *h1_E_tot_43pto1p = new TH1F("h1_E_tot_43pto1p","",n_bins,x_values);
  TH1F *h1_E_rec_43pto1p = new TH1F("h1_E_rec_43pto1p","",n_bins,x_values);
  TH1F *h1_E_tot_3pto2p =new TH1F("h1_E_tot_3pto2p","",n_bins,x_values);
  TH1F *h1_E_rec_3pto2p = new TH1F("h1_E_rec_3pto2p","",n_bins,x_values);
  TH1F *h1_E_tot_4pto1p =new TH1F("h1_E_tot_4pto1p","",n_bins,x_values);
  TH1F *h1_E_rec_4pto1p = new TH1F("h1_E_rec_4pto1p","",n_bins,x_values);
  TH1F *h1_E_tot_4pto3p =new TH1F("h1_E_tot_4pto3p","",n_bins,x_values);
  TH1F *h1_E_rec_4pto3p = new TH1F("h1_E_rec_4pto3p","",n_bins,x_values);
  TH1F *h1_E_tot_4pto2p =new TH1F("h1_E_tot_4pto2p","",n_bins,x_values);
  TH1F *h1_E_rec_4pto2p = new TH1F("h1_E_rec_4pto2p","",n_bins,x_values);
  TH1F *h1_E_rec = new TH1F("h1_E_rec","",n_bins,x_values);
  TH1F *h1_E_rec_0pi = new TH1F("h1_E_rec_0pi","",n_bins,x_values);
  TH1F *h1_E_rec_1pi = new TH1F("h1_E_rec_1pi","",n_bins,x_values);
  TH1F *h1_E_rec_1pi_weight = new TH1F("h1_E_rec_1pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_2pi_weight = new TH1F("h1_E_rec_2pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_3pi_weight = new TH1F("h1_E_rec_3pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_4pi_weight = new TH1F("h1_E_rec_4pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_20pi = new TH1F("h1_E_rec_20pi","",n_bins,x_values);
  TH1F *h1_E_rec_21pi = new TH1F("h1_E_rec_21pi","",n_bins,x_values);
  TH1F *h1_E_rec_30pi = new TH1F("h1_E_rec_30pi","",n_bins,x_values);
  TH1F *h1_E_rec_310pi = new TH1F("h1_E_rec_310pi","",n_bins,x_values);
  TH1F *h1_E_rec_320pi = new TH1F("h1_E_rec_320pi","",n_bins,x_values);
  TH1F *h1_E_rec_3210pi = new TH1F("h1_E_rec_3210pi","",n_bins,x_values);
  TH1F *h1_E_rec_40pi = new TH1F("h1_E_rec_40pi","",n_bins,x_values);
  TH1F *h1_E_rec_410pi = new TH1F("h1_E_rec_410pi","",n_bins,x_values);
  TH1F *h1_E_rec_420pi = new TH1F("h1_E_rec_420pi","",n_bins,x_values);
  TH1F *h1_E_rec_4210pi = new TH1F("h1_E_rec_4210pi","",n_bins,x_values);
  TH1F *h1_E_rec_430pi = new TH1F("h1_E_rec_430pi","",n_bins,x_values);
  TH1F *h1_E_rec_4310pi = new TH1F("h1_E_rec_4310pi","",n_bins,x_values);
  TH1F *h1_E_rec_4320pi = new TH1F("h1_E_rec_4320pi","",n_bins,x_values);
  TH1F *h1_E_rec_43210pi = new TH1F("h1_E_rec_43210pi","",n_bins,x_values);
  TH1F *h1_E_rec_1prot  = new TH1F("h1_E_rec_1prot","",n_bins,x_values);
  TH1F *h1_E_tot_1prot  = new TH1F("h1_E_tot_1prot","",n_bins,x_values);
  TH1F *h1_E_rec_cutpi1_piplpimi = new TH1F("h1_E_rec_cutpi1_piplpimi","",n_bins,x_values);
  TH1F *h1_E_tot_cutpi1_piplpimi = new TH1F("h1_E_tot_cutpi1_piplpimi","",n_bins,x_values);
  TH1F *h1_E_tot = new TH1F("h1_E_tot","",n_bins,x_values);
  TH1F *h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
  TH1F *h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
  TH1F *h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
	TH1F *h1_E_rec_undetfactor  = new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
  TH1F *h1_E_tot_undetfactor  = new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
  TH1F *h1_E_tot_1p2pi  = new TH1F("h1_E_tot_1p2pi","",n_bins,x_values);
  TH1F *h1_E_rec_1p2pi  = new TH1F("h1_E_rec_1p2pi","",n_bins,x_values);
  TH1F *h1_E_tot_1p3pi  = new TH1F("h1_E_tot_1p3pi","",n_bins,x_values);
  TH1F *h1_E_rec_1p3pi  = new TH1F("h1_E_rec_1p3pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p2pi  = new TH1F("h1_E_tot_2p2pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p2pi  = new TH1F("h1_E_rec_2p2pi","",n_bins,x_values);
  TH1F *h1_E_tot_3p1pi  = new TH1F("h1_E_tot_3p1pi","",n_bins,x_values);
  TH1F *h1_E_rec_3p1pi  = new TH1F("h1_E_rec_3p1pi","",n_bins,x_values);
  TH1F *h1_E_tot_1p2pi_1p0pi  = new TH1F("h1_E_tot_1p2pi_1p0pi","",n_bins,x_values);
  TH1F *h1_E_rec_1p2pi_1p0pi  = new TH1F("h1_E_rec_1p2pi_1p0pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p1pi_2p0pi  = new TH1F("h1_E_tot_2p1pi_2p0pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p1pi_2p0pi  = new TH1F("h1_E_rec_2p1pi_2p0pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p1pi_1p1pi  = new TH1F("h1_E_tot_2p1pi_1p1pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p1pi_1p1pi  = new TH1F("h1_E_rec_2p1pi_1p1pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p1pi_1p0pi  = new TH1F("h1_E_tot_2p1pi_1p0pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p1pi_1p0pi  = new TH1F("h1_E_rec_2p1pi_1p0pi","",n_bins,x_values);

  //Defintions of Histogram for each slice
  for(int h = 0; h < n_slice; h++){
   h1_Erec_p_bkgd_slice[h]= new TH1F(Form("h1_Erec_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice[h]= new TH1F(Form("h1_Etot_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3pto1p_slice[h]= new TH1F(Form("h1_Erec_3pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3pto1p_slice[h]= new TH1F(Form("h1_Etot_3pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3pto2p_slice[h]= new TH1F(Form("h1_Erec_3pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3pto2p_slice[h]= new TH1F(Form("h1_Etot_3pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3p1pi_slice[h]= new TH1F(Form("h1_Erec_3p1pi_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3p1pi_slice[h]= new TH1F(Form("h1_Etot_3p1pi_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_43pto1p_slice[h]= new TH1F(Form("h1_Etot_43pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_43pto1p_slice[h]= new TH1F(Form("h1_Erec_43pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto3p_slice[h]= new TH1F(Form("h1_Erec_4pto3p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto3p_slice[h]= new TH1F(Form("h1_Etot_4pto3p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto2p_slice[h]= new TH1F(Form("h1_Erec_4pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto2p_slice[h]= new TH1F(Form("h1_Etot_4pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto1p_slice[h]= new TH1F(Form("h1_Erec_4pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto1p_slice[h]= new TH1F(Form("h1_Etot_4pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_Npi0[h] = new TH1F(Form("h1_Etot_Npi0_%d",h+1),"",n_bins,x_values);
   h1_Erec_Npi0[h] = new TH1F(Form("h1_Erec_Npi0_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_pipl_pimi_fact[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_pipl_pimi_fact_pipl[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pipl_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_pipl_pimi_fact_pimi[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pimi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_pipl_pimi_new_fact[h]= new TH1F(Form("h1_Erec_bkgd_pipl_pimi_new_fact_%d",h+1),"",n_bins,x_values);
   h1_Etot_Npi1[h] = new TH1F(Form("h1_Etot_Npi1_%d",h+1),"",n_bins,x_values);
   h1_Erec_Npi1[h] = new TH1F(Form("h1_Erec_Npi1_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_1p2pi[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_1p2pi[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to2p0pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to2p0pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p2pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p2pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p2pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p2pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_1p2pi_1p0pi[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_1p2pi_1p0pi[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_1p3pi[h] = new TH1F(Form("h1_Etot_bkgd_1p3pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_1p3pi[h] = new TH1F(Form("h1_Erec_bkgd_1p3pi_%d",h+1),"",n_bins,x_values);
  }


  for(int h=0;h<20;h++){
    h2_N_pi_phot[h]=new TH2F(Form("h2_N_pi_phot_%d",h),"",10,0,5,10,0,5);
  }

  TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();


/** Beginning of Event Loop **/
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<200000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    int nb = GetEntry(jentry);
    if( jentry%200000 == 0 )
    {
	         gDirectory->Write("hist_Files", TObject::kOverwrite);
		       cout<<jentry<<endl;
    }

      fTorusCurrent=2250; //only true for 2.2 and 4.4 GeV
//1500 is not used anymore, no runnumber in GENIE simulation files
  /*  if((runnb>18283 && runnb<18289) || (runnb>18300 && runnb<18304) || (runnb>18317 && runnb<18329))        fTorusCurrent=750;    //setting appropriate torrus magnet current
    else if ((runnb>18293 && runnb<18301) || (runnb>18305 && runnb<18317) || (runnb>18328 && runnb<18336))  fTorusCurrent=1500;
    else fTorusCurrent=2250;
*/

    if(jentry == 0){ //was n_evt == 1 before but jentry = n_evnt - 1

          fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
          fiducialcut->SetFiducialCutParameters(fbeam_en);
          std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
          rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
          std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;

    }

    //Resets q vector to (0,0,0)
    rotation->ResetQVector();
		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Outgoing e', Uncorr and corrected are the same.Needs to be changed later to one TLorentzVector

		TLorentzVector V4_el(pxl,pyl,pzl,El);
		TLorentzVector V4_el_uncorr(pxl,pyl,pzl,El);
		TVector3 V3_el(pxl,pyl,pzl);

		//Smearing of Electron Vector from Simulation
		double SmearedPe = gRandom->Gaus(pl,reso_e*pl);
		double SmearedEe = sqrt( SmearedPe*SmearedPe + e_mass * e_mass );
		V3_el.SetXYZ(SmearedPe/pl * pxl,SmearedPe/pl * pyl,SmearedPe/pl * pzl);
		V4_el.SetPxPyPzE(V3_el.X(),V3_el.Y(),V3_el.Z(),SmearedEe);
		double phi_ElectronOut = V3_el.Phi(); //in Radians
// apapadop
		V3_el.SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180),  GENIE coordinate system flipped with respect to CLAS

		//Fiducial Cuts with the smeared values
		if (!EFiducialCut(fbeam_en,V3_el) ) continue; // Electron theta & phi fiducial cuts

//apapadop
		phi_ElectronOut += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
		double el_momentum = V3_el.Mag();
		double el_theta = V3_el.Theta();

		//acceptance_c takes phi in radians and here unmodified by 30 degree.
		double e_acc_ratio = acceptance_c(el_momentum, cos(el_theta), phi_ElectronOut, 11,file_acceptance);


		double el_phi_mod = phi_ElectronOut*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
		if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree

		//Calculated Mott Cross Section and Weights for Inclusive Histograms
		double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(el_theta)+1))/(2*pow(El,2)*pow((1-cos(el_theta)),2.));
		double WeightIncl = wght*e_acc_ratio / Mott_cross_sec;

		//converting theta to degrees
		el_theta = el_theta*TMath::RadToDeg();

      //Energy reconstruction for electron only method
  //double E_rec_old = (2*(m_prot-bind_en[ftarget])*V4_el.E()+m_prot*m_prot-(m_prot-bind_en[ftarget])*(m_prot-bind_en[ftarget]))/(2*(m_prot-bind_en[ftarget]-V4_el.E()+V4_el.Rho()*cz[ind_em]));  //using the same value of single nucleon separation E for Ecal and Eqe
  //FROM AFROS CODE: THIS CORRESPONDS TO THE OLD E_REC
//  double ERecoOnlyePrime = (2*ProtonMass*BE + 2*ProtonMass*El - pow(ElectronMass,2)) / 2 / (ProtonMass - El + pElectronOut*cos(theta_ElectronOut));

    double E_rec = (m_prot*bind_en[ftarget]+m_prot*V4_el.E())/(m_prot-V4_el.E()+V4_el.Rho()*cos(el_theta));  //using the same value of single nucleon separation E Ecal and Eqe

    //Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
    double nu = -(V4_el-V4_beam).E();
    double reco_Q2 = -(V4_el-V4_beam).Mag2();
    double x_bjk = reco_Q2/(2*m_prot*nu);
    TVector3 V3_q = (V4_beam-V4_el).Vect();
    double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

    //Set q vector for the following rotations for the subtraction procedure
    rotation->SetQVector(V3_q);
//    rotation->PrintQVector();

    h1_el_mom->Fill(V4_el_uncorr.Rho(),WeightIncl);
    h1_el_mom_corr->Fill(V4_el.Rho(),WeightIncl);
    h1_el_mom_ratio->Fill(V4_el.Rho()/V4_el_uncorr.Rho(),WeightIncl);
    h2_el_pcorr_puncorr->Fill(V4_el.Rho(),V4_el.Rho()/V4_el_uncorr.Rho(),WeightIncl);
    h2_el_mom_diff->Fill(V4_el.Rho(),V4_el.Rho()-V4_el_uncorr.Rho(),WeightIncl);

    //Filling Histogram for electron kinematics
    h1_xbjk->Fill(x_bjk);
    h1_xbjk_weight->Fill(x_bjk,WeightIncl);
    h1_Q2->Fill(reco_Q2);
    h1_Q2_weight->Fill(reco_Q2,WeightIncl);
    h2_Q2_nu->Fill(nu,reco_Q2);
    h2_Q2_nu_weight->Fill(nu,reco_Q2,WeightIncl);
    h2_Q2_xbjk_weight->Fill(x_bjk,reco_Q2,WeightIncl);
    h1_Wvar->Fill(W_var);
    h1_Wvar_weight->Fill(W_var,WeightIncl);
    h2_Q2_W->Fill(W_var,reco_Q2);
    h2_xB_W->Fill(W_var,x_bjk);
    h2_Q2_W_weight->Fill(W_var,reco_Q2,WeightIncl);

    h1_el_Mott_crosssec->Fill(Mott_cross_sec);
    h2_el_theta_phi->Fill(el_phi_mod,el_theta,WeightIncl);
    h1_el_theta->Fill(el_theta);

    //Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

    //Index variables for hadrons (p and pions)
    int index_p[20]; //index for each proton
    int index_pi[20]; //index for each pion
    int ind_pi_phot[20];
    int index_pipl[20]; //index for each pi plus
    int index_pimi[20]; //index for each pi minus

    //Number of hadrons
    int num_p = 0;
    int num_pi = 0;
    int num_pi_phot = 0; //couting all pions and photons
    int num_pimi = 0;
    int num_pipl = 0;
    int num_pi_phot_nonrad=0; //counting all pions and non-radiation photons
    int num_phot_rad = 0; //counting radiation photons
    //Index and number variables for neutral particles
    int ec_index_n[20];
    int ec_num_n = 0;
    bool ec_radstat_n[20];

    //Array initialize to -1
    for (int i = 0; i<20; i++) {
      index_p[i] = -1;   index_pi[i] = -1;   index_pipl[i] = -1;   index_pimi[i] = -1;   ind_pi_phot[i] = -1;
      ec_index_n[i] = -1;   ec_radstat_n[i] = false;
    }

    const double phot_rad_cut = 40;
    const double phot_e_phidiffcut=30; //electron - photon phi difference cut

    // Creating vectors to store id of particles in the array
    vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
    ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

   //Loop for Hadrons
    for (int i = 0; i < nf; i++)
    {
        //Start of proton selection
       if (pdgf[i] == 2212 && pf[i] > 0.3) {

          num_p = num_p + 1;
          index_p[num_p - 1] = i;
          ProtonID.push_back(i);

       }

       if (pdgf[i] == -211 && pf[i] > 0.15)  {

          num_pimi = num_pimi + 1;
          num_pi = num_pi + 1;
          num_pi_phot = num_pi_phot + 1;
          num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
          index_pimi[num_pi_phot - 1] = i;
          index_pi[num_pi_phot - 1] = i;
          ind_pi_phot[num_pi_phot - 1] = i;
          PiMinusID.push_back(i);

       }

       if ( pdgf[i] == 211 && pf[i] > 0.15)  {

          num_pipl = num_pipl + 1;
          num_pi  = num_pi + 1;
          num_pi_phot = num_pi_phot + 1;
          num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
          index_pipl[num_pi_phot - 1] = i;
          index_pi[num_pi_phot - 1] = i;
          ind_pi_phot[num_pi_phot - 1] = i;
          PiPlusID.push_back(i);

       }

       if (pdgf[i] == 22 && pf[i] > 0.3) {

          ec_num_n = ec_num_n + 1;
          num_pi_phot = num_pi_phot + 1;
          ind_pi_phot[num_pi_phot - 1] = i;  //no events with photons for now F.H. 28.08.19
          PhotonID.push_back(i);
          // WARNING: THe following needs to be implemented for simulation data F.H. 24.08.19
          //Cut on Radiation photon via angle with respect to the electron
          V3_phot_angles.SetXYZ(pxf[i],pyf[i],pzf[i]);
          double neut_phi_mod = V3_phot_angles.Phi()*TMath::RadToDeg() + 30; //Add 30 degree
          if (neut_phi_mod < 0) neut_phi_mod = neut_phi_mod + 360;  //Neutral particle is between 0 and 360 degree

          //within 40 degrees in theta and 30 degrees in phi
          if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg() < phot_rad_cut && fabs(neut_phi_mod-el_phi_mod) < phot_e_phidiffcut ) {
             ec_radstat_n[num_pi_phot - 1] = true; //select radiation photons
             num_phot_rad = num_phot_rad + 1;
          }

       }

    } //end of hadron loop


    //Skip event if there is at least one radiation photon
    if (num_phot_rad > 0) {
      continue;
    }

    //Filling Histograms with multiplicities
    h1_Npi->Fill(num_pi);
    h1_Nprot->Fill(num_p);
    h1_Nphot->Fill(ec_num_n);
    h1_Npipl->Fill(num_pipl);
    h1_Npimi->Fill(num_pimi);
    h1_Npiphot->Fill(num_pi_phot);
    h1_Npiphot_norad->Fill(num_pi_phot_nonrad);
    h2_N_prot_pi->Fill(num_pi,num_p);
    h2_N_prot_pi_phot->Fill(num_pi+ec_num_n,num_p);
    h2_N_prot_pi_phot_nonrad->Fill(num_pi_phot_nonrad,num_p);
    h2_N_pi_phot[num_p]->Fill(ec_num_n,num_pi);

	//Events with exactly 2 protons
	if(num_p == 2) {

          double SmearedPp;
          double SmearedEp;
        //LorentzVectors for protons without momentum smearing
          TLorentzVector V4_prot_uncorr1(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
          TLorentzVector V4_prot_uncorr2(pxf[index_p[1]],pyf[index_p[1]],pzf[index_p[1]],TMath::Sqrt(m_prot*m_prot+pf[index_p[1]]*pf[index_p[1]]));
        //Kinematic for first proton, smearing
    	   	SmearedPp = gRandom->Gaus(pf[index_p[0]],reso_p*pf[index_p[0]]);
      		SmearedEp = sqrt( SmearedPp*SmearedPp + m_prot * m_prot );

	        TVector3 V3_prot_corr1;
	        V3_prot_corr1.SetXYZ(SmearedPp/pf[index_p[0]] * pxf[index_p[0]],SmearedPp/pf[index_p[0]] * pyf[index_p[0]],SmearedPp/pf[index_p[0]] * pzf[index_p[0]]);
          // apapadop
          double phi_prot1 = V3_prot_corr1.Phi();
          V3_prot_corr1.SetPhi(phi_prot1 + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS

          if (!PFiducialCut(fbeam_en, V3_prot_corr1) ) { continue; } // Proton theta & phi fiducial cuts
          phi_prot1 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

	        double p_theta1 = V3_prot_corr1.Theta();
	        double prot_mom_corr1 = V3_prot_corr1.Mag();
	        //Proton 1 weight
	        double p_acc_ratio1 = acceptance_c(prot_mom_corr1, cos(p_theta1), phi_prot1, 2212,file_acceptance_p);
	        p_theta1 = p_theta1*TMath::RadToDeg();

	        //Kinematic for second proton
	        SmearedPp = gRandom->Gaus(pf[index_p[1]],reso_p*pf[index_p[1]]);
	        SmearedEp = sqrt( SmearedPp*SmearedPp + m_prot * m_prot );

	        TVector3 V3_prot_corr2;
	        V3_prot_corr2.SetXYZ(SmearedPp/pf[index_p[1]] * pxf[index_p[1]],SmearedPp/pf[index_p[1]] * pyf[index_p[1]],SmearedPp/pf[index_p[1]] * pzf[index_p[1]]);

          // apapadop
          double phi_prot2 = V3_prot_corr2.Phi();
          V3_prot_corr2.SetPhi(phi_prot2 + TMath::Pi()); // Vec.Phi() is between (-180,180) // GENIE coordinate system flipped with respect to CLAS

          if (!PFiducialCut(fbeam_en, V3_prot_corr2) ) { continue; } // Proton theta & phi fiducial cuts
          phi_prot2 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

          double p_theta2 = V3_prot_corr2.Theta();
          double prot_mom_corr2 = V3_prot_corr2.Mag();
          //Proton 2 weight
          double p_acc_ratio2 = acceptance_c(prot_mom_corr2, cos(p_theta2), phi_prot2, 2212,file_acceptance_p);
          p_theta2 = p_theta2*TMath::RadToDeg();

          //Total proton weight
          double weight_protons = p_acc_ratio1 * p_acc_ratio2;

          TVector3 V3_2prot_uncorr[2];
          V3_2prot_uncorr[0] = V4_prot_uncorr1.Vect();
          V3_2prot_uncorr[1] = V4_prot_uncorr2.Vect();

          TVector3 V3_2prot_corr[2];
          V3_2prot_corr[0] = V3_prot_corr1;
          V3_2prot_corr[1] = V3_prot_corr2;

	//---------------------------------- 2p 0pi->  1p0pi   ----------------------------------------------


          double E_tot_2p[2]={0};
          double p_perp_tot_2p[2]={0};
          N_prot_both = 0;
          double P_N_2p[2]={0};

          rotation->prot2_rot_func( V3_2prot_corr, V3_2prot_uncorr, V4_el, E_tot_2p, p_perp_tot_2p, P_N_2p , &N_prot_both);

          if(num_pi_phot==0 && N_prot_both!=0){
      //    if(num_pi == 0 && N_prot_both != 0){
            double histoweight = weight_protons*e_acc_ratio*wght/Mott_cross_sec; //total weight from 2p acceptance , 1e acceptance, Mott, and GENIE weight

            for(int f = 0; f < num_p; f++){    //looping through two protons
              h1_E_tot_p_bkgd->Fill(E_tot_2p[f],P_N_2p[f]*histoweight);
              h1_E_rec_p_bkgd->Fill(E_rec,P_N_2p[f]*histoweight);
              h2_Erec_pperp_2p->Fill(p_perp_tot_2p[f],E_rec,P_N_2p[f]*histoweight);
              h1_E_tot_p_bkgd_fracfeed->Fill((E_tot_2p[f]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_N_2p[f]*histoweight);
              h1_E_rec_p_bkgd_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_N_2p[f]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[f],-P_N_2p[f]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[f]) *TMath::RadToDeg(),-P_N_2p[f]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[f],-P_N_2p[f]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[f],-P_N_2p[f]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[f],-P_N_2p[f]*histoweight);

              for(int i = 0 ; i < n_slice; i++)
              {
	               if (p_perp_tot_2p[f] < pperp_max[i] && p_perp_tot_2p[f] > pperp_min[i]) {
	                  h1_Etot_p_bkgd_slice[i]->Fill(E_tot_2p[f],P_N_2p[f]*histoweight);
	                  h1_Erec_p_bkgd_slice[i]->Fill(E_rec,P_N_2p[f]*histoweight);
	               }
	            }
            }//looping through two protons

            h1_E_tot_2p_det->Fill(E_tot_2p[0],histoweight);
            h1_E_rec_2p_det->Fill(E_rec,histoweight);
          }//no pions cut and N_prot_both!=0

	//---------------------------------- 2p 1pi   ----------------------------------------------
 //Const int can be placed somewhere up after if for 2 protons F.H. 05.09.19
          const int N_2prot=2;
//Variable might can be placed in a more local context F.H. 05.09.19
          double Ecal_2p1pi_to2p0pi[N_2prot]={0};
          double p_miss_perp_2p1pi_to2p0pi[N_2prot]={0};

          if (num_pi_phot==1) {
//          if (num_pi == 1) {

            int charge = 0;
            if (num_pimi == 1 && num_pipl == 0) {
              charge = -1;
            }
            else if (num_pipl == 1 && num_pimi == 0) {
              charge = 1;
            }
            //radiation status is false if real photon
            else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
              charge = 0;
            }
            //skip radiation photon
            else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
	            charge = 0;
            }
            else {
              std::cout << "2 proton 1 pi loop Problem with Pion/Photon count and identification. Num Pi = " << num_pi << " , Num PiPlus = " << num_pipl << " , Num PiMinus = " << num_pimi << std::endl;
              continue;
            }


            SmearedPp = gRandom->Gaus(pf[index_pi[0]],reso_pi*pf[index_pi[0]]);
            SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

            TVector3 V3_1pi_corr;
            V3_1pi_corr.SetXYZ(SmearedPp/pf[index_pi[0]] * pxf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pyf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pzf[index_pi[0]]);

            // apapadop
            double phi_pion = V3_1pi_corr.Phi();
            V3_1pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
            phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

            // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
            if ( !Pi_phot_fid_united(fbeam_en, V3_1pi_corr, charge) )
            {
              continue;
            }

            double pion_theta = V3_1pi_corr.Theta();
            double pion_mom_corr = V3_1pi_corr.Mag();

            double pion_acc_ratio = 0;
            if (charge == 1) { //acceptance for pi plus
              pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
            }
            else if (charge == -1) {    //acceptance for pi minus. using electron acceptance map
              pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
            }
            else if (charge == 0) {    //acceptance for neutral, setting to 1 for now F.H. 09/24/19
              pion_acc_ratio = 1;
            }
            else { std::cout << "WARNING: 2proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

            double P_2p1pito2p0pi[2] = {0};
            double P_2p1pito1p1pi[2] = {0};
            double P_2p1pito1p0pi[2] = {0};
            double Ptot = 0;

            rotation->prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_1pi_corr, charge, ec_radstat_n[0],V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

            double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19

            for(int z=0; z < N_2prot; z++){ //looping over two protons

	//---------------------------------- 2p 1pi ->2p 0pi   ----------------------------------------------

              h1_E_tot_2p1pi_2p0pi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h1_E_rec_2p1pi_2p0pi->Fill(E_rec,P_2p1pito2p0pi[z]*histoweight);
              h2_Erec_pperp_2p1pi_2p0pi->Fill(p_miss_perp_2p1pi_to2p0pi[z],E_rec,P_2p1pito2p0pi[z]*histoweight);
              h1_E_tot_2p1pi_2p0pi_fracfeed->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
              h1_E_rec_2p1pi_2p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z]) *TMath::RadToDeg(),P_2p1pito2p0pi[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);

              for(int i = 0; i < n_slice; i++){
                if (p_miss_perp_2p1pi_to2p0pi[z]<pperp_max[i] && p_miss_perp_2p1pi_to2p0pi[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
	                 h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(E_rec,P_2p1pito2p0pi[z]*histoweight);
                }
              }
	//---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

              h1_E_tot_2p1pi_1p1pi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h1_E_rec_2p1pi_1p1pi->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
              h2_Erec_pperp_2p1pi_1p1pi->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p1pi[z]*histoweight);
              h1_E_tot_2p1pi_1p1pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
              h1_E_rec_2p1pi_1p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_2p1pito1p1pi[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

              for(int i = 0; i < n_slice; i++){
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
                  h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
                  h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
                }
              }
//---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------

              h1_E_tot_2p1pi_1p0pi->Fill(E_tot_2p[z], P_2p1pito1p0pi[z]*histoweight);
              h1_E_rec_2p1pi_1p0pi->Fill(E_rec,P_2p1pito1p0pi[z]*histoweight);
              h2_Erec_pperp_2p1pi_1p0pi->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p0pi[z]*histoweight);
              h1_E_tot_2p1pi_1p0pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p0pi[z]*histoweight);
              h1_E_rec_2p1pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p0pi[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),-P_2p1pito1p0pi[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

              for(int i = 0; i < n_slice; i++) {
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*histoweight);
	                 h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_rec,P_2p1pito1p0pi[z]*histoweight);
                }
              }

            }//filling the histograms for 2protons
          }//1pi requirement


	//---------------------------------- 2p 2pi   ----------------------------------------------

          const int N_2pi=2;
          int q_pi2[2];
          double Ecal_2p2pi[N_2prot];
          double p_miss_perp_2p2pi[N_2prot];
          double Ptot_2p[2]={0};
          bool ecstat_pi2[N_2pi]={false};

          if (num_pi_phot==2) { //no photons for now F.H. 29.8.19
//          if (num_pi == 2) {

            TVector3 V3_2pi_corr[N_2pi];
            double pion_acc_ratio[N_2pi] = {0};
            for (int i = 0; i < num_pi_phot; i++) {

              ecstat_pi2[i] = ec_radstat_n[i];
              if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
                q_pi2[i] = 1;
              }
              else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
                q_pi2[i] = -1;
              }
              else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
                q_pi2[i] = 0;
              }
              else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
                q_pi2[i] = 0;
              }
              else {  std::cout << "WARNING: 2p 2pion event: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

              SmearedPp = gRandom->Gaus(pf[index_pi[i]],reso_pi*pf[index_pi[i]]);
              SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

              V3_2pi_corr[i].SetXYZ(SmearedPp/pf[index_pi[i]] * pxf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pyf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pzf[index_pi[i]]);

              // apapadop
              double phi_pion = V3_2pi_corr[i].Phi();
              V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
              phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

              double pion_theta = V3_2pi_corr[i].Theta();
              double pion_mom_corr = V3_2pi_corr[i].Mag();

              if (q_pi2[i] == 1) { //acceptance for pi plus
                pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
              }
              else if (q_pi2[i] == -1) {    //acceptance for pi minus. using electron acceptance map
                pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
              }
              else if (q_pi2[i] == 0) {    //acceptance for photon set to 1 for now F.H. 09/24/19
                pion_acc_ratio[i] = 1;
              }
              else { std::cout << "WARNING: 2proton and 2 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

            }

            // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
            // Pi+Pi_phot_fid_united with 0 is for Photons/Pi0
            //If any fiducial is not fullfilled continue to next event
            if ( ( !Pi_phot_fid_united(fbeam_en, V3_2pi_corr[0],  q_pi2[0]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_2pi_corr[1],  q_pi2[1]) )  )
            {
              continue;
            }

            rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,q_pi2,ecstat_pi2 ,V4_el, Ecal_2p2pi,p_miss_perp_2p2pi,Ptot_2p);

            double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
            double histoweight = weight_pions * weight_protons * e_acc_ratio * wght/Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19


            for(int z = 0; z < N_2prot; z++){ //looping over two protons

		//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

              h1_E_tot_2p2pi->Fill(E_tot_2p[z], Ptot_2p[z]*histoweight);
              h1_E_rec_2p2pi->Fill(E_rec,Ptot_2p[z]*histoweight);
              h2_Erec_pperp_2p2pi->Fill(p_perp_tot_2p[z],E_rec,Ptot_2p[z]*histoweight);
              h1_E_tot_2p2pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],Ptot_2p[z]*histoweight);
              h1_E_rec_2p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],Ptot_2p[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],Ptot_2p[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),Ptot_2p[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],Ptot_2p[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],Ptot_2p[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],Ptot_2p[z]*histoweight);

              for(int i = 0; i < n_slice; i++)
              {
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p2pi[i]->Fill(E_tot_2p[z],Ptot_2p[z]*histoweight);
	                 h1_Erec_p_bkgd_slice_2p2pi[i]->Fill(E_rec,Ptot_2p[z]*histoweight);
                }
              }
            }   //Filling the histogram for two protons
          }//2pi requirement

    } //2prot requirement

  //Events with exactly 3 protons
	  if(num_p == 3){

        double SmearedPp;
        double SmearedEp;

	      const int N_3p = 3;
	      TLorentzVector V4_p_uncorr[N_3p], V4_p_corr[N_3p],V4_prot_el[N_3p];
	      TVector3 V3_prot_uncorr[N_3p],V3_prot_corr[N_3p],V3_3p_rot[N_3p];
	      double E_cal[N_3p],p_miss_perp[N_3p],P_3pto1p[N_3p];
	      double N_p1[N_3p]={0};
        double N_p_three=0;
        double E_cal_3pto1p[3]={0};
        double p_miss_perp_3pto1p[3]={0};
	      int N_comb = 3;
	      const int N_2p = 2;
	      double E_cal_3pto2p[3][N_2p]={0};
        double p_miss_perp_3pto2p[3][N_2p]={0};
        double P_3pto2p[3][N_2p]={0};
	      TVector3 V3_2p_rot[N_2p], V3_prot_el[N_3p][N_2p];
        double p_acc_ratio[N_3p] = {0};

	      for(int i = 0; i < N_3p; i++)
	      {
          N_p1[i] = 0;

          V4_p_uncorr[i].SetPxPyPzE(pxf[index_p[i]],pyf[index_p[i]],pzf[index_p[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
          V3_prot_uncorr[i] = V4_p_uncorr[i].Vect();

         //Smearing of Proton
          SmearedPp = gRandom->Gaus(pf[index_p[i]],reso_p*pf[index_p[i]]);
          SmearedEp = TMath::Sqrt( SmearedPp*SmearedPp + m_prot * m_prot );

          V3_prot_corr[i].SetXYZ(SmearedPp/pf[index_p[i]] * pxf[index_p[i]],SmearedPp/pf[index_p[i]] * pyf[index_p[i]],SmearedPp/pf[index_p[i]] * pzf[index_p[i]]);
          V4_p_corr[i].SetPxPyPzE(SmearedPp/pf[index_p[i]] * pxf[index_p[i]],SmearedPp/pf[index_p[i]] * pyf[index_p[i]],SmearedPp/pf[index_p[i]] * pzf[index_p[i]],SmearedEp);

          // apapadop
          double phi_prot = V3_el.Phi(); //in Radians
          V3_prot_corr[i].SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180)
          phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

          double p_theta = V3_prot_corr[i].Theta();
          double prot_mom_corr = V3_prot_corr[i].Mag();
          //Proton acceptance weight
          p_acc_ratio[i] = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);

	        V4_prot_el[i] = V4_p_corr[i] + V4_el;
	        E_cal[i] = V4_el.E()+ V4_p_corr[i].E() - m_prot + bind_en[ftarget];
	        p_miss_perp[i] = TMath::Sqrt(V4_prot_el[i].Px()*V4_prot_el[i].Px() + V4_prot_el[i].Py()*V4_prot_el[i].Py());
	      }

        //If one of the smeared protons doesn't get through fiducials, skip event
        if ( (!PFiducialCut(fbeam_en, V3_prot_corr[0]) ) || (!PFiducialCut(fbeam_en, V3_prot_corr[1]) ) || (!PFiducialCut(fbeam_en, V3_prot_corr[2]) ) ) { continue; }


        V3_prot_el[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
        V3_prot_el[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
        V3_prot_el[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
        V3_prot_el[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
        V3_prot_el[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
        V3_prot_el[2][1]=V4_el.Vect()+V3_prot_uncorr[2];


        rotation->prot3_rot_func( V3_prot_corr,V3_prot_uncorr,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

        //acceptance weight for all three protons
        double weight_protons =  p_acc_ratio[0] * p_acc_ratio[1] * p_acc_ratio[2];

	      if(num_pi_phot==0 && N_p_three!=0){
//	      if(num_pi == 0 && N_p_three != 0){    //no photons for now F.H. 29.8.19
           double histoweight = weight_protons * e_acc_ratio * wght/Mott_cross_sec; //Weight for 3protons, 1 electron, GENIE weight and Mott cross section
	         for(int count = 0; count < N_comb; count++)    { //Loop over number of combinations
               for(int j = 0; j < N_2p; j++)    { //loop over two protons

		//-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------

                  h1_E_tot_3pto2p->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*histoweight);
		              h1_E_rec_3pto2p->Fill(E_rec, P_3pto2p[count][j]*histoweight);
		              h2_Erec_pperp_321p->Fill(p_miss_perp_3pto2p[count][j],E_rec,P_3pto2p[count][j]*histoweight);
		              h1_E_tot_3pto2p_fracfeed->Fill((E_cal_3pto2p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto2p[count][j]*histoweight);
		              h1_E_rec_3pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_3pto2p[count][j]*histoweight);
		              h2_pperp_W->Fill(W_var,p_miss_perp_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
		              h1_theta0->Fill((V4_beam.Vect()).Angle(V3_prot_el[count][j])*TMath::RadToDeg(),P_3pto2p[count][j]*histoweight);
	                h2_Ecal_Eqe->Fill(E_rec,E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
	                h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
	                h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);

		              for(int i = 0; i < n_slice; i++)
		              {
		                  if (p_miss_perp_3pto2p[count][j]<pperp_max[i] && p_miss_perp_3pto2p[count][j]>pperp_min[i]){
	                         h1_Etot_3pto2p_slice[i]->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*histoweight);
			                     h1_Erec_3pto2p_slice[i]->Fill(E_rec, P_3pto2p[count][j]*histoweight);
		                  }
		              }
               } //end loop over protons
	         } //end loop over combination N_comb

		//-----------------------------------------  3p to 1p  -----------------------------------------------------------------------

	         for(int j = 0; j < N_3p; j++)    {

		             P_3pto1p[j]= N_p1[j]/N_p_three;
		             h1_E_tot_3pto1p->Fill(E_cal[j], P_3pto1p[j]*histoweight);
		             h1_E_rec_3pto1p->Fill(E_rec,P_3pto1p[j]*histoweight);
		             h2_Erec_pperp_31p->Fill(p_miss_perp[j],E_rec,P_3pto1p[j]*histoweight);
		             h1_E_tot_3pto1p_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto1p[j]*histoweight);
		             h1_E_rec_3pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_3pto1p[j]*histoweight);
		             h2_pperp_W->Fill(W_var,p_miss_perp[j],-P_3pto1p[j]*histoweight);
		             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),-P_3pto1p[j]*histoweight);
		             h2_Ecal_Eqe->Fill(E_rec,E_cal[j],-P_3pto1p[j]*histoweight);
		             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],-P_3pto1p[j]*histoweight);
		             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],-P_3pto1p[j]*histoweight);

		             for(int i = 0; i < n_slice; i++)
		             {
		                 if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
		                     h1_Etot_3pto1p_slice[i]->Fill(E_cal[j],P_3pto1p[j]*histoweight);
		                     h1_Erec_3pto1p_slice[i]->Fill(E_rec,P_3pto1p[j]*histoweight);
		                 }
		             }
	         } //end loop over N_3p
	      }//end if num_pi_phot==0 && N_p_three!=0, no pions

	//----------------------------------3p 1pi ----------------------------------------------------------

        if (num_pi_phot==1) { //no photons for now F.H. 29.8.19
//        if (num_pi == 1) { //number of pions  = 1

             double P_tot_3p[N_3p]={0};
             double Ecal_3p1pi[N_3p]={0};
             double p_miss_perp_3p1pi[N_3p]={0};
             double charge = 0;
             TVector3 V3_pi_corr;

             if (num_pimi == 1 && num_pipl == 0) {
               charge = -1;
             }
             else if (num_pipl == 1 && num_pimi == 0) {
               charge = 1;
             }
             //radiation status is false if real photon
             else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
               charge = 0;
             }
             //skip radiation photon
             else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
 	            charge = 0;
             }
             else {  std::cout << "WARNING: 3p 1pion event: No charge for one pion/photon could be assigned. Pion number " << index_pi[0] << std::endl; continue; }

             SmearedPp = gRandom->Gaus(pf[index_pi[0]],reso_pi*pf[index_pi[0]]);
             SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

             V3_pi_corr.SetXYZ(SmearedPp/pf[index_pi[0]] * pxf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pyf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pzf[index_pi[0]]);

             // apapadop
             double phi_pion = V3_pi_corr.Phi(); //in Radians
             V3_pi_corr.SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180)
             phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

             // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
             if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, charge) )
             {
               continue;
             }

             double pion_theta = V3_pi_corr.Theta();
             double pion_mom_corr = V3_pi_corr.Mag();

             double pion_acc_ratio = 0;
             if (charge == 1) { //acceptance for pi plus
               pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
             }
             else if (charge == -1) {    //acceptance for pi minus. using electron acceptance map
               pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
             }
             else if (charge == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
               pion_acc_ratio = 1;
             }
             else { std::cout << "WARNING: 3proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }


             rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot_uncorr, V3_pi_corr, charge ,ec_radstat_n[0], V4_el,  Ecal_3p1pi,p_miss_perp_3p1pi, P_tot_3p);

             double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec; //Weight for 3protons, 1 pion, 1 electron, GENIE weight and Mott cross section

	           for(int j = 0; j < N_3p; j++)    { //loop over 3 protons

                 h1_E_tot_3p1pi->Fill(E_cal[j], P_tot_3p[j]*histoweight);
		             h1_E_rec_3p1pi->Fill(E_rec,P_tot_3p[j]*histoweight);
		             h2_Erec_pperp_3p1pi->Fill(p_miss_perp[j],E_rec,P_tot_3p[j]*histoweight);
		             h1_E_tot_3p1pi_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_tot_3p[j]*histoweight);
		             h1_E_rec_3p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_3p[j]*histoweight);
		             h2_pperp_W->Fill(W_var,p_miss_perp[j],P_tot_3p[j]*histoweight);
		             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_tot_3p[j]*histoweight);
		             h2_Ecal_Eqe->Fill(E_rec,E_cal[j],P_tot_3p[j]*histoweight);
		             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],P_tot_3p[j]*histoweight);
		             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],P_tot_3p[j]*histoweight);

		             for(int i = 0; i < n_slice; i++)
		             {
		                 if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
                         h1_Etot_3p1pi_slice[i]->Fill(E_cal[j],P_tot_3p[j]*histoweight);
		                     h1_Erec_3p1pi_slice[i]->Fill(E_rec,P_tot_3p[j]*histoweight);
		                 }
		             }
             } //end loop over N_3p
        }// 1 pi requirement ends

	  }//end if num_p == 3  3proton requiremnet

//Events with exactly 4 protons
	 	 if(num_p==4){

        double SmearedPp;
        double SmearedEp;

	      const int N_p4=4;
	      TLorentzVector V4_p4_uncorr[N_p4], V4_p4_corr[N_p4],V4_prot4_el[N_p4];
	      TVector3 V3_prot4_uncorr[N_p4],V3_prot4_corr[N_p4];
	      double E_cal_p4[N_p4]={0};
        double p_miss_perp_p4[N_p4]={0};
        double P_4pto1p[N_p4]={0};
	      TVector3 V3_p4_rot[N_p4];
	      bool prot4_stat[N_p4]={false};
	      const int Ncomb_4to1 = 4,Ncomb_4to2 = 6, Ncomb_4to3 = 4;
	      double N_p4_p1[Ncomb_4to1]={0};
        double N_p4_p2[Ncomb_4to2]={0};
        double N_p4_p3[Ncomb_4to3]={0};
        double N_p_four = 0;
        double p_acc_ratio[N_p4];

	      for(int i = 0; i < N_p4; i++)  //loop over 4 protons
	      {

          V4_p4_uncorr[i].SetPxPyPzE(pxf[index_p[i]],pyf[index_p[i]],pzf[index_p[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
          V3_prot4_uncorr[i] = V4_p4_uncorr[i].Vect();

         //Smearing of Proton
          SmearedPp = gRandom->Gaus(pf[index_p[i]],reso_p*pf[index_p[i]]);
          SmearedEp = TMath::Sqrt( SmearedPp*SmearedPp + m_prot * m_prot );

          V3_prot4_corr[i].SetXYZ(SmearedPp/pf[index_p[i]] * pxf[index_p[i]],SmearedPp/pf[index_p[i]] * pyf[index_p[i]],SmearedPp/pf[index_p[i]] * pzf[index_p[i]]);
          V4_p4_corr[i].SetPxPyPzE(SmearedPp/pf[index_p[i]] * pxf[index_p[i]],SmearedPp/pf[index_p[i]] * pyf[index_p[i]],SmearedPp/pf[index_p[i]] * pzf[index_p[i]],SmearedEp);

          // apapadop
          double phi_prot = V3_prot4_corr[i].Phi(); //in Radians
          V3_prot4_corr[i].SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
          phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

          double p_theta = V3_prot4_corr[i].Theta();
          double prot_mom_corr = V3_prot4_corr[i].Mag();
          //Proton acceptance weight
	        p_acc_ratio[i] = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);

          V4_prot4_el[i] = V4_p4_corr[i] + V4_el;
          E_cal_p4[i] = V4_el.E() + V4_p4_corr[i].E() - m_prot + bind_en[ftarget];
          p_miss_perp_p4[i] = TMath::Sqrt(V4_prot4_el[i].Px()*V4_prot4_el[i].Px() + V4_prot4_el[i].Py()*V4_prot4_el[i].Py());
        }

        //If one of the smeared protons doesn't get through fiducials, skip event
        if ( (!PFiducialCut(fbeam_en, V3_prot4_corr[0]) ) || (!PFiducialCut(fbeam_en, V3_prot4_corr[1]) ) ||
             (!PFiducialCut(fbeam_en, V3_prot4_corr[2]) ) || (!PFiducialCut(fbeam_en, V3_prot4_corr[3]) ) )
        {
          continue;
        }

        //acceptance weight for all three protons
        double weight_protons =  p_acc_ratio[0] * p_acc_ratio[1] * p_acc_ratio[2] * p_acc_ratio[3];


        if ( num_pi_phot == 0){ //no pion or photon
	            for(int g = 0; g < N_tot; g++){ //this looks like a 4-proton rotation function -> could be placed maybe in an extra function

		              double rot_angle = gRandom->Uniform(0,2*TMath::Pi());
		              for(int i = 0; i < N_p4;i++) {
		                  V3_p4_rot[i]= V3_prot4_uncorr[i];
		                  V3_p4_rot[i].Rotate(rot_angle,V3_q);
		              }
		              for(int i_p = 0; i_p < N_p4; i_p++) {
                    prot4_stat[i_p] = PFiducialCut(fbeam_en, V3_p4_rot[i_p]);
                  }

		              if( prot4_stat[0]  && !prot4_stat[1]   && !prot4_stat[2] && !prot4_stat[3])  N_p4_p1[0]=N_p4_p1[0]+1;//Detecting 1p out of 4p
		              if(!prot4_stat[0]  &&   prot4_stat[1]  && !prot4_stat[2] && !prot4_stat[3])  N_p4_p1[1]=N_p4_p1[1]+1;
		              if(!prot4_stat[0]  &&  !prot4_stat[1]  &&  prot4_stat[2] && !prot4_stat[3])  N_p4_p1[2]=N_p4_p1[2]+1;
		              if(!prot4_stat[0]  &&  !prot4_stat[1]  && !prot4_stat[2] &&  prot4_stat[3])  N_p4_p1[3]=N_p4_p1[3]+1;
		              if( prot4_stat[0]  &&  prot4_stat[1]   &&  prot4_stat[2] &&  prot4_stat[3])  N_p_four=N_p_four+1;   //Detecting 4p out of 4p

		              if( prot4_stat[0]  &&  prot4_stat[1]  && !prot4_stat[2]  && !prot4_stat[3])  N_p4_p2[0]=N_p4_p2[0]+1;//Detecting 2p out of 4p
		              if( prot4_stat[0]  && !prot4_stat[1]  &&  prot4_stat[2]  && !prot4_stat[3])  N_p4_p2[1]=N_p4_p2[1]+1;
		              if( prot4_stat[0]  && !prot4_stat[1]  && !prot4_stat[2]  &&  prot4_stat[3])  N_p4_p2[2]=N_p4_p2[2]+1;
		              if(!prot4_stat[0]  &&  prot4_stat[1]  &&  prot4_stat[2]  && !prot4_stat[3])  N_p4_p2[3]=N_p4_p2[3]+1;
		              if(!prot4_stat[0]  &&  prot4_stat[1]  && !prot4_stat[2]  &&  prot4_stat[3])  N_p4_p2[4]=N_p4_p2[4]+1;
		              if(!prot4_stat[0]  && !prot4_stat[1]  &&  prot4_stat[2]  &&  prot4_stat[3])  N_p4_p2[5]=N_p4_p2[5]+1;


		              if( prot4_stat[0]  &&  prot4_stat[1]  &&  prot4_stat[2]  && !prot4_stat[3])  N_p4_p3[0]=N_p4_p3[0]+1;//Detecting 3p out of 4p
		              if( prot4_stat[0]  &&  prot4_stat[1]  &&  !prot4_stat[2] &&  prot4_stat[3])  N_p4_p3[1]=N_p4_p3[1]+1;
		              if( prot4_stat[0]  && !prot4_stat[1]  &&  prot4_stat[2]  &&  prot4_stat[3])  N_p4_p3[2]=N_p4_p3[2]+1;
		              if(!prot4_stat[0]  &&  prot4_stat[1]  &&  prot4_stat[2]  &&  prot4_stat[3])  N_p4_p3[3]=N_p4_p3[3]+1;
	            }//for loop of 4p rotations ends

                // still no pions
              int N_comb=3;    //number of 2 proton combination out of three
	            const int N_2p=2, N_3p=3;
	            double E_cal_4pto3p[3][N_2p]={0};
              double p_miss_perp_4pto3p[3][N_2p]={0};
              double P_4pto3p[3][N_2p]={0};
	            TVector3 V3_prot_corr[N_3p],V3_prot_uncorr[N_3p],V3_prot_el_4pto3p[N_3p][N_2p],V3_el_prot[N_comb][N_2p];
	            double N_p_three=0,N_p1[N_3p]={0};
              double E_cal_43pto1p[N_3p];
              double p_miss_perp_43pto1p[N_3p];
              double P_43pto1p[3]={0};

              double histoweight = weight_protons * e_acc_ratio * wght/Mott_cross_sec; //Weight for 3protons, 1 electron, GENIE weight and Mott cross section


	            for(int g = 0; g < Ncomb_4to3; g++){   //estimating the undetected 4p contribution to  3p
		             if(g==0) {
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[0]; V3_prot_uncorr[1]=V3_prot4_uncorr[1]; V3_prot_uncorr[2]=V3_prot4_uncorr[2];
		                 V3_prot_corr[0]=V3_prot4_corr[0]; V3_prot_corr[1]=V3_prot4_corr[1]; V3_prot_corr[2]=V3_prot4_corr[2];
		             }
		             if(g==1){
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[0]; V3_prot_uncorr[1]=V3_prot4_uncorr[1]; V3_prot_uncorr[2]=V3_prot4_uncorr[3];
		                 V3_prot_corr[0]=V3_prot4_corr[0]; V3_prot_corr[1]=V3_prot4_corr[1]; V3_prot_corr[2]=V3_prot4_corr[3];
		             }
		             if(g==2){
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[0]; V3_prot_uncorr[1]=V3_prot4_uncorr[2]; V3_prot_uncorr[2]=V3_prot4_uncorr[3];
		                 V3_prot_corr[0]=V3_prot4_corr[0]; V3_prot_corr[1]=V3_prot4_corr[2]; V3_prot_corr[2]=V3_prot4_corr[3];
		             }
		             if(g==3){
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[1]; V3_prot_uncorr[1]=V3_prot4_uncorr[2]; V3_prot_uncorr[2]=V3_prot4_uncorr[3];
		                 V3_prot_corr[0]=V3_prot4_corr[1]; V3_prot_corr[1]=V3_prot4_corr[2]; V3_prot_corr[2]=V3_prot4_corr[3];
		             }

		             rotation->prot3_rot_func(V3_prot_corr, V3_prot_uncorr,V4_el,E_cal_4pto3p,p_miss_perp_4pto3p, P_4pto3p,N_p1,E_cal_43pto1p,p_miss_perp_43pto1p,&N_p_three);

		             V3_el_prot[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
		             V3_el_prot[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
		             V3_el_prot[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
		             V3_el_prot[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
		             V3_el_prot[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
		             V3_el_prot[2][1]=V4_el.Vect()+V3_prot_uncorr[2];

		             if( N_p_three!=0 && N_p_four!=0){
		   	             for(int count = 0; count < N_comb; count++)    { //looping through number of 2 proton combination out of 3 protons
	                      for(int j = 0;j < N_2p; j++)    {  //looping through number of 1 proton combination out of 2 protons

	//-----------------------------------------  4p to 3p->2->1  -----------------------------------------------------------------------

		                      h1_E_tot_4pto3p->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_E_rec_4pto3p->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_Erec_pperp_4321p->Fill(p_miss_perp_4pto3p[count][j],E_rec,P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_E_tot_4pto3p_fracfeed->Fill((E_cal_4pto3p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_E_rec_4pto3p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_pperp_W->Fill(W_var,p_miss_perp_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_theta0->Fill((V4_beam.Vect()).Angle(V3_el_prot[count][j])*TMath::RadToDeg(),-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);

		                      for(int i = 0; i < n_slice; i++)
			                    {
			                         if (p_miss_perp_4pto3p[count][j]<pperp_max[i] && p_miss_perp_4pto3p[count][j]>pperp_min[i]){
			                              h1_Etot_4pto3p_slice[i]->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
			                              h1_Erec_4pto3p_slice[i]->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
			                         }
			                    }
		                    } //end loop over N_2p
		                 } //end loop over N_comb : number of 2 proton combination out of 3 protons

		   //-----------------------------------------  4p to 3p->1p  -----------------------------------------------------------------------
		   	             for(int j = 0; j < N_3p; j++)    { //4p to 3p->1, looping through 1p out of 3p
                        //P_43pto1p doesnt have to be an array, one local variable here
		                     P_43pto1p[j]= N_p1[j]/N_p_three*(N_p4_p3[g]/N_p_four);
		                     h1_E_tot_43pto1p->Fill(E_cal_43pto1p[j], P_43pto1p[j]*histoweight);
		                     h1_E_rec_43pto1p->Fill(E_rec,P_43pto1p[j]*histoweight);
		                     h2_Erec_pperp_431p->Fill(p_miss_perp_43pto1p[j],E_rec,P_43pto1p[j]*histoweight);
		                     h1_E_tot_43pto1p_fracfeed->Fill((E_cal_43pto1p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_43pto1p[j]*histoweight);
		                     h1_E_rec_43pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_43pto1p[j]*histoweight);
		                     h2_pperp_W->Fill(W_var,p_miss_perp_43pto1p[j],P_43pto1p[j]*histoweight);
  		                   h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_43pto1p[j]*histoweight);
		                     h2_Ecal_Eqe->Fill(E_rec,E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
		                     h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
		                     h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_43pto1p[j],P_43pto1p[j]*histoweight);

                         for(int i = 0; i <n_slice; i++)
		                     {
			                        if (p_miss_perp_43pto1p[j]<pperp_max[i] && p_miss_perp_43pto1p[j]>pperp_min[i]){
			                             h1_Etot_43pto1p_slice[i]->Fill(E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
			                             h1_Erec_43pto1p_slice[i]->Fill(E_rec,P_43pto1p[j]*histoweight);
			                        }
		                     }
                    } // end loop over N_3p
		             }//end of N_p_three and N_p_four !=0
	            }//end of the loop through 3p combinations out of 4, g < Ncomb_4to3

              //still no pions or photons num_pi_phot == 0
	            int N_4to2=0;
      	      TVector3 V3p2[2],V3p2_uncorr[2];
	            double E_cal_4pto2p[2]={0};
              double p_miss_perp_4pto2p[2]={0};
              double P_4pto2p[2]={0};
              double N_two=0;

  //-----------------------------------------  4p to 2p->1  -----------------------------------------------------------------------
	            for(int ind1 = 0; ind1 < N_p4; ind1++){          //estimating the undetected 4p contribution to  2p
		              for(int ind2 = 0; ind2 < N_p4; ind2++){
		                  if(ind1!=ind2 && ind1 < ind2){

		                      V3p2[0]=V3_prot4_corr[ind1];
		                      V3p2[1]=V3_prot4_corr[ind2];
		                      V3p2_uncorr[0]=V3_prot4_uncorr[ind1];
		                      V3p2_uncorr[1]=V3_prot4_uncorr[ind2];

		                      rotation->prot2_rot_func( V3p2,V3p2_uncorr, V4_el,E_cal_4pto2p,p_miss_perp_4pto2p,  P_4pto2p, &N_two);

		                      if( N_two!=0  && N_p_four!=0){
		                          for(int j = 0; j < N_2p; j++)  {  //looping through  1 proton combination out of 2 protons

	                               h1_E_tot_4pto2p->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_E_rec_4pto2p->Fill(E_rec, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h2_Erec_pperp_421p->Fill( p_miss_perp_4pto2p[j],E_rec,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_E_tot_4pto2p_fracfeed->Fill((E_cal_4pto2p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_E_rec_4pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h2_pperp_W->Fill(W_var,p_miss_perp_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3p2_uncorr[j])*TMath::RadToDeg(),P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
                                 h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
		                             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);

	                               for(int i = 0; i < n_slice; i++)
	                               {
	                                  if (p_miss_perp_4pto2p[j]<pperp_max[i] && p_miss_perp_4pto2p[j]>pperp_min[i]){
	                                     h1_Etot_4pto2p_slice[i]->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
		                                   h1_Erec_4pto2p_slice[i]->Fill(E_rec, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
		                                }
		                             }
		                          } //end loop over N_2p
                          } //end if N_two!=0  && N_p_four!=0
		                      N_4to2= N_4to2+1;
                      } //end if ind1!=ind2 && ind1 < ind2
                  } //end loop over ind2
	            } //end loop over ind1

 //-----------------------------------------  4p to 1p  -----------------------------------------------------------------------
	            if( N_p_four!=0){
		              for(int j = 0; j < N_p4; j++)    {       //estimating the undetected 4p contribution to  1p
                      //P_4pto1p[j] doesnt have to be an array since it is only used here as a local variable
		                  P_4pto1p[j]= N_p4_p1[j]/N_p_four;
		                  h1_E_tot_4pto1p->Fill(E_cal_p4[j], P_4pto1p[j]*histoweight);
		                  h1_E_rec_4pto1p->Fill(E_rec,P_4pto1p[j]*histoweight);
		                  h2_Erec_pperp_41p->Fill(p_miss_perp_p4[j],E_rec, P_4pto1p[j]*histoweight);
		                  h1_E_tot_4pto1p_fracfeed->Fill((E_cal_p4[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto1p[j]*histoweight);
		                  h1_E_rec_4pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_4pto1p[j]*histoweight);
		                  h2_pperp_W->Fill(W_var,p_miss_perp_p4[j],-P_4pto1p[j]*histoweight);
		                  h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot4_uncorr[j])*TMath::RadToDeg(),-P_4pto1p[j]*histoweight);
		                  h2_Ecal_Eqe->Fill(E_rec,E_cal_p4[j],-P_4pto1p[j]*histoweight);
		                  h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_p4[j],-P_4pto1p[j]*histoweight);
		                  h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_p4[j],-P_4pto1p[j]*histoweight);

		                  for(int i = 0; i < n_slice; i++)
		                  {
		                      if (p_miss_perp_p4[j]<pperp_max[i] && p_miss_perp_p4[j]>pperp_min[i]){
                              h1_Etot_4pto1p_slice[i]->Fill(E_cal_p4[j],P_4pto1p[j]*histoweight);
			                        h1_Erec_4pto1p_slice[i]->Fill(E_rec,P_4pto1p[j]*histoweight);
		                      }
		                  }
                  } //end loop over N_p4
	            } // end if N_p_four!=0
        }//no pion statment ends

	   }//4 proton requirement (num_p == 4)

    //We are not looking for 4 Proton and 1 Pion events!

//No Protons here, Next 150 lines are for the inclusive events

     h1_E_rec->Fill(E_rec,WeightIncl);

     if(num_pi_phot ==0){
//     if(num_pi == 0){    //no photons for now F.H. 29.8.19
       h1_E_rec_0pi->Fill(E_rec,WeightIncl);
	     h1_E_rec_0pi_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],WeightIncl);
     }

	//----------------------------- e- ,1pi  -----------------------------------------

     if(num_pi_phot == 1){
//     if(num_pi == 1){    //no photons for now F.H. 29.8.19

        TVector3 V3_pi_corr;
        double P_undet=0;
        int charge = 0;
        double pion_acc_ratio = 0;

        if (num_pimi == 1 && num_pipl == 0) {
          charge = -1;
        }
        else if (num_pipl == 1 && num_pimi == 0) {
          charge = 1;
        }
        //radiation status is false if real photon
        else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
          charge = 0;
        }
        //skip radiation photon
        else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
         charge = 0;
        }
        else {  std::cout << "WARNING: 1pion events: No charge for one pion/photon could be assigned.  "  << std::endl; continue; }

        double SmearedPp = gRandom->Gaus(pf[index_pi[0]],reso_pi*pf[index_pi[0]]);
        double SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

        V3_pi_corr.SetXYZ(SmearedPp/pf[index_pi[0]] * pxf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pyf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pzf[index_pi[0]]);

        // apapadop
        double phi_pion = V3_pi_corr.Phi(); //in Radians
        V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
        phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

        double pion_theta = V3_pi_corr.Theta();
        double pion_mom_corr = V3_pi_corr.Mag();

        if (charge == 1) { //acceptance for pi plus
          pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
        }
        else if (charge == -1) {    //acceptance for pi minus. using electron acceptance map
          pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
        }
        else if (charge == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
          pion_acc_ratio = 1;
        }
        else { std::cout << "WARNING: 1 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

        // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
        // Pi_phot_fid_united with 0 is for Photon/Pi0
        //If any fiducial is not fullfilled continue to next event
        if ( ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr,  charge) ) )
        {
            continue;
        }


        rotation->pi1_rot_func( V3_pi_corr, charge, ec_radstat_n[0], &P_undet);

        double histoweight = pion_acc_ratio * WeightIncl;

        h1_E_rec_1pi_weight->Fill(E_rec,P_undet*histoweight);
	      h1_E_rec_1pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_undet*histoweight);

	      if(!ec_radstat_n[0])  h1_E_rec_1pi->Fill(E_rec,histoweight);
	      if(ec_num_n==1)       h2_phot_e_angle_Erec->Fill(E_rec,V3_pi_corr.Angle(V4_el.Vect())*TMath::RadToDeg());
	   }

	//----------------------------- e- ,2pi  -----------------------------------------

     if(num_pi_phot == 2){
  //   if(num_pi == 2){    //no photons for now F.H. 29.8.19

	      const int N_2pi = 2;
	      TVector3 V3_2pi_corr[N_2pi];
	      int q_pi2[N_2pi];
        bool radstat_pi2[N_2pi] = {false};
	      double P_1pi[N_2pi] = {0};
        double P_0pi = 0;
        double pion_acc_ratio[N_2pi] = {0};

        for (int i = 0; i < num_pi_phot; i++) {

            radstat_pi2[i] = ec_radstat_n[i];
            if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
              q_pi2[i] = 1;
            }
            else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
              q_pi2[i] = -1;
            }
            else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
              q_pi2[i] = 0;
            }
            else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
              q_pi2[i] = 0;
            }
            else {  std::cout << "WARNING: 2pion events: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

            double SmearedPp = gRandom->Gaus(pf[index_pi[i]],reso_pi*pf[index_pi[i]]);
            double SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

            V3_2pi_corr[i].SetXYZ(SmearedPp/pf[index_pi[i]] * pxf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pyf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pzf[index_pi[i]]);

            // apapadop
            double phi_pion = V3_2pi_corr[i].Phi(); //in Radians
            V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
            phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

            double pion_theta = V3_2pi_corr[i].Theta();
            double pion_mom_corr = V3_2pi_corr[i].Mag();

            if (q_pi2[i] == 1) { //acceptance for pi plus
                pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
            }
            else if (q_pi2[i] == -1) {    //acceptance for pi minus. using electron acceptance map
                pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
            }
            else if (q_pi2[i] == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
                pion_acc_ratio[i] = 1;
            }
            else { std::cout << "WARNING: 2 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

        }

      // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
      // Pi_phot_fid_united with 0 is for Photon/Pi0
      //If any fiducial is not fullfilled continue to next event
        if ( ( !Pi_phot_fid_united(fbeam_en, V3_2pi_corr[0],  q_pi2[0]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_2pi_corr[1],  q_pi2[1]) )  )
        {
            continue;
        }

        rotation->pi2_rot_func(V3_2pi_corr, q_pi2,radstat_pi2, &P_0pi,P_1pi);

        double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
        double histoweight = weight_pions * e_acc_ratio * wght/Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19

//----------------------------- e- ,2pi->0pi (-) -----------------------------------------
        h1_E_rec_2pi_weight->Fill(E_rec,(-P_0pi)*histoweight);
        h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*histoweight);
        h1_E_rec_20pi->Fill(E_rec,(P_0pi)*histoweight);
//----------------------------- e- ,2pi->1pi->0pi (+)  -----------------------------------------
        for(int k = 0; k < N_2pi; k++){ //loop over two pions
          h1_E_rec_2pi_weight->Fill(E_rec,P_1pi[k]*histoweight);
          h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[k]*histoweight);
          h1_E_rec_21pi->Fill(E_rec,(P_1pi[k])*histoweight);
        }
     } //end if for two pion events

//----------------------------- e- ,3pi  -----------------------------------------

	   if(num_pi_phot == 3){
  //   if(num_pi == 3){    //no photons for now F.H. 29.8.19

        const int N_3pi=3;
        const int N_2pi=2;
        TVector3 V3_3pi_corr[N_3pi];
        int q_pi3[N_3pi]={0};
        bool radstat_pi3[N_3pi]={false};
        double P_0pi = 0;
        double P_1pi[N_3pi]={0};
        double P_320pi[N_3pi]={0};
        double P_3210pi[N_3pi][N_2pi]={0};
        double pion_acc_ratio[N_3pi] = {0};

        for (int i = 0; i < num_pi_phot; i++) {

            radstat_pi3[i] = ec_radstat_n[i];
            if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
              q_pi3[i] = 1;
            }
            else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
              q_pi3[i] = -1;
            }
            else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
              q_pi3[i] = 0;
            }
            else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
              q_pi3[i] = 0;
            }
            else {  std::cout << "WARNING: 3pion events: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

            double SmearedPp = gRandom->Gaus(pf[index_pi[i]],reso_pi*pf[index_pi[i]]);
            double SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

            V3_3pi_corr[i].SetXYZ(SmearedPp/pf[index_pi[i]] * pxf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pyf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pzf[index_pi[i]]);

            // apapadop
            double phi_pion = V3_3pi_corr[i].Phi(); //in Radians
            V3_3pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
            phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

            double pion_theta = V3_3pi_corr[i].Theta();
            double pion_mom_corr = V3_3pi_corr[i].Mag();

            if (q_pi3[i] == 1) { //acceptance for pi plus
                pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
            }
            else if (q_pi3[i] == -1) {    //acceptance for pi minus. using electron acceptance map
                pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
            }
            else if (q_pi3[i] == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
                pion_acc_ratio[i] = 1;
            }
            else { std::cout << "WARNING: 3 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

        }

        // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
        // Pi_phot_fid_united with 0 is for Photons/Pi0
        //If any fiducial is not fullfilled continue to next event
        if ( ( !Pi_phot_fid_united(fbeam_en, V3_3pi_corr[0],  q_pi3[0]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_3pi_corr[1],  q_pi3[1]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_3pi_corr[2],  q_pi3[2]) ) )
        {
              continue;
        }

        rotation->pi3_rot_func( V3_3pi_corr, q_pi3,radstat_pi3, &P_0pi, P_1pi, P_320pi,P_3210pi);

        double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2];
        double histoweight = weight_pions * e_acc_ratio * wght/Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19

 //---------------------------3pi->0pi----------------------------------------------
        h1_E_rec_3pi_weight->Fill(E_rec,(-P_0pi)*histoweight);
        h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*histoweight);
        h1_E_rec_30pi->Fill(E_rec,(P_0pi)*histoweight);

        for(int h = 0; h < N_3pi; h++){ //loop over three pions

//---------------------------3pi->1pi->0pi----------------------------------------------
          h1_E_rec_3pi_weight->Fill(E_rec,P_1pi[h]*histoweight);
          h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[h]*histoweight);
          h1_E_rec_310pi->Fill(E_rec,(P_1pi[h])*histoweight);
 //---------------------------3pi->2pi->0pi----------------------------------------------
          h1_E_rec_3pi_weight->Fill(E_rec,P_320pi[h]*histoweight);
          h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_320pi[h]*histoweight);
          h1_E_rec_320pi->Fill(E_rec,(P_320pi[h])*histoweight);
//---------------------------3pi->2pi->1pi->0pi----------------------------------------------
          for(int g = 0; g < N_2pi; g++){ //loop over two pions
            h1_E_rec_3pi_weight->Fill(E_rec,(-P_3210pi[h][g])*histoweight);
            h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_3210pi[h][g])*histoweight);
            h1_E_rec_3210pi->Fill(E_rec,(P_3210pi[h][g])*histoweight);
          }
        }//end of 3pi loop
     }//end of 3pi requirement

//----------------------------- e- ,4pi  -----------------------------------------
	   if(num_pi_phot == 4){
  //   if(num_pi == 4){    //no photons for now F.H. 29.8.19

       const int N_4pi=4;
       TVector3 V3_4pi_corr[N_4pi];
       int q_pi4[N_4pi]={0};
       bool  radstat_pi4[N_4pi]={false};
       double P_0pi=0;
       double P_410pi=0;
       double P_420pi=0;
       double P_4210pi=0;
       double P_430pi=0;
       double P_4310pi=0;
       double P_4320pi=0;
       double P_43210pi=0;
       double pion_acc_ratio[N_4pi] = {0};

       for (int i = 0; i < num_pi_phot; i++) {

           radstat_pi4[i] = ec_radstat_n[i];
           if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
             q_pi4[i] = 1;
           }
           else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
             q_pi4[i] = -1;
           }
           else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
             q_pi4[i] = 0;
           }
           else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
             q_pi4[i] = 0;
           }
           else {  std::cout << "WARNING: 4pion events: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

           double SmearedPp = gRandom->Gaus(pf[index_pi[i]],reso_pi*pf[index_pi[i]]);
    // not used       double SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

           V3_4pi_corr[i].SetXYZ(SmearedPp/pf[index_pi[i]] * pxf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pyf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pzf[index_pi[i]]);

           // apapadop
           double phi_pion = V3_4pi_corr[i].Phi(); //in Radians
           V3_4pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
           phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

           double pion_theta = V3_4pi_corr[i].Theta();
           double pion_mom_corr = V3_4pi_corr[i].Mag();

           if (q_pi4[i] == 1) { //acceptance for pi plus
               pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
           }
           else if (q_pi4[i] == -1) {    //acceptance for pi minus. using electron acceptance map
               pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
           }
           else if (q_pi4[i] == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
               pion_acc_ratio[i] = 1;
           }
           else { std::cout << "WARNING: 4 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

       }

       // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
        // Pi_phot_fid_united with 0 is for Photon/Pi0
       //If any fiducial is not fullfilled continue to next event
       if ( ( !Pi_phot_fid_united(fbeam_en, V3_4pi_corr[0],  q_pi4[0]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_4pi_corr[1],  q_pi4[1]) )
         || ( !Pi_phot_fid_united(fbeam_en, V3_4pi_corr[2],  q_pi4[2]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_4pi_corr[3],  q_pi4[3]) ) )
       {
             continue;
       }

       rotation->pi4_rot_func(V3_4pi_corr, q_pi4, radstat_pi4, &P_0pi,&P_410pi,&P_420pi,&P_4210pi,&P_430pi,&P_4310pi,&P_4320pi,&P_43210pi);

       double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2] * pion_acc_ratio[3];
       double histoweight = weight_pions * e_acc_ratio * wght/Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19


 //---------------------------4pi->0pi----------------------------------------------
//why is it here not split like for 3pi case, sum over all weights is done here F.H 04.08.19
       h1_E_rec_4pi_weight->Fill(E_rec,(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight);
       h1_E_rec_4pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight);
       h1_E_rec_40pi->Fill(E_rec,(P_0pi)*histoweight);
//---------------------------4pi->1pi->0pi----------------------------------------------
       h1_E_rec_410pi->Fill(E_rec,(P_410pi)*histoweight);
//---------------------------4pi->2pi->0pi----------------------------------------------
       h1_E_rec_420pi->Fill(E_rec,(P_420pi)*histoweight);
//---------------------------4pi->2pi->1pi->0pi----------------------------------------------
       h1_E_rec_4210pi->Fill(E_rec,(P_4210pi)*histoweight);
//---------------------------4pi->3pi->0pi----------------------------------------------
       h1_E_rec_430pi->Fill(E_rec,(P_430pi)*histoweight);
//---------------------------4pi->3pi->1pi->0pi----------------------------------------------
       h1_E_rec_4310pi->Fill(E_rec,(P_4310pi)*histoweight);
//---------------------------4pi->3pi->2pi->0pi----------------------------------------------
       h1_E_rec_4320pi->Fill(E_rec,(P_4320pi)*histoweight);
//---------------------------4pi->3pi->2pi->1pi->0pi----------------------------------------------
       h1_E_rec_43210pi->Fill(E_rec,(P_43210pi)*histoweight);

     }//end of 4 pi/photon requirement


   //------------------------------------------requiring there to be a proton -------------------------------------
    //Events with exactly one proton
     if( num_p == 1)
     {

       double SmearedPp;
       double SmearedEp;


     //Vector for proton without momentum smearing
       TLorentzVector V4_prot_uncorr(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
       TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();

     //Smearing of Proton
       SmearedPp = gRandom->Gaus(pf[index_p[0]],reso_p*pf[index_p[0]]);
       SmearedEp = sqrt( SmearedPp*SmearedPp + m_prot * m_prot );

      //Vector for proton with momentum smearing
       TVector3 V3_prot_corr;
       TLorentzVector V4_prot_corr;
       V3_prot_corr.SetXYZ(SmearedPp/pf[index_p[0]] * pxf[index_p[0]],SmearedPp/pf[index_p[0]] * pyf[index_p[0]],SmearedPp/pf[index_p[0]] * pzf[index_p[0]]);
       V4_prot_corr.SetPxPyPzE(SmearedPp/pf[index_p[0]] * pxf[index_p[0]],SmearedPp/pf[index_p[0]] * pyf[index_p[0]],SmearedPp/pf[index_p[0]] * pzf[index_p[0]],SmearedEp);

       // apapadop
       double phi_prot = V3_prot_corr.Phi(); //in Radians
       V3_prot_corr.SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
       phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

       if (!PFiducialCut(fbeam_en, V3_prot_corr) ) { continue; } // Proton theta & phi fiducial cuts

       //Proton kinematic variables
       double p_theta = V3_prot_corr.Theta();
       double prot_mom_corr = V3_prot_corr.Mag();
       //Proton weight
       double p_acc_ratio = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);

       h1_prot_mom->Fill(prot_mom_corr);
       h1_prot_mom_ratio->Fill(prot_mom_corr/pf[index_p[0]]);

       TLorentzVector V4_prot_el_tot = V4_prot_corr + V4_el;

       double p_perp_tot = TMath::Sqrt(V4_prot_el_tot.Px()*V4_prot_el_tot.Px() + V4_prot_el_tot.Py()*V4_prot_el_tot.Py());
       double E_tot = V4_el.E() + V4_prot_corr.E() -m_prot + bind_en[ftarget];

//These Histograms are events with 1 electron and  1 proton and multiple pions
       double histoweight_inc = p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;

       h1_E_tot->Fill(E_tot,histoweight_inc);
       h1_E_rec_1prot->Fill(E_rec,histoweight_inc);
       h1_E_tot_1prot->Fill(E_tot,histoweight_inc);
       h2_Erec_pperp->Fill(p_perp_tot,E_rec,histoweight_inc);

       //---------------------------------- 1p 0pi   ----------------------------------------------
       if(num_pi_phot == 0){
  //     if(num_pi == 0){    //no photons for now F.H. 29.8.19

          double histoweight = p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;

          h2_Erec_pperp_newcut2->Fill(p_perp_tot,E_rec,histoweight);
     	    h1_E_rec_cut2_new->Fill(E_rec,histoweight);
     	    h1_E_tot_cut2->Fill(E_tot,histoweight);
     	    h1_E_tot_cut2_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],histoweight);
     	    h1_E_rec_cut2_new_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],histoweight);
     	    h2_pperp_W->Fill(W_var,p_perp_tot,histoweight);
     	    h1_theta0->Fill((V4_beam.Vect()).Angle(V4_prot_el_tot.Vect()) *TMath::RadToDeg(),histoweight);
     	    h2_Ecal_Eqe->Fill(E_rec,E_tot,histoweight);
     	    h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,histoweight);
     	    h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,histoweight);

          for(int i = 0; i < n_slice; i++)
          {
            if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
              h1_Etot_Npi0[i]->Fill(E_tot,histoweight);
              h1_Erec_Npi0[i]->Fill(E_rec,histoweight);
            }
          }

          if (p_perp_tot < 0.2){
            h1_E_rec_cut005_newcut3->Fill(E_rec,histoweight);
     	      h2_Erec_pperp_cut3->Fill(p_perp_tot,E_rec,histoweight);
     	    }
       }//num pi=0
       //---------------------------------- 1p 1pi   ----------------------------------------------
       if(num_pi_phot == 1){
    //   if(num_pi == 1){    //no photons for now F.H. 29.8.19

          double N_piphot_det;
          double N_piphot_undet;
          TVector3 V3_pi_corr;
          int charge = 0;
          double pion_acc_ratio = 0;

          if (num_pimi == 1 && num_pipl == 0) {
            charge = -1;
          }
          else if (num_pipl == 1 && num_pimi == 0) {
            charge = 1;
          }
          //radiation status is false if real photon
          else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
            charge = 0;
          }
          //skip radiation photon
          else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
           charge = 0;
          }
          else {  std::cout << "WARNING: 1pion events: No charge for one pion/photon could be assigned.  "  << std::endl; continue; }

          SmearedPp = gRandom->Gaus(pf[index_pi[0]],reso_pi*pf[index_pi[0]]);
          SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

          V3_pi_corr.SetXYZ(SmearedPp/pf[index_pi[0]] * pxf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pyf[index_pi[0]],SmearedPp/pf[index_pi[0]] * pzf[index_pi[0]]);

          // apapadop
          double phi_pion = V3_pi_corr.Phi(); //in Radians
          V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
          phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

          double pion_theta = V3_pi_corr.Theta();
          double pion_mom_corr = V3_pi_corr.Mag();

          if (charge == 1) { //acceptance for pi plus
            pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
          }
          else if (charge == -1) {    //acceptance for pi minus. using electron acceptance map
            pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
          }
          else if (charge == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
            pion_acc_ratio = 1;
          }
          else { std::cout << "WARNING: 1 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

          // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
          // Pi_phot_fid_united with 0 is for Photon/Pi0
          //If any fiducial is not fullfilled continue to next event
          if ( ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr,  charge) ) )
          {
              continue;
          }

          rotation->prot1_pi1_rot_func(V3_prot_uncorr,V3_pi_corr, charge, ec_radstat_n[0], &N_piphot_det,&N_piphot_undet);

          double histoweight = pion_acc_ratio * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; //1proton, 1 Pion, 1 electron acceptance, GENIE weight and Mott

          if(N_piphot_det!=0){

             h1_E_rec_undetfactor->Fill(E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
             h1_E_tot_undetfactor->Fill(E_tot,(N_piphot_undet/N_piphot_det)*histoweight);
             h2_Erec_pperp_1p1pi->Fill(p_perp_tot,E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
             h1_E_rec_undetfactor_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(N_piphot_undet/N_piphot_det)*histoweight);
             h1_E_tot_undetfactor_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],(N_piphot_undet/N_piphot_det)*histoweight);
             h2_pperp_W->Fill(W_var,p_perp_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-(N_piphot_undet/N_piphot_det)*histoweight);
             h2_Ecal_Eqe->Fill(E_rec,E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);

             for(int i = 0; i < n_slice; i++)
             {
               if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
                 h1_Etot_bkgd_pipl_pimi_fact[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*histoweight);
                 h1_Erec_bkgd_pipl_pimi_new_fact[i]->Fill(E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
                 h1_Etot_bkgd_pipl_pimi_fact_pipl[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*histoweight);
                 h1_Etot_Npi1[i]->Fill(E_tot,histoweight);
                 h1_Erec_Npi1[i]->Fill(E_rec,histoweight);
               }
             }
          } //end of N_piphot_det!=0

          h1_E_rec_cutpi1_piplpimi->Fill(E_rec,histoweight);
          h1_E_tot_cutpi1_piplpimi->Fill(E_tot,histoweight);

       }//end of 1p 1pi requirement

 //---------------------------------- 1p 2pi   ----------------------------------------------
       if(num_pi_phot == 2){
  //     if(num_pi == 2){    //no photons for now F.H. 29.8.19

         const int N_2pi=2;
         TVector3 V3_2pi_corr[N_2pi],V3_2pi_rot[N_2pi],V3_p_rot;
         int q_pi2[N_2pi];
         bool radstat_pi2[N_2pi]={false};
         double P_1p0pi=0;
         double P_1p1pi[N_2pi]={0};

         double pion_acc_ratio[N_2pi] = {0};

         for (int i = 0; i < num_pi_phot; i++) {

             radstat_pi2[i] = ec_radstat_n[i];
             if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
               q_pi2[i] = 1;
             }
             else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
               q_pi2[i] = -1;
             }
             else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
               q_pi2[i] = 0;
             }
             else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
               q_pi2[i] = 0;
             }
             else {  std::cout << "WARNING: 1Proton 2pion events: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

             SmearedPp = gRandom->Gaus(pf[index_pi[i]],reso_pi*pf[index_pi[i]]);
             SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

             V3_2pi_corr[i].SetXYZ(SmearedPp/pf[index_pi[i]] * pxf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pyf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pzf[index_pi[i]]);

             // apapadop
             double phi_pion = V3_2pi_corr[i].Phi(); //in Radians
             V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
	           phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

             double pion_theta = V3_2pi_corr[i].Theta();
             double pion_mom_corr = V3_2pi_corr[i].Mag();

             if (q_pi2[i] == 1) { //acceptance for pi plus
                 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
             }
             else if (q_pi2[i] == -1) {    //acceptance for pi minus. using electron acceptance map
                 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
             }
             else if (q_pi2[i] == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
                 pion_acc_ratio[i] = 1;
             }
             else { std::cout << "WARNING: 1 Proton 2 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

         }

         // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
         // Pi_phot_fid_united with 0 for Photons/Pi0
         //If any fiducial is not fullfilled continue to next event
         if ( ( !Pi_phot_fid_united(fbeam_en, V3_2pi_corr[0],  q_pi2[0]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_2pi_corr[1],  q_pi2[1]) ) )
         {
               continue;
         }

         rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi_corr,q_pi2,radstat_pi2,&P_1p0pi,P_1p1pi);

         double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
         double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; //1proton, 2 Pion, 1 electron acceptance, GENIE weight and Mott

 //---------------------------------- 1p 2pi->1p1pi   ----------------------------------------------

         for(int z = 0; z < N_2pi; z++){  //to consider 2 diff. 1pi states

           h1_E_tot_1p2pi->Fill(E_tot,P_1p1pi[z]*histoweight);
           h1_E_rec_1p2pi->Fill(E_rec,P_1p1pi[z]*histoweight);
           h2_Erec_pperp_1p2pi_1p1pi->Fill(p_perp_tot,E_rec,P_1p1pi[z]*histoweight);
           h1_E_tot_1p2pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
           h1_E_rec_1p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p1pi[z]*histoweight);
           h2_pperp_W->Fill(W_var,p_perp_tot,P_1p1pi[z]*histoweight);
           h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p1pi[z]*histoweight);
           h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p1pi[z]*histoweight);
           h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p1pi[z]*histoweight);
           h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p1pi[z]*histoweight);

           for(int i = 0; i < n_slice; i++){
	            if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	               h1_Etot_bkgd_1p2pi[i]->Fill(E_tot,P_1p1pi[z]*histoweight);
	               h1_Erec_bkgd_1p2pi[i]->Fill(E_rec,P_1p1pi[z]*histoweight);
	            }
           }
         } //end loop over N_2pi

 //---------------------------------- 1p 2pi->1p0pi   ----------------------------------------------

         h1_E_tot_1p2pi_1p0pi->Fill(E_tot,P_1p0pi*histoweight);
         h1_E_rec_1p2pi_1p0pi->Fill(E_rec,P_1p0pi*histoweight);
         h2_Erec_pperp_1p2pi_1p0pi->Fill(p_perp_tot,E_rec,P_1p0pi*histoweight);
         h1_E_tot_1p2pi_1p0pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p0pi*histoweight);
         h1_E_rec_1p2pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p0pi*histoweight);
         h2_pperp_W->Fill(W_var,p_perp_tot,-P_1p0pi*histoweight);
         h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-P_1p0pi*histoweight);
         h2_Ecal_Eqe->Fill(E_rec,E_tot,-P_1p0pi*histoweight);
         h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-P_1p0pi*histoweight);
         h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-P_1p0pi*histoweight);

         for(int i = 0; i < n_slice; i++){
           if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	            h1_Etot_bkgd_1p2pi_1p0pi[i]->Fill(E_tot,P_1p0pi*histoweight);
	            h1_Erec_bkgd_1p2pi_1p0pi[i]->Fill(E_rec,P_1p0pi*histoweight);
           }
         }
       }//1p 2pi statetment ends

 //---------------------------------- 1p 3pi   ----------------------------------------------
       if(num_pi_phot == 3){
    //   if(num_pi == 3){    //no photons for now F.H. 29.8.19

         const int N_3pi=3;
         TVector3 V3_3pi_corr[N_3pi],V3_3pi_rot[N_3pi],V3_p_rot;
         int q_pi3[N_3pi];
         bool radstat_pi3[N_3pi]={false};
         double P_1p3pi = 0;
         double pion_acc_ratio[N_3pi] = {0};

         for (int i = 0; i < num_pi; i++) {

             radstat_pi3[i] = ec_radstat_n[i];
             if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
               q_pi3[i] = 1;
             }
             else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
               q_pi3[i] = -1;
             }
             else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
               q_pi3[i] = 0;
             }
             else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
               q_pi3[i] = 0;
             }
             else {  std::cout << "WARNING: 3pion events: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

             SmearedPp = gRandom->Gaus(pf[index_pi[i]],reso_pi*pf[index_pi[i]]);
             SmearedEp = sqrt( SmearedPp*SmearedPp + m_pion * m_pion );

             V3_3pi_corr[i].SetXYZ(SmearedPp/pf[index_pi[i]] * pxf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pyf[index_pi[i]],SmearedPp/pf[index_pi[i]] * pzf[index_pi[i]]);

             // apapadop
             double phi_pion = V3_3pi_corr[i].Phi(); //in Radians
             V3_3pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
             phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

             double pion_theta = V3_3pi_corr[i].Theta();
             double pion_mom_corr = V3_3pi_corr[i].Mag();

             if (q_pi3[i] == 1) { //acceptance for pi plus
                 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
             }
             else if (q_pi3[i] == -1) {    //acceptance for pi minus. using electron acceptance map
                 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
             }
             else if (q_pi3[i] == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
                 pion_acc_ratio[i] = 1;
             }
             else { std::cout << "WARNING: 3 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

         }

         // Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus.
         // Pi_phot_fid_united with 0 is for Photons/Pi0
         //If any fiducial is not fullfilled continue to next event
         if ( ( !Pi_phot_fid_united(fbeam_en, V3_3pi_corr[0],  q_pi3[0]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_3pi_corr[1],  q_pi3[1]) ) || ( !Pi_phot_fid_united(fbeam_en, V3_3pi_corr[2],  q_pi3[2]) ) )
         {
               continue;
         }

         rotation->prot1_pi3_rot_func(V3_prot_uncorr,V3_3pi_corr,q_pi3,radstat_pi3,&P_1p3pi);

         double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2];
         double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; //1proton, 3 Pion, 1 electron acceptance, GENIE weight and Mott


 //---------------------------------- 1p 3pi->1p 0pi  total ?? F.H. 08/13/19 check logic here compared to 1p 2pi case ----------------------------------------------
         h1_E_tot_1p3pi->Fill(E_tot,P_1p3pi*histoweight);
         h1_E_rec_1p3pi->Fill(E_rec,P_1p3pi*histoweight);
         h2_Erec_pperp_1p3pi->Fill(p_perp_tot,E_rec,P_1p3pi*histoweight);
         h1_E_tot_1p3pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi*histoweight);
         h1_E_rec_1p3pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p3pi*histoweight);
         h2_pperp_W->Fill(W_var,p_perp_tot,P_1p3pi*histoweight);
         h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p3pi*histoweight);
         h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p3pi*histoweight);
         h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p3pi*histoweight);
         h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p3pi*histoweight);

         for(int i = 0; i < n_slice; i++){
           if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
             h1_Etot_bkgd_1p3pi[i]->Fill(E_tot,P_1p3pi*histoweight);
             h1_Erec_bkgd_1p3pi[i]->Fill(E_rec,P_1p3pi*histoweight);
           }
         }
       }//end of 1p 3pi requirement


     } // 1proton ends


  } //end of event loop (jentry)

  gStyle->SetOptFit(1);

  for(int i = 0; i <= n_slice-1; i++)
  {
  //------------------------------------using the ratio of the pi- to pi+  --------------------------------------

	   h1_Etot_piplpimi_subtruct_fact[i] = (TH1F*)  h1_Etot_Npi0[i]->Clone(Form("h1_Etot_piplpimi_subtruct_fact_%d",i+1));
	   h1_Etot_piplpimi_subtruct_fact[i]->Add(h1_Etot_bkgd_pipl_pimi_fact[i],-1);
	   h1_Erec_piplpimi_subtruct_fact[i]=(TH1F*)  h1_Erec_Npi0[i]->Clone(Form("h1_Erec_piplpimi_subtruct_fact_%d",i+1));
	   h1_Erec_piplpimi_subtruct_fact[i]->Add(h1_Erec_bkgd_pipl_pimi_new_fact[i],-1);

 //------------------------------------subtracting 2p contribution from 1p events  --------------------------------------

	   h1_Etot_p_bkgd_slice_sub[i]=(TH1F*) h1_Etot_piplpimi_subtruct_fact[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub[i]->Add(h1_Etot_p_bkgd_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub[i]=(TH1F*) h1_Erec_piplpimi_subtruct_fact[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub[i]->Add(h1_Erec_p_bkgd_slice[i],-1);

 //------------------------------------undetected 3 to 2 proton subtraction --------------------------------------

	   h1_Etot_p_bkgd_slice_sub32[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub32_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub32[i]->Add(h1_Etot_3pto2p_slice[i]);
	   h1_Erec_p_bkgd_slice_sub32[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub32_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub32[i]->Add(h1_Erec_3pto2p_slice[i]);

 //------------------------------------undetected 3 to 1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub31[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub32[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub31_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub31[i]->Add(h1_Etot_3pto1p_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub31[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub32[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub31_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub31[i]->Add(h1_Erec_3pto1p_slice[i],-1);

 //------------------------------------undetected 4 to 3->2->1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub43[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub31[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub43_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub43[i]->Add(h1_Etot_4pto3p_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub43[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub31[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub43_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub43[i]->Add(h1_Erec_4pto3p_slice[i],-1);

 //------------------------------------undetected 4 to 3->1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub431[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub43[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub431_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub431[i]->Add(h1_Etot_43pto1p_slice[i]);
	   h1_Erec_p_bkgd_slice_sub431[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub43[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub431_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub431[i]->Add(h1_Erec_43pto1p_slice[i]);

 //------------------------------------undetected 4 to 2 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub42[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub431[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub42_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub42[i]->Add(h1_Etot_4pto2p_slice[i]);
	   h1_Erec_p_bkgd_slice_sub42[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub431[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub42_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub42[i]->Add(h1_Erec_4pto2p_slice[i]);

 //------------------------------------undetected 4 to 1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub41[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub42[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub41_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub41[i]->Add(h1_Etot_4pto1p_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub41[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub42[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub41_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub41[i]->Add(h1_Erec_4pto1p_slice[i],-1);

//------------------------------------undetected 1p 2pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub1p2pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub41[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub1p2pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub1p2pi[i]->Add(h1_Etot_bkgd_1p2pi[i]);
	   h1_Erec_p_bkgd_slice_sub1p2pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub41[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub1p2pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub1p2pi[i]->Add(h1_Erec_bkgd_1p2pi[i]);

//------------------------------------undetected 1p 2pi-> 1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub1p2pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub1p2pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub1p2pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub1p2pi_0pi[i]->Add(h1_Etot_bkgd_1p2pi_1p0pi[i],-1);
	   h1_Erec_p_bkgd_slice_sub1p2pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub1p2pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub1p2pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub1p2pi_0pi[i]->Add(h1_Erec_bkgd_1p2pi_1p0pi[i],-1);

//------------------------------------undetected 1p 3pi-> 1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub1p3pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub1p2pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub1p3pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub1p3pi_0pi[i]->Add(h1_Etot_bkgd_1p3pi[i]);
	   h1_Erec_p_bkgd_slice_sub1p3pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub1p2pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub1p3pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub1p3pi_0pi[i]->Add(h1_Erec_bkgd_1p3pi[i]);

//------------------------------------undetected 2p 2pi-> 1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p2pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub1p3pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p2pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p2pi_0pi[i]->Add(h1_Etot_p_bkgd_slice_2p2pi[i]);
	   h1_Erec_p_bkgd_slice_sub2p2pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub1p3pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p2pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p2pi_0pi[i]->Add(h1_Erec_p_bkgd_slice_2p2pi[i]);

//------------------------------------undetected 3p 1pi->1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub3p1pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub2p2pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub3p1pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub3p1pi_0pi[i]->Add(h1_Etot_3p1pi_slice[i]);
	   h1_Erec_p_bkgd_slice_sub3p1pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub2p2pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub3p1pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub3p1pi_0pi[i]->Add(h1_Erec_3p1pi_slice[i]);

//------------------------------------undetected 2p 1pi ->2p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p1pi_2p[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub3p1pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p1pi_2p_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p1pi_2p[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[i]);
	   h1_Erec_p_bkgd_slice_sub2p1pi_2p[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub3p1pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p1pi_2p_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p1pi_2p[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[i]);

//------------------------------------undetected 2p 1pi ->1p 1pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub2p1pi_2p[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p1pi_1p_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[i]);
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub2p1pi_2p[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p1pi_1p_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[i]);

//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[i],-1);
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[i],-1);

  }

 //------------------------------------fractional energy reconstruction plots --------------------------------------

  //------------------------------------using the ratio of the pi- to pi+  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_factor =(TH1F*)  h1_E_rec_cut2_new->Clone("h_Erec_subtruct_piplpimi_factor");
  h_Erec_subtruct_piplpimi_factor->Add(h1_E_rec_undetfactor,-1);

  TH1F *h_Etot_subtruct_piplpimi_factor=(TH1F*)  h1_E_tot_cut2->Clone("h_Etot_subtruct_piplpimi_factor");
  h_Etot_subtruct_piplpimi_factor->Add(h1_E_tot_undetfactor,-1);

  TH2F *h2_Erec_pperp_1p1pisub=(TH2F*) h2_Erec_pperp_newcut2->Clone("h2_Erec_pperp_1p1pisub");
  h2_Erec_pperp_1p1pisub->Add(h2_Erec_pperp_1p1pi,-1);

  TH1F *h_Erec_subtruct_piplpimi_factor_fracfeed =(TH1F*)  h1_E_rec_cut2_new_fracfeed->Clone("h_Erec_subtruct_piplpimi_factor_fracfeed");
  h_Erec_subtruct_piplpimi_factor_fracfeed->Add(h1_E_rec_undetfactor_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_factor_fracfeed=(TH1F*)  h1_E_tot_cut2_fracfeed->Clone("h_Etot_subtruct_piplpimi_factor_fracfeed");
  h_Etot_subtruct_piplpimi_factor_fracfeed->Add(h1_E_tot_undetfactor_fracfeed,-1);

 //-----------------------------------undetected 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_prot=(TH1F*)  h_Erec_subtruct_piplpimi_factor->Clone("h_Erec_subtruct_piplpimi_prot");
  h_Erec_subtruct_piplpimi_prot->Add(h1_E_rec_p_bkgd,-1);

  TH1F *h_Etot_subtruct_piplpimi_prot=(TH1F*)  h_Etot_subtruct_piplpimi_factor->Clone("h_Etot_subtruct_piplpimi_prot");
  h_Etot_subtruct_piplpimi_prot->Add(h1_E_tot_p_bkgd,-1);

  TH2F *h2_Erec_pperp_2psub=(TH2F*) h2_Erec_pperp_1p1pisub->Clone("h2_Erec_pperp_2psub");
  h2_Erec_pperp_2psub->Add(h2_Erec_pperp_2p,-1);

  TH1F *h_Erec_subtruct_piplpimi_prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_factor_fracfeed->Clone("h_Erec_subtruct_piplpimi_prot_fracfeed");
  h_Erec_subtruct_piplpimi_prot_fracfeed->Add(h1_E_rec_p_bkgd_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_factor_fracfeed->Clone("h_Etot_subtruct_piplpimi_prot_fracfeed");
  h_Etot_subtruct_piplpimi_prot_fracfeed->Add(h1_E_tot_p_bkgd_fracfeed,-1);

 //-----------------------------------undetected 3 to 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_32prot=(TH1F*)  h_Erec_subtruct_piplpimi_prot->Clone("h_Erec_subtruct_piplpimi_32prot");
  h_Erec_subtruct_piplpimi_32prot->Add(h1_E_rec_3pto2p);

  TH1F *h_Etot_subtruct_piplpimi_32prot=(TH1F*)  h_Etot_subtruct_piplpimi_prot->Clone("h_Etot_subtruct_piplpimi_32prot");
  h_Etot_subtruct_piplpimi_32prot->Add(h1_E_tot_3pto2p);

  TH2F *h2_Erec_pperp_32psub=(TH2F*) h2_Erec_pperp_2psub->Clone("h2_Erec_pperp_32psub");
  h2_Erec_pperp_32psub->Add(h2_Erec_pperp_321p);

  TH1F *h_Erec_subtruct_piplpimi_32prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_32prot_fracfeed");
  h_Erec_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_rec_3pto2p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_32prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_32prot_fracfeed");
  h_Etot_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_tot_3pto2p_fracfeed);

 //-----------------------------------undetected 3 to 1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_31prot=(TH1F*)  h_Erec_subtruct_piplpimi_32prot->Clone("h_Erec_subtruct_piplpimi_31prot");
  h_Erec_subtruct_piplpimi_31prot->Add(h1_E_rec_3pto1p,-1);

  TH1F *h_Etot_subtruct_piplpimi_31prot=(TH1F*)  h_Etot_subtruct_piplpimi_32prot->Clone("h_Etot_subtruct_piplpimi_31prot");
  h_Etot_subtruct_piplpimi_31prot->Add(h1_E_tot_3pto1p,-1);

  TH2F *h2_Erec_pperp_31psub=(TH2F*) h2_Erec_pperp_32psub->Clone("h2_Erec_pperp_31psub");
  h2_Erec_pperp_31psub->Add(h2_Erec_pperp_31p,-1);

  TH1F *h_Erec_subtruct_piplpimi_31prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_32prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_31prot_fracfeed");
  h_Erec_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_rec_3pto1p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_31prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_32prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_31prot_fracfeed");
  h_Etot_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_tot_3pto1p_fracfeed,-1);

 //-----------------------------------undetected 4 to 3->2->1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_43prot=(TH1F*)  h_Erec_subtruct_piplpimi_31prot->Clone("h_Erec_subtruct_piplpimi_43prot");
  h_Erec_subtruct_piplpimi_43prot->Add(h1_E_rec_4pto3p,-1);

  TH1F *h_Etot_subtruct_piplpimi_43prot=(TH1F*)  h_Etot_subtruct_piplpimi_31prot->Clone("h_Etot_subtruct_piplpimi_43prot");
  h_Etot_subtruct_piplpimi_43prot->Add(h1_E_tot_4pto3p,-1);

  TH2F *h2_Erec_pperp_43psub=(TH2F*) h2_Erec_pperp_31psub->Clone("h2_Erec_pperp_43psub");
  h2_Erec_pperp_43psub->Add(h2_Erec_pperp_4321p,-1);

  TH1F *h_Erec_subtruct_piplpimi_43prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_31prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_43prot_fracfeed");
  h_Erec_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_rec_4pto3p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_43prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_31prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_43prot_fracfeed");
  h_Etot_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_tot_4pto3p_fracfeed,-1);

 //-----------------------------------undetected 4 to 3->1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_431prot=(TH1F*)  h_Erec_subtruct_piplpimi_43prot->Clone("h_Erec_subtruct_piplpimi_431prot");
  h_Erec_subtruct_piplpimi_431prot->Add(h1_E_rec_43pto1p);

  TH1F *h_Etot_subtruct_piplpimi_431prot=(TH1F*)  h_Etot_subtruct_piplpimi_43prot->Clone("h_Etot_subtruct_piplpimi_431prot");
  h_Etot_subtruct_piplpimi_431prot->Add(h1_E_tot_43pto1p);

  TH2F *h2_Erec_pperp_431psub=(TH2F*) h2_Erec_pperp_43psub->Clone("h2_Erec_pperp_431psub");
  h2_Erec_pperp_431psub->Add(h2_Erec_pperp_431p);

  TH1F *h_Erec_subtruct_piplpimi_431prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_43prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_431prot_fracfeed");
  h_Erec_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_rec_43pto1p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_431prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_43prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_431prot_fracfeed");
  h_Etot_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_tot_43pto1p_fracfeed);

//-----------------------------------undetected 4 to 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_42prot=(TH1F*)  h_Erec_subtruct_piplpimi_431prot->Clone("h_Erec_subtruct_piplpimi_42prot");
  h_Erec_subtruct_piplpimi_42prot->Add(h1_E_rec_4pto2p);

  TH1F *h_Etot_subtruct_piplpimi_42prot=(TH1F*) h_Etot_subtruct_piplpimi_431prot->Clone("h_Etot_subtruct_piplpimi_42prot");
  h_Etot_subtruct_piplpimi_42prot->Add(h1_E_tot_4pto2p);

  TH2F *h2_Erec_pperp_42psub=(TH2F*) h2_Erec_pperp_431psub->Clone("h2_Erec_pperp_42psub");
  h2_Erec_pperp_42psub->Add(h2_Erec_pperp_421p);

  TH1F *h_Erec_subtruct_piplpimi_42prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_431prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_42prot_fracfeed");
  h_Erec_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_rec_4pto2p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_42prot_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_431prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_42prot_fracfeed");
  h_Etot_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_tot_4pto2p_fracfeed);

 //-----------------------------------undetected 4 to 1 proton subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_41prot=(TH1F*)  h_Erec_subtruct_piplpimi_42prot->Clone("h_Erec_subtruct_piplpimi_41prot");
  h_Erec_subtruct_piplpimi_41prot->Add(h1_E_rec_4pto1p,-1);

  TH1F *h_Etot_subtruct_piplpimi_41prot=(TH1F*)  h_Etot_subtruct_piplpimi_42prot->Clone("h_Etot_subtruct_piplpimi_41prot");
  h_Etot_subtruct_piplpimi_41prot->Add(h1_E_tot_4pto1p,-1);

  TH2F *h2_Erec_pperp_41psub=(TH2F*) h2_Erec_pperp_42psub->Clone("h2_Erec_pperp_41psub");
  h2_Erec_pperp_41psub->Add(h2_Erec_pperp_41p,-1);

  TH1F *h_Erec_subtruct_piplpimi_41prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_42prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_41prot_fracfeed");
  h_Erec_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_rec_4pto1p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_41prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_42prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_41prot_fracfeed");
  h_Etot_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_tot_4pto1p_fracfeed,-1);

//------------------------------------undetected 1p 2pi ->1 p1pi ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_1p2pi=(TH1F*)  h_Erec_subtruct_piplpimi_41prot->Clone("h_Erec_subtruct_piplpimi_1p2pi");
  h_Erec_subtruct_piplpimi_1p2pi->Add(h1_E_rec_1p2pi);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi=(TH1F*)  h_Etot_subtruct_piplpimi_41prot->Clone("h_Etot_subtruct_piplpimi_1p2pi");
  h_Etot_subtruct_piplpimi_1p2pi->Add(h1_E_tot_1p2pi);

  TH2F *h2_Erec_pperp_sub_1p2pi_1p1pi=(TH2F*) h2_Erec_pperp_41psub->Clone("h2_Erec_pperp_sub_1p2pi_1p1pi");
  h2_Erec_pperp_sub_1p2pi_1p1pi->Add(h2_Erec_pperp_1p2pi_1p1pi);

  TH1F *h_Erec_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_41prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_fracfeed");
  h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Add(h1_E_rec_1p2pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_41prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Add(h1_E_tot_1p2pi_fracfeed);

//------------------------------------undetected 1p 2pi-> 1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi");
  h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Add(h1_E_rec_1p2pi_1p0pi,-1);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi");
  h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Add(h1_E_tot_1p2pi_1p0pi,-1);

  TH2F *h2_Erec_pperp_sub_1p2pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p1pi->Clone("h2_Erec_pperp_sub_1p2pi_1p0pi");
  h2_Erec_pperp_sub_1p2pi_1p0pi->Add(h2_Erec_pperp_1p2pi_1p0pi,-1);

  TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
  h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(h1_E_rec_1p2pi_1p0pi_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(h1_E_tot_1p2pi_1p0pi_fracfeed,-1);

//------------------------------------undetected 1p 3pi-> 1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_1p3pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Erec_subtruct_piplpimi_1p3pi");
  h_Erec_subtruct_piplpimi_1p3pi->Add(h1_E_rec_1p3pi);

  TH1F *h_Etot_subtruct_piplpimi_1p3pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Etot_subtruct_piplpimi_1p3pi");
  h_Etot_subtruct_piplpimi_1p3pi->Add(h1_E_tot_1p3pi);

  TH2F *h2_Erec_pperp_sub_1p3pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p0pi->Clone("h2_Erec_pperp_sub_1p3pi");
  h2_Erec_pperp_sub_1p3pi->Add(h2_Erec_pperp_1p3pi);

  TH1F *h_Erec_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p3pi_fracfeed");
  h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Add(h1_E_rec_1p3pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p3pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Add(h1_E_tot_1p3pi_fracfeed);


//------------------------------------undetected 2p 2pi ->1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p2pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p3pi->Clone("h_Erec_subtruct_piplpimi_2p2pi");
  h_Erec_subtruct_piplpimi_2p2pi->Add(h1_E_rec_2p2pi);

  TH1F *h_Etot_subtruct_piplpimi_2p2pi=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi->Clone("h_Etot_subtruct_piplpimi_2p2pi");
  h_Etot_subtruct_piplpimi_2p2pi->Add(h1_E_tot_2p2pi);

  TH2F *h2_Erec_pperp_sub_2p2pi=(TH2F*) h2_Erec_pperp_sub_1p3pi->Clone("h2_Erec_pperp_sub_2p2pi");
  h2_Erec_pperp_sub_2p2pi->Add(h2_Erec_pperp_2p2pi);

  TH1F *h_Erec_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p2pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Add(h1_E_rec_2p2pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p2pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Add(h1_E_tot_2p2pi_fracfeed);


//------------------------------------undetected 3p 1pi ->1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_3p1pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p2pi->Clone("h_Erec_subtruct_piplpimi_3p1pi");
  h_Erec_subtruct_piplpimi_3p1pi->Add(h1_E_rec_3p1pi);

  TH1F *h_Etot_subtruct_piplpimi_3p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi->Clone("h_Etot_subtruct_piplpimi_3p1pi");
  h_Etot_subtruct_piplpimi_3p1pi->Add(h1_E_tot_3p1pi);

  TH2F *h2_Erec_pperp_sub_3p1pi=(TH2F*) h2_Erec_pperp_sub_2p2pi->Clone("h2_Erec_pperp_sub_3p1pi");
  h2_Erec_pperp_sub_3p1pi->Add(h2_Erec_pperp_3p1pi);

  TH1F *h_Erec_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_3p1pi_fracfeed");
  h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Add(h1_E_rec_3p1pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_3p1pi_fracfeed");
  h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Add(h1_E_tot_3p1pi_fracfeed);

//------------------------------------undetected 2p 1pi ->2p 0pi  --------------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_3p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi");
  h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Add(h1_E_rec_2p1pi_2p0pi);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi");
  h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Add(h1_E_tot_2p1pi_2p0pi);

  TH2F *h2_Erec_pperp_sub_2p1pi_2p0pi=(TH2F*) h2_Erec_pperp_sub_3p1pi->Clone("h2_Erec_pperp_sub_2p1pi_2p0pi");
  h2_Erec_pperp_sub_2p1pi_2p0pi->Add(h2_Erec_pperp_2p1pi_2p0pi);

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(h1_E_rec_2p1pi_2p0pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(h1_E_tot_2p1pi_2p0pi_fracfeed);

//------------------------------------undetected 2p 1pi ->1p 1pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi");
  h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Add(h1_E_rec_2p1pi_1p1pi);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi");
  h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Add(h1_E_tot_2p1pi_1p1pi);

  TH2F *h2_Erec_pperp_sub_2p1pi_1p1pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_2p0pi->Clone("h2_Erec_pperp_sub_2p1pi_1p1pi");
  h2_Erec_pperp_sub_2p1pi_1p1pi->Add(h2_Erec_pperp_2p1pi_1p1pi);

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_rec_2p1pi_1p1pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_tot_2p1pi_1p1pi_fracfeed);

//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi");
  h_Erec_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_rec_2p1pi_1p0pi,-1);

// apapadop
//  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi");
  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("epRecoEnergy_slice_0");
  h_Etot_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_tot_2p1pi_1p0pi,-1);

  TH2F *h2_Erec_pperp_sub_2p1pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_1p1pi->Clone("h2_Erec_pperp_sub_2p1pi_1p0pi");
  h2_Erec_pperp_sub_2p1pi_1p0pi->Add(h2_Erec_pperp_2p1pi_1p0pi,-1);

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(h1_E_rec_2p1pi_1p0pi_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(h1_E_tot_2p1pi_1p0pi_fracfeed,-1);


 //-----------------------------------looking only at e-, 1pi, undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot = (TH1F*)  h1_E_rec_0pi->Clone("h_Erec_subtruct_piplpimi_noprot");
  h_Erec_subtruct_piplpimi_noprot->Add(h1_E_rec_1pi_weight,-1);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed = (TH1F*)  h1_E_rec_0pi_frac_feed->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed");
  h_Erec_subtruct_piplpimi_noprot_frac_feed->Add(h1_E_rec_1pi_weight_frac_feed,-1);
 //-----------------------------------looking only at e-, 2pi undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot_2pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot->Clone("h_Erec_subtruct_piplpimi_noprot_2pi");
  h_Erec_subtruct_piplpimi_noprot_2pi->Add(h1_E_rec_2pi_weight);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed2pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_frac_feed->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed2pi");
  h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Add(h1_E_rec_2pi_weight_frac_feed);

 //-----------------------------------looking only at e-, 3pi, undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot_3pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_2pi->Clone("h_Erec_subtruct_piplpimi_noprot_3pi");
  h_Erec_subtruct_piplpimi_noprot_3pi->Add(h1_E_rec_3pi_weight);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed3pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed3pi");
  h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Add(h1_E_rec_3pi_weight_frac_feed);

 //-----------------------------------looking only at e-, 4pi, undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot_4pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_3pi->Clone("h_Erec_subtruct_piplpimi_noprot_4pi");
  h_Erec_subtruct_piplpimi_noprot_4pi->Add(h1_E_rec_4pi_weight);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed4pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed4pi");
  h_Erec_subtruct_piplpimi_noprot_frac_feed4pi->Add(h1_E_rec_4pi_weight_frac_feed);

  gDirectory->Write("hist_Files", TObject::kOverwrite);
  // skim_tree->AutoSave();

}

//End Loop function


double genie_analysis::acceptance_c(double p, double cost, double phi, int particle_id,TFile* file_acceptance) {

	//Redefinition of the phi angle
	// because the acceptance maps are defined between (-30,330)

// Check that phi is between (0,360)

	//int redef = -30;
	int redef = 0;

	TH3D * acc;
	TH3D * gen;

	acc = (TH3D*)file_acceptance->Get("Accepted Particles");
	gen = (TH3D*)file_acceptance->Get("Generated Particles");

  //map 330 till 360 to [-30:0] for the acceptance map histogram
  if(phi > (2*TMath::Pi() - TMath::Pi()/6.) ) { phi -= 2*TMath::Pi(); }
	//Find number of generated events

	double pbin_gen = gen->GetXaxis()->FindBin(p);
	double tbin_gen = gen->GetYaxis()->FindBin(cost);
	double phibin_gen = gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	double num_gen = gen->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

	//Find number of accepted events

	double pbin_acc = acc->GetXaxis()->FindBin(p);
	double tbin_acc = acc->GetYaxis()->FindBin(cost);
	double phibin_acc = acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	double num_acc = acc->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

	double acc_ratio = (double)num_acc / (double)num_gen;
	double acc_err = (double)sqrt(acc_ratio*(1-acc_ratio)) / (double)num_gen;


	return acc_ratio;

}


void genie_analysis::LoopCLAS()
{

  std::map<std::string,double>bind_en;
  std::map<std::string,double>target_mass;
  std::map<std::string,double>residual_target_mass;
  std::map<std::string, double> Ecal_offset; //that might not be necessary for simulation data


  target_name = ftarget;   //std string for target name
  en_beam["1161"]=1.161;
  en_beam["2261"]=2.261;
  en_beam["4461"]=4.461;

  en_beam_Ecal["1161"]=1.161;
  en_beam_Ecal["2261"]=2.261;
  en_beam_Ecal["4461"]=4.461;

  en_beam_Eqe["1161"]=1.161;
  en_beam_Eqe["2261"]=2.261;
  en_beam_Eqe["4461"]=4.461;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //nentries =8000000;

  double N_prot1 = 0, N_prot2 = 0,N_prot_both = 0;
  const int n_slice=3;
  const double pperp_min[n_slice]={0.,0.2,0.4};
  const double pperp_max[n_slice]={0.2,0.4,10.};
  TVector3 V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;
  TVector3 V3_phot_angles;

  TString E_acc_file;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV  Configuration parameters and cuts
  {
      E_acc_file="1_161";
  }


  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV  Configuration parameters and cuts
  {
      E_acc_file="2_261";
  }

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5)  //4.4 GeV  Configuration parameters and cuts
  {
      E_acc_file="4_461";
  }
  //Further constants for binding energies and target masses
  Ecal_offset["3He"]=0.004;
  Ecal_offset["4He"]=0.005;
  Ecal_offset["C12"]=0.005;
  Ecal_offset["56Fe"]=0.011;

  bind_en["3He"] = He3_bind_en-D2_bind_en + Ecal_offset["3He"]; //the offset is used to shift the peak to be at 0
  bind_en["4He"] = He4_bind_en-H3_bind_en + Ecal_offset["4He"];
  bind_en["C12"] = C12_bind_en-B_bind_en  + Ecal_offset["C12"];
  bind_en["56Fe"]= Fe_bind_en-Mn_bind_en  + Ecal_offset["56Fe"];
  bind_en["CH2"] = C12_bind_en-B_bind_en;

  target_mass["3He"] = 2*m_prot+m_neut-He3_bind_en;
  target_mass["4He"] = 2*m_prot+2*m_neut-He4_bind_en;
  target_mass["C12"] = 6*m_prot+6*m_neut-C12_bind_en;
  target_mass["56Fe"]= 26*m_prot+30*m_neut-Fe_bind_en;
  target_mass["CH2"] = 6*m_prot+6*m_neut-C12_bind_en;

  residual_target_mass["3He"] = m_prot+m_neut-D2_bind_en;
  residual_target_mass["4He"] = m_prot+2*m_neut-H3_bind_en;
  residual_target_mass["C12"] = 5*m_prot+6*m_neut-B_bind_en;
  residual_target_mass["56Fe"]= 25*m_prot+30*m_neut-Mn_bind_en;
  residual_target_mass["CH2"] = 25*m_prot+30*m_neut-Mn_bind_en;

  //Definition of Histograms
  TH1F *h1_Etot_p_bkgd_slice[n_slice], *h1_Erec_p_bkgd_slice[n_slice];
  TH1F *h1_Etot_Npi0[n_slice],*h1_Erec_Npi0[n_slice];
  TH1F *h1_Etot_Npi1[n_slice],*h1_Erec_Npi1[n_slice];

  TH1F *h1_Erec_bkgd_pipl_pimi_new_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact[n_slice];
  TH1F *h1_Etot_bkgd_pipl_pimi_fact_pipl[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pimi[n_slice];

  TH1F *h1_Etot_bkgd_1p2pi[n_slice],      *h1_Erec_bkgd_1p2pi[n_slice];
  TH1F *h1_Etot_bkgd_1p2pi_1p0pi[n_slice],*h1_Erec_bkgd_1p2pi_1p0pi[n_slice];
  TH1F *h1_Etot_bkgd_1p3pi[n_slice],      *h1_Erec_bkgd_1p3pi[n_slice];

  TH1F *h1_Etot_p_bkgd_slice_2p2pi[n_slice],        *h1_Erec_p_bkgd_slice_2p2pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[n_slice];

  TH1F *h1_Etot_piplpimi_subtruct_fact[n_slice],*h1_Erec_piplpimi_subtruct_fact[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub[n_slice],*h1_Erec_p_bkgd_slice_sub[n_slice];
  TH1F *h1_Etot_3pto1p_slice[n_slice],*h1_Erec_3pto1p_slice[n_slice];
  TH1F *h1_Etot_3pto2p_slice[n_slice],*h1_Erec_3pto2p_slice[n_slice];
  TH1F *h1_Etot_3p1pi_slice[n_slice], *h1_Erec_3p1pi_slice[n_slice];
  TH1F *h1_Etot_4pto3p_slice[n_slice],*h1_Erec_4pto3p_slice[n_slice];
  TH1F *h1_Etot_4pto1p_slice[n_slice],*h1_Erec_4pto1p_slice[n_slice];
  TH1F *h1_Etot_4pto2p_slice[n_slice],*h1_Erec_4pto2p_slice[n_slice];
  TH1F *h1_Etot_43pto1p_slice[n_slice],*h1_Erec_43pto1p_slice[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub1p2pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p2pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_1p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_2p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_2p[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub1p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p2pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub1p3pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p3pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub2p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p2pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub3p1pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub3p1pi_0pi[n_slice];

  TH1F *h1_Etot_p_bkgd_slice_sub32[n_slice],*h1_Erec_p_bkgd_slice_sub32[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub31[n_slice],*h1_Erec_p_bkgd_slice_sub31[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub43[n_slice],*h1_Erec_p_bkgd_slice_sub43[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub41[n_slice],*h1_Erec_p_bkgd_slice_sub41[n_slice];
  TH1F *h1_Erec_p_bkgd_slice_sub42[n_slice],*h1_Etot_p_bkgd_slice_sub42[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub431[n_slice],*h1_Erec_p_bkgd_slice_sub431[n_slice];
  TH2F *h2_N_pi_phot[20];


  gRandom = new TRandom3();
  gRandom->SetSeed(10);

  TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
  TLorentzVector V4_target(0,0,0,target_mass[ftarget]);

  //Output file definition
  TFile *file_out = new TFile(Form("e2a_ep_%s_%s_neutrino6_united4_radphot_test_data.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");

  //initialize Fiducial functions for EC limits
  fiducialcut->InitEClimits();

  //Definition and initialization of Histograms
  TH1F *h1_el_Mott_crosssec = new TH1F("h1_el_Mott_crosssec","",200,0.,0.01);
  TH1F *h1_Wvar = new TH1F("h1_Wvar","",400,0,3);
  TH1F *h1_Wvar_weight = new TH1F("h1_Wvar_weight","",400,0,3);
  TH1F *h1_xbjk = new TH1F("h1_xbjk","",400,0,3);
  TH1F *h1_xbjk_weight = new TH1F("h1_xbjk_weight","",400,0,3);
  TH1F *h1_Q2 = new TH1F("h1_Q2","",400,0,6);
  TH1F *h1_Q2_weight = new TH1F("h1_Q2_weight","",400,0,6);
  TH1F *h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
  TH1F *h1_Nprot=new TH1F("h1_Nprot","",10,0,5);
  TH1F *h1_Nphot=new TH1F("h1_Nphot","",10,0,5);
  TH1F *h1_Npiphot=new TH1F("h1_Npiphot","",10,0,5);
  TH1F *h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,0,5);
  TH1F *h1_Npi=new TH1F("h1_Npi","",10,0,5);
  TH1F *h1_Npipl=new TH1F("h1_Npipl","",10,0,5);
  TH1F *h1_Npimi=new TH1F("h1_Npimi","",10,0,5);
  TH1F *h1_el_mom = new TH1F("h1_el_mom","",100,1.2,6);
  TH1F *h1_el_mom_corr = new TH1F("h1_el_mom_corr","",100,1.2,6);
  TH1F *h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
  TH1F *h1_prot_mom = new TH1F("h1_prot_mom","",300,0,3);
  TH1F *h1_prot_mom_ratio = new TH1F("h1_prot_mom_ratio","",50,0.97,1.2);

  //Binning for fractional feed down energy histograms
  int N_qe;
  double *x_qe;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
    N_qe=109;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=64;i++)  x_qe[i]=-1+i*0.015;
    for (int i=0;i<=44;i++)  x_qe[i+65]=-0.04+(i+1)*0.01;
  }

  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
    N_qe=109;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=64;i++)  x_qe[i]=-1+i*0.015;
    for (int i=0;i<=44;i++)  x_qe[i+65]=-0.04+(i+1)*0.01;
  }

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
    N_qe=82;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=29;i++)  x_qe[i]=-1+i*0.03;
    for (int i=0;i<=52;i++)  x_qe[i+30]=-0.13+(i+1)*0.01;
  }

 //Definitions of further Histograms
  TH1F *h1_E_rec_1pi_weight_frac_feed=new TH1F("h1_E_rec_1pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_2pi_weight_frac_feed=new TH1F("h1_E_rec_2pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_3pi_weight_frac_feed=new TH1F("h1_E_rec_3pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pi_weight_frac_feed=new TH1F("h1_E_rec_4pi_weight_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_rec_0pi_frac_feed=new TH1F("h1_E_rec_0pi_frac_feed","",N_qe,x_qe);
  TH1F *h1_E_tot_cut2_fracfeed = new TH1F("h1_E_tot_cut2_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_cut2_new_fracfeed = new TH1F("h1_E_rec_cut2_new_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_p_bkgd_fracfeed = new TH1F("h1_E_tot_p_bkgd_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_p_bkgd_fracfeed = new TH1F("h1_E_rec_p_bkgd_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_3pto2p_fracfeed = new TH1F("h1_E_tot_3pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_3pto2p_fracfeed = new TH1F("h1_E_rec_3pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_3pto1p_fracfeed = new TH1F("h1_E_tot_3pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_3pto1p_fracfeed = new TH1F("h1_E_rec_3pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_4pto3p_fracfeed = new TH1F("h1_E_tot_4pto3p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pto3p_fracfeed = new TH1F("h1_E_rec_4pto3p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_43pto1p_fracfeed = new TH1F("h1_E_tot_43pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_43pto1p_fracfeed = new TH1F("h1_E_rec_43pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_4pto2p_fracfeed = new TH1F("h1_E_tot_4pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pto2p_fracfeed = new TH1F("h1_E_rec_4pto2p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_4pto1p_fracfeed = new TH1F("h1_E_tot_4pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_4pto1p_fracfeed = new TH1F("h1_E_rec_4pto1p_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_1p2pi_fracfeed = new TH1F("h1_E_tot_1p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_1p2pi_fracfeed = new TH1F("h1_E_rec_1p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_1p3pi_fracfeed = new TH1F("h1_E_tot_1p3pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_1p3pi_fracfeed = new TH1F("h1_E_rec_1p3pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_2p2pi_fracfeed = new TH1F("h1_E_tot_2p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_2p2pi_fracfeed = new TH1F("h1_E_rec_2p2pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_3p1pi_fracfeed = new TH1F("h1_E_tot_3p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_3p1pi_fracfeed = new TH1F("h1_E_rec_3p1pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_tot_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_rec_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_rec_undetfactor_fracfeed = new TH1F("h1_E_rec_undetfactor_fracfeed","",N_qe,x_qe);
  TH1F *h1_E_tot_undetfactor_fracfeed = new TH1F("h1_E_tot_undetfactor_fracfeed","",N_qe,x_qe);

  TH1F *h1_theta0=new TH1F("h1_theta0","",300,0,180);
  TH2F *h2_Ecal_Eqe=new TH2F("h2_Ecal_Eqe","",600,0,5,600,0,5);
  TH2F *h2_EqeEcalratio_Eqe=new TH2F("h2_EqeEcalratio_Eqe","",600,0,5,300,0,2);
  TH2F *h2_EqeEcaldiff_Eqe=new TH2F("h2_EqeEcaldiff_Eqe","",600,0,5,300,-3,3);
  TH2F *h2_N_prot_pi=new TH2F("h2_N_prot_pi","",10,0,5,10,0,5);
  TH2F *h2_N_prot_pi_phot=new TH2F("h2_N_prot_pi_phot","",10,0,5,10,0,5);
  TH2F *h2_N_prot_pi_phot_nonrad=new TH2F("h2_N_prot_pi_phot_nonrad","",10,0,5,10,0,5);
  TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,10,60);
  TH2F *h2_el_mom_diff = new TH2F("h2_el_mom_diff","",500,0.,1.,500,-0.1,0.1);
  TH2F *h2_Q2_nu = new TH2F("h2_Q2_nu","",200,0,3.5,200,0,5);
  TH2F *h2_Q2_nu_weight = new TH2F("h2_Q2_nu_weight","",200,0,3.5,200,0,5);
  TH2F *h2_Q2_xbjk_weight = new TH2F("h2_Q2_xbjk_weight","",200,0,3,200,0,5);
  TH2F *h2_Q2_W=new TH2F("h2_Q2_W","",200,0,3,200,0,5);
  TH2F *h2_xB_W=new TH2F("h2_xB_W","",200,0,3,200,0,3);
  TH2F *h2_Q2_W_weight=new TH2F("h2_Q2_W_weight","",200,0,3,200,0,5);
  TH2F *h2_el_pcorr_puncorr = new TH2F("h2_el_pcorr_puncorr","",100,0,1,100,0,3);
  TH2F *h2_Erec_pperp = new TH2F("h2_Erec_pperp","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_newcut2 = new TH2F("h2_Erec_pperp_newcut2","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_cut3 = new TH2F("h2_Erec_pperp_cut3","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p = new TH2F("h2_Erec_pperp_2p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_321p = new TH2F("h2_Erec_pperp_321p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_31p = new TH2F("h2_Erec_pperp_31p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_4321p = new TH2F("h2_Erec_pperp_4321p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_431p = new TH2F("h2_Erec_pperp_431p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_421p = new TH2F("h2_Erec_pperp_421p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_41p = new TH2F("h2_Erec_pperp_41p","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p1pi = new TH2F("h2_Erec_pperp_1p1pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p2pi_1p0pi = new TH2F("h2_Erec_pperp_1p2pi_1p0pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p2pi_1p1pi = new TH2F("h2_Erec_pperp_1p2pi_1p1pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p1pi_2p0pi = new TH2F("h2_Erec_pperp_2p1pi_2p0pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p1pi_1p1pi = new TH2F("h2_Erec_pperp_2p1pi_1p1pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p1pi_1p0pi = new TH2F("h2_Erec_pperp_2p1pi_1p0pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_1p3pi = new TH2F("h2_Erec_pperp_1p3pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_2p2pi = new TH2F("h2_Erec_pperp_2p2pi","",400,0,1,400,0,6.);
  TH2F *h2_Erec_pperp_3p1pi = new TH2F("h2_Erec_pperp_3p1pi","",400,0,1,400,0,6.);
  TH2F *h2_pperp_W=new TH2F("h2_pperp_W","",200,0,3,200,0,2);

  TH2F *h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);

  //Binning for energy reconstruction histograms
  int n_bins;
  double *x_values;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
    n_bins=38;
    x_values=new double[n_bins+1];
    for (int i=0;i<=17;i++) x_values[i]=0.4+i*0.04;
    for (int i=0;i<=20;i++) x_values[i+18]=1.08+(i+1)*0.02;
  }

  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
    n_bins=54;
    x_values=new double[n_bins+1];
    for (int i=0;i<=23;i++) x_values[i]=i*0.09;
    for (int i=0;i<=30;i++) x_values[i+24]=2.07+(i+1)*0.03;
  }

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
    n_bins=38;
    x_values=new double[n_bins+1];
    for (int i=0;i<=21;i++)  x_values[i]=i*0.2;
    for (int i=0;i<=16;i++)  x_values[i+22]=4.2+(i+1)*0.05;
  }

 //Definitions of further Histograms
  TH1F *h1_E_rec_2p_det = new TH1F("h1_E_rec_2p_det","",n_bins,x_values);
  TH1F *h1_E_tot_2p_det = new TH1F("h1_E_tot_2p_det","",n_bins,x_values);
  TH1F *h1_E_tot_p_bkgd = new TH1F("h1_E_tot_p_bkgd","",n_bins,x_values);
  TH1F *h1_E_rec_p_bkgd = new TH1F("h1_E_rec_p_bkgd","",n_bins,x_values);
  TH1F *h1_E_tot_3pto1p = new TH1F("h1_E_tot_3pto1p","",n_bins,x_values);
  TH1F *h1_E_rec_3pto1p = new TH1F("h1_E_rec_3pto1p","",n_bins,x_values);
  TH1F *h1_E_tot_43pto1p = new TH1F("h1_E_tot_43pto1p","",n_bins,x_values);
  TH1F *h1_E_rec_43pto1p = new TH1F("h1_E_rec_43pto1p","",n_bins,x_values);
  TH1F *h1_E_tot_3pto2p =new TH1F("h1_E_tot_3pto2p","",n_bins,x_values);
  TH1F *h1_E_rec_3pto2p = new TH1F("h1_E_rec_3pto2p","",n_bins,x_values);
  TH1F *h1_E_tot_4pto1p =new TH1F("h1_E_tot_4pto1p","",n_bins,x_values);
  TH1F *h1_E_rec_4pto1p = new TH1F("h1_E_rec_4pto1p","",n_bins,x_values);
  TH1F *h1_E_tot_4pto3p =new TH1F("h1_E_tot_4pto3p","",n_bins,x_values);
  TH1F *h1_E_rec_4pto3p = new TH1F("h1_E_rec_4pto3p","",n_bins,x_values);
  TH1F *h1_E_tot_4pto2p =new TH1F("h1_E_tot_4pto2p","",n_bins,x_values);
  TH1F *h1_E_rec_4pto2p = new TH1F("h1_E_rec_4pto2p","",n_bins,x_values);
  TH1F *h1_E_rec = new TH1F("h1_E_rec","",n_bins,x_values);
  TH1F *h1_E_rec_0pi = new TH1F("h1_E_rec_0pi","",n_bins,x_values);
  TH1F *h1_E_rec_1pi = new TH1F("h1_E_rec_1pi","",n_bins,x_values);
  TH1F *h1_E_rec_1pi_weight = new TH1F("h1_E_rec_1pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_2pi_weight = new TH1F("h1_E_rec_2pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_3pi_weight = new TH1F("h1_E_rec_3pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_4pi_weight = new TH1F("h1_E_rec_4pi_weight","",n_bins,x_values);
  TH1F *h1_E_rec_20pi = new TH1F("h1_E_rec_20pi","",n_bins,x_values);
  TH1F *h1_E_rec_21pi = new TH1F("h1_E_rec_21pi","",n_bins,x_values);
  TH1F *h1_E_rec_30pi = new TH1F("h1_E_rec_30pi","",n_bins,x_values);
  TH1F *h1_E_rec_310pi = new TH1F("h1_E_rec_310pi","",n_bins,x_values);
  TH1F *h1_E_rec_320pi = new TH1F("h1_E_rec_320pi","",n_bins,x_values);
  TH1F *h1_E_rec_3210pi = new TH1F("h1_E_rec_3210pi","",n_bins,x_values);
  TH1F *h1_E_rec_40pi = new TH1F("h1_E_rec_40pi","",n_bins,x_values);
  TH1F *h1_E_rec_410pi = new TH1F("h1_E_rec_410pi","",n_bins,x_values);
  TH1F *h1_E_rec_420pi = new TH1F("h1_E_rec_420pi","",n_bins,x_values);
  TH1F *h1_E_rec_4210pi = new TH1F("h1_E_rec_4210pi","",n_bins,x_values);
  TH1F *h1_E_rec_430pi = new TH1F("h1_E_rec_430pi","",n_bins,x_values);
  TH1F *h1_E_rec_4310pi = new TH1F("h1_E_rec_4310pi","",n_bins,x_values);
  TH1F *h1_E_rec_4320pi = new TH1F("h1_E_rec_4320pi","",n_bins,x_values);
  TH1F *h1_E_rec_43210pi = new TH1F("h1_E_rec_43210pi","",n_bins,x_values);
  TH1F *h1_E_rec_1prot  = new TH1F("h1_E_rec_1prot","",n_bins,x_values);
  TH1F *h1_E_tot_1prot  = new TH1F("h1_E_tot_1prot","",n_bins,x_values);
  TH1F *h1_E_rec_cutpi1_piplpimi = new TH1F("h1_E_rec_cutpi1_piplpimi","",n_bins,x_values);
  TH1F *h1_E_tot_cutpi1_piplpimi = new TH1F("h1_E_tot_cutpi1_piplpimi","",n_bins,x_values);
  TH1F *h1_E_tot = new TH1F("h1_E_tot","",n_bins,x_values);
  TH1F *h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
  TH1F *h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
  TH1F *h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
	TH1F *h1_E_rec_undetfactor  = new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
  TH1F *h1_E_tot_undetfactor  = new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
  TH1F *h1_E_tot_1p2pi  = new TH1F("h1_E_tot_1p2pi","",n_bins,x_values);
  TH1F *h1_E_rec_1p2pi  = new TH1F("h1_E_rec_1p2pi","",n_bins,x_values);
  TH1F *h1_E_tot_1p3pi  = new TH1F("h1_E_tot_1p3pi","",n_bins,x_values);
  TH1F *h1_E_rec_1p3pi  = new TH1F("h1_E_rec_1p3pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p2pi  = new TH1F("h1_E_tot_2p2pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p2pi  = new TH1F("h1_E_rec_2p2pi","",n_bins,x_values);
  TH1F *h1_E_tot_3p1pi  = new TH1F("h1_E_tot_3p1pi","",n_bins,x_values);
  TH1F *h1_E_rec_3p1pi  = new TH1F("h1_E_rec_3p1pi","",n_bins,x_values);
  TH1F *h1_E_tot_1p2pi_1p0pi  = new TH1F("h1_E_tot_1p2pi_1p0pi","",n_bins,x_values);
  TH1F *h1_E_rec_1p2pi_1p0pi  = new TH1F("h1_E_rec_1p2pi_1p0pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p1pi_2p0pi  = new TH1F("h1_E_tot_2p1pi_2p0pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p1pi_2p0pi  = new TH1F("h1_E_rec_2p1pi_2p0pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p1pi_1p1pi  = new TH1F("h1_E_tot_2p1pi_1p1pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p1pi_1p1pi  = new TH1F("h1_E_rec_2p1pi_1p1pi","",n_bins,x_values);
  TH1F *h1_E_tot_2p1pi_1p0pi  = new TH1F("h1_E_tot_2p1pi_1p0pi","",n_bins,x_values);
  TH1F *h1_E_rec_2p1pi_1p0pi  = new TH1F("h1_E_rec_2p1pi_1p0pi","",n_bins,x_values);

  //Defintions of Histogram for each slice
  for(int h = 0; h < n_slice; h++){
   h1_Erec_p_bkgd_slice[h]= new TH1F(Form("h1_Erec_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice[h]= new TH1F(Form("h1_Etot_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3pto1p_slice[h]= new TH1F(Form("h1_Erec_3pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3pto1p_slice[h]= new TH1F(Form("h1_Etot_3pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3pto2p_slice[h]= new TH1F(Form("h1_Erec_3pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3pto2p_slice[h]= new TH1F(Form("h1_Etot_3pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3p1pi_slice[h]= new TH1F(Form("h1_Erec_3p1pi_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3p1pi_slice[h]= new TH1F(Form("h1_Etot_3p1pi_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_43pto1p_slice[h]= new TH1F(Form("h1_Etot_43pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_43pto1p_slice[h]= new TH1F(Form("h1_Erec_43pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto3p_slice[h]= new TH1F(Form("h1_Erec_4pto3p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto3p_slice[h]= new TH1F(Form("h1_Etot_4pto3p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto2p_slice[h]= new TH1F(Form("h1_Erec_4pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto2p_slice[h]= new TH1F(Form("h1_Etot_4pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto1p_slice[h]= new TH1F(Form("h1_Erec_4pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto1p_slice[h]= new TH1F(Form("h1_Etot_4pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_Npi0[h] = new TH1F(Form("h1_Etot_Npi0_%d",h+1),"",n_bins,x_values);
   h1_Erec_Npi0[h] = new TH1F(Form("h1_Erec_Npi0_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_pipl_pimi_fact[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_pipl_pimi_fact_pipl[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pipl_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_pipl_pimi_fact_pimi[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pimi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_pipl_pimi_new_fact[h]= new TH1F(Form("h1_Erec_bkgd_pipl_pimi_new_fact_%d",h+1),"",n_bins,x_values);
   h1_Etot_Npi1[h] = new TH1F(Form("h1_Etot_Npi1_%d",h+1),"",n_bins,x_values);
   h1_Erec_Npi1[h] = new TH1F(Form("h1_Erec_Npi1_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_1p2pi[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_1p2pi[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to2p0pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to2p0pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p2pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p2pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_p_bkgd_slice_2p2pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p2pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_1p2pi_1p0pi[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_1p2pi_1p0pi[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_1p0pi_%d",h+1),"",n_bins,x_values);
   h1_Etot_bkgd_1p3pi[h] = new TH1F(Form("h1_Etot_bkgd_1p3pi_%d",h+1),"",n_bins,x_values);
   h1_Erec_bkgd_1p3pi[h] = new TH1F(Form("h1_Erec_bkgd_1p3pi_%d",h+1),"",n_bins,x_values);
  }


  for(int h=0;h<20;h++){
    h2_N_pi_phot[h]=new TH2F(Form("h2_N_pi_phot_%d",h),"",10,0,5,10,0,5);
  }

  TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();


/** Beginning of Event Loop **/
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<200000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    int nb = GetEntry(jentry);
    if( jentry%200000 == 0 )
    {
	         gDirectory->Write("hist_Files", TObject::kOverwrite);
		       cout<<jentry<<endl;
    }

    fTorusCurrent=2250; //only true for 2.2 and 4.4 GeV
    //Make new if condition depending on beam energy. F.H. 9/9/19
//1500 is not used anymore, no runnumber in GENIE simulation files
//Keep 1500 with extra warning if it used :)
  /*  if((runnb>18283 && runnb<18289) || (runnb>18300 && runnb<18304) || (runnb>18317 && runnb<18329))        fTorusCurrent=750;    //setting appropriate torrus magnet current
    else if ((runnb>18293 && runnb<18301) || (runnb>18305 && runnb<18317) || (runnb>18328 && runnb<18336))  fTorusCurrent=1500;
    else fTorusCurrent=2250;
*/
    if(jentry == 0){ //was n_evt == 1 before but jentry = n_evnt - 1
          //SetMomCorrParameters(); Functions is missing F.H. 08/01/19
          fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
          fiducialcut->SetFiducialCutParameters(fbeam_en);
          std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
          rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
          std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;

    }

    //Resets q vector to (0,0,0)
    rotation->ResetQVector();

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Outgoing e', Onl corrected  values from GENIE data structure
		TLorentzVector V4_el(pxl,pyl,pzl,El);
    TVector3 V3_el(pxl,pyl,pzl);


    double el_momentum = V3_el.Mag();
    double el_theta = V3_el.Theta();
    //Definition as for data. WARNING: Needs to be checked if this also works for GENIE data F.H. 08/24/19
    double el_phi_mod = V3_el.Phi()*TMath::RadToDeg() + 30; //Add extra 30 degree rotation in phi;
    if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree



    //Calculated Mott Cross Section and Weights for Inclusive Histograms
    double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(el_theta)+1))/(2*pow(El,2)*pow((1-cos(el_theta)),2.));
    double WeightIncl = 1./Mott_cross_sec;
    el_theta = el_theta*TMath::RadToDeg();

    double E_rec = (m_prot*bind_en[ftarget]+m_prot*V4_el.E())/(m_prot-V4_el.E()+V4_el.Rho()*cos(el_theta));  //using the same value of single nucleon separation E Ecal and Eqe

    //Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
    double nu = -(V4_el-V4_beam).E();
    double reco_Q2 = -(V4_el-V4_beam).Mag2();
    double x_bjk = reco_Q2/(2*m_prot*nu);
    TVector3 V3_q = (V4_beam-V4_el).Vect();
    double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

    //Set q vector for the following rotations for the subtraction procedure
    rotation->SetQVector(V3_q);
//    rotation->PrintQVector();

    h1_el_mom_corr->Fill(V4_el.Rho(),WeightIncl);

    //Filling Histogram for electron kinematics
    h1_xbjk->Fill(x_bjk);
    h1_xbjk_weight->Fill(x_bjk,WeightIncl);
    h1_Q2->Fill(reco_Q2);
    h1_Q2_weight->Fill(reco_Q2,WeightIncl);
    h2_Q2_nu->Fill(nu,reco_Q2);
    h2_Q2_nu_weight->Fill(nu,reco_Q2,WeightIncl);
    h2_Q2_xbjk_weight->Fill(x_bjk,reco_Q2,WeightIncl);
    h1_Wvar->Fill(W_var);
    h1_Wvar_weight->Fill(W_var,WeightIncl);
    h2_Q2_W->Fill(W_var,reco_Q2);
    h2_xB_W->Fill(W_var,x_bjk);
    h2_Q2_W_weight->Fill(W_var,reco_Q2,WeightIncl);

    h1_el_Mott_crosssec->Fill(Mott_cross_sec);
    h2_el_theta_phi->Fill(el_phi_mod,el_theta);
    h1_el_theta->Fill(el_theta);

    //Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

    //Index variables for hadrons (p and pions)
    int index_p[20]; //index for each proton
    int index_pi[20]; //index for each pion
    int ind_pi_phot[20]; //index for pions and photons
    int index_pipl[20]; //index for each pi plus
    int index_pimi[20]; //index for each pi minus
    bool ec_radstat_n[20];
    //Setting arrays
    for (int i = 0; i<20; i++) {
      index_p[i] = -1;   index_pi[i] = -1;   index_pipl[i] = -1;   index_pimi[i] = -1;   ind_pi_phot[i] = -1; ec_radstat_n[i] = false;
    }
    //Number of hadrons
    int num_p = 0;
    int num_pi = 0;
    int num_pi_phot = 0; //couting all pions and photons
    int num_pimi = 0;
    int num_pipl = 0;
    int num_pi_phot_nonrad=0; //counting all pions and non-radiation photons
    int num_phot_rad = 0; //counting radiation photons
    //Index and number variables for neutral particles
    int ec_num_n = 0;

    const double phot_rad_cut = 40;
    const double phot_e_phidiffcut=30; //electron - photon phi difference cut

    // Creating vectors to store id of particles in the array
    vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
    ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

   //Loop for Hadrons
    for (int i = 0; i < nf; i++)
    {
        //Start of proton selection
       if (pdgf[i] == 2212) { // && pf[i] > 0.3) {

          num_p = num_p + 1;
          index_p[num_p - 1] = i;
          ProtonID.push_back(i);

       }

       if (pdgf[i] == -211) { // && pf[i] > 0.15)  {

          num_pimi = num_pimi + 1;
          num_pi = num_pi + 1;
          num_pi_phot = num_pi_phot + 1;
          num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
          index_pimi[num_pi_phot - 1] = i;
          index_pi[num_pi_phot - 1] = i;
          ind_pi_phot[num_pi_phot - 1] = i;
          PiMinusID.push_back(i);

       }

       if ( pdgf[i] == 211) { // && pf[i] > 0.15)  {

          num_pipl = num_pipl + 1;
          num_pi  = num_pi + 1;
          num_pi_phot = num_pi_phot + 1;
          num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
          index_pipl[num_pi_phot - 1] = i;
          index_pi[num_pi_phot - 1] = i;
          ind_pi_phot[num_pi_phot - 1] = i;
          PiPlusID.push_back(i);

       }


       if (pdgf[i] == 22) { // && pf[i] > 0.3) {

              ec_num_n = ec_num_n + 1;
              num_pi_phot = num_pi_phot + 1;
              ind_pi_phot[num_pi_phot - 1] = i;  //no events with photons for now F.H. 28.08.19
       				PhotonID.push_back(i);
              // WARNING: THe following needs to be implemented for simulation data F.H. 24.08.19
              //Cut on Radiation photon via angle with respect to the electron
              V3_phot_angles.SetXYZ(pxf[i],pyf[i],pzf[i]);
              double neut_phi_mod = V3_phot_angles.Phi()*TMath::RadToDeg() + 30; //Add 30 degree
              if (neut_phi_mod < 0) neut_phi_mod = neut_phi_mod + 360;  //Neutral particle is between 0 and 360 degree

              //within 40 degrees in theta and 30 degrees in phi
              if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg() < phot_rad_cut && fabs(neut_phi_mod-el_phi_mod) < phot_e_phidiffcut ) {
                  ec_radstat_n[num_pi_phot - 1] = true; //select radiation photons
                  num_phot_rad = num_phot_rad + 1;
              }
              if(!ec_radstat_n[num_pi_phot]) num_pi_phot_nonrad = num_pi_phot_nonrad + 1;

       }

    } //end of hadron loop

    //Skip event if there is at least one radiation photon
    if (num_phot_rad > 0) {
      continue;
    }
//std::cout << "num_p = " << num_p << "  num_pi = " << num_pi << "  num_pipl = " << num_pipl << "  num_pimi = " << num_pimi << std::endl;

    //Filling Histograms with multiplicities
    h1_Npi->Fill(num_pi);
    h1_Nprot->Fill(num_p);
    h1_Nphot->Fill(ec_num_n);
    h1_Npipl->Fill(num_pipl);
    h1_Npimi->Fill(num_pimi);
    h1_Npiphot->Fill(num_pi_phot);
    h1_Npiphot_norad->Fill(num_pi_phot_nonrad);
    h2_N_prot_pi->Fill(num_pi,num_p);
    h2_N_prot_pi_phot->Fill(num_pi+ec_num_n,num_p);
    h2_N_prot_pi_phot_nonrad->Fill(num_pi_phot_nonrad,num_p);
    h2_N_pi_phot[num_p]->Fill(ec_num_n,num_pi);

    //Events with exactly 2 protons
	  if(num_p == 2)
    {

        //LorentzVectors for protons without momentum corrections
          TLorentzVector V4_prot_uncorr1(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
          TLorentzVector V4_prot_uncorr2(pxf[index_p[1]],pyf[index_p[1]],pzf[index_p[1]],TMath::Sqrt(m_prot*m_prot+pf[index_p[1]]*pf[index_p[1]]));


          TVector3 V3_prot_corr1(pxf[index_p[0]+60],pyf[index_p[0]+60],pzf[index_p[0]+60]);
          TVector3 V3_prot_corr2(pxf[index_p[1]+60],pyf[index_p[1]+60],pzf[index_p[1]+60]);


          TVector3 V3_2prot_uncorr[2];
          V3_2prot_uncorr[0] = V4_prot_uncorr1.Vect();
          V3_2prot_uncorr[1] = V4_prot_uncorr2.Vect();

          TVector3 V3_2prot_corr[2];
          V3_2prot_corr[0] = V3_prot_corr1;
          V3_2prot_corr[1] = V3_prot_corr2;

 //---------------------------------- 2p 0pi->  1p0pi   ----------------------------------------------


          double E_tot_2p[2]={0};
          double p_perp_tot_2p[2]={0};
          N_prot_both = 0;
          double P_N_2p[2]={0};

          rotation->prot2_rot_func(V3_2prot_corr, V3_2prot_uncorr, V4_el, E_tot_2p, p_perp_tot_2p, P_N_2p , &N_prot_both);

          if(num_pi_phot==0 && N_prot_both!=0){ // no photons for now F.H. 29.8.19
	  //          if(num_pi == 0 && N_prot_both != 0){
            double histoweight = 1./Mott_cross_sec; //total weight from 2p acceptance , 1e acceptance, Mott, and GENIE weight

            for(int f = 0; f < num_p; f++){    //looping through two protons
              h1_E_tot_p_bkgd->Fill(E_tot_2p[f],P_N_2p[f]*histoweight);
              h1_E_rec_p_bkgd->Fill(E_rec,P_N_2p[f]*histoweight);
              h2_Erec_pperp_2p->Fill(p_perp_tot_2p[f],E_rec,P_N_2p[f]*histoweight);
              h1_E_tot_p_bkgd_fracfeed->Fill((E_tot_2p[f]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_N_2p[f]*histoweight);
              h1_E_rec_p_bkgd_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_N_2p[f]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[f],-P_N_2p[f]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[f]) *TMath::RadToDeg(),-P_N_2p[f]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[f],-P_N_2p[f]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[f],-P_N_2p[f]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[f],-P_N_2p[f]*histoweight);

              for(int i = 0 ; i < n_slice; i++)
              {
	               if (p_perp_tot_2p[f] < pperp_max[i] && p_perp_tot_2p[f] > pperp_min[i]) {
	                  h1_Etot_p_bkgd_slice[i]->Fill(E_tot_2p[f],P_N_2p[f]*histoweight);
	                  h1_Erec_p_bkgd_slice[i]->Fill(E_rec,P_N_2p[f]*histoweight);
	               }
	            }
            }//looping through two protons

            h1_E_tot_2p_det->Fill(E_tot_2p[0],histoweight);
            h1_E_rec_2p_det->Fill(E_rec,histoweight);
          }//no pions cut and N_prot_both!=0


 //---------------------------------- 2p 1pi   ----------------------------------------------
 //Const int can be placed somewhere up after if for 2 protons F.H. 05.09.19
          const int N_2prot=2;
//Variable might can be placed in a more local context F.H. 05.09.19
          double Ecal_2p1pi_to2p0pi[N_2prot]={0};
          double p_miss_perp_2p1pi_to2p0pi[N_2prot]={0};

          if (num_pi_phot==1) {
//          if (num_pi == 1) {

            int charge = 0;
            if (num_pimi == 1 && num_pipl == 0) {
              charge = -1;
            }
            else if (num_pipl == 1 && num_pimi == 0) {
              charge = 1;
            }
            //radiation status is false if real photon
            else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
              charge = 0;
            }
            //skip radiation photon
            else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
	            charge = 0;
            }
            else {
              std::cout << "2 proton 1 pi loop Problem with Pion count and identification. Num Pi = " << num_pi << " , Num PiPlus = " << num_pipl << " , Num PiMinus = " << num_pimi << std::endl;
              continue;
            }

            TVector3 V3_1pi_corr;
            V3_1pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);

            double P_2p1pito2p0pi[2] = {0};
            double P_2p1pito1p1pi[2] = {0};
            double P_2p1pito1p0pi[2] = {0};
            double Ptot = 0;

            rotation->prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_1pi_corr, charge, ec_radstat_n[0],V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

            double histoweight = 1./Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19

            for(int z=0; z < N_2prot; z++){ //looping over two protons

 //---------------------------------- 2p 1pi ->2p 0pi   ----------------------------------------------

              h1_E_tot_2p1pi_2p0pi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h1_E_rec_2p1pi_2p0pi->Fill(E_rec,P_2p1pito2p0pi[z]*histoweight);
              h2_Erec_pperp_2p1pi_2p0pi->Fill(p_miss_perp_2p1pi_to2p0pi[z],E_rec,P_2p1pito2p0pi[z]*histoweight);
              h1_E_tot_2p1pi_2p0pi_fracfeed->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
              h1_E_rec_2p1pi_2p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z]) *TMath::RadToDeg(),P_2p1pito2p0pi[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);

              for(int i = 0; i < n_slice; i++){
                if (p_miss_perp_2p1pi_to2p0pi[z]<pperp_max[i] && p_miss_perp_2p1pi_to2p0pi[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
	                 h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(E_rec,P_2p1pito2p0pi[z]*histoweight);
                }
              }
 //---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

              h1_E_tot_2p1pi_1p1pi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h1_E_rec_2p1pi_1p1pi->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
              h2_Erec_pperp_2p1pi_1p1pi->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p1pi[z]*histoweight);
              h1_E_tot_2p1pi_1p1pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
              h1_E_rec_2p1pi_1p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_2p1pito1p1pi[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

              for(int i = 0; i < n_slice; i++){
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
                  h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
                  h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
                }
              }
//---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------

              h1_E_tot_2p1pi_1p0pi->Fill(E_tot_2p[z], P_2p1pito1p0pi[z]*histoweight);
              h1_E_rec_2p1pi_1p0pi->Fill(E_rec,P_2p1pito1p0pi[z]*histoweight);
              h2_Erec_pperp_2p1pi_1p0pi->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p0pi[z]*histoweight);
              h1_E_tot_2p1pi_1p0pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p0pi[z]*histoweight);
              h1_E_rec_2p1pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p0pi[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),-P_2p1pito1p0pi[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

              for(int i = 0; i < n_slice; i++) {
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*histoweight);
	                 h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_rec,P_2p1pito1p0pi[z]*histoweight);
                }
              }

            }//filling the histograms for 2protons
          }//1pi requirement


//---------------------------------- 2p 2pi   ----------------------------------------------
          const int N_2pi=2;
          int q_pi2[2];
          double Ecal_2p2pi[N_2prot];
          double p_miss_perp_2p2pi[N_2prot];
          double Ptot_2p[2]={0};
          bool ecstat_pi2[N_2pi]={false};

          if (num_pi_phot==2) {
//          if (num_pi == 2) {

            TVector3 V3_2pi_corr[N_2pi];

            for (int i = 0; i < num_pi_phot; i++) {

              ecstat_pi2[i] = ec_radstat_n[i];
              if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
                q_pi2[i] = 1;
              }
              else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
                q_pi2[i] = -1;
              }
              else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
                q_pi2[i] = 0;
              }
              else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
		            q_pi2[i] = 0;
              }
              else {  std::cout << "WARNING: 2p 2pion event: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

              V3_2pi_corr[i].SetXYZ(pxf[ind_pi_phot[i]],pyf[ind_pi_phot[i]],pzf[ind_pi_phot[i]]);


            }

            rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,q_pi2,ecstat_pi2 ,V4_el, Ecal_2p2pi,p_miss_perp_2p2pi,Ptot_2p);

            double histoweight = 1./Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19


            for(int z = 0; z < N_2prot; z++){ //looping over two protons

//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

              h1_E_tot_2p2pi->Fill(E_tot_2p[z], Ptot_2p[z]*histoweight);
              h1_E_rec_2p2pi->Fill(E_rec,Ptot_2p[z]*histoweight);
              h2_Erec_pperp_2p2pi->Fill(p_perp_tot_2p[z],E_rec,Ptot_2p[z]*histoweight);
              h1_E_tot_2p2pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],Ptot_2p[z]*histoweight);
              h1_E_rec_2p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],Ptot_2p[z]*histoweight);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],Ptot_2p[z]*histoweight);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),Ptot_2p[z]*histoweight);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],Ptot_2p[z]*histoweight);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],Ptot_2p[z]*histoweight);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],Ptot_2p[z]*histoweight);

              for(int i = 0; i < n_slice; i++)
              {
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p2pi[i]->Fill(E_tot_2p[z],Ptot_2p[z]*histoweight);
	                 h1_Erec_p_bkgd_slice_2p2pi[i]->Fill(E_rec,Ptot_2p[z]*histoweight);
                }
              }
            }   //Filling the histogram for two protons
          }//2pi requirement

    } //2prot requirement

  //Events with exactly 3 protons
	  if(num_p == 3){


	      const int N_3p = 3;
	      TLorentzVector V4_p_uncorr[N_3p], V4_p_corr[N_3p],V4_prot_el[N_3p];
	      TVector3 V3_prot_uncorr[N_3p],V3_prot_corr[N_3p],V3_3p_rot[N_3p];
	      double E_cal[N_3p],p_miss_perp[N_3p],P_3pto1p[N_3p];
	      double N_p1[N_3p]={0};
        double N_p_three=0;
        double E_cal_3pto1p[3]={0};
        double p_miss_perp_3pto1p[3]={0};
	      int N_comb = 3;
	      const int N_2p = 2;
	      double E_cal_3pto2p[3][N_2p]={0};
        double p_miss_perp_3pto2p[3][N_2p]={0};
        double P_3pto2p[3][N_2p]={0};
	      TVector3 V3_2p_rot[N_2p], V3_prot_el[N_3p][N_2p];

	      for(int i = 0; i < N_3p; i++)
	      {
          N_p1[i] = 0;

          V4_p_uncorr[i].SetPxPyPzE(pxf[index_p[i]],pyf[index_p[i]],pzf[index_p[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
          V3_prot_uncorr[i] = V4_p_uncorr[i].Vect();


          V3_prot_corr[i].SetXYZ(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60]);
          V4_p_corr[i].SetPxPyPzE(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]+60]*pf[index_p[i]+60]));

	        V4_prot_el[i] = V4_p_corr[i] + V4_el;
	        E_cal[i] = V4_el.E()+ V4_p_corr[i].E() - m_prot + bind_en[ftarget];
	        p_miss_perp[i] = TMath::Sqrt(V4_prot_el[i].Px()*V4_prot_el[i].Px() + V4_prot_el[i].Py()*V4_prot_el[i].Py());
	      }

        V3_prot_el[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
        V3_prot_el[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
        V3_prot_el[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
        V3_prot_el[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
        V3_prot_el[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
        V3_prot_el[2][1]=V4_el.Vect()+V3_prot_uncorr[2];


        rotation->prot3_rot_func(V3_prot_corr,V3_prot_uncorr,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

	      if(num_pi_phot==0 && N_p_three!=0){
//	      if(num_pi == 0 && N_p_three != 0){    //no photons for now F.H. 29.8.19
           double histoweight = 1./Mott_cross_sec; //Weight for 3protons, 1 electron, GENIE weight and Mott cross section
	         for(int count = 0; count < N_comb; count++)    { //Loop over number of combinations
               for(int j = 0; j < N_2p; j++)    { //loop over two protons

 //-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------

                  h1_E_tot_3pto2p->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*histoweight);
		              h1_E_rec_3pto2p->Fill(E_rec, P_3pto2p[count][j]*histoweight);
		              h2_Erec_pperp_321p->Fill(p_miss_perp_3pto2p[count][j],E_rec,P_3pto2p[count][j]*histoweight);
		              h1_E_tot_3pto2p_fracfeed->Fill((E_cal_3pto2p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto2p[count][j]*histoweight);
		              h1_E_rec_3pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_3pto2p[count][j]*histoweight);
		              h2_pperp_W->Fill(W_var,p_miss_perp_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
		              h1_theta0->Fill((V4_beam.Vect()).Angle(V3_prot_el[count][j])*TMath::RadToDeg(),P_3pto2p[count][j]*histoweight);
	                h2_Ecal_Eqe->Fill(E_rec,E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
	                h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
	                h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);

		              for(int i = 0; i < n_slice; i++)
		              {
		                  if (p_miss_perp_3pto2p[count][j]<pperp_max[i] && p_miss_perp_3pto2p[count][j]>pperp_min[i]){
	                         h1_Etot_3pto2p_slice[i]->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*histoweight);
			                     h1_Erec_3pto2p_slice[i]->Fill(E_rec, P_3pto2p[count][j]*histoweight);
		                  }
		              }
               } //end loop over protons
	         } //end loop over combination N_comb

 //-----------------------------------------  3p to 1p  -----------------------------------------------------------------------
	         for(int j = 0; j < N_3p; j++)    {

		             P_3pto1p[j]= N_p1[j]/N_p_three;
		             h1_E_tot_3pto1p->Fill(E_cal[j], P_3pto1p[j]*histoweight);
		             h1_E_rec_3pto1p->Fill(E_rec,P_3pto1p[j]*histoweight);
		             h2_Erec_pperp_31p->Fill(p_miss_perp[j],E_rec,P_3pto1p[j]*histoweight);
		             h1_E_tot_3pto1p_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto1p[j]*histoweight);
		             h1_E_rec_3pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_3pto1p[j]*histoweight);
		             h2_pperp_W->Fill(W_var,p_miss_perp[j],-P_3pto1p[j]*histoweight);
		             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),-P_3pto1p[j]*histoweight);
		             h2_Ecal_Eqe->Fill(E_rec,E_cal[j],-P_3pto1p[j]*histoweight);
		             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],-P_3pto1p[j]*histoweight);
		             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],-P_3pto1p[j]*histoweight);

		             for(int i = 0; i < n_slice; i++)
		             {
		                 if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
		                     h1_Etot_3pto1p_slice[i]->Fill(E_cal[j],P_3pto1p[j]*histoweight);
		                     h1_Erec_3pto1p_slice[i]->Fill(E_rec,P_3pto1p[j]*histoweight);
		                 }
		             }
	         } //end loop over N_3p
	      }//end if num_pi_phot==0 && N_p_three!=0, no pions

//----------------------------------3p 1pi ----------------------------------------------------------

        if (num_pi_phot==1) { //no photons for now F.H. 29.8.19
  //      if (num_pi == 1) { //number of pions  = 1

             double P_tot_3p[N_3p]={0};
             double Ecal_3p1pi[N_3p]={0};
             double p_miss_perp_3p1pi[N_3p]={0};
             double charge = 0;
             TVector3 V3_pi_corr;

             if (num_pimi == 1 && num_pipl == 0) {
               charge = -1;
             }
             else if (num_pipl == 1 && num_pimi == 0) {
               charge = 1;
             }
             //radiation status is false if real photon
             else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
               charge = 0;
             }
             //skip radiation photon
             else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
	              charge = 0;
             }
             else {  std::cout << "WARNING: 3p 1pion event: No charge for one pion/photon could be assigned. Pion number " << ind_pi_phot[0] << std::endl; continue; }

             V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);

             rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot_uncorr, V3_pi_corr, charge ,ec_radstat_n[0], V4_el,  Ecal_3p1pi,p_miss_perp_3p1pi, P_tot_3p);

             double histoweight = 1./Mott_cross_sec; //Weight for 3protons, 1 pion, 1 electron, GENIE weight and Mott cross section

	           for(int j = 0; j < N_3p; j++)    { //loop over 3 protons

                 h1_E_tot_3p1pi->Fill(E_cal[j], P_tot_3p[j]*histoweight);
		             h1_E_rec_3p1pi->Fill(E_rec,P_tot_3p[j]*histoweight);
		             h2_Erec_pperp_3p1pi->Fill(p_miss_perp[j],E_rec,P_tot_3p[j]*histoweight);
		             h1_E_tot_3p1pi_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_tot_3p[j]*histoweight);
		             h1_E_rec_3p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_3p[j]*histoweight);
		             h2_pperp_W->Fill(W_var,p_miss_perp[j],P_tot_3p[j]*histoweight);
		             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_tot_3p[j]*histoweight);
		             h2_Ecal_Eqe->Fill(E_rec,E_cal[j],P_tot_3p[j]*histoweight);
		             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],P_tot_3p[j]*histoweight);
		             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],P_tot_3p[j]*histoweight);

		             for(int i = 0; i < n_slice; i++)
		             {
		                 if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
                         h1_Etot_3p1pi_slice[i]->Fill(E_cal[j],P_tot_3p[j]*histoweight);
		                     h1_Erec_3p1pi_slice[i]->Fill(E_rec,P_tot_3p[j]*histoweight);
		                 }
		             }
             } //end loop over N_3p
        }// 1 pi requirement ends

	  }//end if num_p == 3  3proton requiremnet

//Events with exactly 4 protons
	 	 if(num_p==4){

	      const int N_p4=4;
	      TLorentzVector V4_p4_uncorr[N_p4], V4_p4_corr[N_p4],V4_prot4_el[N_p4];
	      TVector3 V3_prot4_uncorr[N_p4],V3_prot4_corr[N_p4];
	      double E_cal_p4[N_p4]={0};
        double p_miss_perp_p4[N_p4]={0};
        double P_4pto1p[N_p4]={0};
	      TVector3 V3_p4_rot[N_p4];
	      bool prot4_stat[N_p4]={false};
	      const int Ncomb_4to1 = 4,Ncomb_4to2 = 6, Ncomb_4to3 = 4;
	      double N_p4_p1[Ncomb_4to1]={0};
        double N_p4_p2[Ncomb_4to2]={0};
        double N_p4_p3[Ncomb_4to3]={0};
        double N_p_four = 0;

	      for(int i = 0; i < N_p4; i++)  //loop over 4 protons
	      {

          V4_p4_uncorr[i].SetPxPyPzE(pxf[index_p[i]],pyf[index_p[i]],pzf[index_p[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
          V3_prot4_uncorr[i] = V4_p4_uncorr[i].Vect();


          V3_prot4_corr[i].SetXYZ(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60]);
          V4_p4_corr[i].SetPxPyPzE(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60], TMath::Sqrt(m_prot*m_prot+pf[index_p[i]+60]*pf[index_p[i]+60]));

          V4_prot4_el[i] = V4_p4_corr[i] + V4_el;
          E_cal_p4[i] = V4_el.E() + V4_p4_corr[i].E() - m_prot + bind_en[ftarget];
          p_miss_perp_p4[i] = TMath::Sqrt(V4_prot4_el[i].Px()*V4_prot4_el[i].Px() + V4_prot4_el[i].Py()*V4_prot4_el[i].Py());
        }

        if ( num_pi_phot == 0){ //no pion or photon
	            for(int g = 0; g < N_tot; g++){ //this looks like a 4-proton rotation function -> could be placed maybe in an extra function

		              double rot_angle = gRandom->Uniform(0,2*TMath::Pi());
		              for(int i = 0; i < N_p4;i++) {
		                  V3_p4_rot[i]= V3_prot4_uncorr[i];
		                  V3_p4_rot[i].Rotate(rot_angle,V3_q);
		              }
		              for(int i_p = 0; i_p < N_p4; i_p++) {
                    prot4_stat[i_p] = PFiducialCut(fbeam_en, V3_p4_rot[i_p]);
                  }

		              if( prot4_stat[0]  && !prot4_stat[1]   && !prot4_stat[2] && !prot4_stat[3])  N_p4_p1[0]=N_p4_p1[0]+1;//Detecting 1p out of 4p
		              if(!prot4_stat[0]  &&   prot4_stat[1]  && !prot4_stat[2] && !prot4_stat[3])  N_p4_p1[1]=N_p4_p1[1]+1;
		              if(!prot4_stat[0]  &&  !prot4_stat[1]  &&  prot4_stat[2] && !prot4_stat[3])  N_p4_p1[2]=N_p4_p1[2]+1;
		              if(!prot4_stat[0]  &&  !prot4_stat[1]  && !prot4_stat[2] &&  prot4_stat[3])  N_p4_p1[3]=N_p4_p1[3]+1;
		              if( prot4_stat[0]  &&  prot4_stat[1]   &&  prot4_stat[2] &&  prot4_stat[3])  N_p_four=N_p_four+1;   //Detecting 4p out of 4p

		              if( prot4_stat[0]  &&  prot4_stat[1]  && !prot4_stat[2]  && !prot4_stat[3])  N_p4_p2[0]=N_p4_p2[0]+1;//Detecting 2p out of 4p
		              if( prot4_stat[0]  && !prot4_stat[1]  &&  prot4_stat[2]  && !prot4_stat[3])  N_p4_p2[1]=N_p4_p2[1]+1;
		              if( prot4_stat[0]  && !prot4_stat[1]  && !prot4_stat[2]  &&  prot4_stat[3])  N_p4_p2[2]=N_p4_p2[2]+1;
		              if(!prot4_stat[0]  &&  prot4_stat[1]  &&  prot4_stat[2]  && !prot4_stat[3])  N_p4_p2[3]=N_p4_p2[3]+1;
		              if(!prot4_stat[0]  &&  prot4_stat[1]  && !prot4_stat[2]  &&  prot4_stat[3])  N_p4_p2[4]=N_p4_p2[4]+1;
		              if(!prot4_stat[0]  && !prot4_stat[1]  &&  prot4_stat[2]  &&  prot4_stat[3])  N_p4_p2[5]=N_p4_p2[5]+1;


		              if( prot4_stat[0]  &&  prot4_stat[1]  &&  prot4_stat[2]  && !prot4_stat[3])  N_p4_p3[0]=N_p4_p3[0]+1;//Detecting 3p out of 4p
		              if( prot4_stat[0]  &&  prot4_stat[1]  &&  !prot4_stat[2] &&  prot4_stat[3])  N_p4_p3[1]=N_p4_p3[1]+1;
		              if( prot4_stat[0]  && !prot4_stat[1]  &&  prot4_stat[2]  &&  prot4_stat[3])  N_p4_p3[2]=N_p4_p3[2]+1;
		              if(!prot4_stat[0]  &&  prot4_stat[1]  &&  prot4_stat[2]  &&  prot4_stat[3])  N_p4_p3[3]=N_p4_p3[3]+1;
	            }//for loop of 4p rotations ends

                // still no pions
              int N_comb=3;    //number of 2 proton combination out of three
	            const int N_2p=2, N_3p=3;
	            double E_cal_4pto3p[3][N_2p]={0};
              double p_miss_perp_4pto3p[3][N_2p]={0};
              double P_4pto3p[3][N_2p]={0};
	            TVector3 V3_prot_corr[N_3p],V3_prot_uncorr[N_3p],V3_prot_el_4pto3p[N_3p][N_2p],V3_el_prot[N_comb][N_2p];
	            double N_p_three=0,N_p1[N_3p]={0};
              double E_cal_43pto1p[N_3p];
              double p_miss_perp_43pto1p[N_3p];
              double P_43pto1p[3]={0};

              double histoweight = 1./Mott_cross_sec; //Weight for 3protons, 1 electron, GENIE weight and Mott cross section


	            for(int g = 0; g < Ncomb_4to3; g++){   //estimating the undetected 4p contribution to  3p
		             if(g==0) {
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[0]; V3_prot_uncorr[1]=V3_prot4_uncorr[1]; V3_prot_uncorr[2]=V3_prot4_uncorr[2];
		                 V3_prot_corr[0]=V3_prot4_corr[0]; V3_prot_corr[1]=V3_prot4_corr[1]; V3_prot_corr[2]=V3_prot4_corr[2];
		             }
		             if(g==1){
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[0]; V3_prot_uncorr[1]=V3_prot4_uncorr[1]; V3_prot_uncorr[2]=V3_prot4_uncorr[3];
		                 V3_prot_corr[0]=V3_prot4_corr[0]; V3_prot_corr[1]=V3_prot4_corr[1]; V3_prot_corr[2]=V3_prot4_corr[3];
		             }
		             if(g==2){
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[0]; V3_prot_uncorr[1]=V3_prot4_uncorr[2]; V3_prot_uncorr[2]=V3_prot4_uncorr[3];
		                 V3_prot_corr[0]=V3_prot4_corr[0]; V3_prot_corr[1]=V3_prot4_corr[2]; V3_prot_corr[2]=V3_prot4_corr[3];
		             }
		             if(g==3){
		                 V3_prot_uncorr[0]=V3_prot4_uncorr[1]; V3_prot_uncorr[1]=V3_prot4_uncorr[2]; V3_prot_uncorr[2]=V3_prot4_uncorr[3];
		                 V3_prot_corr[0]=V3_prot4_corr[1]; V3_prot_corr[1]=V3_prot4_corr[2]; V3_prot_corr[2]=V3_prot4_corr[3];
		             }

		             rotation->prot3_rot_func(V3_prot_corr, V3_prot_uncorr,V4_el,E_cal_4pto3p,p_miss_perp_4pto3p, P_4pto3p,N_p1,E_cal_43pto1p,p_miss_perp_43pto1p,&N_p_three);

		             V3_el_prot[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
		             V3_el_prot[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
		             V3_el_prot[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
		             V3_el_prot[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
		             V3_el_prot[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
		             V3_el_prot[2][1]=V4_el.Vect()+V3_prot_uncorr[2];

		             if( N_p_three!=0 && N_p_four!=0){
		   	             for(int count = 0; count < N_comb; count++)    { //looping through number of 2 proton combination out of 3 protons
	                      for(int j = 0;j < N_2p; j++)    {  //looping through number of 1 proton combination out of 2 protons

	//-----------------------------------------  4p to 3p->2->1  -----------------------------------------------------------------------

		                      h1_E_tot_4pto3p->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_E_rec_4pto3p->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_Erec_pperp_4321p->Fill(p_miss_perp_4pto3p[count][j],E_rec,P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_E_tot_4pto3p_fracfeed->Fill((E_cal_4pto3p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_E_rec_4pto3p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_pperp_W->Fill(W_var,p_miss_perp_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h1_theta0->Fill((V4_beam.Vect()).Angle(V3_el_prot[count][j])*TMath::RadToDeg(),-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
		                      h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);

		                      for(int i = 0; i < n_slice; i++)
			                    {
			                         if (p_miss_perp_4pto3p[count][j]<pperp_max[i] && p_miss_perp_4pto3p[count][j]>pperp_min[i]){
			                              h1_Etot_4pto3p_slice[i]->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
			                              h1_Erec_4pto3p_slice[i]->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
			                         }
			                    }
		                    } //end loop over N_2p
		                 } //end loop over N_comb : number of 2 proton combination out of 3 protons

		   //-----------------------------------------  4p to 3p->1p  -----------------------------------------------------------------------
		   	             for(int j = 0; j < N_3p; j++)    { //4p to 3p->1, looping through 1p out of 3p
                        //P_43pto1p doesnt have to be an array, one local variable here
		                     P_43pto1p[j]= N_p1[j]/N_p_three*(N_p4_p3[g]/N_p_four);
		                     h1_E_tot_43pto1p->Fill(E_cal_43pto1p[j], P_43pto1p[j]*histoweight);
		                     h1_E_rec_43pto1p->Fill(E_rec,P_43pto1p[j]*histoweight);
		                     h2_Erec_pperp_431p->Fill(p_miss_perp_43pto1p[j],E_rec,P_43pto1p[j]*histoweight);
		                     h1_E_tot_43pto1p_fracfeed->Fill((E_cal_43pto1p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_43pto1p[j]*histoweight);
		                     h1_E_rec_43pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_43pto1p[j]*histoweight);
		                     h2_pperp_W->Fill(W_var,p_miss_perp_43pto1p[j],P_43pto1p[j]*histoweight);
  		                   h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_43pto1p[j]*histoweight);
		                     h2_Ecal_Eqe->Fill(E_rec,E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
		                     h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
		                     h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_43pto1p[j],P_43pto1p[j]*histoweight);

                         for(int i = 0; i <n_slice; i++)
		                     {
			                        if (p_miss_perp_43pto1p[j]<pperp_max[i] && p_miss_perp_43pto1p[j]>pperp_min[i]){
			                             h1_Etot_43pto1p_slice[i]->Fill(E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
			                             h1_Erec_43pto1p_slice[i]->Fill(E_rec,P_43pto1p[j]*histoweight);
			                        }
		                     }
                    } // end loop over N_3p
		             }//end of N_p_three and N_p_four !=0
	            }//end of the loop through 3p combinations out of 4, g < Ncomb_4to3

              //still no pions or photons num_pi_phot == 0
	            int N_4to2=0;
      	      TVector3 V3p2[2],V3p2_uncorr[2];
	            double E_cal_4pto2p[2]={0};
              double p_miss_perp_4pto2p[2]={0};
              double P_4pto2p[2]={0};
              double N_two=0;

  //-----------------------------------------  4p to 2p->1  -----------------------------------------------------------------------
	            for(int ind1 = 0; ind1 < N_p4; ind1++){          //estimating the undetected 4p contribution to  2p
		              for(int ind2 = 0; ind2 < N_p4; ind2++){
		                  if(ind1!=ind2 && ind1 < ind2){

		                      V3p2[0]=V3_prot4_corr[ind1];
		                      V3p2[1]=V3_prot4_corr[ind2];
		                      V3p2_uncorr[0]=V3_prot4_uncorr[ind1];
		                      V3p2_uncorr[1]=V3_prot4_uncorr[ind2];

		                      rotation->prot2_rot_func(V3p2,V3p2_uncorr, V4_el,E_cal_4pto2p,p_miss_perp_4pto2p,  P_4pto2p, &N_two);

		                      if( N_two!=0  && N_p_four!=0){
		                          for(int j = 0; j < N_2p; j++)  {  //looping through  1 proton combination out of 2 protons

	                               h1_E_tot_4pto2p->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_E_rec_4pto2p->Fill(E_rec, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h2_Erec_pperp_421p->Fill( p_miss_perp_4pto2p[j],E_rec,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_E_tot_4pto2p_fracfeed->Fill((E_cal_4pto2p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_E_rec_4pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h2_pperp_W->Fill(W_var,p_miss_perp_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3p2_uncorr[j])*TMath::RadToDeg(),P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
                                 h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
	                               h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
		                             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);

	                               for(int i = 0; i < n_slice; i++)
	                               {
	                                  if (p_miss_perp_4pto2p[j]<pperp_max[i] && p_miss_perp_4pto2p[j]>pperp_min[i]){
	                                     h1_Etot_4pto2p_slice[i]->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
		                                   h1_Erec_4pto2p_slice[i]->Fill(E_rec, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
		                                }
		                             }
		                          } //end loop over N_2p
                          } //end if N_two!=0  && N_p_four!=0
		                      N_4to2= N_4to2+1;
                      } //end if ind1!=ind2 && ind1 < ind2
                  } //end loop over ind2
	            } //end loop over ind1

 //-----------------------------------------  4p to 1p  -----------------------------------------------------------------------
	            if( N_p_four!=0){
		              for(int j = 0; j < N_p4; j++)    {       //estimating the undetected 4p contribution to  1p
                      //P_4pto1p[j] doesnt have to be an array since it is only used here as a local variable
		                  P_4pto1p[j]= N_p4_p1[j]/N_p_four;
		                  h1_E_tot_4pto1p->Fill(E_cal_p4[j], P_4pto1p[j]*histoweight);
		                  h1_E_rec_4pto1p->Fill(E_rec,P_4pto1p[j]*histoweight);
		                  h2_Erec_pperp_41p->Fill(p_miss_perp_p4[j],E_rec, P_4pto1p[j]*histoweight);
		                  h1_E_tot_4pto1p_fracfeed->Fill((E_cal_p4[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto1p[j]*histoweight);
		                  h1_E_rec_4pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_4pto1p[j]*histoweight);
		                  h2_pperp_W->Fill(W_var,p_miss_perp_p4[j],-P_4pto1p[j]*histoweight);
		                  h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot4_uncorr[j])*TMath::RadToDeg(),-P_4pto1p[j]*histoweight);
		                  h2_Ecal_Eqe->Fill(E_rec,E_cal_p4[j],-P_4pto1p[j]*histoweight);
		                  h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_p4[j],-P_4pto1p[j]*histoweight);
		                  h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_p4[j],-P_4pto1p[j]*histoweight);

		                  for(int i = 0; i < n_slice; i++)
		                  {
		                      if (p_miss_perp_p4[j]<pperp_max[i] && p_miss_perp_p4[j]>pperp_min[i]){
                              h1_Etot_4pto1p_slice[i]->Fill(E_cal_p4[j],P_4pto1p[j]*histoweight);
			                        h1_Erec_4pto1p_slice[i]->Fill(E_rec,P_4pto1p[j]*histoweight);
		                      }
		                  }
                  } //end loop over N_p4
	            } // end if N_p_four!=0
        }//no pion statment ends

	   }//4 proton requirement (num_p == 4)

    //We are not looking for 4 Proton and 1 Pion events!

//No Protons here, Next 150 lines are for the inclusive events

     h1_E_rec->Fill(E_rec,WeightIncl);
     if(num_pi_phot ==0){
     //if(num_pi == 0){    //no photons for now F.H. 29.8.19
       h1_E_rec_0pi->Fill(E_rec,WeightIncl);
	     h1_E_rec_0pi_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],WeightIncl);
     }

//----------------------------- e- ,1pi  -----------------------------------------

     if(num_pi_phot == 1){
//     if(num_pi == 1){    //no photons for now F.H. 29.8.19

        TVector3 V3_pi_corr;
        double P_undet=0;
        int charge = 0;
        double pion_acc_ratio = 0;

        if (num_pimi == 1 && num_pipl == 0) {
          charge = -1;
        }
        else if (num_pipl == 1 && num_pimi == 0) {
          charge = 1;
        }
        //radiation status is false if real photon
        else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
          charge = 0;
        }
        //skip radiation photon
        else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
	         charge = 0;
        }
        //skip radiation photon
        else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
          continue;
        }
        else {  std::cout << "WARNING: 1pion events: No charge for one pion/photon could be assigned.  "  << std::endl; continue; }

        V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]], pyf[ind_pi_phot[0]], pzf[ind_pi_phot[0]]);

        rotation->pi1_rot_func(V3_pi_corr, charge, ec_radstat_n[0], &P_undet);

        double histoweight = WeightIncl;

        h1_E_rec_1pi_weight->Fill(E_rec,P_undet*histoweight);
	      h1_E_rec_1pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_undet*histoweight);

	      if(!ec_radstat_n[0])  h1_E_rec_1pi->Fill(E_rec,histoweight);
	      if(ec_num_n==1)       h2_phot_e_angle_Erec->Fill(E_rec,V3_pi_corr.Angle(V4_el.Vect())*TMath::RadToDeg());
	   }
//----------------------------- e- ,2pi  -----------------------------------------

     if(num_pi_phot == 2){
//     if(num_pi == 2){    //no photons for now F.H. 29.8.19

	      const int N_2pi = 2;
	      TVector3 V3_2pi_corr[N_2pi];
	      int q_pi2[N_2pi];
        bool radstat_pi2[N_2pi] = {false};
	      double P_1pi[N_2pi] = {0};
        double P_0pi = 0;

        for (int i = 0; i < num_pi_phot; i++) {

            radstat_pi2[i] = ec_radstat_n[i];
            if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
              q_pi2[i] = 1;
            }
            else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
              q_pi2[i] = -1;
            }
            else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
              q_pi2[i] = 0;
            }
            else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
              continue;
            }
            else {  std::cout << "WARNING: 2pion event: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

            V3_2pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);

        }

        rotation->pi2_rot_func(V3_2pi_corr, q_pi2,radstat_pi2, &P_0pi,P_1pi);

        double histoweight = 1./Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19

//----------------------------- e- ,2pi->0pi (-) -----------------------------------------
        h1_E_rec_2pi_weight->Fill(E_rec,(-P_0pi)*histoweight);
        h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*histoweight);
        h1_E_rec_20pi->Fill(E_rec,(P_0pi)*histoweight);
//----------------------------- e- ,2pi->1pi->0pi (+)  -----------------------------------------
        for(int k = 0; k < N_2pi; k++){ //loop over two pions
          h1_E_rec_2pi_weight->Fill(E_rec,P_1pi[k]*histoweight);
          h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[k]*histoweight);
          h1_E_rec_21pi->Fill(E_rec,(P_1pi[k])*histoweight);
        }
     } //end if for two pion events

//----------------------------- e- ,3pi  -----------------------------------------

     if(num_pi_phot == 3){
//     if(num_pi == 3){    //no photons for now F.H. 29.8.19

        const int N_3pi=3;
        const int N_2pi=2;
        TVector3 V3_3pi_corr[N_3pi];
        int q_pi3[N_3pi]={0};
        bool radstat_pi3[N_3pi]={false};
        double P_0pi = 0;
        double P_1pi[N_3pi]={0};
        double P_320pi[N_3pi]={0};
        double P_3210pi[N_3pi][N_2pi]={0};

        for (int i = 0; i < num_pi_phot; i++) {

            radstat_pi3[i] = ec_radstat_n[i];
            if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
              q_pi3[i] = 1;
            }
            else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
              q_pi3[i] = -1;
            }
            else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
              q_pi3[i] = 0;
            }
            else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
              q_pi3[i] = 0;
            }
            else {  std::cout << "WARNING: 3pion event: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

            V3_3pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);

        }

        rotation->pi3_rot_func(V3_3pi_corr, q_pi3,radstat_pi3, &P_0pi, P_1pi, P_320pi,P_3210pi);

        double histoweight = 1./Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19

 //---------------------------3pi->0pi----------------------------------------------
        h1_E_rec_3pi_weight->Fill(E_rec,(-P_0pi)*histoweight);
        h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*histoweight);
        h1_E_rec_30pi->Fill(E_rec,(P_0pi)*histoweight);

        for(int h = 0; h < N_3pi; h++){ //loop over three pions

//---------------------------3pi->1pi->0pi----------------------------------------------
          h1_E_rec_3pi_weight->Fill(E_rec,P_1pi[h]*histoweight);
          h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[h]*histoweight);
          h1_E_rec_310pi->Fill(E_rec,(P_1pi[h])*histoweight);
 //---------------------------3pi->2pi->0pi----------------------------------------------
          h1_E_rec_3pi_weight->Fill(E_rec,P_320pi[h]*histoweight);
          h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_320pi[h]*histoweight);
          h1_E_rec_320pi->Fill(E_rec,(P_320pi[h])*histoweight);
//---------------------------3pi->2pi->1pi->0pi----------------------------------------------
          for(int g = 0; g < N_2pi; g++){ //loop over two pions
            h1_E_rec_3pi_weight->Fill(E_rec,(-P_3210pi[h][g])*histoweight);
            h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_3210pi[h][g])*histoweight);
            h1_E_rec_3210pi->Fill(E_rec,(P_3210pi[h][g])*histoweight);
          }
        }//end of 3pi loop
     }//end of 3pi requirement

//----------------------------- e- ,4pi  -----------------------------------------
     if(num_pi_phot == 4){
//     if(num_pi == 4){    //no photons for now F.H. 29.8.19

       const int N_4pi=4;
       TVector3 V3_4pi_corr[N_4pi];
       int q_pi4[N_4pi]={0};
       bool  radstat_pi4[N_4pi]={false};
       double P_0pi=0;
       double P_410pi=0;
       double P_420pi=0;
       double P_4210pi=0;
       double P_430pi=0;
       double P_4310pi=0;
       double P_4320pi=0;
       double P_43210pi=0;

       for (int i = 0; i < num_pi_phot; i++) {

           radstat_pi4[i] = ec_radstat_n[i];
           if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
             q_pi4[i] = 1;
           }
           else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
             q_pi4[i] = -1;
           }
           else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
             q_pi4[i] = 0;
           }
           else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
             q_pi4[i] = 0;
           }
           else {  std::cout << "WARNING: 4pion event: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }


           V3_4pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);


       }

       rotation->pi4_rot_func( V3_4pi_corr, q_pi4, radstat_pi4, &P_0pi,&P_410pi,&P_420pi,&P_4210pi,&P_430pi,&P_4310pi,&P_4320pi,&P_43210pi);

       double histoweight = 1./Mott_cross_sec; //Is this correct in the following loop? F.H. 09/01/19


 //---------------------------4pi->0pi----------------------------------------------
//why is it here not split like for 3pi case, sum over all weights is done here F.H 04.08.19
       h1_E_rec_4pi_weight->Fill(E_rec,(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight);
       h1_E_rec_4pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight);
       h1_E_rec_40pi->Fill(E_rec,(P_0pi)*histoweight);
//---------------------------4pi->1pi->0pi----------------------------------------------
       h1_E_rec_410pi->Fill(E_rec,(P_410pi)*histoweight);
//---------------------------4pi->2pi->0pi----------------------------------------------
       h1_E_rec_420pi->Fill(E_rec,(P_420pi)*histoweight);
//---------------------------4pi->2pi->1pi->0pi----------------------------------------------
       h1_E_rec_4210pi->Fill(E_rec,(P_4210pi)*histoweight);
//---------------------------4pi->3pi->0pi----------------------------------------------
       h1_E_rec_430pi->Fill(E_rec,(P_430pi)*histoweight);
//---------------------------4pi->3pi->1pi->0pi----------------------------------------------
       h1_E_rec_4310pi->Fill(E_rec,(P_4310pi)*histoweight);
//---------------------------4pi->3pi->2pi->0pi----------------------------------------------
       h1_E_rec_4320pi->Fill(E_rec,(P_4320pi)*histoweight);
//---------------------------4pi->3pi->2pi->1pi->0pi----------------------------------------------
       h1_E_rec_43210pi->Fill(E_rec,(P_43210pi)*histoweight);

     }//end of 4 pi/photon requirement


   //------------------------------------------requiring there to be a proton -------------------------------------
    //Events with exactly one proton
     if( num_p == 1)
     {

     //Vector for proton without momentum smearing
       TLorentzVector V4_prot_uncorr(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
       TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();

      //Vector for proton with momentum smearing
       TVector3 V3_prot_corr;
       TLorentzVector V4_prot_corr;
       V3_prot_corr.SetXYZ(pxf[index_p[0]+60], pyf[index_p[0]+60], pzf[index_p[0]+60]);
       V4_prot_corr.SetPxPyPzE(pxf[index_p[0]+60], pyf[index_p[0]+60], pzf[index_p[0]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]+60]*pf[index_p[0]+60]));

       double prot_mom_corr = V3_prot_corr.Mag()/V3_prot_uncorr.Mag();
       h1_prot_mom->Fill(prot_mom_corr);
       h1_prot_mom_ratio->Fill(prot_mom_corr/pf[index_p[0]]);

       TLorentzVector V4_prot_el_tot = V4_prot_corr + V4_el;

       double p_perp_tot = TMath::Sqrt(V4_prot_el_tot.Px()*V4_prot_el_tot.Px() + V4_prot_el_tot.Py()*V4_prot_el_tot.Py());
       double E_tot = V4_el.E() + V4_prot_corr.E() -m_prot + bind_en[ftarget];

//These Histograms are events with 1 electron and  1 proton and multiple pions
       double histoweight_inc = 1./Mott_cross_sec;

       h1_E_tot->Fill(E_tot,histoweight_inc);
       h1_E_rec_1prot->Fill(E_rec,histoweight_inc);
       h1_E_tot_1prot->Fill(E_tot,histoweight_inc);
       h2_Erec_pperp->Fill(p_perp_tot,E_rec,histoweight_inc);

       //---------------------------------- 1p 0pi   ----------------------------------------------
       if(num_pi_phot == 0){
  //     if(num_pi == 0){    //no photons for now F.H. 29.8.19

          double histoweight = 1./Mott_cross_sec;

          h2_Erec_pperp_newcut2->Fill(p_perp_tot,E_rec,histoweight);
     	    h1_E_rec_cut2_new->Fill(E_rec,histoweight);
     	    h1_E_tot_cut2->Fill(E_tot,histoweight);
     	    h1_E_tot_cut2_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],histoweight);
     	    h1_E_rec_cut2_new_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],histoweight);
     	    h2_pperp_W->Fill(W_var,p_perp_tot,histoweight);
     	    h1_theta0->Fill((V4_beam.Vect()).Angle(V4_prot_el_tot.Vect()) *TMath::RadToDeg(),histoweight);
     	    h2_Ecal_Eqe->Fill(E_rec,E_tot,histoweight);
     	    h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,histoweight);
     	    h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,histoweight);

          for(int i = 0; i < n_slice; i++)
          {
            if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
              h1_Etot_Npi0[i]->Fill(E_tot,histoweight);
              h1_Erec_Npi0[i]->Fill(E_rec,histoweight);
            }
          }
          if (p_perp_tot < 0.2){
            h1_E_rec_cut005_newcut3->Fill(E_rec,histoweight);
     	      h2_Erec_pperp_cut3->Fill(p_perp_tot,E_rec,histoweight);
     	    }
       }//num pi=0
       //---------------------------------- 1p 1pi   ----------------------------------------------
       if(num_pi_phot == 1){
  //     if(num_pi == 1){    //no photons for now F.H. 29.8.19

          double N_piphot_det;
          double N_piphot_undet;
          TVector3 V3_pi_corr;
          int charge = 0;

          if (num_pimi == 1 && num_pipl == 0) {
            charge = -1;
          }
          else if (num_pipl == 1 && num_pimi == 0) {
            charge = 1;
          }
          //radiation status is false if real photon
          else if (num_pipl == 0 && num_pimi == 0 && (!ec_radstat_n[0]) )  {
            charge = 0;
          }
          //skip radiation photon
          else if (num_pipl == 0 && num_pimi == 0 && ec_radstat_n[0] ) {
	           charge = 0;
          }
          else {  std::cout << "WARNING: 1pion events: No charge for one pion/photon could be assigned.  "  << std::endl; continue; }

          V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]], pyf[ind_pi_phot[0]], pzf[ind_pi_phot[0]]);

          rotation->prot1_pi1_rot_func(V3_prot_uncorr,V3_pi_corr, charge, ec_radstat_n[0], &N_piphot_det,&N_piphot_undet);

          double histoweight = 1./Mott_cross_sec; //1proton, 1 Pion, 1 electron acceptance, GENIE weight and Mott

          if(N_piphot_det!=0){

             h1_E_rec_undetfactor->Fill(E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
             h1_E_tot_undetfactor->Fill(E_tot,(N_piphot_undet/N_piphot_det)*histoweight);
             h2_Erec_pperp_1p1pi->Fill(p_perp_tot,E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
             h1_E_rec_undetfactor_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(N_piphot_undet/N_piphot_det)*histoweight);
             h1_E_tot_undetfactor_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],(N_piphot_undet/N_piphot_det)*histoweight);
             h2_pperp_W->Fill(W_var,p_perp_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-(N_piphot_undet/N_piphot_det)*histoweight);
             h2_Ecal_Eqe->Fill(E_rec,E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);

             for(int i = 0; i < n_slice; i++)
             {
               if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
                 h1_Etot_bkgd_pipl_pimi_fact[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*histoweight);
                 h1_Erec_bkgd_pipl_pimi_new_fact[i]->Fill(E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
                 h1_Etot_bkgd_pipl_pimi_fact_pipl[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*histoweight);
                 h1_Etot_Npi1[i]->Fill(E_tot,histoweight);
                 h1_Erec_Npi1[i]->Fill(E_rec,histoweight);
               }
             }
          } //end of N_piphot_det!=0

          h1_E_rec_cutpi1_piplpimi->Fill(E_rec,histoweight);
          h1_E_tot_cutpi1_piplpimi->Fill(E_tot,histoweight);

       }//end of 1p 1pi requirement

 //---------------------------------- 1p 2pi   ----------------------------------------------
       if(num_pi_phot == 2){
  //     if(num_pi == 2){    //no photons for now F.H. 29.8.19

         const int N_2pi=2;
         TVector3 V3_2pi_corr[N_2pi],V3_2pi_rot[N_2pi],V3_p_rot;
         int q_pi2[N_2pi];
         bool radstat_pi2[N_2pi]={false};
         double P_1p0pi=0;
         double P_1p1pi[N_2pi]={0};

         for (int i = 0; i < num_pi; i++) {

             radstat_pi2[i] = ec_radstat_n[i];
             if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
               q_pi2[i] = 1;
             }
             else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
               q_pi2[i] = -1;
             }
             else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
               q_pi2[i] = 0;
             }
             else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
               q_pi2[i] = 0;
             }
             else {  std::cout << "WARNING: 1p 2pion event: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

             V3_2pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);

         }

         rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi_corr,q_pi2,radstat_pi2,&P_1p0pi,P_1p1pi);

         double histoweight = 1./Mott_cross_sec; //1proton, 2 Pion, 1 electron acceptance, GENIE weight and Mott

 //---------------------------------- 1p 2pi->1p1pi   ----------------------------------------------

         for(int z = 0; z < N_2pi; z++){  //to consider 2 diff. 1pi states

           h1_E_tot_1p2pi->Fill(E_tot,P_1p1pi[z]*histoweight);
           h1_E_rec_1p2pi->Fill(E_rec,P_1p1pi[z]*histoweight);
           h2_Erec_pperp_1p2pi_1p1pi->Fill(p_perp_tot,E_rec,P_1p1pi[z]*histoweight);
           h1_E_tot_1p2pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
           h1_E_rec_1p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p1pi[z]*histoweight);
           h2_pperp_W->Fill(W_var,p_perp_tot,P_1p1pi[z]*histoweight);
           h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p1pi[z]*histoweight);
           h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p1pi[z]*histoweight);
           h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p1pi[z]*histoweight);
           h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p1pi[z]*histoweight);

           for(int i = 0; i < n_slice; i++){
	            if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	               h1_Etot_bkgd_1p2pi[i]->Fill(E_tot,P_1p1pi[z]*histoweight);
	               h1_Erec_bkgd_1p2pi[i]->Fill(E_rec,P_1p1pi[z]*histoweight);
	            }
           }
         } //end loop over N_2pi

 //---------------------------------- 1p 2pi->1p0pi   ----------------------------------------------

         h1_E_tot_1p2pi_1p0pi->Fill(E_tot,P_1p0pi*histoweight);
         h1_E_rec_1p2pi_1p0pi->Fill(E_rec,P_1p0pi*histoweight);
         h2_Erec_pperp_1p2pi_1p0pi->Fill(p_perp_tot,E_rec,P_1p0pi*histoweight);
         h1_E_tot_1p2pi_1p0pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p0pi*histoweight);
         h1_E_rec_1p2pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p0pi*histoweight);
         h2_pperp_W->Fill(W_var,p_perp_tot,-P_1p0pi*histoweight);
         h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-P_1p0pi*histoweight);
         h2_Ecal_Eqe->Fill(E_rec,E_tot,-P_1p0pi*histoweight);
         h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-P_1p0pi*histoweight);
         h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-P_1p0pi*histoweight);

         for(int i = 0; i < n_slice; i++){
           if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	            h1_Etot_bkgd_1p2pi_1p0pi[i]->Fill(E_tot,P_1p0pi*histoweight);
	            h1_Erec_bkgd_1p2pi_1p0pi[i]->Fill(E_rec,P_1p0pi*histoweight);
           }
         }
       }//1p 2pi statetment ends

 //---------------------------------- 1p 3pi   ----------------------------------------------
       if(num_pi_phot == 3){
//       if(num_pi == 3){    //no photons for now F.H. 29.8.19

         const int N_3pi=3;
         TVector3 V3_3pi_corr[N_3pi],V3_3pi_rot[N_3pi],V3_p_rot;
         int q_pi3[N_3pi];
         bool radstat_pi3[N_3pi]={false};
         double P_1p3pi = 0;

	 //         for (int i = 0; i < num_pi; i++) {
         for (int i = 0; i < num_pi_phot; i++) {

             radstat_pi3[i] = ec_radstat_n[i];
             if ( index_pipl[i] == ind_pi_phot[i] ) { //i-th pion is a piplus
               q_pi3[i] = 1;
             }
             else if ( index_pimi[i] == ind_pi_phot[i] ) { //i-th pion is a piminus
               q_pi3[i] = -1;
             }
             else if (ind_pi_phot[i]!= -1 && !ec_radstat_n[i] ) { //i-th particle is a neutral
               q_pi3[i] = 0;
             }
             else if (ind_pi_phot[i]!= -1 && ec_radstat_n[i] ) { //i-th particle is a radiation photon
               q_pi3[i] = 0;
             }
             else {  std::cout << "WARNING: 1p 3pion event: No charge for one pion/photon could be assigned. Pion number " << i << std::endl; continue; }

             V3_3pi_corr[i].SetXYZ(pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);

         }

         rotation->prot1_pi3_rot_func(V3_prot_uncorr,V3_3pi_corr,q_pi3,radstat_pi3,&P_1p3pi);

         double histoweight = 1./Mott_cross_sec; //1proton, 3 Pion, 1 electron acceptance, GENIE weight and Mott


 //---------------------------------- 1p 3pi->1p 0pi  total ?? F.H. 08/13/19 check logic here compared to 1p 2pi case ----------------------------------------------
         h1_E_tot_1p3pi->Fill(E_tot,P_1p3pi*histoweight);
         h1_E_rec_1p3pi->Fill(E_rec,P_1p3pi*histoweight);
         h2_Erec_pperp_1p3pi->Fill(p_perp_tot,E_rec,P_1p3pi*histoweight);
         h1_E_tot_1p3pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi*histoweight);
         h1_E_rec_1p3pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p3pi*histoweight);
         h2_pperp_W->Fill(W_var,p_perp_tot,P_1p3pi*histoweight);
         h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p3pi*histoweight);
         h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p3pi*histoweight);
         h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p3pi*histoweight);
         h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p3pi*histoweight);

         for(int i = 0; i < n_slice; i++){
           if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
             h1_Etot_bkgd_1p3pi[i]->Fill(E_tot,P_1p3pi*histoweight);
             h1_Erec_bkgd_1p3pi[i]->Fill(E_rec,P_1p3pi*histoweight);
           }
         }
       }//end of 1p 3pi requirement


     } // 1proton ends


  } //end of event loop (jentry)

  gStyle->SetOptFit(1);

  for(int i = 0; i <= n_slice-1; i++)
  {
  //------------------------------------using the ratio of the pi- to pi+  --------------------------------------

	   h1_Etot_piplpimi_subtruct_fact[i] = (TH1F*)  h1_Etot_Npi0[i]->Clone(Form("h1_Etot_piplpimi_subtruct_fact_%d",i+1));
	   h1_Etot_piplpimi_subtruct_fact[i]->Add(h1_Etot_bkgd_pipl_pimi_fact[i],-1);
	   h1_Erec_piplpimi_subtruct_fact[i]=(TH1F*)  h1_Erec_Npi0[i]->Clone(Form("h1_Erec_piplpimi_subtruct_fact_%d",i+1));
	   h1_Erec_piplpimi_subtruct_fact[i]->Add(h1_Erec_bkgd_pipl_pimi_new_fact[i],-1);

 //------------------------------------subtracting 2p contribution from 1p events  --------------------------------------

	   h1_Etot_p_bkgd_slice_sub[i]=(TH1F*) h1_Etot_piplpimi_subtruct_fact[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub[i]->Add(h1_Etot_p_bkgd_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub[i]=(TH1F*) h1_Erec_piplpimi_subtruct_fact[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub[i]->Add(h1_Erec_p_bkgd_slice[i],-1);

 //------------------------------------undetected 3 to 2 proton subtraction --------------------------------------

	   h1_Etot_p_bkgd_slice_sub32[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub32_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub32[i]->Add(h1_Etot_3pto2p_slice[i]);
	   h1_Erec_p_bkgd_slice_sub32[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub32_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub32[i]->Add(h1_Erec_3pto2p_slice[i]);

 //------------------------------------undetected 3 to 1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub31[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub32[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub31_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub31[i]->Add(h1_Etot_3pto1p_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub31[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub32[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub31_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub31[i]->Add(h1_Erec_3pto1p_slice[i],-1);

 //------------------------------------undetected 4 to 3->2->1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub43[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub31[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub43_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub43[i]->Add(h1_Etot_4pto3p_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub43[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub31[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub43_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub43[i]->Add(h1_Erec_4pto3p_slice[i],-1);

 //------------------------------------undetected 4 to 3->1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub431[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub43[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub431_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub431[i]->Add(h1_Etot_43pto1p_slice[i]);
	   h1_Erec_p_bkgd_slice_sub431[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub43[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub431_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub431[i]->Add(h1_Erec_43pto1p_slice[i]);

 //------------------------------------undetected 4 to 2 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub42[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub431[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub42_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub42[i]->Add(h1_Etot_4pto2p_slice[i]);
	   h1_Erec_p_bkgd_slice_sub42[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub431[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub42_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub42[i]->Add(h1_Erec_4pto2p_slice[i]);

 //------------------------------------undetected 4 to 1 proton addition --------------------------------------

	   h1_Etot_p_bkgd_slice_sub41[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub42[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub41_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub41[i]->Add(h1_Etot_4pto1p_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub41[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub42[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub41_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub41[i]->Add(h1_Erec_4pto1p_slice[i],-1);

//------------------------------------undetected 1p 2pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub1p2pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub41[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub1p2pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub1p2pi[i]->Add(h1_Etot_bkgd_1p2pi[i]);
	   h1_Erec_p_bkgd_slice_sub1p2pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub41[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub1p2pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub1p2pi[i]->Add(h1_Erec_bkgd_1p2pi[i]);

//------------------------------------undetected 1p 2pi-> 1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub1p2pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub1p2pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub1p2pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub1p2pi_0pi[i]->Add(h1_Etot_bkgd_1p2pi_1p0pi[i],-1);
	   h1_Erec_p_bkgd_slice_sub1p2pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub1p2pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub1p2pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub1p2pi_0pi[i]->Add(h1_Erec_bkgd_1p2pi_1p0pi[i],-1);

//------------------------------------undetected 1p 3pi-> 1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub1p3pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub1p2pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub1p3pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub1p3pi_0pi[i]->Add(h1_Etot_bkgd_1p3pi[i]);
	   h1_Erec_p_bkgd_slice_sub1p3pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub1p2pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub1p3pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub1p3pi_0pi[i]->Add(h1_Erec_bkgd_1p3pi[i]);

//------------------------------------undetected 2p 2pi-> 1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p2pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub1p3pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p2pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p2pi_0pi[i]->Add(h1_Etot_p_bkgd_slice_2p2pi[i]);
	   h1_Erec_p_bkgd_slice_sub2p2pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub1p3pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p2pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p2pi_0pi[i]->Add(h1_Erec_p_bkgd_slice_2p2pi[i]);

//------------------------------------undetected 3p 1pi->1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub3p1pi_0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub2p2pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub3p1pi_0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub3p1pi_0pi[i]->Add(h1_Etot_3p1pi_slice[i]);
	   h1_Erec_p_bkgd_slice_sub3p1pi_0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub2p2pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub3p1pi_0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub3p1pi_0pi[i]->Add(h1_Erec_3p1pi_slice[i]);

//------------------------------------undetected 2p 1pi ->2p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p1pi_2p[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub3p1pi_0pi[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p1pi_2p_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p1pi_2p[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[i]);
	   h1_Erec_p_bkgd_slice_sub2p1pi_2p[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub3p1pi_0pi[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p1pi_2p_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p1pi_2p[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[i]);

//------------------------------------undetected 2p 1pi ->1p 1pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub2p1pi_2p[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p1pi_1p_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[i]);
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub2p1pi_2p[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p1pi_1p_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[i]);

//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

	   h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[i],-1);
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi_%d",i+1));
	   h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[i],-1);

  }

 //------------------------------------fractional energy reconstruction plots --------------------------------------

  //------------------------------------using the ratio of the pi- to pi+  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_factor =(TH1F*)  h1_E_rec_cut2_new->Clone("h_Erec_subtruct_piplpimi_factor");
  h_Erec_subtruct_piplpimi_factor->Add(h1_E_rec_undetfactor,-1);

  TH1F *h_Etot_subtruct_piplpimi_factor=(TH1F*)  h1_E_tot_cut2->Clone("h_Etot_subtruct_piplpimi_factor");
  h_Etot_subtruct_piplpimi_factor->Add(h1_E_tot_undetfactor,-1);

  TH2F *h2_Erec_pperp_1p1pisub=(TH2F*) h2_Erec_pperp_newcut2->Clone("h2_Erec_pperp_1p1pisub");
  h2_Erec_pperp_1p1pisub->Add(h2_Erec_pperp_1p1pi,-1);

  TH1F *h_Erec_subtruct_piplpimi_factor_fracfeed =(TH1F*)  h1_E_rec_cut2_new_fracfeed->Clone("h_Erec_subtruct_piplpimi_factor_fracfeed");
  h_Erec_subtruct_piplpimi_factor_fracfeed->Add(h1_E_rec_undetfactor_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_factor_fracfeed=(TH1F*)  h1_E_tot_cut2_fracfeed->Clone("h_Etot_subtruct_piplpimi_factor_fracfeed");
  h_Etot_subtruct_piplpimi_factor_fracfeed->Add(h1_E_tot_undetfactor_fracfeed,-1);

 //-----------------------------------undetected 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_prot=(TH1F*)  h_Erec_subtruct_piplpimi_factor->Clone("h_Erec_subtruct_piplpimi_prot");
  h_Erec_subtruct_piplpimi_prot->Add(h1_E_rec_p_bkgd,-1);

  TH1F *h_Etot_subtruct_piplpimi_prot=(TH1F*)  h_Etot_subtruct_piplpimi_factor->Clone("h_Etot_subtruct_piplpimi_prot");
  h_Etot_subtruct_piplpimi_prot->Add(h1_E_tot_p_bkgd,-1);

  TH2F *h2_Erec_pperp_2psub=(TH2F*) h2_Erec_pperp_1p1pisub->Clone("h2_Erec_pperp_2psub");
  h2_Erec_pperp_2psub->Add(h2_Erec_pperp_2p,-1);

  TH1F *h_Erec_subtruct_piplpimi_prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_factor_fracfeed->Clone("h_Erec_subtruct_piplpimi_prot_fracfeed");
  h_Erec_subtruct_piplpimi_prot_fracfeed->Add(h1_E_rec_p_bkgd_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_factor_fracfeed->Clone("h_Etot_subtruct_piplpimi_prot_fracfeed");
  h_Etot_subtruct_piplpimi_prot_fracfeed->Add(h1_E_tot_p_bkgd_fracfeed,-1);

 //-----------------------------------undetected 3 to 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_32prot=(TH1F*)  h_Erec_subtruct_piplpimi_prot->Clone("h_Erec_subtruct_piplpimi_32prot");
  h_Erec_subtruct_piplpimi_32prot->Add(h1_E_rec_3pto2p);

  TH1F *h_Etot_subtruct_piplpimi_32prot=(TH1F*)  h_Etot_subtruct_piplpimi_prot->Clone("h_Etot_subtruct_piplpimi_32prot");
  h_Etot_subtruct_piplpimi_32prot->Add(h1_E_tot_3pto2p);

  TH2F *h2_Erec_pperp_32psub=(TH2F*) h2_Erec_pperp_2psub->Clone("h2_Erec_pperp_32psub");
  h2_Erec_pperp_32psub->Add(h2_Erec_pperp_321p);

  TH1F *h_Erec_subtruct_piplpimi_32prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_32prot_fracfeed");
  h_Erec_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_rec_3pto2p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_32prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_32prot_fracfeed");
  h_Etot_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_tot_3pto2p_fracfeed);

 //-----------------------------------undetected 3 to 1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_31prot=(TH1F*)  h_Erec_subtruct_piplpimi_32prot->Clone("h_Erec_subtruct_piplpimi_31prot");
  h_Erec_subtruct_piplpimi_31prot->Add(h1_E_rec_3pto1p,-1);

  TH1F *h_Etot_subtruct_piplpimi_31prot=(TH1F*)  h_Etot_subtruct_piplpimi_32prot->Clone("h_Etot_subtruct_piplpimi_31prot");
  h_Etot_subtruct_piplpimi_31prot->Add(h1_E_tot_3pto1p,-1);

  TH2F *h2_Erec_pperp_31psub=(TH2F*) h2_Erec_pperp_32psub->Clone("h2_Erec_pperp_31psub");
  h2_Erec_pperp_31psub->Add(h2_Erec_pperp_31p,-1);

  TH1F *h_Erec_subtruct_piplpimi_31prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_32prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_31prot_fracfeed");
  h_Erec_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_rec_3pto1p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_31prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_32prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_31prot_fracfeed");
  h_Etot_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_tot_3pto1p_fracfeed,-1);

 //-----------------------------------undetected 4 to 3->2->1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_43prot=(TH1F*)  h_Erec_subtruct_piplpimi_31prot->Clone("h_Erec_subtruct_piplpimi_43prot");
  h_Erec_subtruct_piplpimi_43prot->Add(h1_E_rec_4pto3p,-1);

  TH1F *h_Etot_subtruct_piplpimi_43prot=(TH1F*)  h_Etot_subtruct_piplpimi_31prot->Clone("h_Etot_subtruct_piplpimi_43prot");
  h_Etot_subtruct_piplpimi_43prot->Add(h1_E_tot_4pto3p,-1);

  TH2F *h2_Erec_pperp_43psub=(TH2F*) h2_Erec_pperp_31psub->Clone("h2_Erec_pperp_43psub");
  h2_Erec_pperp_43psub->Add(h2_Erec_pperp_4321p,-1);

  TH1F *h_Erec_subtruct_piplpimi_43prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_31prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_43prot_fracfeed");
  h_Erec_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_rec_4pto3p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_43prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_31prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_43prot_fracfeed");
  h_Etot_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_tot_4pto3p_fracfeed,-1);

 //-----------------------------------undetected 4 to 3->1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_431prot=(TH1F*)  h_Erec_subtruct_piplpimi_43prot->Clone("h_Erec_subtruct_piplpimi_431prot");
  h_Erec_subtruct_piplpimi_431prot->Add(h1_E_rec_43pto1p);

  TH1F *h_Etot_subtruct_piplpimi_431prot=(TH1F*)  h_Etot_subtruct_piplpimi_43prot->Clone("h_Etot_subtruct_piplpimi_431prot");
  h_Etot_subtruct_piplpimi_431prot->Add(h1_E_tot_43pto1p);

  TH2F *h2_Erec_pperp_431psub=(TH2F*) h2_Erec_pperp_43psub->Clone("h2_Erec_pperp_431psub");
  h2_Erec_pperp_431psub->Add(h2_Erec_pperp_431p);

  TH1F *h_Erec_subtruct_piplpimi_431prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_43prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_431prot_fracfeed");
  h_Erec_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_rec_43pto1p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_431prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_43prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_431prot_fracfeed");
  h_Etot_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_tot_43pto1p_fracfeed);

//-----------------------------------undetected 4 to 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_42prot=(TH1F*)  h_Erec_subtruct_piplpimi_431prot->Clone("h_Erec_subtruct_piplpimi_42prot");
  h_Erec_subtruct_piplpimi_42prot->Add(h1_E_rec_4pto2p);

  TH1F *h_Etot_subtruct_piplpimi_42prot=(TH1F*) h_Etot_subtruct_piplpimi_431prot->Clone("h_Etot_subtruct_piplpimi_42prot");
  h_Etot_subtruct_piplpimi_42prot->Add(h1_E_tot_4pto2p);

  TH2F *h2_Erec_pperp_42psub=(TH2F*) h2_Erec_pperp_431psub->Clone("h2_Erec_pperp_42psub");
  h2_Erec_pperp_42psub->Add(h2_Erec_pperp_421p);

  TH1F *h_Erec_subtruct_piplpimi_42prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_431prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_42prot_fracfeed");
  h_Erec_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_rec_4pto2p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_42prot_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_431prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_42prot_fracfeed");
  h_Etot_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_tot_4pto2p_fracfeed);

 //-----------------------------------undetected 4 to 1 proton subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_41prot=(TH1F*)  h_Erec_subtruct_piplpimi_42prot->Clone("h_Erec_subtruct_piplpimi_41prot");
  h_Erec_subtruct_piplpimi_41prot->Add(h1_E_rec_4pto1p,-1);

  TH1F *h_Etot_subtruct_piplpimi_41prot=(TH1F*)  h_Etot_subtruct_piplpimi_42prot->Clone("h_Etot_subtruct_piplpimi_41prot");
  h_Etot_subtruct_piplpimi_41prot->Add(h1_E_tot_4pto1p,-1);

  TH2F *h2_Erec_pperp_41psub=(TH2F*) h2_Erec_pperp_42psub->Clone("h2_Erec_pperp_41psub");
  h2_Erec_pperp_41psub->Add(h2_Erec_pperp_41p,-1);

  TH1F *h_Erec_subtruct_piplpimi_41prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_42prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_41prot_fracfeed");
  h_Erec_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_rec_4pto1p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_41prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_42prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_41prot_fracfeed");
  h_Etot_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_tot_4pto1p_fracfeed,-1);

//------------------------------------undetected 1p 2pi ->1 p1pi ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_1p2pi=(TH1F*)  h_Erec_subtruct_piplpimi_41prot->Clone("h_Erec_subtruct_piplpimi_1p2pi");
  h_Erec_subtruct_piplpimi_1p2pi->Add(h1_E_rec_1p2pi);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi=(TH1F*)  h_Etot_subtruct_piplpimi_41prot->Clone("h_Etot_subtruct_piplpimi_1p2pi");
  h_Etot_subtruct_piplpimi_1p2pi->Add(h1_E_tot_1p2pi);

  TH2F *h2_Erec_pperp_sub_1p2pi_1p1pi=(TH2F*) h2_Erec_pperp_41psub->Clone("h2_Erec_pperp_sub_1p2pi_1p1pi");
  h2_Erec_pperp_sub_1p2pi_1p1pi->Add(h2_Erec_pperp_1p2pi_1p1pi);

  TH1F *h_Erec_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_41prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_fracfeed");
  h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Add(h1_E_rec_1p2pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_41prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Add(h1_E_tot_1p2pi_fracfeed);

//------------------------------------undetected 1p 2pi-> 1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi");
  h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Add(h1_E_rec_1p2pi_1p0pi,-1);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi");
  h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Add(h1_E_tot_1p2pi_1p0pi,-1);

  TH2F *h2_Erec_pperp_sub_1p2pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p1pi->Clone("h2_Erec_pperp_sub_1p2pi_1p0pi");
  h2_Erec_pperp_sub_1p2pi_1p0pi->Add(h2_Erec_pperp_1p2pi_1p0pi,-1);

  TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
  h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(h1_E_rec_1p2pi_1p0pi_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(h1_E_tot_1p2pi_1p0pi_fracfeed,-1);

//------------------------------------undetected 1p 3pi-> 1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_1p3pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Erec_subtruct_piplpimi_1p3pi");
  h_Erec_subtruct_piplpimi_1p3pi->Add(h1_E_rec_1p3pi);

  TH1F *h_Etot_subtruct_piplpimi_1p3pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Etot_subtruct_piplpimi_1p3pi");
  h_Etot_subtruct_piplpimi_1p3pi->Add(h1_E_tot_1p3pi);

  TH2F *h2_Erec_pperp_sub_1p3pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p0pi->Clone("h2_Erec_pperp_sub_1p3pi");
  h2_Erec_pperp_sub_1p3pi->Add(h2_Erec_pperp_1p3pi);

  TH1F *h_Erec_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p3pi_fracfeed");
  h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Add(h1_E_rec_1p3pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p3pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Add(h1_E_tot_1p3pi_fracfeed);


//------------------------------------undetected 2p 2pi ->1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p2pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p3pi->Clone("h_Erec_subtruct_piplpimi_2p2pi");
  h_Erec_subtruct_piplpimi_2p2pi->Add(h1_E_rec_2p2pi);

  TH1F *h_Etot_subtruct_piplpimi_2p2pi=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi->Clone("h_Etot_subtruct_piplpimi_2p2pi");
  h_Etot_subtruct_piplpimi_2p2pi->Add(h1_E_tot_2p2pi);

  TH2F *h2_Erec_pperp_sub_2p2pi=(TH2F*) h2_Erec_pperp_sub_1p3pi->Clone("h2_Erec_pperp_sub_2p2pi");
  h2_Erec_pperp_sub_2p2pi->Add(h2_Erec_pperp_2p2pi);

  TH1F *h_Erec_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p2pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Add(h1_E_rec_2p2pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p2pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Add(h1_E_tot_2p2pi_fracfeed);


//------------------------------------undetected 3p 1pi ->1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_3p1pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p2pi->Clone("h_Erec_subtruct_piplpimi_3p1pi");
  h_Erec_subtruct_piplpimi_3p1pi->Add(h1_E_rec_3p1pi);

  TH1F *h_Etot_subtruct_piplpimi_3p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi->Clone("h_Etot_subtruct_piplpimi_3p1pi");
  h_Etot_subtruct_piplpimi_3p1pi->Add(h1_E_tot_3p1pi);

  TH2F *h2_Erec_pperp_sub_3p1pi=(TH2F*) h2_Erec_pperp_sub_2p2pi->Clone("h2_Erec_pperp_sub_3p1pi");
  h2_Erec_pperp_sub_3p1pi->Add(h2_Erec_pperp_3p1pi);

  TH1F *h_Erec_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_3p1pi_fracfeed");
  h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Add(h1_E_rec_3p1pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_3p1pi_fracfeed");
  h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Add(h1_E_tot_3p1pi_fracfeed);

//------------------------------------undetected 2p 1pi ->2p 0pi  --------------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_3p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi");
  h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Add(h1_E_rec_2p1pi_2p0pi);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi");
  h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Add(h1_E_tot_2p1pi_2p0pi);

  TH2F *h2_Erec_pperp_sub_2p1pi_2p0pi=(TH2F*) h2_Erec_pperp_sub_3p1pi->Clone("h2_Erec_pperp_sub_2p1pi_2p0pi");
  h2_Erec_pperp_sub_2p1pi_2p0pi->Add(h2_Erec_pperp_2p1pi_2p0pi);

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(h1_E_rec_2p1pi_2p0pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(h1_E_tot_2p1pi_2p0pi_fracfeed);

//------------------------------------undetected 2p 1pi ->1p 1pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi");
  h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Add(h1_E_rec_2p1pi_1p1pi);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi");
  h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Add(h1_E_tot_2p1pi_1p1pi);

  TH2F *h2_Erec_pperp_sub_2p1pi_1p1pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_2p0pi->Clone("h2_Erec_pperp_sub_2p1pi_1p1pi");
  h2_Erec_pperp_sub_2p1pi_1p1pi->Add(h2_Erec_pperp_2p1pi_1p1pi);

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_rec_2p1pi_1p1pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_tot_2p1pi_1p1pi_fracfeed);

//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi");
  h_Erec_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_rec_2p1pi_1p0pi,-1);

// apapadop
//  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi");
  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("epRecoEnergy_slice_0");
  h_Etot_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_tot_2p1pi_1p0pi,-1);

  TH2F *h2_Erec_pperp_sub_2p1pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_1p1pi->Clone("h2_Erec_pperp_sub_2p1pi_1p0pi");
  h2_Erec_pperp_sub_2p1pi_1p0pi->Add(h2_Erec_pperp_2p1pi_1p0pi,-1);

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(h1_E_rec_2p1pi_1p0pi_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(h1_E_tot_2p1pi_1p0pi_fracfeed,-1);


 //-----------------------------------looking only at e-, 1pi, undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot = (TH1F*)  h1_E_rec_0pi->Clone("h_Erec_subtruct_piplpimi_noprot");
  h_Erec_subtruct_piplpimi_noprot->Add(h1_E_rec_1pi_weight,-1);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed = (TH1F*)  h1_E_rec_0pi_frac_feed->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed");
  h_Erec_subtruct_piplpimi_noprot_frac_feed->Add(h1_E_rec_1pi_weight_frac_feed,-1);
 //-----------------------------------looking only at e-, 2pi undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot_2pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot->Clone("h_Erec_subtruct_piplpimi_noprot_2pi");
  h_Erec_subtruct_piplpimi_noprot_2pi->Add(h1_E_rec_2pi_weight);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed2pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_frac_feed->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed2pi");
  h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Add(h1_E_rec_2pi_weight_frac_feed);

 //-----------------------------------looking only at e-, 3pi, undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot_3pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_2pi->Clone("h_Erec_subtruct_piplpimi_noprot_3pi");
  h_Erec_subtruct_piplpimi_noprot_3pi->Add(h1_E_rec_3pi_weight);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed3pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed3pi");
  h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Add(h1_E_rec_3pi_weight_frac_feed);

 //-----------------------------------looking only at e-, 4pi, undetected pion subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_noprot_4pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_3pi->Clone("h_Erec_subtruct_piplpimi_noprot_4pi");
  h_Erec_subtruct_piplpimi_noprot_4pi->Add(h1_E_rec_4pi_weight);

  TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed4pi = (TH1F*)  h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed4pi");
  h_Erec_subtruct_piplpimi_noprot_frac_feed4pi->Add(h1_E_rec_4pi_weight_frac_feed);

  gDirectory->Write("hist_Files", TObject::kOverwrite);
  // skim_tree->AutoSave();

}

//End LoopCLAS function
