
#define E2A_EP_C

#include "e2a_ep_neutrino6_united4_radphot.h"
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
#include <TGraph.h>

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en

using namespace std;

void SetFiducialCutParameters(std::string beam_en); // Load Fidicual Parameters for 1.1 and 4.4 GeV from file
//void SetMomCorrParameters();

// Also used by FilterData.{C,h}

TF1 *pipl_deltat_sig,*pipl_deltat_mean,*pimi_deltat_sig,*pimi_deltat_mean, *prot_deltat_sig, *prot_deltat_mean,*el_Epratio_sig,*el_Epratio_mean;

double vz_corr(TF1 *vz_corr_func, double phi,double theta);
TVector3 FindUVW(TVector3 xyz);
Bool_t CutUVW(TVector3 ecxyz);
Float_t ProtonMomCorrection_He3_4Cell(std::string atarget, TLorentzVector V4Pr, Float_t vertex_p);


std::map<std::string,double>vert_min;
std::map<std::string,double>vert_max;
std::map<std::string,double>vertdiff_min;
std::map<std::string,double>vertdiff_max;
std::map<std::string,double>bind_en;
std::map<std::string,double>target_mass;
std::map<std::string,double>residual_target_mass;
std::map<std::pair<std::string, int>, double> EC_time_offset;
std::map<std::string,double>EC_photon_beta;
std::map<std::string,double>LEC_photon_beta;
std::map<std::string, double> Ecal_offset;



//e- E_tot/p vs p PID cut
Double_t FSum_e(Double_t *x,Double_t *par){
  if(x[0]<par[1])       return el_Epratio_mean->EvalPar(x)+par[0]*el_Epratio_sig->EvalPar(x);
  else if(x[0]>=par[1]) return el_Epratio_mean->Eval(par[1])+par[0]*el_Epratio_sig->Eval(par[1]);
  else return -1;
}

Double_t FSub_e(Double_t *x,Double_t *par){
  if(x[0]<par[1])       return el_Epratio_mean->EvalPar(x)-par[0]*el_Epratio_sig->EvalPar(x);
  else if(x[0]>=par[1]) return el_Epratio_mean->Eval(par[1])-par[0]*el_Epratio_sig->Eval(par[1]);
  else return -1;
}

//proton Delta_t vs momentum PID cut

Double_t FSum_prot(Double_t *x, Double_t *par){   //the 2 parameters are the cut range and momentum limit
  if(x[0] < par[1])         return prot_deltat_mean->EvalPar(x)+par[0]*prot_deltat_sig->EvalPar(x);
  else if(x[0] >= par[1])   return prot_deltat_mean->Eval(par[1])+par[0]*prot_deltat_sig->Eval(par[1]);
  else return -1;
}
Double_t FSub_prot(Double_t *x,Double_t *par){
  if(x[0] < par[1])         return prot_deltat_mean->EvalPar(x)-par[0]*prot_deltat_sig->EvalPar(x);
  else if(x[0] >= par[1])   return prot_deltat_mean->Eval(par[1])-par[0]*prot_deltat_sig->Eval(par[1]);
  else return -1;
}


//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions

Double_t FSum_pimi(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pimi_deltat_mean->EvalPar(x)+par[0]*pimi_deltat_sig->EvalPar(x);
  else if(x[0]>=par[1])return pimi_deltat_mean->Eval(par[1])+par[0]*pimi_deltat_sig->Eval(par[1]);
  else return -1;
}
Double_t FSub_pimi(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pimi_deltat_mean->EvalPar(x)-par[0]*pimi_deltat_sig->EvalPar(x);
  else if(x[0]>=par[1])return pimi_deltat_mean->Eval(par[1])-par[0]*pimi_deltat_sig->Eval(par[1]);
  else return -1;
}



//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions
Double_t FSum_pipl(Double_t *x,Double_t *par){

  if(x[0]<par[1])  return pipl_deltat_mean->EvalPar(x)+par[0]*pipl_deltat_sig->EvalPar(x);
  else if(x[0]>=par[1]) return pipl_deltat_mean->Eval(par[1])+par[0]*pipl_deltat_sig->Eval(par[1]);
  else return -1;
}
Double_t FSub_pipl(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pipl_deltat_mean->EvalPar(x)-par[0]*pipl_deltat_sig->EvalPar(x);
  else if(x[0]>=par[1])return pipl_deltat_mean->Eval(par[1])-par[0]*pipl_deltat_sig->Eval(par[1]);
  else return -1;
}


 const int N_tot=100;




void e2a_ep_neutrino6_united4_radphot::Loop()
{

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
  //      nentries =8000000;
  //     nentries =1000000;

  double N_prot1 = 0, N_prot2 = 0,N_prot_both = 0;
  double eps;
  Float_t pimi_phimin = 0, pimi_phimax = 0;
  Float_t pipl_phimin = 0, pipl_phimax = 0;
  const int n_slice=3,nsect=6;
  double beta,delta;
  const double pperp_min[n_slice]={0.,0.2,0.4};
  const double pperp_max[n_slice]={0.2,0.4,10.};
  TVector3 V3_pimi,V3_pipl,V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;
  TVector3 V3_phot_angles;
  double sum_val,sub_val;
  double epratio_sig_cutrange=3.;
  double prot_delt_cutrange=3.;
  int 	el_segment, el_cc_sector;
  double delt_uplim,delt_lowlim;
  double prot_accept_mom_lim;
  double prot_mom_lim;
  double min_good_mom;
  double max_mom;
  Double_t el_sccc_timediff;
  Double_t sc_cc_delt_cut_sect[nsect]={-2,-5,-8,-8,-2,2};
  Double_t el_cc_nphe;
  Double_t elmom_corr_fact[nsect];
  double pipl_maxmom, pimi_maxmom,pimi_delt_cutrange,pipl_delt_cutrange;
  int N_pperp,N_Ecal;
  double *pperp_cut,*Ecal_lowlim,*Ecal_uplim;
  TF1 *vz_corr_func;


  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV  Configuration parameters and cuts
  {
      eps=0.02;
      prot_accept_mom_lim=0.3;
      prot_mom_lim=0.95;
      min_good_mom=0.4;
      max_mom=1.1;
      pipl_maxmom=0.65;
      pimi_maxmom=0.6;
      pimi_delt_cutrange=3.;
      pipl_delt_cutrange=3.;
      N_pperp=2;
      N_Ecal=6;
      pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;
      pperp_cut[1]=0.2;
      for (int i=0;i<N_Ecal;i++){
         Ecal_lowlim[i]=0.45+i*0.18;
	       Ecal_uplim[i]=0.63+i*0.18;
      }
      Ecal_lowlim[5]=0.;
      Ecal_uplim[5]=1.35;
      for (int i=0;i<N_Ecal;i++)	cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<endl;

      vert_min["3He"]=-3.05;
      vert_min["C12"]=4.95;
      vert_min["CH2"]=4.85;
      vert_max["3He"]=-0.18;
      vert_max["C12"]=5.76;
      vert_max["CH2"]=5.62;

      vertdiff_min["3He"]=-1.;
      vertdiff_min["C12"]=-1.;
      vertdiff_min["CH2"]=-1.;

      vertdiff_max["3He"]=1.;
      vertdiff_max["C12"]=1.;
      vertdiff_max["CH2"]=1.;

      EC_photon_beta["3He"]=0.89;
      EC_photon_beta["C12"]=0.89;
      EC_photon_beta["CH2"]=0.91;

      LEC_photon_beta["3He"]=0.97;
      LEC_photon_beta["C12"]=0.97;
      LEC_photon_beta["CH2"]=0.97;

      EC_time_offset[std::make_pair("3He",1)]=-0.73;  EC_time_offset[std::make_pair("3He",2)]=-0.81; EC_time_offset[std::make_pair("3He",3)]=-0.91;
      EC_time_offset[std::make_pair("3He",4)]=-0.94;  EC_time_offset[std::make_pair("3He",5)]=-0.92; EC_time_offset[std::make_pair("3He",6)]=-0.81;

      EC_time_offset[std::make_pair("C12",1)]=-0.71;  EC_time_offset[std::make_pair("C12",2)]=-0.77; EC_time_offset[std::make_pair("C12",3)]=-0.87;
      EC_time_offset[std::make_pair("C12",4)]=-0.91;  EC_time_offset[std::make_pair("C12",5)]=-0.89; EC_time_offset[std::make_pair("C12",6)]=-0.79;

      EC_time_offset[std::make_pair("CH2",1)]=-0.70;  EC_time_offset[std::make_pair("CH2",2)]=-0.80; EC_time_offset[std::make_pair("CH2",3)]=-0.91;
      EC_time_offset[std::make_pair("CH2",4)]=-0.92;  EC_time_offset[std::make_pair("CH2",5)]=-0.91; EC_time_offset[std::make_pair("CH2",6)]=-0.80;

      elmom_corr_fact[0]=1.007;
      elmom_corr_fact[1]=0.988;
      elmom_corr_fact[2]=1.008;
      elmom_corr_fact[3]=1.011;
      elmom_corr_fact[4]=1.014;
      elmom_corr_fact[5]=1.013;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
  }


  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV  Configuration parameters and cuts
  {
      eps=0.02;//GeV
      prot_accept_mom_lim=0.3;
      prot_mom_lim=2.15;
      min_good_mom=0.55;
      max_mom=2.1;
      pipl_maxmom=1.4;
      pimi_maxmom=1.3;
      pimi_delt_cutrange=3.;
      pipl_delt_cutrange=3.;
      N_pperp=2;
      N_Ecal=6;
      pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;
      pperp_cut[1]=0.2;
      for (int i=0;i<N_Ecal;i++){
	       Ecal_lowlim[i]=0.75+i*0.25;
	        Ecal_uplim[i]=1.+i*0.25;
      }
      Ecal_lowlim[5]=0.;
      Ecal_uplim[5]=2.;
      for (int i=0;i<N_Ecal;i++)	cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<endl;

      vert_min["3He"]=-3.29;
      vert_min["4He"]=-2.53;
      vert_min["C12"]=4.8;
      vert_min["56Fe"]=4.6;
      vert_max["3He"]=-0.23;
      vert_max["4He"]=1.73;
      vert_max["C12"]=5.5;
      vert_max["56Fe"]=5.3;

      vertdiff_min["3He"]=-1.;
      vertdiff_min["4He"]=-1.;
      vertdiff_min["C12"]=-1.;
      vertdiff_min["56Fe"]=-1.;

      vertdiff_max["3He"]=1.;
      vertdiff_max["4He"]=1.;
      vertdiff_max["C12"]=1.;
      vertdiff_max["56Fe"]=1.;

      EC_photon_beta["3He"]=0.93;
      EC_photon_beta["4He"]=0.92;
      EC_photon_beta["C12"]=0.92;
      EC_photon_beta["56Fe"]=0.90;

      LEC_photon_beta["3He"]=0.96;
      LEC_photon_beta["4He"]=0.94;
      LEC_photon_beta["C12"]=0.94;
      LEC_photon_beta["56Fe"]=0.95;

      EC_time_offset[std::make_pair("3He",1)]=-1.37;  EC_time_offset[std::make_pair("3He",2)]=-1.42; EC_time_offset[std::make_pair("3He",3)]=-1.55;
      EC_time_offset[std::make_pair("3He",4)]=-1.53;  EC_time_offset[std::make_pair("3He",5)]=-1.49; EC_time_offset[std::make_pair("3He",6)]=-1.44;

      EC_time_offset[std::make_pair("4He",1)]=0.72;  EC_time_offset[std::make_pair("4He",2)]=0.27; EC_time_offset[std::make_pair("4He",3)]=0.16;
      EC_time_offset[std::make_pair("4He",4)]=0.21;  EC_time_offset[std::make_pair("4He",5)]=0.22; EC_time_offset[std::make_pair("4He",6)]=0.21;

      EC_time_offset[std::make_pair("C12",1)]=0.50;  EC_time_offset[std::make_pair("C12",2)]=0.39; EC_time_offset[std::make_pair("C12",3)]=0.29;
      EC_time_offset[std::make_pair("C12",4)]=0.29;  EC_time_offset[std::make_pair("C12",5)]=0.32; EC_time_offset[std::make_pair("C12",6)]=0.33;

      EC_time_offset[std::make_pair("56Fe",1)]=0.75;  EC_time_offset[std::make_pair("56Fe",2)]=0.49; EC_time_offset[std::make_pair("56Fe",3)]=0.37;
      EC_time_offset[std::make_pair("56Fe",4)]=0.39;  EC_time_offset[std::make_pair("56Fe",5)]=0.43; EC_time_offset[std::make_pair("56Fe",6)]=0.44;

      elmom_corr_fact[0]=1.001;
      elmom_corr_fact[1]=0.991;
      elmom_corr_fact[2]=1.005;
      elmom_corr_fact[3]=1.004;
      elmom_corr_fact[4]=1.006;
      elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
  }

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5)  //4.4 GeV  Configuration parameters and cuts
  {
      eps=0.02;//GeV
      prot_accept_mom_lim=0.3;
      prot_mom_lim=2.7;
      min_good_mom=1.1;
      if(ftarget=="3He") min_good_mom=1.3;//the EC threshold was different for He3
      max_mom=3.7;
      pipl_maxmom=1.9;
      pipl_delt_cutrange=3.;
      pimi_maxmom=1.6;
      pimi_delt_cutrange=3.;
      N_pperp=2;
      N_Ecal=6;
      pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;
      pperp_cut[1]=0.2;
      for (int i=0;i<N_Ecal;i++){
	        Ecal_lowlim[i]=1.5+i*0.5;
	        Ecal_uplim[i]=2.+i*0.5;
      }
      Ecal_lowlim[5]=0.;
      Ecal_uplim[5]=4.;
      for (int i=0;i<N_Ecal;i++)	cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<endl;;

      vert_min["3He"]=-3.27;
      vert_min["4He"]=-2.51;
      vert_min["C12"]=4.7;
      vert_min["56Fe"]=4.6;
      vert_max["3He"]=0.07;
      vert_max["4He"]=1.71;
      vert_max["C12"]=5.3;
      vert_max["56Fe"]=5.4;

      vertdiff_min["3He"]=-1.;
      vertdiff_min["4He"]=-1;
      vertdiff_min["C12"]=-1;
      vertdiff_min["56Fe"]=-1;

      vertdiff_max["3He"]=1.;
      vertdiff_max["4He"]=1.;
      vertdiff_max["C12"]=1;
      vertdiff_max["56Fe"]=1;

      EC_photon_beta["3He"]=0.92;
      EC_photon_beta["4He"]=0.91;
      EC_photon_beta["C12"]=0.92;
      EC_photon_beta["56Fe"]=0.91;

      LEC_photon_beta["3He"]=0.97;
      LEC_photon_beta["4He"]=0.97;
      LEC_photon_beta["C12"]=0.95;
      LEC_photon_beta["56Fe"]=0.96;

      EC_time_offset[std::make_pair("3He",1)]=-0.15;  EC_time_offset[std::make_pair("3He",2)]=-0.26; EC_time_offset[std::make_pair("3He",3)]=-0.41;
      EC_time_offset[std::make_pair("3He",4)]=-0.29;  EC_time_offset[std::make_pair("3He",5)]=-0.25; EC_time_offset[std::make_pair("3He",6)]=-0.23;

      EC_time_offset[std::make_pair("4He",1)]=-0.01;  EC_time_offset[std::make_pair("4He",2)]=-0.11; EC_time_offset[std::make_pair("4He",3)]=-0.23;
      EC_time_offset[std::make_pair("4He",4)]=-0.26;  EC_time_offset[std::make_pair("4He",5)]=-0.21; EC_time_offset[std::make_pair("4He",6)]=-0.09;

      EC_time_offset[std::make_pair("C12",1)]=-0.01;  EC_time_offset[std::make_pair("C12",2)]=-0.11; EC_time_offset[std::make_pair("C12",3)]=-0.23;
      EC_time_offset[std::make_pair("C12",4)]=-0.27;  EC_time_offset[std::make_pair("C12",5)]=-0.21; EC_time_offset[std::make_pair("C12",6)]=-0.08;

      EC_time_offset[std::make_pair("56Fe",1)]=-0.49;  EC_time_offset[std::make_pair("56Fe",2)]=-0.14; EC_time_offset[std::make_pair("56Fe",3)]=-0.32;
      EC_time_offset[std::make_pair("56Fe",4)]=-0.25;  EC_time_offset[std::make_pair("56Fe",5)]=-0.17; EC_time_offset[std::make_pair("56Fe",6)]=-0.35;

      elmom_corr_fact[0]=1.001;
      elmom_corr_fact[1]=0.991;
      elmom_corr_fact[2]=1.005;
      elmom_corr_fact[3]=1.004;
      elmom_corr_fact[4]=1.006;
      elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
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
  TH1F *h1_el_cc_deltat[nsect], *h1_el_cc_deltat_cut[nsect],*h1_Erec_p_bkgd_slice[n_slice],*h1_Etot_p_bkgd_slice[n_slice],*h1_Etot_Npi0[n_slice],*h1_Erec_Npi0_new[n_slice],* h1_Erec_bkgd_pipl_pimi_new_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pipl[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pimi[n_slice],*h1_Etot_Npi1[n_slice],*h1_Erec_Npi1[n_slice],*h1_Etot_bkgd_1p2pi[n_slice],*h1_Erec_bkgd_1p2pi[n_slice],*h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[n_slice],*h1_Etot_p_bkgd_slice_2p2pi[n_slice],*h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p2pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[n_slice],*h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[n_slice],*h1_Etot_bkgd_1p2pi_1p0pi[n_slice],*h1_Erec_bkgd_1p2pi_1p0pi[n_slice],*h1_Etot_bkgd_1p3pi[n_slice],*h1_Erec_bkgd_1p3pi[n_slice],*h1_e_mom_corrfuct[nsect];
  TH1F *h1_Etot_piplpimi_subtruct_fact[n_slice],*h1_Erec_piplpimi_subtruct_new_fact[n_slice],*h1_Etot_p_bkgd_slice_sub[n_slice],*h1_Erec_p_bkgd_slice_sub[n_slice],*h1_Etot_3pto1p_slice[n_slice],*h1_Erec_3pto1p_slice[n_slice],*h1_Etot_3pto2p_slice[n_slice],*h1_Erec_3pto2p_slice[n_slice],*h1_Etot_3p1pi_slice[n_slice],*h1_Erec_3p1pi_slice[n_slice],*h1_Etot_4pto3p_slice[n_slice],*h1_Erec_4pto3p_slice[n_slice],*h1_Etot_4pto1p_slice[n_slice],*h1_Erec_4pto1p_slice[n_slice],*h1_Etot_4pto2p_slice[n_slice],*h1_Erec_4pto2p_slice[n_slice],*h1_Etot_43pto1p_slice[n_slice],*h1_Erec_43pto1p_slice[n_slice],*h1_Etot_p_bkgd_slice_sub1p2pi[n_slice], *h1_Erec_p_bkgd_slice_sub1p2pi[n_slice],*h1_Etot_p_bkgd_slice_sub2p1pi_1p[n_slice],*h1_Etot_p_bkgd_slice_sub3p1pi_0p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p[n_slice],*h1_Erec_p_bkgd_slice_sub3p1pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub2p1pi_2p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_2p[n_slice],*h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[n_slice],*h1_Etot_p_bkgd_slice_sub1p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p2pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub1p3pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub3p1pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p3pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub2p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p2pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub32[n_slice],*h1_Erec_p_bkgd_slice_sub32[n_slice],*h1_Etot_p_bkgd_slice_sub31[n_slice],*h1_Erec_p_bkgd_slice_sub31[n_slice],*h1_Etot_p_bkgd_slice_sub43[n_slice],*h1_Erec_p_bkgd_slice_sub43[n_slice],*h1_Etot_p_bkgd_slice_sub41[n_slice],*h1_Erec_p_bkgd_slice_sub41[n_slice],*h1_Erec_p_bkgd_slice_sub42[n_slice],*h1_Etot_p_bkgd_slice_sub42[n_slice],*h1_Etot_p_bkgd_slice_sub431[n_slice],*h1_Erec_p_bkgd_slice_sub431[n_slice];
  TH1F *h1_el_ec_sc_timediff_sect_corr[nsect], *h1_el_ec_sc_timediff_sect[nsect],*h1_beta_ec_corr_sect[nsect],*h1_el_SCpdfidcut[6],*h1_el_SCpd[6];
  TH2F *h2_el_theta_phi_p_beffidcut[6],*h2_el_theta_phi_p_fidcut[6],*h2_el_ec_sc_timediff_ecu[6],*h2_el_ec_sc_timediff_ecv[6],*h2_el_ec_sc_timediff_ecw[6], *h2_el_ec_sc_timediff_SCpd[6],*h2_prot_theta_phi_p_beffidcut[6],*h2_prot_theta_phi_p_fidcut[6],*h2_el_theta_phi_p_beffidcut2[6],*h2_el_theta_phi_p_fidcut2[6],*h2_N_pi_phot[20],*h2_pimi_theta_phi_p_beffidcut[6],*h2_pimi_theta_phi_p_fidcut[6],*h2_pipl_theta_phi_p_beffidcut[6],*h2_pipl_theta_phi_p_fidcut[6],*h2_prot_theta_p[6],*h2_prot_theta_p_cut[6],*h2_pimi_theta_p[6],*h2_pimi_theta_p_cut[6],*h2_pipl_theta_p[6],*h2_pipl_theta_p_cut[6];


  gRandom = new TRandom3();
  gRandom->SetSeed(10);

  TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
  TLorentzVector V4_target(0,0,0,target_mass[ftarget]);

  //Definition of other input files with calibration data
  TFile *file_in1=new TFile(Form("FiducialsCorrections/protdeltat_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in=new TFile(Form("FiducialsCorrections/el_Epratio_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in2=new TFile(Form("FiducialsCorrections/pimideltat_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in3= new TFile(Form("FiducialsCorrections/vz_%s_%s.root",ftarget.c_str(),fbeam_en.c_str()));;
  TFile *file_in4=new TFile(Form("FiducialsCorrections/pipldeltat_mom_%s.root",fbeam_en.c_str()));


  TFile *file_in5;
  TFile *file_in6;
  TFile *file_in7;

  if (en_beam[fbeam_en] < 2.300 && en_beam[fbeam_en] > 2.100 ) {
    file_in5 = new TFile("FiducialsCorrections/vz_56Fe_2261_badruns.root");//vertex correction for 56Fe runs with exploded liquid target cell
    file_in6 = new TFile("FiducialsCorrections/vz_3He_2261_2ndrungroup.root");//vertx correction for 3He 2nd group runs
    file_in7 = new TFile("FiducialsCorrections/vz_4He_2261_2ndrungroup.root");//vertex correction for 4He 2nd group runs
  }

  //Output file definition
  TFile *file_out = new TFile(Form("e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");

  double pars[3];
  if (en_beam[fbeam_en] < 2.300 && en_beam[fbeam_en] > 2.100 ) {
    if(ftarget=="56Fe"){
      pars[0]=((TF1 *)file_in5->Get("f_vz"))->GetParameter(0);
      pars[1]=((TF1 *)file_in5->Get("f_vz"))->GetParameter(1);
      pars[2]=((TF1 *)file_in5->Get("f_vz"))->GetParameter(2);
    }
    if(ftarget=="3He"){
      pars[0]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(0);
      pars[1]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(1);
      pars[2]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(2);
    }
    if(ftarget=="4He"){
     pars[0]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(0);
     pars[1]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(1);
     pars[2]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(2);
   }
  }
  else {
    pars[0] = 0;
    pars[1] = 0;
    pars[2] = 0;
  }

 //Reading of input functions for calibrations
  vz_corr_func = (TF1 *)file_in3->Get("f_vz");
  el_Epratio_mean=(TF1*)file_in->Get("f_mean");
  el_Epratio_sig=(TF1*)file_in->Get("f_sig");
  prot_deltat_sig=(TF1*)file_in1->Get("sig_pol9");
  prot_deltat_mean=(TF1*)file_in1->Get("mean_pol9");
  pipl_deltat_sig=(TF1*)file_in4->Get("sig_pol9");
  pipl_deltat_mean=(TF1*)file_in4->Get("mean_pol9");
  pimi_deltat_sig=(TF1*)file_in2->Get("sig_pol9");
  pimi_deltat_mean=(TF1*)file_in2->Get("mean_pol9");

  //More functions
  TF1 *fsum_e,*fsub_e; // electron TF1 functions for E/p
  TF1 *fsum_prot,*fsub_prot; // proton TF1 functions for E/p
  TF1 *fsum_pimi,*fsub_pimi; // Pi minus TF1 functions for E/p
  TF1 *fsum_pipl,*fsub_pipl; // Pi Plus TF1 functions for E/p

  fsum_pimi=new TF1("fsum_pimi",FSum_pimi,0.,5.,2);
  fsub_pimi=new TF1("fsub_pimi",FSub_pimi,0.,5.,2);
  fsum_pipl=new TF1("fsum_pipl",FSum_pipl,0.,5.,2);
  fsub_pipl=new TF1("fsub_pipl",FSub_pipl,0.,5.,2);
  fsum_prot=new TF1("fsum_prot",FSum_prot,0.,5.,2);
  fsub_prot=new TF1("fsub_pprot",FSub_prot,0.,5.,2);
  fsum_e = new TF1("fsum_e",FSum_e,0.,5.,2);
  fsub_e = new TF1("fsub_e",FSub_e,0.,5.,2);


  //initialize Fiducial functions for EC limits
  fiducialcut->InitEClimits();

  //Definition and initialization of Histograms
  TH1F *h1_el_vertuncorr=new TH1F("h1_el_vertuncorr","",200,-10,10);
  TH1F *h1_el_vertcorr=new TH1F("h1_el_vertcorr","",200,-10,10);
  TH1F *h1_el_Mott_crosssec = new TH1F("h1_el_Mott_crosssec","",200,0.,0.01);
  TH1F *h1_el_Etot = new TH1F("h1_el_Etot","",600,0,4);
  TH1F *h1_el_Ein = new TH1F("h1_el_Ein","",800,0,1.5);
  TH1F *h1_el_Etot_cut = new TH1F("h1_el_Etot_cut","",600,0,4);
  TH1F *h1_el_Ein_cut = new TH1F("h1_el_Ein_cut","",800,0,1.5);
  TH1F *h1_el_cc_nphe = new TH1F("h1_el_cc_nphe","",300,0,30);
  TH1F *h1_el_cc_nphe_cut = new TH1F("h1_el_cc_nphe_cut","",300,0,30);
  TH1F *h1_el_cc_nphe_cut2 = new TH1F("h1_el_cc_nphe_cut2","",300,0,30);
  TH1F *h1_el_cc_chi2 = new TH1F("h1_el_cc_chi2","",200,0,0.7);
  TH1F *h1_Wvar = new TH1F("h1_Wvar","",400,0,3);
  TH1F *h1_Wvar_weight = new TH1F("h1_Wvar_weight","",400,0,3);
  TH1F *h1_Wepp = new TH1F("h1_Wepp","",400,0,3);
  TH1F *h1_Wepp_uncorr = new TH1F("h1_Wepp_uncorr","",400,0,3);
  TH1F *h1_xbjk = new TH1F("h1_xbjk","",400,0,3);
  TH1F *h1_xbjk_weight = new TH1F("h1_xbjk_weight","",400,0,3);
  TH1F *h1_Q2 = new TH1F("h1_Q2","",400,0,6);
  TH1F *h1_Q2_weight = new TH1F("h1_Q2_weight","",400,0,6);
  TH1F *h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
  TH1F *h1_Nprot=new TH1F("h1_Nprot","",10,0,5);
  TH1F *h1_Nprot_NonZeroProt=new TH1F("h1_Nprot_NonZeroProt","",8,1,5);
  TH1F *h1_Nphot=new TH1F("h1_Nphot","",10,0,5);
  TH1F *h1_Npiphot=new TH1F("h1_Npiphot","",10,0,5);
  TH1F *h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,0,5);
  TH1F *h1_photon_E=new TH1F("h1_photon_E","",200,0,4.5);
  TH1F *h1_photon_EC_E=new TH1F("h1_photon_EC_E","",200,0,4.5);
  TH1F *h1_phot_e_angle= new TH1F("h1_phot_e_angle","",300,0,180);
  TH1F *h1_time_ec=new TH1F("h1_time_ec","",200,-20,20);
  TH1F *h1_Npi=new TH1F("h1_Npi","",10,0,5);
  TH1F *h1_Npi_NonZeroProt=new TH1F("h1_Npi_NonZeroProt","",10,0,5);
  TH1F *h1_Npipl=new TH1F("h1_Npipl","",10,0,5);
  TH1F *h1_Npimi=new TH1F("h1_Npimi","",10,0,5);
  TH1F *h1_el_mom = new TH1F("h1_el_mom","",100,1.2,6);
  TH1F *h1_el_mom_corr = new TH1F("h1_el_mom_corr","",100,1.2,6);
  TH1F *h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
  TH1F *h1_el_prot_vertdiff_all= new TH1F("h1_el_prot_vertdiff_all","",300,-10,10);
  TH1F *h1_prot_vertdiff= new TH1F("h1_prot_vertdiff","",300,-10,10);
  TH1F *h1_el_prot_vertdiff= new TH1F("h1_el_prot_vertdiff","",300,-10,10);
  TH1F *h1_el_prot_vertdiff1= new TH1F("h1_el_prot_vertdiff1","",300,-10,10);
  TH1F *h1_el_prot_vertdiff2= new TH1F("h1_el_prot_vertdiff2","",300,-10,10);
  TH1F *h1_el_3prot_vertdiff1= new TH1F("h1_el_3prot_vertdiff1","",300,-10,10);
  TH1F *h1_el_3prot_vertdiff2= new TH1F("h1_el_3prot_vertdiff2","",300,-10,10);
  TH1F *h1_el_3prot_vertdiff3= new TH1F("h1_el_3prot_vertdiff3","",300,-10,10);
  TH1F *h1_el_4prot_vertdiff1= new TH1F("h1_el_4prot_vertdiff1","",300,-10,10);
  TH1F *h1_el_4prot_vertdiff2= new TH1F("h1_el_4prot_vertdiff2","",300,-10,10);
  TH1F *h1_el_4prot_vertdiff3= new TH1F("h1_el_4prot_vertdiff3","",300,-10,10);
  TH1F *h1_el_4prot_vertdiff4= new TH1F("h1_el_4prot_vertdiff4","",300,-10,10);
  TH1F *h1_pipl_prot_vertdiff= new TH1F("h1_pipl_prot_vertdiff","",300,-10,10);
  TH1F *h1_pimi_prot_vertdiff= new TH1F("h1_pimi_prot_vertdiff","",300,-10,10);
  TH1F *h1_prot_mom = new TH1F("h1_prot_mom","",300,0,3);
  TH1F *h1_prot_mom_ratio = new TH1F("h1_prot_mom_ratio","",50,0.97,1.2);
  TH1F *h1_neg_m=new TH1F("h1_neg_m","",300,0,2);
  TH1F *h1_pos_m=new TH1F("h1_pos_m","",300,0,2);


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

  TH1F *h1_beta_ec = new TH1F("h1_beta_ec","",300,0,2);
  TH1F *h1_beta_ec_corr = new TH1F("h1_beta_ec_corr","",300,0,2);
  TH1F *h1_beta_ec_corr_cut = new TH1F("h1_beta_ec_corr_cut","",300,0,2);
  TH1F *h1_el_ec_sc_timediff = new TH1F("h1_el_ec_sc_timediff","",600,-30,30);
  TH1F *h1_el_ec_sc_timediff_corr = new TH1F("h1_el_ec_sc_timediff_corr","",600,-30,30);
  TH1F *h1_el_ec_sc_timediff_allSCpd = new TH1F("h1_el_ec_sc_timediff_allSCpd","",600,-30,30);
  TH1F *h1_el_ec_sc_timediff_corr_allSCpd = new TH1F("h1_el_ec_sc_timediff_corr_allSCpd","",600,-30,30);
  TH1F *h1_theta0=new TH1F("h1_theta0","",300,0,180);
  TH2F *h2_Ecal_Eqe=new TH2F("h2_Ecal_Eqe","",600,0,5,600,0,5);
  TH2F *h2_EqeEcalratio_Eqe=new TH2F("h2_EqeEcalratio_Eqe","",600,0,5,300,0,2);
  TH2F *h2_EqeEcaldiff_Eqe=new TH2F("h2_EqeEcaldiff_Eqe","",600,0,5,300,-3,3);
  TH2F *h2_N_prot_pi=new TH2F("h2_N_prot_pi","",10,0,5,10,0,5);
  TH2F *h2_N_prot_pi_phot=new TH2F("h2_N_prot_pi_phot","",10,0,5,10,0,5);
  TH2F *h2_N_prot_pi_phot_nonrad=new TH2F("h2_N_prot_pi_phot_nonrad","",10,0,5,10,0,5);
  TH2F *h2_el_E_p_ratio_cut = new TH2F("h2_el_E_p_ratio_cut","",200,0,4.5,200,0,0.5);
  TH2F *h2_el_E_p_ratio = new TH2F("h2_el_E_p_ratio","",200,0,4.5,200,0,0.5);
  TH2F *h2_el_E_p_ratio_withoutCC = new TH2F("h2_el_E_p_ratio_withoutCC","",200,0,4.5,200,0,0.5);
  TH2F *h2_el_E_p_ratio_withCC = new TH2F("h2_el_E_p_ratio_withCC","",200,0,4.5,200,0,0.5);
  TH2D *h2_el_Ein_Eout=new TH2D("h2_el_Ein_Eout","",400,0,2,400,0,1.5);
  TH2D *h2_el_Einout_Etot=new TH2D("h2_el_Einout_Etot","",400,0,4,400,0,4);
  TH2F *h2_el_ec_xy = new TH2F("h2_el_ec_xy","",100,-600,600,100,-600,600);
  TH2F *h2_el_ec_xy_fidcut = new TH2F("h2_el_ec_xy_fidcut","",100,-600,600,100,-600,600);
  TH2F *h2_el_phi_vert = new TH2F("h2_el_phi_vert","",500,-6,7,120,0,380);
  TH2F *h2_el_phi_vert_uncorr = new TH2F("h2_el_phi_vert_uncorr","",500,-6,7,120,0,380);
  TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,0,180);
  TH2F *h2_neutral_costheta_phi_EC_all = new TH2F("h2_neutral_costheta_phi_EC_all","",200,0,360,200,0,1.1);
  TH2F *h2_neutral_theta_phi_EC_all = new TH2F("h2_neutral_theta_phi_EC_all","",200,0,360,200,0,180);
  TH2F *h2_neutral_theta_phi_EC_all_fidcut = new TH2F("h2_neutral_theta_phi_EC_all_fidcut","",200,0,360,200,0,180);
  TH2F *h2_pimi_theta_phi = new TH2F("h2_pimi_theta_phi","",200,0,360,200,0,180);
  TH2F *h2_pipl_theta_phi = new TH2F("h2_pipl_theta_phi","",200,0,360,200,0,180);
  TH2F *h2_pimi_theta_phi_beffid = new TH2F("h2_pimi_theta_phi_beffid","",200,0,360,200,0,180);
  TH2F *h2_pipl_theta_phi_beffid = new TH2F("h2_pipl_theta_phi_beffid","",200,0,360,200,0,180);
  TH2F *h2_prot_theta_phi = new TH2F("h2_prot_theta_phi","",200,0,360,200,0,180);
  TH2F *h2_prot_px_py_p = new TH2F("h2_prot_px_py_p","",100,-1,1,100,-1,1);
  TH2F *h2_prot_px_py_p_fidcut = new TH2F("h2_prot_px_py_p_fidcut","",100,-1,1,100,-1,1);
  TH2F *h2_pipl_theta_phi_p = new TH2F("h2_pipl_theta_phi_p","",200,0,360,200,0,180);
  TH2F *h2_pipl_theta_phi_fidcut = new TH2F("h2_pipl_theta_phi_fidcut","",200,0,360,200,0,180);
  TH2F *h2_pimi_theta_phi_p = new TH2F("h2_pimi_theta_phi_p","",200,0,360,200,0,180);
  TH2F *h2_pimi_theta_phi_fidcut = new TH2F("h2_pimi_theta_phi_fidcut","",200,0,360,200,0,180);
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

  TH2F *h2_pipl_delt_p= new TH2F("h2_pipl_delt_p","",300,0.,2.5,300,-15,15);
  TH2F *h2_pimi_delt_p= new TH2F("h2_pimi_delt_p","",300,0.,2.5,300,-15,15);
  TH2F *h2_pos_delt_p= new TH2F("h2_pos_delt_p","",300,0.,2.5,300,-15,15);
  TH2F *h2_neg_delt_p= new TH2F("h2_neg_delt_p","",300,0.,2.5,300,-15,15);
  TH2F *h2_pipl_beta_p = new TH2F("h2_pipl_beta_p","",200,0,4,200,0,1.2);
  TH2F *h2_pimi_beta_p = new TH2F("h2_pimi_beta_p","",200,0,4,200,0,1.2);
  TH2F *h2_neg_beta_p = new TH2F("h2_neg_beta_p","",200,0,4,200,0,1.2);
  TH2F *h2_pos_beta_p = new TH2F("h2_pos_beta_p","",200,0,4,200,0,1.2);
  TH2F *h2_prot_beta_p = new TH2F("h2_prot_beta_p","",200,0,4,200,0,1.2);
  TH2F *h2_neg_E_p = new TH2F("h2_neg_E_p","",200,0,3,200,0,200);
  TH2F *h2_pos_E_p = new TH2F("h2_pos_E_p","",200,0,3,200,0,200);
  TH2F *h2_pimi_E_p = new TH2F("h2_pimi_E_p","",200,0,3,200,0,200);
  TH2F *h2_pipl_E_p = new TH2F("h2_pipl_E_p","",200,0,3,200,0,200);
  TH2F *h2_prot_E_p = new TH2F("h2_prot_E_p","",200,0,3,200,0,200);
  TH2F *h2_prot_Deltat_p= new TH2F("h2_prot_Deltat_p","",300,0.,2.7,300,-15,15);
  TH2F *h2_el_vertcorr_runN= new TH2F("h2_el_vertcorr_runN","",30,18370,18436,300,-7,7);
  TH2F *h2_phot_e_angle_vsphotE= new TH2F("h2_phot_e_angle_vsphotE","",300,0,3,300,0,180);
  TH2F *h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);
  TH2F *h2_Wepp_ephi= new TH2F("h2_Wepp_ephi","",720,0,360,200,0.85,1.05);
  TH2F *h2_Wepp_ephi_corr= new TH2F("h2_Wepp_ephi_corr","",720,0,360,200,0.85,1.05);
  TH2F *h2_Wepp_ephi_uncorrprot= new TH2F("h2_Wepp_ephi_uncorrprot","",720,0,360,200,0.85,1.05);
  TH2F *h2_Wepp_ephi_corr_uncorrprot= new TH2F("h2_Wepp_ephi_corr_uncorrprot","",720,0,360,200,0.85,1.05);

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
  TH1F *h1_Etot = new TH1F("h1_Etot","",n_bins,x_values);
  TH1F *h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
  TH1F *h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
  TH1F *h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
	TH1F *h1_E_rec_undetfactor  = new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
  TH1F *h1_E_tot_undetfactor  = new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
  TH1F *h1_E_tot_undetfactor_pipl  = new TH1F("h1_E_tot_undetfactor_pipl","",n_bins,x_values);
  TH1F *h1_E_tot_undetfactor_pimi  = new TH1F("h1_E_tot_undetfactor_pimi","",n_bins,x_values);
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

  //Defintions of Histogram for each sector
  for(int h=0;h<nsect;h++) {
   h1_el_SCpdfidcut[h] = new TH1F(Form("h1_el_SCpdfidcut_%d",h+1),"",31,-0.5,30.5);
   h1_el_SCpd[h] = new TH1F(Form("h1_el_SCpd_%d",h+1),"",31,-0.5,30.5);
   h1_el_cc_deltat[h] = new TH1F(Form("h1_el_cc_deltat_%d",h+1),"",500,-100,100);
   h1_el_cc_deltat_cut[h] = new TH1F(Form("h1_el_cc_deltat_cut_%d",h+1),"",500,-100,100);
   h1_el_ec_sc_timediff_sect[h] = new TH1F(Form("h1_el_ec_sc_timediff_sect_%d",h+1),"",600,-30,30);
   h1_el_ec_sc_timediff_sect_corr[h] = new TH1F(Form("h1_el_ec_sc_timediff_sect_corr_%d",h+1),"",600,-30,30);
   h2_el_ec_sc_timediff_ecu[h] = new TH2F(Form("h2_el_ec_sc_timediff_ecu_%d",h+1),"",800,0,400,600,-30,30);
   h2_el_ec_sc_timediff_ecv[h] = new TH2F(Form("h2_el_ec_sc_timediff_ecv_%d",h+1),"",800,0,400,600,-30,30);
   h2_el_ec_sc_timediff_ecw[h] = new TH2F(Form("h2_el_ec_sc_timediff_ecw_%d",h+1),"",800,0,400,600,-30,30);
   h2_el_ec_sc_timediff_SCpd[h] = new TH2F(Form("h2_el_ec_sc_timediff_SCpd_%d",h+1),"",75,0,25,600,-30,30);
   h1_beta_ec_corr_sect[h] = new TH1F(Form("h1_beta_ec_corr_sect_%d",h+1),"",300,0,1.3);
   h1_e_mom_corrfuct[h] = new TH1F(Form("h1_e_mom_corrfuct_%d",h+1),"",400,0.7,1.3);
   h2_el_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_el_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_el_theta_phi_p_fidcut[h] = new TH2F(Form("h2_el_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_el_theta_phi_p_beffidcut2[h] = new TH2F(Form("h2_el_theta_phi_p_beffidcut2_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_el_theta_phi_p_fidcut2[h] = new TH2F(Form("h2_el_theta_phi_p_fidcut2_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_prot_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_prot_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_prot_theta_phi_p_fidcut[h] = new TH2F(Form("h2_prot_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pipl_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_pipl_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pipl_theta_phi_p_fidcut[h] = new TH2F(Form("h2_pipl_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pimi_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_pimi_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pimi_theta_phi_p_fidcut[h] = new TH2F(Form("h2_pimi_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_prot_theta_p[h] = new TH2F(Form("h2_prot_theta_p_%d",h+1),"",400,0,4,240,0,120);
   h2_prot_theta_p_cut[h] = new TH2F(Form("h2_prot_theta_p_cut_%d",h+1),"",400,0,4,240,0,120);
   h2_pimi_theta_p[h] = new TH2F(Form("h2_pimi_theta_p_%d",h+1),"",800,0,4,480,0,120);
   h2_pimi_theta_p_cut[h] = new TH2F(Form("h2_pimi_theta_p_cut_%d",h+1),"",400,0,4,240,0,120);
   h2_pipl_theta_p[h] = new TH2F(Form("h2_pipl_theta_p_%d",h+1),"",800,0,4,480,0,120);
   h2_pipl_theta_p_cut[h] = new TH2F(Form("h2_pipl_theta_p_cut_%d",h+1),"",400,0,4,240,0,120);
  }


  for(int h=0;h<n_slice;h++){
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
   h1_Erec_Npi0_new[h] = new TH1F(Form("h1_Erec_Npi0_new_%d",h+1),"",n_bins,x_values);
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
  //Decleration of 2-dim histograms
  TH1F *h1_Etot_Npi0_Ecalcut[2][6], *h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[2][6],*h1_Etot_p_bkgd_slice_Ecalcut[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_3p1pi_1p0pi[2][6];
  TH1F *h1_Etot_piplpimi_subtruct_fact_Ecalcut[2][6], *h1_Etot_p_bkgd_slice_sub_Ecalcut[2][6];
  TH1F *h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio3p[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio4p[2][6], *h1_Etot_p_bkgd_slice_Ecalcut41[2][6], *h1_Etot_p_bkgd_slice_Ecalcut421[2][6], *h1_Etot_p_bkgd_slice_Ecalcut431[2][6], *h1_Etot_p_bkgd_slice_Ecalcut4321[2][6], *h1_Etot_p_bkgd_slice_Ecalcut31[2][6],*h1_Etot_p_bkgd_slice_Ecalcut321[2][6], *h1_Etot_p_bkgd_slice_sub_Ecalcut41[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut42[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut431[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_1p3pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut43[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut31[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut32[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_2p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p2pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[2][6];
  TH1F *h1_E_tot_undetfactor09  = new TH1F("h1_E_tot_undetfactor09","",1,0,6);
  TH1F *h1_E_tot_cut2_09 = new TH1F("h1_E_tot_cut2_09","",1,0,6);
  TH1F *h1_E_tot_p_bkgd09  = new TH1F("h1_E_tot_p_bkgd09","",1,0,6);
  TH1F *h1_Etot_p321_bkgd09  = new TH1F("h1_Etot_p321_bkgd09","",1,0,6);
  TH1F *h1_Etot_p31_bkgd09  = new TH1F("h1_Etot_p31_bkgd09","",1,0,6);
  TH1F *h1_Etot_p4321_bkgd09  = new TH1F("h1_Etot_p4321_bkgd09","",1,0,6);
  TH1F *h1_Etot_p431_bkgd09  = new TH1F("h1_Etot_p431_bkgd09","",1,0,6);
  TH1F *h1_Etot_p421_bkgd09  = new TH1F("h1_Etot_p421_bkgd09","",1,0,6);
  TH1F *h1_Etot_p41_bkgd09  = new TH1F("h1_Etot_p41_bkgd09","",1,0,6);
  TH1F *h1_Etot_bkgd09_2p1pi_2p0pi  = new TH1F("h1_Etot_bkgd09_2p1pi_2p0pi","",1,0,6);
  TH1F *h1_Etot_bkgd09_2p1pi_1p1pi  = new TH1F("h1_Etot_bkgd09_2p1pi_1p1pi","",1,0,6);
  TH1F *h1_Etot_bkgd09_2p1pi_1p0pi  = new TH1F("h1_Etot_bkgd09_2p1pi_1p0pi","",1,0,6);
  TH1F *h1_Etot_bkgd09_1p2pi_1p0pi  = new TH1F("h1_Etot_bkgd09_1p2pi_1p0pi","",1,0,6);
  TH1F *h1_Etot_bkgd09_1p2pi_1p1pi  = new TH1F("h1_Etot_bkgd09_1p2pi_1p1pi","",1,0,6);
  TH1F *h1_Etot_bkgd09_1p3pi  = new TH1F("h1_Etot_bkgd09_1p3pi","",1,0,6);
  TH1F *h1_Etot_bkgd09_2p2pi  = new TH1F("h1_Etot_bkgd09_2p2pi","",1,0,6);
  TH1F *h1_Etot_bkgd09_3p1pi  = new TH1F("h1_Etot_bkgd09_3p1pi","",1,0,6);


  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){

      h1_Etot_Npi0_Ecalcut[i][j]=  new TH1F(Form("h1_Etot_Npi0_Ecalcut_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j]=  new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_Ecalcut_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut41[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut41_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut421[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut421_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut431[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut431_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut4321[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut4321_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut31[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut31_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut321[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut321_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut3_2p1pito2p0pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_1p3pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_2p2pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_3p1pi_%d_%d",i+1,j+1),"",1,0,6);
    }
  }

  //initialization of correct error arrays for the histograms
  h1_E_tot_cut2_09->Sumw2();
  h1_E_tot_undetfactor09->Sumw2();
  h1_E_tot_p_bkgd09->Sumw2();
  h1_Etot_p321_bkgd09->Sumw2();
  h1_Etot_p31_bkgd09->Sumw2();
  h1_Etot_p4321_bkgd09->Sumw2();
  h1_Etot_p431_bkgd09->Sumw2();
  h1_Etot_p421_bkgd09->Sumw2();
  h1_Etot_p41_bkgd09->Sumw2();
  h1_Etot_bkgd09_2p1pi_2p0pi->Sumw2();
  h1_Etot_bkgd09_2p1pi_1p1pi->Sumw2();
  h1_Etot_bkgd09_2p1pi_1p0pi->Sumw2();
  h1_Etot_bkgd09_1p2pi_1p0pi->Sumw2();
  h1_Etot_bkgd09_1p2pi_1p1pi->Sumw2();
  h1_Etot_bkgd09_1p3pi->Sumw2();
  h1_Etot_bkgd09_2p2pi->Sumw2();
  h1_Etot_bkgd09_3p1pi->Sumw2();

  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){

      h1_Etot_Npi0_Ecalcut[i][j]->Sumw2();
      h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut41[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut421[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut431[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut4321[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut31[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut321[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[i][j]->Sumw2();
    }
  }




  for(int h=0;h<n_slice;h++){
   h1_Erec_p_bkgd_slice[h]->Sumw2();
   h1_Etot_p_bkgd_slice[h]->Sumw2();
   h1_Erec_3pto1p_slice[h]->Sumw2();
   h1_Etot_3pto1p_slice[h]->Sumw2();
   h1_Erec_3pto2p_slice[h]->Sumw2();
   h1_Etot_3pto2p_slice[h]->Sumw2();
   h1_Erec_3p1pi_slice[h]->Sumw2();
   h1_Etot_3p1pi_slice[h]->Sumw2();
   h1_Etot_43pto1p_slice[h]->Sumw2();
   h1_Erec_43pto1p_slice[h]->Sumw2();
   h1_Erec_4pto3p_slice[h]->Sumw2();
   h1_Etot_4pto3p_slice[h]->Sumw2();
   h1_Erec_4pto2p_slice[h]->Sumw2();
   h1_Etot_4pto2p_slice[h]->Sumw2();
   h1_Erec_4pto1p_slice[h]->Sumw2();
   h1_Etot_4pto1p_slice[h]->Sumw2();
   h1_Etot_Npi0[h]->Sumw2();
   h1_Erec_Npi0_new[h]->Sumw2();
   h1_Etot_bkgd_pipl_pimi_fact[h]->Sumw2();
   h1_Etot_bkgd_pipl_pimi_fact_pipl[h]->Sumw2();
   h1_Etot_bkgd_pipl_pimi_fact_pimi[h]->Sumw2();
   h1_Erec_bkgd_pipl_pimi_new_fact[h]->Sumw2();
   h1_Etot_Npi1[h]->Sumw2();
   h1_Erec_Npi1[h]->Sumw2();
   h1_Etot_bkgd_1p2pi[h]->Sumw2();
   h1_Erec_bkgd_1p2pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p2pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p2pi[h]->Sumw2();
   h1_Etot_bkgd_1p2pi_1p0pi[h]->Sumw2();
   h1_Erec_bkgd_1p2pi_1p0pi[h]->Sumw2();
   h1_Erec_bkgd_1p3pi[h]->Sumw2();
   h1_Etot_bkgd_1p3pi[h]->Sumw2();
  }

  h1_E_rec_1pi_weight_frac_feed->Sumw2();
  h1_E_rec_2pi_weight_frac_feed->Sumw2();
  h1_E_rec_3pi_weight_frac_feed->Sumw2();
  h1_E_rec_4pi_weight_frac_feed->Sumw2();
  h1_E_rec_0pi_frac_feed->Sumw2();
  h1_E_tot_cut2_fracfeed ->Sumw2();
  h1_E_rec_cut2_new_fracfeed ->Sumw2();
  h1_E_tot_p_bkgd_fracfeed ->Sumw2();
  h1_E_rec_p_bkgd_fracfeed ->Sumw2();
  h1_E_tot_2p1pi_2p0pi_fracfeed ->Sumw2();
  h1_E_rec_2p1pi_2p0pi_fracfeed ->Sumw2();
  h1_E_tot_2p1pi_1p1pi_fracfeed ->Sumw2();
  h1_E_rec_2p1pi_1p1pi_fracfeed ->Sumw2();
  h1_E_tot_2p1pi_1p0pi_fracfeed ->Sumw2();
  h1_E_rec_2p1pi_1p0pi_fracfeed ->Sumw2();
  h1_E_tot_1p3pi_fracfeed ->Sumw2();
  h1_E_rec_1p3pi_fracfeed ->Sumw2();
  h1_E_tot_2p2pi_fracfeed ->Sumw2();
  h1_E_rec_2p2pi_fracfeed ->Sumw2();
  h1_E_tot_3p1pi_fracfeed ->Sumw2();
  h1_E_rec_3p1pi_fracfeed ->Sumw2();
  h1_E_tot_3pto2p_fracfeed->Sumw2();
  h1_E_rec_3pto2p_fracfeed->Sumw2();
  h1_E_tot_3pto1p_fracfeed->Sumw2();
  h1_E_rec_3pto1p_fracfeed->Sumw2();
  h1_E_tot_4pto3p_fracfeed->Sumw2();
  h1_E_rec_4pto3p_fracfeed->Sumw2();
  h1_E_tot_43pto1p_fracfeed->Sumw2();
  h1_E_rec_43pto1p_fracfeed->Sumw2();
  h1_E_tot_4pto2p_fracfeed->Sumw2();
  h1_E_rec_4pto2p_fracfeed->Sumw2();
  h1_E_tot_4pto1p_fracfeed->Sumw2();
  h1_E_rec_4pto1p_fracfeed->Sumw2();
  h1_E_tot_1p2pi_fracfeed->Sumw2();
  h1_E_rec_1p2pi_fracfeed->Sumw2();
  h1_E_tot_1p2pi_1p0pi_fracfeed->Sumw2();
  h1_E_rec_1p2pi_1p0pi_fracfeed->Sumw2();
  h1_E_rec_undetfactor_fracfeed->Sumw2();
  h1_E_tot_undetfactor_fracfeed->Sumw2();
  h1_E_rec_2p_det->Sumw2();
  h1_E_tot_2p_det->Sumw2();
  h1_E_tot_p_bkgd->Sumw2();
  h1_E_rec_p_bkgd->Sumw2();
  h1_E_tot_3pto1p->Sumw2();
  h1_E_rec_3pto1p->Sumw2();
  h1_E_tot_43pto1p->Sumw2();
  h1_E_rec_43pto1p->Sumw2();
  h1_E_tot_3pto2p->Sumw2();
  h1_E_rec_3pto2p->Sumw2();
  h1_E_tot_4pto1p->Sumw2();
  h1_E_rec_4pto1p->Sumw2();
  h1_E_tot_4pto3p->Sumw2();
  h1_E_rec_4pto3p->Sumw2();
  h1_E_tot_4pto2p->Sumw2();
  h1_E_rec_4pto2p->Sumw2();
  h1_E_rec->Sumw2();
  h1_E_rec_0pi->Sumw2();
  h1_E_rec_1pi->Sumw2();
  h1_E_rec_4pi_weight->Sumw2();
  h1_E_rec_40pi->Sumw2();
  h1_E_rec_410pi->Sumw2();
  h1_E_rec_420pi->Sumw2();
  h1_E_rec_4210pi->Sumw2();
  h1_E_rec_430pi->Sumw2();
  h1_E_rec_4310pi->Sumw2();
  h1_E_rec_4320pi->Sumw2();
  h1_E_rec_43210pi->Sumw2();
  h1_E_rec_1pi_weight->Sumw2();
  h1_E_rec_2pi_weight->Sumw2();
  h1_E_rec_3pi_weight->Sumw2();
  h1_E_rec_20pi->Sumw2();
  h1_E_rec_21pi->Sumw2();
  h1_E_rec_30pi->Sumw2();
  h1_E_rec_310pi->Sumw2();
  h1_E_rec_320pi->Sumw2();
  h1_E_rec_3210pi->Sumw2();
  h1_E_rec_1prot->Sumw2();
  h1_E_tot_1prot->Sumw2();
  h1_E_rec_cutpi1_piplpimi->Sumw2();
  h1_E_tot_cutpi1_piplpimi->Sumw2();
  h1_Etot->Sumw2();
  h1_E_rec_cut2_new->Sumw2();
  h1_E_tot_cut2->Sumw2();
  h1_E_rec_cut005_newcut3->Sumw2();
  h1_E_rec_undetfactor->Sumw2();
  h1_E_tot_undetfactor->Sumw2();
  h1_E_tot_undetfactor_pipl->Sumw2();
  h1_E_tot_undetfactor_pimi->Sumw2();
  h1_E_tot_1p2pi->Sumw2();
  h1_E_rec_1p2pi->Sumw2();
  h1_E_tot_1p3pi->Sumw2();
  h1_E_rec_1p3pi->Sumw2();
  h1_E_tot_2p2pi->Sumw2();
  h1_E_rec_2p2pi->Sumw2();
  h1_E_tot_3p1pi->Sumw2();
  h1_E_rec_3p1pi->Sumw2();
  h1_E_tot_1p2pi_1p0pi->Sumw2();
  h1_E_rec_1p2pi_1p0pi->Sumw2();
  h1_E_tot_2p1pi_2p0pi->Sumw2();
  h1_E_rec_2p1pi_2p0pi->Sumw2();
  h1_E_tot_2p1pi_1p1pi->Sumw2();
  h1_E_rec_2p1pi_1p1pi->Sumw2();
  h1_E_tot_2p1pi_1p0pi->Sumw2();
  h1_E_rec_2p1pi_1p0pi->Sumw2();

  int CounterEvents = 0;

/** Beginning of Event Loop **/
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<200000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    int nb = GetEntry(jentry);
    if (ientry < 0) break;

    if( jentry%200000 == 0 )
    {
	         gDirectory->Write("hist_Files", TObject::kOverwrite);
		       cout<<jentry<<endl;
    }

    if (runnb==18258 || runnb==18259 || (runnb>18382 && runnb<18438) || (runnb>18220 && runnb<18253)) {
      vz_corr_func->SetParameters(pars); //setting appropriate parameters of the vertex correction function for the runs with the same target and beam energy, but different vertex correction
    }

    if(runnb==18258 || runnb==18259)   {   //setting appropriate e- vertex cut range for the runs with the same target and beam energy, but different vertex correction
      vert_max["56Fe"] = 6.;
      vert_min["56Fe"] = 5.2;//runs with exploaded cell
    }
    if(runnb>18382 && runnb<18438){
      vert_max["3He"] = 0.01;
      vert_min["3He"] = -3.31; //runs with thin exit window
    }
    if(runnb>18220 && runnb<18253){
      vert_max["4He"] = 0.77;
      vert_min["4He"] = -2.27;  //runs with 5cm liquid target cell
    }

    if((runnb>18283 && runnb<18289) || (runnb>18300 && runnb<18304) || (runnb>18317 && runnb<18329))        fTorusCurrent=750;    //setting appropriate torrus magnet current
    else if ((runnb>18293 && runnb<18301) || (runnb>18305 && runnb<18317) || (runnb>18328 && runnb<18336))  fTorusCurrent=1500;
    else fTorusCurrent=2250;

    if(jentry == 0){ //was n_evt == 1 before but jentry = n_evnt - 1
        fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
        fiducialcut->SetFiducialCutParameters(fbeam_en);
        std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
        rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
        std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;
    }

    //Resets q vector to (0,0,0)
    rotation->ResetQVector();

    int n_elec = 0;
    const int ind_em=0; //Index for electron
    if (ec[ind_em] <=0) {
      //      std::cout << "Possible problem with making electron ec vector. EC index below/equal Zero: ec[ind_em] =  " << ec[ind_em] << std::endl;
      continue;
    }
    if (sc[ind_em] <=0) {
/*      std::cout << "Possible problem with making electron ec vector. SC index below/equal zero: sc[ind_em] =  " << sc[ind_em] << std::endl;*/
      continue;
    }
    //Define electron vectors, angles amd other Information
    TVector3 e_ec_xyz1(ech_x[ec[ind_em]-1],ech_y[ec[ind_em]-1],ech_z[ec[ind_em]-1]);
    TVector3 el_mom1(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em] ,p[ind_em]*cz[ind_em]);
    double sc_time = sc_t[sc[ind_em] - 1];
    double sc_path = sc_r[sc[ind_em] - 1];
    int sc_paddle = sc_pd[sc[ind_em] - 1];
    int sc_sector = sc_sect[sc[ind_em] - 1];
    float el_vert = vz[ind_em];
    double ec_x = ech_x[ec[ind_em]-1];
    double ec_y = ech_y[ec[ind_em]-1];
    double ec_z = ech_z[ec[ind_em]-1];
    double el_theta =  TMath::ACos(cz[ind_em])*TMath::RadToDeg();
    double el_phi_mod = TMath::ATan2(cy[ind_em],cx[ind_em])*TMath::RadToDeg()+30; //Add extra 30 degree rotation in phi
    if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree
    int el_ec_sector = ec_sect[ec[ind_em] - 1];
    double el_vert_corr = el_vert+vz_corr(vz_corr_func,el_phi_mod,el_theta);

    //Variables for electron cuts
    double ece = TMath::Max( ec_ei[ec[ind_em] - 1] + ec_eo[ec[ind_em] - 1],   etot[ec[ind_em] - 1]);
    el_segment = int((cc_segm[cc[ind_em]-1]-int(cc_segm[cc[ind_em]-1]/1000)*1000)/10); //does this work in all cases?? F.H. 08/07/19
    el_cc_sector = cc_sect[cc[ind_em]-1];
    el_sccc_timediff = sc_t[cc[ind_em]-1]-cc_t[cc[ind_em]-1]-(sc_r[cc[ind_em]-1]-cc_r[cc[ind_em]-1])/(c*ns_to_s);
    el_cc_nphe = nphe[cc[ind_em]-1]/10.;
    double ec_SC_timediff_uncorr = ec_t[ec[ind_em]-1]-sc_t[sc[ind_em]-1]-(ec_r[ec[ind_em]-1]-sc_r[sc[ind_em]-1])/(c*ns_to_s);

    //fsum_e and fsub_p are TF1 Functions for electron E/p cuts
    fsum_e->SetParameters(epratio_sig_cutrange, max_mom);
    fsub_e->SetParameters(epratio_sig_cutrange, max_mom);

   //Cuts for 1.1 GeV on energy deposition, momenta, tof and cherenkov
    if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2. && ec[ind_em] > 0.5 && sc[ind_em] > 0.5 && cc[ind_em] > 0.5 && q[ind_em] < 0 &&
      ec_ei[ec[ind_em] - 1] >= 0.03 && ece/p[ind_em] >= fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) && p[ind_em] >= min_good_mom  &&
      cc_c2[cc[ind_em]-1] <= 0.1 && el_sccc_timediff >= sc_cc_delt_cut_sect[el_cc_sector-1] )
    {
	     h1_el_SCpd[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);
	     for(int k=1;k<=6;k++){  //k is sector number
	           if(abs(p[ind_em]-0.45)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut[k-1]->Fill(el_phi_mod,el_theta);
	           if(abs(p[ind_em]-1)<0.05 && sc_sect[sc[ind_em]-1]==k)	  h2_el_theta_phi_p_beffidcut2[k-1]->Fill(el_phi_mod,el_theta);
	      }
    }

   //Cuts for 2.2 GeV on energy deposition, momenta, tof and cherenkov
    if(en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 3. && ec[ind_em] > 0.5 && cc[ind_em] > 0.5 &&  sc[ind_em] > 0.5 && q[ind_em] < 0  &&
      ec_ei[ec[ind_em] - 1] >= 0.06 && ece/p[ind_em] >= fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) && p[ind_em] >= min_good_mom &&
      cc_c2[cc[ind_em]-1] < 0.1 && el_sccc_timediff >= sc_cc_delt_cut_sect[el_cc_sector-1]  &&
      TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass) <= en_beam[fbeam_en]) //this cut only exists for 2.2 GeV to cut some very rare events with p_e > beam mom" F.H. 08/08/19
      {
	       h1_el_SCpd[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);
	       for(int k=1;k<=6;k++){ //k is sector number
	          if(abs(p[ind_em]-1.)<0.05 && sc_sect[sc[ind_em]-1]==k)	 h2_el_theta_phi_p_beffidcut[k-1]->Fill(el_phi_mod,el_theta);
	          if(abs(p[ind_em]-1.65)<0.05 && sc_sect[sc[ind_em]-1]==k) h2_el_theta_phi_p_beffidcut2[k-1]->Fill(el_phi_mod,el_theta);
	       }
      }

    //Cuts for 4.4 GeV on energy deposition, momenta, tof and cherenkov.
    // It is the only one with ece cut 0.33, keeping approved cuts
    if(en_beam[fbeam_en] > 4. && en_beam[fbeam_en] < 5. && ec[ind_em] > 0.5 &&  sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] < 0 &&
       ec_ei[ec[ind_em] - 1] >= 0.055 && ece >= 0.33  && ece/p[ind_em] >= fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) && p[ind_em] >= min_good_mom  &&
       cc_c2[cc[ind_em]-1] < 0.1 && el_sccc_timediff >= sc_cc_delt_cut_sect[el_cc_sector-1] )
      {
	       h1_el_SCpd[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);
	       for(int k=1;k<=6;k++){ //k is sector number
	          if(abs(p[ind_em]-2.5)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut[k-1]->Fill(el_phi_mod,el_theta);
	          if(abs(p[ind_em]-1.4)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut2[k-1]->Fill(el_phi_mod,el_theta);
	       }
      }

    h2_el_ec_xy->Fill(ec_x,ec_y);
    //Main Fiducial Cut for Electron
    if( !EFiducialCut(fbeam_en, el_mom1) ) continue;//theta, phi cuts
    if( !CutUVW(e_ec_xyz1) )               continue; //u>60, v<360, w<400

    h2_el_ec_xy_fidcut->Fill(ec_x,ec_y);

    /** FROM HERE AFTER ELECTRON FIDUCIAL CUTS **/
    //Cuts for 4.4 GeV, no cc and tof cuts (only done for 4.4 GeV)
    if(en_beam[fbeam_en] > 4. && en_beam[fbeam_en] < 5. && ec[ind_em] > 0.5 &&  sc[ind_em] > 0.5  &&  q[ind_em] < 0 &&
      ec_ei[ec[ind_em] - 1] >= 0.055  && ece >= 0.33  && p[ind_em] >= min_good_mom)
    {
      h2_el_E_p_ratio_withoutCC->Fill(p[ind_em], ece/p[ind_em]);
      if(cc[ind_em] > 0.5) h2_el_E_p_ratio_withCC->Fill(p[ind_em], ece/p[ind_em]);
    }

    //Cuts for 1.1 GeV. No cut on minimum momentum and ec_ei cut
    if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2. && ec[ind_em] > 0.5 && sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] < 0 &&
       ece/p[ind_em] >= fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) &&
       el_sccc_timediff >= sc_cc_delt_cut_sect[el_cc_sector-1] &&   cc_c2[cc[ind_em]-1] <= 0.1)
    {
      // if(ec_ei[ec[ind_em] - 1] >= 0.05)
      h1_el_Etot_cut->Fill(ece);
      //  if (p[ind_em]>=min_good_mom)
      h1_el_Ein_cut->Fill(ec_ei[ec[ind_em] - 1]);
    }

    //Cuts for 2.2 GeV. No cut on minimum momentum and ec_ei cut
    if(en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 3. && ec[ind_em] > 0.5 && sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] < 0 &&
      ece/p[ind_em] >= fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) &&
      cc_c2[cc[ind_em]-1] < 0.1 && el_sccc_timediff >= sc_cc_delt_cut_sect[el_cc_sector-1]  &&
      TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass) <= en_beam[fbeam_en]) //this cut only exists for 2.2 GeV to cut some very rare events with p_e > beam mom" F.H. 08/08/19
    {
      // if(ec_ei[ec[ind_em] - 1] >= 0.06)
      h1_el_Etot_cut->Fill(ece);
      // if (p[ind_em]>=min_good_mom)
      h1_el_Ein_cut->Fill(ec_ei[ec[ind_em] - 1]);
    }

  //Cuts for 4.4 GeV. No cut on minimum momentum and ec_ei cut
    if(en_beam[fbeam_en] > 4. && en_beam[fbeam_en] < 5. && ec[ind_em] > 0.5 &&  sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] < 0 &&
      ece/p[ind_em] >= fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em])  &&
      cc_c2[cc[ind_em]-1] < 0.1 && el_sccc_timediff >= sc_cc_delt_cut_sect[el_cc_sector-1])
    {
      // if(ec_ei[ec[ind_em] - 1] >= 0.055)
      h1_el_Etot_cut->Fill(ece);
   // if (p[ind_em]>=min_good_mom)
      h1_el_Ein_cut->Fill(ec_ei[ec[ind_em] - 1]);
    }

    //General cut on EC, SC, CC hit and q (charge) for all events
    if( ec[ind_em] < 0.5 ||  sc[ind_em] < 0.5 ||  cc[ind_em] < 0.5 || q[ind_em] >=0)
    {
      continue;
    }

    h1_el_Etot->Fill(ece);
    h1_el_Ein->Fill(ec_ei[ec[ind_em] - 1]);
    h2_el_E_p_ratio->Fill(p[ind_em], ece/p[ind_em]);
    h1_el_cc_chi2->Fill(cc_c2[cc[ind_em]-1]);
    h1_el_cc_deltat[el_cc_sector-1]->Fill(el_sccc_timediff);
    if(el_cc_nphe>2.5)      h1_el_cc_deltat_cut[el_cc_sector-1]->Fill(el_sccc_timediff);
    if(el_vert_corr<vert_max[ftarget] && el_vert_corr>vert_min[ftarget])  h1_el_cc_nphe->Fill(el_cc_nphe);
    if(el_vert_corr<vert_max[ftarget] && el_vert_corr>vert_min[ftarget] && el_sccc_timediff>sc_cc_delt_cut_sect[el_cc_sector-1] &&  cc_c2[cc[ind_em]-1]<0.1)
    {
      h1_el_cc_nphe_cut->Fill(el_cc_nphe);
    }

    h2_el_Ein_Eout->Fill(ec_eo[ec[ind_em]-1],ec_ei[ec[ind_em]-1]);
    h2_el_Einout_Etot->Fill(ece,ec_ei[ec[ind_em]-1]+ec_eo[ec[ind_em]-1]);

    //Cut on 1.1 GeV events (E/p, energy deposit, TOF and cherenkov)
    if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] <2. &&
      ( ec_ei[ec[ind_em] - 1] < 0.03 || ece/p[ind_em] < fsub_e->Eval(p[ind_em]) || ece/p[ind_em] > fsum_e->Eval(p[ind_em]) ||
        p[ind_em] < min_good_mom || el_sccc_timediff < sc_cc_delt_cut_sect[el_cc_sector-1] ||   cc_c2[cc[ind_em]-1] > 0.1 ) )
    {
        continue;
    }

    //Cut on 2.2 GeV events (E/p, energy deposit, TOF and cherenkov)
    if(en_beam[fbeam_en] < 3.  && en_beam[fbeam_en] > 2 &&
      ( ec_ei[ec[ind_em] - 1] < 0.06 || ece/p[ind_em] < fsub_e->Eval(p[ind_em]) || ece/p[ind_em] > fsum_e->Eval(p[ind_em]) ||
        p[ind_em] < min_good_mom || cc_c2[cc[ind_em]-1] >= 0.1 || el_sccc_timediff < sc_cc_delt_cut_sect[el_cc_sector-1] ||
        TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass)>en_beam[fbeam_en] ) ) //only here a cut on electron momentum to cut some very scarse events where p_e > beam energy (see Mariana's anaysis note)
    {
        continue;
    }

    //Cut on 4.4 GeV events (E/p, energy deposit, TOF and cherenkov)
    if(en_beam[fbeam_en] > 4. && en_beam[fbeam_en] < 5. &&
      ( ec_ei[ec[ind_em] - 1] < 0.055 || ece < 0.33 || ece/p[ind_em] < fsub_e->Eval(p[ind_em]) || ece/p[ind_em] > fsum_e->Eval(p[ind_em]) ||
        p[ind_em] < min_good_mom  || cc_c2[cc[ind_em]-1] >= 0.1 || el_sccc_timediff<sc_cc_delt_cut_sect[el_cc_sector-1] ) )
    {
        continue;
    }

    h1_el_SCpdfidcut[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);

    //Plotting electron phi-theta for each sector for each beam energy. Phi is modified (phi_mod = phi + 30)
    //4.4 GeV
    if(en_beam[fbeam_en] > 4. && en_beam[fbeam_en] < 5. ){
      for(int k=1;k<=6;k++){  //k is sector number
	        if(abs(p[ind_em]-2.5)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut[k-1]->Fill(el_phi_mod,el_theta);
	        if(abs(p[ind_em]-1.4)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut2[k-1]->Fill(el_phi_mod,el_theta);
      }
    }
    //2.2 GeV
    else if (en_beam[fbeam_en] < 3. && en_beam[fbeam_en] > 2){
      for(int k=1;k<=6;k++){  //k is sector number
	        if(abs(p[ind_em]-1.)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut[k-1]->Fill(el_phi_mod,el_theta);
	        if(abs(p[ind_em]-1.65)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut2[k-1]->Fill(el_phi_mod,el_theta);
      }
    }
    else { //1.1 GeV
      for(int k=1;k<=6;k++){  //k is sector number
	        if(abs(p[ind_em]-0.45)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut[k-1]->Fill(el_phi_mod,el_theta);
	        if(abs(p[ind_em]-1)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut2[k-1]->Fill(el_phi_mod,el_theta);
      }
    }

    h2_el_E_p_ratio_cut->Fill(p[ind_em], ece/p[ind_em]);
    //Fill histogram for cherenkov photo electrons with vertex cut
    if( el_vert_corr < vert_max[ftarget] && el_vert_corr > vert_min[ftarget]) h1_el_cc_nphe_cut2->Fill(el_cc_nphe);

    if(sc_paddle == 5) {
      h1_el_ec_sc_timediff->Fill(ec_SC_timediff_uncorr);
      h1_el_ec_sc_timediff_corr->Fill(ec_SC_timediff_uncorr-EC_time_offset[make_pair(ftarget,el_ec_sector)]);
    }
    h1_el_ec_sc_timediff_allSCpd->Fill(ec_SC_timediff_uncorr);
    h1_el_ec_sc_timediff_corr_allSCpd->Fill(ec_SC_timediff_uncorr-EC_time_offset[make_pair(ftarget,el_ec_sector)]);
    h1_el_ec_sc_timediff_sect[el_cc_sector-1]->Fill(ec_SC_timediff_uncorr);
    h1_el_ec_sc_timediff_sect_corr[el_cc_sector-1]->Fill(ec_SC_timediff_uncorr-EC_time_offset[make_pair(ftarget,el_ec_sector)]);
    TVector3 v3_el_ec_uvw = FindUVW(e_ec_xyz1);
    h2_el_ec_sc_timediff_ecu[el_cc_sector-1]->Fill(v3_el_ec_uvw.X(),ec_SC_timediff_uncorr);
    h2_el_ec_sc_timediff_ecv[el_cc_sector-1]->Fill(v3_el_ec_uvw.Y(),ec_SC_timediff_uncorr);
    h2_el_ec_sc_timediff_ecw[el_cc_sector-1]->Fill(v3_el_ec_uvw.Z(),ec_SC_timediff_uncorr);
    h2_el_ec_sc_timediff_SCpd[el_cc_sector-1]->Fill(sc_paddle,ec_SC_timediff_uncorr);

    //Electron vertex cut
    if( !(el_vert_corr < vert_max[ftarget] && el_vert_corr > vert_min[ftarget]) ) continue;

    //Main Electron 4-Vectors with and without momentum correction. Units are GeV
    TLorentzVector V4_el_uncorr(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em],p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass));
    TLorentzVector V4_el(elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cx[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cy[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]*elmom_corr_fact[el_ec_sector-1]*elmom_corr_fact[el_ec_sector-1]+e_mass*e_mass));

    h1_el_mom->Fill(V4_el_uncorr.Rho());
    h1_el_mom_corr->Fill(V4_el.Rho());
    h1_el_mom_ratio->Fill(V4_el.Rho()/V4_el_uncorr.Rho());
    h2_el_pcorr_puncorr->Fill(V4_el.Rho(),V4_el.Rho()/V4_el_uncorr.Rho());
    h2_el_mom_diff->Fill(V4_el.Rho(),V4_el.Rho()-V4_el_uncorr.Rho());

    //Electron momentum vector corrected
    TVector3 V3_el=V4_el.Vect();
    double Mott_cross_sec = (fine_struc_const*fine_struc_const*(cz[ind_em]+1)) / (2*V4_el.E()*V4_el.E()*(1-cz[ind_em])*(1-cz[ind_em]));
      //Energy reconstruction for electron only method
  //double E_rec_old = (2*(m_prot-bind_en[ftarget])*V4_el.E()+m_prot*m_prot-(m_prot-bind_en[ftarget])*(m_prot-bind_en[ftarget]))/(2*(m_prot-bind_en[ftarget]-V4_el.E()+V4_el.Rho()*cz[ind_em]));  //using the same value of single nucleon separation E for Ecal and Eqe
    double E_rec= (m_prot*bind_en[ftarget]+m_prot*V4_el.E())/(m_prot-V4_el.E()+V4_el.Rho()*cz[ind_em]);  //using the same value of single nucleon separation E Ecal and Eqe

    //Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
    double nu = -(V4_el-V4_beam).E();
    double Q2 = -(V4_el-V4_beam).Mag2();
    double x_bjk = Q2/(2*m_prot*nu);
    TVector3 V3_q = (V4_beam-V4_el).Vect();
    double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

    //Set q vector for the following rotations for the subtraction procedure
    rotation->SetQVector(V3_q);
//    rotation->PrintQVector();

    //Filling Histogram for electron kinematics
    h1_xbjk->Fill(x_bjk);
    h1_xbjk_weight->Fill(x_bjk,1/Mott_cross_sec);
    h1_Q2->Fill(Q2);
    h1_Q2_weight->Fill(Q2,1/Mott_cross_sec);
    h2_Q2_nu->Fill(nu,Q2);
    h2_Q2_nu_weight->Fill(nu,Q2,1/Mott_cross_sec);
    h2_Q2_xbjk_weight->Fill(x_bjk,Q2,1/Mott_cross_sec);
    h1_Wvar->Fill(W_var);
    h1_Wvar_weight->Fill(W_var,1/Mott_cross_sec);
    h2_Q2_W->Fill(W_var,Q2);
    h2_xB_W->Fill(W_var,x_bjk);
    h2_Q2_W_weight->Fill(W_var,Q2,1/Mott_cross_sec);

    h1_el_Mott_crosssec->Fill(Mott_cross_sec);
    h2_el_theta_phi->Fill(el_phi_mod,el_theta);
    h1_el_theta->Fill(el_theta);
    h2_el_phi_vert_uncorr->Fill(el_vert,el_phi_mod);
    h2_el_phi_vert->Fill(el_vert_corr,el_phi_mod);
    h1_el_vertuncorr->Fill(el_vert);
    h1_el_vertcorr->Fill(el_vert_corr);
    h2_el_vertcorr_runN->Fill(runnb,el_vert_corr);

    //Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

    //Index variables for hadrons (p and pions)
    int index_p[20]; //index for each proton
    int ind_p;      //temp proton index in the loop
    int index_pi[20]; //index for each pion
    int ind_pi_phot[20];
    int ind_pi;       //temp pion index in the loop
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
    int ec_num_n = 0;
    bool ec_radstat_n[20] = {false};

    //Setting arrays
    for (int i = 0; i<20; i++) {
      index_p[i] = -1;   index_pi[i] = -1;   index_pipl[i] = -1;   index_pimi[i] = -1;   ind_pi_phot[i] = -1;
    }

    double pimi_phi, pimi_phi_mod, pimi_theta; //Pi Minus
    double pipl_phi, pipl_phi_mod, pipl_theta; //Pi Plus
    double prot_phi, prot_phi_mod, prot_theta; //Proton
    double pipl_vert_corr, pipl_vert; //Pi Plus Vertex and correction
    double pimi_vert_corr, pimi_vert; //Pi Minus Vertex and correction
    double p_vert_corr; //Proton Vertex Corrected
    double ecpath_corr;
    double neut_phi, neut_phi_mod, neut_theta; //Neutral
    double neut_ecx, neut_ecy, neut_ecz; //neutrals EC hit pos
    double neut_xvert,neut_yvert,neut_zvert; //neutrals Vertex position
    double neut_ecpath_corr, neut_ectime_corr, neut_beta_corr; //Neutrals Corrected values
    double ec_delta = 0; //for neutrals

    const double pimi_vertcut = 2.5; //Vertexcut Pi minus
    const double pipl_vertcut = 2.5; //Vertexcut Pi plus
    const double phot_rad_cut = 40;
    const double phot_e_phidiffcut=30; //electron - photon phi difference cut
    double photon_ece;
    const double EC_sampling_frac = 0.31; //for photons

    CounterEvents ++;

   //Loop for Hadrons, electrons have i=0
    for( int i = 1; i < TMath::Min(gpart, 20); i++ )
    {
        //Start of proton selection
	     if( sc[i] > 0 && stat[i] > 0 &&  id[i] == 2212 ) //Particle i is a proton and has a sc hit
	     {
      	    ind_p = i;
            beta  = p[ind_p]/TMath::Sqrt(p[ind_p]*p[ind_p]+m_prot*m_prot);
	          delta = sc_t[sc[ind_p]-1] - sc_r[sc[ind_p]-1] / (beta*c*ns_to_s) - tr_time;
            prot_phi = TMath::ATan2(cy[ind_p],cx[ind_p])*TMath::RadToDeg();
	          prot_phi_mod = TMath::ATan2(cy[ind_p],cx[ind_p])*TMath::RadToDeg() + 30; //Add extra 30 degree rotation in phi
	          if(prot_phi_mod<0)   prot_phi_mod = prot_phi_mod + 360;   //Proton will be between 0 and 360
            prot_theta = TMath::ACos(cz[ind_p])*TMath::RadToDeg();

            //fsum_e and fsub_p are TF1 Functions for proton delta cuts
	          fsum_prot->SetParameters(prot_delt_cutrange,prot_mom_lim);
	          fsub_prot->SetParameters(prot_delt_cutrange,prot_mom_lim);
	          h2_prot_Deltat_p->Fill(p[ind_p],delta);

            //proton pid cut and momentum > 0.3GeV cut to get rid of low momentum protons that have a high energy loss and we don't know the efficiency precisely for
	          if(delta < fsum_prot->Eval(p[ind_p]) && delta > fsub_prot->Eval(p[ind_p]) && p[ind_p] >= prot_accept_mom_lim){

              TLorentzVector V4_uncorrprot(p[ind_p]*cx[ind_p],p[ind_p]*cy[ind_p],p[ind_p]*cz[ind_p],TMath::Sqrt(p[ind_p]*p[ind_p]+ m_prot*m_prot ) );
              p_vert_corr = vz[ind_p]+vz_corr(vz_corr_func,prot_phi_mod,prot_theta);

	            h2_prot_px_py_p->Fill(cx[ind_p],cy[ind_p]);
	            for(int k=1;k<=6;k++){ // k is sector number
	             if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2. && abs(p[ind_p]-0.6) < 0.025 && sc_sect[sc[ind_p]-1]==k) { //1.1 GeV
                 h2_prot_theta_phi_p_beffidcut[k-1]->Fill(prot_phi_mod,prot_theta);
               }
	             if(en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 5. && abs(p[ind_p]-0.975) < 0.025 && sc_sect[sc[ind_p]-1]==k)	 { //2.2 and 4.4 GeV
                 h2_prot_theta_phi_p_beffidcut[k-1]->Fill(prot_phi_mod,prot_theta);
               }
	            }
	            h2_prot_theta_p[sc_sect[sc[ind_p]-1]-1]->Fill(p[i],prot_theta);

	            if(PFiducialCut(fbeam_en, V4_uncorrprot.Vect())){ //proton fiducial cuts

	              h2_prot_theta_p_cut[sc_sect[sc[ind_p]-1]-1]->Fill(p[i],prot_theta);
	              for(int k = 1; k <= 6; k++){ //k is sector number
		                if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2. && abs(p[ind_p]-0.6) < 0.025 && sc_sect[sc[ind_p]-1]==k) { //1.1 GeV
                      h2_prot_theta_phi_p_fidcut[k-1]->Fill(prot_phi_mod,prot_theta);
                    }
		                if(en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 5. && abs(p[ind_p]-0.975) < 0.025 && sc_sect[sc[ind_p]-1]==k) { //2.2 and 4.4 GeV
                      h2_prot_theta_phi_p_fidcut[k-1]->Fill(prot_phi_mod,prot_theta);
                    }
	               }
	               h1_el_prot_vertdiff_all->Fill(el_vert_corr - p_vert_corr);
	               h2_prot_px_py_p_fidcut->Fill(cx[ind_p],cy[ind_p]);
	               h2_prot_E_p->Fill(p[ind_p],edep[sc[ind_p]-1]);
	               h2_prot_beta_p->Fill(p[ind_p],b[ind_p]);
	               h2_prot_theta_phi->Fill(prot_phi_mod,prot_theta);

                //main vertex cut for protons
	               if( (el_vert_corr-p_vert_corr) > vertdiff_min[ftarget] && (el_vert_corr-p_vert_corr) < vertdiff_max[ftarget] ){
	                   num_p = num_p + 1;
	                   index_p[num_p-1] = i;
	               }
	             } //end of fiducial cuts
             } //end of if delta condition
	        } //end of if loop to check for proton id


	        if(q[i] < 0 && sc[i] > 0 && dc[i] > 0 && stat[i] > 0 ) //negative particle, possibly pi minus
	        {
	             V3_pimi.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	             beta = p[i]/TMath::Sqrt(p[i]*p[i]+m_pimi*m_pimi);
	             delta = sc_t[sc[i]-1]-sc_r[sc[i]-1]/(beta*c*ns_to_s) - tr_time;

               //fsum_pimi and fsub_pimi are TF1 Functions for piminus delta cuts
      	       fsub_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);
	             fsum_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);

	             h1_neg_m->Fill(TMath::Sqrt(p[i]*p[i]/(b[i]*b[i])-p[i]*p[i]));
	             h2_neg_delt_p->Fill(p[i],delta);
	             h2_neg_E_p->Fill(p[i],edep[sc[i]-1]);
	             h2_neg_beta_p->Fill(p[i],b[i]);
               //Pion pid delta cut
          	   if(delta < fsum_pimi->Eval(p[i]) && delta > fsub_pimi->Eval(p[i])){

	                pimi_vert=vz[i];
	                pimi_phi = TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
	                pimi_phi_mod = pimi_phi + 30;  //Add extra 30 degree rotation in phi
	                if (pimi_phi_mod<0)  pimi_phi_mod = pimi_phi_mod + 360;  //Pi minus is between 0 and 360 degree
	                pimi_theta = TMath::ACos(cz[i])*TMath::RadToDeg();
	                pimi_vert_corr = pimi_vert+vz_corr(vz_corr_func, pimi_phi_mod,pimi_theta);

                  //Some if conditions for histograms
                  if(abs(p[i]-1.) < 0.02 && sc_sect[sc[i]-1]==1 && abs(en_beam[fbeam_en]-4.)<1)   h2_pimi_theta_phi_p->Fill(pimi_phi_mod,pimi_theta); //4.4 GeV, why is scsect ==1 here F.H. 12/08/19
	                if(abs(p[i]-1.) < 0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2.)<1)   h2_pimi_theta_phi_p->Fill(pimi_phi_mod,pimi_theta); //2.2 GeV
	                if(abs(p[i]-0.5) <0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1.)<1)   h2_pimi_theta_phi_p->Fill(pimi_phi_mod,pimi_theta); //1.1 GeV
                  for(int k=1;k<=6;k++){ //k is sector number
	                  if(en_beam[fbeam_en]> 1. && en_beam[fbeam_en]< 2. &&  abs(p[i]-0.3)<0.025 && sc_sect[sc[i]-1]==k)  { //1.1 GeV
                      h2_pimi_theta_phi_p_beffidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
                    }
	                  if(en_beam[fbeam_en]> 2. && en_beam[fbeam_en]< 5. && abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k) { //2.2 and 4.4 GeV
                      h2_pimi_theta_phi_p_beffidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
                    }
	                }
	                h2_pimi_theta_phi_beffid->Fill(pimi_phi_mod,pimi_theta);
	                h2_pimi_theta_p[sc_sect[sc[i]-1]-1]->Fill(p[i],pimi_theta);

	                if(PimiFiducialCut(fbeam_en, V3_pimi, &pimi_phimin, &pimi_phimax)){  //Pi minus fiducial cuts

		                h2_pimi_theta_p_cut[sc_sect[sc[i]-1]-1]->Fill(p[i],pimi_theta);
		                h1_pimi_prot_vertdiff->Fill(el_vert_corr-pimi_vert_corr);
                    //some conditions for histogram
		                if(abs(p[i]-1.)< 0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-4.)<1)	h2_pimi_theta_phi_fidcut->Fill(pimi_phi_mod,pimi_theta); //4.4 GeV
		                if(abs(p[i]-1.)< 0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2.)<1)	h2_pimi_theta_phi_fidcut->Fill(pimi_phi_mod,pimi_theta); //2.2 GeV
		                if(abs(p[i]-0.5)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1.)<1)	h2_pimi_theta_phi_fidcut->Fill(pimi_phi_mod,pimi_theta); //1.1 GeV
                    for(int k=1;k<=6;k++){ //k is sector number
		                    if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2. &&  abs(p[i]-0.3)<0.025 && sc_sect[sc[i]-1]==k) { //1.1 GeV
                          h2_pimi_theta_phi_p_fidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
                        }
		                    if(en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 5. &&  abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k) { //2.2 GeV and 4.4 GeV
                           h2_pimi_theta_phi_p_fidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
                        }
	                  }
		                h2_pimi_theta_phi->Fill(pimi_phi_mod,pimi_theta);
                    //main vertex cut for pi minus
		                if(abs(el_vert_corr-pimi_vert_corr) < pimi_vertcut){

		                    num_pimi = num_pimi + 1;
		                    num_pi = num_pi + 1;
		                    num_pi_phot = num_pi_phot + 1;
                        num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
		                    index_pimi[num_pimi - 1] = i;
		                    index_pi[num_pi - 1] = i;
		                    ind_pi_phot[num_pi_phot - 1] = i;

		                    h2_pimi_beta_p->Fill(p[i],b[i]);
		                    h2_pimi_E_p->Fill(p[i],edep[sc[i]-1]);
		                    h2_pimi_delt_p->Fill(p[i],delta);

		                } //if piminus vertex cut
	               } //if Piminus fiducials
	            } //if piminus delta cut
	        } //if negative particle


			    if(q[i] > 0 &&  sc[i] > 0 && dc[i] > 0 && stat[i] > 0) //positive particle. possible pi plus
	        {
	             V3_pipl.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	             beta = p[i]/TMath::Sqrt(p[i]*p[i]+m_pipl*m_pipl);
	             delta= sc_t[sc[i]-1]-sc_r[sc[i]-1]/(beta*c*ns_to_s) - tr_time;

               //fsum_piplus and fsub_piplus are TF1 Functions for piplus delta cuts
	             fsub_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);
	             fsum_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);

	             h1_pos_m->Fill(TMath::Sqrt(p[i]*p[i]/(b[i]*b[i])-p[i]*p[i]));
	             h2_pos_delt_p->Fill(p[i],delta);
	             h2_pos_beta_p->Fill(p[i],b[i]);
	             h2_pos_E_p->Fill(p[i],edep[sc[i]-1]);
               //Pion pid delta cut
               if(delta < fsum_pipl->Eval(p[i]) && delta > fsub_pipl->Eval(p[i])){

	                pipl_vert=vz[i];
	                pipl_phi = TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
            	    pipl_phi_mod = pipl_phi + 30; //Add 30 degrees
	                if (pipl_phi_mod < 0)  pipl_phi_mod = pipl_phi_mod + 360;  //Pi plus is between 0 and 360 degree
	                pipl_theta = TMath::ACos(cz[i])*TMath::RadToDeg();
	                pipl_vert_corr = pipl_vert + vz_corr(vz_corr_func,pipl_phi_mod,pipl_theta);

                  //Some if conditions for histograms
	                if(abs(p[i]-1) < 0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-4)<1)   h2_pipl_theta_phi_p->Fill(pipl_phi_mod,pipl_theta); //4.4 GeV
	                if(abs(p[i]-1) < 0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2)<1)   h2_pipl_theta_phi_p->Fill(pipl_phi_mod,pipl_theta); //2.2 GeV
	                if(abs(p[i]-0.5)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1)<1)   h2_pipl_theta_phi_p->Fill(pipl_phi_mod,pipl_theta); //1.1 GeV
                  for(int k = 1; k <= 6; k++){ //k is sector number
	                   if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2. &&  abs(p[i]-0.3)<0.025 && sc_sect[sc[i]-1]==k) { //1.1 GeV
                       h2_pipl_theta_phi_p_beffidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
                     }
	                   if(en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 5. && abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k) { //2.2 and 4.4 GeV
                       h2_pipl_theta_phi_p_beffidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
                     }
	                }
	                h2_pipl_theta_phi_beffid->Fill(pipl_phi_mod,pipl_theta);
	                h2_pipl_theta_p[sc_sect[sc[i]-1]-1]->Fill(p[i],pipl_theta);

	                if (PiplFiducialCut(fbeam_en, V3_pipl, &pipl_phimin, &pipl_phimax)){ //Pi Plus fiducial cut

		                h2_pipl_theta_p_cut[sc_sect[sc[i]-1]-1]->Fill(p[i],pipl_theta);
		                h1_pipl_prot_vertdiff->Fill(el_vert_corr-pipl_vert_corr);
	       	          if(abs(p[i]-1) < 0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-4)<1)	h2_pipl_theta_phi_fidcut->Fill(pipl_phi_mod,pipl_theta); //4.4 GeV
		                if(abs(p[i]-1) < 0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2)<1)	h2_pipl_theta_phi_fidcut->Fill(pipl_phi_mod,pipl_theta); //2.2 GeV
		                if(abs(p[i]-0.5)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1)<1)	h2_pipl_theta_phi_fidcut->Fill(pipl_phi_mod,pipl_theta); //1.1 GeV
                    for(int k = 1; k <= 6; k++){ //k is sector number
		                    if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] < 2. && abs(p[i]-0.3)<0.025   && sc_sect[sc[i]-1]==k) { //1.1 GeV
                          h2_pipl_theta_phi_p_fidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
                        }
		                    if(en_beam[fbeam_en] > 2. && en_beam[fbeam_en] < 5. && abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k) {//2.2 and 4.4 GeV
                          h2_pipl_theta_phi_p_fidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
                        }
	                  }

		                h2_pipl_theta_phi->Fill(pipl_phi_mod,pipl_theta);

		                if (abs(el_vert_corr-pipl_vert_corr) < pipl_vertcut){ //pi plus vertex cut
		                    num_pipl = num_pipl + 1;
		                    num_pi  = num_pi + 1;
		                    num_pi_phot = num_pi_phot + 1;
                        num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
		                    index_pipl[num_pipl - 1] = i;
		                    index_pi[num_pi - 1] = i;
		                    ind_pi_phot[num_pi_phot - 1] = i;

		                    h2_pipl_beta_p->Fill(p[i],b[i]);
		                    h2_pipl_E_p->Fill(p[i],edep[sc[i]-1]);
		                    h2_pipl_delt_p->Fill(p[i],delta);


		                 }	 //vert cut ends
	                 } //fidcut ends
	              }//delta cut ends
	         }//pipl ends


           if(ec[i] > 0 && dc[i] <= 0  && sc[i] <= 0  && stat[i] > 0 && q[i] == 0) //neutral particles, only EC
	         {

	              neut_zvert = vz[i];
	              neut_yvert = vy[i];
	              neut_xvert = vx[i];
	              neut_ecx = ech_x[ec[i]-1];
	              neut_ecy = ech_y[ec[i]-1];
	              neut_ecz = ech_z[ec[i]-1];
	              TVector3 V3_phot_ec_xyz;
                V3_phot_ec_xyz.SetXYZ(ech_x[ec[i]-1],ech_y[ec[i]-1],ech_z[ec[i]-1]);
	              TVector3 V3_phot_ec_uvw = FindUVW(V3_phot_ec_xyz);

	              neut_ecpath_corr = TMath::Sqrt((neut_ecx-neut_xvert)*(neut_ecx-neut_xvert)+(neut_ecy-neut_yvert)*(neut_ecy-neut_yvert)+(neut_ecz-neut_zvert)*(neut_ecz-neut_zvert));
	              neut_ectime_corr = neut_ecpath_corr/(b[i]*c*ns_to_s) - EC_time_offset[make_pair(ftarget,ec_sect[ec[i]-1])];
	              neut_beta_corr = neut_ecpath_corr/(neut_ectime_corr*c*ns_to_s);
                neut_phi = TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
	              neut_phi_mod = neut_phi + 30; //Add 30 degree
	              if (neut_phi_mod < 0) neut_phi_mod = neut_phi_mod + 360;  //Neutral particle is between 0 and 360 degree
                neut_theta = TMath::ACos(cz[i])*TMath::RadToDeg();

	              V3_phot_angles.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	              ec_delta = ec_t[ec[i]-1] - neut_ecpath_corr/(c*ns_to_s) + EC_time_offset[make_pair(ftarget,ec_sect[ec[i]-1])] - tr_time;

	              h1_beta_ec->Fill(b[i]);
	              h1_beta_ec_corr->Fill(neut_beta_corr);
	              h1_beta_ec_corr_sect[ec_sect[ec[i]-1]-1]->Fill(neut_beta_corr);

	              if(neut_beta_corr > EC_photon_beta[ftarget]){   //photon identification

	                 h2_neutral_theta_phi_EC_all->Fill(neut_phi_mod,neut_theta);
	                 h2_neutral_costheta_phi_EC_all->Fill(neut_phi_mod,cz[i]);

	                 if(Phot_fid(V3_phot_angles)){ //photon fiducial function

                     ec_num_n = ec_num_n + 1;
	                   num_pi_phot = num_pi_phot + 1;
	                   ind_pi_phot[num_pi_phot - 1] = i;
                     //Photon EC energy deposit
 	                   photon_ece = TMath::Max( ec_ei[ec[i] - 1] + ec_eo[ec[i] - 1],etot[ec[i] - 1]);

                     //Cut on Radiation photon via angle with respect to the electron
                     //within 40 degrees in theta and 30 degrees in phi
	                   if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg() < phot_rad_cut && abs(neut_phi_mod-el_phi_mod) < phot_e_phidiffcut ) {
		                     ec_radstat_n[num_pi_phot - 1] = true; //select radiation photons
                         h1_photon_EC_E->Fill(photon_ece/EC_sampling_frac);
                         num_phot_rad = num_phot_rad + 1;
	                   }
	                   if(!ec_radstat_n[num_pi_phot - 1]) num_pi_phot_nonrad = num_pi_phot_nonrad + 1;

	                   h1_time_ec->Fill(ec_delta);
	                   h2_neutral_theta_phi_EC_all_fidcut->Fill(neut_phi_mod,neut_theta);
	                   h1_beta_ec_corr_cut->Fill(neut_beta_corr);
	                   h1_photon_E->Fill(photon_ece/EC_sampling_frac);
	                   h1_phot_e_angle->Fill(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg());
	                   h2_phot_e_angle_vsphotE->Fill(photon_ece/EC_sampling_frac,V3_phot_angles.Angle(V3_el)*TMath::RadToDeg());
	                }//end if photon fiducial
	             }//n beta
	         } //if neutral particles

    } //end of hadron loop

    //Skip event if there is at least one radiation photon
    if (num_phot_rad > 0) {
      continue;
    }

    //Filling Histograms with multiplicities
    h1_Npi->Fill(num_pi);
    h1_Nprot->Fill(num_p);

    if (num_p > 0) {
	h1_Nprot_NonZeroProt->Fill(num_p);
	h1_Npi_NonZeroProt->Fill(num_pi);
    }

    h1_Nphot->Fill(ec_num_n);
    h1_Npipl->Fill(num_pipl);
    h1_Npimi->Fill(num_pimi);
    h1_Npiphot->Fill(num_pi_phot);
    h1_Npiphot_norad->Fill(num_pi_phot_nonrad);
    h2_N_prot_pi->Fill(num_pi,num_p);
    h2_N_prot_pi_phot->Fill(num_pi+ec_num_n,num_p);
    h2_N_prot_pi_phot_nonrad->Fill(num_pi_phot_nonrad,num_p);
    h2_N_pi_phot[num_p]->Fill(ec_num_n,num_pi);

    int ind_prot1,ind_prot2,ind_phot;

    //Events with exactly 2 protons
	  if(num_p == 2)
    {

        ind_prot1 = index_p[0];
        ind_prot2 = index_p[1];
        //LorentzVectors for protons without momentum correction
        TLorentzVector V4_prot_uncorr1(p[ind_prot1]*cx[ind_prot1],p[ind_prot1]*cy[ind_prot1],p[ind_prot1]*cz[ind_prot1],TMath::Sqrt(m_prot*m_prot+p[ind_prot1]*p[ind_prot1]));
        TLorentzVector V4_prot_uncorr2(p[ind_prot2]*cx[ind_prot2],p[ind_prot2]*cy[ind_prot2],p[ind_prot2]*cz[ind_prot2],TMath::Sqrt(m_prot*m_prot+p[ind_prot2]*p[ind_prot2]));
        //Kinematic for first proton
        double prot_vert1 = vz[ind_prot1];
        double p_phi1=TMath::ATan2(cy[ind_prot1],cx[ind_prot1])*TMath::RadToDeg();
        double p_phi_mod1=p_phi1+30;
        if (p_phi_mod1<0)p_phi_mod1=p_phi_mod1+360;
        double p_theta1=TMath::ACos(cz[ind_prot1])*TMath::RadToDeg();
        double prot_vert_corr1=prot_vert1+vz_corr(vz_corr_func,p_phi_mod1,p_theta1);
        double prot_mom_corr1;
        if(ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr1,prot_vert_corr1) != -1)
        {
          prot_mom_corr1=ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr1,prot_vert_corr1);
        }
        else
        {
          prot_mom_corr1=p[ind_prot1];
        }

     //Kinematic for second proton
        double prot_vert2=vz[ind_prot2];
        double p_phi2=TMath::ATan2(cy[ind_prot2],cx[ind_prot2])*TMath::RadToDeg();
        double p_phi_mod2=p_phi2+30;
        if (p_phi_mod2<0)p_phi_mod2=p_phi_mod2+360;
        double p_theta2=TMath::ACos(cz[ind_prot2])*TMath::RadToDeg();
        double prot_vert_corr2=prot_vert2+vz_corr(vz_corr_func,p_phi_mod2,p_theta2);
        double prot_mom_corr2;
        if(ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr2,prot_vert_corr2) != -1)
        {
          prot_mom_corr2=ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr2,prot_vert_corr2);
        }
        else
        {
          prot_mom_corr2=p[ind_prot2];
        }

        h1_el_prot_vertdiff1->Fill(el_vert_corr-prot_vert_corr1);
        h1_el_prot_vertdiff2->Fill(el_vert_corr-prot_vert_corr2);
        //Electron and Proton Vertex Differnce cut
        if( (el_vert_corr-prot_vert_corr1) > vertdiff_min[ftarget] && (el_vert_corr-prot_vert_corr1) < vertdiff_max[ftarget] &&
            (el_vert_corr-prot_vert_corr2) > vertdiff_min[ftarget] && (el_vert_corr-prot_vert_corr2) < vertdiff_max[ftarget])
        {

          TVector3 V3_2prot_uncorr[2];
          V3_2prot_uncorr[0].SetXYZ(p[ind_prot1]*cx[ind_prot1],p[ind_prot1]*cy[ind_prot1],p[ind_prot1]*cz[ind_prot1]);
          V3_2prot_uncorr[1].SetXYZ(p[ind_prot2]*cx[ind_prot2],p[ind_prot2]*cy[ind_prot2],p[ind_prot2]*cz[ind_prot2]);

          TVector3 V3_2prot_corr[2];
          V3_2prot_corr[0].SetXYZ(prot_mom_corr1*cx[ind_prot1],prot_mom_corr1*cy[ind_prot1],prot_mom_corr1*cz[ind_prot1]);
          V3_2prot_corr[1].SetXYZ(prot_mom_corr2*cx[ind_prot2],prot_mom_corr2*cy[ind_prot2],prot_mom_corr2*cz[ind_prot2]);

          TLorentzVector V4_prot_corr1(V3_2prot_corr[0],TMath::Sqrt(m_prot*m_prot+prot_mom_corr1*prot_mom_corr1));
          TLorentzVector V4_prot_corr2(V3_2prot_corr[1],TMath::Sqrt(m_prot*m_prot+prot_mom_corr2*prot_mom_corr2));
          TLorentzVector V4_prot_el_tot1=V4_prot_corr1+V4_el;
          TLorentzVector V4_prot_el_tot2=V4_prot_corr2+V4_el;

          h1_Wepp->Fill((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2-V4_el).M());
          h1_Wepp_uncorr->Fill((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2-V4_el_uncorr).M());

          h2_Wepp_ephi_corr->Fill(el_phi_mod,(V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2-V4_el).M());
          h2_Wepp_ephi->Fill(el_phi_mod,(V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2-V4_el_uncorr).M());

          h2_Wepp_ephi_corr_uncorrprot->Fill(el_phi_mod,(V4_target+V4_beam-V4_prot_uncorr1-V4_prot_uncorr2-V4_el).M());
          h2_Wepp_ephi_uncorrprot->Fill(el_phi_mod,(V4_target+V4_beam-V4_prot_uncorr1-V4_prot_uncorr2-V4_el_uncorr).M());

          double mult = 0.5*((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2).Mag2()-m_neut*m_neut)/((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2).E()*V4_el_uncorr.Rho()-((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2).Vect()).Dot(V4_el_uncorr.Vect()));

          h1_e_mom_corrfuct[sc_sect[sc[ind_em]-1]-1]->Fill(mult);

 //---------------------------------- 2p 0pi->  1p0pi   ----------------------------------------------

          double E_tot_2p[2]={0};
          double p_perp_tot_2p[2]={0};
          N_prot_both = 0;
          double P_N_2p[2]={0};
          V3_q=(V4_beam-V4_el).Vect();
          rotation->prot2_rot_func(V3_2prot_corr, V3_2prot_uncorr, V4_el, E_tot_2p, p_perp_tot_2p, P_N_2p , &N_prot_both);

          if(num_pi_phot==0 && N_prot_both!=0){
            for(int f = 0; f < num_p; f++){    //looping through two protons

              h1_E_tot_p_bkgd->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);
              h1_E_rec_p_bkgd->Fill(E_rec,P_N_2p[f]*1/Mott_cross_sec);
              h1_E_tot_p_bkgd09->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);
              h2_Erec_pperp_2p->Fill(p_perp_tot_2p[f],E_rec,P_N_2p[f]*1/Mott_cross_sec);
              h1_E_tot_p_bkgd_fracfeed->Fill((E_tot_2p[f]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_N_2p[f]*1/Mott_cross_sec);
              h1_E_rec_p_bkgd_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_N_2p[f]*1/Mott_cross_sec);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[f]) *TMath::RadToDeg(),-P_N_2p[f]*1/Mott_cross_sec);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);

              for (int i = 0; i < N_pperp; i++) {
	               for(int j = 0; j < N_Ecal; j++) {
	                  if(E_tot_2p[f] > Ecal_lowlim[j] && E_tot_2p[f] < Ecal_uplim[j] && p_perp_tot_2p[f] > pperp_cut[i])
                    {
                       h1_Etot_p_bkgd_slice_Ecalcut[i][j]->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);
                    }
	               }
              }
              for(int i = 0 ; i < n_slice; i++)
              {
	               if (p_perp_tot_2p[f] < pperp_max[i] && p_perp_tot_2p[f] > pperp_min[i]) {
	                  h1_Etot_p_bkgd_slice[i]->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);
	                  h1_Erec_p_bkgd_slice[i]->Fill(E_rec,P_N_2p[f]*1/Mott_cross_sec);
	               }
	            }
            }//looping through two protons

            h1_E_tot_2p_det->Fill(E_tot_2p[0],1/Mott_cross_sec);
            h1_E_rec_2p_det->Fill(E_rec,1/Mott_cross_sec);
          }//no pions cut and N_prot_both!=0


 //---------------------------------- 2p 1pi   ----------------------------------------------
          const int N_2prot=2;
          TVector3 V3_1pi,V3_2p_rotated[N_2prot],V3_1pirot;
          bool pi1_stat=false;
          double N_2p_0pi=0,N_all=0,N_1p_1pi[N_2prot]={0},N_1p_0pi[N_2prot]={0};
          double Ecal_2p1pi_to2p0pi[N_2prot]={0},p_miss_perp_2p1pi_to2p0pi[N_2prot]={0}, P_2p1pi_to2p0pi[N_2prot]={0},N_2p_det=0;
          double N_pidet=0,N_piundet=0;

          if (num_pi_phot==1) {

            V3_1pi.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);

            double P_2p1pito2p0pi[2]={0};
            double P_2p1pito1p1pi[2]={0};
            double P_2p1pito1p0pi[2]={0};
            double Ptot=0;

            rotation->prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_1pi, q[ind_pi_phot[0]],V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

            for(int z=0; z < N_2prot; z++){ //looping over two protons

 //---------------------------------- 2p 1pi ->2p 0pi   ----------------------------------------------

              h1_E_tot_2p1pi_2p0pi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h1_E_rec_2p1pi_2p0pi->Fill(E_rec,P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h2_Erec_pperp_2p1pi_2p0pi->Fill(p_miss_perp_2p1pi_to2p0pi[z],E_rec,P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h1_Etot_bkgd09_2p1pi_2p0pi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h1_E_tot_2p1pi_2p0pi_fracfeed->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h1_E_rec_2p1pi_2p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h2_pperp_W->Fill(W_var,p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z]) *TMath::RadToDeg(),P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h2_Ecal_Eqe->Fill(E_rec,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);

              for(int i = 0; i < n_slice; i++){
                if (p_miss_perp_2p1pi_to2p0pi[z]<pperp_max[i] && p_miss_perp_2p1pi_to2p0pi[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
	                 h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(E_rec,P_2p1pito2p0pi[z]*1/Mott_cross_sec);
                }
              }
              for (int i = 0; i < N_pperp; i++){
                for(int j = 0; j < N_Ecal; j++){
	                 if(Ecal_2p1pi_to2p0pi[z] > Ecal_lowlim[j] && Ecal_2p1pi_to2p0pi[z] < Ecal_uplim[j] && p_miss_perp_2p1pi_to2p0pi[z] > pperp_cut[i]) {
                       h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[i][j]->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
                   }
                }
              }
 //---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

              h1_E_tot_2p1pi_1p1pi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h1_E_rec_2p1pi_1p1pi->Fill(E_rec,P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h2_Erec_pperp_2p1pi_1p1pi->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h1_Etot_bkgd09_2p1pi_1p1pi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h1_E_tot_2p1pi_1p1pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h1_E_rec_2p1pi_1p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);

              for(int i = 0; i < n_slice; i++){
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
                  h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
                  h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_rec,P_2p1pito1p1pi[z]*1/Mott_cross_sec);
                }
              }
              for (int i = 0; i < N_pperp; i++){
                for(int j = 0; j < N_Ecal; j++){
                  if(E_tot_2p[z] > Ecal_lowlim[j] && E_tot_2p[z] < Ecal_uplim[j] && p_perp_tot_2p[z] > pperp_cut[i])  {
                     h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[i][j]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
                  }
                }
              }
//---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------

              h1_E_tot_2p1pi_1p0pi->Fill(E_tot_2p[z], P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h1_E_rec_2p1pi_1p0pi->Fill(E_rec,P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h2_Erec_pperp_2p1pi_1p0pi->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h1_Etot_bkgd09_2p1pi_1p0pi->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h1_E_tot_2p1pi_1p0pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h1_E_rec_2p1pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);

              for(int i = 0; i < n_slice; i++) {
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
	                 h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_rec,P_2p1pito1p0pi[z]*1/Mott_cross_sec);
                }
              }
              for (int i = 0;i < N_pperp; i++){
                for(int j = 0; j < N_Ecal; j++){
                  if(E_tot_2p[z] > Ecal_lowlim[j] && E_tot_2p[z] < Ecal_uplim[j] && p_perp_tot_2p[z] > pperp_cut[i]) {
                      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[i][j]->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
                  }
                }
              }

            }//filling the histograms for 2protons
          }//1pi requirement


//---------------------------------- 2p 2pi   ----------------------------------------------
          const int N_2pi=2;
          int q_pi2[2];
          TVector3 V3_2pi[N_2pi];
          double Ecal_2p2pi[N_2prot],p_miss_perp_2p2pi[N_2prot],Ptot_2p[2]={0};
          bool ecstat_pi2[N_2pi]={false};

          if ( num_pi_phot == 2) {

            V3_2pi[0].SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
            V3_2pi[1].SetXYZ(p[ind_pi_phot[1]]*cx[ind_pi_phot[1]],p[ind_pi_phot[1]]*cy[ind_pi_phot[1]],p[ind_pi_phot[1]]*cz[ind_pi_phot[1]]);
            q_pi2[0] = q[ind_pi_phot[0]];
            q_pi2[1] = q[ind_pi_phot[1]];
            ecstat_pi2[0] = ec_radstat_n[0];
            ecstat_pi2[1] = ec_radstat_n[1];

            rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi,q_pi2 ,V4_el, Ecal_2p2pi,p_miss_perp_2p2pi,Ptot_2p);

            for(int z = 0; z < N_2prot; z++){ //looping over two protons

//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

              h1_E_tot_2p2pi->Fill(E_tot_2p[z], Ptot_2p[z]*1/Mott_cross_sec);
              h1_E_rec_2p2pi->Fill(E_rec,Ptot_2p[z]*1/Mott_cross_sec);
              h2_Erec_pperp_2p2pi->Fill(p_perp_tot_2p[z],E_rec,Ptot_2p[z]*1/Mott_cross_sec);
              h1_Etot_bkgd09_2p2pi->Fill(E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
              h1_E_tot_2p2pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],Ptot_2p[z]*1/Mott_cross_sec);
              h1_E_rec_2p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],Ptot_2p[z]*1/Mott_cross_sec);
              h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
              h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),Ptot_2p[z]*1/Mott_cross_sec);
              h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
              h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
              h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);

              for(int i = 0; i < n_slice; i++)
              {
                if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	                 h1_Etot_p_bkgd_slice_2p2pi[i]->Fill(E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
	                 h1_Erec_p_bkgd_slice_2p2pi[i]->Fill(E_rec,Ptot_2p[z]*1/Mott_cross_sec);
                }
              }
              for (int i = 0;i < N_pperp; i++){
                for(int j = 0; j < N_Ecal; j++){
                  if(E_tot_2p[z]>Ecal_lowlim[j] && E_tot_2p[z]<Ecal_uplim[j] && p_perp_tot_2p[z]>pperp_cut[i]) {
                      h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[i][j]->Fill(E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
                  }
                }
              }
            }   //Filling the histogram for two protons
          }//2pi requirement
        }//2prot vert cut
    } //2prot requirement

  //Events with exactly 3 protons
	  if(num_p == 3){

	      const int N_3p=3;
	      TLorentzVector V4_p_uncorr[N_3p], V4_p_corr[N_3p],V4_prot_el[N_3p];
	      float prot_vz[N_3p];
	      double prot_phi[N_3p],prot_theta[N_3p],prot_vz_corr[N_3p],prot_p_corr[N_3p];
	      TVector3 V3_prot[N_3p],V3_prot_corr[N_3p],V3_3p_rot[N_3p];
	      double E_cal[N_3p],p_miss_perp[N_3p],P_3pto1p[N_3p];
	      double N_p1[N_3p]={0};
        double N_p_three=0;
	      double N_p12[N_3p]={0};
        double N_p13[N_3p]={0};
        double N_p23[N_3p]={0};
        double N_p_two=0;
        double E_cal_3pto1p[3]={0};
        double p_miss_perp_3pto1p[3]={0};
	      int N_comb=3;
	      const int N_2p=2;
	      double E_cal_3pto2p[3][N_2p]={0};
        double p_miss_perp_3pto2p[3][N_2p]={0};
        double P_3pto2p[3][N_2p]={0};
	      TVector3 V3_prot_el_3pto2p[N_3p][N_2p];
	      TVector3 V3_2p_rot[N_2p], V3_prot_el[N_3p][N_2p];
 	      bool prot_stat[N_3p]={false};

	      for(int i = 0; i < N_3p; i++)
	      {
	        V4_p_uncorr[i].SetPxPyPzE(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+p[index_p[i]]*p[index_p[i]]));
	        prot_vz[i]=vz[index_p[i]];
	        prot_phi[i]=TMath::ATan2(cy[index_p[i]],cx[index_p[i]])*TMath::RadToDeg()+30;
	        if(prot_phi[i]<0)prot_phi[i]=prot_phi[i]+360;
	        prot_theta[i]=TMath::ACos(cz[index_p[i]])*TMath::RadToDeg();
	        prot_vz_corr[i]=prot_vz[i]+vz_corr(vz_corr_func,prot_phi[i],prot_theta[i]);

	        if(ProtonMomCorrection_He3_4Cell(ftarget,V4_p_uncorr[i],prot_vz_corr[i]) != -1){
		            prot_p_corr[i]=ProtonMomCorrection_He3_4Cell(ftarget,V4_p_uncorr[i],prot_vz_corr[i]);
          }
	        else
		      {
            prot_p_corr[i]=p[index_p[i]];
          }

	        V4_p_corr[i].SetPxPyPzE(prot_p_corr[i]*cx[index_p[i]],prot_p_corr[i]*cy[index_p[i]],prot_p_corr[i]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+prot_p_corr[i]*prot_p_corr[i]));
	        V3_prot[i].SetXYZ(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]]);
	        V3_prot_corr[i].SetXYZ(prot_p_corr[i]*cx[index_p[i]],prot_p_corr[i]*cy[index_p[i]],prot_p_corr[i]*cz[index_p[i]]);

	        V4_prot_el[i]=V4_p_corr[i]+V4_el;
	        E_cal[i]=V4_el.E()+ V4_p_corr[i].E()-m_prot+bind_en[ftarget];
	        p_miss_perp[i]=TMath::Sqrt(V4_prot_el[i].Px()*V4_prot_el[i].Px()+V4_prot_el[i].Py()*V4_prot_el[i].Py());
	      }
        V3_prot_el[0][0]=V4_el.Vect()+V3_prot[0];
        V3_prot_el[0][1]=V4_el.Vect()+V3_prot[1];
        V3_prot_el[1][0]=V4_el.Vect()+V3_prot[0];
        V3_prot_el[1][1]=V4_el.Vect()+V3_prot[2];
        V3_prot_el[2][0]=V4_el.Vect()+V3_prot[1];
        V3_prot_el[2][1]=V4_el.Vect()+V3_prot[2];

	      h1_el_3prot_vertdiff1->Fill(el_vert_corr- prot_vz_corr[0]);
	      h1_el_3prot_vertdiff2->Fill(el_vert_corr- prot_vz_corr[1]);
	      h1_el_3prot_vertdiff3->Fill(el_vert_corr- prot_vz_corr[2]);

        //Electron and Proton vertex difference cut
	      if(   (el_vert_corr- prot_vz_corr[0])>vertdiff_min[ftarget] && (el_vert_corr- prot_vz_corr[0])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot_vz_corr[1])>vertdiff_min[ftarget] && (el_vert_corr- prot_vz_corr[1])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot_vz_corr[2])>vertdiff_min[ftarget] && (el_vert_corr- prot_vz_corr[2])<vertdiff_max[ftarget])
        {

	         for(int i = 0; i < N_3p; i++){
	            N_p1[i]=0;
	         }
	         V3_q=(V4_beam-V4_el).Vect();

           rotation->prot3_rot_func( V3_prot_corr,V3_prot,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

	         if(num_pi_phot==0 && N_p_three!=0){
	           for(int count = 0; count < N_comb;count++)    { //Loop over number of combinations
               for(int j = 0; j < N_2p; j++)    { //loop over two protons

 //-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------

                  h1_E_tot_3pto2p->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*1/Mott_cross_sec);
		              h1_E_rec_3pto2p->Fill(E_rec, P_3pto2p[count][j]*1/Mott_cross_sec);
		              h1_Etot_p321_bkgd09->Fill(E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
		              h2_Erec_pperp_321p->Fill(p_miss_perp_3pto2p[count][j],E_rec,P_3pto2p[count][j]*1/Mott_cross_sec);
		              h1_E_tot_3pto2p_fracfeed->Fill((E_cal_3pto2p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto2p[count][j]*1/Mott_cross_sec);
		              h1_E_rec_3pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_3pto2p[count][j]*1/Mott_cross_sec);
		              h2_pperp_W->Fill(W_var,p_miss_perp_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
		              h1_theta0->Fill((V4_beam.Vect()).Angle(V3_prot_el[count][j])*TMath::RadToDeg(),P_3pto2p[count][j]*1/Mott_cross_sec);
	                h2_Ecal_Eqe->Fill(E_rec,E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
	                h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
	                h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);

		              for (int n = 0; n < N_pperp;n++){
		                  for(int z = 0; z < N_Ecal; z++){
		                      if(E_cal_3pto2p[count][j]>Ecal_lowlim[z] && E_cal_3pto2p[count][j]<Ecal_uplim[z] && p_miss_perp_3pto2p[count][j]>pperp_cut[n]) {
                              h1_Etot_p_bkgd_slice_Ecalcut321[n][z]->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*1/Mott_cross_sec);
                          }
		                  }
		              }

		              for(int i = 0; i < n_slice; i++)
		              {
		                  if (p_miss_perp_3pto2p[count][j]<pperp_max[i] && p_miss_perp_3pto2p[count][j]>pperp_min[i]){
	                         h1_Etot_3pto2p_slice[i]->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*1/Mott_cross_sec);
			                     h1_Erec_3pto2p_slice[i]->Fill(E_rec, P_3pto2p[count][j]*1/Mott_cross_sec);
		                  }
		              }
		            } //end loop over protons
	           } //end loop over combination N_comb

 //-----------------------------------------  3p to 1p  -----------------------------------------------------------------------
	           for(int j = 0; j < N_3p; j++)    {

		             P_3pto1p[j]= N_p1[j]/N_p_three;
		             h1_E_tot_3pto1p->Fill(E_cal[j], P_3pto1p[j]*1/Mott_cross_sec);
		             h1_E_rec_3pto1p->Fill(E_rec,P_3pto1p[j]*1/Mott_cross_sec);
		             h1_Etot_p31_bkgd09->Fill(E_cal[j],P_3pto1p[j]*1/Mott_cross_sec);
		             h2_Erec_pperp_31p->Fill(p_miss_perp[j],E_rec,P_3pto1p[j]*1/Mott_cross_sec);
		             h1_E_tot_3pto1p_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto1p[j]*1/Mott_cross_sec);
		             h1_E_rec_3pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_3pto1p[j]*1/Mott_cross_sec);
		             h2_pperp_W->Fill(W_var,p_miss_perp[j],-P_3pto1p[j]*1/Mott_cross_sec);
		             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot[j])*TMath::RadToDeg(),-P_3pto1p[j]*1/Mott_cross_sec);
		             h2_Ecal_Eqe->Fill(E_rec,E_cal[j],-P_3pto1p[j]*1/Mott_cross_sec);
		             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],-P_3pto1p[j]*1/Mott_cross_sec);
		             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],-P_3pto1p[j]*1/Mott_cross_sec);

		             for (int n = 0; n < N_pperp; n++){
		                 for(int z = 0; z < N_Ecal; z++){
		                     if(E_cal[j] > Ecal_lowlim[z] && E_cal[j] < Ecal_uplim[z] && p_miss_perp[j]>pperp_cut[n]) {
                             h1_Etot_p_bkgd_slice_Ecalcut31[n][z]->Fill(E_cal[j], P_3pto1p[j]*1/Mott_cross_sec);
                         }
		                 }
		             }

		             for(int i = 0; i < n_slice; i++)
		             {
		                 if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
		                     h1_Etot_3pto1p_slice[i]->Fill(E_cal[j],P_3pto1p[j]*1/Mott_cross_sec);
		                     h1_Erec_3pto1p_slice[i]->Fill(E_rec,P_3pto1p[j]*1/Mott_cross_sec);
		                 }
		             }
	           } //end loop over N_3p
	         }//end if num_pi_phot==0 && N_p_three!=0, no pions

//----------------------------------3p 1pi ----------------------------------------------------------

           if (num_pi_phot==1) { //number of pions or photons = 1

             double P_tot_3p[N_3p]={0};
             double Ecal_3p1pi[N_3p]={0};
             double p_miss_perp_3p1pi[N_3p]={0};
             TVector3 V3_pi_phot;
             V3_q=(V4_beam-V4_el).Vect();
             V3_pi_phot.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);

             rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot, V3_pi_phot,q[ind_pi_phot[0]], V4_el,  Ecal_3p1pi,p_miss_perp_3p1pi, P_tot_3p);

	           for(int j = 0; j < N_3p; j++)    { //loop over 3 protons

                 h1_E_tot_3p1pi->Fill(E_cal[j], P_tot_3p[j]*1/Mott_cross_sec);
		             h1_E_rec_3p1pi->Fill(E_rec,P_tot_3p[j]*1/Mott_cross_sec);
		             h1_Etot_bkgd09_3p1pi->Fill(E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		             h2_Erec_pperp_3p1pi->Fill(p_miss_perp[j],E_rec,P_tot_3p[j]*1/Mott_cross_sec);
		             h1_E_tot_3p1pi_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_tot_3p[j]*1/Mott_cross_sec);
		             h1_E_rec_3p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_3p[j]*1/Mott_cross_sec);
		             h2_pperp_W->Fill(W_var,p_miss_perp[j],P_tot_3p[j]*1/Mott_cross_sec);
		             h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot[j])*TMath::RadToDeg(),P_tot_3p[j]*1/Mott_cross_sec);
		             h2_Ecal_Eqe->Fill(E_rec,E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		             h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);

		             for (int n = 0; n < N_pperp; n++){
		                 for(int z = 0; z < N_Ecal; z++){
		                     if(E_cal[j]>Ecal_lowlim[z] && E_cal[j]<Ecal_uplim[z] && p_miss_perp[j]>pperp_cut[n]) {
                             h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[n][z]->Fill(E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
                         }
		                 }
		             }
		             for(int i = 0; i < n_slice; i++)
		             {
		                 if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
                         h1_Etot_3p1pi_slice[i]->Fill(E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		                     h1_Erec_3p1pi_slice[i]->Fill(E_rec,P_tot_3p[j]*1/Mott_cross_sec);
		                 }
		             }
             } //end loop over N_3p
           }// 1 pi requirement ends
	      }//vertex cut for all three protons
	  }//end if num_p == 3  3proton requiremnet

//Events with exactly 4 protons
	 	 if(num_p==4){

	      const int N_p4=4;
	      TLorentzVector V4_p4_uncorr[N_p4], V4_p4_corr[N_p4],V4_prot4_el[N_p4];
	      float prot4_vz[N_p4];
	      double prot4_phi[N_p4],prot4_theta[N_p4],prot4_vz_corr[N_p4],prot4_p_corr[N_p4];
	      TVector3 V3_prot4[N_p4],V3_prot4_corr[N_p4];
	      double E_cal_p4[N_p4]={0};
        double p_miss_perp_p4[N_p4]={0};
        double P_4pto1p[N_p4]={0};
	      TVector3 V3_p4_rot[N_p4];
	      bool prot4_stat[N_p4]={false};
	      const int Ncomb_4to1 = 4,Ncomb_4to2 = 6, Ncomb_4to3 = 4;
	      double N_p4_p1[Ncomb_4to1]={0};
        double N_p4_p2[Ncomb_4to2]={0};
        double N_p4_p3[Ncomb_4to3]={0};
        double N_p_four=0;

	      for(int i = 0; i < N_p4; i++)  //loop over 4 protons
	      {
	         V4_p4_uncorr[i].SetPxPyPzE(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+p[index_p[i]]*p[index_p[i]]));
	         prot4_vz[i]=vz[index_p[i]];
	         prot4_phi[i]=TMath::ATan2(cy[index_p[i]],cx[index_p[i]])*TMath::RadToDeg()+30;
	         if(prot4_phi[i]<0)prot4_phi[i]=prot4_phi[i]+360;
	         prot4_theta[i]=TMath::ACos(cz[index_p[i]])*TMath::RadToDeg();
	         prot4_vz_corr[i]=prot4_vz[i]+vz_corr(vz_corr_func,prot4_phi[i],prot4_theta[i]);

	         if(ProtonMomCorrection_He3_4Cell(ftarget,V4_p4_uncorr[i],prot4_vz_corr[i]) != -1)
           {
		           prot4_p_corr[i]=ProtonMomCorrection_He3_4Cell(ftarget,V4_p4_uncorr[i],prot4_vz_corr[i]);
           }
	         else
		       {
               prot4_p_corr[i]=p[index_p[i]];
           }

	         V3_prot4[i].SetXYZ(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]]);
	         V3_prot4_corr[i].SetXYZ(prot4_p_corr[i]*cx[index_p[i]],prot4_p_corr[i]*cy[index_p[i]],prot4_p_corr[i]*cz[index_p[i]]);
	         V4_p4_corr[i].SetPxPyPzE(prot4_p_corr[i]*cx[index_p[i]],prot4_p_corr[i]*cy[index_p[i]],prot4_p_corr[i]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+prot4_p_corr[i]*prot4_p_corr[i]));
	         V4_prot4_el[i]=V4_p4_corr[i]+V4_el;
	         E_cal_p4[i]=V4_el.E()+ V4_p4_corr[i].E()-m_prot+bind_en[ftarget];
	         p_miss_perp_p4[i]=TMath::Sqrt(V4_prot4_el[i].Px()*V4_prot4_el[i].Px()+V4_prot4_el[i].Py()*V4_prot4_el[i].Py());
	      } //end loop over 4 protons

	      h1_el_4prot_vertdiff1->Fill(el_vert_corr- prot4_vz_corr[0]);
	      h1_el_4prot_vertdiff2->Fill(el_vert_corr- prot4_vz_corr[1]);
	      h1_el_4prot_vertdiff3->Fill(el_vert_corr- prot4_vz_corr[2]);
	      h1_el_4prot_vertdiff4->Fill(el_vert_corr- prot4_vz_corr[3]);

        //Electron and Proton vertex difference cut
	      if(   (el_vert_corr- prot4_vz_corr[0])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[0])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot4_vz_corr[1])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[1])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot4_vz_corr[2])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[2])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot4_vz_corr[3])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[3])<vertdiff_max[ftarget])
        {

	         V3_q=(V4_beam-V4_el).Vect();

           if ( num_pi_phot == 0){ //no pion or photon
	            for(int g = 0; g < N_tot; g++){ //this looks like a 4-proton rotation function -> could be placed maybe in an extra function

		              double rot_angle = gRandom->Uniform(0,2*TMath::Pi());
		              for(int i = 0; i < N_p4;i++) {
		                  V3_p4_rot[i]= V3_prot4[i];
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
	            TVector3 V3_prot[N_3p],V3_prot_uncorr[N_3p],V3_prot_el_4pto3p[N_3p][N_2p],V3_el_prot[N_comb][N_2p];
	            double N_p_three=0,N_p1[N_3p]={0};
              double E_cal_43pto1p[N_3p];
              double p_miss_perp_43pto1p[N_3p];
              double P_43pto1p[3]={0};

	            for(int g = 0; g < Ncomb_4to3; g++){   //estimating the undetected 4p contribution to  3p
		             if(g==0) {
		                 V3_prot_uncorr[0]=V3_prot4[0]; V3_prot_uncorr[1]=V3_prot4[1]; V3_prot_uncorr[2]=V3_prot4[2];
		                 V3_prot[0]=V3_prot4_corr[0]; V3_prot[1]=V3_prot4_corr[1]; V3_prot[2]=V3_prot4_corr[2];
		             }
		             if(g==1){
		                 V3_prot_uncorr[0]=V3_prot4[0]; V3_prot_uncorr[1]=V3_prot4[1]; V3_prot_uncorr[2]=V3_prot4[3];
		                 V3_prot[0]=V3_prot4_corr[0]; V3_prot[1]=V3_prot4_corr[1]; V3_prot[2]=V3_prot4_corr[3];
		             }
		             if(g==2){
		                 V3_prot_uncorr[0]=V3_prot4[0]; V3_prot_uncorr[1]=V3_prot4[2]; V3_prot_uncorr[2]=V3_prot4[3];
		                 V3_prot[0]=V3_prot4_corr[0]; V3_prot[1]=V3_prot4_corr[2]; V3_prot[2]=V3_prot4_corr[3];
		             }
		             if(g==3){
		                 V3_prot_uncorr[0]=V3_prot4[1]; V3_prot_uncorr[1]=V3_prot4[2]; V3_prot_uncorr[2]=V3_prot4[3];
		                 V3_prot[0]=V3_prot4_corr[1]; V3_prot[1]=V3_prot4_corr[2]; V3_prot[2]=V3_prot4_corr[3];
		             }

		             rotation->prot3_rot_func( V3_prot, V3_prot_uncorr,V4_el,E_cal_4pto3p,p_miss_perp_4pto3p, P_4pto3p,N_p1,E_cal_43pto1p,p_miss_perp_43pto1p,&N_p_three);

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

		                      h1_E_tot_4pto3p->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h1_E_rec_4pto3p->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h1_Etot_p4321_bkgd09->Fill(E_cal_4pto3p[count][j],P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h2_Erec_pperp_4321p->Fill(p_miss_perp_4pto3p[count][j],E_rec,P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h1_E_tot_4pto3p_fracfeed->Fill((E_cal_4pto3p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h1_E_rec_4pto3p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h2_pperp_W->Fill(W_var,p_miss_perp_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h1_theta0->Fill((V4_beam.Vect()).Angle(V3_el_prot[count][j])*TMath::RadToDeg(),-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		                      h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);

		                      for (int n = 0; n < N_pperp; n++){
			                       for(int z = 0; z < N_Ecal; z++){
			                            if(E_cal_4pto3p[count][j] > Ecal_lowlim[z] && E_cal_4pto3p[count][j] < Ecal_uplim[z] && p_miss_perp_4pto3p[count][j]>pperp_cut[n])  {
                                     h1_Etot_p_bkgd_slice_Ecalcut4321[n][z]->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
                                  }
                             }
		                      }
		                      for(int i = 0; i < n_slice; i++)
			                    {
			                         if (p_miss_perp_4pto3p[count][j]<pperp_max[i] && p_miss_perp_4pto3p[count][j]>pperp_min[i]){
			                              h1_Etot_4pto3p_slice[i]->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
			                              h1_Erec_4pto3p_slice[i]->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
			                         }
			                    }
		                    } //end loop over N_2p
		                 } //end loop over N_comb : number of 2 proton combination out of 3 protons

		   //-----------------------------------------  4p to 3p->1p  -----------------------------------------------------------------------
		   	             for(int j = 0; j < N_3p; j++)    { //4p to 3p->1, looping through 1p out of 3p
                        //P_43pto1p doesnt have to be an array, one local variable here
		                     P_43pto1p[j]= N_p1[j]/N_p_three*(N_p4_p3[g]/N_p_four);
		                     h1_E_tot_43pto1p->Fill(E_cal_43pto1p[j], P_43pto1p[j]*1/Mott_cross_sec);
		                     h1_E_rec_43pto1p->Fill(E_rec,P_43pto1p[j]*1/Mott_cross_sec);
		                     h1_Etot_p431_bkgd09->Fill(E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
		                     h2_Erec_pperp_431p->Fill(p_miss_perp_43pto1p[j],E_rec,P_43pto1p[j]*1/Mott_cross_sec);
		                     h1_E_tot_43pto1p_fracfeed->Fill((E_cal_43pto1p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_43pto1p[j]*1/Mott_cross_sec);
		                     h1_E_rec_43pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_43pto1p[j]*1/Mott_cross_sec);
		                     h2_pperp_W->Fill(W_var,p_miss_perp_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
  		                   h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_43pto1p[j]*1/Mott_cross_sec);
		                     h2_Ecal_Eqe->Fill(E_rec,E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
		                     h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
		                     h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);

		                     for (int n = 0; n < N_pperp; n++){
		                         for(int z = 0; z < N_Ecal; z++){
			                            if(E_cal_43pto1p[j] > Ecal_lowlim[z] && E_cal_43pto1p[j] < Ecal_uplim[z] && p_miss_perp_43pto1p[j]>pperp_cut[n]) {
                                      h1_Etot_p_bkgd_slice_Ecalcut431[n][z]->Fill(E_cal_43pto1p[j], P_43pto1p[j]*1/Mott_cross_sec);
                                  }
		                         }
		                     }
                         for(int i = 0; i <n_slice; i++)
		                     {
			                        if (p_miss_perp_43pto1p[j]<pperp_max[i] && p_miss_perp_43pto1p[j]>pperp_min[i]){
			                             h1_Etot_43pto1p_slice[i]->Fill(E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
			                             h1_Erec_43pto1p_slice[i]->Fill(E_rec,P_43pto1p[j]*1/Mott_cross_sec);
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
		                      V3p2_uncorr[0]=V3_prot4[ind1];
		                      V3p2_uncorr[1]=V3_prot4[ind2];

		                      rotation->prot2_rot_func( V3p2,V3p2_uncorr, V4_el,E_cal_4pto2p,p_miss_perp_4pto2p,  P_4pto2p, &N_two);

		                      if( N_two!=0  && N_p_four!=0){
		                          for(int j = 0; j < N_2p; j++)  {  //looping through  1 proton combination out of 2 protons

	                               h1_E_tot_4pto2p->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h1_E_rec_4pto2p->Fill(E_rec, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h1_Etot_p421_bkgd09->Fill(E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h2_Erec_pperp_421p->Fill( p_miss_perp_4pto2p[j],E_rec,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h1_E_tot_4pto2p_fracfeed->Fill((E_cal_4pto2p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h1_E_rec_4pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h2_pperp_W->Fill(W_var,p_miss_perp_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3p2_uncorr[j])*TMath::RadToDeg(),P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
                                 h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
	                               h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
		                             h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);

	                               for (int n = 0; n < N_pperp; n++){
	                                  for(int z = 0; z < N_Ecal; z++){
	                                     if(E_cal_4pto2p[j] > Ecal_lowlim[z] && E_cal_4pto2p[j] < Ecal_uplim[z] && p_miss_perp_4pto2p[j]>pperp_cut[n]) {
                                          h1_Etot_p_bkgd_slice_Ecalcut421[n][z]->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
                                       }
	                                  }
	                               }
	                               for(int i = 0; i < n_slice; i++)
	                               {
	                                  if (p_miss_perp_4pto2p[j]<pperp_max[i] && p_miss_perp_4pto2p[j]>pperp_min[i]){
	                                     h1_Etot_4pto2p_slice[i]->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
		                                   h1_Erec_4pto2p_slice[i]->Fill(E_rec, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
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
		                  h1_E_tot_4pto1p->Fill(E_cal_p4[j], P_4pto1p[j]*1/Mott_cross_sec);
		                  h1_E_rec_4pto1p->Fill(E_rec,P_4pto1p[j]*1/Mott_cross_sec);
		                  h1_Etot_p41_bkgd09->Fill(E_cal_p4[j],P_4pto1p[j]*1/Mott_cross_sec);
		                  h2_Erec_pperp_41p->Fill(p_miss_perp_p4[j],E_rec, P_4pto1p[j]*1/Mott_cross_sec);
		                  h1_E_tot_4pto1p_fracfeed->Fill((E_cal_p4[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto1p[j]*1/Mott_cross_sec);
		                  h1_E_rec_4pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_4pto1p[j]*1/Mott_cross_sec);
		                  h2_pperp_W->Fill(W_var,p_miss_perp_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);
		                  h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot4[j])*TMath::RadToDeg(),-P_4pto1p[j]*1/Mott_cross_sec);
		                  h2_Ecal_Eqe->Fill(E_rec,E_cal_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);
		                  h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);
		                  h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);

		                  for (int n = 0; n < N_pperp; n++){
		                      for(int z = 0; z < N_Ecal; z++){
		                          if(E_cal_p4[j] > Ecal_lowlim[z] && E_cal_p4[j] < Ecal_uplim[z] && p_miss_perp_p4[j]>pperp_cut[n]) {
                                  h1_Etot_p_bkgd_slice_Ecalcut41[n][z]->Fill(E_cal_p4[j], P_4pto1p[j]*1/Mott_cross_sec);
                              }
		                      }
		                  }
		                  for(int i = 0; i < n_slice; i++)
		                  {
		                      if (p_miss_perp_p4[j]<pperp_max[i] && p_miss_perp_p4[j]>pperp_min[i]){
                              h1_Etot_4pto1p_slice[i]->Fill(E_cal_p4[j],P_4pto1p[j]*1/Mott_cross_sec);
			                        h1_Erec_4pto1p_slice[i]->Fill(E_rec,P_4pto1p[j]*1/Mott_cross_sec);
		                      }
		                  }
                  } //end loop over N_p4
	            } // end if N_p_four!=0
           }//no pion statment ends
	      }//4 proton vertex cut
	   }//4 proton requirement (num_p == 4)

     double P_undet=0;
     TVector3 V3_pi;

	   V3_q=(V4_beam-V4_el).Vect();
	   h1_E_rec->Fill(E_rec,1/Mott_cross_sec);

	   if (num_pi_phot==0){
       h1_E_rec_0pi->Fill(E_rec,1/Mott_cross_sec);
	     h1_E_rec_0pi_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],1/Mott_cross_sec);
     }

//----------------------------- e- ,1pi  -----------------------------------------

	   if (num_pi_phot==1) {

	      V3_pi.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
	      rotation->pi1_rot_func(V3_pi,q[ind_pi_phot[0]],&P_undet);

	      h1_E_rec_1pi_weight->Fill(E_rec,P_undet*1/Mott_cross_sec);
	      h1_E_rec_1pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_undet*1/Mott_cross_sec);

	      if(!ec_radstat_n[0])  h1_E_rec_1pi->Fill(E_rec,1/Mott_cross_sec);
	      if(ec_num_n==1)       h2_phot_e_angle_Erec->Fill(E_rec,V3_pi.Angle(V4_el.Vect())*TMath::RadToDeg());
	   }
//----------------------------- e- ,2pi  -----------------------------------------

     if (num_pi_phot==2) {

	      const int N_2pi=2;
	      TVector3 V3_2pi[N_2pi];
	      int q_pi2[N_2pi];
        bool radstat_pi2[N_2pi]={false};
	      double P_1pi[N_2pi]={0};
        double P_0pi=0;

        V3_2pi[0].SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
        V3_2pi[1].SetXYZ(p[ind_pi_phot[1]]*cx[ind_pi_phot[1]],p[ind_pi_phot[1]]*cy[ind_pi_phot[1]],p[ind_pi_phot[1]]*cz[ind_pi_phot[1]]);
        q_pi2[0]=q[ind_pi_phot[0]];
        q_pi2[1]=q[ind_pi_phot[1]];
        radstat_pi2[0]=ec_radstat_n[0];
        radstat_pi2[1]=ec_radstat_n[1];

        rotation->pi2_rot_func( V3_2pi, q_pi2,&P_0pi,P_1pi);

//----------------------------- e- ,2pi->0pi (-) -----------------------------------------
        h1_E_rec_2pi_weight->Fill(E_rec,(-P_0pi)*1/Mott_cross_sec);
        h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*1/Mott_cross_sec);
        h1_E_rec_20pi->Fill(E_rec,(P_0pi)*1/Mott_cross_sec);
//----------------------------- e- ,2pi->1pi->0pi (+)  -----------------------------------------
        for(int k = 0; k < N_2pi; k++){ //loop over two pions
          h1_E_rec_2pi_weight->Fill(E_rec,P_1pi[k]*1/Mott_cross_sec);
          h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[k]*1/Mott_cross_sec);
          h1_E_rec_21pi->Fill(E_rec,(P_1pi[k])*1/Mott_cross_sec);
        }
     } //end if for two pion events

//----------------------------- e- ,3pi  -----------------------------------------

     if (num_pi_phot==3) {

        const int N_3pi=3;
        const int N_2pi=2;
        TVector3 V3_3pi[N_3pi];
        int q_pi3[N_3pi]={0};
        bool radstat_pi3[N_3pi]={false};
        double P_0pi=0;
        double P_1pi[N_3pi]={0};
        double P_320pi[N_3pi]={0};
        double P_3210pi[N_3pi][N_2pi]={0};

        for (int h = 0; h < N_3pi; h++){ //loop over three pions

          V3_3pi[h].SetXYZ(p[ind_pi_phot[h]]*cx[ind_pi_phot[h]],p[ind_pi_phot[h]]*cy[ind_pi_phot[h]],p[ind_pi_phot[h]]*cz[ind_pi_phot[h]]);
          q_pi3[h]=q[ind_pi_phot[h]];
          radstat_pi3[h]=ec_radstat_n[h];
        }

        rotation->pi3_rot_func( V3_3pi, q_pi3,&P_0pi, P_1pi, P_320pi,P_3210pi);

 //---------------------------3pi->0pi----------------------------------------------
        h1_E_rec_3pi_weight->Fill(E_rec,(-P_0pi)*1/Mott_cross_sec);
        h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*1/Mott_cross_sec);
        h1_E_rec_30pi->Fill(E_rec,(P_0pi)*1/Mott_cross_sec);

        for(int h = 0; h < N_3pi; h++){ //loop over three pions

//---------------------------3pi->1pi->0pi----------------------------------------------
          h1_E_rec_3pi_weight->Fill(E_rec,P_1pi[h]*1/Mott_cross_sec);
          h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[h]*1/Mott_cross_sec);
          h1_E_rec_310pi->Fill(E_rec,(P_1pi[h])*1/Mott_cross_sec);
 //---------------------------3pi->2pi->0pi----------------------------------------------
          h1_E_rec_3pi_weight->Fill(E_rec,P_320pi[h]*1/Mott_cross_sec);
          h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_320pi[h]*1/Mott_cross_sec);
          h1_E_rec_320pi->Fill(E_rec,(P_320pi[h])*1/Mott_cross_sec);
//---------------------------3pi->2pi->1pi->0pi----------------------------------------------
          for(int g = 0; g < N_2pi; g++){ //loop over two pions
            h1_E_rec_3pi_weight->Fill(E_rec,(-P_3210pi[h][g])*1/Mott_cross_sec);
            h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_3210pi[h][g])*1/Mott_cross_sec);
            h1_E_rec_3210pi->Fill(E_rec,(P_3210pi[h][g])*1/Mott_cross_sec);
          }
        }//end of 3pi loop
     }//end of 3pi requirement

//----------------------------- e- ,4pi  -----------------------------------------
     if (num_pi_phot==4) {

       const int N_4pi=4;
       TVector3 V3_4pi[N_4pi];
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

       for (int h = 0; h < N_4pi;h++){ //loop over four pions

         V3_4pi[h].SetXYZ(p[ind_pi_phot[h]]*cx[ind_pi_phot[h]],p[ind_pi_phot[h]]*cy[ind_pi_phot[h]],p[ind_pi_phot[h]]*cz[ind_pi_phot[h]]);
         q_pi4[h]=q[ind_pi_phot[h]];
         radstat_pi4[h]=ec_radstat_n[h];
       }

       rotation->pi4_rot_func( V3_4pi, q_pi4,&P_0pi,&P_410pi,&P_420pi,&P_4210pi,&P_430pi,&P_4310pi,&P_4320pi,&P_43210pi);

 //---------------------------4pi->0pi----------------------------------------------
//why is it here not split like for 3pi case, sum over all weights is done here
       h1_E_rec_4pi_weight->Fill(E_rec,(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*1/Mott_cross_sec);
       h1_E_rec_4pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*1/Mott_cross_sec);
       h1_E_rec_40pi->Fill(E_rec,(P_0pi)*1/Mott_cross_sec);
//---------------------------4pi->1pi->0pi----------------------------------------------
       h1_E_rec_410pi->Fill(E_rec,(P_410pi)*1/Mott_cross_sec);
//---------------------------4pi->2pi->0pi----------------------------------------------
       h1_E_rec_420pi->Fill(E_rec,(P_420pi)*1/Mott_cross_sec);
//---------------------------4pi->2pi->1pi->0pi----------------------------------------------
       h1_E_rec_4210pi->Fill(E_rec,(P_4210pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->0pi----------------------------------------------
       h1_E_rec_430pi->Fill(E_rec,(P_430pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->1pi->0pi----------------------------------------------
       h1_E_rec_4310pi->Fill(E_rec,(P_4310pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->2pi->0pi----------------------------------------------
       h1_E_rec_4320pi->Fill(E_rec,(P_4320pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->2pi->1pi->0pi----------------------------------------------
       h1_E_rec_43210pi->Fill(E_rec,(P_43210pi)*1/Mott_cross_sec);

     }//end of 4 pi/photon requirement


   //------------------------------------------requiring there to be a proton -------------------------------------
    //Events with exactly one proton
     if( num_p==1)
     {

       ind_p = index_p[0];

       TLorentzVector V4_prot_uncorr(p[ind_p]*cx[ind_p],p[ind_p]*cy[ind_p],p[ind_p]*cz[ind_p],TMath::Sqrt(m_prot*m_prot+p[ind_p]*p[ind_p]));
      //Proton kinematics
       float prot_vert=vz[ind_p];
       double p_phi=TMath::ATan2(cy[ind_p],cx[ind_p])*TMath::RadToDeg();
       double p_phi_mod=p_phi+30;
       if (p_phi_mod<0)p_phi_mod=p_phi_mod+360;
       double	p_theta=TMath::ACos(cz[ind_p])*TMath::RadToDeg();
       double prot_vert_corr=prot_vert+vz_corr(vz_corr_func,p_phi_mod,p_theta);
       double prot_mom_corr;

       h1_el_prot_vertdiff->Fill(el_vert_corr-prot_vert_corr);

       if(ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr,prot_vert_corr) != -1)
       {
	        prot_mom_corr = ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr,prot_vert_corr);
       }
       else
	     {
          prot_mom_corr=p[ind_p];
       }

       h1_prot_mom->Fill(prot_mom_corr);
       h1_prot_mom_ratio->Fill(prot_mom_corr/p[ind_p]);

       TLorentzVector V4_prot_corr(prot_mom_corr*cx[ind_p],prot_mom_corr*cy[ind_p],prot_mom_corr*cz[ind_p],TMath::Sqrt(m_prot*m_prot+prot_mom_corr*prot_mom_corr));
       TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();
       TLorentzVector V4_prot_el_tot = V4_prot_corr + V4_el;
       double p_perp_tot=TMath::Sqrt(V4_prot_el_tot.Px()*V4_prot_el_tot.Px()+V4_prot_el_tot.Py()*V4_prot_el_tot.Py());
       double p_z_tot=V4_prot_el_tot.Pz();
       double p_tot=V4_prot_el_tot.Rho();
       double E_tot=V4_el.E()+V4_prot_corr.E()-m_prot+bind_en[ftarget];

       //Vertex cut was removed here since it is already in the event selection loop F.H 08/13/19

 //---------------------------------- 1p 2pi   ----------------------------------------------
       if (num_pi_phot==2) {

         const int N_2pi=2;
         TVector3 V3_2pi[N_2pi],V3_2pi_rot[N_2pi],V3_p_rot;
         bool pi2_stat[N_2pi]={false};
         int q_pi2[N_2pi];
         bool radstat_pi2[N_2pi]={false};
         double P_1p0pi=0;
         double P_1p1pi[N_2pi]={0};
         V3_q=(V4_beam-V4_el).Vect();

         V3_2pi[0].SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
         V3_2pi[1].SetXYZ(p[ind_pi_phot[1]]*cx[ind_pi_phot[1]],p[ind_pi_phot[1]]*cy[ind_pi_phot[1]],p[ind_pi_phot[1]]*cz[ind_pi_phot[1]]);
         q_pi2[0]=q[ind_pi_phot[0]];
         q_pi2[1]=q[ind_pi_phot[1]];
         radstat_pi2[0]=ec_radstat_n[0];
         radstat_pi2[1]=ec_radstat_n[1];

         rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi,q_pi2,&P_1p0pi,P_1p1pi);

 //---------------------------------- 1p 2pi->1p1pi   ----------------------------------------------

         for(int z = 0; z < N_2pi; z++){  //to consider 2 diff. 1pi states

           h1_E_tot_1p2pi->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
           h1_E_rec_1p2pi->Fill(E_rec,P_1p1pi[z]*1/Mott_cross_sec);
           h2_Erec_pperp_1p2pi_1p1pi->Fill(p_perp_tot,E_rec,P_1p1pi[z]*1/Mott_cross_sec);
           h1_Etot_bkgd09_1p2pi_1p1pi->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
           h1_E_tot_1p2pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*1/Mott_cross_sec);
           h1_E_rec_1p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p1pi[z]*1/Mott_cross_sec);
           h2_pperp_W->Fill(W_var,p_perp_tot,P_1p1pi[z]*1/Mott_cross_sec);
           h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p1pi[z]*1/Mott_cross_sec);
           h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p1pi[z]*1/Mott_cross_sec);
           h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p1pi[z]*1/Mott_cross_sec);
           h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p1pi[z]*1/Mott_cross_sec);

           for(int i = 0; i < n_slice; i++){
	            if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	               h1_Etot_bkgd_1p2pi[i]->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
	               h1_Erec_bkgd_1p2pi[i]->Fill(E_rec,P_1p1pi[z]*1/Mott_cross_sec);
	            }
           }
           for (int i = 0; i < N_pperp; i++){
	            for(int j = 0; j < N_Ecal; j++){
	               if(E_tot > Ecal_lowlim[j] && E_tot < Ecal_uplim[j] && p_perp_tot>pperp_cut[i])  {
                    h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[i][j]->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
                 }
	            }
           }
         } //end loop over N_2pi

 //---------------------------------- 1p 2pi->1p0pi   ----------------------------------------------

         h1_E_tot_1p2pi_1p0pi->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
         h1_E_rec_1p2pi_1p0pi->Fill(E_rec,P_1p0pi*1/Mott_cross_sec);
         h2_Erec_pperp_1p2pi_1p0pi->Fill(p_perp_tot,E_rec,P_1p0pi*1/Mott_cross_sec);
         h1_Etot_bkgd09_1p2pi_1p0pi->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
         h1_E_tot_1p2pi_1p0pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p0pi*1/Mott_cross_sec);
         h1_E_rec_1p2pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p0pi*1/Mott_cross_sec);
         h2_pperp_W->Fill(W_var,p_perp_tot,-P_1p0pi*1/Mott_cross_sec);
         h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-P_1p0pi*1/Mott_cross_sec);
         h2_Ecal_Eqe->Fill(E_rec,E_tot,-P_1p0pi*1/Mott_cross_sec);
         h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-P_1p0pi*1/Mott_cross_sec);
         h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-P_1p0pi*1/Mott_cross_sec);

         for(int i = 0; i < n_slice; i++){
           if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	            h1_Etot_bkgd_1p2pi_1p0pi[i]->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
	            h1_Erec_bkgd_1p2pi_1p0pi[i]->Fill(E_rec,P_1p0pi*1/Mott_cross_sec);
           }
         }
         for (int i = 0;i < N_pperp; i++){
           for(int j = 0; j < N_Ecal; j++){
	            if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i]) {
                  h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[i][j]->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
              }
           }
         }
       }//1p 2pi statetment ends

 //---------------------------------- 1p 3pi   ----------------------------------------------
       if (num_pi_phot==3) {

         const int N_3pi=3;
         TVector3 V3_3pi[N_3pi],V3_3pi_rot[N_3pi],V3_p_rot;
         bool pi3_stat[N_3pi]={false};
         int q_pi3[N_3pi];
         bool radstat_pi3[N_3pi]={false};
         double N_3pi_p=0;
         double N_1pi_1p[N_3pi]={0};
         double N_pi=0;
         double N_nopi=0;
         double N_0pi_1p=0;
         double P_1p3pi=0;
         V3_q=(V4_beam-V4_el).Vect();

         V3_3pi[0].SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
         V3_3pi[1].SetXYZ(p[ind_pi_phot[1]]*cx[ind_pi_phot[1]],p[ind_pi_phot[1]]*cy[ind_pi_phot[1]],p[ind_pi_phot[1]]*cz[ind_pi_phot[1]]);
         V3_3pi[2].SetXYZ(p[ind_pi_phot[2]]*cx[ind_pi_phot[2]],p[ind_pi_phot[2]]*cy[ind_pi_phot[2]],p[ind_pi_phot[2]]*cz[ind_pi_phot[2]]);
         q_pi3[0]=q[ind_pi_phot[0]];
         q_pi3[1]=q[ind_pi_phot[1]];
         q_pi3[2]=q[ind_pi_phot[2]];
         radstat_pi3[0]=ec_radstat_n[0];
         radstat_pi3[1]=ec_radstat_n[1];
         radstat_pi3[2]=ec_radstat_n[2];

         rotation->prot1_pi3_rot_func(V3_prot_uncorr,V3_3pi,q_pi3,&P_1p3pi);

 //---------------------------------- 1p 3pi->1p 0pi  total ?? F.H. 08/13/19 check logic here compared to 1p 2pi case ----------------------------------------------
         h1_E_tot_1p3pi->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
         h1_E_rec_1p3pi->Fill(E_rec,P_1p3pi*1/Mott_cross_sec);
         h2_Erec_pperp_1p3pi->Fill(p_perp_tot,E_rec,P_1p3pi*1/Mott_cross_sec);
         h1_Etot_bkgd09_1p3pi->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
         h1_E_tot_1p3pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi*1/Mott_cross_sec);
         h1_E_rec_1p3pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p3pi*1/Mott_cross_sec);
         h2_pperp_W->Fill(W_var,p_perp_tot,P_1p3pi*1/Mott_cross_sec);
         h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p3pi*1/Mott_cross_sec);
         h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p3pi*1/Mott_cross_sec);
         h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p3pi*1/Mott_cross_sec);
         h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p3pi*1/Mott_cross_sec);

         for(int i = 0; i < n_slice; i++){
           if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
             h1_Etot_bkgd_1p3pi[i]->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
             h1_Erec_bkgd_1p3pi[i]->Fill(E_rec,P_1p3pi*1/Mott_cross_sec);
           }
         }
         for (int i = 0; i < N_pperp; i++){
           for(int j = 0; j < N_Ecal; j++){
             if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i])  {
                h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[i][j]->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
             }
           }
         }
       }//end of 1p 3pi requirement

  //---------------------------------- 1p 1pi   ----------------------------------------------
       if (num_pi_phot==1) {
         double N_piphot_det,N_piphot_undet;
         TVector3 V3_pi_phot;
         V3_q=(V4_beam-V4_el).Vect();
         V3_pi_phot.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);	N_piphot_det=N_piphot_undet=0;

         rotation->prot1_pi1_rot_func(V3_prot_uncorr,V3_pi_phot, q[ind_pi_phot[0]], &N_piphot_det,&N_piphot_undet);

		     if(N_piphot_det!=0){

	          h1_E_rec_undetfactor->Fill(E_rec,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
	          h1_E_tot_undetfactor->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
	          h1_E_tot_undetfactor_pipl->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
	          h1_E_tot_undetfactor09->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
	          h2_Erec_pperp_1p1pi->Fill(p_perp_tot,E_rec,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		        h1_E_rec_undetfactor_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		        h1_E_tot_undetfactor_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		        h2_pperp_W->Fill(W_var,p_perp_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		        h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		        h2_Ecal_Eqe->Fill(E_rec,E_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		        h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
	          h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);

            for(int i = 0; i < n_slice; i++)
            {
              if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
                h1_Etot_bkgd_pipl_pimi_fact[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
                h1_Erec_bkgd_pipl_pimi_new_fact[i]->Fill(E_rec,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
                h1_Etot_bkgd_pipl_pimi_fact_pipl[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
                h1_Etot_Npi1[i]->Fill(E_tot,1/Mott_cross_sec);
                h1_Erec_Npi1[i]->Fill(E_rec,1/Mott_cross_sec);
              }
            }
            for (int i = 0; i < N_pperp; i++){
              for(int j = 0; j < N_Ecal; j++){
                if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i]) {
                    h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
                }
              }
		        }
		     } //end of N_piphot_det!=0

	       h1_E_rec_cutpi1_piplpimi->Fill(E_rec,1/Mott_cross_sec);
	       h1_E_tot_cutpi1_piplpimi->Fill(E_tot,1/Mott_cross_sec);

       }//end of 1p 1pi requirement

  //---------------------------------- 1p 0pi   ----------------------------------------------
       if (num_pi_phot==0){

          h2_Erec_pperp_newcut2->Fill(p_perp_tot,E_rec,1/Mott_cross_sec);
	        h1_E_rec_cut2_new->Fill(E_rec,1/Mott_cross_sec);
	        h1_E_tot_cut2->Fill(E_tot,1/Mott_cross_sec);
	        h1_E_tot_cut2_09->Fill(E_tot,1/Mott_cross_sec);
	        h1_E_tot_cut2_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],1/Mott_cross_sec);
	        h1_E_rec_cut2_new_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],1/Mott_cross_sec);
	        h2_pperp_W->Fill(W_var,p_perp_tot,1/Mott_cross_sec);
	        h1_theta0->Fill((V4_beam.Vect()).Angle(V4_prot_el_tot.Vect()) *TMath::RadToDeg(),1/Mott_cross_sec);
	        h2_Ecal_Eqe->Fill(E_rec,E_tot,1/Mott_cross_sec);
	        h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,1/Mott_cross_sec);
	        h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,1/Mott_cross_sec);

          for(int i = 0; i < n_slice; i++)
          {
            if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
              h1_Etot_Npi0[i]->Fill(E_tot,1/Mott_cross_sec);
              h1_Erec_Npi0_new[i]->Fill(E_rec,1/Mott_cross_sec);
            }
          }
          for (int i = 0; i < N_pperp; i++){
	           for(int j = 0; j < N_Ecal; j++){
		             if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i])  {
                     h1_Etot_Npi0_Ecalcut[i][j]->Fill(E_tot,1/Mott_cross_sec);
                 }
             }
	        }
          if (p_perp_tot < 0.2){
            h1_E_rec_cut005_newcut3->Fill(E_rec,1/Mott_cross_sec);
	          h2_Erec_pperp_cut3->Fill(p_perp_tot,E_rec,1/Mott_cross_sec);
	        }
       }//num pi=0

       h1_Etot->Fill(E_tot,1/Mott_cross_sec);
       h1_E_rec_1prot->Fill(E_rec,1/Mott_cross_sec);
       h1_E_tot_1prot->Fill(E_tot,1/Mott_cross_sec);
       h2_Erec_pperp->Fill(p_perp_tot,E_rec,1/Mott_cross_sec);

     } // 1proton ends


  } //end of event loop (jentry)

  gStyle->SetOptFit(1);

  for(int i = 0; i <= n_slice-1; i++)
  {
  //------------------------------------using the ratio of the pi- to pi+  --------------------------------------

	   h1_Etot_piplpimi_subtruct_fact[i] = (TH1F*)  h1_Etot_Npi0[i]->Clone(Form("h1_Etot_piplpimi_subtruct_fact_%d",i+1));
	   h1_Etot_piplpimi_subtruct_fact[i]->Add(h1_Etot_bkgd_pipl_pimi_fact[i],-1);
	   h1_Erec_piplpimi_subtruct_new_fact[i]=(TH1F*)  h1_Erec_Npi0_new[i]->Clone(Form("h1_Erec_piplpimi_subtruct_new_fact_%d",i+1));
	   h1_Erec_piplpimi_subtruct_new_fact[i]->Add(h1_Erec_bkgd_pipl_pimi_new_fact[i],-1);

 //------------------------------------subtracting 2p contribution from 1p events  --------------------------------------

	   h1_Etot_p_bkgd_slice_sub[i]=(TH1F*) h1_Etot_piplpimi_subtruct_fact[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub_%d",i+1));
	   h1_Etot_p_bkgd_slice_sub[i]->Add(h1_Etot_p_bkgd_slice[i],-1);
	   h1_Erec_p_bkgd_slice_sub[i]=(TH1F*) h1_Erec_piplpimi_subtruct_new_fact[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub_%d",i+1));
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

  for (int i = 0; i < N_pperp; i++) {
    for(int j = 0;j < N_Ecal; j++) {

      //------------------------------------using the ratio of the pi- to pi+  --------------------------------------
      h1_Etot_piplpimi_subtruct_fact_Ecalcut[i][j]=(TH1F*)  h1_Etot_Npi0_Ecalcut[i][j]->Clone(Form("h1_Etot_piplpimi_subtruct_fact_Ecalcut_%d_%d",i+1,i+j));
      h1_Etot_piplpimi_subtruct_fact_Ecalcut[i][j]->Add(h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j],-1);

      //------------------------------------subtracting 2p contribution from 1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut[i][j]=(TH1F*) h1_Etot_piplpimi_subtruct_fact_Ecalcut[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut[i][j],-1);

 //------------------------------------subtracting 3p to  2p->1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut32[i][j]=(TH1F*) h1_Etot_p_bkgd_slice_sub_Ecalcut[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut32_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut32[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut321[i][j]);

 //------------------------------------subtracting 3p to  1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut31[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut32[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut31_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut31[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut31[i][j],-1);

 //------------------------------------subtracting 4p to  3p->2p->1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut43[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut31[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut43_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut43[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut4321[i][j],-1);

 //------------------------------------subtracting 4p to  3p->1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut431[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut43[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut431_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut431[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut431[i][j]);

//------------------------------------subtracting 4p to  2p->1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut42[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut431[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut42_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut42[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut421[i][j]);

//------------------------------------subtracting 4p to  1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut41[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut42[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut41_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut41[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut41[i][j],-1);

//------------------------------------undetected 1p 2pi ->1p1pi ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p1pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut41[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p1pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p1pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[i][j]);

//------------------------------------undetected 1p 2pi-> 1p 0pi  ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p0pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p1pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p0pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p0pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[i][j],-1);

//------------------------------------undetected 1p 3pi-> 1p 0pi  ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_1p3pi_1p0pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p0pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_1p3pi_1p0pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_1p3pi_1p0pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[i][j]);

//------------------------------------undetected 2p 2pi ->1p 0pi  ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p2pi_1p0pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut_1p3pi_1p0pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_2p2pi_1p0pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p2pi_1p0pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[i][j]);

//------------------------------------undetected 3p 1pi ->1p 0pi  ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_3p1pi_1p0pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut_2p2pi_1p0pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_3p1pi_1p0pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_3p1pi_1p0pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[i][j]);

//------------------------------------undetected 2p 1pi ->2p 0pi  ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_2p0pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut_3p1pi_1p0pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_2p0pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_2p0pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[i][j]);

//------------------------------------undetected 2p 1pi ->1p 1pi  ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p1pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_2p0pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p1pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p1pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[i][j]);

//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p0pi[i][j]=(TH1F*)  h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p1pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p0pi_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p0pi[i][j]->Add(h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[i][j],-1);

    } //end of loop over N_Ecal
  } // end of loop over N_pperp

  //------------------------------------using the ratio of the pi- to pi+  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_factor =(TH1F*)  h1_E_rec_cut2_new->Clone("h_Erec_subtruct_piplpimi_factor");
  h_Erec_subtruct_piplpimi_factor->Add(h1_E_rec_undetfactor,-1);

  TH1F *h_Etot_subtruct_piplpimi_factor=(TH1F*)  h1_E_tot_cut2->Clone("h_Etot_subtruct_piplpimi_factor");
  h_Etot_subtruct_piplpimi_factor->Add(h1_E_tot_undetfactor,-1);

  TH1F *h_Etot_subtruct_piplpimi_factor09=(TH1F*)  h1_E_tot_cut2_09->Clone("h_Etot_subtruct_piplpimi_factor09");
  h_Etot_subtruct_piplpimi_factor09->Add(h1_E_tot_undetfactor09,-1);

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

  TH1F *h_Etot_subtruct_piplpimi_prot09=(TH1F*)  h_Etot_subtruct_piplpimi_factor09->Clone("h_Etot_subtruct_piplpimi_prot09");
  h_Etot_subtruct_piplpimi_prot09->Add(h1_E_tot_p_bkgd09,-1);

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

  TH1F *h_Etot_subtruct_piplpimi_32prot09=(TH1F*)  h_Etot_subtruct_piplpimi_prot09->Clone("h_Etot_subtruct_piplpimi_32prot09");
  h_Etot_subtruct_piplpimi_32prot09->Add(h1_Etot_p321_bkgd09);

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

  TH1F *h_Etot_subtruct_piplpimi_31prot09=(TH1F*)  h_Etot_subtruct_piplpimi_32prot09->Clone("h_Etot_subtruct_piplpimi_31prot09");
  h_Etot_subtruct_piplpimi_31prot09->Add(h1_Etot_p31_bkgd09,-1);

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

  TH1F *h_Etot_subtruct_piplpimi_43prot09=(TH1F*)  h_Etot_subtruct_piplpimi_31prot09->Clone("h_Etot_subtruct_piplpimi_43prot09");
  h_Etot_subtruct_piplpimi_43prot09->Add(h1_Etot_p4321_bkgd09,-1);

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

  TH1F *h_Etot_subtruct_piplpimi_431prot09=(TH1F*)  h_Etot_subtruct_piplpimi_43prot09->Clone("h_Etot_subtruct_piplpimi_431prot09");
  h_Etot_subtruct_piplpimi_431prot09->Add(h1_Etot_p431_bkgd09);

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

  TH1F *h_Etot_subtruct_piplpimi_42prot09=(TH1F*)  h_Etot_subtruct_piplpimi_431prot09->Clone("h_Etot_subtruct_piplpimi_42prot09");
  h_Etot_subtruct_piplpimi_42prot09->Add(h1_Etot_p421_bkgd09);

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

  TH1F *h_Etot_subtruct_piplpimi_41prot09=(TH1F*)  h_Etot_subtruct_piplpimi_42prot09->Clone("h_Etot_subtruct_piplpimi_41prot09");
  h_Etot_subtruct_piplpimi_41prot09->Add(h1_Etot_p41_bkgd09,-1);

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

  TH1F *h_Etot_subtruct_piplpimi09_1p2pi_1p1pi=(TH1F*)  h_Etot_subtruct_piplpimi_41prot09->Clone("h_Etot_subtruct_piplpimi09_1p2pi_1p1pi");
  h_Etot_subtruct_piplpimi09_1p2pi_1p1pi->Add(h1_Etot_bkgd09_1p2pi_1p1pi);

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

  TH1F *h_Etot_subtruct_piplpimi09_1p2pi_1p0pi=(TH1F*)  h_Etot_subtruct_piplpimi09_1p2pi_1p1pi->Clone("h_Etot_subtruct_piplpimi09_1p2pi_1p0pi");
  h_Etot_subtruct_piplpimi09_1p2pi_1p0pi->Add(h1_Etot_bkgd09_1p2pi_1p0pi,-1);

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

  TH1F *h_Etot_subtruct_piplpimi09_1p3pi=(TH1F*)  h_Etot_subtruct_piplpimi09_1p2pi_1p0pi->Clone("h_Etot_subtruct_piplpimi09_1p3pi");
  h_Etot_subtruct_piplpimi09_1p3pi->Add(h1_Etot_bkgd09_1p3pi);

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

  TH1F *h_Etot_subtruct_piplpimi09_2p2pi=(TH1F*)  h_Etot_subtruct_piplpimi09_1p3pi->Clone("h_Etot_subtruct_piplpimi09_2p2pi");
  h_Etot_subtruct_piplpimi09_2p2pi->Add(h1_Etot_bkgd09_2p2pi);

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

  TH1F *h_Etot_subtruct_piplpimi09_3p1pi=(TH1F*)  h_Etot_subtruct_piplpimi09_2p2pi->Clone("h_Etot_subtruct_piplpimi09_3p1pi");
  h_Etot_subtruct_piplpimi09_3p1pi->Add(h1_Etot_bkgd09_3p1pi);

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

  TH1F *h_Etot_subtruct_piplpimi09_2p1pi_2p0pi=(TH1F*)  h_Etot_subtruct_piplpimi09_3p1pi->Clone("h_Etot_subtruct_piplpimi09_2p1pi_2p0pi");
  h_Etot_subtruct_piplpimi09_2p1pi_2p0pi->Add(h1_Etot_bkgd09_2p1pi_2p0pi);

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

  TH1F *h_Etot_subtruct_piplpimi09_2p1pi_1p1pi=(TH1F*)  h_Etot_subtruct_piplpimi09_2p1pi_2p0pi->Clone("h_Etot_subtruct_piplpimi09_2p1pi_1p1pi");
  h_Etot_subtruct_piplpimi09_2p1pi_1p1pi->Add(h1_Etot_bkgd09_2p1pi_1p1pi);

  TH2F *h2_Erec_pperp_sub_2p1pi_1p1pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_2p0pi->Clone("h2_Erec_pperp_sub_2p1pi_1p1pi");
  h2_Erec_pperp_sub_2p1pi_1p1pi->Add(h2_Erec_pperp_2p1pi_1p1pi);

  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
  h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_rec_2p1pi_1p1pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_tot_2p1pi_1p1pi_fracfeed);

//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

// apapadop
//  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi");
  TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Clone("eRecoEnergy_slice_0");
  h_Erec_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_rec_2p1pi_1p0pi,-1);

// apapadop
  //  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi");
  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("epRecoEnergy_slice_0");
  h_Etot_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_tot_2p1pi_1p0pi,-1);

  TH1F *h_Etot_subtruct_piplpimi09_2p1pi_1p0pi=(TH1F*)  h_Etot_subtruct_piplpimi09_2p1pi_1p1pi->Clone("h_Etot_subtruct_piplpimi09_2p1pi_1p0pi");
  h_Etot_subtruct_piplpimi09_2p1pi_1p0pi->Add(h1_Etot_bkgd09_2p1pi_1p0pi,-1);

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

  cout<<" ---------------------------2p subtracted -----------------------"<<endl;

  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){
      h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio[i][j]=   (TH1F *)  h1_Etot_p_bkgd_slice_sub_Ecalcut[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio[i][j]->Divide(h_Etot_subtruct_piplpimi_prot09);
      cout<<200+i*200<<"   "<<j+1<<"   "<<h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio[i][j]->GetBinContent(1)<<"  Error  "<<h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio[i][j]->GetBinError(1)<<endl;
    }
  }

  cout<<" ---------------------------3p subtracted -----------------------"<<endl;

  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){
      h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio3p[i][j]=   (TH1F *)   h1_Etot_p_bkgd_slice_sub_Ecalcut31[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio3p_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio3p[i][j]->Divide(h_Etot_subtruct_piplpimi_31prot09);
      cout<<200+i*200<<"   "<<j+1<<"   "<<h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio3p[i][j]->GetBinContent(1)<<"  Error  "<<h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio3p[i][j]->GetBinError(1)<<endl;
    }
  }

  cout<<" ---------------------------final subtracted -----------------------"<<endl;
  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){
      h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio4p[i][j]=   (TH1F *)   h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p0pi[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio4p_%d_%d",i+1,j+1));
      h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio4p[i][j]->Divide(h_Etot_subtruct_piplpimi09_2p1pi_1p0pi);
      cout<<200+i*200<<"   "<<j+1<<"   "<<h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio4p[i][j]->GetBinContent(1)<<"  Error  "<<h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio4p[i][j]->GetBinError(1)<<endl;


    }
  }


  fsub_pipl->Write();
  fsum_pipl->Write();
  fsub_pimi->Write();
  fsum_pimi->Write();
  fsum_e->Write();
  fsub_e->Write();
  fsum_prot->Write();
  fsub_prot->Write();
  gDirectory->Write("hist_Files", TObject::kOverwrite);
  // skim_tree->AutoSave();


  std::cout << "CounterEvents = " << CounterEvents << std::endl;

  delete[]  pperp_cut;
  delete[] Ecal_lowlim;
  delete[] Ecal_uplim;
  delete pipl_deltat_sig;delete pipl_deltat_mean;delete pimi_deltat_sig;delete pimi_deltat_mean;delete fsum_pimi;delete fsub_pimi;delete fsum_pipl;delete fsub_pipl;delete prot_deltat_sig;delete prot_deltat_mean;delete fsum_prot;delete fsub_prot;delete el_Epratio_sig;delete el_Epratio_mean;delete fsum_e;delete fsub_e;

}

//End Loop function



double vz_corr(TF1 *vz_corr_func, double phi,double theta)            //correction function for vertex , takes the arguments in Deg.
{
  //  return (0.2)*cos((phi+47.9)*TMath::DegToRad())/tan(theta*TMath::DegToRad());
   // vertex correction function obtained for the empty runs 18393 and 18394, works fine for 3He runs at 2.261[GeV] beam energy

  return (-(vz_corr_func->GetParameter(1)))*cos((phi-(vz_corr_func->GetParameter(2)))*TMath::DegToRad())/tan(theta*TMath::DegToRad());
  //vertex correction function for 4He runs at 2.261[GeV] beam energy obtained for the empty run18283

}

TVector3 FindUVW(TVector3 xyz)
{
  // get the U V W distance to EC edge for the purpose of geometry cut
  // ported from Stepan's function ec_xyz_duvw. the input is lab coordinates
  // of the EC hit.
  Float_t x = xyz.X(); Float_t y = xyz.Y(); Float_t z = xyz.Z();
  Float_t xi,yi,zi,u,v,w;
  Float_t ec_the = 0.4363323;
  Float_t ylow = -182.974; Float_t yhi = 189.956;
  Float_t tgrho=1.95325; Float_t sinrho=0.8901256; Float_t cosrho=0.455715;
  Float_t phi=xyz.Phi()*180./TMath::Pi(); if(phi<-30) phi+=360;
  Int_t ec_sect = (phi+30)/60.; if(ec_sect<0)ec_sect=0; if(ec_sect>5)ec_sect=5;
  Float_t ec_phi = ec_sect*TMath::Pi()/3.;
  xi = -x*sin(ec_phi) + y*cos(ec_phi);
  yi = x*cos(ec_the)*cos(ec_phi) + y*cos(ec_the)*sin(ec_phi) - z*sin(ec_the);
  zi = x*sin(ec_the)*cos(ec_phi) + y*sin(ec_the)*sin(ec_phi) + z*cos(ec_the);
  zi -= 510.32;
  u = (yi-ylow)/sinrho;
  v = (yhi-ylow)/tgrho - xi + (yhi-yi)/tgrho;
  w = ((yhi-ylow)/tgrho + xi + (yhi-yi)/tgrho)/2./cosrho;
  TVector3 uvw(u,v,w);
  return uvw;
}

Bool_t CutUVW(TVector3 ecxyz)
{
  // Cut the edges of EC according to UVW distance threshold defined by par_EcUVW array.
  // If it passes the cut, return true, if not return false
  TVector3 ecuvw = FindUVW(ecxyz);
  Float_t phi=ecxyz.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
  Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
  return (ecuvw.X()>par_EcUVW[sector][0] && ecuvw.Y()<par_EcUVW[sector][1] && ecuvw.Z()<par_EcUVW[sector][2]);
}


Float_t ProtonMomCorrection_He3_4Cell(std::string atarget, TLorentzVector V4Pr, Float_t vertex_p ){

  // Low energy proton momentum correction function
  // to be used with He3 target (4th target cell) (RUN # 18338-18438)
  // Input: Proton momentum 4 vector, and Z coord of proton vertex.
  // Returns the corrected MAGNITUDE of the proton momentum,

  Float_t  up_parm[6]   = {2.001,  -14.94,  47.2,   -77.59,  65.73,  -22.85};
  Float_t  down_parm[6] = {1.4165, -13.004, 48.897, -92.443, 86.984, -32.424};

  Float_t  proton_p     = V4Pr.Vect().Mag();
  Float_t  thetta_p     = V4Pr.Vect().Theta()*57.3;
  Float_t  polinom_up   = (((((up_parm[5]*proton_p+up_parm[4])*proton_p+up_parm[3])
			  *proton_p+up_parm[2])*proton_p+up_parm[1])*proton_p+up_parm[0]);

  Float_t polinom_down = (((((down_parm[5]*proton_p+down_parm[4])*proton_p+down_parm[3])
			  *proton_p+down_parm[2])*proton_p+down_parm[1])*proton_p+down_parm[0]);


  if(polinom_up<0.  ) polinom_up   = 0;
  if(polinom_down<0.) polinom_down = 0;

  Float_t  p_corr_up   = proton_p + proton_p*polinom_up;
  Float_t  p_corr_down = proton_p + proton_p*polinom_down;

  p_corr_down=p_corr_down*4./3-proton_p/3;//artificial cut to match with Bins distribution
  p_corr_up=p_corr_down;//artificial cut to match with Bins distribution

  if((thetta_p>=70.)) return p_corr_up;

  if((thetta_p < 30.)||(vertex_p>=(1/20.*thetta_p-5/2))||
     (thetta_p<=(-200*proton_p+86))){
    return p_corr_down;
  }



  if((thetta_p<=70.)&&(thetta_p>=30)&&
     (thetta_p>(20*vertex_p+50))){
    return p_corr_up;
  }else if(proton_p<0.57){
    return p_corr_down;
  } else { return p_corr_up;}

  return -1.;
}
