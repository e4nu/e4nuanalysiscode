
#define E2A_EP_C

#include "e2a_ep_neutrino6_united4_radphot.h"
#include "FilterData.h"
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

using namespace std;
using namespace Constants;

void SetFiducialCutParameters(std::string beam_en); // Load Fidicual Parameters for 1.1 and 4.4 GeV from file
//void SetMomCorrParameters();


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

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------

void FilterData::Loop()
{

  target_name = ftarget;   //std string for target name

// -----------------------------------------------------------------------------------------------------

// apapadop

int TargetPdgCode, TargetZ, TargetA;

if (ftarget=="3He") { TargetPdgCode = 1000020030; TargetZ = 2; TargetA = 3; }
if (ftarget=="4He") { TargetPdgCode = 1000020040; TargetZ = 2; TargetA = 4; }
if (ftarget=="C12") { TargetPdgCode = 1000060120; TargetZ = 6; TargetA = 12; }
if (ftarget=="56Fe") { TargetPdgCode = 1000260560; TargetZ = 26; TargetA = 56; }

// -----------------------------------------------------------------------------------------------------

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
  TFile *file_out = new TFile(Form("genie_filtered_data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");

// __________________________________________________________________________________________________________________________________________________________________________________________________

	// apapadop
	TTree* mytree = new TTree("gst","gst");

	Int_t           genie_iev;
	Int_t           genie_neu;
	Int_t           genie_fspl;
	Int_t           genie_tgt;
	Int_t           genie_Z;
	Int_t           genie_A;
	Int_t           genie_hitnuc;
	Int_t           genie_hitqrk;
	Int_t           genie_resid;
	Bool_t          genie_sea;
	Bool_t          genie_qel;
	Bool_t          genie_mec;
	Bool_t          genie_res;
	Bool_t          genie_dis;
	Bool_t          genie_coh;
	Bool_t          genie_dfr;
	Bool_t          genie_imd;
	Bool_t          genie_imdanh;
	Bool_t          genie_singlek;
	Bool_t          genie_nuel;
	Bool_t          genie_em;
	Bool_t          genie_CC;
	Bool_t          genie_nc;
	Bool_t          genie_charm;
	Int_t           genie_neut_code;
	Int_t           genie_nuance_code;
	Double_t        genie_wght;
	Double_t        genie_xs;
	Double_t        genie_ys;
	Double_t        genie_ts;
	Double_t        genie_Q2s;
	Double_t        genie_Ws;
	Double_t        genie_x;
	Double_t        genie_y;
	Double_t        genie_t;
	Double_t        genie_Q2reco;
	Double_t        genie_W;
	Double_t        genie_EvRF;
	Double_t        genie_Ev;
	Double_t        genie_pxv;
	Double_t        genie_pyv;
	Double_t        genie_pzv;
	Double_t        genie_En;
	Double_t        genie_pxn;
	Double_t        genie_pyn;
	Double_t        genie_pzn;
	Double_t        genie_El;
	Double_t        genie_pxl;
	Double_t        genie_pyl;
	Double_t        genie_pzl;
	Double_t        genie_pl;
	Double_t        genie_cthl;
	Int_t           genie_nfp;
	Int_t           genie_nfn;
	Int_t           genie_nfpip;
	Int_t           genie_nfpim;
	Int_t           genie_nfpi0;
	Int_t           genie_nfkp;
	Int_t           genie_nfkm;
	Int_t           genie_nfk0;
	Int_t           genie_nfem;
	Int_t           genie_nfother;
	Int_t           genie_nip;
	Int_t           genie_nin;
	Int_t           genie_nipip;
	Int_t           genie_nipim;
	Int_t           genie_nipi0;
	Int_t           genie_nikp;
	Int_t           genie_nikm;
	Int_t           genie_nik0;
	Int_t           genie_niem;
	Int_t           genie_niother;
	Int_t           genie_ni;
	int InitialStateParticles = 2;
	Int_t           genie_pdgi[InitialStateParticles];   //[ni]
	Int_t           genie_resc[InitialStateParticles];   //[ni]
	Double_t        genie_Ei[InitialStateParticles];   //[ni]
	Double_t        genie_pxi[InitialStateParticles];   //[ni]
	Double_t        genie_pyi[InitialStateParticles];   //[ni]
	Double_t        genie_pzi[InitialStateParticles];   //[ni]
	Int_t           genie_nf;
	const int FinalStateParticles = 120;
	Int_t           genie_pdgf[FinalStateParticles];   //[nf]
	Double_t        genie_Ef[FinalStateParticles];   //[nf]
	Double_t        genie_pxf[FinalStateParticles];   //[nf]
	Double_t        genie_pyf[FinalStateParticles];   //[nf]
	Double_t        genie_pzf[FinalStateParticles];   //[nf]
	Double_t        genie_pf[FinalStateParticles];   //[nf]
	Double_t        genie_cthf[FinalStateParticles];   //[nf]
	Double_t        genie_vtxx;
	Double_t        genie_vtxy;
	Double_t        genie_vtxz;
	Double_t        genie_vtxt;
	Double_t        genie_sumKEf;
	Double_t        genie_calresp0;

	mytree->Branch("iev", &genie_iev, "iev/I");
	mytree->Branch("neu", &genie_neu, "neu/I");
	mytree->Branch("fspl", &genie_fspl, "fspl/I");
	mytree->Branch("tgt", &genie_tgt, "tgt/I");
	mytree->Branch("Z", &genie_Z, "Z/I");
	mytree->Branch("A", &genie_A, "A/I");
	mytree->Branch("hitnuc", &genie_hitnuc, "hitnuc/I");
	mytree->Branch("hitqrk", &genie_hitqrk, "hitqrk/I");
	mytree->Branch("resid", &genie_resid, "resid/I");
	mytree->Branch("sea", &genie_sea, "sea/O");
	mytree->Branch("qel", &genie_qel, "qel/O");
	mytree->Branch("mec", &genie_mec, "mec/O");
	mytree->Branch("res", &genie_res, "res/O");
	mytree->Branch("dis", &genie_dis, "dis/O");
	mytree->Branch("coh", &genie_coh, "coh/O");
	mytree->Branch("dfr", &genie_dfr, "dfr/O");
	mytree->Branch("imd", &genie_imd, "imd/O");
	mytree->Branch("imdanh", &genie_imdanh, "imdanh/O");
	mytree->Branch("singlek", &genie_singlek, "singlek/O");
	mytree->Branch("nuel", &genie_nuel, "nuel/O");
	mytree->Branch("em", &genie_em, "em/O");
	mytree->Branch("cc", &genie_CC, "cc/O");
	mytree->Branch("nc", &genie_nc, "nc/O");
	mytree->Branch("charm", &genie_charm, "charm/O");
	mytree->Branch("neut_code", &genie_neut_code, "neut_code/I");
	mytree->Branch("nuance_code", &genie_nuance_code, "nuance_code/I");
	mytree->Branch("wght", &genie_wght, "wght/D");
	mytree->Branch("xs", &genie_xs, "xs/D");
	mytree->Branch("ys", &genie_ys, "ys/D");
	mytree->Branch("ts", &genie_ts, "ts/D");
	mytree->Branch("Q2s", &genie_Q2s, "Q2s/D");
	mytree->Branch("Ws", &genie_Ws, "Ws/D");
	mytree->Branch("x", &genie_x, "x/D");
	mytree->Branch("y", &genie_y, "x/D");
	mytree->Branch("t", &genie_t, "t/D");
	mytree->Branch("Q2", &genie_Q2reco, "Q2/D");
	mytree->Branch("W", &genie_W, "W/D");
	mytree->Branch("EvRF", &genie_EvRF, "EvRF/D");
	mytree->Branch("Ev", &genie_Ev, "Ev/D");
	mytree->Branch("pxv", &genie_pxv, "pxv/D");
	mytree->Branch("pyv", &genie_pyv, "pyv/D");
	mytree->Branch("pzv", &genie_pzv, "pzv/D");
	mytree->Branch("En", &genie_En, "En/D");
	mytree->Branch("pxn", &genie_pxn, "pxn/D");
	mytree->Branch("pyn", &genie_pyn, "pyn/D");
	mytree->Branch("pzn", &genie_pzn, "pzn/D");
	mytree->Branch("El", &genie_El, "El/D");
	mytree->Branch("pxl", &genie_pxl, "pxl/D");
	mytree->Branch("pyl", &genie_pyl, "pyl/D");
	mytree->Branch("pzl", &genie_pzl, "pzl/D");
	mytree->Branch("pl", &genie_pl, "pl/D");
	mytree->Branch("cthl", &genie_cthl, "cthl/D");
	mytree->Branch("nfp", &genie_nfp, "nfp/I");
	mytree->Branch("nfn", &genie_nfn, "nfn/I");
	mytree->Branch("nfpip", &genie_nfpip, "nfpip/I");
	mytree->Branch("nfpim", &genie_nfpim, "nfpim/I");
	mytree->Branch("nfpi0", &genie_nfpi0, "nfpi0/I");
	mytree->Branch("nfkp", &genie_nfkp, "nfkp/I");
	mytree->Branch("nfkm", &genie_nfkm, "nfkm/I");
	mytree->Branch("nfk0", &genie_nfk0, "nfk0/I");
	mytree->Branch("nfem", &genie_nfem, "nfem/I");
	mytree->Branch("nfother", &genie_nfother, "nfother/I");
	mytree->Branch("nip", &genie_nip, "nip/I");
	mytree->Branch("nin", &genie_nin, "nin/I");
	mytree->Branch("nipip", &genie_nipip, "nipip/I");
	mytree->Branch("nipim", &genie_nipim, "nipim/I");
	mytree->Branch("nipi0", &genie_nipi0, "nipi0/I");
	mytree->Branch("nikp", &genie_nikp, "nikp/I");
	mytree->Branch("nikm", &genie_nikm, "nikm/I");
	mytree->Branch("nik0", &genie_nik0, "nik0/I");
	mytree->Branch("niem", &genie_niem, "niem/I");
	mytree->Branch("niother", &genie_niother, "niother/I");
	mytree->Branch("ni", &genie_ni, "ni/I");
	mytree->Branch("pdgi", &genie_pdgi,"pdgi[2]/I");
	mytree->Branch("resc", &genie_resc, "resc[2]/I");
	mytree->Branch("Ei", &genie_Ei, "Ei[2]/D");
	mytree->Branch("pxi", &genie_pxi, "pxi[2]/D");
	mytree->Branch("pyi", &genie_pyi, "pyi[2]/D");
	mytree->Branch("pzi", &genie_pzi, "pzi[2]/D");
	mytree->Branch("nf", &genie_nf, "nf/I");

	mytree->Branch("pdgf", &genie_pdgf, "pdgf[120]/I");
	mytree->Branch("Ef", &genie_Ef, "Ef[120]/D");
	mytree->Branch("pxf", &genie_pxf, "pxf[120]/D");
	mytree->Branch("pyf", &genie_pyf, "pyf[120]/D");
	mytree->Branch("pzf", &genie_pzf, "pzf[120]/D");
	mytree->Branch("pf", &genie_pf, "pf[120]/D");
	mytree->Branch("cthf", &genie_cthf, "cthf[120]/D");
	mytree->Branch("vtxx", &genie_vtxx, "vtxx/D");
	mytree->Branch("vtxy", &genie_vtxy, "vtxy/D");
	mytree->Branch("vtxz", &genie_vtxz, "vtxz/D");
	mytree->Branch("vtxt", &genie_vtxt, "vtxt/D");
	mytree->Branch("sumKEf", &genie_sumKEf, "sumKEf/D");
	mytree->Branch("calresp0", &genie_calresp0, "calresp0/D");

	int shift = 60;

// __________________________________________________________________________________________________________________________________________________________________________________________________






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
  fiducialcut-InitEClimits();

/** Beginning of Event Loop **/
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<200000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
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
          //SetMomCorrParameters(); Functions is missing F.H. 08/01/19
          fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam, en_beam_Ecal, en_beam_Eqe);
          SetFiducialCutParameters(fbeam_en);
    }

    int n_elec = 0;
    const int ind_em=0; //Index for electron
    if (ec[ind_em] <=0) {
      std::cout << "Possible problem with making electron ec vector. EC index below/equal Zero: ec[ind_em] =  " << ec[ind_em] << std::endl;
      continue;
    }
    if (sc[ind_em] <=0) {
      std::cout << "Possible problem with making electron ec vector. SC index below/equal zero: sc[ind_em] =  " << sc[ind_em] << std::endl;
      continue;
    }

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

genie_neu = 11;
genie_fspl = 11;
genie_tgt = TargetPdgCode;
genie_Z = TargetZ;
genie_A = TargetA;
genie_hitnuc = -99.;
genie_hitqrk = -99.;
genie_resid = -99.;
genie_sea = false;
genie_qel = false;
genie_mec = false;
genie_res = false;
genie_dis = false;
genie_coh = false;
genie_dfr = false;
genie_imd = false;
genie_imdanh = false;
genie_singlek = false;
genie_nuel = false;
genie_em = true;
genie_CC = false;
genie_nc = false;
genie_charm = false;
genie_neut_code = -99;
genie_nuance_code = -99;
genie_wght = 1.;
genie_xs = -99;
genie_ys = -99;
genie_ts = -99;
genie_Q2s = -99;
genie_Ws = -99;
genie_Ev = en_beam[fbeam_en];
genie_EvRF = -99.;
genie_pxv = 0.;
genie_pyv = 0.;
genie_pzv = sqrt(en_beam[fbeam_en] * en_beam[fbeam_en] - e_mass*e_mass); 
genie_En = -99.;
genie_pxn = -99.;
genie_pyn = -99.;
genie_pzn = -99.;
genie_nfkp = -99;
genie_nfkm = -99;
genie_nfk0 = -99;
genie_nfem = -99;
genie_nfother = -99;
genie_nip = -99;
genie_nin = -99;
genie_nipip = -99;
genie_nipim = -99;
genie_nipi0 = -99;
genie_nikp = -99;
genie_nikm = -99;
genie_nik0 = -99;
genie_niem = -99;
genie_niother = -99;
genie_ni = -99;
genie_calresp0 = -99;
genie_y = -99.;
genie_t = -99;

int NEventsTotal = 0;

for (int WhichInitialStateParticle = 0; WhichInitialStateParticle < InitialStateParticles; WhichInitialStateParticle ++) {

	genie_pdgi[WhichInitialStateParticle] = -99;
	genie_resc[WhichInitialStateParticle] = -99;
	genie_Ei[WhichInitialStateParticle] = -99.;
	genie_pxi[WhichInitialStateParticle] = -99.;
	genie_pyi[WhichInitialStateParticle] = -99.;
	genie_pzi[WhichInitialStateParticle] = -99.;

}

for (int WhichFinalStateParticle = 0; WhichFinalStateParticle < FinalStateParticles; WhichFinalStateParticle ++) {

	genie_pdgf[WhichFinalStateParticle] = -99;
	genie_Ef[WhichFinalStateParticle] = -99.;
	genie_pxf[WhichFinalStateParticle] = -99.;
	genie_pyf[WhichFinalStateParticle] = -99.;
	genie_pzf[WhichFinalStateParticle] = -99.;
	genie_pf[WhichFinalStateParticle] = -99.;
	genie_cthf[WhichFinalStateParticle] = -99.;

}

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

    //Main Fiducial Cut for Electron
    if( !EFiducialCut(fbeam_en, el_mom1) ) continue;//theta, phi cuts
    if( !CutUVW(e_ec_xyz1) )               continue; //u>60, v<360, w<400

    //General cut on EC, SC, CC hit and q (charge) for all events
    if( ec[ind_em] < 0.5 ||  sc[ind_em] < 0.5 ||  cc[ind_em] < 0.5 || q[ind_em] >=0)
    {
      continue;
    }

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

    // --------------------------------------------------------------------------------------------

    // apapadop

    //Electron vertex cut
    if( !(el_vert_corr < vert_max[ftarget] && el_vert_corr > vert_min[ftarget]) ) continue;

    // --------------------------------------------------------------------------------------------

    //Main Electron 4-Vectors with and without momentum correction. Units are GeV
    TLorentzVector V4_el_uncorr(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em],p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass));
    TLorentzVector V4_el(elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cx[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cy[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]*elmom_corr_fact[el_ec_sector-1]*elmom_corr_fact[el_ec_sector-1]+e_mass*e_mass));

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

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

genie_x = x_berk;
genie_Q2reco = Q2;
genie_W = W_var;
genie_El = V4_el.E();
genie_pxl = V4_el.Px();
genie_pyl = V4_el.Py();
genie_pzl = V4_el.Pz();
genie_pl = V4_el.Rho();
genie_cthl = cos(V4_el.Theta());
genie_vtxx = 0;
genie_vtxy = 0;
genie_vtxz = 0;

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    //Index and number variables for neutral particles
    int ec_index_n[20];
    int ec_num_n = 0;
    bool ec_radstat_n[20] = {false};

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

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

int IndexInFinalStateParticleArray = 0;

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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

            //proton pid cut and momentum > 0.3GeV cut to get rid of low momentum protons that have a high energy loss and we don't know the efficiency precisely for
	          if(delta < fsum_prot->Eval(p[ind_p]) && delta > fsub_prot->Eval(p[ind_p]) && p[ind_p] >= prot_accept_mom_lim){

              TLorentzVector V4_uncorrprot(p[ind_p]*cx[ind_p],p[ind_p]*cy[ind_p],p[ind_p]*cz[ind_p],TMath::Sqrt(p[ind_p]*p[ind_p]+ m_prot*m_prot ) );
              p_vert_corr = vz[ind_p]+vz_corr(vz_corr_func,prot_phi_mod,prot_theta);


	            if(PFiducialCut(fbeam_en, V4_uncorrprot.Vect())){ //proton fiducial cuts


                //main vertex cut for protons
	               if( (el_vert_corr-p_vert_corr) > vertdiff_min[ftarget] && (el_vert_corr-p_vert_corr) < vertdiff_max[ftarget] ){
	                   num_p = num_p + 1;
	                   index_p[num_p-1] = i;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

genie_pdgf[IndexInFinalStateParticleArray] = id[ind_p];
genie_Ef[IndexInFinalStateParticleArray] = V4_uncorrprot.E();
genie_pxf[IndexInFinalStateParticleArray] = V4_uncorrprot.Px();
genie_pyf[IndexInFinalStateParticleArray] = V4_uncorrprot.Py();
genie_pzf[IndexInFinalStateParticleArray] = V4_uncorrprot.Pz();
genie_pf[IndexInFinalStateParticleArray] = V4_uncorrprot.Rho();
genie_cthf[IndexInFinalStateParticleArray] = V4_uncorrprot.CosTheta();

double CorrectedProtonMomentum = ProtonMomCorrection_He3_4Cell(ftarget,V4_uncorrprot,p_vert_corr);
V4_uncorrprot.SetRho(CorrectedProtonMomentum);
double CorrectedProtonE = sqrt(CorrectedProtonMomentum*CorrectedProtonMomentum + m_prot*m_prot);
V4_uncorrprot.SetE(CorrectedProtonE);


genie_pdgf[IndexInFinalStateParticleArray+shift] = id[ind_p];
genie_Ef[IndexInFinalStateParticleArray+shift] = V4_uncorrprot.E();
genie_pxf[IndexInFinalStateParticleArray+shift] = V4_uncorrprot.Px();
genie_pyf[IndexInFinalStateParticleArray+shift] = V4_uncorrprot.Py();
genie_pzf[IndexInFinalStateParticleArray+shift] = V4_uncorrprot.Pz();
genie_pf[IndexInFinalStateParticleArray+shift] = V4_uncorrprot.Rho();
genie_cthf[IndexInFinalStateParticleArray+shift] = V4_uncorrprot.CosTheta();
IndexInFinalStateParticleArray++;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------


	               } // end of vertex cuts for protons

	             } //end of fiducial cuts
             } //end of if delta condition
	        } //end of if loop to check for proton id

		// -----------------------------------------------------------------------------------------------------------------------------------------------------

		// PiMinus


	        if(q[i] < 0 && sc[i] > 0 && dc[i] > 0 && stat[i] > 0 ) //negative particle, possibly pi minus
	        {
	             V3_pimi.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	             beta = p[i]/TMath::Sqrt(p[i]*p[i]+m_pimi*m_pimi);
	             delta = sc_t[sc[i]-1]-sc_r[sc[i]-1]/(beta*c*ns_to_s) - tr_time;

               //fsum_pimi and fsub_pimi are TF1 Functions for piminus delta cuts
      	       fsub_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);
	             fsum_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);

               //Pion pid delta cut
          	   if(delta < fsum_pimi->Eval(p[i]) && delta > fsub_pimi->Eval(p[i])){

	                pimi_vert=vz[i];
	                pimi_phi = TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
	                pimi_phi_mod = pimi_phi + 30;  //Add extra 30 degree rotation in phi
	                if (pimi_phi_mod<0)  pimi_phi_mod = pimi_phi_mod + 360;  //Pi minus is between 0 and 360 degree
	                pimi_theta = TMath::ACos(cz[i])*TMath::RadToDeg();
	                pimi_vert_corr = pimi_vert+vz_corr(vz_corr_func, pimi_phi_mod,pimi_theta);


	                if(PimiFiducialCut(fbeam_en, V3_pimi, &pimi_phimin, &pimi_phimax)){  //Pi minus fiducial cuts

                    for(int k=1;k<=6;k++){ //k is sector number

                    //main vertex cut for pi minus
		                if(abs(el_vert_corr-pimi_vert_corr) < pimi_vertcut){

		                    num_pimi = num_pimi + 1;
		                    num_pi = num_pi + 1;
		                    num_pi_phot = num_pi_phot + 1;
                        num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
		                    index_pimi[num_pimi - 1] = i;
		                    index_pi[num_pi - 1] = i;
		                    ind_pi_phot[num_pi_phot - 1] = i;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

double PiMinusP = V3_pimi.Mag();
double PiMinusE = TMath::Sqrt(PiMinusP*PiMinusP + m_pimi * m_pimi);

genie_pdgf[IndexInFinalStateParticleArray] = id[i];
genie_Ef[IndexInFinalStateParticleArray] = PiMinusE;
genie_pxf[IndexInFinalStateParticleArray] = V3_pimi.X();
genie_pyf[IndexInFinalStateParticleArray] = V3_pimi.Y();
genie_pzf[IndexInFinalStateParticleArray] = V3_pimi.Z();
genie_pf[IndexInFinalStateParticleArray] = PiMinusP;
genie_cthf[IndexInFinalStateParticleArray] = V3_pimi.CosTheta();
IndexInFinalStateParticleArray++;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------


		                } //if piminus vertex cut
	               } //if Piminus fiducials
	            } //if piminus delta cut
	        } //if negative particle


// --------------------------------------------------------------------------------------------------------------------------------------------------------------------


		// PiPlus

			    if(q[i] > 0 &&  sc[i] > 0 && dc[i] > 0 && stat[i] > 0) //positive particle. possible pi plus
	        {
	             V3_pipl.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	             beta = p[i]/TMath::Sqrt(p[i]*p[i]+m_pipl*m_pipl);
	             delta= sc_t[sc[i]-1]-sc_r[sc[i]-1]/(beta*c*ns_to_s) - tr_time;

               //fsum_piplus and fsub_piplus are TF1 Functions for piplus delta cuts
	             fsub_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);
	             fsum_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);

               //Pion pid delta cut
               if(delta < fsum_pipl->Eval(p[i]) && delta > fsub_pipl->Eval(p[i])){

	                pipl_vert=vz[i];
	                pipl_phi = TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
            	    pipl_phi_mod = pipl_phi + 30; //Add 30 degrees
	                if (pipl_phi_mod < 0)  pipl_phi_mod = pipl_phi_mod + 360;  //Pi plus is between 0 and 360 degree
	                pipl_theta = TMath::ACos(cz[i])*TMath::RadToDeg();
	                pipl_vert_corr = pipl_vert + vz_corr(vz_corr_func,pipl_phi_mod,pipl_theta);

	                if (PiplFiducialCut(fbeam_en, V3_pipl, &pipl_phimin, &pipl_phimax)){ //Pi Plus fiducial cut


		                if (abs(el_vert_corr-pipl_vert_corr) < pipl_vertcut){ //pi plus vertex cut
		                    num_pipl = num_pipl + 1;
		                    num_pi  = num_pi + 1;
		                    num_pi_phot = num_pi_phot + 1;
                        num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
		                    index_pipl[num_pipl - 1] = i;
		                    index_pi[num_pi - 1] = i;
		                    ind_pi_phot[num_pi_phot - 1] = i;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

double PiPlusP = V3_pipl.Mag();
double PiPlusE = TMath::Sqrt(PiPlusP*PiPlusP + m_pipl * m_pipl);

genie_pdgf[IndexInFinalStateParticleArray] = id[i];
genie_Ef[IndexInFinalStateParticleArray] = PiPlusE;
genie_pxf[IndexInFinalStateParticleArray] = V3_pipl.X();
genie_pyf[IndexInFinalStateParticleArray] = V3_pipl.Y();
genie_pzf[IndexInFinalStateParticleArray] = V3_pipl.Z();
genie_pf[IndexInFinalStateParticleArray] = PiPlusP;
genie_cthf[IndexInFinalStateParticleArray] = V3_pipl.CosTheta();
IndexInFinalStateParticleArray++;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------


		                 }	 //vert cut ends
	                 } //fidcut ends
	              }//delta cut ends
	         }//pipl ends


		// -----------------------------------------------------------------------------------------------------------------------------------------------------


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

	              if(neut_beta_corr > EC_photon_beta[ftarget]){   //photon identification


	                 if(Phot_fid(V3_phot_angles)){ //photon fiducial function

                     ec_num_n = ec_num_n + 1;
	                   ec_index_n[ec_num_n - 1] = i;
	                   num_pi_phot = num_pi_phot + 1;
	                   ind_pi_phot[num_pi_phot - 1] = i;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

genie_pdgf[IndexInFinalStateParticleArray] = id[i];
genie_Ef[IndexInFinalStateParticleArray] = V3_phot_angles.Mag();
genie_pxf[IndexInFinalStateParticleArray] = V3_phot_angles.X();
genie_pyf[IndexInFinalStateParticleArray] = V3_phot_angles.Y();
genie_pzf[IndexInFinalStateParticleArray] = V3_phot_angles.Z();
genie_pf[IndexInFinalStateParticleArray] = V3_phot_angles.Mag();
genie_cthf[IndexInFinalStateParticleArray] = V3_phot_angles.CosTheta();
IndexInFinalStateParticleArray++;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------



                     //Photon EC energy deposit
 	                   photon_ece = TMath::Max( ec_ei[ec[i] - 1] + ec_eo[ec[i] - 1],etot[ec[i] - 1]);

                     //Cut on Radiation photon via angle with respect to the electron
                     //within 40 degrees in theta and 30 degrees in phi
	                   if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg() < phot_rad_cut && abs(neut_phi_mod-el_phi_mod) < phot_e_phidiffcut ) {
		                     ec_radstat_n[num_pi_phot - 1] = true; //select radiation photons

	                   }
	                   if(!ec_radstat_n[num_pi_phot - 1]) num_pi_phot_nonrad = num_pi_phot_nonrad + 1;


	                }//end if photon fiducial
	             }//n beta
	         } //if neutral particles

    } //end of hadron loop

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

// apapadop

genie_nfp = num_p;
genie_nfn = 0;
genie_nfpip = num_pipl;
genie_nfpim = num_pimi;
genie_nfpi0 = ec_num_n;
genie_nf = genie_nfp + genie_nfn + genie_nfpip + genie_nfpim + genie_nfpi0;

genie_iev = NEventsTotal;
NEventsTotal++;

mytree->Fill();

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


	// This is the end of the filter

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  } //end of event loop (jentry)

  gStyle->SetOptFit(1);


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

void  FilterData::prot3_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot[3],TVector3  V3prot_uncorr[3],TLorentzVector V4el,double Ecal_3pto2p[][2],double  pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3],double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det){

  std::string fbeam_en = beam_en;
  double m_prot=0.9382720813;
  const int N_3p=3, N_2p=2;
  double N_p2[N_3p]={0},N_p1det[3][N_3p]={0};
  TVector3 V3_3p_rot[N_3p],V3_2p_rot[N_3p],V3_prot_el_3pto2p[N_2p],V3_prot_el_3pto1p[N_3p];
  double rot_angle;
  double N_pthree=0;
  bool prot_stat[N_3p];
  int count =0;



     for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());

       for(int i=0;i<N_3p;i++) {
	 V3_3p_rot[i]= V3prot_uncorr[i];
	 V3_3p_rot[i].Rotate(rot_angle,V3q);
       }


       for(int ind_p=0;ind_p<N_3p;ind_p++) prot_stat[ind_p]=PFiducialCut(fbeam_en, V3_3p_rot[ind_p]);

       if( prot_stat[0]  && !prot_stat[1]  && !prot_stat[2])  N_p1[0]=N_p1[0]+1;
       if(!prot_stat[0] &&   prot_stat[1]  && !prot_stat[2])  N_p1[1]=N_p1[1]+1;
       if(!prot_stat[0] &&  !prot_stat[1]  &&  prot_stat[2])  N_p1[2]=N_p1[2]+1;
       if( prot_stat[0]  &&  prot_stat[1]  &&  prot_stat[2])  N_pthree=N_pthree+1;

       if(prot_stat[0] && prot_stat[1] && !prot_stat[2])   N_p2[0]=N_p2[0]+1;
       if(prot_stat[0] && !prot_stat[1] && prot_stat[2])   N_p2[1]=N_p2[1]+1;
       if(!prot_stat[0] && prot_stat[1] && prot_stat[2])   N_p2[2]=N_p2[2]+1;


     }//for loop of 3p rotations ends

 //-----------------------------------------  3p to 1p  -----------------------------------------------------------------------
     for(int j=0;j<N_3p;j++)    { //looping through 1p combinations out of 3protons

       V3_prot_el_3pto1p[j]=V4el.Vect()+ V3prot[j];
       Ecal_3pto1p[j]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[j].Mag()*V3prot[j].Mag())-m_prot+bind_en[target_name];
       pmiss_perp_3pto1p[j]=TMath::Sqrt(V3_prot_el_3pto1p[j].Px()*V3_prot_el_3pto1p[j].Px()+V3_prot_el_3pto1p[j].Py()*V3_prot_el_3pto1p[j].Py());
     }


 //-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------
     for(int ind1=0;ind1<N_3p;ind1++){   //looping through 2p combinations  out of 3p
     for(int ind2=0;ind2<N_3p;ind2++){
	 if(ind1!=ind2 && ind1<ind2){

	     for(int g1=0; g1<N_tot; g1++){

	       rot_angle=gRandom->Uniform(0,2*TMath::Pi());

	       V3_2p_rot[ind1]=V3prot_uncorr[ind1];
	       V3_2p_rot[ind2]=V3prot_uncorr[ind2];
	       V3_2p_rot[ind1].Rotate(rot_angle,V3q);
	       V3_2p_rot[ind2].Rotate(rot_angle,V3q);

	       if(PFiducialCut(fbeam_en, V3_2p_rot[ind1])  && !PFiducialCut(fbeam_en, V3_2p_rot[ind2])) N_p1det[count][0]=N_p1det[count][0]+1;
	       if(!PFiducialCut(fbeam_en, V3_2p_rot[ind1]) && PFiducialCut(fbeam_en, V3_2p_rot[ind2]))  N_p1det[count][1]=N_p1det[count][1]+1;
	       if(PFiducialCut(fbeam_en, V3_2p_rot[ind1])  && PFiducialCut(fbeam_en, V3_2p_rot[ind2]))  N_p1det[count][2]=N_p1det[count][2]+1;
	     }


	     if( N_p1det[count][2]!=0   && N_pthree!=0){


	       V3_prot_el_3pto2p[0]=V4el.Vect()+ V3prot[ind1];
	       Ecal_3pto2p[count][0]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[ind1].Mag()*V3prot[ind1].Mag())-m_prot+bind_en[target_name];
	       pmiss_perp_3pto2p[count][0]=TMath::Sqrt(V3_prot_el_3pto2p[0].Px()*V3_prot_el_3pto2p[0].Px()+V3_prot_el_3pto2p[0].Py()*V3_prot_el_3pto2p[0].Py());
	       P3pto2p[count][0]=N_p1det[count][0]/N_p1det[count][2]*(N_p2[count]/N_pthree);

	       V3_prot_el_3pto2p[1]=V4el.Vect()+ V3prot[ind2];
	       Ecal_3pto2p[count][1]=V4el.E()+TMath::Sqrt(m_prot*m_prot+V3prot[ind2].Mag()*V3prot[ind2].Mag())-m_prot+bind_en[target_name];
	       pmiss_perp_3pto2p[count][1]=TMath::Sqrt(V3_prot_el_3pto2p[1].Px()*V3_prot_el_3pto2p[1].Px()+V3_prot_el_3pto2p[1].Py()*V3_prot_el_3pto2p[1].Py());
	       P3pto2p[count][1]=N_p1det[count][1]/N_p1det[count][2]*(N_p2[count]/N_pthree);
	     }


    count=count +1;

	 }
     }
   }
  *N_p3det=N_pthree;
}






void  FilterData::prot2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot[2],TVector3  V3prot_uncorr[2],TLorentzVector V4el,double Ecal_2pto1p[2],double  pmiss_perp_2pto1p[2],double  P2pto1p[2], double *Nboth){

  std::string fbeam_en = beam_en;
  const int N2=2;
  double rot_angle, N_p2to1[N2]={0},m_prot=0.9382720813;
  TVector3 V3_2prot[N2], V3_prot_el_2pto1p[N2];
  double N_2=0;



  for(int g1=0; g1<N_tot; g1++){


    rot_angle=gRandom->Uniform(0,2*TMath::Pi());

    V3_2prot[0]=V3prot_uncorr[0];
    V3_2prot[1]=V3prot_uncorr[1];
    V3_2prot[0].Rotate(rot_angle,V3q);
    V3_2prot[1].Rotate(rot_angle,V3q);


    if(PFiducialCut(fbeam_en, V3_2prot[0])  && !PFiducialCut(fbeam_en, V3_2prot[1])) N_p2to1[0]=N_p2to1[0]+1;
    if(!PFiducialCut(fbeam_en, V3_2prot[0]) && PFiducialCut(fbeam_en, V3_2prot[1]))  N_p2to1[1]=N_p2to1[1]+1;
    if(PFiducialCut(fbeam_en, V3_2prot[0])  && PFiducialCut(fbeam_en, V3_2prot[1]))  N_2=N_2+1;
  }

   //-----------------------------------------  2p to 1p  -----------------------------------------------------------------------
  V3_prot_el_2pto1p[0]=V4el.Vect()+ V3prot[0];
    Ecal_2pto1p[0]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[0].Mag()*V3prot[0].Mag())-m_prot+bind_en[target_name];
    pmiss_perp_2pto1p[0]=TMath::Sqrt(V3_prot_el_2pto1p[0].Px()*V3_prot_el_2pto1p[0].Px()+V3_prot_el_2pto1p[0].Py()*V3_prot_el_2pto1p[0].Py());

    V3_prot_el_2pto1p[1]=V4el.Vect()+ V3prot[1];
    Ecal_2pto1p[1]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[1].Mag()*V3prot[1].Mag())-m_prot+bind_en[target_name];
    pmiss_perp_2pto1p[1]=TMath::Sqrt(V3_prot_el_2pto1p[1].Px()*V3_prot_el_2pto1p[1].Px()+V3_prot_el_2pto1p[1].Py()*V3_prot_el_2pto1p[1].Py());


  if( N_2!=0){
  P2pto1p[0]=N_p2to1[0]/N_2;
    P2pto1p[1]=N_p2to1[1]/N_2;
  }
  else{
  P2pto1p[0]=0;
    P2pto1p[1]=0;
  }

  *Nboth=N_2;
}






void FilterData::prot1_pi1_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi, int q_pi,bool radstat, double *N_pi_p,double *N_nopi_p){

  std::string fbeam_en = beam_en;
  double rotation_ang;
  TVector3 V3_pi_rot, V3_p_rot;
  Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;
  bool pi_stat=true;

     double Npi_p = 0;
     double Nnopi_p = 0;

     for(int g=0; g<N_tot; g++){

       rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
       V3_p_rot= V3prot;
       V3_p_rot.Rotate(rotation_ang,V3q);

       if(!radstat){
       V3_pi_rot=V3pi;
       V3_pi_rot.Rotate(rotation_ang,V3q);
       pi_stat=Pi_phot_fid_united(fbeam_en, V3_pi_rot,q_pi);
}

       if(PFiducialCut(fbeam_en, V3_p_rot)  && pi_stat) Npi_p=Npi_p+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && !pi_stat) Nnopi_p=Nnopi_p+1;

     }
     *N_pi_p=Npi_p;
     *N_nopi_p=Nnopi_p;
}






void FilterData::prot1_pi2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[2], int q_pi[2],bool radstat[2], double *P_1p0pi,double P_1p1pi[2]){

  std::string fbeam_en = beam_en;
  const int N_pi=2;
  double rotation_ang;
  TVector3 V3_rot_pi[2], V3_p_rot;
  bool status_pi[2]={true};

     double N_all = 0;
     double Nnopi = 0,N_1p1pi[2]={0};

     for(int g=0; g<N_tot; g++){

       rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
       V3_p_rot= V3prot;

       V3_p_rot.Rotate(rotation_ang,V3q);

       for(int i=0;i<N_pi;i++){
	 if(!radstat[i]){
	 V3_rot_pi[i]=V3pi[i];
	 V3_rot_pi[i].Rotate(rotation_ang,V3q);
	 status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
	 }
       }

       if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] ) N_1p1pi[0]=N_1p1pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] ) N_1p1pi[1]=N_1p1pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] ) N_all=N_all+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && !status_pi[1]) Nnopi=Nnopi+1;

     }


     if(N_all!=0){
       //----------------------1p2pi->1p0pi
       *P_1p0pi=Nnopi/N_all;

       //----------------------1p2pi->1p1pi->1p0pi
       double N_nopi_p = 0,N_pi_p=0;

       for(int h=0;h<N_pi;h++){

	 prot1_pi1_rot_func(fbeam_en, V3q,V3prot,V3pi[h],q_pi[h],radstat[h],&N_pi_p,&N_nopi_p);
	 if(N_pi_p!=0) P_1p1pi[h]=(N_1p1pi[h]/N_all)*(N_nopi_p/N_pi_p);
	 else  P_1p1pi[h]=0;
       }
     }   //N_all!=0 statement

     else{
       *P_1p0pi=0;
       P_1p1pi[0]=0;
       P_1p1pi[1]=0;
     }
}




void FilterData::prot1_pi3_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[3], int q_pi[3],bool radstat[3],double *P_tot){

  std::string fbeam_en = beam_en;
  *P_tot=0;
  const int N_pi=3;
  double rotation_ang;
  TVector3 V3_rot_pi[N_pi], V3_p_rot;
  bool status_pi[N_pi]={true};
     double N_all = 0;
     double Nnopi = 0,N_1p1pi[3]={0},N_1p2pi[3]={0};

     for(int g=0; g<N_tot; g++){

       rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
       V3_p_rot= V3prot;

       V3_p_rot.Rotate(rotation_ang,V3q);

       for(int i=0;i<N_pi;i++){
	 if(!radstat[i]){
	 V3_rot_pi[i]=V3pi[i];
	 V3_rot_pi[i].Rotate(rotation_ang,V3q);
	 status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
	 }
       }

       if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] && !status_pi[2]) N_1p1pi[0]=N_1p1pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] && !status_pi[2]) N_1p1pi[1]=N_1p1pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && !status_pi[1] && status_pi[2]) N_1p1pi[2]=N_1p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] && !status_pi[2]) N_1p2pi[0]=N_1p2pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && !status_pi[1] && status_pi[2]) N_1p2pi[1]=N_1p2pi[1]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && status_pi[1] && status_pi[2]) N_1p2pi[2]=N_1p2pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && status_pi[0]  && status_pi[1] && status_pi[2])  N_all=N_all+1;
       if(PFiducialCut(fbeam_en, V3_p_rot)  && !status_pi[0]  && !status_pi[1] && !status_pi[2])Nnopi=Nnopi+1;
     }

     if(N_all!=0){
       //----------------------1p3pi->1p0pi
       double P_1p0pi=0;

       P_1p0pi=Nnopi/N_all;

       //----------------------1p3pi->1p1pi->1p0pi
       double N_nopi_p = 0,N_pi_p=0;
       const int N2pi=2;
       TVector3 V3_pion[N2pi];
       int q_pion[N2pi];
       bool rad_stat[N2pi];
       double  P_1p1pi=0,P_1p2pi=0;

     for(int h=0;h<N_pi;h++){

       prot1_pi1_rot_func(fbeam_en, V3q,V3prot,V3pi[h],q_pi[h],radstat[h],&N_pi_p,&N_nopi_p);
	 if(N_pi_p!=0) P_1p1pi=P_1p1pi+(N_1p1pi[h]/N_all)*(N_nopi_p/N_pi_p);

  //----------------------1p3pi->1p2pi->1p0pi
	 double P_1p1pion[N2pi]={0},P_1p0pion=0;

	 if(h==0)   {
	   V3_pion[0]=   V3pi[0];   V3_pion[1]=V3pi[1];
	   q_pion[0]=    q_pi[0];    q_pion[1]=q_pi[1];
	   rad_stat[0]=  radstat[0]; rad_stat[1]=radstat[1];
	 }
	 if(h==1)   {
	   V3_pion[0]=   V3pi[0];    V3_pion[1]=V3pi[2];
	   q_pion[0]=    q_pi[0];     q_pion[1]=q_pi[2];
	   rad_stat[0]=  radstat[0]; rad_stat[1]=radstat[2];
	 }
	 if(h==2)   {
	   V3_pion[0]=   V3pi[1];   V3_pion[1]=V3pi[2];
	   q_pion[0]=    q_pi[1];    q_pion[1]=q_pi[2];
	   rad_stat[0]=  radstat[1]; rad_stat[1]=radstat[2];
	 }
 //----------------------1p3pi->1p2pi->1p1pi->1p0pi
	 prot1_pi2_rot_func(fbeam_en, V3q,V3prot,V3_pion, q_pion,rad_stat,&P_1p0pion,P_1p1pion);
	 P_1p2pi=P_1p2pi+(N_1p2pi[h]/N_all)*(P_1p0pion-P_1p1pion[0]-P_1p1pion[1]);

     }//for loop ends

     *P_tot=P_1p2pi+P_1p1pi-P_1p0pi;
     }   //N_all!=0 statement

     else{
       *P_tot=0;
     }

}




void FilterData::prot2_pi1_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_2p1pi_to2p0pi[2],double p_miss_perp_2p1pi_to2p0pi[2],double P_2p1pito2p0pi[2],double P_2p1pito1p1pi[2],double P_2p1pito1p0pi[2],double *P_tot){

  std::string fbeam_en = beam_en;
  const int N_2prot=2;
  TVector3 V3_2p_rotated[N_2prot],V3_1pirot;
  bool pi1_stat=true;
  double N_2p_0pi=0,N_all=0,N_1p_1pi[N_2prot]={0},N_1p_0pi[N_2prot]={0};
  double P_2pto1p[N_2prot]={0},N_2p_det=0;
  double   N_pidet=0,N_piundet=0,rot_angle;
  *P_tot=0;


     for(int g=0; g<N_tot; g++){

     rot_angle=gRandom->Uniform(0,2*TMath::Pi());


     V3_2p_rotated[0]=V3_2prot_uncorr[0];
     V3_2p_rotated[1]=V3_2prot_uncorr[1];
     V3_2p_rotated[0].Rotate(rot_angle,V3_q);
     V3_2p_rotated[1].Rotate(rot_angle,V3_q);

     if(!radstat){
       V3_1pirot=V3_1pi;
       V3_1pirot.Rotate(rot_angle,V3_q);
       pi1_stat=Pi_phot_fid_united(fbeam_en, V3_1pirot, q_pi);
     }

     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi1_stat) N_2p_0pi=N_2p_0pi+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[0]=N_1p_1pi[0]+1;
     if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_1p_1pi[1]=N_1p_1pi[1]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi1_stat) N_1p_0pi[0]=N_1p_0pi[0]+1;
     if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi1_stat) N_1p_0pi[1]=N_1p_0pi[1]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi1_stat) N_all=N_all+1;
   }

 //---------------------------------- 2p 1pi ->2p 0pi   ----------------------------------------------
  if(N_all!=0){

  prot2_rot_func(fbeam_en, V3_q, V3_2prot_corr,V3_2prot_uncorr, V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2pto1p ,&N_2p_det);

    for(int z=0;z<N_2prot;z++){

	P_2p1pito2p0pi[z]=(N_2p_0pi/N_all)*P_2pto1p[z];

      //---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

	prot1_pi1_rot_func(fbeam_en, V3_q,V3_2prot_uncorr[z],V3_1pi, q_pi,radstat, &N_pidet,&N_piundet);
      if(N_pidet!=0) P_2p1pito1p1pi[z]=(N_1p_1pi[z]/N_all)*(N_piundet/N_pidet);
      else P_2p1pito1p1pi[z]=0;

      //---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------
      P_2p1pito1p0pi[z]=(N_1p_0pi[z]/N_all);

      *P_tot=*P_tot+P_2p1pito2p0pi[z]+P_2p1pito1p1pi[z]-P_2p1pito1p0pi[z];

    }//looping through 2p
    }

if(N_all==0){
    P_2p1pito2p0pi[0]=P_2p1pito2p0pi[1]=0;
    P_2p1pito1p1pi[0]= P_2p1pito1p1pi[1]=0;
    P_2p1pito1p0pi[0]=P_2p1pito1p0pi[1]=0;
    *P_tot=0;
  }

}


void FilterData::prot2_pi2_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], int q_pi[2],bool radstat[2],TLorentzVector V4_el, double Ecal_2p2pi[2],double p_miss_perp_2p2pi[2],double P_tot_2p[2]){

  std::string fbeam_en = beam_en;
  const int N_2prot=2,N_2pi=2;
  TVector3 V3_2p_rotated[N_2prot],V3_2pirot[N_2pi];
  bool pi2_stat[N_2pi]={true};
  double   rot_angle;
  double N_2p_0pi=0,N_2p_1pi[N_2pi]={0},N_1p_2pi[N_2prot]={0},N_all=0,N_1p_1pi[N_2prot][N_2pi]={0},N_1p_0pi[N_2prot]={0};
  double   N_pidet=0,N_piundet=0;
  double P_2pto1p[N_2prot]={0},N_2p_det=0;
  double P_1p0pi=0,P_1p1pi[N_2pi]={0};
  double P_2p1pito2p0pi[2]={0}, P_2p1pito1p1pi[2]={0},P_2p1pito1p0pi[2]={0},Ptot=0;
  double P_2p2pito1p0pi[N_2prot]={0},P_2p2pito1p1pi[N_2prot]={0},P_2p2pito1p2pi[N_2prot]={0},P_2p2pito2p1pi[N_2prot]={0};
  P_tot_2p[0]=P_tot_2p[1]=0;

     for(int g=0; g<N_tot; g++){

     rot_angle=gRandom->Uniform(0,2*TMath::Pi());

     for(int k=0; k<N_2pi; k++){

       V3_2p_rotated[k]=V3_2prot_uncorr[k];
       V3_2p_rotated[k].Rotate(rot_angle,V3_q);

       if(!radstat[k]){
	 V3_2pirot[k]=V3_2pi[k];
	 V3_2pirot[k].Rotate(rot_angle,V3_q);
	 pi2_stat[k]=Pi_phot_fid_united(fbeam_en, V3_2pirot[k], q_pi[k]);
       }
     }

     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_2p_1pi[0]=N_2p_1pi[0]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_2p_1pi[1]=N_2p_1pi[1]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && !pi2_stat[1]) N_2p_0pi=N_2p_0pi+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[0]=N_1p_2pi[0]+1;
     if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[1]=N_1p_2pi[1]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[0][0]=N_1p_1pi[0][0]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[0][1]=N_1p_1pi[0][1]+1;
     if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[1][0]=N_1p_1pi[1][0]+1;
     if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[1][1]=N_1p_1pi[1][1]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && !pi2_stat[1])  N_1p_0pi[0]=N_1p_0pi[0]+1;
     if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && !pi2_stat[1])  N_1p_0pi[1]=N_1p_0pi[1]+1;
     if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_all=N_all+1;
   }


  if(N_all!=0){

 prot2_rot_func(fbeam_en, V3_q, V3_2prot_corr,V3_2prot_uncorr, V4_el,Ecal_2p2pi,p_miss_perp_2p2pi,P_2pto1p ,&N_2p_det);

 for(int z=0;z<N_2prot;z++){
 //---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

   P_2p2pito1p0pi[z]= N_1p_0pi[z]/N_all;

 //---------------------------------- 2p 2pi ->1p 1pi   ----------------------------------------------

for(int k=0;k<N_2pi;k++){
  N_pidet=N_piundet=0;
  prot1_pi1_rot_func(fbeam_en, V3_q,V3_2prot_uncorr[z],V3_2pi[k], q_pi[k],radstat[k],&N_pidet,&N_piundet);
   if(N_pidet!=0) P_2p2pito1p1pi[z]=P_2p2pito1p1pi[z]+(N_1p_1pi[z][k]/N_all)*(N_piundet/N_pidet);
 }

 //---------------------------------- 2p 2pi ->2p 0pi   ----------------------------------------------

 if(N_2p_det!=0)P_2p1pito2p0pi[z]=(N_2p_0pi/N_all)*P_2pto1p[z];

      //---------------------------------- 2p 2pi ->1p 2pi   ----------------------------------------------

 P_1p0pi=P_1p1pi[0]=P_1p1pi[1]=0;
 prot1_pi2_rot_func(fbeam_en, V3_q,V3_2prot_uncorr[z],V3_2pi,q_pi,radstat,&P_1p0pi,P_1p1pi);
 P_2p2pito1p2pi[z]=(N_1p_2pi[z]/N_all)*(P_1p0pi-P_1p1pi[0]-P_1p1pi[1]);

 //---------------------------------- 2p 2pi ->2p 1pi   ----------------------------------------------

 P_2p1pito2p0pi[0]=P_2p1pito2p0pi[1]=0; P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0; P_2p1pito1p0pi[0]=P_2p1pito1p0pi[1]=0;Ptot=0;
 prot2_pi1_rot_func(fbeam_en, V3_q,V3_2prot_corr,V3_2prot_uncorr,V3_2pi[z], q_pi[z],radstat[z],V4_el,Ecal_2p2pi,p_miss_perp_2p2pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

  // P_2p2pito2p1pi[z]=(N_2p_1pi[0]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z])+(N_2p_1pi[1]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z]);
 //P_tot_2p[z]=-P_2p2pito1p0pi[z]+P_2p2pito1p1pi[z]+P_2p1pito2p0pi[z]+P_2p2pito1p2pi[z]+P_2p2pito2p1pi[z];

 //P_2p2pito2p1pi[z]=(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[0]- P_2p1pito1p1pi[0]+P_2p1pito1p0pi[0])+(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[1]- P_2p1pito1p1pi[1]+P_2p1pito1p0pi[1]);

 P_2p2pito2p1pi[0]= P_2p2pito2p1pi[0]+(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[0]- P_2p1pito1p1pi[0]+P_2p1pito1p0pi[0]);
 P_2p2pito2p1pi[1]= P_2p2pito2p1pi[1]+(N_2p_1pi[z]/N_all)*(-P_2p1pito2p0pi[1]- P_2p1pito1p1pi[1]+P_2p1pito1p0pi[1]);

    }//looping through 2p

P_tot_2p[0]=-P_2p2pito1p0pi[0]+P_2p2pito1p1pi[0]+P_2p1pito2p0pi[0]+P_2p2pito1p2pi[0]+P_2p2pito2p1pi[0];
P_tot_2p[1]=-P_2p2pito1p0pi[1]+P_2p2pito1p1pi[1]+P_2p1pito2p0pi[1]+P_2p2pito1p2pi[1]+P_2p2pito2p1pi[1];

  }//N_all!=0

if(N_all==0){
  P_tot_2p[0]= P_tot_2p[1]=0;
  }

}




void FilterData::prot3_pi1_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_3p1pi[3],double p_miss_perp_3p1pi[3],double P_tot_3p[3]){

  std::string fbeam_en = beam_en;
  const int N_3prot=3;
  TVector3 V3_3p_rotated[N_3prot],V3_pirot;
  bool pi_stat=true;
  double   rot_angle;

  double N_1p0pi[N_3prot]={0},N_all=0,N_1p1pi[N_3prot]={0},N_2p0pi[N_3prot]={0},N_2p1pi[N_3prot]={0},N_3p0pi=0;
  double  P_3p1pito1p0pi[N_3prot]={0},P_3p1pito1p1pi[N_3prot]={0};
  double   N_pidet=0,N_piundet=0;
  TVector3 V3_2p_corr[N_3prot],V3_2p_uncorr[N_3prot];
  double Ecal_2p[2],p_miss_perp_2p[2],P_2pto1p[2]={0},N_2p_det=0,P_3p1pito2p0pi[N_3prot]={0};
  int count=0;
  double N_p1[N_3prot]={0},N_p_three=0;
  double E_cal_3pto2p[N_3prot][2]={0},p_miss_perp_3pto2p[N_3prot][2]={0},P_3pto2p[N_3prot][2]={0},P_3p1pito3p0pi[N_3prot]={0};
  double Ecal_2p1pi[2],p_miss_perp_2p1pi[2],P_3p1pito2p1pi[N_3prot]={0};
  double P_2p1pito2p0pi[2]={0},P_2p1pito1p1pi[2]={0},P_2p1pito1p0pi[2]={0},Ptot=0;

     for(int g=0; g<N_tot; g++){

     rot_angle=gRandom->Uniform(0,2*TMath::Pi());
     for(int k=0; k<N_3prot; k++){

     V3_3p_rotated[k]=V3_3prot_uncorr[k];
     V3_3p_rotated[k].Rotate(rot_angle,V3_q);
     }

     if(!radstat){ V3_pirot=V3_pi;
     V3_pirot.Rotate(rot_angle,V3_q);
     pi_stat=Pi_phot_fid_united(fbeam_en, V3_pirot, q_pi);
     }

     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_1p0pi[0]=N_1p0pi[0]+1;
     if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_1p0pi[1]=N_1p0pi[1]+1;
     if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_1p0pi[2]=N_1p0pi[2]+1;
     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[0]=N_1p1pi[0]+1;
     if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[1]=N_1p1pi[1]+1;
     if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[2]=N_1p1pi[2]+1;
     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_2p0pi[0]=N_2p0pi[0]+1;
     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_2p0pi[1]=N_2p0pi[1]+1;
     if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_2p0pi[2]=N_2p0pi[2]+1;
     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[0]=N_2p1pi[0]+1;
     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[1]=N_2p1pi[1]+1;
     if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[2]=N_2p1pi[2]+1;
     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && !pi_stat)  N_3p0pi=N_3p0pi+1;
     if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_all=N_all+1;
   }


  if(N_all!=0){

//----------------------------------3p 1pi ->3p 0pi->1p0pi   ----------------------------------------------
          prot3_rot_func(fbeam_en, V3_q, V3_3prot_uncorr,V3_3prot_corr,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, Ecal_3p1pi,p_miss_perp_3p1pi,&N_p_three);

 if(N_p_three!=0){
	       P_3p1pito3p0pi[0]= (N_3p0pi/N_all)*(N_p1[0]/N_p_three);
	       P_3p1pito3p0pi[1]= (N_3p0pi/N_all)*(N_p1[1]/N_p_three);
	       P_3p1pito3p0pi[2]= (N_3p0pi/N_all)*(N_p1[2]/N_p_three);
	     }


 for(int z=0;z<N_3prot;z++){
 //---------------------------------- 3p 1pi ->1p 0pi   ----------------------------------------------

   P_3p1pito1p0pi[z]= -N_1p0pi[z]/N_all;

 //---------------------------------- 3p 1pi ->1p 1pi   ----------------------------------------------

   N_pidet=N_piundet=0;
   prot1_pi1_rot_func(fbeam_en, V3_q,V3_3prot_uncorr[z],V3_pi, q_pi,radstat, &N_pidet,&N_piundet);
   if(N_pidet!=0) P_3p1pito1p1pi[z]=(N_1p1pi[z]/N_all)*(N_piundet/N_pidet);

 //---------------------------------- 3p 1pi ->2p 0pi   ----------------------------------------------

		      for(int i=0;i<N_3prot;i++){       //looping through 2p combinations  out of 3p
			if(z!=i && z<i){               // 3 pairs of 2proton combinations with z, i indexes(z<i)

			  V3_2p_corr[0]=V3_3prot_corr[z];V3_2p_corr[1]=V3_3prot_corr[i];
			  V3_2p_uncorr[0]=V3_3prot_uncorr[z];V3_2p_uncorr[1]=V3_3prot_uncorr[i];

			  P_2pto1p[0]=0;P_2pto1p[1]=0;N_2p_det=0;
			  prot2_rot_func(fbeam_en, V3_q, V3_2p_corr,V3_2p_uncorr, V4_el,Ecal_2p,p_miss_perp_2p,P_2pto1p ,&N_2p_det);
			  if(N_2p_det!=0){
			    P_3p1pito2p0pi[z]=P_3p1pito2p0pi[z]+(N_2p0pi[count]/N_all)*P_2pto1p[0];
			    P_3p1pito2p0pi[i]=P_3p1pito2p0pi[i]+(N_2p0pi[count]/N_all)*P_2pto1p[1];
			  }


      //---------------------------------- 3p 1pi ->3p 0pi->2p0pi   ----------------------------------------------

	     if(N_p_three!=0){
	       P_3p1pito3p0pi[z]= P_3p1pito3p0pi[z]+(N_3p0pi/N_all)*(-P_3pto2p[count][0]);
	       P_3p1pito3p0pi[i]= P_3p1pito3p0pi[i]+(N_3p0pi/N_all)*(-P_3pto2p[count][1]);
	     }

 //---------------------------------- 3p 1pi ->2p 1pi   ----------------------------------------------
 P_2p1pito2p0pi[0]=P_2p1pito2p0pi[1]=0; P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0; P_2p1pito1p0pi[0]=P_2p1pito1p0pi[1]=0;Ptot=0;

 prot2_pi1_rot_func(fbeam_en, V3_q,V3_2p_corr,V3_2p_uncorr,V3_pi, q_pi,radstat,V4_el,Ecal_2p1pi,p_miss_perp_2p1pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

 // P_3p1pito2p1pi[z]=(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[z]- P_2p1pito1p1pi[z]+P_2p1pito1p0pi[z]);
 // P_3p1pito2p1pi[i]=(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[i]- P_2p1pito1p1pi[i]+P_2p1pito1p0pi[i]);

  P_3p1pito2p1pi[z]= P_3p1pito2p1pi[z]+(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[0]- P_2p1pito1p1pi[0]+P_2p1pito1p0pi[0]);
  P_3p1pito2p1pi[i]= P_3p1pito2p1pi[i]+(N_2p1pi[count]/N_all)*(-P_2p1pito2p0pi[1]- P_2p1pito1p1pi[1]+P_2p1pito1p0pi[1]);

	     count=count+1;
			}
		      }

		      // P_tot_3p[z]=P_3p1pito2p1pi[z]+P_3p1pito3p0pi[z]+P_3p1pito2p0pi[z]+P_3p1pito1p1pi[z]+ P_3p1pito1p0pi[z];

    }//looping through 3p

 P_tot_3p[0]=P_3p1pito2p1pi[0]+P_3p1pito3p0pi[0]+P_3p1pito2p0pi[0]+P_3p1pito1p1pi[0]+ P_3p1pito1p0pi[0];
 P_tot_3p[1]=P_3p1pito2p1pi[1]+P_3p1pito3p0pi[1]+P_3p1pito2p0pi[1]+P_3p1pito1p1pi[1]+ P_3p1pito1p0pi[1];
 P_tot_3p[2]=P_3p1pito2p1pi[2]+P_3p1pito3p0pi[2]+P_3p1pito2p0pi[2]+P_3p1pito1p1pi[2]+ P_3p1pito1p0pi[2];

    }

if(N_all==0){
  P_tot_3p[0]= P_tot_3p[1]=P_tot_3p[2]=0;
  }

}









void FilterData::pi1_rot_func(std::string beam_en, TVector3 V3_pi, int q_pi,bool radstat,TVector3 V3_q,double *P_pi){

  std::string fbeam_en = beam_en;
  double N_pion=0;
  double rot_angle;
  TVector3 V3_rot_pi;
  Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;

   if(!radstat){
  for(int g=0; g<N_tot; g++){

    V3_rot_pi=V3_pi;
    rot_angle=gRandom->Uniform(0,2*TMath::Pi());
    V3_rot_pi.Rotate(rot_angle,V3_q);
    if(Pi_phot_fid_united(fbeam_en, V3_rot_pi,q_pi)) N_pion=N_pion+1;
  }
    }
    else {
     N_pion=N_tot;
   }

  if(N_pion!=0)     *P_pi=(N_tot-N_pion)/N_pion;
  else *P_pi=0;

  // if(!radstat) cout<<"nereqev     "<<N_pion<<"   radstat"<<"     "<<radstat<<endl;
}




void FilterData::pi2_rot_func(std::string beam_en, TVector3 V3_pi[2], int q_pi[2],bool radstat[2], TVector3 V3_q,double *P_0pi,double P_1pi[2]){

  std::string fbeam_en = beam_en;
  const int N_pi=2;
  TVector3 V3_rot_pi[N_pi];
 double rot_angle;
 bool status_pi[N_pi]={true};
 Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;
 double N_bothpi=0,N_nopi=0,N_1pi[N_pi]={0},P_pi1[N_pi]={0};

    if(!radstat[0] || !radstat[1]){
     for(int g=0; g<N_tot; g++){

       rot_angle=gRandom->Uniform(0,2*TMath::Pi());
       for(int i=0;i<N_pi;i++){
	 if(!radstat[i]){
	   V3_rot_pi[i]=V3_pi[i];
	   V3_rot_pi[i].Rotate(rot_angle,V3_q);
	   status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
	 }
       }
       if( status_pi[0] && !status_pi[1]) N_1pi[0]=N_1pi[0]+1;
       if(!status_pi[0] &&  status_pi[1]) N_1pi[1]=N_1pi[1]+1;
       if(!status_pi[0] && !status_pi[1]) N_nopi=N_nopi+1;
       if( status_pi[0] &&  status_pi[1]) N_bothpi=N_bothpi+1;
     }
    }
    else {N_bothpi=N_tot;}

    pi1_rot_func(fbeam_en, V3_pi[0],  q_pi[0],radstat[0], V3_q,P_pi1);
    pi1_rot_func(fbeam_en, V3_pi[1],  q_pi[1],radstat[1], V3_q,P_pi1+1);

      if(N_bothpi!=0){
	*P_0pi=N_nopi/N_bothpi;
	P_1pi[0]=N_1pi[0]/N_bothpi*P_pi1[0];
	P_1pi[1]=N_1pi[1]/N_bothpi*P_pi1[1];
      }
      else{
	*P_0pi=0;
	P_1pi[0]=0;
	P_1pi[1]=0;
      }
}





void FilterData::pi3_rot_func(std::string beam_en, TVector3 V3_pi[3], int q_pi[3],bool radstat[3], TVector3 V3_q,double *P_0pi,double P_1pi[3],double P_320[3],double P_3210[][2]){

 std::string fbeam_en = beam_en;
 const int N_pi=3;
  TVector3 V3_rot_pi[N_pi];
 double rot_angle;
 bool status_pi[N_pi]={true};
 Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;
 double N_1pi[N_pi]={0},N_allpi=0,N_nopi=0,N_2pi[N_pi]={0};


 if(!radstat[0] || !radstat[1]  || !radstat[2]){
  for(int g=0; g<N_tot; g++){

         rot_angle=gRandom->Uniform(0,2*TMath::Pi());

       for(int i=0;i<N_pi;i++){

	 if(!radstat[i]){
	 V3_rot_pi[i]=V3_pi[i];
	 V3_rot_pi[i].Rotate(rot_angle,V3_q);
	 status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
	 }
       }

       if( status_pi[0]  && !status_pi[1] &&  !status_pi[2]) N_1pi[0]=N_1pi[0]+1;
       if(!status_pi[0] &&   status_pi[1]  && !status_pi[2]) N_1pi[1]=N_1pi[1]+1;
       if(!status_pi[0] &&  !status_pi[1] &&   status_pi[2]) N_1pi[2]=N_1pi[2]+1;
       if( status_pi[0]  &&  status_pi[1] &&  !status_pi[2]) N_2pi[0]=N_2pi[0]+1;
       if( status_pi[0] &&  !status_pi[1]  &&  status_pi[2]) N_2pi[1]=N_2pi[1]+1;
       if(!status_pi[0] &&   status_pi[1] &&   status_pi[2]) N_2pi[2]=N_2pi[2]+1;
       if(!status_pi[0] &&  !status_pi[1] &&  !status_pi[2]) N_nopi=N_nopi+1;
       if( status_pi[0]  &&  status_pi[1]  &&  status_pi[2]) N_allpi=N_allpi+1;
     }
 }
 else N_allpi=N_tot;

  const int N_pi2=2;
  double P_pi=0;
  double P_1pion[N_pi2]={0},P_0pion=0;
  TVector3 V3_pion[N_pi2];
  int q_pion[N_pi2];
  bool rad_stat[N_pi2];

  if(N_allpi!=0){
 //---------------------------3pi->0pi----------------------------------------------
    *P_0pi=N_nopi/N_allpi;
 //---------------------------3pi->1pi->0pi----------------------------------------------
    for(int h=0;h<N_pi;h++){
      pi1_rot_func(fbeam_en, V3_pi[h],q_pi[h],radstat[h],V3_q,&P_pi);
      P_1pi[h]=P_pi*(N_1pi[h]/N_allpi);
    //---------------------------3pi->2pi->0pi----------------------------------------------

      if(h==0)   {
	V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[1];
	q_pion[0]=q_pi[0];q_pion[1]=q_pi[1];
	rad_stat[0]=radstat[0];rad_stat[1]=radstat[1];
      }
      if(h==1)   {
	V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[2];
	q_pion[0]=q_pi[0];q_pion[1]=q_pi[2];
	rad_stat[0]=radstat[0];rad_stat[1]=radstat[2];
      }
      if(h==2)   {
	V3_pion[0]=V3_pi[1];V3_pion[1]=V3_pi[2];
	q_pion[0]=q_pi[1];q_pion[1]=q_pi[2];
	rad_stat[0]=radstat[1];rad_stat[1]=radstat[2];
      }
      pi2_rot_func(fbeam_en, V3_pion,q_pion,rad_stat,V3_q,&P_0pion, P_1pion);
      P_320[h]=P_0pion*(N_2pi[h]/N_allpi);

//---------------------------3pi->2pi->1pi->0pi----------------------------------------------

    P_3210[h][0]=P_1pion[0]*(N_2pi[h]/N_allpi);
    P_3210[h][1]=P_1pion[1]*(N_2pi[h]/N_allpi);

    }//end of 3p loop


  }// end of N_allpi!=0 statement


  else{
    for(int h=0;h<3;h++){
      P_320[h]=0;
      P_3210[h][0]=0;
      P_3210[h][1]=0;
    }
  }


}







void FilterData::pi4_rot_func(std::string beam_en, TVector3 V3_pi[4], int q_pi[4],bool radstat[4], TVector3 V3_q,double *P_0pi,double *P_410,double *P_420,double *P_4210,double *P_430,double *P_4310,double *P_4320,double *P_43210){

 std::string fbeam_en = beam_en;
 const int N_pi=4;
  TVector3 V3_rot_pi[N_pi];
 double rot_angle;
 bool status_pi[N_pi]={true};
 double N_1pi[N_pi]={0},N_allpi=0,N_nopi=0,N_2pi[6]={0},N_3pi[4]={0};


if(!radstat[0] || !radstat[1]  || !radstat[2]|| !radstat[3]){
  for(int g=0; g<N_tot; g++){

         rot_angle=gRandom->Uniform(0,2*TMath::Pi());

       for(int i=0;i<N_pi;i++){

	 if(!radstat[i]){
	   V3_rot_pi[i]=V3_pi[i];
	   V3_rot_pi[i].Rotate(rot_angle,V3_q);
	   status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
	 }
       }

       if( status_pi[0]  && !status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_1pi[0]=N_1pi[0]+1; //1pi or phot
       if( !status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_1pi[1]=N_1pi[1]+1;
       if( !status_pi[0]  && !status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_1pi[2]=N_1pi[2]+1;
       if( !status_pi[0]  && !status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_1pi[3]=N_1pi[3]+1;

       if( status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_2pi[0]=N_2pi[0]+1;//2pi or phot
       if( status_pi[0]  &&! status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_2pi[1]=N_2pi[1]+1;
       if( status_pi[0]  &&! status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_2pi[2]=N_2pi[2]+1;
       if( !status_pi[0]  && status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_2pi[3]=N_2pi[3]+1;
       if( !status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_2pi[4]=N_2pi[4]+1;
       if( !status_pi[0]  && !status_pi[1] &&  status_pi[2]  &&  status_pi[3]) N_2pi[5]=N_2pi[5]+1;

       if( status_pi[0]  && status_pi[1] &&  status_pi[2]  &&  !status_pi[3]) N_3pi[0]=N_3pi[0]+1;//3pi or phot
       if( status_pi[0]  && status_pi[1] &&  !status_pi[2]  &&  status_pi[3]) N_3pi[1]=N_3pi[1]+1;
       if( status_pi[0]  && !status_pi[1] &&  status_pi[2]  &&  status_pi[3]) N_3pi[2]=N_3pi[2]+1;
       if( !status_pi[0]  && status_pi[1] &&  status_pi[2]  &&  status_pi[3]) N_3pi[3]=N_3pi[3]+1;

       if( !status_pi[0]  && !status_pi[1] &&  !status_pi[2]  &&  !status_pi[3]) N_nopi=N_nopi+1; //0 pi or phot
       if( status_pi[0]  && status_pi[1] &&  status_pi[2]  && status_pi[3]) N_allpi=N_allpi+1; //4pi or phot
     }
 }
  else N_allpi=N_tot;

   double P_pi=0;
  const int N_pi3=3;
  TVector3 V3_pion[N_pi3];
  double P_1pion[N_pi3]={0},P_0pion=0, P_320_pion[3]={0}, P_3210_pion[3][2]={0};
  int q_pion[N_pi3];
  bool rad_stat[N_pi3];

  if(N_allpi!=0){
 //---------------------------4pi->0pi----------------------------------------------
    *P_0pi=N_nopi/N_allpi;
 //---------------------------4pi->1pi->0pi----------------------------------------------
    for(int h=0;h<N_pi;h++){
      pi1_rot_func(fbeam_en, V3_pi[h],q_pi[h],radstat[h],V3_q,&P_pi);
      *P_410=*P_410+P_pi*(N_1pi[h]/N_allpi);

    //---------------------------4pi->3pi->0pi----------------------------------------------
      if(h==0)   {
	V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[1];V3_pion[2]=V3_pi[2];
	q_pion[0]=q_pi[0];q_pion[1]=q_pi[1];q_pion[2]=q_pi[2];
	rad_stat[0]=radstat[0];rad_stat[1]=radstat[1];rad_stat[2]=radstat[2];
      }
      if(h==1)   {
	V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[1];V3_pion[2]=V3_pi[3];
	q_pion[0]=q_pi[0];q_pion[1]=q_pi[1];q_pion[2]=q_pi[3];
	rad_stat[0]=radstat[0];rad_stat[1]=radstat[1];rad_stat[2]=radstat[3];
      }
      if(h==2)   {
	V3_pion[0]=V3_pi[0];V3_pion[1]=V3_pi[2];V3_pion[2]=V3_pi[3];
	q_pion[0]=q_pi[0];q_pion[1]=q_pi[2];q_pion[2]=q_pi[3];
	rad_stat[0]=radstat[0];rad_stat[1]=radstat[2];rad_stat[2]=radstat[3];
      }
      if(h==3)   {
	V3_pion[0]=V3_pi[1];V3_pion[1]=V3_pi[2];V3_pion[2]=V3_pi[3];
	q_pion[0]=q_pi[1];q_pion[1]=q_pi[2];q_pion[2]=q_pi[3];
	rad_stat[0]=radstat[1];rad_stat[1]=radstat[2];rad_stat[2]=radstat[3];
      }

      pi3_rot_func(fbeam_en, V3_pion,q_pion,rad_stat,V3_q,&P_0pion, P_1pion,P_320_pion,P_3210_pion);
      *P_430=*P_430+P_0pion*(N_3pi[h]/N_allpi);

//---------------------------4pi->3pi->1pi->0pi----------------------------------------------

      *P_4310= *P_4310+(P_1pion[0]+P_1pion[1]+P_1pion[2])*(N_3pi[h]/N_allpi);

//---------------------------4pi->3pi->2pi->0pi----------------------------------------------

      *P_4320=*P_4320+(P_320_pion[0]+P_320_pion[1]+P_320_pion[2])*(N_3pi[h]/N_allpi);

//---------------------------4pi->3pi->2pi->1pi->0pi----------------------------------------------

      *P_43210=*P_43210+((P_3210_pion[0][0]+P_3210_pion[0][1])+(P_3210_pion[1][0]+P_3210_pion[1][1])+
			 (P_3210_pion[2][0]+P_3210_pion[2][1]))*(N_3pi[h]/N_allpi);

    }//end of 4pi loop


//---------------------------4pi->2pi->0pi----------------------------------------------
  const int N2pi=2;
  TVector3 V3_2pi[N2pi];
  int q_2pi[N2pi];
  double P_0pi=0, P_1pi[N2pi]={0};
  bool radstat_2pi[N2pi];

    for(int h=0;h<6;h++){

 if(h==0)   {
       V3_2pi[0]=V3_pi[0];      V3_2pi[1]=V3_pi[1];
         q_2pi[0]=q_pi[0];        q_2pi[1]=q_pi[1];
	 radstat_2pi[0]=radstat[0];  radstat_2pi[1]=radstat[1];
      }
 if(h==1)   {
       V3_2pi[0]=V3_pi[0];      V3_2pi[1]=V3_pi[2];
         q_2pi[0]=q_pi[0];        q_2pi[1]=q_pi[2];
	 radstat_2pi[0]=radstat[0];  radstat_2pi[1]=radstat[2];
      }
 if(h==2)   {
       V3_2pi[0]=V3_pi[0];      V3_2pi[1]=V3_pi[3];
         q_2pi[0]=q_pi[0];        q_2pi[1]=q_pi[3];
	 radstat_2pi[0]=radstat[0];  radstat_2pi[1]=radstat[3];
      }
 if(h==3)   {
       V3_2pi[0]=V3_pi[1];      V3_2pi[1]=V3_pi[2];
         q_2pi[0]=q_pi[1];        q_2pi[1]=q_pi[2];
	 radstat_2pi[0]=radstat[1];  radstat_2pi[1]=radstat[2];
      }
 if(h==4)   {
       V3_2pi[0]=V3_pi[1];      V3_2pi[1]=V3_pi[3];
         q_2pi[0]=q_pi[1];        q_2pi[1]=q_pi[3];
	 radstat_2pi[0]=radstat[1];  radstat_2pi[1]=radstat[3];
      }
 if(h==5)   {
       V3_2pi[0]=V3_pi[2];      V3_2pi[1]=V3_pi[3];
         q_2pi[0]=q_pi[2];        q_2pi[1]=q_pi[3];
	 radstat_2pi[0]=radstat[2];  radstat_2pi[1]=radstat[3];
      }

 pi2_rot_func(fbeam_en, V3_2pi,q_2pi,radstat_2pi,V3_q,&P_0pi, P_1pi);

 *P_420=*P_420+P_0pi*(N_2pi[h]/N_allpi);

//---------------------------4pi->2pi->1pi->0pi----------------------------------------------

 *P_4210=*P_4210+(P_1pi[0]+P_1pi[0])*(N_2pi[h]/N_allpi);

    }

  }// end of N_allpi!=0 statement

 else{
   *P_0pi=0;
   *P_410=0;
   *P_430=0;
   *P_4310=0;
   *P_4320=0;
   *P_43210=0;
   *P_420=0;
   *P_4210=0;
 }



}
