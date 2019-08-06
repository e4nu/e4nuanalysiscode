
#define E2A_EP_C

#include "e2a_ep_neutrino6_united4_radphot.h"
#include "Constants.h"
#include "Fiducial.h"
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
using namespace Fiducial;

 double e_mass=0.000510998;
int fTorusCurrent;
bool SCpdcut=true;

std::string target_name;
std::map<std::string,double> en_beam;
std::map<std::string,double> en_beam_Ecal;
std::map<std::string,double> en_beam_Eqe;


void SetFiducialCutParameters(std::string beam_en); // Load Fidicual Parameters for 1.1 and 4.4 GeV from file
//void SetMomCorrParameters();


TLorentzVector EMomentumCorrection(TLorentzVector V4el);

TF1 *vz_corr_func;
double vz_corr(double phi,double theta);
TVector3 FindUVW(TVector3 xyz);
Bool_t CutUVW(TVector3 ecxyz);
Bool_t EFiducialCut(std::string beam_en, TVector3 momentum);
Bool_t PFiducialCut(std::string beam_en, TVector3 momentum);
Float_t ProtonMomCorrection_He3_4Cell(std::string atarget, TLorentzVector V4Pr, Float_t vertex_p);
Bool_t PiplFiducialCut(std::string beam_en, TVector3 momentum,Float_t *philow,Float_t *phiup);
Bool_t PimiFiducialCut(std::string beam_en, TVector3 momentum, Float_t *pimi_philow, Float_t *pimi_phiup);
bool Phot_fid(TVector3 V3_phot);
bool Pi_phot_fid_united(std::string beam_en, TVector3 V3_pi_phot, int q_pi_phot);
Bool_t GetEPhiLimits(std::string beam_en, Float_t momentum, Float_t theta, Int_t sector,
		     Float_t *EPhiMin, Float_t *EPhiMax);
void prot3_rot_func(std::string beam_en, TVector3 V3q, TVector3 V3prot[3],TVector3  V3prot_uncorr[3],TLorentzVector V4el,double Ecal_3pto2p[][2],double  pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3],double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det);
void prot2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot[2],TVector3  V3prot_uncorr[2],TLorentzVector V4el,double Ecal_2pto1p[2],double  pmiss_perp_2pto1p[2],double  P2pto1p[2], double *Nboth);
void prot1_pi1_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi, int q_pi,bool radstat, double *N_pi_p,double *N_nopi_p);
void prot1_pi2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[2], int q_pi[2],bool radstat[2], double *P_1p0pi,double P_1p1pi[2]);
void prot1_pi3_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[3], int q_pi[3],bool radstat[3],double *P_tot);
void prot2_pi1_rot_func(std::string beam_en, TVector3 V3q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_2p1pi_to2p0pi[2],double p_miss_perp_2p1pi_to2p0pi[2],double P_2p1pito2p0pi[2],double P_2p1pito1p1pi[2],double P_2p1pito1p0pi[2],double *P_tot);
void prot2_pi2_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], int q_pi[2],bool radstat[2],TLorentzVector V4_el, double Ecal_2p2pi[2],double p_miss_perp_2p2pi[2],double P_tot_2p[2]);
void prot3_pi1_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_3p1pi[3],double p_miss_perp_3p1pi[3],double P_tot_3p[3]);
void pi1_rot_func(std::string beam_en, TVector3 V3_pi, int q_pi,bool radstat,TVector3 V3_q,double *P_pi);
void pi2_rot_func(std::string beam_en, TVector3 V3_pi[2], int q_pi[2],bool radstat[2], TVector3 V3_q,double *P_0pi,double P_1pi[2]);
void pi3_rot_func(std::string beam_en, TVector3 V3_pi[3], int q_pi[3],bool radstat[3], TVector3 V3_q,double *P_0pi,double P_1pi[3],double P_320[3],double P_3210[][2]);
void pi4_rot_func(std::string beam_en, TVector3 V3_pi[4], int q_pi[4],bool radstat[4], TVector3 V3_q,double *P_0pi,double *P_410,double *P_420,double *P_4210,double *P_430,double *P_4310,double *P_4320,double *P_43210);

TF1 *pipl_deltat_sig,*pipl_deltat_mean,*pimi_deltat_sig,*pimi_deltat_mean, *fsum_pimi,*fsub_pimi, *fsum_pipl,*fsub_pipl, *prot_deltat_sig, *prot_deltat_mean,*fsum_prot,*fsub_prot,*el_Epratio_sig,*el_Epratio_mean,*fsum_e,*fsub_e;
TF1 *up_lim1_ec, *up_lim2_ec,*up_lim3_ec,*up_lim4_ec, *up_lim5_ec,*up_lim6_ec,*low_lim1_ec,*low_lim2_ec,*low_lim3_ec, *low_lim4_ec,*low_lim5_ec,*low_lim6_ec;
TF1 *leftside_lim1_ec, *leftside_lim2_ec,*leftside_lim3_ec, *leftside_lim4_ec,*leftside_lim5_ec, *leftside_lim6_ec,*rightside_lim1_ec, *rightside_lim2_ec,*rightside_lim3_ec, *rightside_lim4_ec,*rightside_lim5_ec, *rightside_lim6_ec;


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
  if(x[0]<par[1]) return el_Epratio_mean->EvalPar(x)+par[0]*el_Epratio_sig->EvalPar(x);
if(x[0]>=par[1]) return el_Epratio_mean->Eval(par[1])+par[0]*el_Epratio_sig->Eval(par[1]);
}

Double_t FSub_e(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return el_Epratio_mean->EvalPar(x)-par[0]*el_Epratio_sig->EvalPar(x);
if(x[0]>=par[1]) return el_Epratio_mean->Eval(par[1])-par[0]*el_Epratio_sig->Eval(par[1]);
}




//proton Delta_t vs momentum PID cut

Double_t FSum_prot(Double_t *x, Double_t *par){   //the 2 parameters are the cut range and momentum limit
  if(x[0] < par[1])    return prot_deltat_mean->EvalPar(x)+par[0]*prot_deltat_sig->EvalPar(x);
  if(x[0] >= par[1])   return prot_deltat_mean->Eval(par[1])+par[0]*prot_deltat_sig->Eval(par[1]);
}
Double_t FSub_prot(Double_t *x,Double_t *par){
  if(x[0] < par[1])    return prot_deltat_mean->EvalPar(x)-par[0]*prot_deltat_sig->EvalPar(x);
  if(x[0] >= par[1])   return prot_deltat_mean->Eval(par[1])-par[0]*prot_deltat_sig->Eval(par[1]);
}


//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions

Double_t FSum_pimi(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pimi_deltat_mean->EvalPar(x)+par[0]*pimi_deltat_sig->EvalPar(x);
  if(x[0]>=par[1])return pimi_deltat_mean->Eval(par[1])+par[0]*pimi_deltat_sig->Eval(par[1]);
}
Double_t FSub_pimi(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pimi_deltat_mean->EvalPar(x)-par[0]*pimi_deltat_sig->EvalPar(x);
  if(x[0]>=par[1])return pimi_deltat_mean->Eval(par[1])-par[0]*pimi_deltat_sig->Eval(par[1]);
}



//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions
Double_t FSum_pipl(Double_t *x,Double_t *par){

  if(x[0]<par[1])  return pipl_deltat_mean->EvalPar(x)+par[0]*pipl_deltat_sig->EvalPar(x);
  if(x[0]>=par[1]) return pipl_deltat_mean->Eval(par[1])+par[0]*pipl_deltat_sig->Eval(par[1]);
}
Double_t FSub_pipl(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pipl_deltat_mean->EvalPar(x)-par[0]*pipl_deltat_sig->EvalPar(x);
  if(x[0]>=par[1])return pipl_deltat_mean->Eval(par[1])-par[0]*pipl_deltat_sig->Eval(par[1]);
}


 const int N_tot=10;




void e2a_ep_neutrino6_united4_radphot::Loop()
{

  target_name = ftarget;   //std string for target name
  cout << fbeam_en << endl;
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
  nentries =8000000;

  const double m_pimi= 0.139570,m_pipl= 0.139570,Mpi = 0.139570;
  const double m_prot= 0.9382720813,m_neut=0.939565,H3_bind_en=0.008481,He4_bind_en=0.0283,C12_bind_en=0.09215, B_bind_en=0.0762,He3_bind_en=0.0077,D2_bind_en=0.00222,Fe_bind_en=0.49226,Mn_bind_en=0.4820764;
  double N_prot1=0,N_prot2=0,N_prot_both=0;
  double eps;
  Float_t cphil = 0, pimi_phimin, pimi_phimax;
  Float_t cphir = 0;
  double ece;
  const Double_t c = 2.99792E+10;
  Double_t ns_to_s = 1.0E-9;
  const int n_slice=3,nsect=6;
  double bett,deltt;
 const double pperp_min[n_slice]={0.,0.2,0.4};
 const double pperp_max[n_slice]={0.2,0.4,10.};
  TVector3 V3_pimi,V3_pipl,V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;
  int N_comb=3;
  double sum_val,sub_val;
  double epratio_sig_cutrange=3.;
  double prot_delt_cutrange=3.;
  int 	el_segment, el_cc_sector;
  double delt_uplim,delt_lowlim;
 double prot_accept_mom_lim;
  double prot_mom_lim;
  double min_good_mom;
  double max_mom;
  Double_t el_sccc_timediff, sc_cc_delt_cut_sect[nsect]={-2,-5,-8,-8,-2,2}, sc_cc_delt_cut_sect_max[nsect]={15,15,10,15,20,20}, el_cc_nphe,elmom_corr_fact[nsect];
  double pipl_maxmom, pimi_maxmom,pimi_delt_cutrange,pipl_delt_cutrange;
 int N_pperp,N_Ecal;
  double *pperp_cut,*Ecal_lowlim,*Ecal_uplim;




if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.)
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
      N_pperp=2,N_Ecal=6;
      pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;  pperp_cut[1]=0.2;
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

 elmom_corr_fact[0]=1.007;elmom_corr_fact[1]=0.988;elmom_corr_fact[2]=1.008;elmom_corr_fact[3]=1.011;elmom_corr_fact[4]=1.014;elmom_corr_fact[5]=1.013;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
    }




  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.)
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
      N_pperp=2,N_Ecal=6;
      pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;  pperp_cut[1]=0.2;
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


    //map<std::pair<std::string, int>, double> EC_time_offset;
 EC_time_offset[std::make_pair("3He",1)]=-1.37;  EC_time_offset[std::make_pair("3He",2)]=-1.42; EC_time_offset[std::make_pair("3He",3)]=-1.55;
 EC_time_offset[std::make_pair("3He",4)]=-1.53;  EC_time_offset[std::make_pair("3He",5)]=-1.49; EC_time_offset[std::make_pair("3He",6)]=-1.44;

 EC_time_offset[std::make_pair("4He",1)]=0.72;  EC_time_offset[std::make_pair("4He",2)]=0.27; EC_time_offset[std::make_pair("4He",3)]=0.16;
 EC_time_offset[std::make_pair("4He",4)]=0.21;  EC_time_offset[std::make_pair("4He",5)]=0.22; EC_time_offset[std::make_pair("4He",6)]=0.21;

 EC_time_offset[std::make_pair("C12",1)]=0.50;  EC_time_offset[std::make_pair("C12",2)]=0.39; EC_time_offset[std::make_pair("C12",3)]=0.29;
 EC_time_offset[std::make_pair("C12",4)]=0.29;  EC_time_offset[std::make_pair("C12",5)]=0.32; EC_time_offset[std::make_pair("C12",6)]=0.33;

 EC_time_offset[std::make_pair("56Fe",1)]=0.75;  EC_time_offset[std::make_pair("56Fe",2)]=0.49; EC_time_offset[std::make_pair("56Fe",3)]=0.37;
 EC_time_offset[std::make_pair("56Fe",4)]=0.39;  EC_time_offset[std::make_pair("56Fe",5)]=0.43; EC_time_offset[std::make_pair("56Fe",6)]=0.44;

 elmom_corr_fact[0]=1.001;elmom_corr_fact[1]=0.991;elmom_corr_fact[2]=1.005;elmom_corr_fact[3]=1.004;elmom_corr_fact[4]=1.006;elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
    }


  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5)
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
  N_pperp=2,N_Ecal=6;
  pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;  pperp_cut[1]=0.2;
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

 elmom_corr_fact[0]=1.001;elmom_corr_fact[1]=0.991;elmom_corr_fact[2]=1.005;elmom_corr_fact[3]=1.004;elmom_corr_fact[4]=1.006;elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
    }

  Ecal_offset["3He"]=0.004;
  Ecal_offset["4He"]=0.005;
  Ecal_offset["C12"]=0.005;
  Ecal_offset["56Fe"]=0.011;

  bind_en["3He"]=He3_bind_en-D2_bind_en+Ecal_offset["3He"]; //the offset is used to shift the peak to be at 0
  bind_en["4He"]=He4_bind_en-H3_bind_en+Ecal_offset["4He"];
  bind_en["C12"]=C12_bind_en-B_bind_en+Ecal_offset["C12"];
  bind_en["56Fe"]=Fe_bind_en-Mn_bind_en+Ecal_offset["56Fe"];
  bind_en["CH2"]=C12_bind_en-B_bind_en;

  target_mass["3He"]=2*m_prot+m_neut-He3_bind_en;
  target_mass["4He"]=2*m_prot+2*m_neut-He4_bind_en;
  target_mass["C12"]=6*m_prot+6*m_neut-C12_bind_en;
  target_mass["56Fe"]=26*m_prot+30*m_neut-Fe_bind_en;
  target_mass["CH2"]=6*m_prot+6*m_neut-C12_bind_en;

  residual_target_mass["3He"]=m_prot+m_neut-D2_bind_en;
  residual_target_mass["4He"]=m_prot+2*m_neut-H3_bind_en;
  residual_target_mass["C12"]=5*m_prot+6*m_neut-B_bind_en;
  residual_target_mass["56Fe"]=25*m_prot+30*m_neut-Mn_bind_en;
  residual_target_mass["CH2"]=25*m_prot+30*m_neut-Mn_bind_en;


  TH1F *h1_el_cc_deltat[nsect], *h1_el_cc_deltat_cut[nsect],*h1_Erec_p_bkgd_slice[n_slice],*h1_Etot_p_bkgd_slice[n_slice],*h1_Etot_Npi0[n_slice],*h1_Erec_Npi0_new[n_slice],* h1_Erec_bkgd_pipl_pimi_new_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pipl[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pimi[n_slice],*h1_Etot_Npi1[n_slice],*h1_Erec_Npi1[n_slice],*h1_Etot_bkgd_1p2pi[n_slice],*h1_Erec_bkgd_1p2pi[n_slice],*h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[n_slice],*h1_Etot_p_bkgd_slice_2p2pi[n_slice],*h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[n_slice],*h1_Erec_p_bkgd_slice_2p2pi[n_slice],*h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[n_slice],*h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[n_slice],*h1_Etot_bkgd_1p2pi_1p0pi[n_slice],*h1_Erec_bkgd_1p2pi_1p0pi[n_slice],*h1_Etot_bkgd_1p3pi[n_slice],*h1_Erec_bkgd_1p3pi[n_slice],*h1_e_mom_corrfuct[nsect];
  TH1F *h_Etot_piplpimi_subtruct_fact[n_slice],*h_Erec_piplpimi_subtruct_new_fact[n_slice],*h1_Etot_p_bkgd_slice_sub[n_slice],*h1_Erec_p_bkgd_slice_sub[n_slice],*h1_Etot_3pto1p_slice[n_slice],*h1_Erec_3pto1p_slice[n_slice],*h1_Etot_3pto2p_slice[n_slice],*h1_Erec_3pto2p_slice[n_slice],*h1_Etot_3p1pi_slice[n_slice],*h1_Erec_3p1pi_slice[n_slice],*h1_Etot_4pto3p_slice[n_slice],*h1_Erec_4pto3p_slice[n_slice],*h1_Etot_4pto1p_slice[n_slice],*h1_Erec_4pto1p_slice[n_slice],*h1_Etot_4pto2p_slice[n_slice],*h1_Erec_4pto2p_slice[n_slice],*h1_Etot_43pto1p_slice[n_slice],*h1_Erec_43pto1p_slice[n_slice],*h1_Etot_p_bkgd_slice_sub1p2pi[n_slice], *h1_Erec_p_bkgd_slice_sub1p2pi[n_slice],*h1_Etot_p_bkgd_slice_sub2p1pi_1p[n_slice],*h1_Etot_p_bkgd_slice_sub3p1pi_0p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p[n_slice],*h1_Erec_p_bkgd_slice_sub3p1pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub2p1pi_2p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_2p[n_slice],*h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p0pi[n_slice],*h1_Etot_p_bkgd_slice_sub1p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p2pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub1p3pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub3p1pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p3pi_0pi[n_slice],*h1_Etot_p_bkgd_slice_sub2p2pi_0pi[n_slice],*h1_Erec_p_bkgd_slice_sub2p2pi_0pi[n_slice];
  TH1F *h1_Etot_p_bkgd_slice_sub32[n_slice],*h1_Erec_p_bkgd_slice_sub32[n_slice],*h1_Etot_p_bkgd_slice_sub31[n_slice],*h1_Erec_p_bkgd_slice_sub31[n_slice],*h1_Etot_p_bkgd_slice_sub43[n_slice],*h1_Erec_p_bkgd_slice_sub43[n_slice],*h1_Etot_p_bkgd_slice_sub41[n_slice],*h1_Erec_p_bkgd_slice_sub41[n_slice],*h1_Erec_p_bkgd_slice_sub42[n_slice],*h1_Etot_p_bkgd_slice_sub42[n_slice],*h1_Etot_p_bkgd_slice_sub431[n_slice],*h1_Erec_p_bkgd_slice_sub431[n_slice];
  TH1F *h1_el_ec_sc_timediff_sect_corr[nsect], *h1_el_ec_sc_timediff_sect[nsect],*h1_beta_ec_corr_sect[nsect],*h1_el_SCpdfidcut[6],*h1_el_SCpd[6],*pos_m_slices[4],*prot_m_slices[4],*pipl_m_slices[4];
  TH2F *h2_el_theta_phi_p_beffidcut[6],*h2_el_theta_phi_p_fidcut[6],*h2_el_ec_sc_timediff_ecu[6],*h2_el_ec_sc_timediff_ecv[6],*h2_el_ec_sc_timediff_ecw[6], *h2_el_ec_sc_timediff_SCpd[6],*h2_prot_theta_phi_p_beffidcut[6],*h2_prot_theta_phi_p_fidcut[6],*h2_el_theta_phi_p_beffidcut2[6],*h2_el_theta_phi_p_fidcut2[6],*h2_N_pi_phot[20],*h2_pimi_theta_phi_p_beffidcut[6],*h2_pimi_theta_phi_p_fidcut[6],*h2_pipl_theta_phi_p_beffidcut[6],*h2_pipl_theta_phi_p_fidcut[6],*h2_prot_theta_p[6],*h2_prot_theta_p_cut[6],*h2_pimi_theta_p[6],*h2_pimi_theta_p_cut[6],*h2_pipl_theta_p[6],*h2_pipl_theta_p_cut[6];


  gRandom = new TRandom3();
  //  gRandom->SetSeed(0);

  TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
  TLorentzVector V4_target(0,0,0,target_mass[ftarget]);

  TFile *file_in1=new TFile(Form("FiducialsCorrections/protdeltat_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in=new TFile(Form("FiducialsCorrections/el_Epratio_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in2=new TFile(Form("FiducialsCorrections/pimideltat_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in3= new TFile(Form("FiducialsCorrections/vz_%s_%s.root",ftarget.c_str(),fbeam_en.c_str()));;
  TFile *file_in4=new TFile(Form("FiducialsCorrections/pipldeltat_mom_%s.root",fbeam_en.c_str()));


  TFile *file_in5;
  TFile *file_in6;
  TFile *file_in7;
  std::cout << "fbeam_en " << fbeam_en.c_str() << std::endl;
  if (en_beam[fbeam_en] < 2.300 && en_beam[fbeam_en] > 2.100 ) {
    file_in5 = new TFile("FiducialsCorrections/vz_56Fe_2261_badruns.root");//vertex correction for 56Fe runs with exploded liquid target cell
    file_in6 = new TFile("FiducialsCorrections/vz_3He_2261_2ndrungroup.root");//vertx correction for 3He 2nd group runs
    file_in7 = new TFile("FiducialsCorrections/vz_4He_2261_2ndrungroup.root");//vertex correction for 4He 2nd group runs
  }


  TFile *file_out = new TFile(Form("e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");
  // TFile *file_out = new TFile("e2a_ep_neutrino6_united4_radphot.root", "Recreate");//to submit jobs
  // TFile *file_out = new TFile(Form("/volatile/clas/clase2/Mariana/Skim_Mariana/Files_forTaofeng/e2a_ep_%s_%s_eNp.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");


  //  TTree *skim_tree=(TTree*)fChain->CloneTree(0);


  double pars[3];
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



  vz_corr_func = (TF1 *)file_in3->Get("f_vz");
  el_Epratio_mean=(TF1*)file_in->Get("f_mean");
  el_Epratio_sig=(TF1*)file_in->Get("f_sig");
  prot_deltat_sig=(TF1*)file_in1->Get("sig_pol9");
  prot_deltat_mean=(TF1*)file_in1->Get("mean_pol9");
  pipl_deltat_sig=(TF1*)file_in4->Get("sig_pol9");
  pipl_deltat_mean=(TF1*)file_in4->Get("mean_pol9");
  pimi_deltat_sig=(TF1*)file_in2->Get("sig_pol9");
  pimi_deltat_mean=(TF1*)file_in2->Get("mean_pol9");

  fsum_pimi=new TF1("fsum_pimi",FSum_pimi,0.,5.,2);
  fsub_pimi=new TF1("fsub_pimi",FSub_pimi,0.,5.,2);
  fsum_pipl=new TF1("fsum_pipl",FSum_pipl,0.,5.,2);
  fsub_pipl=new TF1("fsub_pipl",FSub_pipl,0.,5.,2);
  fsum_prot=new TF1("fsum_prot",FSum_prot,0.,5.,2);
  fsub_prot=new TF1("fsub_pprot",FSub_prot,0.,5.,2);
  fsum_e=new TF1("fsum_e",FSum_e,0.,5.,2);
  fsub_e=new TF1("fsub_e",FSub_e,0.,5.,2);


  //Defining EC limits
  up_lim1_ec=new TF1("up_lim1_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim2_ec=new TF1("up_lim2_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim3_ec=new TF1("up_lim3_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim4_ec=new TF1("up_lim4_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim5_ec=new TF1("up_lim5_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim6_ec=new TF1("up_lim6_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim1_ec=new TF1("low_lim1_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim2_ec=new TF1("low_lim2_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim3_ec=new TF1("low_lim3_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim4_ec=new TF1("low_lim4_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim5_ec=new TF1("low_lim5_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim6_ec=new TF1("low_lim6_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  rightside_lim1_ec=new TF1("rightside_lim1_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim1_ec=new TF1("leftside_lim1_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim2_ec=new TF1("rightside_lim2_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim2_ec=new TF1("leftside_lim2_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim3_ec=new TF1("rightside_lim3_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim3_ec=new TF1("leftside_lim3_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim4_ec=new TF1("rightside_lim4_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim4_ec=new TF1("leftside_lim4_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim5_ec=new TF1("rightside_lim5_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim5_ec=new TF1("leftside_lim5_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim6_ec=new TF1("rightside_lim6_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim6_ec=new TF1("leftside_lim6_ec","[0]*(x+[1])+[2]",0,360);


  up_lim1_ec->SetParameters(0.995,30,-0.0001);
  up_lim2_ec->SetParameters(0.995,90,-0.0001);
  up_lim3_ec->SetParameters(0.995,150,-0.0001);
  up_lim4_ec->SetParameters(0.995,210,-0.0001);
  up_lim5_ec->SetParameters(0.995,270,-0.0001);
  up_lim6_ec->SetParameters(0.995,330,-0.0001);
  low_lim1_ec->SetParameters(0.7,30,-0.00005);
  low_lim2_ec->SetParameters(0.7,90,-0.00005);
  low_lim3_ec->SetParameters(0.7,150,-0.00005);
  low_lim4_ec->SetParameters(0.7,210,-0.00005);
  low_lim5_ec->SetParameters(0.7,270,-0.00005);
  low_lim6_ec->SetParameters(0.7,330,-0.00005);
  leftside_lim1_ec->SetParameters(0.11,0,0.03);
  rightside_lim1_ec->SetParameters(-0.11,-60,0.03);
  leftside_lim2_ec->SetParameters(0.11,-60,0.03);
  rightside_lim2_ec->SetParameters(-0.11,-120,0.03);
  leftside_lim3_ec->SetParameters(0.11,-120,0.03);
  rightside_lim3_ec->SetParameters(-0.11,-180,0.03);
  leftside_lim4_ec->SetParameters(0.11,-180,0.03);
  rightside_lim4_ec->SetParameters(-0.11,-240,0.03);
  leftside_lim5_ec->SetParameters(0.11,-240,0.03);
  rightside_lim5_ec->SetParameters(-0.11,-300,0.03);
  leftside_lim6_ec->SetParameters(0.11,-300,0.03);
  rightside_lim6_ec->SetParameters(-0.11,-360,0.03);

 TH1F *h1_ephotphidiff=new TH1F("h1_ephotphidiff","",720,-180,180);
 TH1F *h1_ephotphidiff_cut=new TH1F("h1_ephotphidiff_cut","",720,-180,180);
 TH1F *h1_el_vertuncorr=new TH1F("h1_el_vertuncorr","",200,-10,10);
  TH1F *h1_el_vertcorr=new TH1F("h1_el_vertcorr","",200,-10,10);
  TH1F *h1_el_vertcorr_cut=new TH1F("h1_el_vertcorr_cut","",200,-10,10);
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
 TH1F *h1_xberk = new TH1F("h1_xberk","",400,0,3);
 TH1F *h1_xberk_weight = new TH1F("h1_xberk_weight","",400,0,3);
  TH1F *h1_Q2 = new TH1F("h1_Q2","",400,0,6);
 TH1F *h1_Q2_weight = new TH1F("h1_Q2_weight","",400,0,6);
  TH1F *h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
  TH1F *h1_necpath=new TH1F("h1_necpath","",400,480,600);
TH1F *h1_p_neutron=new TH1F("h1_p_neutron","",400,0,4);
TH1F *h1_p_neutron_uncorr=new TH1F("h1_p_neutron_uncorr","",400,0,4);
  TH1F *h1_Nprot=new TH1F("h1_Nprot","",10,0,5);
  TH1F *h1_Nphot=new TH1F("h1_Nphot","",10,0,5);
  TH1F *h1_Npiphot=new TH1F("h1_Npiphot","",10,0,5);
  TH1F *h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,0,5);
  TH1F *h1_Nphot_EC=new TH1F("h1_Nphot_EC","",10,0,5);
  TH1F *h1_Nphot_LEC=new TH1F("h1_Nphot_LEC","",10,0,5);
  TH1F *h1_photon_E=new TH1F("h1_photon_E","",200,0,4.5);
 TH1F *h1_photon_EC_E=new TH1F("h1_photon_EC_E","",200,0,4.5);
 TH1F *h1_photon_energy=new TH1F("h1_photon_energy","",200,0,4.5);
 TH1F *h1_photon_LAC_E=new TH1F("h1_photon_LAC_E","",200,0,4.5);
TH1F *h1_phot_e_angle= new TH1F("h1_phot_e_angle","",300,0,180);
 TH1F *h1_time_ec=new TH1F("h1_time_ec","",200,-20,20);
  TH1F *h1_Npi=new TH1F("h1_Npi","",10,0,5);
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
 TH1F *neg_m=new TH1F("neg_m","",300,0,2);
 TH1F *pos_m=new TH1F("pos_m","",300,0,2);
TH1F *pos_m_ep=new TH1F("pos_m_ep","",300,0,2);
TH1F *neg_m_ep=new TH1F("neg_m_ep","",300,0,2);

 int N_qe;
  double *x_qe;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
 N_qe=109;
 x_qe=new double[N_qe+1];
   for (int i=0;i<=64;i++)x_qe[i]=-1+i*0.015;
   for (int i=0;i<=44;i++)x_qe[i+65]=-0.04+(i+1)*0.01;
    }

 if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
   N_qe=109;
   x_qe=new double[N_qe+1];
   for (int i=0;i<=64;i++)x_qe[i]=-1+i*0.015;
   for (int i=0;i<=44;i++)x_qe[i+65]=-0.04+(i+1)*0.01;
 }


if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){

   N_qe=82;
   x_qe=new double[N_qe+1];
   for (int i=0;i<=29;i++)x_qe[i]=-1+i*0.03;
   for (int i=0;i<=52;i++)x_qe[i+30]=-0.13+(i+1)*0.01;
 }


TH1F *h1_E_rec_1pi_weight_frac_feed=new TH1F("h1_E_rec_1pi_weight_frac_feed","",N_qe,x_qe);
 TH1F *h1_E_rec_2pi_weight_frac_feed=new TH1F("h1_E_rec_2pi_weight_frac_feed","",N_qe,x_qe);
TH1F *h1_E_rec_3pi_weight_frac_feed=new TH1F("h1_E_rec_3pi_weight_frac_feed","",N_qe,x_qe);
TH1F *h1_E_rec_4pi_weight_frac_feed=new TH1F("h1_E_rec_4pi_weight_frac_feed","",N_qe,x_qe);
TH1F *h1_E_rec_0pi_frac_feed=new TH1F("h1_E_rec_0pi_frac_feed","",N_qe,x_qe);
 TH1F *h1_E_tot_cut2_fracfeed = new TH1F("h1_E_tot_cut2_fracfeed","",N_qe,x_qe);
TH1F *h1_E_rec_cut2_new_fracfeed = new TH1F("h1_E_rec_cut2_new_fracfeed","",N_qe,x_qe);
TH1F *Etot_p_bkgd_fracfeed = new TH1F("Etot_p_bkgd_fracfeed","",N_qe,x_qe);
TH1F *Erec_p_bkgd_fracfeed = new TH1F("Erec_p_bkgd_fracfeed","",N_qe,x_qe);
TH1F *Etot_2p1pi_2p0pi_fracfeed = new TH1F("Etot_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_2p1pi_2p0pi_fracfeed = new TH1F("Erec_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
TH1F *Etot_2p1pi_1p1pi_fracfeed = new TH1F("Etot_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_2p1pi_1p1pi_fracfeed = new TH1F("Erec_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
TH1F *Etot_2p1pi_1p0pi_fracfeed = new TH1F("Etot_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_2p1pi_1p0pi_fracfeed = new TH1F("Erec_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
TH1F *Etot_3pto2p_fracfeed = new TH1F("Etot_3pto2p_fracfeed","",N_qe,x_qe);
TH1F *Erec_3pto2p_fracfeed = new TH1F("Erec_3pto2p_fracfeed","",N_qe,x_qe);
TH1F *Etot_3pto1p_fracfeed = new TH1F("Etot_3pto1p_fracfeed","",N_qe,x_qe);
TH1F *Erec_3pto1p_fracfeed = new TH1F("Erec_3pto1p_fracfeed","",N_qe,x_qe);
TH1F *Etot_4pto3p_fracfeed = new TH1F("Etot_4pto3p_fracfeed","",N_qe,x_qe);
TH1F *Erec_4pto3p_fracfeed = new TH1F("Erec_4pto3p_fracfeed","",N_qe,x_qe);
TH1F *Etot_43pto1p_fracfeed = new TH1F("Etot_43pto1p_fracfeed","",N_qe,x_qe);
TH1F *Erec_43pto1p_fracfeed = new TH1F("Erec_43pto1p_fracfeed","",N_qe,x_qe);
TH1F *Etot_4pto2p_fracfeed = new TH1F("Etot_4pto2p_fracfeed","",N_qe,x_qe);
TH1F *Erec_4pto2p_fracfeed = new TH1F("Erec_4pto2p_fracfeed","",N_qe,x_qe);
TH1F *Etot_4pto1p_fracfeed = new TH1F("Etot_4pto1p_fracfeed","",N_qe,x_qe);
TH1F *Erec_4pto1p_fracfeed = new TH1F("Erec_4pto1p_fracfeed","",N_qe,x_qe);
TH1F *Etot_1p2pi_fracfeed = new TH1F("Etot_1p2pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_1p2pi_fracfeed = new TH1F("Erec_1p2pi_fracfeed","",N_qe,x_qe);
TH1F *Etot_1p3pi_fracfeed = new TH1F("Etot_1p3pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_1p3pi_fracfeed = new TH1F("Erec_1p3pi_fracfeed","",N_qe,x_qe);
TH1F *Etot_2p2pi_fracfeed = new TH1F("Etot_2p2pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_2p2pi_fracfeed = new TH1F("Erec_2p2pi_fracfeed","",N_qe,x_qe);
TH1F *Etot_3p1pi_fracfeed = new TH1F("Etot_3p1pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_3p1pi_fracfeed = new TH1F("Erec_3p1pi_fracfeed","",N_qe,x_qe);
TH1F *Etot_1p2pi_1p0pi_fracfeed = new TH1F("Etot_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
TH1F *Erec_1p2pi_1p0pi_fracfeed = new TH1F("Erec_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
TH1F *h1_E_rec_undetfactor_fracfeed = new TH1F("h1_E_rec_undetfactor_fracfeed","",N_qe,x_qe);
TH1F *h1_E_tot_undetfactor_fracfeed = new TH1F("h1_E_tot_undetfactor_fracfeed","",N_qe,x_qe);
 TH1F *h1_beta_ec = new TH1F("h1_beta_ec","",300,0,2);
 TH1F *h1_beta_lec = new TH1F("h1_beta_lec","",300,0,2);
TH1F *h1_beta_ec_corr = new TH1F("h1_beta_ec_corr","",300,0,2);
TH1F *h1_beta_ec_corr_cut = new TH1F("h1_beta_ec_corr_cut","",300,0,2);
 TH1F *h1_beta_lec_cut = new TH1F("h1_beta_lec_cut","",300,0,2);
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
TH2F *h2_el_E_p_ratio_LAC = new TH2F("h2_el_E_p_ratio_LAC","",200,0,4.5,200,0,0.5);
 TH2D *h2_e_Ein_Eout=new TH2D("h2_e_Ein_Eout","",400,0,2,400,0,1.5);
 TH2D *h2_e_Einout_Etot=new TH2D("h2_e_Einout_Etot","",400,0,4,400,0,4);
 TH2F *h2_e_ec_xy = new TH2F("h2_e_ec_xy","",100,-600,600,100,-600,600);
  TH2F *h2_e_ec_xy_fidcut = new TH2F("h2_e_ec_xy_fidcut","",100,-600,600,100,-600,600);
 TH2F *h2_el_phi_vert = new TH2F("h2_el_phi_vert","",500,-6,7,120,0,380);
  TH2F *h2_el_phi_vert_uncorr = new TH2F("h2_el_phi_vert_uncorr","",500,-6,7,120,0,380);
  TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,0,180);
  TH2F *h2_neutral_theta_phi_LAC = new TH2F("h2_neutral_theta_phi_LAC","",200,0,360,200,0,180);
 TH2F *h2_neutral_theta_phi_EC = new TH2F("h2_neutral_theta_phi_EC","",200,0,360,200,0,180);
 TH2F *h2_neutral_costheta_phi_EC_all = new TH2F("h2_neutral_costheta_phi_EC_all","",200,0,360,200,0,1.1);
 TH2F *h2_neutral_costheta_phi_EC = new TH2F("h2_neutral_costheta_phi_EC","",200,0,360,200,0,1.1);
 TH2F *h2_neutral_theta_phi_LAC_all = new TH2F("h2_neutral_theta_phi_LAC_all","",200,0,360,200,0,180);
 TH2F *h2_neutral_theta_phi_EC_all = new TH2F("h2_neutral_theta_phi_EC_all","",200,0,360,200,0,180);
 TH2F *h2_neutral_theta_phi_LAC_all_fidcut = new TH2F("h2_neutral_theta_phi_LAC_all_fidcut","",200,0,360,200,0,180);
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
 TH2F *h2_Q2_xberk_weight = new TH2F("h2_Q2_xberk_weight","",200,0,3,200,0,5);
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

 TH2F *pipl_delt_p= new TH2F("pipl_delt_p","",300,0.,2.5,300,-15,15);
 TH2F *pimi_delt_p= new TH2F("pimi_delt_p","",300,0.,2.5,300,-15,15);
 TH2F *pos_delt_p= new TH2F("pos_delt_p","",300,0.,2.5,300,-15,15);
 TH2F *neg_delt_p= new TH2F("neg_delt_p","",300,0.,2.5,300,-15,15);
 TH2F *pipl_betta_p = new TH2F("pipl_betta_p","",200,0,4,200,0,1.2);
 TH2F *pimi_betta_p = new TH2F("pimi_betta_p","",200,0,4,200,0,1.2);
 TH2F *neg_betta_p = new TH2F("neg_betta_p","",200,0,4,200,0,1.2);
 TH2F *pos_betta_p = new TH2F("pos_betta_p","",200,0,4,200,0,1.2);
 TH2F *prot_betta_p = new TH2F("prot_betta_p","",200,0,4,200,0,1.2);
TH2F *neg_E_p = new TH2F("neg_E_p","",200,0,3,200,0,200);
TH2F *pos_E_p = new TH2F("pos_E_p","",200,0,3,200,0,200);
TH2F *pimi_E_p = new TH2F("pimi_E_p","",200,0,3,200,0,200);
TH2F *pipl_E_p = new TH2F("pipl_E_p","",200,0,3,200,0,200);
TH2F *pipl_delt_p_ep= new TH2F("pipl_delt_p_ep","",300,0.,2.5,300,-15,15);
 TH2F *pimi_delt_p_ep= new TH2F("pimi_delt_p_ep","",300,0.,2.5,300,-15,15);
TH2F *pipl_delt_p_epfidcut= new TH2F("pipl_delt_p_epfidcut","",300,0.,2.5,300,-15,15);
 TH2F *pimi_delt_p_epfidcut= new TH2F("pimi_delt_p_epfidcut","",300,0.,2.5,300,-15,15);
TH2F *prot_E_p = new TH2F("prot_E_p","",200,0,3,200,0,200);
 TH2F *prot_Deltat_p= new TH2F("prot_Deltat_p","",300,0.,2.7,300,-15,15);
 TH2F *h2_el_vertcorr_runN= new TH2F("h2_el_vertcorr_runN","",30,18370,18436,300,-7,7);
TH2F *h2_phot_e_angle_vsphotE= new TH2F("h2_phot_e_angle_vsphotE","",300,0,3,300,0,180);
TH2F *h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);
TH2F *h2_photE_ephotangle_allsect= new TH2F("h2_photE_ephotangle_allsect","",300,0,100,300,0,3);
TH2F *h2_photE_ephotangle_sect_all= new TH2F("h2_photE_ephotangle_sect_all","",300,0,100,300,0,3);
TH2F *h2_ephotangle_ephotphidiff_cut= new TH2F("h2_ephotangle_ephotphidiff_cut","",360,-180,180,300,0,100);
TH2F *h2_ephotangle_ephotphidiff= new TH2F("h2_ephotangle_ephotphidiff","",360,-180,180,300,0,100);
TH2F *h2_Wepp_ephi= new TH2F("h2_Wepp_ephi","",720,0,360,200,0.85,1.05);
TH2F *h2_Wepp_ephi_corr= new TH2F("h2_Wepp_ephi_corr","",720,0,360,200,0.85,1.05);
TH2F *h2_Wepp_ephi_uncorrprot= new TH2F("h2_Wepp_ephi_uncorrprot","",720,0,360,200,0.85,1.05);
TH2F *h2_Wepp_ephi_corr_uncorrprot= new TH2F("h2_Wepp_ephi_corr_uncorrprot","",720,0,360,200,0.85,1.05);


int n_bins;
double *x_values;

 if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
 n_bins=38;
 x_values=new double[n_bins+1];

   for (int i=0;i<=17;i++)x_values[i]=0.4+i*0.04;
   for (int i=0;i<=20;i++)x_values[i+18]=1.08+(i+1)*0.02;
 }


 if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
 n_bins=54;
 x_values=new double[n_bins+1];

   for (int i=0;i<=23;i++)x_values[i]=i*0.09;
   for (int i=0;i<=30;i++)x_values[i+24]=2.07+(i+1)*0.03;
 }


 if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
 n_bins=38;
 x_values=new double[n_bins+1];

   for (int i=0;i<=21;i++) x_values[i]=i*0.2;
   for (int i=0;i<=16;i++)  x_values[i+22]=4.2+(i+1)*0.05;
   }

/*
 if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
 n_bins=98;
 x_values=new double[n_bins+1];

   for (int i=0;i<=17;i++)x_values[i]=0.4+i*0.04;
   for (int i=0;i<=80;i++)x_values[i+18]=1.08+(i+1)*0.005;
   }


 if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
 n_bins=204;
 x_values=new double[n_bins+1];

   for (int i=0;i<=23;i++)x_values[i]=i*0.09;
   for (int i=0;i<=180;i++)x_values[i+24]=2.07+(i+1)*0.005;
 }


 if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
 n_bins=100;
 x_values=new double[n_bins+1];

   for (int i=0;i<=35;i++) x_values[i]=i*0.12;
   for (int i=0;i<=64;i++)  x_values[i+36]=4.2+(i+1)*0.0125;
   }
*/
 TH1F *Erec_2p_det=new TH1F("Erec_2p_det","",n_bins,x_values);
 TH1F *Etot_2p_det=new TH1F("Etot_2p_det","",n_bins,x_values);
 TH1F *Etot_p_bkgd=new TH1F("Etot_p_bkgd","",n_bins,x_values);
 TH1F *Erec_p_bkgd = new TH1F("Erec_p_bkgd","",n_bins,x_values);
 TH1F *Etot_3pto1p=new TH1F("Etot_3pto1p","",n_bins,x_values);
 TH1F *Erec_3pto1p = new TH1F("Erec_3pto1p","",n_bins,x_values);
 TH1F *Etot_43pto1p=new TH1F("Etot_43pto1p","",n_bins,x_values);
 TH1F *Erec_43pto1p = new TH1F("Erec_43pto1p","",n_bins,x_values);
 TH1F *Etot_3pto2p=new TH1F("Etot_3pto2p","",n_bins,x_values);
 TH1F *Erec_3pto2p = new TH1F("Erec_3pto2p","",n_bins,x_values);
TH1F *Etot_4pto1p=new TH1F("Etot_4pto1p","",n_bins,x_values);
 TH1F *Erec_4pto1p = new TH1F("Erec_4pto1p","",n_bins,x_values);
 TH1F *Etot_4pto3p=new TH1F("Etot_4pto3p","",n_bins,x_values);
 TH1F *Erec_4pto3p = new TH1F("Erec_4pto3p","",n_bins,x_values);
 TH1F *Etot_4pto2p=new TH1F("Etot_4pto2p","",n_bins,x_values);
 TH1F *Erec_4pto2p = new TH1F("Erec_4pto2p","",n_bins,x_values);
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
 TH1F *Erec_1prot  = new TH1F("Erec_1prot","",n_bins,x_values);
 TH1F *Etot_1prot  = new TH1F("Etot_1prot","",n_bins,x_values);
  TH1F *h1_E_rec_cutpi1_piplpimi = new TH1F("h1_E_rec_cutpi1_piplpimi","",n_bins,x_values);
  TH1F * h1_E_tot_cutpi1_piplpimi = new TH1F("h1_E_tot_cutpi1_piplpimi","",n_bins,x_values);
  TH1F *h1_Etot = new TH1F("h1_Etot","",n_bins,x_values);
  TH1F *h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
  TH1F *h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
 TH1F *h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
	TH1F *h1_E_rec_undetfactor  = new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
 TH1F *h1_E_tot_undetfactor  = new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
 TH1F *h1_E_tot_undetfactor_pipl  = new TH1F("h1_E_tot_undetfactor_pipl","",n_bins,x_values);
 TH1F *h1_E_tot_undetfactor_pimi  = new TH1F("h1_E_tot_undetfactor_pimi","",n_bins,x_values);
 TH1F *Etot_1p2pi  = new TH1F("Etot_1p2pi","",n_bins,x_values);
 TH1F *Erec_1p2pi  = new TH1F("Erec_1p2pi","",n_bins,x_values);
 TH1F *Etot_1p3pi  = new TH1F("Etot_1p3pi","",n_bins,x_values);
 TH1F *Erec_1p3pi  = new TH1F("Erec_1p3pi","",n_bins,x_values);
TH1F *Etot_2p2pi  = new TH1F("Etot_2p2pi","",n_bins,x_values);
 TH1F *Erec_2p2pi  = new TH1F("Erec_2p2pi","",n_bins,x_values);
 TH1F *Etot_3p1pi  = new TH1F("Etot_3p1pi","",n_bins,x_values);
 TH1F *Erec_3p1pi  = new TH1F("Erec_3p1pi","",n_bins,x_values);
  TH1F *Etot_1p2pi_1p0pi  = new TH1F("Etot_1p2pi_1p0pi","",n_bins,x_values);
 TH1F *Erec_1p2pi_1p0pi  = new TH1F("Erec_1p2pi_1p0pi","",n_bins,x_values);
 TH1F *Etot_2p1pi_2p0pi  = new TH1F("Etot_2p1pi_2p0pi","",n_bins,x_values);
 TH1F *Erec_2p1pi_2p0pi  = new TH1F("Erec_2p1pi_2p0pi","",n_bins,x_values);
 TH1F *Etot_2p1pi_1p1pi  = new TH1F("Etot_2p1pi_1p1pi","",n_bins,x_values);
 TH1F *Erec_2p1pi_1p1pi  = new TH1F("Erec_2p1pi_1p1pi","",n_bins,x_values);
 TH1F *Etot_2p1pi_1p0pi  = new TH1F("Etot_2p1pi_1p0pi","",n_bins,x_values);
 TH1F *Erec_2p1pi_1p0pi  = new TH1F("Erec_2p1pi_1p0pi","",n_bins,x_values);

 for(int k=0;k<4;k++){

   pos_m_slices[k]=new TH1F(Form("pos_m_slices_%d",k+1),"",300,0,2);
   prot_m_slices[k]=new TH1F(Form("prot_m_slices_%d",k+1),"",300,0,2);
   pipl_m_slices[k]=new TH1F(Form("pipl_m_slices_%d",k+1),"",300,0,2);
 }


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
   h1_Etot_p_bkgd_slice[h]=new TH1F(Form("h1_Etot_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3pto1p_slice[h]= new TH1F(Form("h1_Erec_3pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3pto1p_slice[h]=new TH1F(Form("h1_Etot_3pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3pto2p_slice[h]= new TH1F(Form("h1_Erec_3pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3pto2p_slice[h]=new TH1F(Form("h1_Etot_3pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_3p1pi_slice[h]= new TH1F(Form("h1_Erec_3p1pi_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_3p1pi_slice[h]=new TH1F(Form("h1_Etot_3p1pi_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_43pto1p_slice[h]=new TH1F(Form("h1_Etot_43pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_43pto1p_slice[h]= new TH1F(Form("h1_Erec_43pto1p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto3p_slice[h]= new TH1F(Form("h1_Erec_4pto3p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto3p_slice[h]=new TH1F(Form("h1_Etot_4pto3p_slice_%d",h+1),"",n_bins,x_values);
   h1_Erec_4pto2p_slice[h]= new TH1F(Form("h1_Erec_4pto2p_slice_%d",h+1),"",n_bins,x_values);
   h1_Etot_4pto2p_slice[h]=new TH1F(Form("h1_Etot_4pto2p_slice_%d",h+1),"",n_bins,x_values);
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

 TH1F *h1_Etot_Npi0_Ecalcut[2][6], *h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[2][6],*h1_Etot_p_bkgd_slice_Ecalcut[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_3p1pi_1p0pi[2][6];
 TH1F *h_Etot_piplpimi_subtruct_fact_Ecalcut[2][6], *h1_Etot_p_bkgd_slice_sub_Ecalcut[2][6];
 TH1F *h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio3p[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_ratio4p[2][6], *h1_Etot_p_bkgd_slice_Ecalcut41[2][6], *h1_Etot_p_bkgd_slice_Ecalcut421[2][6], *h1_Etot_p_bkgd_slice_Ecalcut431[2][6], *h1_Etot_p_bkgd_slice_Ecalcut4321[2][6], *h1_Etot_p_bkgd_slice_Ecalcut31[2][6],*h1_Etot_p_bkgd_slice_Ecalcut321[2][6], *h1_Etot_p_bkgd_slice_sub_Ecalcut41[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut42[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut431[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_1p3pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut43[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut31[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut32[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_2p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_1p2pi_1p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p1pi_1p1pi[2][6],*h1_Etot_p_bkgd_slice_sub_Ecalcut_2p2pi_1p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[2][6],*h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[2][6];
 TH1F *h1_E_tot_undetfactor09  = new TH1F("h1_E_tot_undetfactor09","",1,0,6);
 TH1F * h1_E_tot_cut2_09 = new TH1F("h1_E_tot_cut2_09","",1,0,6);
 TH1F *Etot_p_bkgd09  = new TH1F("Etot_p_bkgd09","",1,0,6);
 TH1F *Etot_p321_bkgd09  = new TH1F("Etot_p321_bkgd09","",1,0,6);
 TH1F *Etot_p31_bkgd09  = new TH1F("Etot_p31_bkgd09","",1,0,6);
 TH1F *Etot_p4321_bkgd09  = new TH1F("Etot_p4321_bkgd09","",1,0,6);
 TH1F *Etot_p431_bkgd09  = new TH1F("Etot_p431_bkgd09","",1,0,6);
 TH1F *Etot_p421_bkgd09  = new TH1F("Etot_p421_bkgd09","",1,0,6);
 TH1F *Etot_p41_bkgd09  = new TH1F("Etot_p41_bkgd09","",1,0,6);
 TH1F *Etot_bkgd09_2p1pi_2p0pi  = new TH1F("Etot_bkgd09_2p1pi_2p0pi","",1,0,6);
TH1F *Etot_bkgd09_2p1pi_1p1pi  = new TH1F("Etot_bkgd09_2p1pi_1p1pi","",1,0,6);
TH1F *Etot_bkgd09_2p1pi_1p0pi  = new TH1F("Etot_bkgd09_2p1pi_1p0pi","",1,0,6);
TH1F *Etot_bkgd09_1p2pi_1p0pi  = new TH1F("Etot_bkgd09_1p2pi_1p0pi","",1,0,6);
TH1F *Etot_bkgd09_1p2pi_1p1pi  = new TH1F("Etot_bkgd09_1p2pi_1p1pi","",1,0,6);
TH1F *Etot_bkgd09_1p3pi  = new TH1F("Etot_bkgd09_1p3pi","",1,0,6);
TH1F *Etot_bkgd09_2p2pi  = new TH1F("Etot_bkgd09_2p2pi","",1,0,6);
TH1F *Etot_bkgd09_3p1pi  = new TH1F("Etot_bkgd09_3p1pi","",1,0,6);


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


  h1_E_tot_cut2_09->Sumw2();
  h1_E_tot_undetfactor09->Sumw2();
  Etot_p_bkgd09->Sumw2();
  Etot_p321_bkgd09->Sumw2();
  Etot_p31_bkgd09->Sumw2();
  Etot_p4321_bkgd09->Sumw2();
  Etot_p431_bkgd09->Sumw2();
  Etot_p421_bkgd09->Sumw2();
  Etot_p41_bkgd09->Sumw2();
  Etot_bkgd09_2p1pi_2p0pi->Sumw2();
  Etot_bkgd09_2p1pi_1p1pi->Sumw2();
  Etot_bkgd09_2p1pi_1p0pi->Sumw2();
  Etot_bkgd09_1p2pi_1p0pi->Sumw2();
  Etot_bkgd09_1p2pi_1p1pi->Sumw2();
  Etot_bkgd09_1p3pi->Sumw2();
  Etot_bkgd09_2p2pi->Sumw2();
  Etot_bkgd09_3p1pi->Sumw2();

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
Etot_p_bkgd_fracfeed ->Sumw2();
Erec_p_bkgd_fracfeed ->Sumw2();
Etot_2p1pi_2p0pi_fracfeed ->Sumw2();
Erec_2p1pi_2p0pi_fracfeed ->Sumw2();
Etot_2p1pi_1p1pi_fracfeed ->Sumw2();
Erec_2p1pi_1p1pi_fracfeed ->Sumw2();
Etot_2p1pi_1p0pi_fracfeed ->Sumw2();
Erec_2p1pi_1p0pi_fracfeed ->Sumw2();
Etot_1p3pi_fracfeed ->Sumw2();
Erec_1p3pi_fracfeed ->Sumw2();
Etot_2p2pi_fracfeed ->Sumw2();
Erec_2p2pi_fracfeed ->Sumw2();
Etot_3p1pi_fracfeed ->Sumw2();
Erec_3p1pi_fracfeed ->Sumw2();
Etot_3pto2p_fracfeed->Sumw2();
Erec_3pto2p_fracfeed->Sumw2();
Etot_3pto1p_fracfeed->Sumw2();
Erec_3pto1p_fracfeed->Sumw2();
Etot_4pto3p_fracfeed->Sumw2();
Erec_4pto3p_fracfeed->Sumw2();
Etot_43pto1p_fracfeed->Sumw2();
Erec_43pto1p_fracfeed->Sumw2();
Etot_4pto2p_fracfeed->Sumw2();
Erec_4pto2p_fracfeed->Sumw2();
Etot_4pto1p_fracfeed->Sumw2();
Erec_4pto1p_fracfeed->Sumw2();
Etot_1p2pi_fracfeed->Sumw2();
Erec_1p2pi_fracfeed->Sumw2();
Etot_1p2pi_1p0pi_fracfeed->Sumw2();
Erec_1p2pi_1p0pi_fracfeed->Sumw2();
h1_E_rec_undetfactor_fracfeed->Sumw2();
h1_E_tot_undetfactor_fracfeed->Sumw2();
Erec_2p_det->Sumw2();
 Etot_2p_det->Sumw2();
 Etot_p_bkgd->Sumw2();
 Erec_p_bkgd->Sumw2();
 Etot_3pto1p->Sumw2();
 Erec_3pto1p->Sumw2();
 Etot_43pto1p->Sumw2();
 Erec_43pto1p->Sumw2();
 Etot_3pto2p->Sumw2();
 Erec_3pto2p->Sumw2();
Etot_4pto1p->Sumw2();
 Erec_4pto1p->Sumw2();
 Etot_4pto3p->Sumw2();
 Erec_4pto3p->Sumw2();
 Etot_4pto2p->Sumw2();
 Erec_4pto2p->Sumw2();
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
 Erec_1prot->Sumw2();
 Etot_1prot->Sumw2();
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
 Etot_1p2pi->Sumw2();
 Erec_1p2pi->Sumw2();
 Etot_1p3pi->Sumw2();
 Erec_1p3pi->Sumw2();
 Etot_2p2pi->Sumw2();
 Erec_2p2pi->Sumw2();
 Etot_3p1pi->Sumw2();
 Erec_3p1pi->Sumw2();
 Etot_1p2pi_1p0pi->Sumw2();
 Erec_1p2pi_1p0pi->Sumw2();
 Etot_2p1pi_2p0pi->Sumw2();
 Erec_2p1pi_2p0pi->Sumw2();
 Etot_2p1pi_1p1pi->Sumw2();
 Erec_2p1pi_1p1pi->Sumw2();
 Etot_2p1pi_1p0pi->Sumw2();
 Erec_2p1pi_1p0pi->Sumw2();




 Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<200000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if( jentry%200000 == 0 )
    {
	         gDirectory->Write("hist_Files", TObject::kOverwrite);
		       cout<<jentry<<endl;
    }





    if (runnb==18258 || runnb==18259 || (runnb>18382 && runnb<18438) || (runnb>18220 && runnb<18253)) vz_corr_func->SetParameters(pars); //setting appropriate parameters of the vertex correction function for the runs with the same target and bema energy, but different vertex correction

    if(runnb==18258 || runnb==18259)   {   //setting appropriate e- vertex cut range for the runs with the same target and bema energy, but different vertex correction
      vert_max["56Fe"]=6.;
      vert_min["56Fe"]=5.2;//runs with exploaded cell
    }
    if(runnb>18382 && runnb<18438){
      vert_max["3He"]=0.01;
      vert_min["3He"]=-3.31; //runs with thin exit window
    }
    if(runnb>18220 && runnb<18253){
      vert_max["4He"]=0.77;
      vert_min["4He"]=-2.27;  //runs with 5cm liquid target cell
    }



    if((runnb>18283 && runnb<18289) || (runnb>18300 && runnb<18304) || (runnb>18317 && runnb<18329))fTorusCurrent=750;    //setting appropriate torrus magnet current
    else if ((runnb>18293 && runnb<18301) || (runnb>18305 && runnb<18317) || (runnb>18328 && runnb<18336))fTorusCurrent=1500;
    else fTorusCurrent=2250;


    if(jentry == 0){ //was n_evt == 1 before but jentry = n_evnt - 1
          //SetMomCorrParameters(); Functions is missing F.H. 08/01/19
          SetFiducialCutParameters(fbeam_en);
    }


    int n_elec = 0;
    const int ind_em=0;

    TVector3 e_ec_xyz1(ech_x[ec[ind_em]-1],ech_y[ec[ind_em]-1],ech_z[ec[ind_em]-1]);
    TVector3 el_mom1(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em] ,p[ind_em]*cz[ind_em]);
    double sc_time = sc_t[sc[ind_em] - 1];
    double sc_path = sc_r[sc[ind_em] - 1];
    int sc_paddle = sc_pd[sc[ind_em] - 1];
    int sc_sector = sc_sect[sc[ind_em] - 1];
    float el_vert=vz[ind_em];
    double ec_x=ech_x[ec[ind_em]-1];
    double ec_y=ech_y[ec[ind_em]-1];
    double ec_z=ech_z[ec[ind_em]-1];
    double el_theta=TMath::ACos(cz[ind_em])*TMath::RadToDeg();
    double el_phi_mod=TMath::ATan2(cy[ind_em],cx[ind_em])*TMath::RadToDeg()+30;
    if(el_phi_mod<0)el_phi_mod=el_phi_mod+360;
    int el_ec_sector = ec_sect[ec[ind_em] - 1];
    double el_vert_corr=el_vert+vz_corr(el_phi_mod,el_theta);


    ece = TMath::Max( ec_ei[ec[ind_em] - 1] + ec_eo[ec[ind_em] - 1],
		      etot[ec[ind_em] - 1]);
    el_segment=int((cc_segm[cc[ind_em]-1]-int(cc_segm[cc[ind_em]-1]/1000)*1000)/10);
    el_cc_sector=cc_sect[cc[ind_em]-1];
    el_sccc_timediff=sc_t[cc[ind_em]-1]-cc_t[cc[ind_em]-1]-(sc_r[cc[ind_em]-1]-cc_r[cc[ind_em]-1])/(c*ns_to_s);
    el_cc_nphe=nphe[cc[ind_em]-1]/10.;
    double ec_SC_timediff_uncorr=ec_t[ec[ind_em]-1]-sc_t[sc[ind_em]-1]-(ec_r[ec[ind_em]-1]-sc_r[sc[ind_em]-1])/(c*ns_to_s);


    fsum_e->SetParameters(epratio_sig_cutrange, max_mom);
    fsub_e->SetParameters(epratio_sig_cutrange, max_mom);


  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. && ec[ind_em] > 0.5 && sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] <0 &&  ec_ei[ec[ind_em] - 1] >= 0.03 && ece/p[ind_em]>=fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) && p[ind_em]>=min_good_mom  &&  el_sccc_timediff>=sc_cc_delt_cut_sect[el_cc_sector-1] &&   cc_c2[cc[ind_em]-1]<=0.1)
   {
	h1_el_SCpd[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);
	for(int k=1;k<=6;k++){
	  if(abs(p[ind_em]-0.45)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut[k-1]->Fill(el_phi_mod,el_theta);
	  if(abs(p[ind_em]-1)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut2[k-1]->Fill(el_phi_mod,el_theta);
	}
}


  if(en_beam[fbeam_en]<3.&& en_beam[fbeam_en]>2    && ec[ind_em] > 0.5 && cc[ind_em] > 0.5 &&  sc[ind_em] > 0.5 && q[ind_em] <0  &&  ec_ei[ec[ind_em] - 1] >= 0.06 && ece/p[ind_em]>=fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <=fsum_e->Eval(p[ind_em]) && p[ind_em]>=min_good_mom && cc_c2[cc[ind_em]-1]<0.1 && el_sccc_timediff>=sc_cc_delt_cut_sect[el_cc_sector-1]  && TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass)<=en_beam[fbeam_en]) //electron pid cuts
      {
	h1_el_SCpd[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);
	for(int k=1;k<=6;k++){
	  if(abs(p[ind_em]-1.)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut[k-1]->Fill(el_phi_mod,el_theta);
	  if(abs(p[ind_em]-1.65)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut2[k-1]->Fill(el_phi_mod,el_theta);
	}
}


  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5. && ec[ind_em] > 0.5 &&  sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] <0 && ec_ei[ec[ind_em] - 1] >= 0.055 && ece>=0.33  && ece/p[ind_em]>=fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) && p[ind_em]>=min_good_mom  && cc_c2[cc[ind_em]-1]<0.1 && el_sccc_timediff>=sc_cc_delt_cut_sect[el_cc_sector-1])
      {
	h1_el_SCpd[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);
	for(int k=1;k<=6;k++){
	  if(abs(p[ind_em]-2.5)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut[k-1]->Fill(el_phi_mod,el_theta);
	  if(abs(p[ind_em]-1.4)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_beffidcut2[k-1]->Fill(el_phi_mod,el_theta);
	}
}



    h2_e_ec_xy->Fill(ec_x,ec_y);
    if(!EFiducialCut(fbeam_en, el_mom1) )continue;//theta, phi cuts
    if( !CutUVW(e_ec_xyz1))continue; //u>60, v<360, w<400
    h2_e_ec_xy_fidcut->Fill(ec_x,ec_y);

if(en_beam[fbeam_en]>4. && ec[ind_em] > 0.5 &&  sc[ind_em] > 0.5  &&  q[ind_em] <0 && ec[ind_em] >0.5 &&  sc[ind_em] > 0.5 && q[ind_em] <0 && ec_ei[ec[ind_em] - 1] >= 0.055 && ece>=0.33  && p[ind_em]>=min_good_mom)
    {
      h2_el_E_p_ratio_withoutCC->Fill(p[ind_em], ece/p[ind_em]);
      if(cc[ind_em] > 0.5) h2_el_E_p_ratio_withCC->Fill(p[ind_em], ece/p[ind_em]);
    }


if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. && ec[ind_em] > 0.5 && sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] <0 && ece/p[ind_em]>=fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em]) &&  el_sccc_timediff>=sc_cc_delt_cut_sect[el_cc_sector-1] &&   cc_c2[cc[ind_em]-1]<=0.1)
   {
     // if(ec_ei[ec[ind_em] - 1] >= 0.05)
   h1_el_Etot_cut->Fill(ece);
   //  if (p[ind_em]>=min_good_mom)
 h1_el_Ein_cut->Fill(ec_ei[ec[ind_em] - 1]);
   }

  if(en_beam[fbeam_en]<3.&& en_beam[fbeam_en]>2    && ec[ind_em] > 0.5 && cc[ind_em] > 0.5 &&  sc[ind_em] > 0.5 && q[ind_em] <0  && ece/p[ind_em]>=fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <=fsum_e->Eval(p[ind_em])  && cc_c2[cc[ind_em]-1]<0.1 && el_sccc_timediff>=sc_cc_delt_cut_sect[el_cc_sector-1]  && TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass)<=en_beam[fbeam_en]) //electron pid cuts
    {
      // if(ec_ei[ec[ind_em] - 1] >= 0.06)
  h1_el_Etot_cut->Fill(ece);
  // if (p[ind_em]>=min_good_mom)
      h1_el_Ein_cut->Fill(ec_ei[ec[ind_em] - 1]);
    }


  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5. && ec[ind_em] > 0.5 &&  sc[ind_em] > 0.5 &&  cc[ind_em] > 0.5 && q[ind_em] <0 && ece/p[ind_em]>=fsub_e->Eval(p[ind_em]) && ece/p[ind_em] <= fsum_e->Eval(p[ind_em])  && cc_c2[cc[ind_em]-1]<0.1 && el_sccc_timediff>=sc_cc_delt_cut_sect[el_cc_sector-1])
    {
      // if(ec_ei[ec[ind_em] - 1] >= 0.055)
   h1_el_Etot_cut->Fill(ece);
   // if (p[ind_em]>=min_good_mom)
 h1_el_Ein_cut->Fill(ec_ei[ec[ind_em] - 1]);
    }





 if( ec[ind_em] < 0.5 ||  sc[ind_em] < 0.5 ||  cc[ind_em] < 0.5 || q[ind_em] >=0)
      {continue;}

    h1_el_Etot->Fill(ece);
    h1_el_Ein->Fill(ec_ei[ec[ind_em] - 1]);
    h2_el_E_p_ratio->Fill(p[ind_em], ece/p[ind_em]);
    h2_el_E_p_ratio_LAC->Fill(p[ind_em], lec_etot[lec[ind_em]-1]/p[ind_em]);
    h1_el_cc_chi2->Fill(cc_c2[cc[ind_em]-1]);
    h1_el_cc_deltat[el_cc_sector-1]->Fill(el_sccc_timediff);
    if(el_cc_nphe>2.5) h1_el_cc_deltat_cut[el_cc_sector-1]->Fill(el_sccc_timediff);
   if(el_vert_corr<vert_max[ftarget] && el_vert_corr>vert_min[ftarget])  h1_el_cc_nphe->Fill(el_cc_nphe);
   if(el_vert_corr<vert_max[ftarget] && el_vert_corr>vert_min[ftarget] && el_sccc_timediff>sc_cc_delt_cut_sect[el_cc_sector-1] &&  cc_c2[cc[ind_em]-1]<0.1) h1_el_cc_nphe_cut->Fill(el_cc_nphe);

    h2_e_Ein_Eout->Fill(ec_eo[ec[ind_em]-1],ec_ei[ec[ind_em]-1]);
    h2_e_Einout_Etot->Fill(ece,ec_ei[ec[ind_em]-1]+ec_eo[ec[ind_em]-1]);


  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. && (ec[ind_em] < 0.5 || sc[ind_em] < 0.5 ||  cc[ind_em] < 0.5 || q[ind_em] >=0 || ec_ei[ec[ind_em] - 1] < 0.03 || ece/p[ind_em]<fsub_e->Eval(p[ind_em]) || ece/p[ind_em] >fsum_e->Eval(p[ind_em]) || p[ind_em]<min_good_mom || el_sccc_timediff<sc_cc_delt_cut_sect[el_cc_sector-1] ||   cc_c2[cc[ind_em]-1]>0.1))
   {continue;}

  if(en_beam[fbeam_en]<3.&& en_beam[fbeam_en]>2    &&(ec[ind_em] < 0.5 || cc[ind_em] < 0.5 ||  sc[ind_em] < 0.5 || q[ind_em] >=0 || ec_ei[ec[ind_em] - 1] < 0.06 || ece/p[ind_em]<fsub_e->Eval(p[ind_em]) || ece/p[ind_em] >fsum_e->Eval(p[ind_em]) || p[ind_em]<min_good_mom || cc_c2[cc[ind_em]-1]>=0.1 || el_sccc_timediff<sc_cc_delt_cut_sect[el_cc_sector-1]  || TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass)>en_beam[fbeam_en])) //electron pid cuts
      {continue;}

    if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5. &&(ec[ind_em] < 0.5 ||  sc[ind_em] < 0.5 ||  cc[ind_em] < 0.5 || q[ind_em] >=0 || ec_ei[ec[ind_em] - 1] < 0.055|| ece<0.33  || ece/p[ind_em]<fsub_e->Eval(p[ind_em]) || ece/p[ind_em] >fsum_e->Eval(p[ind_em]) || p[ind_em]<min_good_mom  || cc_c2[cc[ind_em]-1]>=0.1 || el_sccc_timediff<sc_cc_delt_cut_sect[el_cc_sector-1]))
      {continue;}


    h1_el_SCpdfidcut[ec_sect[ec[ind_em]-1]-1]->Fill(sc_pd[sc[ind_em]-1]);
    if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
      for(int k=1;k<=6;k++){
	if(abs(p[ind_em]-2.5)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut[k-1]->Fill(el_phi_mod,el_theta);
	if(abs(p[ind_em]-1.4)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut2[k-1]->Fill(el_phi_mod,el_theta);
      }
    }
    else if (en_beam[fbeam_en]<3.&& en_beam[fbeam_en]>2){
      for(int k=1;k<=6;k++){
	if(abs(p[ind_em]-1.)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut[k-1]->Fill(el_phi_mod,el_theta);
	if(abs(p[ind_em]-1.65)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut2[k-1]->Fill(el_phi_mod,el_theta);
      }
    }
    else {
      for(int k=1;k<=6;k++){
	if(abs(p[ind_em]-0.45)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut[k-1]->Fill(el_phi_mod,el_theta);
	if(abs(p[ind_em]-1)<0.05 && sc_sect[sc[ind_em]-1]==k)	h2_el_theta_phi_p_fidcut2[k-1]->Fill(el_phi_mod,el_theta);
      }
    }

    h2_el_E_p_ratio_cut->Fill(p[ind_em], ece/p[ind_em]);
      if(el_vert_corr<vert_max[ftarget] && el_vert_corr>vert_min[ftarget]) h1_el_cc_nphe_cut2->Fill(el_cc_nphe);

    if(sc_paddle==5) {
      h1_el_ec_sc_timediff->Fill(ec_SC_timediff_uncorr);
      h1_el_ec_sc_timediff_corr->Fill(ec_SC_timediff_uncorr-EC_time_offset[make_pair(ftarget,el_ec_sector)]);
    }
    h1_el_ec_sc_timediff_allSCpd->Fill(ec_SC_timediff_uncorr);
     h1_el_ec_sc_timediff_corr_allSCpd->Fill(ec_SC_timediff_uncorr-EC_time_offset[make_pair(ftarget,el_ec_sector)]);
    h1_el_ec_sc_timediff_sect[el_cc_sector-1]->Fill(ec_SC_timediff_uncorr);
     h1_el_ec_sc_timediff_sect_corr[el_cc_sector-1]->Fill(ec_SC_timediff_uncorr-EC_time_offset[make_pair(ftarget,el_ec_sector)]);
    TVector3 v3_el_ec_uvw=FindUVW(e_ec_xyz1);
    h2_el_ec_sc_timediff_ecu[el_cc_sector-1]->Fill(v3_el_ec_uvw.X(),ec_SC_timediff_uncorr);
    h2_el_ec_sc_timediff_ecv[el_cc_sector-1]->Fill(v3_el_ec_uvw.Y(),ec_SC_timediff_uncorr);
    h2_el_ec_sc_timediff_ecw[el_cc_sector-1]->Fill(v3_el_ec_uvw.Z(),ec_SC_timediff_uncorr);
    h2_el_ec_sc_timediff_SCpd[el_cc_sector-1]->Fill(sc_paddle,ec_SC_timediff_uncorr);





    TLorentzVector V4_el_uncorr(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em],p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass));
    TLorentzVector V4_el(elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cx[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cy[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]*elmom_corr_fact[el_ec_sector-1]*elmom_corr_fact[el_ec_sector-1]+e_mass*e_mass));

    h1_el_mom->Fill(V4_el_uncorr.Rho());
    h1_el_mom_corr->Fill(V4_el.Rho());
    h1_el_mom_ratio->Fill(V4_el.Rho()/V4_el_uncorr.Rho());
    h2_el_pcorr_puncorr->Fill(V4_el.Rho(),V4_el.Rho()/V4_el_uncorr.Rho());
    h2_el_mom_diff->Fill(V4_el.Rho(),V4_el.Rho()-V4_el_uncorr.Rho());

    TVector3 V3_el=V4_el.Vect();
    double fine_struc_const=0.007297;
    double Mott_cross_sec=(fine_struc_const*fine_struc_const*(cz[ind_em]+1))/(2*V4_el.E()*V4_el.E()*(1-cz[ind_em])*(1-cz[ind_em]));
    //  double E_rec=(2*(m_prot-eps)*V4_el.E()+m_prot*m_prot-(m_prot-eps)*(m_prot-eps))/(2*(m_prot-eps-V4_el.E()+V4_el.Rho()*cz[ind_em]));  //using the same value of single nucleon separation E for all targets
    // double E_rec_new= (m_prot*eps+m_prot*V4_el.E())/(m_prot-V4_el.E()+V4_el.Rho()*cz[ind_em]);  //using the same value of single nucleon separation E for all targets
   double E_rec=(2*(m_prot-bind_en[ftarget])*V4_el.E()+m_prot*m_prot-(m_prot-bind_en[ftarget])*(m_prot-bind_en[ftarget]))/(2*(m_prot-bind_en[ftarget]-V4_el.E()+V4_el.Rho()*cz[ind_em]));  //using the same value of single nucleon separation E for Ecal and Eqe
     double E_rec_new= (m_prot*bind_en[ftarget]+m_prot*V4_el.E())/(m_prot-V4_el.E()+V4_el.Rho()*cz[ind_em]);  //using the same value of single nucleon separation E Ecal and Eqe
    E_rec=E_rec_new;
    double nu=-(V4_el-V4_beam).E();
    double Q2=-(V4_el-V4_beam).Mag2();
    double x_berk = Q2/(2*m_prot*nu);
    TVector3 v3q=(V4_beam-V4_el).Vect();
    double W_var=TMath::Sqrt((m_prot+nu)*(m_prot+nu)-v3q*v3q);





    h1_xberk->Fill(x_berk);
    h1_xberk_weight->Fill(x_berk,1/Mott_cross_sec);
    h1_Q2->Fill(Q2);
    h1_Q2_weight->Fill(Q2,1/Mott_cross_sec);
    h2_Q2_nu->Fill(nu,Q2);
    h2_Q2_nu_weight->Fill(nu,Q2,1/Mott_cross_sec);
    h2_Q2_xberk_weight->Fill(x_berk,Q2,1/Mott_cross_sec);
    h1_Wvar->Fill(W_var);
    h1_Wvar_weight->Fill(W_var,1/Mott_cross_sec);
    h2_Q2_W->Fill(W_var,Q2);
    h2_xB_W->Fill(W_var,x_berk);
    h2_Q2_W_weight->Fill(W_var,Q2,1/Mott_cross_sec);

    h1_el_Mott_crosssec->Fill(Mott_cross_sec);
    h2_el_theta_phi->Fill(el_phi_mod,el_theta);
    h1_el_theta->Fill(el_theta);
    h2_el_phi_vert_uncorr->Fill(el_vert,el_phi_mod);
    h2_el_phi_vert->Fill(el_vert_corr,el_phi_mod);
    h1_el_vertuncorr->Fill(el_vert);
    h1_el_vertcorr->Fill(el_vert_corr);
    h2_el_vertcorr_runN->Fill(runnb,el_vert_corr);




    int index_p[20],ind_p,index_pi[20],ind_pi_phot[20],ind_pi;
    int num_p = 0,num_pi=0,num_pi_phot=0,num_pimi=0,num_pipl=0,num_pi_phot_nonrad=0;
    int index_n[20],ec_index_n[20],ec_index_neutrons[20],index_pipl[20],index_pimi[20];
    bool ec_radstat_n[20]={false}, econly_radstat[20]={false};
    int lac_num_n = 0,ec_num_n = 0,ec_neutrons=0;
    double pi_phi,pi_phi_mod, pi_theta,pimi_phi,pimi_phi_mod,pimi_theta,pipl_phi,pipl_phi_mod, pipl_theta,prot_phi,pipl_vert_corr,pipl_vert,pimi_vert_corr,pimi_vert;
    double nectime,necpath,necin,necout,nectot,nbetta_corr,ecpath_corr,n_mom_corr,n_mom,neutr_phi_mod;
    double nech_x, nech_y, nech_z;
    const double n_track_corr_in=3.72607,n_track_corr_out=-22.7486,n_track_corr_both=-0.238097,Mn=0.939565378,ns_to_s=1.0E-9;
    TLorentzVector V4_neutron;
    bool LEC_status=false;
    const double pimi_vertcut=2.5,pipl_vertcut=2.5,phot_rad_cut=40,phot_e_phidiffcut=30;
    double neut_ecx,neut_ecy,neut_ecz, neut_xvert,neut_yvert,neut_zvert, neut_ecpath_corr,neut_ectime_corr,neut_beta_corr;
    TLorentzVector V4_uncorrprot;
    double p_vert_corr,photon_ece,EC_sampling_frac=0.31,ec_deltt=0;
    double pos_m_slices_min[]={0,0.5,1,1.5}, pos_m_slices_max[]={0.5,1,1.5,8};
    TVector3 V3_phot_angles,V3_el_angles_SC,phot_ec_xyz,v3_phot_ec_uvw;

    for( int i = 1; i < TMath::Min(gpart, 20); i++ )
      {
	if( sc[i] > 0 && stat[i] > 0 &&  id[i] == 2212 )
	  {
	    ind_p=i;

	    bett=p[ind_p]/TMath::Sqrt(p[ind_p]*p[ind_p]+m_prot*m_prot);
	    deltt=sc_t[sc[ind_p]-1]-sc_r[sc[ind_p]-1]/(bett*c*ns_to_s) - tr_time;
	    prot_phi=TMath::ATan2(cy[ind_p],cx[ind_p])*TMath::RadToDeg()+30;
	    if(prot_phi<0)prot_phi=prot_phi+360;

	    fsum_prot->SetParameters(prot_delt_cutrange,prot_mom_lim);
	    fsub_prot->SetParameters(prot_delt_cutrange,prot_mom_lim);

	    prot_Deltat_p->Fill(p[ind_p],deltt);


	     if(deltt<fsum_prot->Eval(p[ind_p]) && deltt>fsub_prot->Eval(p[ind_p]) && p[ind_p]>=prot_accept_mom_lim){ //proton pid cut and momentum >0.3GeV cut to get rid of low momentum protons that have a high energy loss and we don't know the efficiency precisely for

	       V4_uncorrprot.SetPxPyPzE(p[ind_p]*cx[ind_p],p[ind_p]*cy[ind_p],p[ind_p]*cz[ind_p],TMath::Sqrt(m_prot*m_prot+p[ind_p]*p[ind_p]));
	       p_vert_corr=vz[ind_p]+vz_corr(prot_phi,TMath::ACos(cz[ind_p])*TMath::RadToDeg());

	    h2_prot_px_py_p->Fill(cx[ind_p],cy[ind_p]);
	    for(int k=1;k<=6;k++){
	      if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. &&  abs(p[ind_p]-0.6)<0.025 && sc_sect[sc[ind_p]-1]==k)h2_prot_theta_phi_p_beffidcut[k-1]->Fill(prot_phi,TMath::ACos(cz[ind_p])*TMath::RadToDeg());
	      if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<5. && abs(p[ind_p]-0.975)<0.025 && sc_sect[sc[ind_p]-1]==k)	h2_prot_theta_phi_p_beffidcut[k-1]->Fill(prot_phi,TMath::ACos(cz[ind_p])*TMath::RadToDeg());
	    }
	    h2_prot_theta_p[sc_sect[sc[ind_p]-1]-1]->Fill(p[i],TMath::ACos(cz[ind_p])*TMath::RadToDeg());

	    if(PFiducialCut(fbeam_en, V4_uncorrprot.Vect())){ //proton fiducial cuts

	      h2_prot_theta_p_cut[sc_sect[sc[ind_p]-1]-1]->Fill(p[i],TMath::ACos(cz[ind_p])*TMath::RadToDeg());
	      for(int k=1;k<=6;k++){
		if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. &&  abs(p[ind_p]-0.6)<0.025 && sc_sect[sc[ind_p]-1]==k) h2_prot_theta_phi_p_fidcut[k-1]->Fill(prot_phi,TMath::ACos(cz[ind_p])*TMath::RadToDeg());
		if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<5. && abs(p[ind_p]-0.975)<0.025 && sc_sect[sc[ind_p]-1]==k)   h2_prot_theta_phi_p_fidcut[k-1]->Fill(prot_phi,TMath::ACos(cz[ind_p])*TMath::RadToDeg());
	      }
	      h1_el_prot_vertdiff_all->Fill(el_vert_corr-p_vert_corr);
	      h2_prot_px_py_p_fidcut->Fill(cx[ind_p],cy[ind_p]);
	      prot_E_p->Fill(p[ind_p],edep[sc[ind_p]-1]);
	      prot_betta_p->Fill(p[ind_p],b[ind_p]);
	      h2_prot_theta_phi->Fill(prot_phi,TMath::ACos(cz[ind_p])*TMath::RadToDeg());
	      for(int k=0;k<4;k++){
		if(p[ind_p]>pos_m_slices_min[k] && p[ind_p]<pos_m_slices_max[k])prot_m_slices[k]->Fill(TMath::Sqrt(p[ind_p]*p[ind_p]/(b[ind_p]*b[ind_p])-p[ind_p]*p[ind_p]));}
	      if((el_vert_corr-p_vert_corr)>vertdiff_min[ftarget] && (el_vert_corr-p_vert_corr)<vertdiff_max[ftarget]){
	    num_p = num_p + 1;
	    index_p[num_p-1]=i;

	      }//vertex cut
	    }
	    }


	  }

	//	if(q[i]<0 &&  sc[i] > 0 && dc[i]>0 && stat[i] > 0 &&  id[i] == -211)
	if(q[i]<0 &&  sc[i] > 0 && dc[i]>0 && stat[i] > 0 )
	  {
	    V3_pimi.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	    bett=p[i]/TMath::Sqrt(p[i]*p[i]+m_pimi*m_pimi);
	    deltt=sc_t[sc[i]-1]-sc_r[sc[i]-1]/(bett*c*ns_to_s) - tr_time;

      	    fsub_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);
	    fsum_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);

	    neg_m->Fill(TMath::Sqrt(p[i]*p[i]/(b[i]*b[i])-p[i]*p[i]));
	    neg_delt_p->Fill(p[i],deltt);
	    neg_E_p->Fill(p[i],edep[sc[i]-1]);
	    neg_betta_p->Fill(p[i],b[i]);

	    if(deltt<fsum_pimi->Eval(p[i]) && deltt>fsub_pimi->Eval(p[i])){

	      pimi_vert=vz[i];
	      pimi_phi=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
	      pimi_phi_mod=pimi_phi+30;
	      if (pimi_phi_mod<0)pimi_phi_mod=pimi_phi_mod+360;
	      pimi_theta=TMath::ACos(cz[i])*TMath::RadToDeg();
	      pimi_vert_corr=pimi_vert+vz_corr(pimi_phi_mod,pimi_theta);


	      if(abs(p[i]-1.)<0.02 && sc_sect[sc[i]-1]==1 && abs(en_beam[fbeam_en]-4.)<1)   h2_pimi_theta_phi_p->Fill(pimi_phi_mod,pimi_theta);
	      if(abs(p[i]-1.)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2.)<1)   h2_pimi_theta_phi_p->Fill(pimi_phi_mod,pimi_theta);
	      if(abs(p[i]-0.5)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1.)<1)   h2_pimi_theta_phi_p->Fill(pimi_phi_mod,pimi_theta);
  for(int k=1;k<=6;k++){
	      if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. &&  abs(p[i]-0.3)<0.025 && sc_sect[sc[i]-1]==k)    h2_pimi_theta_phi_p_beffidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
	      if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<5. && abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k)   h2_pimi_theta_phi_p_beffidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
	    }
	      h2_pimi_theta_phi_beffid->Fill(pimi_phi_mod,pimi_theta);
	      h2_pimi_theta_p[sc_sect[sc[i]-1]-1]->Fill(p[i],TMath::ACos(cz[i])*TMath::RadToDeg());


	      if(PimiFiducialCut(fbeam_en, V3_pimi, &pimi_phimin, &pimi_phimax)){

		h2_pimi_theta_p_cut[sc_sect[sc[i]-1]-1]->Fill(p[i],TMath::ACos(cz[i])*TMath::RadToDeg());
		h1_pimi_prot_vertdiff->Fill(el_vert_corr-pimi_vert_corr);
		if(abs(p[i]-1.)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-4.)<1)	h2_pimi_theta_phi_fidcut->Fill(pimi_phi_mod,pimi_theta);
		if(abs(p[i]-1.)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2.)<1)	h2_pimi_theta_phi_fidcut->Fill(pimi_phi_mod,pimi_theta);
		if(abs(p[i]-0.5)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1.)<1)	h2_pimi_theta_phi_fidcut->Fill(pimi_phi_mod,pimi_theta);
   for(int k=1;k<=6;k++){
		if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. &&  abs(p[i]-0.3)<0.025 && sc_sect[sc[i]-1]==k) h2_pimi_theta_phi_p_fidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
		if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<5. && abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k)   h2_pimi_theta_phi_p_fidcut[k-1]->Fill(pimi_phi_mod,pimi_theta);
	      }
		h2_pimi_theta_phi->Fill(pimi_phi_mod,pimi_theta);

		if(abs(el_vert_corr-pimi_vert_corr)<pimi_vertcut){

		  num_pimi = num_pimi+1;
		  num_pi=num_pi+1;
		  num_pi_phot=num_pi_phot+1;
		  index_pimi[num_pimi - 1]=i;
		  index_pi[num_pi - 1]=i;
		  ind_pi_phot[num_pi_phot - 1]=i;
		  pimi_betta_p->Fill(p[i],b[i]);
		  pimi_E_p->Fill(p[i],edep[sc[i]-1]);
		  pimi_delt_p->Fill(p[i],deltt);
		  num_pi_phot_nonrad=num_pi_phot_nonrad+1;
		}
	      }
	    }
	  }
	//	if(q[i]>0 &&  sc[i] > 0 && dc[i]>0 && stat[i] > 0 &&  id[i] == 211)
			if(q[i]>0 &&  sc[i] > 0 && dc[i]>0 && stat[i] > 0)
	  {
	    V3_pipl.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	    bett=p[i]/TMath::Sqrt(p[i]*p[i]+m_pipl*m_pipl);
	    deltt=sc_t[sc[i]-1]-sc_r[sc[i]-1]/(bett*c*ns_to_s) - tr_time;

	    fsub_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);
	    fsum_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);

	    pos_m->Fill(TMath::Sqrt(p[i]*p[i]/(b[i]*b[i])-p[i]*p[i]));
	    pos_delt_p->Fill(p[i],deltt);
	    pos_betta_p->Fill(p[i],b[i]);
	    pos_E_p->Fill(p[i],edep[sc[i]-1]);
	     for(int k=0;k<4;k++){
	    if(p[i]>pos_m_slices_min[k] && p[i]<pos_m_slices_max[k])pos_m_slices[k]->Fill(TMath::Sqrt(p[i]*p[i]/(b[i]*b[i])-p[i]*p[i]));}

	    if(deltt<fsum_pipl->Eval(p[i]) && deltt>fsub_pipl->Eval(p[i])){

	    pipl_vert=vz[i];
	    pipl_phi=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
	    pipl_phi_mod=pipl_phi+30;
	    if (pipl_phi_mod<0)pipl_phi_mod=pipl_phi_mod+360;
	    pipl_theta=TMath::ACos(cz[i])*TMath::RadToDeg();
	    pipl_vert_corr=pipl_vert+vz_corr(pipl_phi_mod,pipl_theta);


	    if(abs(p[i]-1)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-4)<1)  h2_pipl_theta_phi_p->Fill(pipl_phi_mod,pipl_theta);
	    if(abs(p[i]-1)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2)<1)   h2_pipl_theta_phi_p->Fill(pipl_phi_mod,pipl_theta);
	    if(abs(p[i]-0.5)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1)<1)   h2_pipl_theta_phi_p->Fill(pipl_phi_mod,pipl_theta);
for(int k=1;k<=6;k++){
	      if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. &&  abs(p[i]-0.3)<0.025 && sc_sect[sc[i]-1]==k)h2_pipl_theta_phi_p_beffidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
	      if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<5. && abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k)   h2_pipl_theta_phi_p_beffidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
	    }
	    h2_pipl_theta_phi_beffid->Fill(pipl_phi_mod,pipl_theta);
	    h2_pipl_theta_p[sc_sect[sc[i]-1]-1]->Fill(p[i],TMath::ACos(cz[i])*TMath::RadToDeg());

	      if (PiplFiducialCut(fbeam_en, V3_pipl, &cphil, &cphir)){

		h2_pipl_theta_p_cut[sc_sect[sc[i]-1]-1]->Fill(p[i],TMath::ACos(cz[i])*TMath::RadToDeg());
		h1_pipl_prot_vertdiff->Fill(el_vert_corr-pipl_vert_corr);
	       	if(abs(p[i]-1)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-4)<1)	h2_pipl_theta_phi_fidcut->Fill(pipl_phi_mod,pipl_theta);
		if(abs(p[i]-1)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-2)<1)	h2_pipl_theta_phi_fidcut->Fill(pipl_phi_mod,pipl_theta);
		if(abs(p[i]-0.5)<0.02 && sc_sect[sc[i]-1]==5 && abs(en_beam[fbeam_en]-1)<1)	h2_pipl_theta_phi_fidcut->Fill(pipl_phi_mod,pipl_theta);
   for(int k=1;k<=6;k++){
		if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. &&  abs(p[i]-0.3)<0.025 && sc_sect[sc[i]-1]==k) h2_pipl_theta_phi_p_fidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
		if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<5. && abs(p[i]-0.975)<0.025 && sc_sect[sc[i]-1]==k)   h2_pipl_theta_phi_p_fidcut[k-1]->Fill(pipl_phi_mod,pipl_theta);
	      }

		h2_pipl_theta_phi->Fill(pipl_phi_mod,pipl_theta);

		if (abs(el_vert_corr-pipl_vert_corr)<pipl_vertcut){
		num_pipl = num_pipl+1;
		num_pi=num_pi+1;
		num_pi_phot=num_pi_phot+1;
		index_pipl[num_pipl - 1]=i;
		index_pi[num_pi - 1]=i;
		ind_pi_phot[num_pi_phot - 1]=i;
		pipl_betta_p->Fill(p[i],b[i]);
		pipl_E_p->Fill(p[i],edep[sc[i]-1]);
		pipl_delt_p->Fill(p[i],deltt);
		num_pi_phot_nonrad=num_pi_phot_nonrad+1;
		for(int k=0;k<4;k++){
		  if(p[i]>pos_m_slices_min[k] && p[i]<pos_m_slices_max[k])pipl_m_slices[k]->Fill(TMath::Sqrt(p[i]*p[i]/(b[i]*b[i])-p[i]*p[i]));}
		}	 //vert cut ends
	      }//fidcut ends
	      }//delt cut ends
	  }//pipl ends

			if( lec[i] > 0 && dc[i]<=0  && sc[i]<=0  && stat[i] > 0 &&  q[i] == 0)
	  {
	    neutr_phi_mod=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg()+30;
	    if (neutr_phi_mod<0)neutr_phi_mod=neutr_phi_mod+360;
	    V3_phot_angles.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	    h1_beta_lec->Fill(b[i]);

	    if(b[i] > LEC_photon_beta[ftarget]){  //photon identification

	      h2_neutral_theta_phi_LAC_all->Fill(neutr_phi_mod,TMath::ACos(cz[i])*TMath::RadToDeg());

	      lac_num_n = lac_num_n + 1;
	      index_n[lac_num_n - 1]=i;
	      photon_ece = TMath::Max( ec_ei[ec[i] - 1] + ec_eo[ec[i] - 1],etot[ec[i] - 1]);
	      h2_neutral_theta_phi_LAC_all_fidcut->Fill(neutr_phi_mod,TMath::ACos(cz[i])*TMath::RadToDeg());
	      h1_beta_lec_cut->Fill(b[i]);
	      h1_photon_E->Fill(photon_ece/EC_sampling_frac);
	      h1_photon_LAC_E->Fill(photon_ece/EC_sampling_frac);
	      h2_neutral_theta_phi_LAC->Fill(neutr_phi_mod,TMath::ACos(cz[i])*TMath::RadToDeg());

	    }

	  }
       	if(ec[i] > 0 && dc[i]<=0  && sc[i]<=0  && stat[i] > 0 &&  q[i] == 0)
	  {

	    neut_zvert=vz[i];
	    neut_yvert=vy[i];
	    neut_xvert=vx[i];
	    neut_ecx=ech_x[ec[i]-1];
	    neut_ecy=ech_y[ec[i]-1];
	    neut_ecz=ech_z[ec[i]-1];
	    phot_ec_xyz.SetXYZ(ech_x[ec[i]-1],ech_y[ec[i]-1],ech_z[ec[i]-1]);
	    v3_phot_ec_uvw=FindUVW(phot_ec_xyz);

	    neut_ecpath_corr=TMath::Sqrt((neut_ecx-neut_xvert)*(neut_ecx-neut_xvert)+(neut_ecy-neut_yvert)*(neut_ecy-neut_yvert)+(neut_ecz-neut_zvert)*(neut_ecz-neut_zvert));
	    neut_ectime_corr=neut_ecpath_corr/(b[i]*c*ns_to_s)-EC_time_offset[make_pair(ftarget,ec_sect[ec[i]-1])];
	    neut_beta_corr=neut_ecpath_corr/(neut_ectime_corr*c*ns_to_s);

	    neutr_phi_mod=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg()+30;
	    if (neutr_phi_mod<0)neutr_phi_mod=neutr_phi_mod+360;
	    V3_phot_angles.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
	    V3_el_angles_SC.SetXYZ(dc_cxsc[dc[ind_em]-1],dc_cysc[dc[ind_em]-1],dc_czsc[dc[ind_em]-1]);
	    ec_deltt=ec_t[ec[i]-1]-neut_ecpath_corr/(c*ns_to_s)+EC_time_offset[make_pair(ftarget,ec_sect[ec[i]-1])] - tr_time;
	    //ec_deltt=ec_t[ec[i]-1]-ec_r[ec[i]-1]/(c*ns_to_s) - tr_time;


	    h1_beta_ec->Fill(b[i]);
	    h1_beta_ec_corr->Fill(neut_beta_corr);
	    h1_beta_ec_corr_sect[ec_sect[ec[i]-1]-1]->Fill(neut_beta_corr);

	    if(neut_beta_corr >EC_photon_beta[ftarget]){   //photon identification

	      h2_neutral_theta_phi_EC_all->Fill(neutr_phi_mod,TMath::ACos(cz[i])*TMath::RadToDeg());
	      h2_neutral_costheta_phi_EC_all->Fill(neutr_phi_mod,cz[i]);


	      if(Phot_fid(V3_phot_angles)){//photon fiducial function

	      ec_num_n = ec_num_n + 1;
	      ec_index_n[ec_num_n - 1]=i;
	      num_pi_phot=num_pi_phot+1;
	      ind_pi_phot[num_pi_phot - 1]=i;
	      if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg()<phot_rad_cut && abs(neutr_phi_mod-el_phi_mod)<phot_e_phidiffcut) {
		ec_radstat_n[num_pi_phot - 1]=true;//select radiation photons
		econly_radstat[ec_num_n-1]=true;
	      }
	      if(!ec_radstat_n[num_pi_phot - 1])num_pi_phot_nonrad=num_pi_phot_nonrad+1;
	      photon_ece = TMath::Max( ec_ei[ec[i] - 1] + ec_eo[ec[i] - 1],
		      etot[ec[i] - 1]);
	      h1_time_ec->Fill(ec_deltt);
	      h2_neutral_theta_phi_EC_all_fidcut->Fill(neutr_phi_mod,TMath::ACos(cz[i])*TMath::RadToDeg());
	      h1_beta_ec_corr_cut->Fill(neut_beta_corr);
	      h1_photon_E->Fill(photon_ece/EC_sampling_frac);
	      if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg()<phot_rad_cut && abs(neutr_phi_mod-el_phi_mod)<phot_e_phidiffcut) h1_photon_EC_E->Fill(photon_ece/EC_sampling_frac);
	      h1_phot_e_angle->Fill(V3_phot_angles.Angle(V4_el.Vect())*TMath::RadToDeg());
	      h2_phot_e_angle_vsphotE->Fill(photon_ece/EC_sampling_frac,V3_phot_angles.Angle(V4_el.Vect())*TMath::RadToDeg());
	        if(ec_num_n==1) {
		if(ec_sect[ec[i]-1]!= ec_sect[ec[ind_em]-1]){

	      	  h2_photE_ephotangle_allsect->Fill(V3_phot_angles.Angle(V4_el.Vect())*TMath::RadToDeg(),photon_ece/EC_sampling_frac);
	       	  h2_ephotangle_ephotphidiff_cut->Fill(neutr_phi_mod-el_phi_mod,V3_phot_angles.Angle(V4_el.Vect())*TMath::RadToDeg());
		  h1_ephotphidiff_cut->Fill(neutr_phi_mod-el_phi_mod);
		}
	       	h2_photE_ephotangle_sect_all->Fill(V3_phot_angles.Angle(V4_el.Vect())*TMath::RadToDeg(),photon_ece/EC_sampling_frac);
	       	h2_ephotangle_ephotphidiff->Fill(neutr_phi_mod-el_phi_mod,V3_phot_angles.Angle(V4_el.Vect())*TMath::RadToDeg());
		h1_ephotphidiff->Fill(neutr_phi_mod-el_phi_mod);

	      h2_neutral_theta_phi_EC->Fill(neutr_phi_mod,TMath::ACos(cz[i])*TMath::RadToDeg());
	      h2_neutral_costheta_phi_EC->Fill(neutr_phi_mod,cz[i]);
	      }

	    }//ec phot

	    }//n beta
	  }




      }

    h1_Npi->Fill(num_pi);
    h1_Nprot->Fill(num_p);
    h1_Nphot->Fill(ec_num_n+lac_num_n);
    h1_Nphot_EC->Fill(ec_num_n);
    h1_Nphot_LEC->Fill(lac_num_n);
    h1_Npipl->Fill(num_pipl);
    h1_Npimi->Fill(num_pimi);
    h1_Npiphot->Fill(num_pi_phot);
    h1_Npiphot_norad->Fill(num_pi_phot_nonrad);
    h2_N_prot_pi->Fill(num_pi,num_p);
    h2_N_prot_pi_phot->Fill(num_pi+ec_num_n,num_p);
    h2_N_prot_pi_phot_nonrad->Fill(num_pi_phot_nonrad,num_p);
    h2_N_pi_phot[num_p]->Fill(ec_num_n,num_pi);

    int ind_prot1,ind_prot2,ind_phot;
    TVector3 V3_q,V3_photon;
    double rot_angle,P_N1, P_N2,photon_E;



    if(num_p==1){
      for(int i=1;i<=ec_num_n;i++){

	ind_phot= ec_index_n[i-1];
	V3_photon.SetXYZ(p[ind_phot]*cx[ind_phot],p[ind_phot]*cy[ind_phot],p[ind_phot]*cz[ind_phot]);
	//	if(econly_radstat[i-1]){
	if(V3_photon.Angle(V3_el)*TMath::RadToDeg()<20){

	  photon_E = TMath::Max( ec_ei[ec[ind_phot] - 1] + ec_eo[ec[ind_phot] - 1],
				 etot[ec[ind_phot] - 1]);
	  h1_photon_energy->Fill(photon_E/EC_sampling_frac);
	}

      }
    }






   ec_num_n=0;


    if(el_vert_corr<vert_max[ftarget] && el_vert_corr>vert_min[ftarget])
      {

	h1_el_vertcorr_cut->Fill(el_vert_corr);
	//	if(ec_num_n==0 && lac_num_n==0) skim_tree->Fill();

	// if( num_p>=1 && num_pi>=1){ //for Lucas

	   //	  p[ind_p]=prot_mom_corr;
	//	 skim_tree->Fill();  }


	 //	 if( num_p==1 && ec_num_neutron==1){
	 //   p[ind_p]=prot_mom_corr;
	 //  b[ec_index_neutron[ec_num_neutron - 1]]=ec_neutron_beta[ec_num_neutron - 1];
	 //  ec_t[ec[ec_index_neutron[ec_num_neutron - 1]]-1]=ec_neutron_time[ec_num_neutron - 1];
	  // skim_tree->Fill(); }




	/*	if( num_p>=1){

	  TLorentzVector V4_p[20];
	  double p_vz[20], p_phi[20], p_theta[20],p_vz_corr[20],p_p_corr[20];
	  for(int i=0;i<num_p;i++)
	    {
	      V4_p[i].SetPxPyPzE(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+p[index_p[i]]*p[index_p[i]]));
	      p_vz[i]=vz[index_p[i]];
	      p_phi[i]=TMath::ATan2(cy[index_p[i]],cx[index_p[i]])*TMath::RadToDeg()+30;
	      if(p_phi[i]<0)p_phi[i]=p_phi[i]+360;
	      p_theta[i]=TMath::ACos(cz[index_p[i]])*TMath::RadToDeg();
	      p_vz_corr[i]=p_vz[i]+vz_corr(p_phi[i],p_theta[i]);

	      if(ProtonMomCorrection_He3_4Cell(ftarget,V4_p[i],p_vz_corr[i]) != -1){
		p_p_corr[i]=ProtonMomCorrection_He3_4Cell(ftarget,V4_p[i],p_vz_corr[i]);}
	      else
		{p_p_corr[i]=p[index_p[i]];}

	      p[index_p[i]]=p_p_corr[i];

	    }

	  p[ind_em]=elmom_corr_fact[el_ec_sector-1]*p[ind_em];
	  skim_tree->Fill();
	  }*/




	 if(num_p==2){

   ind_prot1=index_p[0];
   ind_prot2=index_p[1];

   TLorentzVector V4_prot_uncorr1(p[ind_prot1]*cx[ind_prot1],p[ind_prot1]*cy[ind_prot1],p[ind_prot1]*cz[ind_prot1],TMath::Sqrt(m_prot*m_prot+p[ind_prot1]*p[ind_prot1]));
   TLorentzVector V4_prot_uncorr2(p[ind_prot2]*cx[ind_prot2],p[ind_prot2]*cy[ind_prot2],p[ind_prot2]*cz[ind_prot2],TMath::Sqrt(m_prot*m_prot+p[ind_prot2]*p[ind_prot2]));

   float prot_vert1=vz[ind_prot1];
   double p_phi1=TMath::ATan2(cy[ind_prot1],cx[ind_prot1])*TMath::RadToDeg();
   double p_phi_mod1=p_phi1+30;
   if (p_phi_mod1<0)p_phi_mod1=p_phi_mod1+360;
   double p_theta1=TMath::ACos(cz[ind_prot1])*TMath::RadToDeg();
   double prot_vert_corr1=prot_vert1+vz_corr(p_phi_mod1,p_theta1);
   double prot_mom_corr1;
   if(ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr1,prot_vert_corr1) != -1){
     prot_mom_corr1=ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr1,prot_vert_corr1);}
   else
     {prot_mom_corr1=p[ind_prot1];}

   float prot_vert2=vz[ind_prot2];
   double p_phi2=TMath::ATan2(cy[ind_prot2],cx[ind_prot2])*TMath::RadToDeg();
   double p_phi_mod2=p_phi2+30;
   if (p_phi_mod2<0)p_phi_mod2=p_phi_mod2+360;
   double p_theta2=TMath::ACos(cz[ind_prot2])*TMath::RadToDeg();
   double prot_vert_corr2=prot_vert2+vz_corr(p_phi_mod2,p_theta2);
   double prot_mom_corr2;
   if(ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr2,prot_vert_corr2) != -1){
     prot_mom_corr2=ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr2,prot_vert_corr2);}
   else
     {prot_mom_corr2=p[ind_prot2];}

   h1_el_prot_vertdiff1->Fill(el_vert_corr-prot_vert_corr1);
   h1_el_prot_vertdiff2->Fill(el_vert_corr-prot_vert_corr2);

   if((el_vert_corr-prot_vert_corr1)>vertdiff_min[ftarget] && (el_vert_corr-prot_vert_corr1)<vertdiff_max[ftarget] && (el_vert_corr-prot_vert_corr2)>vertdiff_min[ftarget] && (el_vert_corr-prot_vert_corr2)<vertdiff_max[ftarget]){


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

     double mult=0.5*((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2).Mag2()-m_neut*m_neut)/((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2).E()*V4_el_uncorr.Rho()-((V4_target+V4_beam-V4_prot_corr1-V4_prot_corr2).Vect()).Dot(V4_el_uncorr.Vect()));

       h1_e_mom_corrfuct[sc_sect[sc[ind_em]-1]-1]->Fill(mult);

 //---------------------------------- 2p 0pi->  1p0pi   ----------------------------------------------

    double E_tot_2p[2]={0},p_perp_tot_2p[2]={0};
     N_prot_both=0;
     double P_N_2p[2]={0};
   V3_q=(V4_beam-V4_el).Vect();
   prot2_rot_func(fbeam_en, V3_q, V3_2prot_corr,V3_2prot_uncorr, V4_el,E_tot_2p,p_perp_tot_2p,P_N_2p ,&N_prot_both);



 if(ec_num_n==0 && num_pi_phot==0 && N_prot_both!=0){

       for(int f=0;f<num_p;f++){    //looping through two protons

       Etot_p_bkgd->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);
       Erec_p_bkgd->Fill(E_rec_new,P_N_2p[f]*1/Mott_cross_sec);
       Etot_p_bkgd09->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);
       h2_Erec_pperp_2p->Fill(p_perp_tot_2p[f],E_rec_new,P_N_2p[f]*1/Mott_cross_sec);
       Etot_p_bkgd_fracfeed->Fill((E_tot_2p[f]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_N_2p[f]*1/Mott_cross_sec);
       Erec_p_bkgd_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_N_2p[f]*1/Mott_cross_sec);
       h2_pperp_W->Fill(W_var,p_perp_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);
       h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[f]) *TMath::RadToDeg(),-P_N_2p[f]*1/Mott_cross_sec);
       h2_Ecal_Eqe->Fill(E_rec_new,E_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);
       h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);
       h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot_2p[f],-P_N_2p[f]*1/Mott_cross_sec);
       for (int i=0;i<N_pperp;i++){
	 for(int j=0;j<N_Ecal;j++){
	   if(E_tot_2p[f]>Ecal_lowlim[j] && E_tot_2p[f]<Ecal_uplim[j] && p_perp_tot_2p[f]>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut[i][j]->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);

	 }
       }
       for(int i=0;i<n_slice;i++)
	 {
	   if (p_perp_tot_2p[f]<pperp_max[i] && p_perp_tot_2p[f]>pperp_min[i]){
	     h1_Etot_p_bkgd_slice[i]->Fill(E_tot_2p[f],P_N_2p[f]*1/Mott_cross_sec);
	     h1_Erec_p_bkgd_slice[i]->Fill(E_rec_new,P_N_2p[f]*1/Mott_cross_sec);
	   }
	 }

       }//looping through two protons

       Etot_2p_det->Fill(E_tot_2p[0],1/Mott_cross_sec);
       Erec_2p_det->Fill(E_rec_new,1/Mott_cross_sec);
     }//no pions cut


 //---------------------------------- 2p 1pi   ----------------------------------------------
     const int N_2prot=2;
     TVector3 V3_1pi,V3_2p_rotated[N_2prot],V3_1pirot;
     bool pi1_stat=false;
     double N_2p_0pi=0,N_all=0,N_1p_1pi[N_2prot]={0},N_1p_0pi[N_2prot]={0};
 double Ecal_2p1pi_to2p0pi[N_2prot]={0},p_miss_perp_2p1pi_to2p0pi[N_2prot]={0}, P_2p1pi_to2p0pi[N_2prot]={0},N_2p_det=0;
 double   N_pidet=0,N_piundet=0;

 if (ec_num_n==0 && num_pi_phot==1) {

   V3_1pi.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);

 //---------------------------------- 2p 1pi ->2p 0pi   ----------------------------------------------

   double P_2p1pito2p0pi[2]={0}, P_2p1pito1p1pi[2]={0},P_2p1pito1p0pi[2]={0},Ptot=0;
   prot2_pi1_rot_func(fbeam_en, V3_q,V3_2prot_corr,V3_2prot_uncorr,V3_1pi, q[ind_pi_phot[0]],ec_radstat_n[0],V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

 for(int z=0;z<N_2prot;z++){

     Etot_2p1pi_2p0pi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
     Erec_2p1pi_2p0pi->Fill(E_rec_new,P_2p1pito2p0pi[z]*1/Mott_cross_sec);
  h2_Erec_pperp_2p1pi_2p0pi->Fill(p_miss_perp_2p1pi_to2p0pi[z],E_rec_new,P_2p1pito2p0pi[z]*1/Mott_cross_sec);
     Etot_bkgd09_2p1pi_2p0pi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
     Etot_2p1pi_2p0pi_fracfeed->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
     Erec_2p1pi_2p0pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
     h2_pperp_W->Fill(W_var,p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
     h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z]) *TMath::RadToDeg(),P_2p1pito2p0pi[z]*1/Mott_cross_sec);
       h2_Ecal_Eqe->Fill(E_rec_new,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
       h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
       h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);

     for(int i=0;i<n_slice;i++){
       if (p_miss_perp_2p1pi_to2p0pi[z]<pperp_max[i] && p_miss_perp_2p1pi_to2p0pi[z]>pperp_min[i]){
	 h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
	 h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[i]->Fill(E_rec_new,P_2p1pito2p0pi[z]*1/Mott_cross_sec);
       }
     }
     for (int i=0;i<N_pperp;i++){
       for(int j=0;j<N_Ecal;j++){
	 if(Ecal_2p1pi_to2p0pi[z]>Ecal_lowlim[j] && Ecal_2p1pi_to2p0pi[z]<Ecal_uplim[j] && p_miss_perp_2p1pi_to2p0pi[z]>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[i][j]->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*1/Mott_cross_sec);
       }
     }
 //---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

   Etot_2p1pi_1p1pi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
   Erec_2p1pi_1p1pi->Fill(E_rec_new,P_2p1pito1p1pi[z]*1/Mott_cross_sec);
 h2_Erec_pperp_2p1pi_1p1pi->Fill(p_perp_tot_2p[z],E_rec_new,P_2p1pito1p1pi[z]*1/Mott_cross_sec);
   Etot_bkgd09_2p1pi_1p1pi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
   Etot_2p1pi_1p1pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
   Erec_2p1pi_1p1pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
  h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
  h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_2p1pito1p1pi[z]*1/Mott_cross_sec);
   h2_Ecal_Eqe->Fill(E_rec_new,E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
       h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
       h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);

   for(int i=0;i<n_slice;i++){
     if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
       h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
       h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[i]->Fill(E_rec_new,P_2p1pito1p1pi[z]*1/Mott_cross_sec);
     }
   }
   for (int i=0;i<N_pperp;i++){
     for(int j=0;j<N_Ecal;j++){
       if(E_tot_2p[z]>Ecal_lowlim[j] && E_tot_2p[z]<Ecal_uplim[j] && p_perp_tot_2p[z]>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[i][j]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*1/Mott_cross_sec);
     }
   }
//---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------
  Etot_2p1pi_1p0pi->Fill(E_tot_2p[z], P_2p1pito1p0pi[z]*1/Mott_cross_sec);
   Erec_2p1pi_1p0pi->Fill(E_rec_new,P_2p1pito1p0pi[z]*1/Mott_cross_sec);
   h2_Erec_pperp_2p1pi_1p0pi->Fill(p_perp_tot_2p[z],E_rec_new,P_2p1pito1p0pi[z]*1/Mott_cross_sec);
   Etot_bkgd09_2p1pi_1p0pi->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
   Etot_2p1pi_1p0pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
   Erec_2p1pi_1p0pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
   h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
 h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
  h2_Ecal_Eqe->Fill(E_rec_new,E_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
       h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);
       h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot_2p[z],-P_2p1pito1p0pi[z]*1/Mott_cross_sec);

   for(int i=0;i<n_slice;i++)
     {
       if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	 h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
	 h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_rec_new,P_2p1pito1p0pi[z]*1/Mott_cross_sec);
       }
     }
   for (int i=0;i<N_pperp;i++){
     for(int j=0;j<N_Ecal;j++){
       if(E_tot_2p[z]>Ecal_lowlim[j] && E_tot_2p[z]<Ecal_uplim[j] && p_perp_tot_2p[z]>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[i][j]->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*1/Mott_cross_sec);
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

 if (ec_num_n==0 && num_pi_phot==2) {

   V3_2pi[0].SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
   V3_2pi[1].SetXYZ(p[ind_pi_phot[1]]*cx[ind_pi_phot[1]],p[ind_pi_phot[1]]*cy[ind_pi_phot[1]],p[ind_pi_phot[1]]*cz[ind_pi_phot[1]]);
   q_pi2[0]=q[ind_pi_phot[0]];
   q_pi2[1]=q[ind_pi_phot[1]];
   ecstat_pi2[0]=ec_radstat_n[0];
   ecstat_pi2[1]=ec_radstat_n[1];

   prot2_pi2_rot_func(fbeam_en, V3_q,V3_2prot_corr,V3_2prot_uncorr,V3_2pi,q_pi2,ecstat_pi2 ,V4_el, Ecal_2p2pi,p_miss_perp_2p2pi,Ptot_2p);

//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

for(int z=0;z<N_2prot;z++){

  Etot_2p2pi->Fill(E_tot_2p[z], Ptot_2p[z]*1/Mott_cross_sec);
  Erec_2p2pi->Fill(E_rec_new,Ptot_2p[z]*1/Mott_cross_sec);
  h2_Erec_pperp_2p2pi->Fill(p_perp_tot_2p[z],E_rec_new,Ptot_2p[z]*1/Mott_cross_sec);
  Etot_bkgd09_2p2pi->Fill(E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
  Etot_2p2pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],Ptot_2p[z]*1/Mott_cross_sec);
  Erec_2p2pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],Ptot_2p[z]*1/Mott_cross_sec);
  h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
  h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),Ptot_2p[z]*1/Mott_cross_sec);
  h2_Ecal_Eqe->Fill(E_rec_new,E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
  h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
  h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);

   for(int i=0;i<n_slice;i++)
     {
       if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
	 h1_Etot_p_bkgd_slice_2p2pi[i]->Fill(E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
	 h1_Erec_p_bkgd_slice_2p2pi[i]->Fill(E_rec_new,Ptot_2p[z]*1/Mott_cross_sec);
       }
     }
   for (int i=0;i<N_pperp;i++){
     for(int j=0;j<N_Ecal;j++){
       if(E_tot_2p[z]>Ecal_lowlim[j] && E_tot_2p[z]<Ecal_uplim[j] && p_perp_tot_2p[z]>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[i][j]->Fill(E_tot_2p[z],Ptot_2p[z]*1/Mott_cross_sec);
     }
   }

 }   //Filling the histogram for two protons




  }//2pi requirement




   }//2prot vert cut


   }//2prot requirement








	 if(num_p==3){


	   const int N_3p=3;
	   TLorentzVector V4_p_uncorr[N_3p], V4_p_corr[N_3p],V4_prot_el[N_3p];
	   float prot_vz[N_3p];
	   double prot_phi[N_3p],prot_theta[N_3p],prot_vz_corr[N_3p],prot_p_corr[N_3p];
	   TVector3 V3_prot[N_3p],V3_prot_corr[N_3p],V3_3p_rot[N_3p];
	   double E_cal[N_3p],p_miss_perp[N_3p],P_3pto1p[N_3p];
	   double N_p1[N_3p]={0},N_p_three=0;
	   double N_p12[N_3p]={0},N_p13[N_3p]={0},N_p23[N_3p]={0},N_p_two=0, E_cal_3pto1p[3]={0},p_miss_perp_3pto1p[3]={0};
	   N_comb=3;
	   const int N_2p=2;
	   double E_cal_3pto2p[3][N_2p]={0},p_miss_perp_3pto2p[3][N_2p]={0},P_3pto2p[3][N_2p]={0};
	   TVector3 V3_prot_el_3pto2p[N_3p][N_2p];
	   TVector3 V3_2p_rot[N_2p], V3_prot_el[N_3p][N_2p];
 	   bool prot_stat[N_3p]={false};
	   // prot3_rot_func( V3_q, V3_prot_corr,V3_prot,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

	   for(int i=0;i<N_3p;i++)
	     {
	       V4_p_uncorr[i].SetPxPyPzE(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+p[index_p[i]]*p[index_p[i]]));
	       prot_vz[i]=vz[index_p[i]];
	       prot_phi[i]=TMath::ATan2(cy[index_p[i]],cx[index_p[i]])*TMath::RadToDeg()+30;
	       if(prot_phi[i]<0)prot_phi[i]=prot_phi[i]+360;
	       prot_theta[i]=TMath::ACos(cz[index_p[i]])*TMath::RadToDeg();
	       prot_vz_corr[i]=prot_vz[i]+vz_corr(prot_phi[i],prot_theta[i]);

	       if(ProtonMomCorrection_He3_4Cell(ftarget,V4_p_uncorr[i],prot_vz_corr[i]) != -1){
		 prot_p_corr[i]=ProtonMomCorrection_He3_4Cell(ftarget,V4_p_uncorr[i],prot_vz_corr[i]);}
	       else
		 {prot_p_corr[i]=p[index_p[i]];}

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

	   if((el_vert_corr- prot_vz_corr[0])>vertdiff_min[ftarget] && (el_vert_corr- prot_vz_corr[0])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot_vz_corr[1])>vertdiff_min[ftarget] && (el_vert_corr- prot_vz_corr[1])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot_vz_corr[2])>vertdiff_min[ftarget] && (el_vert_corr- prot_vz_corr[2])<vertdiff_max[ftarget]){

	     for(int i=0;i<N_3p;i++){
	       N_p1[i]=0;
	     }
	     V3_q=(V4_beam-V4_el).Vect();


	     prot3_rot_func(fbeam_en, V3_q, V3_prot_corr,V3_prot,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

	     if(ec_num_n==0 && num_pi_phot==0 && N_p_three!=0){
 //-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------
	       for(int count=0;count<N_comb;count++)    {

		 for(int j=0;j<N_2p;j++)    {

		   Etot_3pto2p->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*1/Mott_cross_sec);
		   Erec_3pto2p->Fill(E_rec_new, P_3pto2p[count][j]*1/Mott_cross_sec);
		   Etot_p321_bkgd09->Fill(E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
		   h2_Erec_pperp_321p->Fill(p_miss_perp_3pto2p[count][j],E_rec_new,P_3pto2p[count][j]*1/Mott_cross_sec);
		   Etot_3pto2p_fracfeed->Fill((E_cal_3pto2p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto2p[count][j]*1/Mott_cross_sec);
		   Erec_3pto2p_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_3pto2p[count][j]*1/Mott_cross_sec);
		   h2_pperp_W->Fill(W_var,p_miss_perp_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
		   h1_theta0->Fill((V4_beam.Vect()).Angle(V3_prot_el[count][j])*TMath::RadToDeg(),P_3pto2p[count][j]*1/Mott_cross_sec);
		   h2_Ecal_Eqe->Fill(E_rec_new,E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
		   h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);
		   h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_cal_3pto2p[count][j],P_3pto2p[count][j]*1/Mott_cross_sec);

		   for (int n=0;n<N_pperp;n++){
		     for(int m=0;m<N_Ecal;m++){
		       if(E_cal_3pto2p[count][j]>Ecal_lowlim[m] && E_cal_3pto2p[count][j]<Ecal_uplim[m] && p_miss_perp_3pto2p[count][j]>pperp_cut[n])   h1_Etot_p_bkgd_slice_Ecalcut321[n][m]->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*1/Mott_cross_sec);

		     }
		   }


		   for(int i=0;i<n_slice;i++)
		     {
		       if (p_miss_perp_3pto2p[count][j]<pperp_max[i] && p_miss_perp_3pto2p[count][j]>pperp_min[i]){

			 h1_Etot_3pto2p_slice[i]->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*1/Mott_cross_sec);
			 h1_Erec_3pto2p_slice[i]->Fill(E_rec_new, P_3pto2p[count][j]*1/Mott_cross_sec);
		       }
		     }
		 }
	       }
 //-----------------------------------------  3p to 1p  -----------------------------------------------------------------------
	       for(int j=0;j<N_3p;j++)    {

		 P_3pto1p[j]= N_p1[j]/N_p_three;
		 Etot_3pto1p->Fill(E_cal[j], P_3pto1p[j]*1/Mott_cross_sec);
		 Erec_3pto1p->Fill(E_rec_new,P_3pto1p[j]*1/Mott_cross_sec);
		 Etot_p31_bkgd09->Fill(E_cal[j],P_3pto1p[j]*1/Mott_cross_sec);
		 h2_Erec_pperp_31p->Fill(p_miss_perp[j],E_rec_new,P_3pto1p[j]*1/Mott_cross_sec);
		 Etot_3pto1p_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto1p[j]*1/Mott_cross_sec);
		 Erec_3pto1p_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_3pto1p[j]*1/Mott_cross_sec);
		 h2_pperp_W->Fill(W_var,p_miss_perp[j],-P_3pto1p[j]*1/Mott_cross_sec);
		 h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot[j])*TMath::RadToDeg(),-P_3pto1p[j]*1/Mott_cross_sec);
		 h2_Ecal_Eqe->Fill(E_rec_new,E_cal[j],-P_3pto1p[j]*1/Mott_cross_sec);
		 h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_cal[j],-P_3pto1p[j]*1/Mott_cross_sec);
		 h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_cal[j],-P_3pto1p[j]*1/Mott_cross_sec);

		 for (int n=0;n<N_pperp;n++){
		   for(int m=0;m<N_Ecal;m++){
		     if(E_cal[j]>Ecal_lowlim[m] && E_cal[j]<Ecal_uplim[m] && p_miss_perp[j]>pperp_cut[n])   h1_Etot_p_bkgd_slice_Ecalcut31[n][m]->Fill(E_cal[j], P_3pto1p[j]*1/Mott_cross_sec);

		   }
		 }

		 for(int i=0;i<n_slice;i++)
		   {
		     if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){

		       h1_Etot_3pto1p_slice[i]->Fill(E_cal[j],P_3pto1p[j]*1/Mott_cross_sec);
		       h1_Erec_3pto1p_slice[i]->Fill(E_rec_new,P_3pto1p[j]*1/Mott_cross_sec);
		     }
		   }

	       }


	     }//no pions cut


	     //----------------------------------3p 1pi ----------------------------------------------------------


 if (ec_num_n==0 && num_pi_phot==1) {

   double P_tot_3p[N_3p]={0},Ecal_3p1pi[N_3p]={0},p_miss_perp_3p1pi[N_3p]={0};
 TVector3 V3_pi_phot;
 V3_q=(V4_beam-V4_el).Vect();

   V3_pi_phot.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);

   prot3_pi1_rot_func(fbeam_en,  V3_q,V3_prot_corr,V3_prot, V3_pi_phot,q[ind_pi_phot[0]],ec_radstat_n[0], V4_el,  Ecal_3p1pi,p_miss_perp_3p1pi, P_tot_3p);

	       for(int j=0;j<N_3p;j++)    {

		 Etot_3p1pi->Fill(E_cal[j], P_tot_3p[j]*1/Mott_cross_sec);
		 Erec_3p1pi->Fill(E_rec_new,P_tot_3p[j]*1/Mott_cross_sec);
		 Etot_bkgd09_3p1pi->Fill(E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		 h2_Erec_pperp_3p1pi->Fill(p_miss_perp[j],E_rec_new,P_tot_3p[j]*1/Mott_cross_sec);
		 Etot_3p1pi_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_tot_3p[j]*1/Mott_cross_sec);
		 Erec_3p1pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_3p[j]*1/Mott_cross_sec);
		 h2_pperp_W->Fill(W_var,p_miss_perp[j],P_tot_3p[j]*1/Mott_cross_sec);
		 h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot[j])*TMath::RadToDeg(),P_tot_3p[j]*1/Mott_cross_sec);
		 h2_Ecal_Eqe->Fill(E_rec_new,E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		 h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		 h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);

		 for (int n=0;n<N_pperp;n++){
		   for(int m=0;m<N_Ecal;m++){
		     if(E_cal[j]>Ecal_lowlim[m] && E_cal[j]<Ecal_uplim[m] && p_miss_perp[j]>pperp_cut[n])   h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[n][m]->Fill(E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);

		   }
		 }
		 for(int i=0;i<n_slice;i++)
		   {
		     if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){

		       h1_Etot_3p1pi_slice[i]->Fill(E_cal[j],P_tot_3p[j]*1/Mott_cross_sec);
		       h1_Erec_3p1pi_slice[i]->Fill(E_rec_new,P_tot_3p[j]*1/Mott_cross_sec);
		     }
		   }

	       }

 }// 1 pi requirement ends


	   }//vertex cut
	 }//N_p=3  3proton requiremnet











	 	 if(num_p==4){

	   const int N_p4=4;
	   TLorentzVector V4_p4_uncorr[N_p4], V4_p4_corr[N_p4],V4_prot4_el[N_p4];
	   float prot4_vz[N_p4];
	   double prot4_phi[N_p4],prot4_theta[N_p4],prot4_vz_corr[N_p4],prot4_p_corr[N_p4];
	   TVector3 V3_prot4[N_p4],V3_prot4_corr[N_p4];
	   double E_cal_p4[N_p4]={0},p_miss_perp_p4[N_p4]={0},P_4pto1p[N_p4]={0};

	   TVector3 V3_p4_rot[N_p4];
	   bool prot4_stat[N_p4]={false};
	   const int Ncomb_4to1=4,Ncomb_4to2=6, Ncomb_4to3=4;
	   double N_p4_p1[Ncomb_4to1]={0},N_p4_p2[Ncomb_4to2]={0},N_p4_p3[Ncomb_4to3]={0},N_p_four=0;



	   for(int i=0;i<N_p4;i++)
	     {
	       V4_p4_uncorr[i].SetPxPyPzE(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+p[index_p[i]]*p[index_p[i]]));
	       prot4_vz[i]=vz[index_p[i]];
	       prot4_phi[i]=TMath::ATan2(cy[index_p[i]],cx[index_p[i]])*TMath::RadToDeg()+30;
	       if(prot4_phi[i]<0)prot4_phi[i]=prot4_phi[i]+360;
	       prot4_theta[i]=TMath::ACos(cz[index_p[i]])*TMath::RadToDeg();
	       prot4_vz_corr[i]=prot4_vz[i]+vz_corr(prot4_phi[i],prot4_theta[i]);

	       if(ProtonMomCorrection_He3_4Cell(ftarget,V4_p4_uncorr[i],prot4_vz_corr[i]) != -1){
		 prot4_p_corr[i]=ProtonMomCorrection_He3_4Cell(ftarget,V4_p4_uncorr[i],prot4_vz_corr[i]);}
	       else
		 {prot4_p_corr[i]=p[index_p[i]];}

	       V3_prot4[i].SetXYZ(p[index_p[i]]*cx[index_p[i]],p[index_p[i]]*cy[index_p[i]],p[index_p[i]]*cz[index_p[i]]);
	       V3_prot4_corr[i].SetXYZ(prot4_p_corr[i]*cx[index_p[i]],prot4_p_corr[i]*cy[index_p[i]],prot4_p_corr[i]*cz[index_p[i]]);
	       V4_p4_corr[i].SetPxPyPzE(prot4_p_corr[i]*cx[index_p[i]],prot4_p_corr[i]*cy[index_p[i]],prot4_p_corr[i]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+prot4_p_corr[i]*prot4_p_corr[i]));
	       V4_prot4_el[i]=V4_p4_corr[i]+V4_el;
	       E_cal_p4[i]=V4_el.E()+ V4_p4_corr[i].E()-m_prot+bind_en[ftarget];
	       p_miss_perp_p4[i]=TMath::Sqrt(V4_prot4_el[i].Px()*V4_prot4_el[i].Px()+V4_prot4_el[i].Py()*V4_prot4_el[i].Py());
	     }


	   h1_el_4prot_vertdiff1->Fill(el_vert_corr- prot4_vz_corr[0]);
	   h1_el_4prot_vertdiff2->Fill(el_vert_corr- prot4_vz_corr[1]);
	   h1_el_4prot_vertdiff3->Fill(el_vert_corr- prot4_vz_corr[2]);
	   h1_el_4prot_vertdiff4->Fill(el_vert_corr- prot4_vz_corr[3]);



	   if((el_vert_corr- prot4_vz_corr[0])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[0])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot4_vz_corr[1])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[1])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot4_vz_corr[2])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[2])<vertdiff_max[ftarget] &&
              (el_vert_corr- prot4_vz_corr[3])>vertdiff_min[ftarget] && (el_vert_corr- prot4_vz_corr[3])<vertdiff_max[ftarget]){

	     V3_q=(V4_beam-V4_el).Vect();


	     if (ec_num_n==0 && num_pi_phot==0){

	       for(int g=0; g<N_tot; g++){

		 rot_angle=gRandom->Uniform(0,2*TMath::Pi());

		 for(int i=0;i<N_p4;i++) {
		   V3_p4_rot[i]= V3_prot4[i];
		   V3_p4_rot[i].Rotate(rot_angle,V3_q);
		 }


		 for(int ind_p=0;ind_p<N_p4;ind_p++) prot4_stat[ind_p]=PFiducialCut(fbeam_en, V3_p4_rot[ind_p]);

		 if( prot4_stat[0]  && !prot4_stat[1]   && !prot4_stat[2] && !prot4_stat[3])  N_p4_p1[0]=N_p4_p1[0]+1;//Detecting 1p out of 4p
		 if(!prot4_stat[0]  &&   prot4_stat[1]  && !prot4_stat[2] && !prot4_stat[3])  N_p4_p1[1]=N_p4_p1[1]+1;
		 if(!prot4_stat[0]  &&  !prot4_stat[1]  &&  prot4_stat[2] && !prot4_stat[3])  N_p4_p1[2]=N_p4_p1[2]+1;
		 if(!prot4_stat[0]  &&  !prot4_stat[1]  && !prot4_stat[2] &&  prot4_stat[3])  N_p4_p1[3]=N_p4_p1[3]+1;
		 if( prot4_stat[0]  &&  prot4_stat[1]   &&  prot4_stat[2] &&  prot4_stat[3])   N_p_four=N_p_four+1;   //Detecting 4p out of 4p

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




	       N_comb=3;    //number of 2 proton combination out of three
	       const int N_2p=2, N_3p=3;
	       double E_cal_4pto3p[3][N_2p]={0},p_miss_perp_4pto3p[3][N_2p]={0},P_4pto3p[3][N_2p]={0};
	       TVector3 V3_prot[N_3p],V3_prot_uncorr[N_3p],V3_prot_el_4pto3p[N_3p][N_2p],V3_el_prot[N_comb][N_2p];
	       double N_p_three=0,N_p1[N_3p]={0};
	       double P_4pto1p[N_p4]={0},E_cal_43pto1p[N_3p],p_miss_perp_43pto1p[N_3p], P_43pto1p[3]={0};



	       for(int g=0;g<Ncomb_4to3;g++){   //estimating the undetected 4p contribution to  3p
//SHouldnt that be g == 0,1,2,3 F.H. 08/05/19
		 if(g=0) {
		   V3_prot_uncorr[0]=V3_prot4[0]; V3_prot_uncorr[1]=V3_prot4[1]; V3_prot_uncorr[2]=V3_prot4[2];
		   V3_prot[0]=V3_prot4_corr[0]; V3_prot[1]=V3_prot4_corr[1]; V3_prot[2]=V3_prot4_corr[2];
		 }
		 if(g=1){
		   V3_prot_uncorr[0]=V3_prot4[0]; V3_prot_uncorr[1]=V3_prot4[1]; V3_prot_uncorr[2]=V3_prot4[3];
		   V3_prot[0]=V3_prot4_corr[0]; V3_prot[1]=V3_prot4_corr[1]; V3_prot[2]=V3_prot4_corr[3];
		 }
		 if(g=2){
		   V3_prot_uncorr[0]=V3_prot4[0]; V3_prot_uncorr[1]=V3_prot4[2]; V3_prot_uncorr[2]=V3_prot4[3];
		   V3_prot[0]=V3_prot4_corr[0]; V3_prot[1]=V3_prot4_corr[2]; V3_prot[2]=V3_prot4_corr[3];
		 }
		 if(g=3){
		   V3_prot_uncorr[0]=V3_prot4[1]; V3_prot_uncorr[1]=V3_prot4[2]; V3_prot_uncorr[2]=V3_prot4[3];
		   V3_prot[0]=V3_prot4_corr[1]; V3_prot[1]=V3_prot4_corr[2]; V3_prot[2]=V3_prot4_corr[3];
		 }

		 prot3_rot_func(fbeam_en,  V3_q, V3_prot, V3_prot_uncorr,V4_el,E_cal_4pto3p,p_miss_perp_4pto3p, P_4pto3p,N_p1,E_cal_43pto1p,p_miss_perp_43pto1p,&N_p_three);

		 V3_el_prot[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
		 V3_el_prot[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
		 V3_el_prot[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
		 V3_el_prot[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
		 V3_el_prot[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
		 V3_el_prot[2][1]=V4_el.Vect()+V3_prot_uncorr[2];

		 if( N_p_three!=0 && N_p_four!=0){
	   //-----------------------------------------  4p to 3p->2->1  -----------------------------------------------------------------------
		   	   for(int count=0;count<N_comb;count++)    { //looping through number of 2 proton combination out of 3 protons

		     for(int j=0;j<N_2p;j++)    {  //looping through number of 1 proton combination out of 2 protons

		       Etot_4pto3p->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       Erec_4pto3p->Fill(E_rec_new, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       Etot_p4321_bkgd09->Fill(E_cal_4pto3p[count][j],P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       h2_Erec_pperp_4321p->Fill(p_miss_perp_4pto3p[count][j],E_rec_new,P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       Etot_4pto3p_fracfeed->Fill((E_cal_4pto3p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       Erec_4pto3p_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       h2_pperp_W->Fill(W_var,p_miss_perp_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       h1_theta0->Fill((V4_beam.Vect()).Angle(V3_el_prot[count][j])*TMath::RadToDeg(),-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       h2_Ecal_Eqe->Fill(E_rec_new,E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
		       h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);

		       for (int n=0;n<N_pperp;n++){
			 for(int m=0;m<N_Ecal;m++){
			   if(E_cal_4pto3p[count][j]>Ecal_lowlim[m] && E_cal_4pto3p[count][j]<Ecal_uplim[m] && p_miss_perp_4pto3p[count][j]>pperp_cut[n])   h1_Etot_p_bkgd_slice_Ecalcut4321[n][m]->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);

			 }
		       }


		       for(int i=0;i<n_slice;i++)
			 {
			   if (p_miss_perp_4pto3p[count][j]<pperp_max[i] && p_miss_perp_4pto3p[count][j]>pperp_min[i]){

			     h1_Etot_4pto3p_slice[i]->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
			     h1_Erec_4pto3p_slice[i]->Fill(E_rec_new, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*1/Mott_cross_sec);
			   }
			 }
		     }
		     }

		   //-----------------------------------------  4p to 3p->1  -----------------------------------------------------------------------
		   	   for(int j=0;j<N_3p;j++)    { //4p to 3p->1, looping through 1p out of 3p

		     P_43pto1p[j]= N_p1[j]/N_p_three*(N_p4_p3[g]/N_p_four);
		     Etot_43pto1p->Fill(E_cal_43pto1p[j], P_43pto1p[j]*1/Mott_cross_sec);
		     Erec_43pto1p->Fill(E_rec_new,P_43pto1p[j]*1/Mott_cross_sec);
		     Etot_p431_bkgd09->Fill(E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
		     h2_Erec_pperp_431p->Fill(p_miss_perp_43pto1p[j],E_rec_new,P_43pto1p[j]*1/Mott_cross_sec);
		     Etot_43pto1p_fracfeed->Fill((E_cal_43pto1p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_43pto1p[j]*1/Mott_cross_sec);
		     Erec_43pto1p_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_43pto1p[j]*1/Mott_cross_sec);
		     h2_pperp_W->Fill(W_var,p_miss_perp_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
		     h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_43pto1p[j]*1/Mott_cross_sec);
		     h2_Ecal_Eqe->Fill(E_rec_new,E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
		     h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
		     h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);

		     for (int n=0;n<N_pperp;n++){
		       for(int m=0;m<N_Ecal;m++){
			 if(E_cal_43pto1p[j]>Ecal_lowlim[m] && E_cal_43pto1p[j]<Ecal_uplim[m] && p_miss_perp_43pto1p[j]>pperp_cut[n])   h1_Etot_p_bkgd_slice_Ecalcut431[n][m]->Fill(E_cal_43pto1p[j], P_43pto1p[j]*1/Mott_cross_sec);

		       }
		     }

		     for(int i=0;i<n_slice;i++)
		       {
			 if (p_miss_perp_43pto1p[j]<pperp_max[i] && p_miss_perp_43pto1p[j]>pperp_min[i]){

			   h1_Etot_43pto1p_slice[i]->Fill(E_cal_43pto1p[j],P_43pto1p[j]*1/Mott_cross_sec);
			   h1_Erec_43pto1p_slice[i]->Fill(E_rec_new,P_43pto1p[j]*1/Mott_cross_sec);
			 }
		       }

		   }

		 }//end of N_p_three and N_p_four !=0
	       }//end of the loop through 3p combinations out of 4





	       int N_4to2=0;
	       TVector3 V3p2[2],V3p2_uncorr[2];
	       double E_cal_4pto2p[2]={0}, p_miss_perp_4pto2p[2]={0},P_4pto2p[2]={0}, N_two=0;

  //-----------------------------------------  4p to 2p->1  -----------------------------------------------------------------------
	       for(int ind1=0;ind1<N_p4;ind1++){          //estimating the undetected 4p contribution to  2p
		 for(int ind2=0;ind2<N_p4;ind2++){
		   if(ind1!=ind2 && ind1<ind2){

		     V3p2[0]=V3_prot4_corr[ind1];
		     V3p2[1]=V3_prot4_corr[ind2];
		     V3p2_uncorr[0]=V3_prot4[ind1];
		     V3p2_uncorr[1]=V3_prot4[ind2];

		     prot2_rot_func(fbeam_en, V3_q, V3p2,V3p2_uncorr, V4_el,E_cal_4pto2p,p_miss_perp_4pto2p,  P_4pto2p, &N_two);

		     if( N_two!=0  && N_p_four!=0){

		       for(int j=0;j<N_2p;j++)    {  //looping through  1 proton combination out of 2 protons

			 Etot_4pto2p->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 Erec_4pto2p->Fill(E_rec_new, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 Etot_p421_bkgd09->Fill(E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 h2_Erec_pperp_421p->Fill( p_miss_perp_4pto2p[j],E_rec_new,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 Etot_4pto2p_fracfeed->Fill((E_cal_4pto2p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 Erec_4pto2p_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 h2_pperp_W->Fill(W_var,p_miss_perp_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3p2_uncorr[j])*TMath::RadToDeg(),P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 h2_Ecal_Eqe->Fill(E_rec_new,E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			 h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);


			 for (int n=0;n<N_pperp;n++){
			   for(int m=0;m<N_Ecal;m++){
			     if(E_cal_4pto2p[j]>Ecal_lowlim[m] && E_cal_4pto2p[j]<Ecal_uplim[m] && p_miss_perp_4pto2p[j]>pperp_cut[n])   h1_Etot_p_bkgd_slice_Ecalcut421[n][m]->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);

			   }
			 }

			 for(int i=0;i<n_slice;i++)
			   {
			     if (p_miss_perp_4pto2p[j]<pperp_max[i] && p_miss_perp_4pto2p[j]>pperp_min[i]){

			       h1_Etot_4pto2p_slice[i]->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			       h1_Erec_4pto2p_slice[i]->Fill(E_rec_new, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*1/Mott_cross_sec);
			     }
			   }
		       }



		     }
		     N_4to2= N_4to2+1;

		   }
		 }
	       }

 //-----------------------------------------  4p to 1p  -----------------------------------------------------------------------
	       if(N_p_four!=0){
		 for(int j=0;j<N_p4;j++)    {       //estimating the undetected 4p contribution to  1p

		   P_4pto1p[j]= N_p4_p1[j]/N_p_four;
		   Etot_4pto1p->Fill(E_cal_p4[j], P_4pto1p[j]*1/Mott_cross_sec);
		   Erec_4pto1p->Fill(E_rec_new,P_4pto1p[j]*1/Mott_cross_sec);
		   Etot_p41_bkgd09->Fill(E_cal_p4[j],P_4pto1p[j]*1/Mott_cross_sec);
		   h2_Erec_pperp_41p->Fill(p_miss_perp_p4[j],E_rec_new, P_4pto1p[j]*1/Mott_cross_sec);
		   Etot_4pto1p_fracfeed->Fill((E_cal_p4[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto1p[j]*1/Mott_cross_sec);
		   Erec_4pto1p_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_4pto1p[j]*1/Mott_cross_sec);
		   h2_pperp_W->Fill(W_var,p_miss_perp_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);
		   h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot4[j])*TMath::RadToDeg(),-P_4pto1p[j]*1/Mott_cross_sec);
		   h2_Ecal_Eqe->Fill(E_rec_new,E_cal_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);
		   h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_cal_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);
		   h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_cal_p4[j],-P_4pto1p[j]*1/Mott_cross_sec);

		   for (int n=0;n<N_pperp;n++){
		     for(int m=0;m<N_Ecal;m++){
		       if(E_cal_p4[j]>Ecal_lowlim[m] && E_cal_p4[j]<Ecal_uplim[m] && p_miss_perp_p4[j]>pperp_cut[n])   h1_Etot_p_bkgd_slice_Ecalcut41[n][m]->Fill(E_cal_p4[j], P_4pto1p[j]*1/Mott_cross_sec);

		     }
		   }

		   for(int i=0;i<n_slice;i++)
		     {
		       if (p_miss_perp_p4[j]<pperp_max[i] && p_miss_perp_p4[j]>pperp_min[i]){

			 h1_Etot_4pto1p_slice[i]->Fill(E_cal_p4[j],P_4pto1p[j]*1/Mott_cross_sec);
			 h1_Erec_4pto1p_slice[i]->Fill(E_rec_new,P_4pto1p[j]*1/Mott_cross_sec);
		       }
		     }

		 }
	       }




	     }//no pion statment ends


	   }//4 proton vertex cut
	 }//4 proton requirement



  double P_undet=0;
 TVector3 V3_pi;

	 V3_q=(V4_beam-V4_el).Vect();
	h1_E_rec->Fill(E_rec_new,1/Mott_cross_sec);

	if (ec_num_n==0 && num_pi_phot==0){

	   h1_E_rec_0pi->Fill(E_rec_new,1/Mott_cross_sec);
	   h1_E_rec_0pi_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],1/Mott_cross_sec);

	 }

//----------------------------- e- ,1pi  -----------------------------------------

	 if (ec_num_n==0 && num_pi_phot==1) {

	   V3_pi.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
	   pi1_rot_func(fbeam_en, V3_pi,q[ind_pi_phot[0]],ec_radstat_n[0],V3_q, &P_undet);

	   h1_E_rec_1pi_weight->Fill(E_rec_new,P_undet*1/Mott_cross_sec);
	   h1_E_rec_1pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_undet*1/Mott_cross_sec);

	  if(!ec_radstat_n[0])  h1_E_rec_1pi->Fill(E_rec_new,1/Mott_cross_sec);
	   if(ec_num_n==1) h2_phot_e_angle_Erec->Fill(E_rec_new,V3_pi.Angle(V4_el.Vect())*TMath::RadToDeg());
	 }
//----------------------------- e- ,2pi  -----------------------------------------

	 const int N_2pi=2;
	 TVector3 V3_2pi[N_2pi];
	 int pi2_q[N_2pi];
	 double P_1pi[N_2pi]={0},P_0pi=0;
	 bool pi2_radstat[N_2pi]={false};
if (ec_num_n==0 && num_pi_phot==2) {

  V3_2pi[0].SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
  V3_2pi[1].SetXYZ(p[ind_pi_phot[1]]*cx[ind_pi_phot[1]],p[ind_pi_phot[1]]*cy[ind_pi_phot[1]],p[ind_pi_phot[1]]*cz[ind_pi_phot[1]]);
  pi2_q[0]=q[ind_pi_phot[0]];
  pi2_q[1]=q[ind_pi_phot[1]];
  pi2_radstat[0]=ec_radstat_n[0];
  pi2_radstat[1]=ec_radstat_n[1];

  pi2_rot_func(fbeam_en, V3_2pi, pi2_q,pi2_radstat, V3_q,&P_0pi,P_1pi);
 //----------------------------- e- ,2pi->0pi (-) -----------------------------------------
 h1_E_rec_2pi_weight->Fill(E_rec_new,(-P_0pi)*1/Mott_cross_sec);
 h1_E_rec_2pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*1/Mott_cross_sec);
h1_E_rec_20pi->Fill(E_rec_new,(P_0pi)*1/Mott_cross_sec);
//----------------------------- e- ,2pi->1pi->0pi (+)  -----------------------------------------
 for(int k=0;k<N_2pi;k++){
 h1_E_rec_2pi_weight->Fill(E_rec_new,P_1pi[k]*1/Mott_cross_sec);
 h1_E_rec_2pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[k]*1/Mott_cross_sec);
h1_E_rec_21pi->Fill(E_rec_new,(P_1pi[k])*1/Mott_cross_sec);
 }
}


//----------------------------- e- ,3pi  -----------------------------------------
 const int N_3pi=3;
TVector3 V3_3pi[N_3pi];
 int q_pi3[N_3pi]={0};
 bool radstat_pi3[N_3pi]={false};
 double P_0pion=0,P_1pion[N_3pi]={0},P_320pi[N_3pi]={0},P_3210pi[N_3pi][N_2pi]={0};

if (ec_num_n==0 && num_pi_phot==3) {
  for (int h=0;h<N_3pi;h++){

    V3_3pi[h].SetXYZ(p[ind_pi_phot[h]]*cx[ind_pi_phot[h]],p[ind_pi_phot[h]]*cy[ind_pi_phot[h]],p[ind_pi_phot[h]]*cz[ind_pi_phot[h]]);
    q_pi3[h]=q[ind_pi_phot[h]];
    radstat_pi3[h]=ec_radstat_n[h];
  }
  pi3_rot_func(fbeam_en, V3_3pi, q_pi3,radstat_pi3,V3_q,&P_0pion, P_1pion, P_320pi,P_3210pi);


 //---------------------------3pi->0pi----------------------------------------------
    h1_E_rec_3pi_weight->Fill(E_rec_new,(-P_0pion)*1/Mott_cross_sec);
    h1_E_rec_3pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pion)*1/Mott_cross_sec);
h1_E_rec_30pi->Fill(E_rec_new,(P_0pion)*1/Mott_cross_sec);

 //---------------------------3pi->1pi->0pi----------------------------------------------
    for(int h=0;h<N_3pi;h++){

  h1_E_rec_3pi_weight->Fill(E_rec_new,P_1pion[h]*1/Mott_cross_sec);
    h1_E_rec_3pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pion[h]*1/Mott_cross_sec);
h1_E_rec_310pi->Fill(E_rec_new,(P_1pion[h])*1/Mott_cross_sec);

 //---------------------------3pi->2pi->0pi----------------------------------------------
    h1_E_rec_3pi_weight->Fill(E_rec_new,P_320pi[h]*1/Mott_cross_sec);
    h1_E_rec_3pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_320pi[h]*1/Mott_cross_sec);
    h1_E_rec_320pi->Fill(E_rec_new,(P_320pi[h])*1/Mott_cross_sec);

//---------------------------3pi->2pi->1pi->0pi----------------------------------------------
    for(int g=0;g<N_2pi;g++){
    h1_E_rec_3pi_weight->Fill(E_rec_new,(-P_3210pi[h][g])*1/Mott_cross_sec);
    h1_E_rec_3pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_3210pi[h][g])*1/Mott_cross_sec);
    h1_E_rec_3210pi->Fill(E_rec_new,(P_3210pi[h][g])*1/Mott_cross_sec);
    }

    }//end of 3pi loop

 }//end of 3pi requirement



//----------------------------- e- ,4pi  -----------------------------------------
if (ec_num_n==0 && num_pi_phot==4) {

 const int N_4pi=4;
TVector3 V3_4pi[N_4pi];
 int q_pi4[N_4pi]={0};
 bool  radstat_pi4[N_4pi]={false};
 double P_0pion=0,P_410pi=0,P_420pi=0,P_4210pi=0,P_430pi=0,P_4310pi=0,P_4320pi=0,P_43210pi=0;

  for (int h=0;h<N_4pi;h++){

    V3_4pi[h].SetXYZ(p[ind_pi_phot[h]]*cx[ind_pi_phot[h]],p[ind_pi_phot[h]]*cy[ind_pi_phot[h]],p[ind_pi_phot[h]]*cz[ind_pi_phot[h]]);
    q_pi4[h]=q[ind_pi_phot[h]];
    radstat_pi4[h]=ec_radstat_n[h];
  }

  pi4_rot_func(fbeam_en, V3_4pi, q_pi4,radstat_pi4,V3_q,&P_0pion,&P_410pi,&P_420pi,&P_4210pi,&P_430pi,&P_4310pi,&P_4320pi,&P_43210pi);



  //  cout<<"verev "<<P_0pion<<"  "<<P_410pi<<"  "<<P_420pi<<"  "<<P_4210pi<<"  "<<P_430pi<<"  "<<P_4310pi<<"  "<<P_4320pi<<"  "<<P_43210pi<<"  "<<endl;
 //---------------------------4pi->0pi----------------------------------------------

    h1_E_rec_4pi_weight->Fill(E_rec_new,(-P_0pion+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*1/Mott_cross_sec);
    h1_E_rec_4pi_weight_frac_feed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pion+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*1/Mott_cross_sec);
h1_E_rec_40pi->Fill(E_rec_new,(P_0pion)*1/Mott_cross_sec);

//---------------------------4pi->1pi->0pi----------------------------------------------

h1_E_rec_410pi->Fill(E_rec_new,(P_410pi)*1/Mott_cross_sec);
//---------------------------4pi->2pi->0pi----------------------------------------------
h1_E_rec_420pi->Fill(E_rec_new,(P_420pi)*1/Mott_cross_sec);
//---------------------------4pi->2pi->1pi->0pi----------------------------------------------
h1_E_rec_4210pi->Fill(E_rec_new,(P_4210pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->0pi----------------------------------------------
h1_E_rec_430pi->Fill(E_rec_new,(P_430pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->1pi->0pi----------------------------------------------
h1_E_rec_4310pi->Fill(E_rec_new,(P_4310pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->2pi->0pi----------------------------------------------
h1_E_rec_4320pi->Fill(E_rec_new,(P_4320pi)*1/Mott_cross_sec);
//---------------------------4pi->3pi->2pi->1pi->0pi----------------------------------------------
h1_E_rec_43210pi->Fill(E_rec_new,(P_43210pi)*1/Mott_cross_sec);

 }//end of 4 pi/photon requirement




   //------------------------------------------requiring there to be a proton -------------------------------------

   if( num_p==1)
     {


       ind_p = index_p[0];

       TLorentzVector V4_prot_uncorr(p[ind_p]*cx[ind_p],p[ind_p]*cy[ind_p],p[ind_p]*cz[ind_p],TMath::Sqrt(m_prot*m_prot+p[ind_p]*p[ind_p]));


       float prot_vert=vz[ind_p];
       double p_phi=TMath::ATan2(cy[ind_p],cx[ind_p])*TMath::RadToDeg();
       double p_phi_mod=p_phi+30;
       if (p_phi_mod<0)p_phi_mod=p_phi_mod+360;
       double	p_theta=TMath::ACos(cz[ind_p])*TMath::RadToDeg();
       double prot_vert_corr=prot_vert+vz_corr(p_phi_mod,p_theta);
       double prot_mom_corr;

       h1_el_prot_vertdiff->Fill(el_vert_corr-prot_vert_corr);


       if(ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr,prot_vert_corr) != -1){
	 prot_mom_corr=ProtonMomCorrection_He3_4Cell(ftarget,V4_prot_uncorr,prot_vert_corr);}
       else
	 {prot_mom_corr=p[ind_p];}

       h1_prot_mom->Fill(prot_mom_corr);
       h1_prot_mom_ratio->Fill(prot_mom_corr/p[ind_p]);

       TLorentzVector V4_prot_corr(prot_mom_corr*cx[ind_p],prot_mom_corr*cy[ind_p],prot_mom_corr*cz[ind_p],TMath::Sqrt(m_prot*m_prot+prot_mom_corr*prot_mom_corr));
       TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();


       TLorentzVector V4_prot_el_tot=V4_prot_corr+V4_el;
       double p_perp_tot=TMath::Sqrt(V4_prot_el_tot.Px()*V4_prot_el_tot.Px()+V4_prot_el_tot.Py()*V4_prot_el_tot.Py());
       double p_z_tot=V4_prot_el_tot.Pz();
       double p_tot=V4_prot_el_tot.Rho();
       double E_tot=V4_el.E()+V4_prot_corr.E()-m_prot+bind_en[ftarget];


       if((el_vert_corr-prot_vert_corr)>vertdiff_min[ftarget] && (el_vert_corr-prot_vert_corr)<vertdiff_max[ftarget]){

 //---------------------------------- 1p 2pi   ----------------------------------------------
    if (ec_num_n==0 && num_pi_phot==2) {

      const int N_2pi=2;
      TVector3 V3_2pi[N_2pi],V3_2pi_rot[N_2pi],V3_p_rot;
      bool pi2_stat[N_2pi]={false};
      int q_pi2[N_2pi];
      bool radstat_pi2[N_2pi]={false};
      V3_q=(V4_beam-V4_el).Vect();

      V3_2pi[0].SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);
      V3_2pi[1].SetXYZ(p[ind_pi_phot[1]]*cx[ind_pi_phot[1]],p[ind_pi_phot[1]]*cy[ind_pi_phot[1]],p[ind_pi_phot[1]]*cz[ind_pi_phot[1]]);
      q_pi2[0]=q[ind_pi_phot[0]];
      q_pi2[1]=q[ind_pi_phot[1]];
      radstat_pi2[0]=ec_radstat_n[0];
      radstat_pi2[1]=ec_radstat_n[1];

      double P_1p0pi=0,P_1p1pi[N_2pi]={0};
      prot1_pi2_rot_func(fbeam_en, V3_q,V3_prot_uncorr,V3_2pi,q_pi2,radstat_pi2,&P_1p0pi,P_1p1pi);

 //---------------------------------- 1p 2pi->1p1pi   ----------------------------------------------

   for(int z=0;z<N_2pi;z++){  //to consider 2 diff. 1pi states

       Etot_1p2pi->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
       Erec_1p2pi->Fill(E_rec_new,P_1p1pi[z]*1/Mott_cross_sec);
       h2_Erec_pperp_1p2pi_1p1pi->Fill(p_perp_tot,E_rec_new,P_1p1pi[z]*1/Mott_cross_sec);
       Etot_bkgd09_1p2pi_1p1pi->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
       Etot_1p2pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*1/Mott_cross_sec);
       Erec_1p2pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p1pi[z]*1/Mott_cross_sec);
       h2_pperp_W->Fill(W_var,p_perp_tot,P_1p1pi[z]*1/Mott_cross_sec);
       h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p1pi[z]*1/Mott_cross_sec);
       h2_Ecal_Eqe->Fill(E_rec_new,E_tot,P_1p1pi[z]*1/Mott_cross_sec);
       h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot,P_1p1pi[z]*1/Mott_cross_sec);
       h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot,P_1p1pi[z]*1/Mott_cross_sec);

       for(int i=0;i<n_slice;i++){
	 if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	   h1_Etot_bkgd_1p2pi[i]->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
	   h1_Erec_bkgd_1p2pi[i]->Fill(E_rec_new,P_1p1pi[z]*1/Mott_cross_sec);
	 }
       }
       for (int i=0;i<N_pperp;i++){
	 for(int j=0;j<N_Ecal;j++){
	   if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[i][j]->Fill(E_tot,P_1p1pi[z]*1/Mott_cross_sec);
	 }
       }
   }

 //---------------------------------- 1p 2pi->1p0pi   ----------------------------------------------

   Etot_1p2pi_1p0pi->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
      Erec_1p2pi_1p0pi->Fill(E_rec_new,P_1p0pi*1/Mott_cross_sec);
      h2_Erec_pperp_1p2pi_1p0pi->Fill(p_perp_tot,E_rec_new,P_1p0pi*1/Mott_cross_sec);
      Etot_bkgd09_1p2pi_1p0pi->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
      Etot_1p2pi_1p0pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p0pi*1/Mott_cross_sec);
      Erec_1p2pi_1p0pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p0pi*1/Mott_cross_sec);
      h2_pperp_W->Fill(W_var,p_perp_tot,-P_1p0pi*1/Mott_cross_sec);
      h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-P_1p0pi*1/Mott_cross_sec);
      h2_Ecal_Eqe->Fill(E_rec_new,E_tot,-P_1p0pi*1/Mott_cross_sec);
      h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot,-P_1p0pi*1/Mott_cross_sec);
      h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot,-P_1p0pi*1/Mott_cross_sec);

      for(int i=0;i<n_slice;i++){
      if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
	h1_Etot_bkgd_1p2pi_1p0pi[i]->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
	h1_Erec_bkgd_1p2pi_1p0pi[i]->Fill(E_rec_new,P_1p0pi*1/Mott_cross_sec);
      }
    }
    for (int i=0;i<N_pperp;i++){
      for(int j=0;j<N_Ecal;j++){
	if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[i][j]->Fill(E_tot,P_1p0pi*1/Mott_cross_sec);
      }
    }

 }//1p 2pi statetment ends






 //---------------------------------- 1p 3pi   ----------------------------------------------

 if (ec_num_n==0 && num_pi_phot==3) {

   const int N_3pi=3;
   TVector3 V3_3pi[N_3pi],V3_3pi_rot[N_3pi],V3_p_rot;
   bool pi3_stat[N_3pi]={false};
   int q_pi3[N_3pi];
   bool radstat_pi3[N_3pi]={false};
   double N_3pi_p=0,N_1pi_1p[N_3pi]={0},N_pi=0,N_nopi=0,N_0pi_1p=0;
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

   double P_1p3pi=0;

   prot1_pi3_rot_func(fbeam_en, V3_q,V3_prot_uncorr,V3_3pi,q_pi3,radstat_pi3,&P_1p3pi);
   Etot_1p3pi->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
   Erec_1p3pi->Fill(E_rec_new,P_1p3pi*1/Mott_cross_sec);
   h2_Erec_pperp_1p3pi->Fill(p_perp_tot,E_rec_new,P_1p3pi*1/Mott_cross_sec);
   Etot_bkgd09_1p3pi->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
   Etot_1p3pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi*1/Mott_cross_sec);
   Erec_1p3pi_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p3pi*1/Mott_cross_sec);
   h2_pperp_W->Fill(W_var,p_perp_tot,P_1p3pi*1/Mott_cross_sec);
   h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p3pi*1/Mott_cross_sec);
   h2_Ecal_Eqe->Fill(E_rec_new,E_tot,P_1p3pi*1/Mott_cross_sec);
   h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot,P_1p3pi*1/Mott_cross_sec);
   h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot,P_1p3pi*1/Mott_cross_sec);

   for(int i=0;i<n_slice;i++){
     if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){
       h1_Etot_bkgd_1p3pi[i]->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
       h1_Erec_bkgd_1p3pi[i]->Fill(E_rec_new,P_1p3pi*1/Mott_cross_sec);
     }
   }
   for (int i=0;i<N_pperp;i++){
     for(int j=0;j<N_Ecal;j++){
       if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i])   h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[i][j]->Fill(E_tot,P_1p3pi*1/Mott_cross_sec);
     }
   }

 }//end of 1p 3pi requirement



  //---------------------------------- 1p 1pi   ----------------------------------------------
 double N_piphot_det,N_piphot_undet;
 TVector3 V3_pi_phot;
 V3_q=(V4_beam-V4_el).Vect();

 if (ec_num_n==0 && num_pi_phot==1) {
   V3_pi_phot.SetXYZ(p[ind_pi_phot[0]]*cx[ind_pi_phot[0]],p[ind_pi_phot[0]]*cy[ind_pi_phot[0]],p[ind_pi_phot[0]]*cz[ind_pi_phot[0]]);	N_piphot_det=N_piphot_undet=0;
   prot1_pi1_rot_func(fbeam_en, V3_q,V3_prot_uncorr,V3_pi_phot, q[ind_pi_phot[0]],ec_radstat_n[0], &N_piphot_det,&N_piphot_undet);
 }

	 Erec_1prot->Fill(E_rec_new,1/Mott_cross_sec);
	 Etot_1prot->Fill(E_tot,1/Mott_cross_sec);
	 h2_Erec_pperp->Fill(p_perp_tot,E_rec_new,1/Mott_cross_sec);

	 for(int i=0;i<n_slice;i++)
	   {
	     if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){

	       if (ec_num_n==0 && num_pi_phot==0)	   {

		 h1_Etot_Npi0[i]->Fill(E_tot,1/Mott_cross_sec);
		 h1_Erec_Npi0_new[i]->Fill(E_rec_new,1/Mott_cross_sec);
	       }

	       if (ec_num_n==0 && num_pi_phot==1) {

		 if(N_piphot_det!=0){

		   h1_E_rec_undetfactor->Fill(E_rec_new,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_E_tot_undetfactor->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_E_tot_undetfactor_pipl->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_Etot_bkgd_pipl_pimi_fact[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_Erec_bkgd_pipl_pimi_new_fact[i]->Fill(E_rec_new,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_Etot_bkgd_pipl_pimi_fact_pipl[i]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_E_tot_undetfactor09->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h2_Erec_pperp_1p1pi->Fill(p_perp_tot,E_rec_new,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_E_rec_undetfactor_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_E_tot_undetfactor_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h2_pperp_W->Fill(W_var,p_perp_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h2_Ecal_Eqe->Fill(E_rec_new,E_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);
		   h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot,-(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);

		   for (int i=0;i<N_pperp;i++){
		     for(int j=0;j<N_Ecal;j++){
		       if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i])   h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j]->Fill(E_tot,(N_piphot_undet/N_piphot_det)*1/Mott_cross_sec);

		     }
		   }
		 }

	       h1_E_rec_cutpi1_piplpimi->Fill(E_rec,1/Mott_cross_sec);
	       h1_E_tot_cutpi1_piplpimi->Fill(E_tot,1/Mott_cross_sec);
	       h1_Etot_Npi1[i]->Fill(E_tot,1/Mott_cross_sec);
	       h1_Erec_Npi1[i]->Fill(E_rec_new,1/Mott_cross_sec);

	     }


	     }//pperp slices
	   }//for loop for pperp slices

	 h1_Etot->Fill(E_tot,1/Mott_cross_sec);


   if (ec_num_n==0 && num_pi_phot==0){
	     h2_Erec_pperp_newcut2->Fill(p_perp_tot,E_rec_new,1/Mott_cross_sec);
	     h1_E_rec_cut2_new->Fill(E_rec_new,1/Mott_cross_sec);
	     h1_E_tot_cut2->Fill(E_tot,1/Mott_cross_sec);
	     h1_E_tot_cut2_09->Fill(E_tot,1/Mott_cross_sec);
	     h1_E_tot_cut2_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],1/Mott_cross_sec);
	     h1_E_rec_cut2_new_fracfeed->Fill((E_rec_new-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],1/Mott_cross_sec);
	     h2_pperp_W->Fill(W_var,p_perp_tot,1/Mott_cross_sec);
	     h1_theta0->Fill((V4_beam.Vect()).Angle(V4_prot_el_tot.Vect()) *TMath::RadToDeg(),1/Mott_cross_sec);
	     h2_Ecal_Eqe->Fill(E_rec_new,E_tot,1/Mott_cross_sec);
	     h2_EqeEcalratio_Eqe->Fill(E_rec_new,E_rec_new/E_tot,1/Mott_cross_sec);
	     h2_EqeEcaldiff_Eqe->Fill(E_rec_new,E_rec_new-E_tot,1/Mott_cross_sec);

	     for (int i=0;i<N_pperp;i++){
	       for(int j=0;j<N_Ecal;j++){
		 if(E_tot>Ecal_lowlim[j] && E_tot<Ecal_uplim[j] && p_perp_tot>pperp_cut[i])    h1_Etot_Npi0_Ecalcut[i][j]->Fill(E_tot,1/Mott_cross_sec);

	       }
	     }

	     if (p_perp_tot<0.2){

	       h1_E_rec_cut005_newcut3->Fill(E_rec_new,1/Mott_cross_sec);
	       h2_Erec_pperp_cut3->Fill(p_perp_tot,E_rec,1/Mott_cross_sec);
	     }


	 }//num pi=0

       }//ep vert

     }//1proton ends




      } //vertex cut



   }


 gStyle->SetOptFit(1);


 for(int i=0;i<=n_slice-1;i++)
      {


  //------------------------------------using the ratio of the pi- to pi+  --------------------------------------

	h_Etot_piplpimi_subtruct_fact[i]=(TH1F*)  h1_Etot_Npi0[i]->Clone(Form("h_Etot_piplpimi_subtruct_fact_%d",i+1));
	h_Etot_piplpimi_subtruct_fact[i]->Add(h1_Etot_bkgd_pipl_pimi_fact[i],-1);
	h_Erec_piplpimi_subtruct_new_fact[i]=(TH1F*)  h1_Erec_Npi0_new[i]->Clone(Form("h_Erec_piplpimi_subtruct_new_fact_%d",i+1));
	h_Erec_piplpimi_subtruct_new_fact[i]->Add(h1_Erec_bkgd_pipl_pimi_new_fact[i],-1);

 //------------------------------------subtracting 2p contribution from 1p events  --------------------------------------

	h1_Etot_p_bkgd_slice_sub[i]=(TH1F*) h_Etot_piplpimi_subtruct_fact[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub_%d",i+1));
	h1_Etot_p_bkgd_slice_sub[i]->Add(h1_Etot_p_bkgd_slice[i],-1);
	h1_Erec_p_bkgd_slice_sub[i]=(TH1F*) h_Erec_piplpimi_subtruct_new_fact[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub_%d",i+1));
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

  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){


      //------------------------------------using the ratio of the pi- to pi+  --------------------------------------
      h_Etot_piplpimi_subtruct_fact_Ecalcut[i][j]=(TH1F*)  h1_Etot_Npi0_Ecalcut[i][j]->Clone(Form("h_Etot_piplpimi_subtruct_fact_Ecalcut_%d_%d",i+1,i+j));
          h_Etot_piplpimi_subtruct_fact_Ecalcut[i][j]->Add(h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j],-1);

      //------------------------------------subtracting 2p contribution from 1p events  --------------------------------------

      h1_Etot_p_bkgd_slice_sub_Ecalcut[i][j]=(TH1F*) h_Etot_piplpimi_subtruct_fact_Ecalcut[i][j]->Clone(Form("h1_Etot_p_bkgd_slice_sub_Ecalcut_%d_%d",i+1,j+1));
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

    }
  }






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
  h_Erec_subtruct_piplpimi_prot->Add(Erec_p_bkgd,-1);

TH1F *h_Etot_subtruct_piplpimi_prot=(TH1F*)  h_Etot_subtruct_piplpimi_factor->Clone("h_Etot_subtruct_piplpimi_prot");
  h_Etot_subtruct_piplpimi_prot->Add(Etot_p_bkgd,-1);

TH1F *h_Etot_subtruct_piplpimi_prot09=(TH1F*)  h_Etot_subtruct_piplpimi_factor09->Clone("h_Etot_subtruct_piplpimi_prot09");
  h_Etot_subtruct_piplpimi_prot09->Add(Etot_p_bkgd09,-1);

 TH2F *h2_Erec_pperp_2psub=(TH2F*) h2_Erec_pperp_1p1pisub->Clone("h2_Erec_pperp_2psub");
 h2_Erec_pperp_2psub->Add(h2_Erec_pperp_2p,-1);

TH1F *h_Erec_subtruct_piplpimi_prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_factor_fracfeed->Clone("h_Erec_subtruct_piplpimi_prot_fracfeed");
  h_Erec_subtruct_piplpimi_prot_fracfeed->Add(Erec_p_bkgd_fracfeed,-1);

TH1F *h_Etot_subtruct_piplpimi_prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_factor_fracfeed->Clone("h_Etot_subtruct_piplpimi_prot_fracfeed");
  h_Etot_subtruct_piplpimi_prot_fracfeed->Add(Etot_p_bkgd_fracfeed,-1);

 //-----------------------------------undetected 3 to 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_32prot=(TH1F*)  h_Erec_subtruct_piplpimi_prot->Clone("h_Erec_subtruct_piplpimi_32prot");
  h_Erec_subtruct_piplpimi_32prot->Add(Erec_3pto2p);

  TH1F *h_Etot_subtruct_piplpimi_32prot=(TH1F*)  h_Etot_subtruct_piplpimi_prot->Clone("h_Etot_subtruct_piplpimi_32prot");
  h_Etot_subtruct_piplpimi_32prot->Add(Etot_3pto2p);

TH1F *h_Etot_subtruct_piplpimi_32prot09=(TH1F*)  h_Etot_subtruct_piplpimi_prot09->Clone("h_Etot_subtruct_piplpimi_32prot09");
  h_Etot_subtruct_piplpimi_32prot09->Add(Etot_p321_bkgd09);

 TH2F *h2_Erec_pperp_32psub=(TH2F*) h2_Erec_pperp_2psub->Clone("h2_Erec_pperp_32psub");
 h2_Erec_pperp_32psub->Add(h2_Erec_pperp_321p);

  TH1F *h_Erec_subtruct_piplpimi_32prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_32prot_fracfeed");
  h_Erec_subtruct_piplpimi_32prot_fracfeed->Add(Erec_3pto2p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_32prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_32prot_fracfeed");
  h_Etot_subtruct_piplpimi_32prot_fracfeed->Add(Etot_3pto2p_fracfeed);

 //-----------------------------------undetected 3 to 1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_31prot=(TH1F*)  h_Erec_subtruct_piplpimi_32prot->Clone("h_Erec_subtruct_piplpimi_31prot");
  h_Erec_subtruct_piplpimi_31prot->Add(Erec_3pto1p,-1);

  TH1F *h_Etot_subtruct_piplpimi_31prot=(TH1F*)  h_Etot_subtruct_piplpimi_32prot->Clone("h_Etot_subtruct_piplpimi_31prot");
  h_Etot_subtruct_piplpimi_31prot->Add(Etot_3pto1p,-1);

TH1F *h_Etot_subtruct_piplpimi_31prot09=(TH1F*)  h_Etot_subtruct_piplpimi_32prot09->Clone("h_Etot_subtruct_piplpimi_31prot09");
  h_Etot_subtruct_piplpimi_31prot09->Add(Etot_p31_bkgd09,-1);

 TH2F *h2_Erec_pperp_31psub=(TH2F*) h2_Erec_pperp_32psub->Clone("h2_Erec_pperp_31psub");
 h2_Erec_pperp_31psub->Add(h2_Erec_pperp_31p,-1);

 TH1F *h_Erec_subtruct_piplpimi_31prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_32prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_31prot_fracfeed");
  h_Erec_subtruct_piplpimi_31prot_fracfeed->Add(Erec_3pto1p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_31prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_32prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_31prot_fracfeed");
  h_Etot_subtruct_piplpimi_31prot_fracfeed->Add(Etot_3pto1p_fracfeed,-1);

 //-----------------------------------undetected 4 to 3->2->1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_43prot=(TH1F*)  h_Erec_subtruct_piplpimi_31prot->Clone("h_Erec_subtruct_piplpimi_43prot");
  h_Erec_subtruct_piplpimi_43prot->Add(Erec_4pto3p,-1);

  TH1F *h_Etot_subtruct_piplpimi_43prot=(TH1F*)  h_Etot_subtruct_piplpimi_31prot->Clone("h_Etot_subtruct_piplpimi_43prot");
  h_Etot_subtruct_piplpimi_43prot->Add(Etot_4pto3p,-1);

TH1F *h_Etot_subtruct_piplpimi_43prot09=(TH1F*)  h_Etot_subtruct_piplpimi_31prot09->Clone("h_Etot_subtruct_piplpimi_43prot09");
  h_Etot_subtruct_piplpimi_43prot09->Add(Etot_p4321_bkgd09,-1);

 TH2F *h2_Erec_pperp_43psub=(TH2F*) h2_Erec_pperp_31psub->Clone("h2_Erec_pperp_43psub");
 h2_Erec_pperp_43psub->Add(h2_Erec_pperp_4321p,-1);

 TH1F *h_Erec_subtruct_piplpimi_43prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_31prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_43prot_fracfeed");
  h_Erec_subtruct_piplpimi_43prot_fracfeed->Add(Erec_4pto3p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_43prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_31prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_43prot_fracfeed");
  h_Etot_subtruct_piplpimi_43prot_fracfeed->Add(Etot_4pto3p_fracfeed,-1);

 //-----------------------------------undetected 4 to 3->1 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_431prot=(TH1F*)  h_Erec_subtruct_piplpimi_43prot->Clone("h_Erec_subtruct_piplpimi_431prot");
  h_Erec_subtruct_piplpimi_431prot->Add(Erec_43pto1p);

  TH1F *h_Etot_subtruct_piplpimi_431prot=(TH1F*)  h_Etot_subtruct_piplpimi_43prot->Clone("h_Etot_subtruct_piplpimi_431prot");
  h_Etot_subtruct_piplpimi_431prot->Add(Etot_43pto1p);

  TH1F *h_Etot_subtruct_piplpimi_431prot09=(TH1F*)  h_Etot_subtruct_piplpimi_43prot09->Clone("h_Etot_subtruct_piplpimi_431prot09");
  h_Etot_subtruct_piplpimi_431prot09->Add(Etot_p431_bkgd09);

TH2F *h2_Erec_pperp_431psub=(TH2F*) h2_Erec_pperp_43psub->Clone("h2_Erec_pperp_431psub");
 h2_Erec_pperp_431psub->Add(h2_Erec_pperp_431p);

 TH1F *h_Erec_subtruct_piplpimi_431prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_43prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_431prot_fracfeed");
  h_Erec_subtruct_piplpimi_431prot_fracfeed->Add(Erec_43pto1p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_431prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_43prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_431prot_fracfeed");
  h_Etot_subtruct_piplpimi_431prot_fracfeed->Add(Etot_43pto1p_fracfeed);

//-----------------------------------undetected 4 to 2 proton subtraction  ---------------------------------------
  TH1F *h_Erec_subtruct_piplpimi_42prot=(TH1F*)  h_Erec_subtruct_piplpimi_431prot->Clone("h_Erec_subtruct_piplpimi_42prot");
  h_Erec_subtruct_piplpimi_42prot->Add(Erec_4pto2p);

  TH1F *h_Etot_subtruct_piplpimi_42prot=(TH1F*) h_Etot_subtruct_piplpimi_431prot->Clone("h_Etot_subtruct_piplpimi_42prot");
  h_Etot_subtruct_piplpimi_42prot->Add(Etot_4pto2p);

  TH1F *h_Etot_subtruct_piplpimi_42prot09=(TH1F*)  h_Etot_subtruct_piplpimi_431prot09->Clone("h_Etot_subtruct_piplpimi_42prot09");
  h_Etot_subtruct_piplpimi_42prot09->Add(Etot_p421_bkgd09);

TH2F *h2_Erec_pperp_42psub=(TH2F*) h2_Erec_pperp_431psub->Clone("h2_Erec_pperp_42psub");
 h2_Erec_pperp_42psub->Add(h2_Erec_pperp_421p);

  TH1F *h_Erec_subtruct_piplpimi_42prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_431prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_42prot_fracfeed");
  h_Erec_subtruct_piplpimi_42prot_fracfeed->Add(Erec_4pto2p_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_42prot_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_431prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_42prot_fracfeed");
  h_Etot_subtruct_piplpimi_42prot_fracfeed->Add(Etot_4pto2p_fracfeed);

 //-----------------------------------undetected 4 to 1 proton subtraction  ---------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_41prot=(TH1F*)  h_Erec_subtruct_piplpimi_42prot->Clone("h_Erec_subtruct_piplpimi_41prot");
  h_Erec_subtruct_piplpimi_41prot->Add(Erec_4pto1p,-1);

  TH1F *h_Etot_subtruct_piplpimi_41prot=(TH1F*)  h_Etot_subtruct_piplpimi_42prot->Clone("h_Etot_subtruct_piplpimi_41prot");
  h_Etot_subtruct_piplpimi_41prot->Add(Etot_4pto1p,-1);

  TH1F *h_Etot_subtruct_piplpimi_41prot09=(TH1F*)  h_Etot_subtruct_piplpimi_42prot09->Clone("h_Etot_subtruct_piplpimi_41prot09");
  h_Etot_subtruct_piplpimi_41prot09->Add(Etot_p41_bkgd09,-1);

TH2F *h2_Erec_pperp_41psub=(TH2F*) h2_Erec_pperp_42psub->Clone("h2_Erec_pperp_41psub");
 h2_Erec_pperp_41psub->Add(h2_Erec_pperp_41p,-1);

  TH1F *h_Erec_subtruct_piplpimi_41prot_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_42prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_41prot_fracfeed");
  h_Erec_subtruct_piplpimi_41prot_fracfeed->Add(Erec_4pto1p_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_41prot_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_42prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_41prot_fracfeed");
  h_Etot_subtruct_piplpimi_41prot_fracfeed->Add(Etot_4pto1p_fracfeed,-1);

//------------------------------------undetected 1p 2pi ->1 p1pi ------ --------------------------------------

  TH1F *h_Erec_subtruct_piplpimi_1p2pi=(TH1F*)  h_Erec_subtruct_piplpimi_41prot->Clone("h_Erec_subtruct_piplpimi_1p2pi");
  h_Erec_subtruct_piplpimi_1p2pi->Add(Erec_1p2pi);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi=(TH1F*)  h_Etot_subtruct_piplpimi_41prot->Clone("h_Etot_subtruct_piplpimi_1p2pi");
  h_Etot_subtruct_piplpimi_1p2pi->Add(Etot_1p2pi);

 TH1F *h_Etot_subtruct_piplpimi09_1p2pi_1p1pi=(TH1F*)  h_Etot_subtruct_piplpimi_41prot09->Clone("h_Etot_subtruct_piplpimi09_1p2pi_1p1pi");
  h_Etot_subtruct_piplpimi09_1p2pi_1p1pi->Add(Etot_bkgd09_1p2pi_1p1pi);

TH2F *h2_Erec_pperp_sub_1p2pi_1p1pi=(TH2F*) h2_Erec_pperp_41psub->Clone("h2_Erec_pperp_sub_1p2pi_1p1pi");
 h2_Erec_pperp_sub_1p2pi_1p1pi->Add(h2_Erec_pperp_1p2pi_1p1pi);

  TH1F *h_Erec_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_41prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_fracfeed");
  h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Add(Erec_1p2pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)  h_Etot_subtruct_piplpimi_41prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Add(Etot_1p2pi_fracfeed);

//------------------------------------undetected 1p 2pi-> 1p 0pi  ------ --------------------------------------

 TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi");
 h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Add(Erec_1p2pi_1p0pi,-1);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi");
  h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Add(Etot_1p2pi_1p0pi,-1);

 TH1F *h_Etot_subtruct_piplpimi09_1p2pi_1p0pi=(TH1F*)  h_Etot_subtruct_piplpimi09_1p2pi_1p1pi->Clone("h_Etot_subtruct_piplpimi09_1p2pi_1p0pi");
 h_Etot_subtruct_piplpimi09_1p2pi_1p0pi->Add(Etot_bkgd09_1p2pi_1p0pi,-1);

TH2F *h2_Erec_pperp_sub_1p2pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p1pi->Clone("h2_Erec_pperp_sub_1p2pi_1p0pi");
 h2_Erec_pperp_sub_1p2pi_1p0pi->Add(h2_Erec_pperp_1p2pi_1p0pi,-1);

 TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
 h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(Erec_1p2pi_1p0pi_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(Etot_1p2pi_1p0pi_fracfeed,-1);

//------------------------------------undetected 1p 3pi-> 1p 0pi  ------ --------------------------------------

 TH1F *h_Erec_subtruct_piplpimi_1p3pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Erec_subtruct_piplpimi_1p3pi");
 h_Erec_subtruct_piplpimi_1p3pi->Add(Erec_1p3pi);

  TH1F *h_Etot_subtruct_piplpimi_1p3pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Etot_subtruct_piplpimi_1p3pi");
  h_Etot_subtruct_piplpimi_1p3pi->Add(Etot_1p3pi);

 TH1F *h_Etot_subtruct_piplpimi09_1p3pi=(TH1F*)  h_Etot_subtruct_piplpimi09_1p2pi_1p0pi->Clone("h_Etot_subtruct_piplpimi09_1p3pi");
 h_Etot_subtruct_piplpimi09_1p3pi->Add(Etot_bkgd09_1p3pi);

TH2F *h2_Erec_pperp_sub_1p3pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p0pi->Clone("h2_Erec_pperp_sub_1p3pi");
 h2_Erec_pperp_sub_1p3pi->Add(h2_Erec_pperp_1p3pi);

 TH1F *h_Erec_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p3pi_fracfeed");
 h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Add(Erec_1p3pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p3pi_fracfeed");
  h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Add(Etot_1p3pi_fracfeed);


//------------------------------------undetected 2p 2pi ->1p 0pi  ------ --------------------------------------

 TH1F *h_Erec_subtruct_piplpimi_2p2pi=(TH1F*)  h_Erec_subtruct_piplpimi_1p3pi->Clone("h_Erec_subtruct_piplpimi_2p2pi");
 h_Erec_subtruct_piplpimi_2p2pi->Add(Erec_2p2pi);

  TH1F *h_Etot_subtruct_piplpimi_2p2pi=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi->Clone("h_Etot_subtruct_piplpimi_2p2pi");
  h_Etot_subtruct_piplpimi_2p2pi->Add(Etot_2p2pi);

 TH1F *h_Etot_subtruct_piplpimi09_2p2pi=(TH1F*)  h_Etot_subtruct_piplpimi09_1p3pi->Clone("h_Etot_subtruct_piplpimi09_2p2pi");
 h_Etot_subtruct_piplpimi09_2p2pi->Add(Etot_bkgd09_2p2pi);

TH2F *h2_Erec_pperp_sub_2p2pi=(TH2F*) h2_Erec_pperp_sub_1p3pi->Clone("h2_Erec_pperp_sub_2p2pi");
 h2_Erec_pperp_sub_2p2pi->Add(h2_Erec_pperp_2p2pi);

 TH1F *h_Erec_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p2pi_fracfeed");
 h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Add(Erec_2p2pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p2pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Add(Etot_2p2pi_fracfeed);


//------------------------------------undetected 3p 1pi ->1p 0pi  ------ --------------------------------------

 TH1F *h_Erec_subtruct_piplpimi_3p1pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p2pi->Clone("h_Erec_subtruct_piplpimi_3p1pi");
 h_Erec_subtruct_piplpimi_3p1pi->Add(Erec_3p1pi);

  TH1F *h_Etot_subtruct_piplpimi_3p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi->Clone("h_Etot_subtruct_piplpimi_3p1pi");
  h_Etot_subtruct_piplpimi_3p1pi->Add(Etot_3p1pi);

 TH1F *h_Etot_subtruct_piplpimi09_3p1pi=(TH1F*)  h_Etot_subtruct_piplpimi09_2p2pi->Clone("h_Etot_subtruct_piplpimi09_3p1pi");
 h_Etot_subtruct_piplpimi09_3p1pi->Add(Etot_bkgd09_3p1pi);

TH2F *h2_Erec_pperp_sub_3p1pi=(TH2F*) h2_Erec_pperp_sub_2p2pi->Clone("h2_Erec_pperp_sub_3p1pi");
 h2_Erec_pperp_sub_3p1pi->Add(h2_Erec_pperp_3p1pi);

 TH1F *h_Erec_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_3p1pi_fracfeed");
 h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Add(Erec_3p1pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_3p1pi_fracfeed");
  h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Add(Etot_3p1pi_fracfeed);



//------------------------------------undetected 2p 1pi ->2p 0pi  --------------------------------------------

 TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_3p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi");
 h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Add(Erec_2p1pi_2p0pi);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi");
  h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Add(Etot_2p1pi_2p0pi);

 TH1F *h_Etot_subtruct_piplpimi09_2p1pi_2p0pi=(TH1F*)  h_Etot_subtruct_piplpimi09_3p1pi->Clone("h_Etot_subtruct_piplpimi09_2p1pi_2p0pi");
 h_Etot_subtruct_piplpimi09_2p1pi_2p0pi->Add(Etot_bkgd09_2p1pi_2p0pi);

TH2F *h2_Erec_pperp_sub_2p1pi_2p0pi=(TH2F*) h2_Erec_pperp_sub_3p1pi->Clone("h2_Erec_pperp_sub_2p1pi_2p0pi");
 h2_Erec_pperp_sub_2p1pi_2p0pi->Add(h2_Erec_pperp_2p1pi_2p0pi);

 TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
 h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(Erec_2p1pi_2p0pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(Etot_2p1pi_2p0pi_fracfeed);


//------------------------------------undetected 2p 1pi ->1p 1pi  ------ --------------------------------------

 TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi");
 h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Add(Erec_2p1pi_1p1pi);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi");
  h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Add(Etot_2p1pi_1p1pi);

 TH1F *h_Etot_subtruct_piplpimi09_2p1pi_1p1pi=(TH1F*)  h_Etot_subtruct_piplpimi09_2p1pi_2p0pi->Clone("h_Etot_subtruct_piplpimi09_2p1pi_1p1pi");
 h_Etot_subtruct_piplpimi09_2p1pi_1p1pi->Add(Etot_bkgd09_2p1pi_1p1pi);

TH2F *h2_Erec_pperp_sub_2p1pi_1p1pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_2p0pi->Clone("h2_Erec_pperp_sub_2p1pi_1p1pi");
 h2_Erec_pperp_sub_2p1pi_1p1pi->Add(h2_Erec_pperp_2p1pi_1p1pi);

 TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
 h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(Erec_2p1pi_1p1pi_fracfeed);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(Etot_2p1pi_1p1pi_fracfeed);

//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

 TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi");
 h_Erec_subtruct_piplpimi_2p1pi_1p0pi->Add(Erec_2p1pi_1p0pi,-1);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi");
  h_Etot_subtruct_piplpimi_2p1pi_1p0pi->Add(Etot_2p1pi_1p0pi,-1);

 TH1F *h_Etot_subtruct_piplpimi09_2p1pi_1p0pi=(TH1F*)  h_Etot_subtruct_piplpimi09_2p1pi_1p1pi->Clone("h_Etot_subtruct_piplpimi09_2p1pi_1p0pi");
 h_Etot_subtruct_piplpimi09_2p1pi_1p0pi->Add(Etot_bkgd09_2p1pi_1p0pi,-1);

TH2F *h2_Erec_pperp_sub_2p1pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_1p1pi->Clone("h2_Erec_pperp_sub_2p1pi_1p0pi");
 h2_Erec_pperp_sub_2p1pi_1p0pi->Add(h2_Erec_pperp_2p1pi_1p0pi,-1);

 TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*)  h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
 h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(Erec_2p1pi_1p0pi_fracfeed,-1);

  TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
  h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(Etot_2p1pi_1p0pi_fracfeed,-1);








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


     delete[]  pperp_cut;
     delete[] Ecal_lowlim;
     delete[] Ecal_uplim;
     delete pipl_deltat_sig;delete pipl_deltat_mean;delete pimi_deltat_sig;delete pimi_deltat_mean;delete fsum_pimi;delete fsub_pimi;delete fsum_pipl;delete fsub_pipl;delete prot_deltat_sig;delete prot_deltat_mean;delete fsum_prot;delete fsub_prot;delete el_Epratio_sig;delete el_Epratio_mean;delete fsum_e;delete fsub_e;
delete up_lim1_ec;delete up_lim2_ec;delete up_lim3_ec;delete up_lim4_ec;delete up_lim5_ec;delete up_lim6_ec;delete low_lim1_ec;delete low_lim2_ec;delete low_lim3_ec;delete low_lim4_ec;delete low_lim5_ec;delete low_lim6_ec;
 delete  rightside_lim1_ec;delete rightside_lim2_ec;delete rightside_lim3_ec;delete rightside_lim4_ec; delete rightside_lim5_ec;delete rightside_lim6_ec;delete leftside_lim1_ec;delete leftside_lim2_ec; delete leftside_lim3_ec;delete leftside_lim4_ec;delete leftside_lim5_ec;delete leftside_lim6_ec;

}





double vz_corr(double phi,double theta)            //correction function for vertex , takes the arguments in Deg.
{
  //  return (0.2)*cos((phi+47.9)*TMath::DegToRad())/tan(theta*TMath::DegToRad()); // vertex correction function obtained for the empty runs 18393 and 18394, works fine for 3He runs at 2.261[GeV] beam energy
  return (-(vz_corr_func->GetParameter(1)))*cos((phi-(vz_corr_func->GetParameter(2)))*TMath::DegToRad())/tan(theta*TMath::DegToRad()); //vertex correction function for 4He runs at 2.261[GeV] beam energy obtained for the empty run18283

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
  return (ecuvw.X()>Constants::par_EcUVW[sector][0] && ecuvw.Y()<Constants::par_EcUVW[sector][1] && ecuvw.Z()<Constants::par_EcUVW[sector][2]);
}









void SetFiducialCutParameters(std::string beam_en){
// reads from a file the parameters of the fiducial cut functions
// Please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu(UNH)
  std::string fbeam_en = beam_en;


 if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){    //
   // reads FC parameters for 4.4GeV , e- and p fiducial cut parameters at 4GeV
   //

  std::ifstream param_file2(Form("/FiducialsCorrections/PFID_%s_%d.dat",fbeam_en.c_str(),fTorusCurrent));//reading the proton fiducial cut parameters at 4GeV
  std::ifstream param_file(Form("/FiducialsCorrections/FCP_%s_%d.dat",fbeam_en.c_str(),fTorusCurrent));


   //	std::ifstream param_file("./FCP_4461_2250.dat");
   int param_type, sector;
   double data[6];
   while ( (sector!=6 || param_type!=21))
     {
       param_file >> param_type;
       param_file >> sector >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];
       // Test the type of parameter and assign it to the proper data array
       //  std::cout << param_type << " " << sector << std::endl;
       switch (param_type)
	 {
	 case  0:
	   for(int k=0; k<2; k++) fgPar_4Gev_2250_Efid_t0_p[sector-1][k] = data[k];
	   break;
	 case  1:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_t1_p[sector-1][k] = data[k];
	   break;
	 case 10:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[sector-1][0][k] = data[k];
	   break;
	 case 11:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[sector-1][1][k] = data[k];
	   break;
	 case 20:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[sector-1][0][k] = data[k];
	   break;
	 case 21:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[sector-1][1][k] = data[k];
	   break;
	 default:
	   printf("Error in Efid parameter file!\nReceived parameter type %d, which is not found.\nAborting!\n\n\n",param_type);
	   break;
	 }
     } // Done reading in Fiducial Region Parameters
  throw "4.4GeV reading";
// ---
 for(int i = 0 ; i < 4 ; i++){
   for(int j = 0 ; j < 8 ; j++){
     param_file >> fgPar_4Gev_2250_Efid_Theta_S3[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 2 ; i++){
   for(int j = 0 ; j < 8 ; j++){
     param_file >> fgPar_4Gev_2250_Efid_Theta_S4[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 8 ; i++){
   for(int j = 0 ; j < 8 ; j++){
     param_file >> fgPar_4Gev_2250_Efid_Theta_S5[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 4 ; i++){
   for(int j = 0 ; j < 4 ; j++){
     param_file >> fgPar_4Gev_2250_Efid_Theta_S3_extra[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 2 ; i++){
   for(int j = 0 ; j < 4 ; j++){
     param_file >> fgPar_4Gev_2250_Efid_Theta_S4_extra[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 8 ; i++){
   for(int j = 0 ; j < 4 ; j++){
     param_file >> fgPar_4Gev_2250_Efid_Theta_S5_extra[i][j];
   }
 }
	param_file.close();



   for(int i = 0 ; i < 6 ; i++){
     for(int j = 0 ; j < 6 ; j++){
       param_file2 >> fgPar_4Gev_2250_Pfidft1l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidft1r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidft2l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidft2r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt1l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt1r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt2l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt2r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbl  [i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbr  [i][j];
     }
   }

for(int i = 0 ; i < 2 ; i++){//reading the proton bad TOF cuts at 4GeV
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS2[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS3[i][j];
  }
 }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS4[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS5[i][j];
  }
 }
for(int i = 0 ; i < 2 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][j];
  }
 }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][j];
  }
 }



 }




 else if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){

  std::ifstream param_file2(Form("./PFID_%s_%d.dat",fbeam_en.c_str(),fTorusCurrent));//reading the proton fiducial cut parameters at 4GeV
  std::ifstream param_file(Form("./FCP_%s_%d.dat",fbeam_en.c_str(),fTorusCurrent));
  std::ifstream param_file3(Form("./PIPFID_%s_%d.dat",fbeam_en.c_str(),fTorusCurrent));
  std::ifstream param_file4(Form("./PIMFID_%s_%d.dat",fbeam_en.c_str(),fTorusCurrent));
  //
   // reads FC parameters for 1.1GeV , e- fiducial cut parameters at 1GeV
   //

 if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {

for(Int_t sector=0;sector<6;sector++)
  {
    for(Int_t thetapar=0;thetapar<5;thetapar++)
      {
        for(Int_t mompar=0;mompar<6;mompar++)
          {
            param_file >> fgPar_1gev_1500_Efid[sector][thetapar][mompar];
          }
      }
  }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_1500_Efid_Theta_S3[i][j];
  }
 }
// ---
for(int i = 0 ; i < 2 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_1500_Efid_Theta_S4[i][j];
  }
 }
// ---
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_1500_Efid_Theta_S5[i][j];
  }
 }

    }


 if ( fTorusCurrent < 760 && fTorusCurrent > 740){

for(Int_t sector=0;sector<6;sector++)
  {
    for(Int_t thetapar=0;thetapar<5;thetapar++)
      {
        for(Int_t mompar=0;mompar<6;mompar++)
          {
            param_file >> fgPar_1gev_750_Efid[sector][thetapar][mompar];
          }
      }
  }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_750_Efid_Theta_S3[i][j];
  }
 }
// ---
for(int i = 0 ; i < 2 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_750_Efid_Theta_S4[i][j];
  }
 }
// ---
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_750_Efid_Theta_S5[i][j];
  }
 }
 }

	param_file.close();



  //
   // reads FC parameters for 1.1GeV , p fiducial cut parameters at 1GeV
   //

 if ( fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file2 >> fgPar_1gev_1500_Pfid[sector][phipar][mompar];
                  //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
                  //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
                }
            }
        }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS5[i][j];
        }
      }
    }
  if (fTorusCurrent< 760 && fTorusCurrent > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file2 >> fgPar_1gev_750_Pfid[sector][phipar][mompar];
                  //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
                  //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
                }
            }
        }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS5[i][j];
        }
      }
    }


  param_file2.close();



  //reads pimi fiducial cut parameters at 1GeV




  if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t thetapar=0;thetapar<5;thetapar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file4 >> fgPar_1gev_1500_Pimfid[sector][thetapar][mompar];
                }
            }
        }
      // ---
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S3[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S4[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S5[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4>> fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][j];
        }
      }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][j];
        }
      }
    }
  if (fTorusCurrent< 760 && fTorusCurrent > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t thetapar=0;thetapar<5;thetapar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file4 >> fgPar_1gev_750_Pimfid[sector][thetapar][mompar];
                }
            }
        }
      // ---
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S3[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S4[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S5[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S3_extra[i][j];
        }
      }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S4_extra[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S5_extra[i][j];
	  //	  std::cout << fgPar_1gev_750_Pimfid_Theta_S5_extra[i][j] << std::endl;
	}
      }
    }
	param_file4.close();


	//reads fiducial cut parameters for pi+

 if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file3 >> fgPar_1gev_1500_Piplfid[sector][phipar][mompar];
		  //  std::cout << "PFID " << fgPar_1gev_1500_Pfid[sector][phipar][mompar]  << std::endl;
                }
            }
        }


      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS5[i][j];
        }
      }


    }


  if ( fTorusCurrent < 760 && fTorusCurrent > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file3 >> fgPar_1gev_750_Piplfid[sector][phipar][mompar];
		  //  std::cout << "PFID " << fgPar_1gev_750_Pfid[sector][phipar][mompar]  << std::endl;
                }
            }
        }

   for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS5[i][j];
        }
      }




    }
  param_file3.close();



 }
else printf("There are no fiducial cut parameters to be read at %3.1f GeV!\n", en_beam[fbeam_en]);
 //	param_file2.close();



}



Bool_t GetEPhiLimits(std::string beam_en, Float_t momentum, Float_t theta, Int_t sector,Float_t *EPhiMin, Float_t *EPhiMax){
//Begin_Html
/*</pre>
 Information for electron fiducial cut,
    returns the minimum and maximum phi accepted for a given momentum, theta and sector
    momentum is in GeV/c, theta is in degrees, 0 <= sector <= 5
    EPhiMin and EPhiMax are in degrees
    Function returns False if inputs are out of bounds
    1.1 GeV not implemented yet
 tested against EFiducialCut to make sure the limits are identical
    2.2 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.1 < Ef < 2.261
             2 inconsistent events out of 10^6
    4.4 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.3 < Ef < 4.461
             0 inconsistent events out of 10^6
 Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
 For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
<pre>
*/
//End_Html
  std::string fbeam_en = beam_en;
  if (sector < 0 || sector > 5) return kFALSE;    // bad input

  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){// 4.4GeV fiducial cuts by protopop@jlab.org
    if ((theta < 15.) || (momentum < 0.9)) return kFALSE;         // out of range
    Float_t t0, t1, b[2], a[2];

    if (momentum > 3.7) momentum = 3.7; // don't extrapolate past the data


    // uncomment this if you want 100MeV energy bins
    //Enrgy = 0.100*int(Enrgy/0.100);


    // calculates parameters of cut functions for this energy
    t0 = fgPar_4Gev_2250_Efid_t0_p[sector][0]/pow(momentum, fgPar_4Gev_2250_Efid_t0_p[sector][1]);
    t1 = 0.; for(int k=0; k<6; k++) t1 += (fgPar_4Gev_2250_Efid_t1_p[sector][k]*pow(momentum, k));
    for(int l=0; l<2; l++){
      b[l] = 0.; for(int k=0; k<6; k++) b[l] += (fgPar_4Gev_2250_Efid_b_p[sector][l][k]*pow(momentum, k));
      a[l] = 0.; for(int k=0; k<6; k++) a[l] += (fgPar_4Gev_2250_Efid_a_p[sector][l][k]*pow(momentum, k));
    }



    // adjust upper limit according to hardware
    if(t1 < 45.) t1 = 45.;
    if(t0 < theta && theta < t1){

      *EPhiMin = 60.*sector - b[0]*(1. - 1/((theta - t0)/(b[0]/a[0]) + 1.));
      *EPhiMax = 60.*sector + b[1]*(1. - 1/((theta - t0)/(b[1]/a[1]) + 1.));
      // if(momentum<1.65 && momentum>1.60)cout<<sector<<"  "<<a[0]<<"    "<<a[1]<<"    "<<a[2]<<endl;
    }
    else {
      *EPhiMin = 60.*sector;
      *EPhiMax = 60.*sector;
    }


  }   // 4.4 GeV e2a
  else {
    return kFALSE;     // wrong beam energy/torus
  }
  return kTRUE;
}










Bool_t EFiducialCut(std::string beam_en, TVector3 momentum)
{

// Electron fiducial cut, return kTRUE if pass or kFALSE if not
  Bool_t status = kTRUE;
  std::string fbeam_en = beam_en;

 if(en_beam[fbeam_en]>1. &&  en_beam[fbeam_en]<2. && fTorusCurrent>740 && fTorusCurrent<1510) {

  Float_t phiMin, phiMax;
  Float_t mom = momentum.Mag();
  Float_t phi = momentum.Phi()*180./TMath::Pi();
  if(phi<-30.) phi += 360.;
  Float_t theta = momentum.Theta()*180./TMath::Pi();
  Int_t  sector = (Int_t)((phi+30.)/60.);
  if(sector < 0) sector = 0;
  if(sector > 5) sector = 5; // to match array index


    phi -= sector*60;
    Double_t elmom = (momentum.Mag())*1000;
    Double_t thetapars[5]={0,0,0,0,0};

    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t thetapar=0;thetapar<5;thetapar++) {
	if((fTorusCurrent>1490) && (fTorusCurrent<1510)) {
	  // 1500A torus current
	  thetapars[thetapar]+=fgPar_1gev_1500_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
	}
	if((fTorusCurrent>740) && (fTorusCurrent<760)) {
	  // 750A torus current
	  thetapars[thetapar]+=fgPar_1gev_750_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
	}
      }
    }

    Int_t uplow;
    Double_t thetacutoff;
    if(phi<=0) {
      uplow=1;
      thetacutoff=((phi*(thetapars[0]-(thetapars[1]/thetapars[2])))+
		   (double(uplow)*thetapars[2]*thetapars[0]))/(phi+(double(uplow)*thetapars[2]));
    }
    else {
      uplow=-1;
      thetacutoff=( (phi*(thetapars[0]-(thetapars[3]/thetapars[4]))) +
		   (double(uplow)*thetapars[4]*thetapars[0]))/(phi+(double(uplow)*thetapars[4]) );
    }

    status = (theta>thetacutoff) && (thetacutoff>=thetapars[0]) && (elmom>300) && (elmom<=1100);



		bool SCpdcut = true;
		if (SCpdcut && (fTorusCurrent>1490) && (fTorusCurrent<1510) ){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
		  if (status){
		    int tsector = sector + 1;
		    // sector 3 has two bad paddles
		    if (tsector == 3){
		      float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
		      for (int i=0; i<4; i++){
			badpar3[i] = 0;
			// calculate the parameters using pol7
			for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom + fgPar_1gev_1500_Efid_Theta_S3[i][d];}
		      }
		      for(int ipar=0;ipar<2;ipar++)
			status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
		    }
		    // sector 4 has one bad paddle
		    else if (tsector == 4){
		      float badpar4[2];     // 2 parameters to determine the position of the theta gap
		      for (int i=0; i<2; i++){
			badpar4[i] = 0;
			// calculate the parameters using pol7
			for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom + fgPar_1gev_1500_Efid_Theta_S4[i][d];}
		      }
		      status = !(theta>badpar4[0] && theta<badpar4[1]);
		    }
		    // sector 5 has four bad paddles
		    else if (tsector == 5){
		      Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
		      for (Int_t i=0; i<8; i++){
			badpar5[i] = 0;
			// calculate the parameters using pol7
			for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom + fgPar_1gev_1500_Efid_Theta_S5[i][d];}
		      }
		      if (mom<1.25) badpar5[0] = 23.4*1500/2250;
		      if (mom<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.

		      for(Int_t ip=0;ip<4;ip++)status = status && !(theta>badpar5[2*ip] && theta<badpar5[2*ip+1]);
		    }
		  }
		}


	if (SCpdcut && (fTorusCurrent>740) && (fTorusCurrent<760) ){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

		  if (status){

		    int tsector = sector + 1;
		    // sector 3 has two bad paddles
		    if (tsector == 3){
		      float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
		      for (int i=0; i<4; i++){
			badpar3[i] = 0;
			// calculate the parameters using pol7
			for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom + fgPar_1gev_750_Efid_Theta_S3[i][d];}
		      }
		      for(int ipar=0;ipar<2;ipar++)
			status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);

		    }
		    // sector 4 has one bad paddle
		    else if (tsector == 4){
		      float badpar4[2];     // 2 parameters to determine the position of the theta gap
		      for (int i=0; i<2; i++){
			badpar4[i] = 0;
			// calculate the parameters using pol7
			for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom + fgPar_1gev_750_Efid_Theta_S4[i][d];}
		      }
		      status = !(theta>badpar4[0] && theta<badpar4[1]);
		    }
		    // sector 5 has four bad paddles
		    else if (tsector == 5){
		      Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
		      for (Int_t i=0; i<8; i++){
			badpar5[i] = 0;
			// calculate the parameters using pol7
			for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom + fgPar_1gev_750_Efid_Theta_S5[i][d];}
		      }
		      if (mom<1.25) badpar5[0] = 23.4*750/2250;
		      if (mom<1.27) badpar5[1] = 24.0*750/2250; // some dummy constants. see fiducial cuts webpage.

		      for(Int_t ipar=0;ipar<4;ipar++){
			status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
		      }

		    }
		  }
		}


    return status;
  }





  if ( en_beam[fbeam_en]>2. &&  en_beam[fbeam_en]<3. && fTorusCurrent>2240 && fTorusCurrent<2260){
    Float_t phi=momentum.Phi()*180./TMath::Pi();
    if(phi<-30.) phi+=360.;
    Int_t sector = (Int_t)((phi+30.)/60.);
    if(sector<0)sector=0;
    if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180./TMath::Pi();
    Float_t mom = momentum.Mag();
    Float_t par[6];               // six parameters to determine the outline of Theta vs Phi
    for (Int_t i=0; i<6; i++){
      par[i] = 0;
      for (Int_t d=8; d>=0; d--){
	par[i] = par[i]*mom + Constants::fgPar_2GeV_2250_Efid[sector][i][d];
      }                          // calculate the parameters using pol8
    }
    if (phi < 0) {
      Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
      status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);
    }
    else {
      Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
      status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);
    }
    // by now, we have checked if the electron is within the outline of theta vs phi plot
    if (SCpdcut){  // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
      if (status){
	Int_t tsector = sector + 1;
	if (tsector == 3){               // sector 3 has two bad paddles
	  Float_t badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	  for (Int_t i=0; i<4; i++){
	    badpar3[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar3[i] = badpar3[i]*mom + Constants::fgPar_2GeV_2250_EfidTheta_S3[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  for(Int_t ipar=0;ipar<2;ipar++)
	    status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
	else if (tsector == 4){         // sector 4 has one bad paddle
	  Float_t badpar4[2];           // 2 parameters to determine the position of the theta gap
	  for (Int_t i=0; i<2; i++){
	    badpar4[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar4[i] = badpar4[i]*mom + Constants::fgPar_2GeV_2250_EfidTheta_S4[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  status = !(theta>badpar4[0] && theta<badpar4[1]);
	}
	else if (tsector == 5){         // sector 5 has four bad paddles
	  Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar5[i] = badpar5[i]*mom + Constants::fgPar_2GeV_2250_EfidTheta_S5[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  if (mom<1.25) badpar5[0] = 23.4;
	  if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	  for(Int_t ipar=0;ipar<4;ipar++)
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }
  }


  if ( en_beam[fbeam_en]>4. &&  en_beam[fbeam_en]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){


//Begin_Html
/*</pre>
  Electron fiducial cut, return kTRUE if the electron is in the fiducial volume
  modified 14 May 2001 lbw
  Now calls GetEPhiLimits for 2.2 and 4.4 GeV
  tested against EFiducialCut for both 2.2 (with and without bad scintillator cuts) and 4.4 GeV
  discrepancy less than 2 in 10^6 events
  Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
  For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
  Please refer to <a href="http://www.jlab.org/Hall-B/secure/e2/stevenmc/FiducialCuts/index.html">1.1 GeV fiducial cuts</a> -- Steven McLauchlan (GU).
<pre>
*/
//End_Html

  Float_t phiMin, phiMax;
  Float_t mom = momentum.Mag();
  Float_t phi = momentum.Phi()*180./TMath::Pi();
  if(phi<-30.) phi += 360.;
  Float_t theta = momentum.Theta()*180./TMath::Pi();
  Int_t  sector = (Int_t)((phi+30.)/60.);
  if(sector < 0) sector = 0;
  if(sector > 5) sector = 5; // to match array index
  // all the work is now done in GetEPhiLimits

  status = GetEPhiLimits(fbeam_en,mom, theta, sector, &phiMin, &phiMax);

  if (status) {
    status = status && (phi > phiMin) && (phi < phiMax);
  }



  if(mom <= 2.0)
    {
      bool SCpdcut = true;
      if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
	if (status){
	  int tsector = sector + 1;
	  // sector 3 has two bad paddles
	  if (tsector == 3){
	    float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	    for (int i=0; i<4; i++){
	      badpar3[i] = 0;
	      // calculate the parameters using pol7
	      for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom + fgPar_4Gev_2250_Efid_Theta_S3[i][d];}
	    }
	    for(int ipar=0;ipar<2;ipar++)
	      status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	  }
	  // sector 4 has one bad paddle
	  else if (tsector == 4){
	    float badpar4[2];     // 2 parameters to determine the position of the theta gap
	    for (int i=0; i<2; i++){
	      badpar4[i] = 0;
	      // calculate the parameters using pol7
	      for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom + fgPar_4Gev_2250_Efid_Theta_S4[i][d];}
	    }
	    status = !(theta>badpar4[0] && theta<badpar4[1]);
	  }
	  // sector 5 has four bad paddles
	  else if (tsector == 5){
	    Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	    for (Int_t i=0; i<8; i++){
	      badpar5[i] = 0;
	      // calculate the parameters using pol7
	      for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom + fgPar_4Gev_2250_Efid_Theta_S5[i][d];}
	    }
	    if (mom<1.25) badpar5[0] = 23.4;
	    if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	    for(Int_t ipar=0;ipar<4;ipar++)
	      status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	  }
	}
      }
      return (status && (phi < phiMax) && (phi>phiMin));
    }
  else{
    bool SCpdcut = true;
    if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
      if (status){
	int tsector = sector + 1;
	// sector 3 has two bad paddles
	if (tsector == 3){
	  float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	  for (int i=0; i<4; i++){
	    badpar3[i] = 0;
	    // calculate the parameters using 1/p
	    badpar3[i] = fgPar_4Gev_2250_Efid_Theta_S3_extra[i][0] + fgPar_4Gev_2250_Efid_Theta_S3_extra[i][1]/mom + fgPar_4Gev_2250_Efid_Theta_S3_extra[i][2]/(mom*mom) + fgPar_4Gev_2250_Efid_Theta_S3_extra[i][3]/(mom*mom*mom);
	  }
	  for(int ipar=0;ipar<2;ipar++)
	    status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
	// sector 4 has one bad paddle
	else if (tsector == 4){
	  float badpar4[2];     // 2 parameters to determine the position of the theta gap
	  for (int i=0; i<2; i++){
	    badpar4[i] = 0;
	    // calculate the parameters using 1/p
	    badpar4[i] = fgPar_4Gev_2250_Efid_Theta_S4_extra[i][0] + fgPar_4Gev_2250_Efid_Theta_S4_extra[i][1]/mom + fgPar_4Gev_2250_Efid_Theta_S4_extra[i][2]/(mom*mom) + fgPar_4Gev_2250_Efid_Theta_S4_extra[i][3]/(mom*mom*mom);
	  }
	  status = !(theta>badpar4[0] && theta<badpar4[1]);
	}
	// sector 5 has four bad paddles
	else if (tsector == 5){
	  Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    // calculate the parameters using 1/p
	    badpar5[i] = fgPar_4Gev_2250_Efid_Theta_S5_extra[i][0] + fgPar_4Gev_2250_Efid_Theta_S5_extra[i][1]/mom + fgPar_4Gev_2250_Efid_Theta_S5_extra[i][2]/(mom*mom) + fgPar_4Gev_2250_Efid_Theta_S5_extra[i][3]/(mom*mom*mom);
	  }
	  if (mom<1.25) badpar5[0] = 23.4;
	  if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	  for(Int_t ipar=0;ipar<4;ipar++)
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }
    return (status && (phi < phiMax) && (phi>phiMin));
  }





  }

  return status;
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







Bool_t PFiducialCut(std::string beam_en, TVector3 momentum){
  //Positive Hadron Fiducial Cut
  //Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).

   Bool_t status = kTRUE;
   std::string fbeam_en = beam_en;


  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){

	Float_t theta = momentum.Theta()*180/M_PI;
	Float_t phi   = momentum.Phi()  *180/M_PI;
	if(phi<-30) phi+=360;
	Int_t sector = Int_t ((phi+30)/60);
	if(sector<0) sector=0;
	if(sector>5) sector=5;
	phi -= sector*60;
	Float_t p = momentum.Mag();


if ( fTorusCurrent < 1510 && fTorusCurrent > 1490){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < .3)
      return false;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_1500_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = p;          // Momentum for bad sc paddles cuts
			if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
			if(tsector==2){      // sector 2 has one bad paddle
				Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
				for (Int_t i=0; i<2; i++){
					badpar2[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS2[i][d];
					}                // calculate the parameters using pol5
				}
				status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
			}
			else if(tsector==3){ // sector 3 has four bad paddles
				Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar3[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS3[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
				}
			}
			else if(tsector==4){ // sector 4 has two bad paddles
				Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<4; i++){
          if (i==0 || i==1)
            if (mom_scpd > .65)
              mom_scpd = .65;
					badpar4[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS4[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<2;ipar++){
					status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
				}
			}
			else if(tsector==5){ // sector 5 has four bad paddles
				Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar5[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS5[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }
  if (fTorusCurrent < 760 && fTorusCurrent > 740){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < .3)
      return false;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_750_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
			if (mom_scpd<0.15)mom_scpd=0.15; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
			if(tsector==2){      // sector 2 has one bad paddle
				Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
				for (Int_t i=0; i<2; i++){
					badpar2[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS2[i][d];
					}                // calculate the parameters using pol5
				}
				status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
			}
			else if(tsector==3){ // sector 3 has four bad paddles
				Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar3[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS3[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
				}
			}
			else if(tsector==4){ // sector 4 has two bad paddles
				Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<4; i++){
					badpar4[i] = 0;
          if (i==0 || i==1)
            if (mom_scpd > .65)
              mom_scpd = .65;
					for (Int_t d=5; d>=0; d--){
						badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS4[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<2;ipar++){
					status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
				}
			}
			else if(tsector==5){ // sector 5 has four bad paddles
				Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar5[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS5[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }

}




  if (en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3. && fTorusCurrent>2240 && fTorusCurrent<2260){
    Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();
    Float_t mom_for = p;              // momentum for forward constraints
    if (mom_for<0.3) mom_for = 0.3;   // momentum smaller than 300 MeV/c, use 300 MeV/c
    if (mom_for>1.6) mom_for = 1.6;   // momentum greater than 1.6 GeV/c, use 1.6 GeV/c
    Float_t mom_bak = p;              // momentum for backward constraints
    if (mom_bak<0.2) mom_bak = 0.2;   // momentum smaller than 200 MeV/c, use 200 MeV/c
    if (mom_bak>1.0) mom_bak = 1.0;   // momentum greater than 1.0 GeV/c, use 1.0 GeV/c
    Float_t theta0 = 8.5;
    Float_t phi_lower = -24.0;
    Float_t phi_upper = 24.0;
    Float_t par_for[4], par_bak[4];
    for (Int_t i=0; i<4; i++){
      par_for[i] = 0; par_bak[i] = 0;
      for (Int_t d=6; d>=0; d--){
	par_for[i] = par_for[i]*mom_for + Constants::fgPar_2GeV_2250_Pfid_For[sector][i][d];
	par_bak[i] = par_bak[i]*mom_bak + Constants::fgPar_2GeV_2250_Pfid_Bak[sector][i][d];
      }
    }
    if (phi < 0) {
      Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
      status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
    }
    else {
      Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
      status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);
    }                     // now the forward constrains are checked
    if ( status ) {       // now check the backward constrains
      if(theta>par_bak[0]) status = kFALSE;
      else if(theta>par_bak[1]) status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]);
    }

    if(status && SCpdcut){ // cut bad scintillator paddles

      Int_t tsector = sector + 1;
      Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
      if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if(tsector==2){      // sector 2 has one bad paddle
	Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
	for (Int_t i=0; i<2; i++){
	  badpar2[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar2[i] = badpar2[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS2[i][d];
	  }                // calculate the parameters using pol5
	}
	status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
      }
      else if(tsector==3){ // sector 3 has four bad paddles
	Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar3[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar3[i] = badpar3[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS3[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
      }
      else if(tsector==4){ // sector 4 has two bad paddles
	Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<4; i++){
	  badpar4[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar4[i] = badpar4[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS4[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<2;ipar++){
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	}
      }
      else if(tsector==5){ // sector 5 has four bad paddles
	Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar5[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar5[i] = badpar5[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS5[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }
  }



  if (en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){//4 GeV Fiducial Cut Rustam Niyazov

    Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = Int_t ((phi+30)/60); if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();

    Float_t parfidl[3];for(Int_t i=0; i<3; i++){parfidl[i]=0;}
    Float_t parfidr[3];for(Int_t i=0; i<3; i++){parfidr[i]=0;}
    Float_t parfidbl[2];for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
    Float_t parfidbr[2];for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
    Float_t cphil=0;Float_t cphir=0;
    Float_t phi45l=0; Float_t phi45r=0;
    Float_t phi60l=0; Float_t phi60r=0;
    Float_t theta_min=11;

    bool Forward=kFALSE; //defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
    Int_t thetab=45; //this variable defines the edge point for Forward<->Backward regions
    Float_t p1=0.575; //last bin momentum for region p<0.6 GeV/c
    Float_t theta_max=140;
    if(p<0.2)p=0.2; //momentum less than 0.2 GeV/c, use 0.2 GeV/c
    if(p>4.4)p=4.4; //momentum greater than 4.4 GeV/c, use 4.4 GeV/c

    //get parametrized values of theta_max for p<0.6 GeV/c region
    if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
    //get parametrized values of theta_max for p>0.6 GeV/c region
    else{theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p1;}

    //Get the momentum dependent parameters for Forward Region (theta <45 deg)
    Forward=kTRUE;
    if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
      //parameters for hyperbolic function
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;
      }
    }
    else{//forward2 defines  regions of momenta and p>0.6 GeV/c
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;
      }
    }
    phi45l=parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
    phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
    if(theta>thetab){//backward region defined by theta >45 deg.
      if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
      if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

      //Get the momentum dependent parameters for Backward Region

      Forward=kFALSE;
      if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;
        }
      }
      else{//backward2 defines  regions of momenta p>0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;
        }
      }
    }

    if(Forward){//Forward region
      if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.
        cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
        cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
    }
    else{//Backward region
      phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
      phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

      if(theta<60){
        cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
        cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
      }
      Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0;
      //dr and er are theta_flat and phi_edge parameters for phi>0;
      dl=parfidbl[0];el=parfidbl[1];
      dr=parfidbr[0];er=parfidbr[1];

      if(theta>45&&theta<60){ //BackwardA region
        //try to match parametrized values from Forward region to Backward region parameters
        if(cphil>phi45l)cphil=phi45l;
        if(cphir<phi45r)cphir=phi45r;
      }
      //BackwardB region & phi<0
      else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant
      else if(theta>dl&&theta<=theta_max){
        cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line
      else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
      //BackwardB region & phi>0
      if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant
      else if(theta>dr&&theta<=theta_max){
        cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line
      else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
    }//Backward Region

    if(phi<0) status=(phi>cphil); //check the constrains
    else if(phi>=0) {status=(phi<cphir);
    }

    if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min

    if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

   //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for
   //some range of momentum, where edge does not look good.
    bool s1s4=(theta<11.7&&(sector==0||sector==3));
    bool s5=(theta<12.2&&sector==4);
    bool s6=(theta<11.4&&sector==5);
    if(p>=0.6&&p<1.5&&(s1s4||s5||s6)) status=kFALSE;



bool SCpdcut = true;
    if(status && SCpdcut){ // cut bad scintillator paddles
      if(p < 1.0){
        Int_t tsector = sector + 1;
        Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
        if (mom_scpd<0.3)mom_scpd=0.3; // momentum smaller than 200 MeV/c, use 200 MeV/c
        if(tsector==2){      // sector 2 has one bad paddle
          Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
          for (Int_t i=0; i<2; i++){
            badpar2[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar2[i] = badpar2[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2[i][d];
            }                // calculate the parameters using pol5
          }
          status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
        }
        else if(tsector==3){ // sector 3 has four bad paddles
          Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar3[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar3[i] = badpar3[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
        }
        else if(tsector==4){ // sector 4 has two bad paddles
          Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<4; i++){
            badpar4[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar4[i] = badpar4[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<2;ipar++){
            status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
        }
        else if(tsector==5){ // sector 5 has four bad paddles
          Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar5[i] = badpar5[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
          }
        }
      }
      else{
        int tsector = sector + 1;
        double mom_scpd =p;
        // sector 2 has one bad paddles
        if (tsector == 2){
          float badpar2[2];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<2; i++){
            badpar2[i] = 0;
            // calculate the parameters using 1/p
            badpar2[i] = fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<1;ipar++)
	    status = status && !(theta>badpar2[2*ipar] && theta<badpar2[2*ipar+1]);
	  // status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
        if (tsector == 3){
          float badpar3[8];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<8; i++){
            badpar3[i] = 0;
            // calculate the parameters using 1/p
            badpar3[i] = fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
        }
        // sector 4 has two bad paddle
        else if (tsector == 4){
          float badpar4[4];     // 2 parameters to determine the position of the theta gap
          for (int i=0; i<4; i++){
            badpar4[i] = 0;
            // calculate the parameters using 1/p
            badpar4[i] = fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<2;ipar++)
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
        // sector 5 has four bad paddles
        else if (tsector == 5){
          Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            // calculate the parameters using 1/p
            badpar5[i] = fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(Int_t ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
        }
      }
    }



  }
  return status;
}






Bool_t PiplFiducialCut(std::string beam_en, TVector3 momentum, Float_t *philow, Float_t *phiup){
  //Positive Hadron Fiducial Cut
  //Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
  std::string fbeam_en = beam_en;
  Bool_t status = kTRUE;

 if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){

	Float_t theta = momentum.Theta()*180/M_PI;
	Float_t phi   = momentum.Phi()  *180/M_PI;
	if(phi<-30) phi+=360;
	Int_t sector = Int_t ((phi+30)/60);
	if(sector<0) sector=0;
	if(sector>5) sector=5;
	phi -= sector*60;
	Float_t p = momentum.Mag();


if (fTorusCurrent < 1510 && fTorusCurrent > 1490){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < 0.15)p=0.15;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_1500_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = p;          // Momentum for bad sc paddles cuts
			if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
			if(tsector==2){      // sector 2 has one bad paddle
				Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
				for (Int_t i=0; i<2; i++){
					badpar2[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS2[i][d];
					}                // calculate the parameters using pol5
				}
				status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
			}
			else if(tsector==3){ // sector 3 has four bad paddles
				Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar3[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS3[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
				}
			}
			else if(tsector==4){ // sector 4 has two bad paddles
				Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<4; i++){
          if (i==0 || i==1)
            if (mom_scpd > .65)
              mom_scpd = .65;
					badpar4[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS4[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<2;ipar++){
					status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
				}
			}
			else if(tsector==5){ // sector 5 has four bad paddles
				Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar5[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS5[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }
  if (fTorusCurrent < 760 && fTorusCurrent > 740){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < 0.15)p=0.15;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_750_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
			if (mom_scpd<0.15)mom_scpd=0.15; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
			if(tsector==2){      // sector 2 has one bad paddle
				Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
				for (Int_t i=0; i<2; i++){
					badpar2[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS2[i][d];
					}                // calculate the parameters using pol5
				}
				status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
			}
			else if(tsector==3){ // sector 3 has four bad paddles
				Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar3[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS3[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
				}
			}
			else if(tsector==4){ // sector 4 has two bad paddles
				Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<4; i++){
					badpar4[i] = 0;
          if (i==0 || i==1)
            if (mom_scpd > .65)
              mom_scpd = .65;
					for (Int_t d=5; d>=0; d--){
						badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS4[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<2;ipar++){
					status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
				}
			}
			else if(tsector==5){ // sector 5 has four bad paddles
				Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar5[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_750_Pfid_ScpdS5[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }

}




  if (en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3. && fTorusCurrent>2240 && fTorusCurrent<2260){

    Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();
    Float_t mom_for = p;              // momentum for forward constraints
    if (mom_for<0.3) mom_for = 0.3;   // momentum smaller than 300 MeV/c, use 300 MeV/c
    if (mom_for>1.6) mom_for = 1.6;   // momentum greater than 1.6 GeV/c, use 1.6 GeV/c
    Float_t mom_bak = p;              // momentum for backward constraints
    if (mom_bak<0.2) mom_bak = 0.2;   // momentum smaller than 200 MeV/c, use 200 MeV/c
    if (mom_bak>1.0) mom_bak = 1.0;   // momentum greater than 1.0 GeV/c, use 1.0 GeV/c
    Float_t theta0 = 8.5;
    Float_t phi_lower = -24.0;
    Float_t phi_upper = 24.0;
    Float_t phimin, phimax;
    Float_t par_for[4], par_bak[4];
    for (Int_t i=0; i<4; i++){
      par_for[i] = 0; par_bak[i] = 0;
      for (Int_t d=6; d>=0; d--){
	par_for[i] = par_for[i]*mom_for + Constants::fgPar_2GeV_2250_Pfid_For[sector][i][d];
	par_bak[i] = par_bak[i]*mom_bak + Constants::fgPar_2GeV_2250_Pfid_Bak[sector][i][d];
      }
    }
    if (phi < 0) {
      Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
       phimin = par_for[1]/((theta-theta0)+par_for[1]/par_for[0])-par_for[0];
       phimax = par_for[0]-par_for[1]/((theta-theta0)+par_for[1]/par_for[0]);
       *philow = phimin;
       *phiup = phimax;
      status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
    }
    else {
      Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
      phimin = par_for[3]/(theta-theta0+par_for[3]/par_for[2])-par_for[2];
      phimax = par_for[2]-par_for[3]/(theta-theta0+par_for[3]/par_for[2]);
      *phiup = phimax;
      *philow = phimin;
      status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);
    }                     // now the forward constrains are checked
    if ( status ) {       // now check the backward constrains
      if(theta>par_bak[0]) status = kFALSE;
      else if(theta>par_bak[1]) status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]);
    }

    if(status && SCpdcut){ // cut bad scintillator paddles

      Int_t tsector = sector + 1;
      Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
      if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if(tsector==2){      // sector 2 has one bad paddle
	Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
		for (Int_t i=0; i<2; i++){
	  badpar2[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar2[i] = badpar2[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS2[i][d];
	  }                // calculate the parameters using pol5
	}
	status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
      }
      else if(tsector==3){ // sector 3 has four bad paddles
	Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar3[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar3[i] = badpar3[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS3[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
      }
      else if(tsector==4){ // sector 4 has two bad paddles
	Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<4; i++){
	  badpar4[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar4[i] = badpar4[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS4[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<2;ipar++){
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	}
      }
      else if(tsector==5){ // sector 5 has four bad paddles
	Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar5[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar5[i] = badpar5[i]*mom_scpd + Constants::fgPar_2GeV_2250_Pfid_ScpdS5[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }

  }


  if (en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){//4 GeV Fiducial Cut Rustam Niyazov

   Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = Int_t ((phi+30)/60); if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();

    Float_t parfidl[3];for(Int_t i=0; i<3; i++){parfidl[i]=0;}
    Float_t parfidr[3];for(Int_t i=0; i<3; i++){parfidr[i]=0;}
    Float_t parfidbl[2];for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
    Float_t parfidbr[2];for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
    Float_t cphil=0;Float_t cphir=0;
    Float_t phi45l=0; Float_t phi45r=0;
    Float_t phi60l=0; Float_t phi60r=0;
    Float_t theta_min=11;

    bool Forward=kFALSE; //defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
    Int_t thetab=45; //this variable defines the edge point for Forward<->Backward regions
    Float_t p1=0.575; //last bin momentum for region p<0.6 GeV/c
    Float_t theta_max=140;
    if(p<0.2)p=0.2; //momentum less than 0.2 GeV/c, use 0.2 GeV/c
    if(p>4.4)p=4.4; //momentum greater than 4.4 GeV/c, use 4.4 GeV/c

    //get parametrized values of theta_max for p<0.6 GeV/c region
    if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
    //get parametrized values of theta_max for p>0.6 GeV/c region
    else{theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p1;}

    //Get the momentum dependent parameters for Forward Region (theta <45 deg)
    Forward=kTRUE;
    if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
      //parameters for hyperbolic function
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;
      }
    }
    else{//forward2 defines  regions of momenta and p>0.6 GeV/c
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;
      }
    }
    phi45l=parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
    phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
    if(theta>thetab){//backward region defined by theta >45 deg.
      if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
      if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

      //Get the momentum dependent parameters for Backward Region

      Forward=kFALSE;
      if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;
        }
      }
      else{//backward2 defines  regions of momenta p>0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;
        }
      }
    }

    if(Forward){//Forward region
      if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.
        cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
        cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
    }
    else{//Backward region
      phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
      phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

      if(theta<60){
        cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
        cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
      }
      Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0;
      //dr and er are theta_flat and phi_edge parameters for phi>0;
      dl=parfidbl[0];el=parfidbl[1];
      dr=parfidbr[0];er=parfidbr[1];

      if(theta>45&&theta<60){ //BackwardA region
        //try to match parametrized values from Forward region to Backward region parameters
        if(cphil>phi45l)cphil=phi45l;
        if(cphir<phi45r)cphir=phi45r;
      }
      //BackwardB region & phi<0
      else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant
      else if(theta>dl&&theta<=theta_max){
        cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line
      else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
      //BackwardB region & phi>0
      if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant
      else if(theta>dr&&theta<=theta_max){
        cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line
      else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
    }//Backward Region


    if(phi<0) status=(phi>cphil); //check the constrains
    else if(phi>=0) {status=(phi<cphir);
  }

    if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min

    if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

   //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for
   //some range of momentum, where edge does not look good.
    bool s1s4=(theta<11.7&&(sector==0||sector==3));
    bool s5=(theta<12.2&&sector==4);
    bool s6=(theta<11.4&&sector==5);
    if(p>=0.6&&p<1.5&&(s1s4||s5||s6)) status=kFALSE;

    *philow = cphil;
    *phiup = cphir;



bool SCpdcut = true;
    if(status && SCpdcut){ // cut bad scintillator paddles
      if(p < 1.0){
        Int_t tsector = sector + 1;
        Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
        if (mom_scpd<0.3)mom_scpd=0.3; // momentum smaller than 200 MeV/c, use 200 MeV/c
        if(tsector==2){      // sector 2 has one bad paddle
          Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
          for (Int_t i=0; i<2; i++){
            badpar2[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar2[i] = badpar2[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2[i][d];
            }                // calculate the parameters using pol5
          }
          status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
        }
        else if(tsector==3){ // sector 3 has four bad paddles
          Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar3[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar3[i] = badpar3[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
        }
        else if(tsector==4){ // sector 4 has two bad paddles
          Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<4; i++){
            badpar4[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar4[i] = badpar4[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<2;ipar++){
            status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
        }
        else if(tsector==5){ // sector 5 has four bad paddles
          Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar5[i] = badpar5[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
          }
        }
      }
      else{
        int tsector = sector + 1;
        double mom_scpd =p;
        // sector 2 has one bad paddles
        if (tsector == 2){
          float badpar2[2];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<2; i++){
            badpar2[i] = 0;
            // calculate the parameters using 1/p
            badpar2[i] = fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<1;ipar++)
            status = status && !(theta>badpar2[2*ipar] && theta<badpar2[2*ipar+1]);
        }
        if (tsector == 3){
          float badpar3[8];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<8; i++){
            badpar3[i] = 0;
            // calculate the parameters using 1/p
            badpar3[i] = fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
        }
        // sector 4 has two bad paddle
        else if (tsector == 4){
          float badpar4[4];     // 2 parameters to determine the position of the theta gap
          for (int i=0; i<4; i++){
            badpar4[i] = 0;
            // calculate the parameters using 1/p
            badpar4[i] = fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<2;ipar++)
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
        // sector 5 has four bad paddles
        else if (tsector == 5){
          Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            // calculate the parameters using 1/p
            badpar5[i] = fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(Int_t ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
        }
      }
    }



  }


  return status;
}








//using two GeV pimi fiducial cuts for both 2 and 4 GeV analysis

Bool_t PimiFiducialCut(std::string beam_en, TVector3 momentum, Float_t *pimi_philow, Float_t *pimi_phiup){
  // Electron fiducial cut, return kTRUE if pass or kFALSE if not
  std::string fbeam_en = beam_en;
  Bool_t status = kTRUE;

  if(en_beam[fbeam_en]>1. &&  en_beam[fbeam_en]<2.) {


  TVector3 mom = momentum;
  double phi = mom.Phi();
  if (phi < -M_PI/6.) phi+= 2.*M_PI;
  int sector = (phi+M_PI/6.)/(M_PI/3.);
  sector = sector%6;
  double phi_deg = phi * 180./M_PI;
  phi_deg -= sector*60;

  double theta = mom.Theta();
  double theta_deg = theta * 180./M_PI;
  double mom_e = mom.Mag();


  if( fTorusCurrent>1490 && fTorusCurrent<1510){

    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    if (mom_e < .15)
      mom_e = .15;
    if (mom_e > 1.1)
      mom_e = 1.1;

    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_1500_Pimfid[sector][phipar][mompar]*pow(mom_e,mompar);
        //std::cout << mom_e << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi_deg<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta_deg-phipars[4])/phipars[3]+1.)));
      status = ((phi_deg>phicutoff) && (theta_deg>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta_deg-phipars[4])/phipars[2]+1.)));
      status = ((phi_deg<phicutoff) && (theta_deg>phipars[4]));
    }
    if (mom_e >= .3)
      {
        bool SCpdcut = true;
        if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
          if (status){
            int tsector = sector + 1;
            // sector 3 has two bad paddles
            if (tsector == 3){
              float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
              for (int i=0; i<4; i++){
                badpar3[i] = 0;
                // calculate the parameters using pol7
                for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_e + fgPar_1gev_1500_Pimfid_Theta_S3[i][d];}
              }
              for(int ipar=0;ipar<2;ipar++)
                status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
            }
            // sector 4 has one bad paddle
            else if (tsector == 4){
              float badpar4[2];     // 2 parameters to determine the position of the theta gap
              for (int i=0; i<2; i++){
                badpar4[i] = 0;
                // calculate the parameters using pol7
                for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_e + fgPar_1gev_1500_Pimfid_Theta_S4[i][d];}
              }
              status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
            }
            // sector 5 has four bad paddles
            else if (tsector == 5){
              Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
              for (Int_t i=0; i<8; i++){
                badpar5[i] = 0;
                // calculate the parameters using pol7
                for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_e + fgPar_1gev_1500_Pimfid_Theta_S5[i][d];}
              }
              if (mom_e<1.25) badpar5[0] = 23.4*1500/2250;
              if (mom_e<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.
              for(Int_t ipar=0;ipar<4;ipar++)
                status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
            }
          }
        }
        return status;
      }
    else{
      bool SCpdcut = true;
      if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
        if (status){
          int tsector = sector + 1;
          // sector 3 has two bad paddles
          if (tsector == 3){
            float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
            for (int i=0; i<4; i++){
              badpar3[i] = 0;
              // calculate the parameters using 1/p
              badpar3[i] = fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][1]/mom_e + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            for(int ipar=0;ipar<2;ipar++)
              status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
          }
          // sector 4 has one bad paddle
          else if (tsector == 4){
            float badpar4[2];     // 2 parameters to determine the position of the theta gap
            for (int i=0; i<2; i++){
              badpar4[i] = 0;
              // calculate the parameters using 1/p
              badpar4[i] = fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][1]/mom_e + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
          }
          // sector 5 has four bad paddles
          else if (tsector == 5){
            Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
            for (Int_t i=0; i<8; i++){
              badpar5[i] = 0;
              // calculate the parameters using 1/p
              badpar5[i] = fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][1]/mom_e + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            if (mom_e<1.25) badpar5[0] = 23.4*1500/2250;
            if (mom_e<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.
            for(Int_t ipar=0;ipar<4;ipar++)
              status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
          }
        }
      }
      return (status);
    }
  }
  if ( fTorusCurrent>740 && fTorusCurrent<760){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    if (mom_e < .15)
      mom_e = .15;
    if (mom_e > 1.1)
      mom_e = 1.1;
    //    std::cout << mom_e << std::endl;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_750_Pimfid[sector][phipar][mompar]*pow(mom_e,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi_deg<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta_deg-phipars[4])/phipars[3]+1.)));
      status = ((phi_deg>phicutoff) && (theta_deg>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta_deg-phipars[4])/phipars[2]+1.)));
      status = ((phi_deg<phicutoff) && (theta_deg>phipars[4]));
    }
    if (mom_e >= .3)
      {
        bool SCpdcut = true;
        if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
          if (status){
            int tsector = sector + 1;
            // sector 3 has two bad paddles
            if (tsector == 3){
              float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
              for (int i=0; i<4; i++){
                badpar3[i] = 0;
                // calculate the parameters using pol7
                for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_e + fgPar_1gev_750_Pimfid_Theta_S3[i][d];}
              }
              for(int ipar=0;ipar<2;ipar++)
                status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
            }
            // sector 4 has one bad paddle
            else if (tsector == 4){
              float badpar4[2];     // 2 parameters to determine the position of the theta gap
              for (int i=0; i<2; i++){
                badpar4[i] = 0;
                // calculate the parameters using pol7
                for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_e + fgPar_1gev_750_Pimfid_Theta_S4[i][d];}
              }
              status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
            }
            // sector 5 has four bad paddles
            else if (tsector == 5){
              Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
              for (Int_t i=0; i<8; i++){
                badpar5[i] = 0;
                // calculate the parameters using pol7
                for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_e + fgPar_1gev_750_Pimfid_Theta_S5[i][d];}
              }
              if (mom_e<1.25) badpar5[0] = 23.4*750/2250;
              if (mom_e<1.27) badpar5[1] = 24.0*750/2250; // some dummy constants. see fiducial cuts webpage.
              for(Int_t ipar=0;ipar<4;ipar++)
                status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
            }
          }
        }
        return status;
      }
    else{
      bool SCpdcut = true;
      if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
        if (status){
          int tsector = sector + 1;
          // sector 3 has two bad paddles
          if (tsector == 3){
            float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
            for (int i=0; i<4; i++){
              badpar3[i] = 0;
              // calculate the parameters using 1/p
              badpar3[i] = fgPar_1gev_750_Pimfid_Theta_S3_extra[i][0] + fgPar_1gev_750_Pimfid_Theta_S3_extra[i][1]/mom_e + fgPar_1gev_750_Pimfid_Theta_S3_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_750_Pimfid_Theta_S3_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            for(int ipar=0;ipar<2;ipar++)
              status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
          }
          // sector 4 has one bad paddle
          else if (tsector == 4){
            float badpar4[2];     // 2 parameters to determine the position of the theta gap
            for (int i=0; i<2; i++){
              badpar4[i] = 0;
              // calculate the parameters using 1/p
              badpar4[i] = fgPar_1gev_750_Pimfid_Theta_S4_extra[i][0] + fgPar_1gev_750_Pimfid_Theta_S4_extra[i][1]/mom_e + fgPar_1gev_750_Pimfid_Theta_S4_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_750_Pimfid_Theta_S4_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
          }
          // sector 5 has four bad paddles
          else if (tsector == 5){
            Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
            for (Int_t i=0; i<8; i++){
              badpar5[i] = 0;
              // calculate the parameters using 1/p
              badpar5[i] = fgPar_1gev_750_Pimfid_Theta_S5_extra[i][0] + fgPar_1gev_750_Pimfid_Theta_S5_extra[i][1]/mom_e + fgPar_1gev_750_Pimfid_Theta_S5_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_750_Pimfid_Theta_S5_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            if (mom_e<1.25) badpar5[0] = 23.4*750/2250;
            if (mom_e<1.27) badpar5[1] = 24.0*750/2250; // some dummy constants. see fiducial cuts webpage.
            for(Int_t ipar=0;ipar<4;ipar++)
              status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
          }
        }
      }
      return (status);
    }
  }
 }


   if ( en_beam[fbeam_en]>2. &&  en_beam[fbeam_en]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){

 Float_t phi=momentum.Phi()*180./TMath::Pi();
    if(phi<-30.) phi+=360.;
    Int_t sector = (Int_t)((phi+30.)/60.);
    if(sector<0)sector=0;
    if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180./TMath::Pi();
    if (theta>45)theta=45;   //to extrapolate the cut to higher theta for pimi
   Float_t mom = momentum.Mag(), phimin, phimax;
   if(mom>2.2)mom=2.2; //to extrapolate the cut to higher momenta for pimi
    Float_t par[6];               // six parameters to determine the outline of Theta vs Phi
    for (Int_t i=0; i<6; i++){
      par[i] = 0;
      for (Int_t d=8; d>=0; d--){
	par[i] = par[i]*mom + Constants::fgPar_2GeV_2250_Efid[sector][i][d];
      }                          // calculate the parameters using pol8
    }
    if (phi < 0) {
      Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
      phimin =  par[3]/((theta-par[0])+par[3]/par[2])-par[2];
      phimax =  par[2]-par[3]/((theta-par[0])+par[3]/par[2]);
      *pimi_philow = phimin;
      *pimi_phiup = phimax;
      status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);
    }
    else {
      Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
      phimin =  par[5]/((theta-par[0])+par[5]/par[4])-par[4];
      phimax =  par[4]-par[5]/((theta-par[0])+par[5]/par[4]);
      *pimi_philow = phimin;
      *pimi_phiup = phimax;
      status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);
    }
    // by now, we have checked if the electron is within the outline of theta vs phi plot
    if (SCpdcut){  // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
      if (status){
	Int_t tsector = sector + 1;
	if (tsector == 3){               // sector 3 has two bad paddles
	  Float_t badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	  for (Int_t i=0; i<4; i++){
	    badpar3[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar3[i] = badpar3[i]*mom + Constants::fgPar_2GeV_2250_EfidTheta_S3[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  for(Int_t ipar=0;ipar<2;ipar++)
	    status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
	else if (tsector == 4){         // sector 4 has one bad paddle
	  Float_t badpar4[2];           // 2 parameters to determine the position of the theta gap
	  for (Int_t i=0; i<2; i++){
	    badpar4[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar4[i] = badpar4[i]*mom + Constants::fgPar_2GeV_2250_EfidTheta_S4[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  status = !(theta>badpar4[0] && theta<badpar4[1]);
	}
	else if (tsector == 5){         // sector 5 has four bad paddles
	  Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar5[i] = badpar5[i]*mom + Constants::fgPar_2GeV_2250_EfidTheta_S5[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  if (mom<1.25) badpar5[0] = 23.4;
	  if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	  for(Int_t ipar=0;ipar<4;ipar++)
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }


  }
  return status;
}





bool Phot_fid(TVector3 V3_phot){

  bool status=true;
  //  double costheta=V3_phot.Pz();
  double costheta=V3_phot.CosTheta();
  double theta_deg=TMath::ACos(V3_phot.Pz())*TMath::RadToDeg();
  double phi_deg=TMath::ATan2(V3_phot.Py(),V3_phot.Px())*TMath::RadToDeg()+30;
  if(phi_deg<0)phi_deg=phi_deg+360;
  bool hot_spot=(phi_deg>185 && phi_deg<191 && costheta<0.71 && costheta>0.67) || (phi_deg>221 && phi_deg<236 && costheta<0.73 && costheta>0.67); // used to kill the two hot spots in sector 4


  if((costheta>low_lim1_ec->Eval(phi_deg) && costheta<up_lim1_ec->Eval(phi_deg) && costheta<rightside_lim1_ec->Eval(phi_deg) && costheta<leftside_lim1_ec->Eval(phi_deg))  ||
     (costheta>low_lim2_ec->Eval(phi_deg) && costheta<up_lim2_ec->Eval(phi_deg) && costheta<rightside_lim2_ec->Eval(phi_deg) && costheta<leftside_lim2_ec->Eval(phi_deg))       ||
     (costheta>low_lim3_ec->Eval(phi_deg) && costheta<up_lim3_ec->Eval(phi_deg) && costheta<rightside_lim3_ec->Eval(phi_deg) && costheta<leftside_lim3_ec->Eval(phi_deg))       ||
     (costheta>low_lim4_ec->Eval(phi_deg) && costheta<up_lim4_ec->Eval(phi_deg) && costheta<rightside_lim4_ec->Eval(phi_deg) && costheta<leftside_lim4_ec->Eval(phi_deg) && !hot_spot)       ||
     (costheta>low_lim5_ec->Eval(phi_deg) && costheta<up_lim5_ec->Eval(phi_deg) && costheta<rightside_lim5_ec->Eval(phi_deg) && costheta<leftside_lim5_ec->Eval(phi_deg))       ||
     (costheta>low_lim6_ec->Eval(phi_deg) && costheta<up_lim6_ec->Eval(phi_deg) && costheta<rightside_lim6_ec->Eval(phi_deg) && costheta<leftside_lim6_ec->Eval(phi_deg))
 ){ //EC only
      status=true;}
    else {
      status=false;}

     return status;
}



bool Pi_phot_fid_united(std::string beam_en, TVector3 V3_pi_phot, int q_pi_phot){

  bool status=false;
  std::string fbeam_en = beam_en;
  Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;

  if(q_pi_phot==0) status=Phot_fid(V3_pi_phot);
  if(q_pi_phot>0) status=PiplFiducialCut(fbeam_en, V3_pi_phot, &pi_cphil, &pi_cphir);
  if(q_pi_phot<0) status=PimiFiducialCut(fbeam_en, V3_pi_phot, &pi_phimin, &pi_phimax);
  return status;


}






void  prot3_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot[3],TVector3  V3prot_uncorr[3],TLorentzVector V4el,double Ecal_3pto2p[][2],double  pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3],double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det){

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






void  prot2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot[2],TVector3  V3prot_uncorr[2],TLorentzVector V4el,double Ecal_2pto1p[2],double  pmiss_perp_2pto1p[2],double  P2pto1p[2], double *Nboth){

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






void prot1_pi1_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi, int q_pi,bool radstat, double *N_pi_p,double *N_nopi_p){

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






void prot1_pi2_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[2], int q_pi[2],bool radstat[2], double *P_1p0pi,double P_1p1pi[2]){

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




void prot1_pi3_rot_func(std::string beam_en, TVector3 V3q, TVector3  V3prot,TVector3 V3pi[3], int q_pi[3],bool radstat[3],double *P_tot){

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




void prot2_pi1_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_2p1pi_to2p0pi[2],double p_miss_perp_2p1pi_to2p0pi[2],double P_2p1pito2p0pi[2],double P_2p1pito1p1pi[2],double P_2p1pito1p0pi[2],double *P_tot){

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


void prot2_pi2_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], int q_pi[2],bool radstat[2],TLorentzVector V4_el, double Ecal_2p2pi[2],double p_miss_perp_2p2pi[2],double P_tot_2p[2]){

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




void prot3_pi1_rot_func(std::string beam_en, TVector3 V3_q,TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, int q_pi,bool radstat,TLorentzVector V4_el, double Ecal_3p1pi[3],double p_miss_perp_3p1pi[3],double P_tot_3p[3]){

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









void pi1_rot_func(std::string beam_en, TVector3 V3_pi, int q_pi,bool radstat,TVector3 V3_q,double *P_pi){

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




void pi2_rot_func(std::string beam_en, TVector3 V3_pi[2], int q_pi[2],bool radstat[2], TVector3 V3_q,double *P_0pi,double P_1pi[2]){

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





void pi3_rot_func(std::string beam_en, TVector3 V3_pi[3], int q_pi[3],bool radstat[3], TVector3 V3_q,double *P_0pi,double P_1pi[3],double P_320[3],double P_3210[][2]){

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







void pi4_rot_func(std::string beam_en, TVector3 V3_pi[4], int q_pi[4],bool radstat[4], TVector3 V3_q,double *P_0pi,double *P_410,double *P_420,double *P_4210,double *P_430,double *P_4310,double *P_4320,double *P_43210){

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
