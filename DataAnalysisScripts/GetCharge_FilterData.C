#define GETCHARGE_FILTERDATA_C

//#include "e2a_ep_neutrino6_united4_radphot.h"
#include "GetCharge_FilterData.h"
#include "Constants.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>

#include <fstream>
#include <exception>
#include <iostream>
#include <map>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <utility>

TString ToString(int num) {

	std::ostringstream start;
	start << num;
	std::string start1 = start.str();
	return start1;

}

// Also used by e2a_ep_neutrino6_united4_radphot.{C,h}

extern TF1 *pipl_deltat_sig,*pipl_deltat_mean,*pimi_deltat_sig,*pimi_deltat_mean, *prot_deltat_sig, *prot_deltat_mean,*el_Epratio_sig,*el_Epratio_mean;

void SetFiducialCutParameters(std::string beam_en); // Load Fidicual Parameters for 1.1 and 4.4 GeV from file

extern double vz_corr(TF1 *vz_corr_func, double phi,double theta);
extern TVector3 FindUVW(TVector3 xyz);
extern Bool_t CutUVW(TVector3 ecxyz);
extern Float_t ProtonMomCorrection_He3_4Cell(std::string atarget, TLorentzVector V4Pr, Float_t vertex_p);

//e- E_tot/p vs p PID cut
extern Double_t FSum_e(Double_t *x,Double_t *par);

extern Double_t FSub_e(Double_t *x,Double_t *par);

//proton Delta_t vs momentum PID cut

extern Double_t FSum_prot(Double_t *x, Double_t *par);
extern Double_t FSub_prot(Double_t *x,Double_t *par);

//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions

extern Double_t FSum_pimi(Double_t *x,Double_t *par);

extern Double_t FSub_pimi(Double_t *x,Double_t *par);

//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions
extern Double_t FSum_pipl(Double_t *x,Double_t *par);

extern Double_t FSub_pipl(Double_t *x,Double_t *par);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------

void GetCharge_FilterData::Loop()
{

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

	target_name = ftarget;   //std string for target name

	// -----------------------------------------------------------------------------------------------------

	// apapadop

	int TargetPdgCode, TargetZ, TargetA;

	if (ftarget=="3He") { TargetPdgCode = 1000020030; TargetZ = 2; TargetA = 3; }
	if (ftarget=="4He") { TargetPdgCode = 1000020040; TargetZ = 2; TargetA = 4; }
	if (ftarget=="CH2") { TargetPdgCode = 1000080140; TargetZ = 8; TargetA = 14; }
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

	int GetCharge_CounterEvents = 0;
	int EC_CounterEvents = 0;
	int SC_CounterEvents = 0;
	int ElectronFid_Counter = 0;
	int CutUVW_Counter = 0;
	int EC_SC_CC_Q_Counter = 0;
	int EnergyDep_Counter = 0;
	int Vertex_Counter = 0;

	Long64_t nentries = fChain->GetEntriesFast();
	//	nentries = 8000000;

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


	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV Configuration parameters and cuts
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

		// ------------------------------------

		// Mariana's fitting

//		vert_min["3He"]=-3.05;
//		vert_min["C12"]=4.95;
//		vert_min["CH2"]=4.85;

//		vert_max["3He"]=-0.18;
//		vert_max["C12"]=5.76;
//		vert_max["CH2"]=5.62;

		// apapadop's fitting // Nov 9 2020

//		vert_min["3He"]=-3.05;
//		vert_min["C12"]=4.43;
//		vert_min["CH2"]=4.64;

//		vert_max["3He"]=-0.18;
//		vert_max["C12"]=6.16;
//		vert_max["CH2"]=5.82;

		// Or's eyeballing // Nov 10 2020

		vert_min["3He"]=-3.05;
		vert_min["C12"]=4.3;
		vert_min["CH2"]=4.3;

		vert_max["3He"]=-0.18;
		vert_max["C12"]=6.5;
		vert_max["CH2"]=6.5;

		// ------------------------------------

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

		// -----------------------------------------

		// Mariana's fitting

//		vert_min["3He"]=-3.29;
//		vert_min["4He"]=-2.53;
//		vert_min["C12"]=4.8;
//		vert_min["56Fe"]=4.6;

//		vert_max["3He"]=-0.23;
//		vert_max["4He"]=1.73;
//		vert_max["C12"]=5.5;
//		vert_max["56Fe"]=5.3;

		// apapadop's fitting // Nov 9 2020

//		vert_min["3He"]=-3.29;
//		vert_min["4He"]=-2.53;
//		vert_min["C12"]=4.58;
//		vert_min["56Fe"]=4.46;

//		vert_max["3He"]=-0.23;
//		vert_max["4He"]=1.73;
//		vert_max["C12"]=5.6;
//		vert_max["56Fe"]=5.47;

		// Or's eyeballing // Nov 10 2020

		vert_min["3He"]=-3.29;
		vert_min["4He"]=-2.53;
		vert_min["C12"]=4.3;
		vert_min["56Fe"]=4.3;

		vert_max["3He"]=-0.23;
		vert_max["4He"]=1.73;
		vert_max["C12"]=6.5;
		vert_max["56Fe"]=6.5;

		// -----------------------------------------

		vertdiff_min["3He"]=-1.;
		vertdiff_min["4He"]=-1.;
		vertdiff_min["C12"]=-1.;
		vertdiff_min["CH2"]=-1.;
		vertdiff_min["56Fe"]=-1.;

		vertdiff_max["3He"]=1.;
		vertdiff_max["4He"]=1.;
		vertdiff_max["C12"]=1.;
		vertdiff_max["CH2"]=1.;
		vertdiff_max["56Fe"]=1.;

		EC_photon_beta["3He"]=0.93;
		EC_photon_beta["4He"]=0.92;
		EC_photon_beta["C12"]=0.92;
		EC_photon_beta["CH2"]=0.92;
		EC_photon_beta["56Fe"]=0.90;

		LEC_photon_beta["3He"]=0.96;
		LEC_photon_beta["4He"]=0.94;
		LEC_photon_beta["C12"]=0.94;
		LEC_photon_beta["CH2"]=0.94;
		LEC_photon_beta["56Fe"]=0.95;

		EC_time_offset[std::make_pair("3He",1)]=-1.37;  EC_time_offset[std::make_pair("3He",2)]=-1.42; EC_time_offset[std::make_pair("3He",3)]=-1.55;
		EC_time_offset[std::make_pair("3He",4)]=-1.53;  EC_time_offset[std::make_pair("3He",5)]=-1.49; EC_time_offset[std::make_pair("3He",6)]=-1.44;

		EC_time_offset[std::make_pair("4He",1)]=0.72;  EC_time_offset[std::make_pair("4He",2)]=0.27; EC_time_offset[std::make_pair("4He",3)]=0.16;
		EC_time_offset[std::make_pair("4He",4)]=0.21;  EC_time_offset[std::make_pair("4He",5)]=0.22; EC_time_offset[std::make_pair("4He",6)]=0.21;

		EC_time_offset[std::make_pair("C12",1)]=0.50;  EC_time_offset[std::make_pair("C12",2)]=0.39; EC_time_offset[std::make_pair("C12",3)]=0.29;
		EC_time_offset[std::make_pair("C12",4)]=0.29;  EC_time_offset[std::make_pair("C12",5)]=0.32; EC_time_offset[std::make_pair("C12",6)]=0.33;

		EC_time_offset[std::make_pair("CH2",1)]=0.50;  EC_time_offset[std::make_pair("CH2",2)]=0.39; EC_time_offset[std::make_pair("CH2",3)]=0.29;
		EC_time_offset[std::make_pair("CH2",4)]=0.29;  EC_time_offset[std::make_pair("CH2",5)]=0.32; EC_time_offset[std::make_pair("CH2",6)]=0.33;

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

		// -----------------------------------------

		// Mariana's fitting

//		vert_min["3He"]=-3.27;
//		vert_min["4He"]=-2.51;
//		vert_min["C12"]=4.7;
//		vert_min["56Fe"]=4.6;

//		vert_max["3He"]=0.07;
//		vert_max["4He"]=1.71;
//		vert_max["C12"]=5.3;
//		vert_max["56Fe"]=5.4;

		// apapadop's fitting // Nov 9 2020

//		vert_min["3He"]=-3.27;
//		vert_min["4He"]=-2.51;
//		vert_min["C12"]=4.44;
//		vert_min["56Fe"]=4.41;

//		vert_max["3He"]=0.07;
//		vert_max["4He"]=1.71;
//		vert_max["C12"]=5.56;
//		vert_max["56Fe"]=5.46;

		// Or's eyeballing // Nov 10 2020

		vert_min["3He"]=-3.27;
		vert_min["4He"]=-2.51;
		vert_min["C12"]=4.3;
		vert_min["56Fe"]=4.3;

		vert_max["3He"]=0.07;
		vert_max["4He"]=1.71;
		vert_max["C12"]=6.5;
		vert_max["56Fe"]=6.5;

		// -----------------------------------------

		vertdiff_min["3He"]=-1.;
		vertdiff_min["4He"]=-1;
		vertdiff_min["C12"]=-1;
		vertdiff_min["CH2"]=-1;
		vertdiff_min["56Fe"]=-1;

		vertdiff_max["3He"]=1.;
		vertdiff_max["4He"]=1.;
		vertdiff_max["C12"]=1;
		vertdiff_max["CH2"]=1;
		vertdiff_max["56Fe"]=1;

		EC_photon_beta["3He"]=0.92;
		EC_photon_beta["4He"]=0.91;
		EC_photon_beta["C12"]=0.92;
		EC_photon_beta["CH2"]=0.92;
		EC_photon_beta["56Fe"]=0.91;

		LEC_photon_beta["3He"]=0.97;
		LEC_photon_beta["4He"]=0.97;
		LEC_photon_beta["C12"]=0.95;
		LEC_photon_beta["CH2"]=0.95;
		LEC_photon_beta["56Fe"]=0.96;

		EC_time_offset[std::make_pair("3He",1)]=-0.15;  EC_time_offset[std::make_pair("3He",2)]=-0.26; EC_time_offset[std::make_pair("3He",3)]=-0.41;
		EC_time_offset[std::make_pair("3He",4)]=-0.29;  EC_time_offset[std::make_pair("3He",5)]=-0.25; EC_time_offset[std::make_pair("3He",6)]=-0.23;

		EC_time_offset[std::make_pair("4He",1)]=-0.01;  EC_time_offset[std::make_pair("4He",2)]=-0.11; EC_time_offset[std::make_pair("4He",3)]=-0.23;
		EC_time_offset[std::make_pair("4He",4)]=-0.26;  EC_time_offset[std::make_pair("4He",5)]=-0.21; EC_time_offset[std::make_pair("4He",6)]=-0.09;

		EC_time_offset[std::make_pair("C12",1)]=-0.01;  EC_time_offset[std::make_pair("C12",2)]=-0.11; EC_time_offset[std::make_pair("C12",3)]=-0.23;
		EC_time_offset[std::make_pair("C12",4)]=-0.27;  EC_time_offset[std::make_pair("C12",5)]=-0.21; EC_time_offset[std::make_pair("C12",6)]=-0.08;

		EC_time_offset[std::make_pair("CH2",1)]=-0.01;  EC_time_offset[std::make_pair("CH2",2)]=-0.11; EC_time_offset[std::make_pair("CH2",3)]=-0.23;
		EC_time_offset[std::make_pair("CH2",4)]=-0.27;  EC_time_offset[std::make_pair("CH2",5)]=-0.21; EC_time_offset[std::make_pair("CH2",6)]=-0.08;

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
	bind_en["C12"] = C12_bind_en-B_bind_en	+ Ecal_offset["C12"];
	bind_en["56Fe"]= Fe_bind_en-Mn_bind_en	+ Ecal_offset["56Fe"];
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

	TFile *file_in5 = NULL;
	TFile *file_in6 = NULL;
	TFile *file_in7 = NULL;

	if (en_beam[fbeam_en] < 2.300 && en_beam[fbeam_en] > 2.100 ) {
		file_in5 = new TFile("FiducialsCorrections/vz_56Fe_2261_badruns.root");//vertex correction for 56Fe runs with exploded liquid target cell
		file_in6 = new TFile("FiducialsCorrections/vz_3He_2261_2ndrungroup.root");//vertx correction for 3He 2nd group runs
		file_in7 = new TFile("FiducialsCorrections/vz_4He_2261_2ndrungroup.root");//vertex correction for 4He 2nd group runs
	}

	//Output file definition

	// e4v analysis files
//	TString FileName = Form("/w/hallb-scifs17exp/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_%s_%s_neutrino6_united4_radphot_test_100M.root",ftarget.c_str(),fbeam_en.c_str());

	// e4v workshop files
	TString FileName = Form("/w/hallb-scifs17exp/clas/claseg2/apapadop/e4vWorkshop_%s_%s.root",ftarget.c_str(),fbeam_en.c_str());

	std::cout << "File " << FileName << " to be  created" << std::endl;

	TFile *file_out = new TFile(FileName, "Recreate");

	// -------------------------------------------------------------------------------------------------------------------------------------------------------

	// apapadop
	TTree* mytree = new TTree("gst","gst");

	Double_t        genie_q_l;
	Int_t           genie_RunNumber;

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
	Double_t        genie_thetal;
	Double_t        genie_phil;

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
	Double_t        genie_thetaf[FinalStateParticles];   //[nf]
	Double_t        genie_phif[FinalStateParticles];   //[nf]

	Double_t        genie_vtxx;
	Double_t        genie_vtxy;
	Double_t        genie_vtxz;
	Double_t        genie_vtxt;
	Double_t        genie_sumKEf;
	Double_t        genie_calresp0;

	mytree->Branch("q_l", &genie_q_l, "q_l/D");
	mytree->Branch("RunNumber", &genie_RunNumber, "RunNumber/I");

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
	mytree->Branch("thetal", &genie_thetal, "thetal/D");
	mytree->Branch("phil", &genie_phil, "phil/D");

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
	mytree->Branch("thetaf", &genie_thetaf, "thetaf[120]/D");
	mytree->Branch("phif", &genie_phif, "phif[120]/D");

	mytree->Branch("vtxx", &genie_vtxx, "vtxx/D");
	mytree->Branch("vtxy", &genie_vtxy, "vtxy/D");
	mytree->Branch("vtxz", &genie_vtxz, "vtxz/D");
	mytree->Branch("vtxt", &genie_vtxt, "vtxt/D");
	mytree->Branch("sumKEf", &genie_sumKEf, "sumKEf/D");
	mytree->Branch("calresp0", &genie_calresp0, "calresp0/D");

	int shift = 60;

	// -------------------------------------------------------------------------------------------------------------------------------------------------------

	double pars[3];
	if (en_beam[fbeam_en] < 2.300 && en_beam[fbeam_en] > 2.100 ) {
		if(ftarget=="56Fe") {
			pars[0]=((TF1 *)file_in5->Get("f_vz"))->GetParameter(0);
			pars[1]=((TF1 *)file_in5->Get("f_vz"))->GetParameter(1);
			pars[2]=((TF1 *)file_in5->Get("f_vz"))->GetParameter(2);
		}
		if(ftarget=="3He") {
			pars[0]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(0);
			pars[1]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(1);
			pars[2]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(2);
		}
		if(ftarget=="4He"){
		 pars[0]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(0);
		 pars[1]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(1);
		 pars[2]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(2);
		}

	} else {
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

	fiducialcut->InitPiMinusFit(fbeam_en);

	//initialize Fiducial functions for EC limits

	fiducialcut->InitEClimits();

	// ------------------------------------------------------------------------------------------------------------------------------------------------

	// Histo declaration

	TH2D* h2_Electron_UncorrectedVertex_Phi = new TH2D("h2_Electron_UncorrectedVertex_Phi",";Uncorrected Vertex [cm];#phi_{e'} [deg]",2000,-10,10,360,0,360);
	TH2D* h2_Electron_CorrectedVertex_Phi = new TH2D("h2_Electron_CorrectedVertex_Phi",";Corrected Vertex [cm];#phi_{e'} [deg]",2000,-10,10,360,0,360);

	TH1D* h1_el_timediff = new TH1D("h1_el_timediff",";el_timediff",100,-50,50);
	TH1D* h1_cc_c2 = new TH1D("h1_cc_c2",";cc_c2",2000,0,20);
	TH1D* h1_ece = new TH1D("h1_ece",";ece",2000,0,20);
	TH1D* h1_eceOverP = new TH1D("h1_eceOverP",";ece/p",200,0,2);
	TH1D* h1_ec_ei = new TH1D("h1_ec_ei",";ec_ei",700,0,0.7);

	TH1D* h1_el_timediff_Sector[6]; 

	for (int i = 0; i < 6; i++) {

		h1_el_timediff_Sector[i] = new TH1D("h1_el_timediff_Sector_"+ToString(i),";el_timediff",100,-50,50);

	}

	// ------------------------------------------------------------------------------------------------------------------------------------------------

	std::cout << "Initial number of events = " << fChain->GetEntries() << std::endl;

	// ---------------------------------------------------------------------------------------------------------------

	// Justification for the parameter choice
	// https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing

	TF1 *myElectronFit = new TF1("myElectronFit","[0]+[1]/x",0.,5.);

	if (en_beam[fbeam_en] == 1.161) { myElectronFit->SetParameters(17,7); }
	if (en_beam[fbeam_en] == 2.261) { myElectronFit->SetParameters(16,10.5); }
	if (en_beam[fbeam_en] == 4.461) { myElectronFit->SetParameters(13.5,15); }

	// ------------------------------------------------------------------------------------------------------------------------------------------------

	// Starting the loop

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		int nb = GetEntry(jentry);
		if (ientry < 0) break;

		if (jentry%1000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/fChain->GetEntries()*100. << " %"<< std::endl;}

		if( jentry%200000 == 0 )
		{
			gDirectory->Write("hist_Files", TObject::kOverwrite);
			std::cout<<jentry<<std::endl;
		}

		if (runnb==18258 || runnb==18259 || (runnb>18382 && runnb<18438) || (runnb>18220 && runnb<18253)) {
			vz_corr_func->SetParameters(pars); 
			//setting appropriate parameters of the vertex correction function for the runs with the same target and beam energy, but different vertex correction
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

		if((runnb>18283 && runnb<18289) || (runnb>18300 && runnb<18304) || (runnb>18317 && runnb<18329)) fTorusCurrent=750; //setting appropriate torrus magnet current
		else if ( 
			(runnb>18293 && runnb<18301) || (runnb>18305 && runnb<18317) || (runnb>18328 && runnb<18336) 
			|| runnb == 18334 /*CH2*/
			)  
			{ fTorusCurrent = 1500; }
		else fTorusCurrent=2250;

		// apapadop, CH2 only run @ 1.1 GeV with 1500 torus current, but we don't have fiducials there, thus using the 750 torus current
		if (runnb == 18334) { fTorusCurrent = 750; } 

		// apapadop, we want only the 750 results for the main e4v paper
		// 1500 to be used only for the hydrogen study
		if (fbeam_en == "1161" && fTorusCurrent > 760) { continue; }  
	
		// For the 1500 (high torus current runs) on 12C @ 1.1 GeV, comment in the lines below                                                            
                //if (fbeam_en == "1161" && fTorusCurrent < 760) { continue; }
		//fTorusCurrent = 750;

		if(jentry == 0){ //was n_evt == 1 before but jentry = n_evnt - 1
			//SetMomCorrParameters(); Functions is missing F.H. 08/01/19
			fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
			fiducialcut->SetFiducialCutParameters(fbeam_en);
		}

		int n_elec = 0;
		const int ind_em=0; //Index for electron
		if (ec[ind_em] <=0) {
//			std::cout << "Possible problem with making electron ec vector. EC index below/equal Zero: ec[ind_em] =  " << ec[ind_em] << std::endl;
			continue;
		}

		EC_CounterEvents++;

		if (sc[ind_em] <=0) {
//			std::cout << "Possible problem with making electron ec vector. SC index below/equal zero: sc[ind_em] =  " << sc[ind_em] << std::endl;
			continue;
		}

		SC_CounterEvents++;

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
			genie_thetaf[WhichFinalStateParticle] = -99.;
			genie_phif[WhichFinalStateParticle] = -99.;

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		//Define electron vectors, angles and other Information

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

		h2_Electron_UncorrectedVertex_Phi->Fill(el_vert,el_phi_mod);
		h2_Electron_CorrectedVertex_Phi->Fill(el_vert_corr,el_phi_mod);

		// Variables for electron cuts

		double ece = TMath::Max( ec_ei[ec[ind_em] - 1] + ec_eo[ec[ind_em] - 1],   etot[ec[ind_em] - 1]);
		el_segment = int((cc_segm[cc[ind_em]-1]-int(cc_segm[cc[ind_em]-1]/1000)*1000)/10); //does this work in all cases?? F.H. 08/07/19
		el_cc_sector = cc_sect[cc[ind_em]-1];
		el_sccc_timediff = sc_t[cc[ind_em]-1]-cc_t[cc[ind_em]-1]-(sc_r[cc[ind_em]-1]-cc_r[cc[ind_em]-1])/(c*ns_to_s);
		el_cc_nphe = nphe[cc[ind_em]-1]/10.;
		double ec_SC_timediff_uncorr = ec_t[ec[ind_em]-1]-sc_t[sc[ind_em]-1]-(ec_r[ec[ind_em]-1]-sc_r[sc[ind_em]-1])/(c*ns_to_s);

		// fsum_e and fsub_p are TF1 Functions for electron E/p cuts

		fsum_e->SetParameters(epratio_sig_cutrange, max_mom);
		fsub_e->SetParameters(epratio_sig_cutrange, max_mom);

		// Main Fiducial Cut for Electron

		if( !EFiducialCut(fbeam_en, el_mom1) ) continue;//theta, phi cuts

		// Extra layer of fiducials for electrons

		double theta_min = myElectronFit->Eval(el_mom1.Mag());
		if (el_theta < theta_min) { continue; }

		ElectronFid_Counter++;

		if( !CutUVW(e_ec_xyz1) ) continue; //u>60, v<360, w<400

		CutUVW_Counter++;

		// General cut on EC, SC, CC hit and q (charge) for all events

		if( ec[ind_em] < 0.5 ||  sc[ind_em] < 0.5 ||  cc[ind_em] < 0.5 || q[ind_em] >= 0) {
		  
			continue;

		}

		EC_SC_CC_Q_Counter++;

		h1_el_timediff->Fill(el_sccc_timediff);
		h1_cc_c2->Fill(cc_c2[cc[ind_em]-1]);
		h1_ec_ei->Fill(ec_ei[ec[ind_em] - 1]);
		h1_ece->Fill(ece);
		h1_eceOverP->Fill(ece/p[ind_em]);
		h1_el_timediff_Sector[(int)(el_phi_mod/60.)]->Fill(el_sccc_timediff);

		// Cut on 1.1 GeV events (E/p, energy deposit, TOF and cherenkov)

		if(en_beam[fbeam_en] > 1. && en_beam[fbeam_en] <2. &&
		  ( ec_ei[ec[ind_em] - 1] < 0.03 || ece/p[ind_em] < fsub_e->Eval(p[ind_em]) || ece/p[ind_em] > fsum_e->Eval(p[ind_em]) ||
				p[ind_em] < min_good_mom || el_sccc_timediff < sc_cc_delt_cut_sect[el_cc_sector-1] ||   cc_c2[cc[ind_em]-1] > 0.1 ) )
		{
				continue;
		}

		// Cut on 2.2 GeV events (E/p, energy deposit, TOF and cherenkov)

		if(en_beam[fbeam_en] < 3.  && en_beam[fbeam_en] > 2 &&
		  ( ec_ei[ec[ind_em] - 1] < 0.06 || ece/p[ind_em] < fsub_e->Eval(p[ind_em]) || ece/p[ind_em] > fsum_e->Eval(p[ind_em]) ||
				p[ind_em] < min_good_mom || cc_c2[cc[ind_em]-1] >= 0.1 || el_sccc_timediff < sc_cc_delt_cut_sect[el_cc_sector-1] ||
				TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass)>en_beam[fbeam_en] ) ) 
		//only here a cut on electron momentum to cut some very scarse events where p_e > beam energy (see Mariana's anaysis note)
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

		EnergyDep_Counter++;

		// --------------------------------------------------------------------------------------------

		// apapadop

		//Electron vertex cut
		if( !(el_vert_corr < vert_max[ftarget] && el_vert_corr > vert_min[ftarget]) ) continue;

		Vertex_Counter++;

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

		genie_x = x_bjk;
		genie_Q2reco = Q2;
		genie_W = W_var;
		genie_El = V4_el.E();
		genie_pxl = V4_el.Px();
		genie_pyl = V4_el.Py();
		genie_pzl = V4_el.Pz();
		genie_pl = V4_el.Rho();
		genie_cthl = cos(V4_el.Theta());
		genie_thetal = V4_el.Theta()*180./TMath::Pi();

		double ElectronPhi = V4_el.Phi()*180./TMath::Pi()+30;
		if (ElectronPhi < 0) { ElectronPhi += 360.; }
		if (ElectronPhi >360) { ElectronPhi -= 360.; }
		genie_phil = ElectronPhi;

		genie_vtxx = 0;
		genie_vtxy = 0;
		genie_vtxz = el_vert_corr;

		// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		//Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

		//Index variables for hadrons (p and pions)
		int index_p[20]; //index for each proton
		int ind_p;		  //temp proton index in the loop
		int index_pi[20]; //index for each pion
		int ind_pi_phot[20];
		int ind_pi;		   //temp pion index in the loop
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

		GetCharge_CounterEvents++;

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
				if(delta < fsum_prot->Eval(p[ind_p]) && delta > fsub_prot->Eval(p[ind_p]) && p[ind_p] >= prot_accept_mom_lim) {

					TLorentzVector V4_uncorrprot(p[ind_p]*cx[ind_p],p[ind_p]*cy[ind_p],p[ind_p]*cz[ind_p],TMath::Sqrt(p[ind_p]*p[ind_p]+ m_prot*m_prot ) );
					p_vert_corr = vz[ind_p]+vz_corr(vz_corr_func,prot_phi_mod,prot_theta);

					// Proton fiducial cuts & extra outline

					if(PFiducialCut(fbeam_en, V4_uncorrprot.Vect()) && PFiducialCutExtra(fbeam_en,V4_uncorrprot.Vect()) ) {

						//main vertex cut for protons
						if( (el_vert_corr-p_vert_corr) > vertdiff_min[ftarget] && (el_vert_corr-p_vert_corr) < vertdiff_max[ftarget] ) {

							num_p = num_p + 1;
							index_p[num_p-1] = i;

							// ------------------------------------------------------------------------------------------------------------------------------------------

							// apapadop

							genie_pdgf[IndexInFinalStateParticleArray] = id[ind_p];
							genie_Ef[IndexInFinalStateParticleArray] = V4_uncorrprot.E();
							genie_pxf[IndexInFinalStateParticleArray] = V4_uncorrprot.Px();
							genie_pyf[IndexInFinalStateParticleArray] = V4_uncorrprot.Py();
							genie_pzf[IndexInFinalStateParticleArray] = V4_uncorrprot.Pz();
							genie_pf[IndexInFinalStateParticleArray] = V4_uncorrprot.Rho();
							genie_cthf[IndexInFinalStateParticleArray] = V4_uncorrprot.CosTheta();
							genie_thetaf[IndexInFinalStateParticleArray] = V4_uncorrprot.Theta()*180./TMath::Pi();

							double UncorrPhi = V4_uncorrprot.Phi()*180./TMath::Pi()+30;
							if (UncorrPhi < 0) { UncorrPhi += 360.; }
							if (UncorrPhi > 360) { UncorrPhi -= 360.; }
							genie_phif[IndexInFinalStateParticleArray] = UncorrPhi;

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
							genie_thetaf[IndexInFinalStateParticleArray+shift] = V4_uncorrprot.Theta()*180./TMath::Pi();

							double ShiftUncorrPhi = V4_uncorrprot.Phi()*180./TMath::Pi()+30;
							if (ShiftUncorrPhi < 0) { ShiftUncorrPhi += 360.; }
							if (ShiftUncorrPhi > 360) { ShiftUncorrPhi -= 360.; }
							genie_phif[IndexInFinalStateParticleArray+shift] = ShiftUncorrPhi;

							IndexInFinalStateParticleArray++;

							// ---------------------------------------------------------------------------------------------------------------------------------

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

				if(delta < fsum_pimi->Eval(p[i]) && delta > fsub_pimi->Eval(p[i])) {

					pimi_vert=vz[i];
					pimi_phi = TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
					pimi_phi_mod = pimi_phi + 30;  //Add extra 30 degree rotation in phi
					if (pimi_phi_mod<0)  pimi_phi_mod = pimi_phi_mod + 360;  //Pi minus is between 0 and 360 degree
					pimi_theta = TMath::ACos(cz[i])*TMath::RadToDeg();
					pimi_vert_corr = pimi_vert+vz_corr(vz_corr_func, pimi_phi_mod,pimi_theta);

					// Pi minus fiducial cuts & extra outline 

					if (PimiFiducialCut(fbeam_en, V3_pimi, &pimi_phimin, &pimi_phimax) && PimiFiducialCutExtra(fbeam_en, V3_pimi) ) {  

						//main vertex cut for pi minus
						if(fabs(el_vert_corr-pimi_vert_corr) < pimi_vertcut) {

							num_pimi = num_pimi + 1;
							num_pi = num_pi + 1;
							num_pi_phot = num_pi_phot + 1;
							num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
							index_pimi[num_pimi - 1] = i;
							index_pi[num_pi - 1] = i;
							ind_pi_phot[num_pi_phot - 1] = i;

							// --------------------------------------------------------------------------------------------------------------------------------------

							// apapadop

							double PiMinusP = V3_pimi.Mag();
							double PiMinusE = TMath::Sqrt(PiMinusP*PiMinusP + m_pimi * m_pimi);

							genie_pdgf[IndexInFinalStateParticleArray] = -211;
							genie_Ef[IndexInFinalStateParticleArray] = PiMinusE;
							genie_pxf[IndexInFinalStateParticleArray] = V3_pimi.X();
							genie_pyf[IndexInFinalStateParticleArray] = V3_pimi.Y();
							genie_pzf[IndexInFinalStateParticleArray] = V3_pimi.Z();
							genie_pf[IndexInFinalStateParticleArray] = PiMinusP;
							genie_cthf[IndexInFinalStateParticleArray] = V3_pimi.CosTheta();
							genie_thetaf[IndexInFinalStateParticleArray] = V3_pimi.Theta()*180./TMath::Pi();

							double CorrPhi = V3_pimi.Phi()*180./TMath::Pi() + 30;
							if (CorrPhi < 0) { CorrPhi += 360.; }
							if (CorrPhi > 360) { CorrPhi -= 360.; }

							genie_phif[IndexInFinalStateParticleArray] = CorrPhi;

							IndexInFinalStateParticleArray++;

							// ----------------------------------------------------------------------------------------------------------------------------------------


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
				if(delta < fsum_pipl->Eval(p[i]) && delta > fsub_pipl->Eval(p[i])) {

					pipl_vert=vz[i];
					pipl_phi = TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
	 				pipl_phi_mod = pipl_phi + 30; //Add 30 degrees
					if (pipl_phi_mod < 0)  pipl_phi_mod = pipl_phi_mod + 360;  //Pi plus is between 0 and 360 degree
					pipl_theta = TMath::ACos(cz[i])*TMath::RadToDeg();
					pipl_vert_corr = pipl_vert + vz_corr(vz_corr_func,pipl_phi_mod,pipl_theta);

					// Pi Plus fiducial cut & extra outline

					if ( PiplFiducialCut(fbeam_en, V3_pipl, &pipl_phimin, &pipl_phimax) && PiplFiducialCutExtra(fbeam_en, V3_pipl)){

						if (fabs(el_vert_corr-pipl_vert_corr) < pipl_vertcut) { //pi plus vertex cut
							num_pipl = num_pipl + 1;
							num_pi  = num_pi + 1;
							num_pi_phot = num_pi_phot + 1;
							num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
							index_pipl[num_pipl - 1] = i;
							index_pi[num_pi - 1] = i;
							ind_pi_phot[num_pi_phot - 1] = i;

							// ------------------------------------------------------------------------------------------------------------------------

							// apapadop

							double PiPlusP = V3_pipl.Mag();
							double PiPlusE = TMath::Sqrt(PiPlusP*PiPlusP + m_pipl * m_pipl);

							genie_pdgf[IndexInFinalStateParticleArray] = 211;
							genie_Ef[IndexInFinalStateParticleArray] = PiPlusE;
							genie_pxf[IndexInFinalStateParticleArray] = V3_pipl.X();
							genie_pyf[IndexInFinalStateParticleArray] = V3_pipl.Y();
							genie_pzf[IndexInFinalStateParticleArray] = V3_pipl.Z();
							genie_pf[IndexInFinalStateParticleArray] = PiPlusP;
							genie_cthf[IndexInFinalStateParticleArray] = V3_pipl.CosTheta();
							genie_thetaf[IndexInFinalStateParticleArray] = V3_pipl.Theta()*180./TMath::Pi();

							double CorrPhi = V3_pipl.Phi()*180./TMath::Pi()+30;
							if (CorrPhi < 0) { CorrPhi += 360.; }
							if (CorrPhi > 360) { CorrPhi -= 360.; }

							genie_phif[IndexInFinalStateParticleArray] = CorrPhi;

							IndexInFinalStateParticleArray++;

							// -----------------------------------------------------------------------------------------------------------------------------------

						} //vert cut ends

					} //fidcut ends

				} //delta cut ends

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

					// Photon fiducial function

					if ( Phot_fid(V3_phot_angles) && Phot_fidExtra(V3_phot_angles) ) {

						ec_num_n = ec_num_n + 1;
						ec_index_n[ec_num_n - 1] = i;
						num_pi_phot = num_pi_phot + 1;
						ind_pi_phot[num_pi_phot - 1] = i;

						// -------------------------------------------------------------------------------------------------------------------------

						// apapadop

						genie_pdgf[IndexInFinalStateParticleArray] = 22;
						genie_Ef[IndexInFinalStateParticleArray] = V3_phot_angles.Mag();
						genie_pxf[IndexInFinalStateParticleArray] = V3_phot_angles.X();
						genie_pyf[IndexInFinalStateParticleArray] = V3_phot_angles.Y();
						genie_pzf[IndexInFinalStateParticleArray] = V3_phot_angles.Z();
						genie_pf[IndexInFinalStateParticleArray] = V3_phot_angles.Mag();
						genie_cthf[IndexInFinalStateParticleArray] = V3_phot_angles.CosTheta();
						genie_thetaf[IndexInFinalStateParticleArray] = V3_phot_angles.Theta()*180./TMath::Pi();

						double CorrPhi = V3_phot_angles.Phi()*180./TMath::Pi()+30;
						if (CorrPhi < 0) { CorrPhi += 360.; }
						if (CorrPhi > 180) { CorrPhi -= 360.; }

						genie_phif[IndexInFinalStateParticleArray] = CorrPhi;

						IndexInFinalStateParticleArray++;

						// -----------------------------------------------------------------------------------------------------------------------------

						//Photon EC energy deposit

						photon_ece = TMath::Max( ec_ei[ec[i] - 1] + ec_eo[ec[i] - 1],etot[ec[i] - 1]);

						//Cut on Radiation photon via angle with respect to the electron
						//within 40 degrees in theta and 30 degrees in phi
						if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg() < phot_rad_cut && fabs(neut_phi_mod-el_phi_mod) < phot_e_phidiffcut ) {
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

		genie_q_l = q_l;

		genie_RunNumber = runnb;
		genie_iev = NEventsTotal;
		NEventsTotal++;

		mytree->Fill();

		// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


	// This is the end of the filter

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	} //end of event loop (jentry)

	std::cout << std::endl << "File " << FileName << " has been created" << std::endl << std::endl;

	std::cout << "Initial number of events = " << fChain->GetEntries() << std::endl;

	std::cout << "EC_CounterEvents = " << EC_CounterEvents << std::endl;
	std::cout << "SC_CounterEvents = " << SC_CounterEvents << std::endl;
	std::cout << "ElectronFid_Counter = " << ElectronFid_Counter << std::endl;
	std::cout << "CutUVW_Counter = " << CutUVW_Counter << std::endl;
	std::cout << "EC_SC_CC_Q_Counter = " << EC_SC_CC_Q_Counter << std::endl;
	std::cout << "EnergyDep_Counter = " << EnergyDep_Counter << std::endl;
	std::cout << "Vertex_Counter = " << Vertex_Counter << std::endl;

	std::cout << "GetCharge_CounterEvents = " << GetCharge_CounterEvents << std::endl;

	gStyle->SetOptFit(1);
	gDirectory->Write("hist_Files", TObject::kOverwrite);

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Store number of initial / selected events

	TH1D* InitialNEvents = new TH1D("InitialNEvents","",1,0,1);
	InitialNEvents->SetBinContent(1,fChain->GetEntries());
	file_out->cd();
	InitialNEvents->Write();

	TH1D* SelectedNEvents = new TH1D("SelectedNEvents","",1,0,1);
	SelectedNEvents->SetBinContent(1,fChain->GetEntries());
	file_out->cd();
	SelectedNEvents->Write();

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

} //End Loop function

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
