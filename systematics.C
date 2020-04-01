#define SYSTEMATICS_C

#include "systematics.h"
#include "Constants.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TH1D.h>
#include <TMatrixD.h>
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
#include <TString.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

vector<double> CalculateCalKineVars(double ECal,TLorentzVector FSElectron) {

	vector<double> CalKineVars; CalKineVars.clear();

	TLorentzVector V4_beam_Cal(0,0,ECal,ECal);
	double nu_Cal = -(FSElectron-V4_beam_Cal).E();
	double Q2_Cal = -(FSElectron-V4_beam_Cal).Mag2();
	double x_bjk_Cal = Q2_Cal/(2*m_prot*nu_Cal);
	TVector3 V3_q_Cal = (V4_beam_Cal-FSElectron).Vect();
	double W_var_Cal = TMath::Sqrt((m_prot+nu_Cal)*(m_prot+nu_Cal)-V3_q_Cal*V3_q_Cal);

	CalKineVars.push_back(nu_Cal); // 0-th element: energy transfer using Ecal
	CalKineVars.push_back(Q2_Cal); // 1st element: Q2 using Ecal
	CalKineVars.push_back(x_bjk_Cal); // 2nd element: xB using Ecal
	CalKineVars.push_back(W_var_Cal); // 3rd element: invariant mass using Ecal

	return CalKineVars;

}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en

void systematics::Loop(Int_t choice, std::string tweak, int sigma) {

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

	//Choice = 0 is for analysis of CLAS data while choice = 1 is for the analysis of GENIE Simulation
	if (choice != 1 && choice != 0) {
		std::cout << "This parameter value is not implemented in systematics::Loop(). It should be either 0 or 1. The given value is " << choice << std::endl;
		std::exit(0);
	}

	std::map<std::string,double>bind_en;
	std::map<std::string,double>target_mass;
	std::map<std::string,double>residual_target_mass;
	std::map<std::string, double> Ecal_offset; //that might not be necessary for simulation data

	target_name = ftarget; //std string for target name
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

	//Resolutions for Smearing for GENIE simulation data
	double reso_p = 0.01; // smearing for the proton
	double reso_e = 0.005; // smearing for the electrons
	//double reso_pipl = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
	//double reso_pimi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
	double reso_pi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)

	// Resolution defined above seems to be insufficient at 1.1 GeV -> tripled it for all particles
	if(fbeam_en == "1161") { reso_p = 3*reso_p; reso_e = 3*reso_e; reso_pi = 3*reso_pi; }

	double Wcut = 2; //cut for all beam energies < 2
	double Q2cut = 0; // cut for 1.1 GeV > 0.1, for 2.2 GeV > 0.4 and 4.4 GeV > 0.8

	const int n_slice=3; // Stick to the 3 slices
	const double pperp_min[n_slice]={0.,0.2,0.4};
	const double pperp_max[n_slice]={0.2,0.4,10.};

//	const int n_slice=3; // Stick to the 3 slices
//	const double pperp_min[n_slice]={0.,0.3,10.};
//	const double pperp_max[n_slice]={0.3,10.,100.};

	TVector3 V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;

	TString E_acc_file;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) //1.1 GeV Configuration parameters and cuts
	{
		E_acc_file="1_161";
		Q2cut = 0.1;
	}


	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV Configuration parameters and cuts
	{
		E_acc_file="2_261";
		Q2cut = 0.4;
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) //4.4 GeV Configuration parameters and cuts
	{
		E_acc_file="4_461";
		Q2cut = 0.8;
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

	//Definition of Histograms
	TH1F *h1_Etot_p_bkgd_slice[n_slice], *h1_Erec_p_bkgd_slice[n_slice];
	TH1F *h1_Etot_Npi0[n_slice],*h1_Erec_Npi0[n_slice];
	TH1F *h1_Etot_Npi1[n_slice],*h1_Erec_Npi1[n_slice];

	TH1F *h1_Erec_bkgd_pipl_pimi_new_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact[n_slice];
	TH1F *h1_Etot_bkgd_pipl_pimi_fact_pipl[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pimi[n_slice];

	TH1F *h1_Etot_bkgd_1p2pi[n_slice], *h1_Erec_bkgd_1p2pi[n_slice];
	TH1F *h1_Etot_bkgd_1p2pi_1p0pi[n_slice],*h1_Erec_bkgd_1p2pi_1p0pi[n_slice];
	TH1F *h1_Etot_bkgd_1p3pi[n_slice], *h1_Erec_bkgd_1p3pi[n_slice];

	TH1F *h1_Etot_p_bkgd_slice_2p2pi[n_slice], *h1_Erec_p_bkgd_slice_2p2pi[n_slice];
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
	TFile* file_acceptance;
	TFile* file_acceptance_p;
	TFile* file_acceptance_pip;

	TString Target = "12C";
	if (ftarget.c_str() == "3He") { Target = "3He"; }
	if (ftarget.c_str() == "4He") { Target = "4He"; }

	if (choice == 1) { //Only need acceptance maps for GENIE simulation data
		file_acceptance = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+".root");
		file_acceptance_p = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_p.root");
		file_acceptance_pip = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_pip.root");
	}

	// ---------------------------------------------------------------------------------------------------------------

	// GENIE Systematic Uncertainties

	TFile *fweights = new TFile("/pnfs/uboone/persistent/users/apapadop/GenieProduction/R-3_0_6/Rad-R-3_0_6_Clas_"+TString(fbeam_en.c_str())+"GeV/G18_10a_02_11a/"+TString(ftarget.c_str())+"/weights_"+TString(tweak.c_str())+".root");

	TTree *tweights = (TTree*)fweights->Get(TString(tweak.c_str()));
	int NtweightsEntries = tweights->GetEntries();

	TArrayF* weights = NULL;
	double fArray;
	tweights->SetBranchAddress("weights", &weights);
	int Ntfileentries = fChain->GetEntries();

	// ---------------------------------------------------------------------------------------------------------------

	double XSecScale = 1.;
//	TFile* XSecFile = TFile::Open("/uboone/app/users/apapadop/R-3_0_6/mySplines/xsec_gxspl-FNALbig.root"); 

//	TGraph* gr = NULL;

//	if (XSecFile) {
//		TDirectory* dir = (TDirectory*)(XSecFile->Get("nu_mu_C12"));
//		gr = (TGraph*)dir->Get("tot_cc");
//	}

	// ---------------------------------------------------------------------------------------------------------------

	//Output file definition

	TFile *file_out;
	if (choice == 0) { file_out = new TFile(Form("data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");}
	else { 

		TString Dir = "Up";
		if (sigma == 0) { Dir = "Down"; } 
		file_out = new TFile("mySystematics"+Dir+"/genie_e2a_ep_"+TString(ftarget.c_str())+"_"+TString(fbeam_en.c_str())+"_neutrino6_united4_radphot_test_"+TString(tweak.c_str())+".root", "Recreate");
	}

	// ---------------------------------------------------------------------------------------------------------------

	//initialize Fiducial functions for EC limits
	fiducialcut->InitEClimits();
	std::cout << " Test InitEClimits Loop " << fiducialcut->up_lim1_ec->Eval(60) << std::endl;

	//Definition and initialization of Histograms
	TH1F *h1_el_Mott_crosssec = new TH1F("h1_el_Mott_crosssec","",200,0.,0.01);
	TH1F *h1_Wvar = new TH1F("h1_Wvar","",400,0,3);
	TH1F *h1_xbjk = new TH1F("h1_xbjk","",400,0,3);
	TH1F *h1_Q2 = new TH1F("h1_Q2","",400,0,6);
	TH1F *h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
	TH1F *h1_Nprot=new TH1F("h1_Nprot","",10,-0.5,4.5);
	TH1F *h1_Nprot_NonZeroProt=new TH1F("h1_Nprot_NonZeroProt","",8,0.5,4.5);
	TH1F *h1_Nphot=new TH1F("h1_Nphot","",10,-0.5,4.5);
	TH1F *h1_Npiphot=new TH1F("h1_Npiphot","",10,-0.5,4.5);
	TH1F *h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,-0.5,4.5);
	TH1F *h1_Npi=new TH1F("h1_Npi","",10,-0.5,4.5);
	TH1F *h1_Npi_NonZeroProt=new TH1F("h1_Npi_NonZeroProt","",10,-0.5,4.5);
	TH1F *h1_Npipl=new TH1F("h1_Npipl","",10,-0.5,4.5);
	TH1F *h1_Npimi=new TH1F("h1_Npimi","",10,-0.5,4.5);
	TH1F *h1_MissMomentum = new TH1F("MissMomentum","",80,0.,1.);
	TH1F *h1_el_mom = new TH1F("h1_el_mom","",100,0.2,6);
	TH1F *h1_el_mom_corr = new TH1F("h1_el_mom_corr","",100,0.,5.);
	TH1F *h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
	TH1F *h1_prot_mom = new TH1F("h1_prot_mom","",300,0,3);
	TH1F *h1_prot_mom_ratio = new TH1F("h1_prot_mom_ratio","",50,0.97,1.2);
	TH1F *h1_Wvar_weight = new TH1F("h1_Wvar_weight","",400,0,3);
	TH1F *h1_xbjk_weight = new TH1F("h1_xbjk_weight","",400,0,3);
	TH1F *h1_Q2_weight = new TH1F("h1_Q2_weight","",400,0,6);
	TH1F *h1_nu_weight = new TH1F("h1_nu_weight","",400,0,4);
	TH1F *h1_WvarCal_weight = new TH1F("h1_WvarCal_weight","",400,0,3);
	TH1F *h1_xbjkCal_weight = new TH1F("h1_xbjkCal_weight","",400,0,3);
	TH1F *h1_Q2Cal_weight = new TH1F("h1_Q2Cal_weight","",400,0,6);
	TH1F *h1_nuCal_weight = new TH1F("h1_nuCal_weight","",400,0,3);

	// -------------------------------------------------------------------------------------------------------

	//Binning for energy reconstruction histograms
	int n_bins;
	double *x_values;
	double *x_qe;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=17;i++) { x_values[i]=0.4+i*0.04; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=20;i++) { x_values[i+18]=1.08+(i+1)*0.02; x_qe[i+18] = (x_values[i+18] - en_beam[fbeam_en]) / en_beam[fbeam_en]; }
	}

	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
		n_bins=54;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=23;i++) { x_values[i]=i*0.09; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=30;i++) { x_values[i+24]=2.07+(i+1)*0.03; x_qe[i+24] = (x_values[i+24] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=21;i++)	{ x_values[i]=i*0.2; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=16;i++)	{ x_values[i+22]=4.2+(i+1)*0.05; x_qe[i+22] = (x_values[i+22] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	// -------------------------------------------------------------------------------------------------------

	//Definitions of further Histograms

	TH1F *CosDeltaThetaElectronPhotonAboveThreshold=new TH1F("CosDeltaThetaElectronPhotonAboveThreshold","",100,-1.,1.);
	TH1F *CosDeltaPhiElectronPhotonAboveThreshold=new TH1F("CosDeltaPhiElectronPhotonAboveThreshold","",100,-1.,1.);

	TH2F *RadCosThetaGammaEgamma = new TH2F("RadCosThetaGammaEgamma","",100,-1.,1.,600,0.,6.);
	TH2F *RadCosDeltaThetaGammaEgamma = new TH2F("RadCosDeltaThetaGammaEgamma","",100,-1.,1.,600,0.,6.);
	TH2F *NonRadThetaVsPhiGamma = new TH2F("NonRadThetaVsPhiGamma","",360,0.,360.,180,0.,180.);

	TH1F *h1_E_rec_1pi_weight_frac_feed=new TH1F("h1_E_rec_1pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_2pi_weight_frac_feed=new TH1F("h1_E_rec_2pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_3pi_weight_frac_feed=new TH1F("h1_E_rec_3pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pi_weight_frac_feed=new TH1F("h1_E_rec_4pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_0pi_frac_feed=new TH1F("h1_E_rec_0pi_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_tot_cut2_fracfeed = new TH1F("h1_E_tot_cut2_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_cut2_new_fracfeed = new TH1F("h1_E_rec_cut2_new_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_p_bkgd_fracfeed = new TH1F("h1_E_tot_p_bkgd_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_p_bkgd_fracfeed = new TH1F("h1_E_rec_p_bkgd_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_2p0pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_2p0pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p1pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p1pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p0pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p0pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_3pto2p_fracfeed = new TH1F("h1_E_tot_3pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_3pto2p_fracfeed = new TH1F("h1_E_rec_3pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_3pto1p_fracfeed = new TH1F("h1_E_tot_3pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_3pto1p_fracfeed = new TH1F("h1_E_rec_3pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_4pto3p_fracfeed = new TH1F("h1_E_tot_4pto3p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pto3p_fracfeed = new TH1F("h1_E_rec_4pto3p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_43pto1p_fracfeed = new TH1F("h1_E_tot_43pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_43pto1p_fracfeed = new TH1F("h1_E_rec_43pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_4pto2p_fracfeed = new TH1F("h1_E_tot_4pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pto2p_fracfeed = new TH1F("h1_E_rec_4pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_4pto1p_fracfeed = new TH1F("h1_E_tot_4pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pto1p_fracfeed = new TH1F("h1_E_rec_4pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_1p2pi_fracfeed = new TH1F("h1_E_tot_1p2pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_1p2pi_fracfeed = new TH1F("h1_E_rec_1p2pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_1p3pi_fracfeed = new TH1F("h1_E_tot_1p3pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_1p3pi_fracfeed = new TH1F("h1_E_rec_1p3pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_2p2pi_fracfeed = new TH1F("h1_E_tot_2p2pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_2p2pi_fracfeed = new TH1F("h1_E_rec_2p2pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_3p1pi_fracfeed = new TH1F("h1_E_tot_3p1pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_3p1pi_fracfeed = new TH1F("h1_E_rec_3p1pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_tot_1p2pi_1p0pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_rec_1p2pi_1p0pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_undetfactor_fracfeed = new TH1F("h1_E_rec_undetfactor_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_undetfactor_fracfeed = new TH1F("h1_E_tot_undetfactor_fracfeed","",n_bins,x_qe);

	TH1F *h1_theta0=new TH1F("h1_theta0","",300,0,180);
	TH2F *h2_Ecal_Eqe=new TH2F("h2_Ecal_Eqe","",800,0,8.,800,0,8.);
	TH1F *h1_Ecal=new TH1F("h1_Ecal","",100,0,10.);
	TH1F *h1_Ecal_Reso=new TH1F("h1_Ecal_Reso","",100,-1.,1.);
	TH2F *h2_Ecal_Etrue=new TH2F("h2_Ecal_Etrue","",100,0,10.,100,0,10.);
	TH2F *h2_Etrue_Ecal=new TH2F("h2_Etrue_Ecal","",100,0,10.,100,0,10.);
	TH2F *h2_EqeEcalratio_Eqe=new TH2F("h2_EqeEcalratio_Eqe","",600,0,5,300,0,2);
	TH2F *h2_EqeEcaldiff_Eqe=new TH2F("h2_EqeEcaldiff_Eqe","",600,0,5,300,-3,3);
	TH2F *h2_N_prot_pi=new TH2F("h2_N_prot_pi","",10,0,5,10,0,5);
	TH2F *h2_N_prot_pi_phot=new TH2F("h2_N_prot_pi_phot","",10,0,5,10,0,5);
	TH2F *h2_N_prot_pi_phot_nonrad=new TH2F("h2_N_prot_pi_phot_nonrad","",10,0,5,10,0,5);
//	TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,0,180);
	TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,10,60);
	TH2F *h2_el_mom_diff = new TH2F("h2_el_mom_diff","",500,0.,1.,500,-0.1,0.1);

	int NBinsNu = 300, NBinsQ2 = 300;
	double MinNu = 0., MaxNu = 4.; double MinQ2 = 0., MaxQ2 = 6.;
	TH2F *h2_Q2_nu = new TH2F("h2_Q2_nu","",NBinsNu,MinNu,MaxNu,NBinsQ2,MinQ2,MaxQ2);
	TH2F *h2_Q2_nu_weight = new TH2F("h2_Q2_nu_weight","",NBinsNu,MinNu,MaxNu,NBinsQ2,MinQ2,MaxQ2);
	TH2F *h2_Q2_nu_weight_FirstSector = new TH2F("h2_Q2_nu_weight_FirstSector","",0.7*NBinsNu,MinNu,MaxNu,0.7*NBinsQ2,MinQ2,MaxQ2);

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
	TH2F *h2_Etot_pperp = new TH2F("h2_Etot_pperp","",400,0,1,400,0,6.);

	TH2F *h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);

	TH2F* h2_QVector_theta_phi = new TH2F("h2_QVector_theta_phi","",200,0,360,200,0,80);

//	//Binning for energy reconstruction histograms
//	int n_bins;
//	double *x_values;

//	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
//		n_bins=38;
//		x_values=new double[n_bins+1];
//		for (int i=0;i<=17;i++) x_values[i]=0.4+i*0.04;
//		for (int i=0;i<=20;i++) x_values[i+18]=1.08+(i+1)*0.02;
//	}

//	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
//		n_bins=54;
//		x_values=new double[n_bins+1];
//		for (int i=0;i<=23;i++) x_values[i]=i*0.09;
//		for (int i=0;i<=30;i++) x_values[i+24]=2.07+(i+1)*0.03;
//	}

//	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
//		n_bins=38;
//		x_values=new double[n_bins+1];
//		for (int i=0;i<=21;i++)	x_values[i]=i*0.2;
//		for (int i=0;i<=16;i++)	x_values[i+22]=4.2+(i+1)*0.05;
//	}

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
	TH1F *h1_E_rec_1prot	= new TH1F("h1_E_rec_1prot","",n_bins,x_values);
	TH1F *h1_E_tot_1prot	= new TH1F("h1_E_tot_1prot","",n_bins,x_values);
	TH1F *h1_E_rec_cutpi1_piplpimi = new TH1F("h1_E_rec_cutpi1_piplpimi","",n_bins,x_values);
	TH1F *h1_E_tot_cutpi1_piplpimi = new TH1F("h1_E_tot_cutpi1_piplpimi","",n_bins,x_values);
	TH1F *h1_E_tot = new TH1F("h1_E_tot","",n_bins,x_values);
	TH1F *h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
	TH1F *h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
	TH1F *h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
	TH1F *h1_E_rec_undetfactor	= new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
	TH1F *h1_E_tot_undetfactor	= new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
	TH1F *h1_E_tot_1p2pi	= new TH1F("h1_E_tot_1p2pi","",n_bins,x_values);
	TH1F *h1_E_rec_1p2pi	= new TH1F("h1_E_rec_1p2pi","",n_bins,x_values);
	TH1F *h1_E_tot_1p3pi	= new TH1F("h1_E_tot_1p3pi","",n_bins,x_values);
	TH1F *h1_E_rec_1p3pi	= new TH1F("h1_E_rec_1p3pi","",n_bins,x_values);
	TH1F *h1_E_tot_2p2pi	= new TH1F("h1_E_tot_2p2pi","",n_bins,x_values);
	TH1F *h1_E_rec_2p2pi	= new TH1F("h1_E_rec_2p2pi","",n_bins,x_values);
	TH1F *h1_E_tot_3p1pi	= new TH1F("h1_E_tot_3p1pi","",n_bins,x_values);
	TH1F *h1_E_rec_3p1pi	= new TH1F("h1_E_rec_3p1pi","",n_bins,x_values);
	TH1F *h1_E_tot_1p2pi_1p0pi	= new TH1F("h1_E_tot_1p2pi_1p0pi","",n_bins,x_values);
	TH1F *h1_E_rec_1p2pi_1p0pi	= new TH1F("h1_E_rec_1p2pi_1p0pi","",n_bins,x_values);
	TH1F *h1_E_tot_2p1pi_2p0pi	= new TH1F("h1_E_tot_2p1pi_2p0pi","",n_bins,x_values);
	TH1F *h1_E_rec_2p1pi_2p0pi	= new TH1F("h1_E_rec_2p1pi_2p0pi","",n_bins,x_values);
	TH1F *h1_E_tot_2p1pi_1p1pi	= new TH1F("h1_E_tot_2p1pi_1p1pi","",n_bins,x_values);
	TH1F *h1_E_rec_2p1pi_1p1pi	= new TH1F("h1_E_rec_2p1pi_1p1pi","",n_bins,x_values);
	TH1F *h1_E_tot_2p1pi_1p0pi	= new TH1F("h1_E_tot_2p1pi_1p0pi","",n_bins,x_values);
	TH1F *h1_E_rec_2p1pi_1p0pi	= new TH1F("h1_E_rec_2p1pi_1p0pi","",n_bins,x_values);

	// Unweighted plots for the number of events

	TH1F *h1_MissMomentum_NoWeight = new TH1F("MissMomentum_NoWeight","",100,0.,1.);

	TH1F *h1_ECal_Slice0_NoWeight = new TH1F("epRecoEnergy_slice_0_NoWeight","",n_bins,x_values);
	TH1F *h1_ECal_Slice1_NoWeight = new TH1F("epRecoEnergy_slice_1_NoWeight","",n_bins,x_values);
	TH1F *h1_ECal_Slice2_NoWeight = new TH1F("epRecoEnergy_slice_2_NoWeight","",n_bins,x_values);
	TH1F *h1_ECal_Slice3_NoWeight = new TH1F("epRecoEnergy_slice_3_NoWeight","",n_bins,x_values);

	TH1F *h1_EQE_Slice0_NoWeight = new TH1F("eRecoEnergy_slice_0_NoWeight","",n_bins,x_values);
	TH1F *h1_EQE_Slice1_NoWeight = new TH1F("eRecoEnergy_slice_1_NoWeight","",n_bins,x_values);
	TH1F *h1_EQE_Slice2_NoWeight = new TH1F("eRecoEnergy_slice_2_NoWeight","",n_bins,x_values);
	TH1F *h1_EQE_Slice3_NoWeight = new TH1F("eRecoEnergy_slice_3_NoWeight","",n_bins,x_values);

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

	// Plots for interaction break down for GENIE samples
	const int NInt = 6; // All Interactions = 0, QE = 1, MEC = 2, RES = 3, DIS = 4, Other = 6
	TH1D* ECal_BreakDown[NInt];
	TH1D* EQE_BreakDown[NInt];
	TH1D* InclusiveEQE_BreakDown[NInt];
	TH1D* Pmiss_BreakDown[NInt];
	TH1D* Q2_BreakDown[NInt];
	TH1D* Nu_BreakDown[NInt];
	TH1D* Pe_BreakDown[NInt];

	for (int WhichInt = 0; WhichInt < NInt; WhichInt++) {

		ECal_BreakDown[WhichInt] = new TH1D(Form("ECal_Int_%d",WhichInt),";E^{Cal} (GeV)",n_bins,x_values);
		EQE_BreakDown[WhichInt] = new TH1D(Form("EQE_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);
		InclusiveEQE_BreakDown[WhichInt] = new TH1D(Form("InclusiveEQE_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);
		Pmiss_BreakDown[WhichInt] = new TH1D(Form("Pmiss_Int_%d",WhichInt),";P_{miss}^{#perp} [GeV/c]",80,0.,1.);
		Q2_BreakDown[WhichInt] = new TH1D(Form("Q2_Int_%d",WhichInt),";Q^{2} [GeV^{2}/c^{2}]",400,0,6);
		Nu_BreakDown[WhichInt] = new TH1D(Form("Nu_Int_%d",WhichInt),";Energy Transfer [GeV]",400,0,4);
		Pe_BreakDown[WhichInt] = new TH1D(Form("Pe_Int_%d",WhichInt),";P_{e} [GeV/c]",100,0.,5.);
	}

	// Vector containing kinematic variables using Ecal
	vector<double> CalKineVars{};
	// Weight to fill the plots mentioned above
	double LocalWeight;

	// Signal Event Counter -> 1e1p0pi events (everything lese is bkg)
	int SignalEvents = 0;
	int QESignalEvents = 0;
	int MECSignalEvents = 0;
	int RESSignalEvents = 0;
	int DISSignalEvents = 0;
	int OtherSignalEvents = 0;

	int EQESignalEventsWithin5Perc = 0, EQESignalEventsWithin5Perc_FirstSlice = 0, EQESignalEventsWithin5Perc_SecondSlice = 0, EQESignalEventsWithin5Perc_ThirdSlice = 0;
	int ECalSignalEventsWithin5Perc = 0, ECalSignalEventsWithin5Perc_FirstSlice = 0, ECalSignalEventsWithin5Perc_SecondSlice = 0, ECalSignalEventsWithin5Perc_ThirdSlice = 0;
	int PMiss_FirstSlice = 0, PMiss_SecondSlice = 0, PMiss_ThirdSlice = 0;

	// ---------------------------------------------------------------------------------------------------------------

	// Get the number of events to run overall
	
//	int Nentries = TMath::Min(Ntfileentries,NtweightsEntries);

	// ---------------------------------------------------------------------------------------------------------------

	/** Beginning of Event Loop **/
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//	for (Long64_t jentry=0; jentry<Nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		//Read Entry
		int nb = GetEntry(jentry);
		if (nb == 0) { std::cout <<"Event loop: 0 byte read for entry " << jentry << ". Indicate failure in reading the file" <<	std::endl;}

		if (jentry%1000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/fChain->GetEntries()*100. << " %"<< std::endl;}

		if( jentry%200000 == 0 )
		{
			gDirectory->Write("hist_Files", TObject::kOverwrite);
			//cout<<jentry<<endl;
		}

		if(jentry == 0){ //first entry to initialize TorusCurrent, Fiducials and Subtraction classes

			//The TorusField has to be set before the Fiducialcut parameters are initialized
			if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. ) //1.1 GeV, we are not using the 1.1 GeV data with 1500 current field
			{
				 fTorusCurrent = 750;
			}
			else if( (en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) || (en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) ) //2.2 GeV	or 4.4 GeV
			{
				 fTorusCurrent = 2250;
			}
			else { std::cout << "systematics::Loop(): fTorusCurrent could not be assigned" << std::endl;}

			fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
			fiducialcut->SetFiducialCutParameters(fbeam_en);
			std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
			rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
			std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;
		}

		//Resets q vector to (0,0,0)
		rotation->ResetQVector();

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		double SmearedPe;
		double SmearedEe;
		double e_acc_ratio = 1.;	//will be 1 for CLAS data

		// Outgoing e',	Uncorr and corrected are the same read from root file.
		//V4_el and V3_el will be changed by smearing for GENIE simulation data
		TLorentzVector V4_el(pxl,pyl,pzl,El);
		TLorentzVector V4_el_uncorr(pxl,pyl,pzl,El);
		TVector3 V3_el(pxl,pyl,pzl);

		double el_momentum = V3_el.Mag();
		double el_theta = V3_el.Theta();

		if (choice == 1) { //smearing, fiducials and acceptance ratio for GENIE simulation data

			//Smearing of Electron Vector from Simulation
			SmearedPe = gRandom->Gaus(pl,reso_e*pl);
			SmearedEe = sqrt( SmearedPe*SmearedPe + e_mass * e_mass );
			V3_el.SetXYZ(SmearedPe/pl * pxl,SmearedPe/pl * pyl,SmearedPe/pl * pzl);
			V4_el.SetPxPyPzE(V3_el.X(),V3_el.Y(),V3_el.Z(),SmearedEe);
			double phi_ElectronOut = V3_el.Phi(); //in Radians

			V3_el.SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180), GENIE coordinate system flipped with respect to CLAS

			//Fiducial Cuts with the smeared values
			if (!EFiducialCut(fbeam_en,V3_el) ) continue; // Electron theta & phi fiducial cuts

			phi_ElectronOut += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
			el_momentum = V3_el.Mag(); //Momentum after smearing
			el_theta = V3_el.Theta(); //Angle after smearing

			//acceptance_c takes phi in radians and here unmodified by 30 degree.
			e_acc_ratio = acceptance_c(el_momentum, cos(el_theta), phi_ElectronOut, 11,file_acceptance);
			if ( fabs(e_acc_ratio) != e_acc_ratio ) { continue; }

			// --------------------------------------------------------------------------------------------------

			// GENIE Systematic Uncertainties

			tweights->GetEntry(jentry);
			float* ArrayWeights = weights->GetArray();
			double TuningWeight = ArrayWeights[sigma];
			e_acc_ratio = e_acc_ratio * TuningWeight;

			// --------------------------------------------------------------------------------------------------		

		}

		// Explicit cuts on electron momentum
		if (fbeam_en=="1161" && el_momentum < 0.4) { continue; }
		if (fbeam_en=="2261" && el_momentum < 0.55) { continue; }
		if (fbeam_en=="4461" && el_momentum < 1.1) { continue; }

		//Definition as for data. It is also correct for GENIE simulation data since V3_el is rotated above by 180 degree in phi
		double el_phi_mod = V3_el.Phi()*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
		if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree


		//Calculated Mott Cross Section and Weights for Inclusive Histograms
		//Wght and e_acc_ratio is 1 for CLAS data
		//double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(el_theta)+1))/(2*pow(El,2.)*pow((1-cos(el_theta)),2.));

		double reco_Q2 = -(V4_el-V4_beam).Mag2();
		double Q4 = reco_Q2 * reco_Q2;
		double Mott_cross_sec = (1./Q4) * XSecScale;

		// ---------------------------------------------------------------------------------------------------------------------

		// For neutrino scattering 
		// switch to true for nu scattering to account for the difference in the propagator

//		bool neutrino = false;

//		if (neutrino) { 

//			XSecScale = gr->Eval(Ev);
//			Mott_cross_sec = XSecScale;

//		}

		// ---------------------------------------------------------------------------------------------------------------------

		double WeightIncl = wght*e_acc_ratio / Mott_cross_sec;

		// Securing ourselves against infinities
		if ( fabs(WeightIncl) != WeightIncl ) { continue; }

		//Calculation of Reconstructed Energy from ELectron only
		//using the same value of single nucleon separation E Ecal and Eqe
		double E_rec = (m_prot*bind_en[ftarget]+m_prot*V4_el.E())/(m_prot-V4_el.E()+V4_el.Rho()*cos(el_theta));

		//Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
		double nu = -(V4_el-V4_beam).E();
		double x_bjk = reco_Q2/(2*m_prot*nu);

		// QE selection
		//if ( fabs(x_bjk - 1.) > 0.2) { continue; }

		// ---------------------------------------------------------------------------------------------------------------------

		TVector3 V3_q = (V4_beam-V4_el).Vect();
		double V3_q_theta_deg = V3_q.Theta() * 180. / TMath::Pi();
		double V3_q_phi_deg = V3_q.Phi() * 180. / TMath::Pi() + 30.; 
		if (V3_q_phi_deg > 360) { V3_q_phi_deg = V3_q_phi_deg - 360.; } 
		if (V3_q_phi_deg < 0) { V3_q_phi_deg = V3_q_phi_deg + 360.; }
		double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

		//converting theta to degrees
		el_theta = el_theta*TMath::RadToDeg();

		//Cuts on Q2 and W, only keep events with Q2 > Q2cut and W < Wcut
		if ( reco_Q2 < Q2cut || W_var > Wcut) continue;

		//Set q vector for the following rotations for the subtraction procedure
		rotation->SetQVector(V3_q);
//		rotation->PrintQVector();

		h1_el_mom->Fill(V4_el_uncorr.Rho(),WeightIncl);
		h1_el_mom_ratio->Fill(V4_el.Rho()/V4_el_uncorr.Rho(),WeightIncl);
		h2_el_pcorr_puncorr->Fill(V4_el.Rho(),V4_el.Rho()/V4_el_uncorr.Rho(),WeightIncl);
		h2_el_mom_diff->Fill(V4_el.Rho(),V4_el.Rho()-V4_el_uncorr.Rho(),WeightIncl);

		//Filling Histogram for electron kinematics
		h1_xbjk->Fill(x_bjk);
		h1_Q2->Fill(reco_Q2);
		h1_Wvar->Fill(W_var);
		h1_el_Mott_crosssec->Fill(Mott_cross_sec);
		h2_el_theta_phi->Fill(el_phi_mod,el_theta,WeightIncl);
		h1_el_theta->Fill(el_theta);
		h2_Q2_nu->Fill(nu,reco_Q2);
		h2_Q2_xbjk_weight->Fill(x_bjk,reco_Q2,WeightIncl);
		h2_Q2_W->Fill(W_var,reco_Q2);
		h2_xB_W->Fill(W_var,x_bjk);
		h2_Q2_W_weight->Fill(W_var,reco_Q2,WeightIncl);

		//Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

		//Index variables for hadrons (p and pions)
		int index_p[20]; //index for each proton
		int index_pi[20]; //index for each pion
		int ind_pi_phot[20]; //index for pions and photons
		int index_pipl[20]; //index for each pi plus
		int index_pimi[20]; //index for each pi minus

		int charge_pi[20]; //Charge for the pions and photons
		//Smeared Momentum and Energy values for GENIE (simulation) data
		double Smeared_Pp[20]; //smeared momentum values for protons
		double Smeared_Ep[20]; //smeared energy values for protons
		double Smeared_Ppi[20]; //smeared momentum values for pions
		double Smeared_Epi[20]; //smeared energy values for pions

		//Number of hadrons
		int num_p = 0;
		int num_pi = 0;
		int num_pi_phot = 0; //couting all pions and photons
		int num_pimi = 0;
		int num_pipl = 0;
		int num_pi_phot_nonrad = 0; //counting all pions and non-radiation photons
		int num_phot_rad = 0; //counting radiation photons
		//Index and number variables for neutral particles
		int ec_index_n[20];
		int ec_num_n = 0;
		bool ec_radstat_n[20];

		//Array initialize to -1 or false
		for (int i = 0; i < 20; i++) {
			index_p[i] = -1;   index_pi[i] = -1;   index_pipl[i] = -1;   index_pimi[i] = -1;   ind_pi_phot[i] = -1;
			ec_index_n[i] = -1;   ec_radstat_n[i] = false;
			charge_pi[i] = -2; //default number should be not a possible real charge
			Smeared_Pp[i]  = 0; Smeared_Ep[i]  = 0;  //default 0 momentum and energy after smearing
			Smeared_Ppi[i] = 0; Smeared_Epi[i] = 0;  //default 0 momentum and energy after smearing
		}

		const double phot_rad_cut = 40;
		const double phot_e_phidiffcut = 30; //electron - photon phi difference cut

		// Creating vectors to store id of particles in the array
		vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
		ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

		//Loop for Hadrons
		for (int i = 0; i < nf; i++) {

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			//Start of proton selection
			if (pdgf[i] == 2212  && pf[i] > 0.3) {

				if ( choice == 1) { //GENIE data

					//Smearing of proton
					double temp_smear_P = gRandom->Gaus(pf[i],reso_p*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_prot * m_prot );

					TVector3 V3_prot_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_prot = V3_prot_corr.Phi();
					V3_prot_corr.SetPhi(phi_prot + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS
					if (!PFiducialCut(fbeam_en, V3_prot_corr) ) { continue; } // Proton theta & phi fiducial cuts

					num_p = num_p + 1;
					index_p[num_p - 1] = i;
					ProtonID.push_back(i);
					Smeared_Pp[num_p - 1] = temp_smear_P;
					Smeared_Ep[num_p - 1] = temp_smear_E;
				}
				else { //CLAS data does not need Fiducial Cut again

						num_p = num_p + 1;
						index_p[num_p - 1] = i;
						ProtonID.push_back(i);
				}
			}

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if (pdgf[i] == -211  && pf[i] > 0.15)  { //PI minus

				if ( choice == 1) { //GENIE data

					//Smearing of pi minus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					// Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
					if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, -1) )     {  continue; }

					num_pimi = num_pimi + 1;
					num_pi = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pimi[num_pi_phot - 1] = i;
					index_pi[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiMinusID.push_back(i);
					charge_pi[num_pi_phot - 1] = -1;
					Smeared_Ppi[num_pi_phot - 1] = temp_smear_P;
					Smeared_Epi[num_pi_phot - 1] = temp_smear_E;
				}
				else { //CLAS data does not need Fiducial Cut again
					num_pimi = num_pimi + 1;
					num_pi = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pimi[num_pi_phot - 1] = i;
					index_pi[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiMinusID.push_back(i);
					charge_pi[num_pi_phot - 1] = -1;
				}
			}

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if ( pdgf[i] == 211  && pf[i] > 0.15)  {

				if ( choice == 1) { //GENIE data
					//Smearing of pi plus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					// Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
					if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, 1) )     {  continue; }

					num_pipl = num_pipl + 1;
					num_pi  = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pipl[num_pi_phot - 1] = i;
					index_pi[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiPlusID.push_back(i);
					charge_pi[num_pi_phot - 1] = 1;
					Smeared_Ppi[num_pi_phot - 1] = temp_smear_P;
					Smeared_Epi[num_pi_phot - 1] = temp_smear_E;
				}
				else { //CLAS data does not need Fiducial Cut again
					num_pipl = num_pipl + 1;
					num_pi  = num_pi + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pipl[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiPlusID.push_back(i);
					charge_pi[num_pi_phot - 1] = 1;
				}
			}

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if (pdgf[i] == 22  && pf[i] > 0.3) {

				//Determine photon vector for the cut on radiation photon via angle with respect to the electron
				TVector3 V3_phot_angles(pxf[i],pyf[i],pzf[i]);
				if (choice == 1) { //GENIE data
					//no smearing of GENIE photons
					double phi_photon = V3_phot_angles.Phi();
					V3_phot_angles.SetPhi(phi_photon + TMath::Pi()); // Vec.Phi() is between (-180,180)
					if ( !Pi_phot_fid_united(fbeam_en, V3_phot_angles, 0) )  { continue;}
				}

				double neut_phi_mod = V3_phot_angles.Phi()*TMath::RadToDeg() + 30; //Add 30 degree
				if (neut_phi_mod < 0) neut_phi_mod = neut_phi_mod + 360;  //Neutral particle is between 0 and 360 degree

				ec_num_n = ec_num_n + 1;
				num_pi_phot = num_pi_phot + 1;
				ind_pi_phot[num_pi_phot - 1] = i;
				PhotonID.push_back(i);

				Smeared_Ppi[num_pi_phot - 1] = V3_phot_angles.Mag();
				Smeared_Epi[num_pi_phot - 1] = V3_phot_angles.Mag();

				CosDeltaThetaElectronPhotonAboveThreshold->Fill( cos( V3_phot_angles.Angle(V3_el) ) );
				CosDeltaPhiElectronPhotonAboveThreshold->Fill( cos( neut_phi_mod-el_phi_mod*TMath::Pi()/180. ) );

				 //within 40 degrees in theta and 30 degrees in phi. Electron phi has already added 30 degree and between 0 to 360

				 if(V3_phot_angles.Angle(V3_el)*TMath::RadToDeg() < phot_rad_cut && fabs(neut_phi_mod-el_phi_mod) < phot_e_phidiffcut ) {

					ec_radstat_n[num_pi_phot - 1] = true; //select radiation photons
					num_phot_rad = num_phot_rad + 1;
					RadCosThetaGammaEgamma->Fill(V3_phot_angles.CosTheta(),V3_phot_angles.Mag() ,WeightIncl);
					RadCosDeltaThetaGammaEgamma->Fill( cos( V3_phot_angles.Angle(V3_el) ) ,V3_phot_angles.Mag() ,WeightIncl);

				 }

				 if(!ec_radstat_n[num_pi_phot - 1]) {
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					charge_pi[num_pi_phot - 1] = 0;
					NonRadThetaVsPhiGamma->Fill(neut_phi_mod,V3_phot_angles.Theta()*TMath::RadToDeg(),WeightIncl);
				 }
			}

		} //end of hadron loop

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------

		//Skip event if there is at least one radiation photon
		if (num_phot_rad > 0) {
		  continue;
		}

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------

		// For GENIE samples, identify the interaction type

		int Interaction = -1;
		if (choice == 1) {

			if (qel) { Interaction = 1; }
			if (mec) { Interaction = 2; }
			if (res) { Interaction = 3; }
			if (dis) { Interaction = 4; }

		}

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------

		//Filling Histograms with multiplicities
		h1_Npi->Fill(num_pi);
		h1_Nprot->Fill(num_p);

		if (num_p > 0) {
			h1_Nprot_NonZeroProt->Fill(num_p);
			h1_Npi_NonZeroProt->Fill(num_pi);
			h2_QVector_theta_phi->Fill(V3_q_phi_deg,V3_q_theta_deg,WeightIncl);
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

		//Events with exactly 2 protons
		if(num_p == 2) {

			//LorentzVectors for protons without momentum smearing or corrections
			TLorentzVector V4_prot_uncorr1(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
			TLorentzVector V4_prot_uncorr2(pxf[index_p[1]],pyf[index_p[1]],pzf[index_p[1]],TMath::Sqrt(m_prot*m_prot+pf[index_p[1]]*pf[index_p[1]]));
			//LorentzVectors for protons with momentum smearing or corrections
			TVector3 V3_prot_corr1;
			TVector3 V3_prot_corr2;

			double p_acc_ratio1 = 1; //will be 1 for CLAS data
			double p_acc_ratio2 = 1; //will be 1 for CLAS data

			if (choice == 0) { //CLAS data
				V3_prot_corr1.SetXYZ(pxf[index_p[0]+60],pyf[index_p[0]+60],pzf[index_p[0]+60]);
				V3_prot_corr2.SetXYZ(pxf[index_p[1]+60],pyf[index_p[1]+60],pzf[index_p[1]+60]);
			}

			if (choice == 1) { //GENIE data, fiducials are done in hadron loop

				V3_prot_corr1.SetXYZ(Smeared_Pp[0]/pf[index_p[0]] * pxf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pyf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pzf[index_p[0]]);
				double phi_prot1 = V3_prot_corr1.Phi();
				V3_prot_corr1.SetPhi(phi_prot1 + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS
				phi_prot1 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				double p_theta1 = V3_prot_corr1.Theta();
				double prot_mom_corr1 = V3_prot_corr1.Mag();
				//Proton 1 weight
				p_acc_ratio1 = acceptance_c(prot_mom_corr1, cos(p_theta1), phi_prot1, 2212,file_acceptance_p);
				if ( fabs(p_acc_ratio1) != p_acc_ratio1 ) { continue; }

				V3_prot_corr2.SetXYZ(Smeared_Pp[1]/pf[index_p[1]] * pxf[index_p[1]],Smeared_Pp[1]/pf[index_p[1]] * pyf[index_p[1]],Smeared_Pp[1]/pf[index_p[1]] * pzf[index_p[1]]);
				double phi_prot2 = V3_prot_corr2.Phi();
 				V3_prot_corr2.SetPhi(phi_prot2 + TMath::Pi()); // Vec.Phi() is between (-180,180) // GENIE coordinate system flipped with respect to CLAS
				phi_prot2 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				double p_theta2 = V3_prot_corr2.Theta();
				double prot_mom_corr2 = V3_prot_corr2.Mag();
				//Proton 2 weight
				p_acc_ratio2 = acceptance_c(prot_mom_corr2, cos(p_theta2), phi_prot2, 2212,file_acceptance_p);
				if ( fabs(p_acc_ratio2) != p_acc_ratio2 ) { continue; }

			}

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
			double N_prot_both = 0;
			double P_N_2p[2]={0};

			rotation->prot2_rot_func( V3_2prot_corr, V3_2prot_uncorr, V4_el, E_tot_2p, p_perp_tot_2p, P_N_2p , &N_prot_both);

			if(num_pi_phot==0 && N_prot_both!=0){

				double histoweight = weight_protons*e_acc_ratio*wght/Mott_cross_sec; //total weight from 2p acceptance , 1e acceptance, Mott, and GENIE weight

				for(int f = 0; f < num_p; f++){    //looping through two protons

					h1_E_tot_p_bkgd->Fill(E_tot_2p[f],P_N_2p[f]*histoweight);
					h1_E_rec_p_bkgd->Fill(E_rec,P_N_2p[f]*histoweight);
					h2_Erec_pperp_2p->Fill(p_perp_tot_2p[f],E_rec,P_N_2p[f]*histoweight);
					h2_Etot_pperp->Fill(p_perp_tot_2p[f],E_tot_2p[f],-P_N_2p[f]*histoweight);
					h1_E_tot_p_bkgd_fracfeed->Fill((E_tot_2p[f]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_N_2p[f]*histoweight);
					h1_E_rec_p_bkgd_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_N_2p[f]*histoweight);
					h2_pperp_W->Fill(W_var,p_perp_tot_2p[f],-P_N_2p[f]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[f]) *TMath::RadToDeg(),-P_N_2p[f]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[f],-P_N_2p[f]*histoweight);
					h1_Ecal->Fill(E_tot_2p[f],-P_N_2p[f]*histoweight);
					h1_Ecal_Reso->Fill((E_tot_2p[f]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_N_2p[f]*histoweight);
					h2_Ecal_Etrue->Fill(E_tot_2p[f],Ev,-P_N_2p[f]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_tot_2p[f],-P_N_2p[f]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[f],-P_N_2p[f]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[f],-P_N_2p[f]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,-P_N_2p[f]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,-P_N_2p[f]*histoweight);
					h1_Wvar_weight->Fill(W_var,-P_N_2p[f]*histoweight);
					h1_nu_weight->Fill(nu,-P_N_2p[f]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),-P_N_2p[f]*histoweight);
					h1_prot_mom->Fill(V3_2prot_corr[f].Mag(),-P_N_2p[f]*histoweight);
					h1_MissMomentum->Fill(p_perp_tot_2p[f],-P_N_2p[f]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot_2p[f],V4_el);
					LocalWeight = -P_N_2p[f]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_tot_2p[f],LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_perp_tot_2p[f],LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

 					if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_tot_2p[f],LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_perp_tot_2p[f],LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0 ; i < n_slice; i++) {

						if (p_perp_tot_2p[f] < pperp_max[i] && p_perp_tot_2p[f] > pperp_min[i]) {
							h1_Etot_p_bkgd_slice[i]->Fill(E_tot_2p[f],P_N_2p[f]*histoweight);
							h1_Erec_p_bkgd_slice[i]->Fill(E_rec,P_N_2p[f]*histoweight);
						}

					}

				} //looping through two protons

				h1_E_tot_2p_det->Fill(E_tot_2p[0],histoweight);
				h1_E_rec_2p_det->Fill(E_rec,histoweight);

			}//no pions cut and N_prot_both!=0

			//---------------------------------- 2p 1pi   ----------------------------------------------
			//Const int can be placed somewhere up after if for 2 protons F.H. 05.09.19
			const int N_2prot=2;
			//Variable might/could be placed in a more local context F.H. 05.09.19
			double Ecal_2p1pi_to2p0pi[N_2prot]={0};
			double p_miss_perp_2p1pi_to2p0pi[N_2prot]={0};

			if (num_pi_phot==1) {

				TVector3 V3_1pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 0) { //CLAS data
					V3_1pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);
				}

				if (choice == 1) { //GENIE data
					pion_acc_ratio = 0;//reset to 0 just to be save
					V3_1pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);

					double phi_pion = V3_1pi_corr.Phi();
					V3_1pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_1pi_corr.Theta();
					double pion_mom_corr = V3_1pi_corr.Mag();

					if (charge_pi[0] == 1) { //acceptance for pi plus
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == -1) {    //acceptance for pi minus. using electron acceptance map
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == 0) {    //acceptance for neutral, setting to 1 for now F.H. 09/24/19
						pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 2proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }
				}

				double P_2p1pito2p0pi[2] = {0};
				double P_2p1pito1p1pi[2] = {0};
				double P_2p1pito1p0pi[2] = {0};
				double Ptot = 0;

				rotation->prot2_pi1_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_1pi_corr, charge_pi[0], V4_el,Ecal_2p1pi_to2p0pi,p_miss_perp_2p1pi_to2p0pi,P_2p1pito2p0pi, P_2p1pito1p1pi, P_2p1pito1p0pi,&Ptot);

				double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Is this correct in the following loop? F.H. 09/01/19

				for(int z=0; z < N_2prot; z++){ //looping over two protons

					//---------------------------------- 2p 1pi ->2p 0pi ----------------------------------------------

					h1_E_tot_2p1pi_2p0pi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_E_rec_2p1pi_2p0pi->Fill(E_rec,P_2p1pito2p0pi[z]*histoweight);
					h2_Erec_pperp_2p1pi_2p0pi->Fill(p_miss_perp_2p1pi_to2p0pi[z],E_rec,P_2p1pito2p0pi[z]*histoweight);
					h2_Etot_pperp->Fill(p_miss_perp_2p1pi_to2p0pi[z],Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_E_tot_2p1pi_2p0pi_fracfeed->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
					h1_E_rec_2p1pi_2p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
					h2_pperp_W->Fill(W_var,p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z]) *TMath::RadToDeg(),P_2p1pito2p0pi[z]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_Ecal->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_Ecal_Reso->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
					h2_Ecal_Etrue->Fill(Ecal_2p1pi_to2p0pi[z],Ev,P_2p1pito2p0pi[z]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,P_2p1pito2p0pi[z]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,P_2p1pito2p0pi[z]*histoweight);
					h1_Wvar_weight->Fill(W_var,P_2p1pito2p0pi[z]*histoweight);
					h1_nu_weight->Fill(nu,P_2p1pito2p0pi[z]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),P_2p1pito2p0pi[z]*histoweight);
					h1_prot_mom->Fill(V3_2prot_corr[z].Mag(),P_2p1pito2p0pi[z]*histoweight);
					h1_MissMomentum->Fill(p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(Ecal_2p1pi_to2p0pi[z],V4_el);
					LocalWeight = P_2p1pito2p0pi[z]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(Ecal_2p1pi_to2p0pi[z],LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_miss_perp_2p1pi_to2p0pi[z],LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

 					if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(Ecal_2p1pi_to2p0pi[z],LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_miss_perp_2p1pi_to2p0pi[z],LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

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
					h2_Etot_pperp->Fill(p_perp_tot_2p[z],E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_E_tot_2p1pi_1p1pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
					h1_E_rec_2p1pi_1p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
					h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_2p1pito1p1pi[z]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_Ecal->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_Ecal_Reso->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
					h2_Ecal_Etrue->Fill(E_tot_2p[z],Ev,P_2p1pito1p1pi[z]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,P_2p1pito1p1pi[z]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,P_2p1pito1p1pi[z]*histoweight);
					h1_Wvar_weight->Fill(W_var,P_2p1pito1p1pi[z]*histoweight);
					h1_nu_weight->Fill(nu,P_2p1pito1p1pi[z]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),P_2p1pito1p1pi[z]*histoweight);
					h1_prot_mom->Fill(V3_2prot_corr[z].Mag(),P_2p1pito1p1pi[z]*histoweight);
					h1_MissMomentum->Fill(p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot_2p[z],V4_el);
					LocalWeight = P_2p1pito1p1pi[z]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_tot_2p[z],LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_perp_tot_2p[z],LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

 					if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_tot_2p[z],LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_perp_tot_2p[z],LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

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
					h2_Etot_pperp->Fill(p_perp_tot_2p[z],E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h1_E_tot_2p1pi_1p0pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p0pi[z]*histoweight);
					h1_E_rec_2p1pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p0pi[z]*histoweight);
					h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),-P_2p1pito1p0pi[z]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h1_Ecal->Fill(E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h1_Ecal_Reso->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-P_2p1pito1p0pi[z]*histoweight);
					h2_Ecal_Etrue->Fill(E_tot_2p[z],Ev,-P_2p1pito1p0pi[z]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,-P_2p1pito1p0pi[z]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,-P_2p1pito1p0pi[z]*histoweight);
					h1_Wvar_weight->Fill(W_var,-P_2p1pito1p0pi[z]*histoweight);
					h1_nu_weight->Fill(nu,-P_2p1pito1p0pi[z]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),-P_2p1pito1p0pi[z]*histoweight);
					h1_prot_mom->Fill(V3_2prot_corr[z].Mag(),-P_2p1pito1p0pi[z]*histoweight);
					h1_MissMomentum->Fill(p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot_2p[z],V4_el);
					LocalWeight = -P_2p1pito1p0pi[z]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_tot_2p[z],LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_perp_tot_2p[z],LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

 					if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_tot_2p[z],LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_perp_tot_2p[z],LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0; i < n_slice; i++) {

						if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){

							h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_tot_2p[z],P_2p1pito1p0pi[z]*histoweight);
							h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[i]->Fill(E_rec,P_2p1pito1p0pi[z]*histoweight);
						}
					}

				}//filling the histograms for 2protons

			}//1pi requirement

			//---------------------------------- 2p 2pi ----------------------------------------------

			const int N_2pi=2;
			double Ecal_2p2pi[N_2prot];
			double p_miss_perp_2p2pi[N_2prot];
			double Ptot_2p[2]={0};

			if (num_pi_phot == 2) {

				TVector3 V3_2pi_corr[N_2pi];
				double pion_acc_ratio[N_2pi] = {1};
				for (int i = 0; i < num_pi_phot; i++) {

					if (choice == 0) { //CLAS data
						V3_2pi_corr[i].SetXYZ(pxf[ind_pi_phot[i]],pyf[ind_pi_phot[i]],pzf[ind_pi_phot[i]]);
						pion_acc_ratio[i] = 1; //Acceptance 1 for CLAS data
					}

					if (choice == 1) { //GENIE data
						pion_acc_ratio[i] = 0; //reset to 0 just to be same
						V3_2pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
								Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						double phi_pion = V3_2pi_corr[i].Phi();
						V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_2pi_corr[i].Theta();
						double pion_mom_corr = V3_2pi_corr[i].Mag();

						if (charge_pi[i] == 1) { //acceptance for pi plus
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == -1) {		//acceptance for pi minus. using electron acceptance map
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == 0) {		//acceptance for photon set to 1 for now F.H. 09/24/19
							pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 2proton and 2 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue;
						}
					}
				}

				rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,charge_pi, V4_el, Ecal_2p2pi,p_miss_perp_2p2pi,Ptot_2p);

				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
				double histoweight = weight_pions * weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Is this correct in the following loop? F.H. 09/01/19


				for(int z = 0; z < N_2prot; z++){ //looping over two protons

					//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------

					h1_E_tot_2p2pi->Fill(E_tot_2p[z], Ptot_2p[z]*histoweight);
					h1_E_rec_2p2pi->Fill(E_rec,Ptot_2p[z]*histoweight);
					h2_Erec_pperp_2p2pi->Fill(p_perp_tot_2p[z],E_rec,Ptot_2p[z]*histoweight);
					h2_Etot_pperp->Fill(p_perp_tot_2p[z],E_tot_2p[z],Ptot_2p[z]*histoweight);
					h1_E_tot_2p2pi_fracfeed->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],Ptot_2p[z]*histoweight);
					h1_E_rec_2p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],Ptot_2p[z]*histoweight);
					h2_pperp_W->Fill(W_var,p_perp_tot_2p[z],Ptot_2p[z]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),Ptot_2p[z]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_tot_2p[z],Ptot_2p[z]*histoweight);
					h1_Ecal->Fill(E_tot_2p[z],Ptot_2p[z]*histoweight);
					h1_Ecal_Reso->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],Ptot_2p[z]*histoweight);
					h2_Ecal_Etrue->Fill(E_tot_2p[z],Ev,Ptot_2p[z]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_tot_2p[z],Ptot_2p[z]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot_2p[z],Ptot_2p[z]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot_2p[z],Ptot_2p[z]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,Ptot_2p[z]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,Ptot_2p[z]*histoweight);
					h1_Wvar_weight->Fill(W_var,Ptot_2p[z]*histoweight);
					h1_nu_weight->Fill(nu,Ptot_2p[z]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),Ptot_2p[z]*histoweight);
					h1_prot_mom->Fill(V3_2prot_corr[z].Mag(),Ptot_2p[z]*histoweight);
					h1_MissMomentum->Fill(p_perp_tot_2p[z],Ptot_2p[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot_2p[z],V4_el);
					LocalWeight = Ptot_2p[z]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_tot_2p[z],LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_perp_tot_2p[z],LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

 					if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_tot_2p[z],LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_perp_tot_2p[z],LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0; i < n_slice; i++) {

						if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
							h1_Etot_p_bkgd_slice_2p2pi[i]->Fill(E_tot_2p[z],Ptot_2p[z]*histoweight);
							h1_Erec_p_bkgd_slice_2p2pi[i]->Fill(E_rec,Ptot_2p[z]*histoweight);
						}
					}

				} //Filling the histogram for two protons

			}//2pi requirement

		} //2prot requirement

		// -------------------------------------------------------------------------------------------------------------------------------------

		//Events with exactly 3 protons

		if(num_p == 3) {

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
			double p_acc_ratio[N_3p] = {1};

			for(int i = 0; i < N_3p; i++) {

				N_p1[i] = 0;

				V4_p_uncorr[i].SetPxPyPzE(pxf[index_p[i]],pyf[index_p[i]],pzf[index_p[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
				V3_prot_uncorr[i] = V4_p_uncorr[i].Vect();

				if (choice == 0) { //CLAS data

					V3_prot_corr[i].SetXYZ(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60]);
					V4_p_corr[i].SetPxPyPzE(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]+60]*pf[index_p[i]+60]));
					p_acc_ratio[i] = 1; //Acceptance is 1 for CLAS datafile
				}

				if (choice == 1) { //GENIE data

					p_acc_ratio[i] = 0; //Reset just to be sure
					V3_prot_corr[i].SetXYZ(Smeared_Pp[i]/pf[index_p[i]] * pxf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pyf[index_p[i]],
						Smeared_Pp[i]/pf[index_p[i]] * pzf[index_p[i]]);
					V4_p_corr[i].SetPxPyPzE(Smeared_Pp[i]/pf[index_p[i]] * pxf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pyf[index_p[i]],
						Smeared_Pp[i]/pf[index_p[i]] * pzf[index_p[i]],Smeared_Ep[i]);

					double phi_prot = V3_el.Phi(); //in Radians
					V3_prot_corr[i].SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double p_theta = V3_prot_corr[i].Theta();
					double prot_mom_corr = V3_prot_corr[i].Mag();
					//Proton acceptance weight
					p_acc_ratio[i] = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);
					if ( fabs(p_acc_ratio[i]) != p_acc_ratio[i] ) { continue; }
				}

				V4_prot_el[i] = V4_p_corr[i] + V4_el;
				E_cal[i] = V4_el.E()+ V4_p_corr[i].E() - m_prot + bind_en[ftarget];
				p_miss_perp[i] = TMath::Sqrt(V4_prot_el[i].Px()*V4_prot_el[i].Px() + V4_prot_el[i].Py()*V4_prot_el[i].Py());

			} //end loop over N_3p

			V3_prot_el[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
			V3_prot_el[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
			V3_prot_el[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
			V3_prot_el[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
			V3_prot_el[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
			V3_prot_el[2][1]=V4_el.Vect()+V3_prot_uncorr[2];


			rotation->prot3_rot_func( V3_prot_corr,V3_prot_uncorr,V4_el,E_cal_3pto2p,p_miss_perp_3pto2p, P_3pto2p,N_p1, E_cal_3pto1p,p_miss_perp_3pto1p,&N_p_three);

			//acceptance weight for all three protons ( = 1 for CLAS data)
			double weight_protons =	p_acc_ratio[0] * p_acc_ratio[1] * p_acc_ratio[2];

			if(num_pi_phot==0 && N_p_three!=0){

				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Weight for 3protons, 1 electron, GENIE weight and Mott cross section

				for(int count = 0; count < N_comb; count++) { //Loop over number of combinations

					for(int j = 0; j < N_2p; j++) { //loop over two protons

						//-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------

						h1_E_tot_3pto2p->Fill(E_cal_3pto2p[count][j], P_3pto2p[count][j]*histoweight);
						h1_E_rec_3pto2p->Fill(E_rec, P_3pto2p[count][j]*histoweight);
						h2_Erec_pperp_321p->Fill(p_miss_perp_3pto2p[count][j],E_rec,P_3pto2p[count][j]*histoweight);
						h2_Etot_pperp->Fill(p_miss_perp_3pto2p[count][j],E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
						h1_E_tot_3pto2p_fracfeed->Fill((E_cal_3pto2p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto2p[count][j]*histoweight);
						h1_E_rec_3pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], P_3pto2p[count][j]*histoweight);
						h2_pperp_W->Fill(W_var,p_miss_perp_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
						h1_theta0->Fill((V4_beam.Vect()).Angle(V3_prot_el[count][j])*TMath::RadToDeg(),P_3pto2p[count][j]*histoweight);
						h2_Ecal_Eqe->Fill(E_rec,E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
						h1_Ecal->Fill(E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
						h1_Ecal_Reso->Fill((E_cal_3pto2p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_3pto2p[count][j]*histoweight);
						h2_Ecal_Etrue->Fill(E_cal_3pto2p[count][j],Ev,P_3pto2p[count][j]*histoweight);
						h2_Etrue_Ecal->Fill(Ev,E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
						h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);
						h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_3pto2p[count][j],P_3pto2p[count][j]*histoweight);

						h1_xbjk_weight->Fill(x_bjk,P_3pto2p[count][j]*histoweight);
						h1_Q2_weight->Fill(reco_Q2,P_3pto2p[count][j]*histoweight);
						h1_Wvar_weight->Fill(W_var,P_3pto2p[count][j]*histoweight);
						h1_nu_weight->Fill(nu,P_3pto2p[count][j]*histoweight);
						h1_el_mom_corr->Fill(V4_el.Rho(),P_3pto2p[count][j]*histoweight);
						h1_prot_mom->Fill(V3_prot_corr[j].Mag(),P_3pto2p[count][j]*histoweight);
						h1_MissMomentum->Fill(p_miss_perp_3pto2p[count][j],P_3pto2p[count][j]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(E_cal_3pto2p[count][j],V4_el);
						LocalWeight = P_3pto2p[count][j]*histoweight;

						h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						ECal_BreakDown[0]->Fill(E_cal_3pto2p[count][j],LocalWeight);
						EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[0]->Fill(p_miss_perp_3pto2p[count][j],LocalWeight);
						Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[0]->Fill(nu,LocalWeight);
						Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

	 					if (choice == 1) {
							ECal_BreakDown[Interaction]->Fill(E_cal_3pto2p[count][j],LocalWeight);
							EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
							Pmiss_BreakDown[Interaction]->Fill(p_miss_perp_3pto2p[count][j],LocalWeight);
							Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
							Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
							Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++) {

							if (p_miss_perp_3pto2p[count][j]<pperp_max[i] && p_miss_perp_3pto2p[count][j]>pperp_min[i]) {

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
					h2_Etot_pperp->Fill(p_miss_perp[j],E_cal[j],-P_3pto1p[j]*histoweight);
					h1_E_tot_3pto1p_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_3pto1p[j]*histoweight);
					h1_E_rec_3pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_3pto1p[j]*histoweight);
					h2_pperp_W->Fill(W_var,p_miss_perp[j],-P_3pto1p[j]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),-P_3pto1p[j]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_cal[j],-P_3pto1p[j]*histoweight);
					h1_Ecal->Fill(E_cal[j],-P_3pto1p[j]*histoweight);
					h1_Ecal_Reso->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-P_3pto1p[j]*histoweight);
					h2_Ecal_Etrue->Fill(E_cal[j],Ev,-P_3pto1p[j]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_cal[j],-P_3pto1p[j]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],-P_3pto1p[j]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],-P_3pto1p[j]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,-P_3pto1p[j]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,-P_3pto1p[j]*histoweight);
					h1_Wvar_weight->Fill(W_var,-P_3pto1p[j]*histoweight);
					h1_nu_weight->Fill(nu,-P_3pto1p[j]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),-P_3pto1p[j]*histoweight);
					h1_prot_mom->Fill(V3_prot_corr[j].Mag(),-P_3pto1p[j]*histoweight);
					h1_MissMomentum->Fill(p_miss_perp[j],-P_3pto1p[j]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_cal[j],V4_el);
					LocalWeight = -P_3pto1p[j]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_cal[j],LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_miss_perp[j],LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

	 				if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_cal[j],LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_miss_perp[j],LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0; i < n_slice; i++) {
					
						if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
							h1_Etot_3pto1p_slice[i]->Fill(E_cal[j],P_3pto1p[j]*histoweight);
							h1_Erec_3pto1p_slice[i]->Fill(E_rec,P_3pto1p[j]*histoweight);
						}
					}

				} //end loop over N_3p

			} //end if num_pi_phot==0 && N_p_three!=0, no pions

			//----------------------------------3p 1pi ----------------------------------------------------------

			if (num_pi_phot==1) {

				double P_tot_3p[N_3p]={0};
				double Ecal_3p1pi[N_3p]={0};
				double p_miss_perp_3p1pi[N_3p]={0};
				TVector3 V3_pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 0) { //CLAS data
					V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);
					pion_acc_ratio = 1; //Acceptance is 1 for CLAS datafile
				}

				if (choice == 1){ //GENIE data

					pion_acc_ratio = 0; //Reset to 0 just to be sure
					V3_pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);

					double phi_pion = V3_pi_corr.Phi(); //in Radians
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_pi_corr.Theta();
					double pion_mom_corr = V3_pi_corr.Mag();

					if (charge_pi[0] == 1) { //acceptance for pi plus
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == -1) {		//acceptance for pi minus. using electron acceptance map
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == 0) {		//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						 pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 3proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }

				}

				rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot_uncorr, V3_pi_corr, charge_pi[0] , V4_el,	Ecal_3p1pi,p_miss_perp_3p1pi, P_tot_3p);

				//for CLAS data is histoweight = 1/Mott_cross_sec
				double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Weight for 3protons, 1 pion, 1 electron, GENIE weight and Mott cross section

				for(int j = 0; j < N_3p; j++) { //loop over 3 protons

					h1_E_tot_3p1pi->Fill(E_cal[j], P_tot_3p[j]*histoweight);
					h1_E_rec_3p1pi->Fill(E_rec,P_tot_3p[j]*histoweight);
					h2_Erec_pperp_3p1pi->Fill(p_miss_perp[j],E_rec,P_tot_3p[j]*histoweight);
					h2_Etot_pperp->Fill(p_miss_perp[j],E_cal[j],P_tot_3p[j]*histoweight);
					h1_E_tot_3p1pi_fracfeed->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_tot_3p[j]*histoweight);
					h1_E_rec_3p1pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_3p[j]*histoweight);
					h2_pperp_W->Fill(W_var,p_miss_perp[j],P_tot_3p[j]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_tot_3p[j]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_cal[j],P_tot_3p[j]*histoweight);
					h1_Ecal->Fill(E_cal[j],P_tot_3p[j]*histoweight);
					h1_Ecal_Reso->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_tot_3p[j]*histoweight);
					h2_Ecal_Etrue->Fill(E_cal[j],Ev,P_tot_3p[j]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_cal[j],P_tot_3p[j]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal[j],P_tot_3p[j]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal[j],P_tot_3p[j]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,P_tot_3p[j]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,P_tot_3p[j]*histoweight);
					h1_Wvar_weight->Fill(W_var,P_tot_3p[j]*histoweight);
					h1_nu_weight->Fill(nu,P_tot_3p[j]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),P_tot_3p[j]*histoweight);
					h1_prot_mom->Fill(V3_prot_corr[j].Mag(),P_tot_3p[j]*histoweight);
					h1_MissMomentum->Fill(p_miss_perp[j],P_tot_3p[j]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_cal[j],V4_el);
					LocalWeight = P_tot_3p[j]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_cal[j],LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_miss_perp[j],LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

	 				if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_cal[j],LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_miss_perp[j],LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0; i < n_slice; i++)
					{
						if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
							h1_Etot_3p1pi_slice[i]->Fill(E_cal[j],P_tot_3p[j]*histoweight);
							h1_Erec_3p1pi_slice[i]->Fill(E_rec,P_tot_3p[j]*histoweight);
						}
					
					}

				} //end loop over N_3p

			} // 1 pi requirement ends

		} //end if num_p == 3  3proton requirement

		// ------------------------------------------------------------------------------------------------------------------------------------------

		//Events with exactly 4 protons

	 	if(num_p==4) {

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
			double p_acc_ratio[N_p4] = {1};

			for(int i = 0; i < N_p4; i++)  //loop over 4 protons
			{

				V4_p4_uncorr[i].SetPxPyPzE(pxf[index_p[i]],pyf[index_p[i]],pzf[index_p[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p[i]]*pf[index_p[i]]));
				V3_prot4_uncorr[i] = V4_p4_uncorr[i].Vect();

				if (choice == 0) { //CLAS data

					V3_prot4_corr[i].SetXYZ(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60]);
					V4_p4_corr[i].SetPxPyPzE(pxf[index_p[i]+60], pyf[index_p[i]+60], pzf[index_p[i]+60], TMath::Sqrt(m_prot*m_prot+pf[index_p[i]+60]*pf[index_p[i]+60]));
					p_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
				}

				if (choice == 1) { //GENIE data

					p_acc_ratio[i] = 0; //Reset to 0 just to be sure
					V3_prot4_corr[i].SetXYZ(Smeared_Pp[i]/pf[index_p[i]] * pxf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pyf[index_p[i]],
							Smeared_Pp[i]/pf[index_p[i]] * pzf[index_p[i]]);
					V4_p4_corr[i].SetPxPyPzE(Smeared_Pp[i]/pf[index_p[i]] * pxf[index_p[i]],Smeared_Pp[i]/pf[index_p[i]] * pyf[index_p[i]],
							Smeared_Pp[i]/pf[index_p[i]] * pzf[index_p[i]],Smeared_Ep[i]);

					double phi_prot = V3_prot4_corr[i].Phi(); //in Radians
					V3_prot4_corr[i].SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double p_theta = V3_prot4_corr[i].Theta();
					double prot_mom_corr = V3_prot4_corr[i].Mag();
					//Proton acceptance weight
					p_acc_ratio[i] = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);
					if ( fabs(p_acc_ratio[i]) != p_acc_ratio[i] ) { continue; }

				}

				V4_prot4_el[i] = V4_p4_corr[i] + V4_el;
				E_cal_p4[i] = V4_el.E() + V4_p4_corr[i].E() - m_prot + bind_en[ftarget];
				p_miss_perp_p4[i] = TMath::Sqrt(V4_prot4_el[i].Px()*V4_prot4_el[i].Px() + V4_prot4_el[i].Py()*V4_prot4_el[i].Py());
	
			} //end loop over N_p4

			//acceptance weight for all four protons. It is 1 for CLAS data
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

	        		} //for loop of 4p rotations ends

				// still no pions
				int N_comb = 3;    //number of 2 proton combination out of three
				const int N_2p = 2, N_3p = 3;
				double E_cal_4pto3p[3][N_2p] = {0};
				double p_miss_perp_4pto3p[3][N_2p] = {0};
				double P_4pto3p[3][N_2p] = {0};
				TVector3 V3_prot_corr[N_3p],V3_prot_uncorr[N_3p],V3_prot_el_4pto3p[N_3p][N_2p],V3_el_prot[N_comb][N_2p];
				double N_p_three=0,N_p1[N_3p] = {0};
				double E_cal_43pto1p[N_3p];
				double p_miss_perp_43pto1p[N_3p];
				double P_43pto1p[3] = {0};

				double histoweight = weight_protons * e_acc_ratio * wght/Mott_cross_sec; 
				//Weight for 3protons, 1 electron, GENIE weight and Mott cross section

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

								//-----------------------------------------  4p to 3p->2->1  ----------------------------------------------

								h1_E_tot_4pto3p->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_E_rec_4pto3p->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_Erec_pperp_4321p->Fill(p_miss_perp_4pto3p[count][j],E_rec,P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_Etot_pperp->Fill(p_miss_perp_4pto3p[count][j],E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_E_tot_4pto3p_fracfeed->Fill((E_cal_4pto3p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], 
												P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_E_rec_4pto3p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], 
												P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_pperp_W->Fill(W_var,p_miss_perp_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_theta0->Fill((V4_beam.Vect()).Angle(V3_el_prot[count][j])*TMath::RadToDeg(),-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_Ecal->Fill(E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_Ecal_Reso->Fill((E_cal_4pto3p[count][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_Ecal_Etrue->Fill(E_cal_4pto3p[count][j],Ev,-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_Etrue_Ecal->Fill(Ev,E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);

								h1_xbjk_weight->Fill(x_bjk,-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_Q2_weight->Fill(reco_Q2,-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_Wvar_weight->Fill(W_var,-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_nu_weight->Fill(nu,-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_el_mom_corr->Fill(V4_el.Rho(),-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_prot_mom->Fill(V3_prot_corr[j].Mag(),-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
								h1_MissMomentum->Fill(p_miss_perp_4pto3p[count][j],-P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);

								// -----------------------------------------------------------------------------------------------
								// Reconstruct xB, W, Q2 using Ecal instead of Etrue

								CalKineVars = CalculateCalKineVars(E_cal_4pto3p[count][j],V4_el);
								LocalWeight = -P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight;

								h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
								h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
								h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
								h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

								h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
								if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

								// Fill plots based on underlying interactions

								ECal_BreakDown[0]->Fill(E_cal_4pto3p[count][j],LocalWeight);
								EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
								Pmiss_BreakDown[0]->Fill(p_miss_perp_4pto3p[count][j],LocalWeight);
								Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
								Nu_BreakDown[0]->Fill(nu,LocalWeight);
								Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

				 				if (choice == 1) {
									ECal_BreakDown[Interaction]->Fill(E_cal_4pto3p[count][j],LocalWeight);
									EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
									Pmiss_BreakDown[Interaction]->Fill(p_miss_perp_4pto3p[count][j],LocalWeight);
									Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
									Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
									Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
								}

								// -----------------------------------------------------------------------------------------------

								for(int i = 0; i < n_slice; i++) {

									if (p_miss_perp_4pto3p[count][j]<pperp_max[i] && p_miss_perp_4pto3p[count][j]>pperp_min[i]){

										h1_Etot_4pto3p_slice[i]->Fill(E_cal_4pto3p[count][j], P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
									        h1_Erec_4pto3p_slice[i]->Fill(E_rec, P_4pto3p[count][j]*(N_p4_p3[g]/N_p_four)*histoweight);
									}
								}

							} //end loop over N_2p

						} //end loop over N_comb : number of 2 proton combination out of 3 protons

						//-----------------------------------------  4p to 3p->1p  -----------------------------------------------------------------------

						for(int j = 0; j < N_3p; j++) { //4p to 3p->1, looping through 1p out of 3p

							//P_43pto1p doesnt have to be an array, one local variable here
							P_43pto1p[j]= N_p1[j]/N_p_three*(N_p4_p3[g]/N_p_four);
							h1_E_tot_43pto1p->Fill(E_cal_43pto1p[j], P_43pto1p[j]*histoweight);
							h1_E_rec_43pto1p->Fill(E_rec,P_43pto1p[j]*histoweight);
							h2_Erec_pperp_431p->Fill(p_miss_perp_43pto1p[j],E_rec,P_43pto1p[j]*histoweight);
							h2_Etot_pperp->Fill(p_miss_perp_43pto1p[j],E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
							h1_E_tot_43pto1p_fracfeed->Fill((E_cal_43pto1p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_43pto1p[j]*histoweight);
							h1_E_rec_43pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_43pto1p[j]*histoweight);
							h2_pperp_W->Fill(W_var,p_miss_perp_43pto1p[j],P_43pto1p[j]*histoweight);
							h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_43pto1p[j]*histoweight);
							h2_Ecal_Eqe->Fill(E_rec,E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
							h1_Ecal->Fill(E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
							h1_Ecal_Reso->Fill((E_cal_43pto1p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_43pto1p[j]*histoweight);
							h2_Ecal_Etrue->Fill(E_cal_43pto1p[j],Ev,P_43pto1p[j]*histoweight);
							h2_Etrue_Ecal->Fill(Ev,E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
							h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_43pto1p[j],P_43pto1p[j]*histoweight);
							h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_43pto1p[j],P_43pto1p[j]*histoweight);

							h1_xbjk_weight->Fill(x_bjk,P_43pto1p[j]*histoweight);
							h1_Q2_weight->Fill(reco_Q2,P_43pto1p[j]*histoweight);
							h1_Wvar_weight->Fill(W_var,P_43pto1p[j]*histoweight);
							h1_nu_weight->Fill(nu,P_43pto1p[j]*histoweight);
							h1_el_mom_corr->Fill(V4_el.Rho(),P_43pto1p[j]*histoweight);
							h1_prot_mom->Fill(V3_prot_corr[j].Mag(),P_43pto1p[j]*histoweight);
							h1_MissMomentum->Fill(p_miss_perp_43pto1p[j],P_43pto1p[j]*histoweight);

							// -----------------------------------------------------------------------------------------------
							// Reconstruct xB, W, Q2 using Ecal instead of Etrue

							CalKineVars = CalculateCalKineVars(E_cal_43pto1p[j],V4_el);
							LocalWeight = P_43pto1p[j]*histoweight;

							h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
							h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
							h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
							h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

							h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
							if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

							// Fill plots based on underlying interactions

							ECal_BreakDown[0]->Fill(E_cal_43pto1p[j],LocalWeight);
							EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
							Pmiss_BreakDown[0]->Fill(p_miss_perp_43pto1p[j],LocalWeight);
							Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
							Nu_BreakDown[0]->Fill(nu,LocalWeight);
							Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

				 			if (choice == 1) {
								ECal_BreakDown[Interaction]->Fill(E_cal_43pto1p[j],LocalWeight);
								EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
								Pmiss_BreakDown[Interaction]->Fill(p_miss_perp_43pto1p[j],LocalWeight);
								Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
								Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
								Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
							}

							// -----------------------------------------------------------------------------------------------

							for(int i = 0; i <n_slice; i++) {

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

							rotation->prot2_rot_func( V3p2, V3p2_uncorr, V4_el,E_cal_4pto2p,p_miss_perp_4pto2p,  P_4pto2p, &N_two);

							if( N_two!=0  && N_p_four!=0){

								for(int j = 0; j < N_2p; j++)  {  //looping through  1 proton combination out of 2 protons

									h1_E_tot_4pto2p->Fill(E_cal_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_E_rec_4pto2p->Fill(E_rec, P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_Erec_pperp_421p->Fill( p_miss_perp_4pto2p[j],E_rec,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_Etot_pperp->Fill( p_miss_perp_4pto2p[j],E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_E_tot_4pto2p_fracfeed->Fill((E_cal_4pto2p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], 
													P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_E_rec_4pto2p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en], 
													P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_pperp_W->Fill(W_var,p_miss_perp_4pto2p[j], P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3p2_uncorr[j])*TMath::RadToDeg(),
										P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_Ecal_Eqe->Fill(E_rec,E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_Ecal->Fill(E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_Ecal_Reso->Fill((E_cal_4pto2p[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_Ecal_Etrue->Fill(E_cal_4pto2p[j],Ev,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_Etrue_Ecal->Fill(Ev,E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);

									h1_xbjk_weight->Fill(x_bjk,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_Q2_weight->Fill(reco_Q2,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_Wvar_weight->Fill(W_var,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_nu_weight->Fill(nu,P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_el_mom_corr->Fill(V4_el.Rho(),P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_prot_mom->Fill(V3p2[j].Mag(),P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);
									h1_MissMomentum->Fill(p_miss_perp_4pto2p[j],P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight);

									// -----------------------------------------------------------------------------------------------
									// Reconstruct xB, W, Q2 using Ecal instead of Etrue

									CalKineVars = CalculateCalKineVars(E_cal_4pto2p[j],V4_el);
									LocalWeight = P_4pto2p[j]*(N_p4_p2[N_4to2]/N_p_four)*histoweight;

									h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
									h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
									h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
									h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

									h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
									if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

									// Fill plots based on underlying interactions

									ECal_BreakDown[0]->Fill(E_cal_4pto2p[j],LocalWeight);
									EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
									Pmiss_BreakDown[0]->Fill(p_miss_perp_4pto2p[j],LocalWeight);
									Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
									Nu_BreakDown[0]->Fill(nu,LocalWeight);
									Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

						 			if (choice == 1) {
										ECal_BreakDown[Interaction]->Fill(E_cal_4pto2p[j],LocalWeight);
										EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
										Pmiss_BreakDown[Interaction]->Fill(p_miss_perp_4pto2p[j],LocalWeight);
										Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
										Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
										Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
									}

									// -----------------------------------------------------------------------------------------------

									for(int i = 0; i < n_slice; i++){

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
						h2_Etot_pperp->Fill(p_miss_perp_p4[j],E_cal_p4[j],-P_4pto1p[j]*histoweight);
						h1_E_tot_4pto1p_fracfeed->Fill((E_cal_p4[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_4pto1p[j]*histoweight);
						h1_E_rec_4pto1p_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_4pto1p[j]*histoweight);
						h2_pperp_W->Fill(W_var,p_miss_perp_p4[j],-P_4pto1p[j]*histoweight);
						h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot4_uncorr[j])*TMath::RadToDeg(),-P_4pto1p[j]*histoweight);
						h2_Ecal_Eqe->Fill(E_rec,E_cal_p4[j],-P_4pto1p[j]*histoweight);
						h1_Ecal->Fill(E_cal_p4[j],-P_4pto1p[j]*histoweight);
						h1_Ecal_Reso->Fill((E_cal_p4[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-P_4pto1p[j]*histoweight);
						h2_Ecal_Etrue->Fill(E_cal_p4[j],Ev,-P_4pto1p[j]*histoweight);
						h2_Etrue_Ecal->Fill(Ev,E_cal_p4[j],-P_4pto1p[j]*histoweight);
						h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_cal_p4[j],-P_4pto1p[j]*histoweight);
						h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_cal_p4[j],-P_4pto1p[j]*histoweight);

						h1_xbjk_weight->Fill(x_bjk,-P_4pto1p[j]*histoweight);
						h1_Q2_weight->Fill(reco_Q2,-P_4pto1p[j]*histoweight);
						h1_Wvar_weight->Fill(W_var,-P_4pto1p[j]*histoweight);
						h1_nu_weight->Fill(nu,-P_4pto1p[j]*histoweight);
						h1_el_mom_corr->Fill(V4_el.Rho(),-P_4pto1p[j]*histoweight);
						h1_prot_mom->Fill(V3_prot_corr[j].Mag(),-P_4pto1p[j]*histoweight);
						h1_MissMomentum->Fill(p_miss_perp_p4[j],-P_4pto1p[j]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(E_cal_p4[j],V4_el);
						LocalWeight = -P_4pto1p[j]*histoweight;

						h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						ECal_BreakDown[0]->Fill(E_cal_p4[j],LocalWeight);
						EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[0]->Fill(p_miss_perp_p4[j],LocalWeight);
						Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[0]->Fill(nu,LocalWeight);
						Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

						if (choice == 1) {
							ECal_BreakDown[Interaction]->Fill(E_cal_p4[j],LocalWeight);
							EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
							Pmiss_BreakDown[Interaction]->Fill(p_miss_perp_p4[j],LocalWeight);
							Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
							Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
							Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++) {

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

		// --------------------------------------------------------------------------------------------------------------------------------------------------------

		//No Protons here, Next 150 lines are for the inclusive events

		h1_E_rec->Fill(E_rec,WeightIncl);

		if(num_pi_phot ==0){

			h1_E_rec_0pi->Fill(E_rec,WeightIncl);
			h1_E_rec_0pi_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],WeightIncl);

			// Inclusive Case BreakDown
			InclusiveEQE_BreakDown[0]->Fill(E_rec,WeightIncl);
			if (choice == 1 && WeightIncl > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,WeightIncl); }
		}

		//----------------------------- e- ,1pi  -----------------------------------------

		if(num_pi_phot == 1){

			TVector3 V3_pi_corr;
			double P_undet=0;
			double pion_acc_ratio = 1;

			if (choice == 0) { //CLAS data
				V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]], pyf[ind_pi_phot[0]], pzf[ind_pi_phot[0]]);
				pion_acc_ratio = 1; //acceptance is 1 for CLAS data
			}

			if (choice == 1) { //GENIE data
				pion_acc_ratio = 0; //reset just to be sure
				V3_pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
						Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);
				double phi_pion = V3_pi_corr.Phi(); //in Radians
				V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
				phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				double pion_theta = V3_pi_corr.Theta();
				double pion_mom_corr = V3_pi_corr.Mag();

				if (charge_pi[0] == 1) { //acceptance for pi plus
					pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
					if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
				}
				else if (charge_pi[0] == -1) {	 //acceptance for pi minus. using electron acceptance map
					pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
					if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
				}
				else if (charge_pi[0] == 0) {	 //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
					pion_acc_ratio = 1;
				}
				else { std::cout << "WARNING: 1 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

			}

			rotation->pi1_rot_func( V3_pi_corr, charge_pi[0], &P_undet);

			//histoweight is 1/Mott_cross_sec for CLAS data
			double histoweight = pion_acc_ratio * WeightIncl;

			h1_E_rec_1pi_weight->Fill(E_rec,P_undet*histoweight);
			h1_E_rec_1pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_undet*histoweight);

			// Inclusive Case BreakDown
			InclusiveEQE_BreakDown[0]->Fill(E_rec,-P_undet*histoweight);
			if (choice == 1 && P_undet*histoweight > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,-P_undet*histoweight); }

			if(!ec_radstat_n[0])  h1_E_rec_1pi->Fill(E_rec,histoweight);
			if(ec_num_n==1)     	h2_phot_e_angle_Erec->Fill(E_rec,V3_pi_corr.Angle(V4_el.Vect())*TMath::RadToDeg());
		}

		//----------------------------- e- ,2pi  -----------------------------------------

		if(num_pi_phot == 2) {

			const int N_2pi = 2;
			TVector3 V3_2pi_corr[N_2pi];
			double P_1pi[N_2pi] = {0};
			double P_0pi = 0;
			double pion_acc_ratio[N_2pi] = {1};

			for (int i = 0; i < num_pi_phot; i++) {

				if (choice == 0) { //CLAS data
					V3_2pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
					pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
				}

				if (choice == 1) { //GENIE data
					pion_acc_ratio[i] = 0; //Reset just to be secure
					V3_2pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
							Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
					double phi_pion = V3_2pi_corr[i].Phi(); //in Radians
					V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_2pi_corr[i].Theta();
					double pion_mom_corr = V3_2pi_corr[i].Mag();

					if (charge_pi[i] == 1) { //acceptance for pi plus
						pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
					}
					else if (charge_pi[i] == -1) {	//acceptance for pi minus. using electron acceptance map
						pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
					}
					else if (charge_pi[i] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						pion_acc_ratio[i] = 1;
					}
					else { std::cout << "WARNING: 2 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }
				}
			} //end loop over num_pi_phot

			rotation->pi2_rot_func(V3_2pi_corr, charge_pi, &P_0pi,P_1pi);
			//weight_pions is 1 for CLAS data
			double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
			//histoweight is 1/Mott_cross_sec for CLAS data
			double histoweight = weight_pions * e_acc_ratio * wght/Mott_cross_sec;

			//----------------------------- e- ,2pi->0pi (-) -----------------------------------------

			h1_E_rec_2pi_weight->Fill(E_rec,(-P_0pi)*histoweight);
			h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*histoweight);
			h1_E_rec_20pi->Fill(E_rec,(P_0pi)*histoweight);

			// Inclusive Case BreakDown
			InclusiveEQE_BreakDown[0]->Fill(E_rec,(-P_0pi)*histoweight);
			if (choice == 1 && P_0pi*histoweight > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,(-P_0pi)*histoweight); }

			//----------------------------- e- ,2pi->1pi->0pi (+)  -----------------------------------------

			for(int k = 0; k < N_2pi; k++){ //loop over two pions

				h1_E_rec_2pi_weight->Fill(E_rec,P_1pi[k]*histoweight);
				h1_E_rec_2pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[k]*histoweight);
				h1_E_rec_21pi->Fill(E_rec,(P_1pi[k])*histoweight);

				// Inclusive Case BreakDown
				InclusiveEQE_BreakDown[0]->Fill(E_rec,(P_1pi[k])*histoweight);
				if (choice == 1 && P_1pi[k]*histoweight > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,(P_1pi[k])*histoweight); }
			}

		} //end if for two pion events

		//----------------------------- e- ,3pi  -----------------------------------------

		if(num_pi_phot == 3){

			const int N_3pi=3;
			const int N_2pi=2;
			TVector3 V3_3pi_corr[N_3pi];
			double P_0pi = 0;
			double P_1pi[N_3pi]={0};
			double P_320pi[N_3pi]={0};
			double P_3210pi[N_3pi][N_2pi]={0};
			double pion_acc_ratio[N_3pi] = {1};

			for (int i = 0; i < num_pi_phot; i++) {

				if (choice == 0) { //CLAS data

					V3_3pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
					pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
				}

				if (choice == 1) { //GENIE data

					pion_acc_ratio[i] = 0; //Reset just to be sure
					V3_3pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
								Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
					double phi_pion = V3_3pi_corr[i].Phi(); //in Radians
					V3_3pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_3pi_corr[i].Theta();
					double pion_mom_corr = V3_3pi_corr[i].Mag();

					if (charge_pi[i] == 1) { //acceptance for pi plus
						pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
					}
					else if (charge_pi[i] == -1) {	//acceptance for pi minus. using electron acceptance map
						pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
					}
					else if (charge_pi[i] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						pion_acc_ratio[i] = 1;
					}
					else { std::cout << "WARNING: 3 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

				}

			} //end loop over num_pi_phot

			rotation->pi3_rot_func( V3_3pi_corr, charge_pi, &P_0pi, P_1pi, P_320pi,P_3210pi);
			//weight_pions is 1 for CLAS data
			double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2];
			//histoweight is 1/Mott_cross_sec for CLAS data
			double histoweight = weight_pions * e_acc_ratio * wght/Mott_cross_sec;

			//---------------------------3pi->0pi----------------------------------------------

			h1_E_rec_3pi_weight->Fill(E_rec,(-P_0pi)*histoweight);
			h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi)*histoweight);
			h1_E_rec_30pi->Fill(E_rec,(P_0pi)*histoweight);

			// Inclusive Case BreakDown
			InclusiveEQE_BreakDown[0]->Fill(E_rec,(P_0pi)*histoweight);
			if (choice == 1 && P_0pi*histoweight > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,(P_0pi)*histoweight); }

			for(int h = 0; h < N_3pi; h++){ //loop over three pions

				//---------------------------3pi->1pi->0pi----------------------------------------------

				h1_E_rec_3pi_weight->Fill(E_rec,P_1pi[h]*histoweight);
				h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1pi[h]*histoweight);
				h1_E_rec_310pi->Fill(E_rec,(P_1pi[h])*histoweight);

				// Inclusive Case BreakDown
				InclusiveEQE_BreakDown[0]->Fill(E_rec,(P_1pi[h])*histoweight);
				if (choice == 1 && P_1pi[h]*histoweight > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,(P_1pi[h])*histoweight); }

				//---------------------------3pi->2pi->0pi----------------------------------------------

				h1_E_rec_3pi_weight->Fill(E_rec,P_320pi[h]*histoweight);
				h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_320pi[h]*histoweight);
				h1_E_rec_320pi->Fill(E_rec,(P_320pi[h])*histoweight);

				// Inclusive Case BreakDown
				InclusiveEQE_BreakDown[0]->Fill(E_rec,(P_320pi[h])*histoweight);
				if (choice == 1 && P_320pi[h]*histoweight > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,(P_320pi[h])*histoweight); }

				//---------------------------3pi->2pi->1pi->0pi----------------------------------------------

				for(int g = 0; g < N_2pi; g++){ //loop over two pions

					h1_E_rec_3pi_weight->Fill(E_rec,(-P_3210pi[h][g])*histoweight);
					h1_E_rec_3pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_3210pi[h][g])*histoweight);
					h1_E_rec_3210pi->Fill(E_rec,(P_3210pi[h][g])*histoweight);

					// Inclusive Case BreakDown
					InclusiveEQE_BreakDown[0]->Fill(E_rec,(-P_3210pi[h][g])*histoweight);
					if (choice == 1 && P_3210pi[h][g]*histoweight > 0) { InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,(-P_3210pi[h][g])*histoweight); }

				}

			}//end of 3pi loop

		}//end of 3pi requirement

		//----------------------------- e- ,4pi  -----------------------------------------

		if(num_pi_phot == 4){

			const int N_4pi=4;
			TVector3 V3_4pi_corr[N_4pi];
			double P_0pi=0;
			double P_410pi=0;
			double P_420pi=0;
			double P_4210pi=0;
			double P_430pi=0;
			double P_4310pi=0;
			double P_4320pi=0;
			double P_43210pi=0;
			double pion_acc_ratio[N_4pi] = {1};

			for (int i = 0; i < num_pi_phot; i++) {

				if (choice == 0) { //CLAS data
					V3_4pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
					pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
				}

				if (choice == 1) { //GENIE data

						pion_acc_ratio[i] = 0; //Reset just to be sure
						V3_4pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						// apapadop
						double phi_pion = V3_4pi_corr[i].Phi(); //in Radians
						V3_4pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_4pi_corr[i].Theta();
						double pion_mom_corr = V3_4pi_corr[i].Mag();

						if (charge_pi[i] == 1) { //acceptance for pi plus
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == -1) {    //acceptance for pi minus. using electron acceptance map
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == 0) {    //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
							pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 4 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }
				}

			} //end loop over num_pi_phot

			rotation->pi4_rot_func(V3_4pi_corr, charge_pi, &P_0pi,&P_410pi,&P_420pi,&P_4210pi,&P_430pi,&P_4310pi,&P_4320pi,&P_43210pi);

			//weight_pions is 1 for CLAS data
			double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2] * pion_acc_ratio[3];
			//histoweight is 1/Mott_cross_sec for CLAS data
			double histoweight = weight_pions * e_acc_ratio * wght/Mott_cross_sec;

			//---------------------------4pi->0pi----------------------------------------------

			//why is it here not split like for 3pi case, sum over all weights is done here F.H 04.08.19
			h1_E_rec_4pi_weight->Fill(E_rec,(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight);
			h1_E_rec_4pi_weight_frac_feed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight);
			h1_E_rec_40pi->Fill(E_rec,(P_0pi)*histoweight);

			// Inclusive Case BreakDown
			InclusiveEQE_BreakDown[0]->Fill(E_rec,(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight);
			if (choice == 1 && (-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight != 0) 
				{ InclusiveEQE_BreakDown[Interaction]->Fill(E_rec,(-P_0pi+P_410pi+P_420pi-P_4210pi+P_430pi-P_4310pi-P_4320pi+P_43210pi)*histoweight); }

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

		if( num_p == 1) {

			//Vector for proton without momentum smearing
			TLorentzVector V4_prot_uncorr(pxf[index_p[0]],pyf[index_p[0]],pzf[index_p[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]]*pf[index_p[0]]));
			TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();

			//Vector for proton with momentum smearing or correction (energy loss)
			TVector3 V3_prot_corr;
			TLorentzVector V4_prot_corr;

			double p_acc_ratio = 1; //acceptance is 1 for CLAS data

			if (choice == 0) { //CLAS data
				V3_prot_corr.SetXYZ(pxf[index_p[0]+60], pyf[index_p[0]+60], pzf[index_p[0]+60]);
				V4_prot_corr.SetPxPyPzE(pxf[index_p[0]+60], pyf[index_p[0]+60], pzf[index_p[0]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p[0]+60]*pf[index_p[0]+60]));
			}
			if (choice == 1) { //GENIE data
				p_acc_ratio = 0; //Reset just to be sure
				//Fiducial cuts are done in the hadron loop
				//Vector for proton with momentum smearing
				V3_prot_corr.SetXYZ(Smeared_Pp[0]/pf[index_p[0]] * pxf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pyf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pzf[index_p[0]]);
				V4_prot_corr.SetPxPyPzE(Smeared_Pp[0]/pf[index_p[0]] * pxf[index_p[0]],Smeared_Pp[0]/pf[index_p[0]] * pyf[index_p[0]],
							Smeared_Pp[0]/pf[index_p[0]] * pzf[index_p[0]],Smeared_Ep[0]);

				double phi_prot = V3_prot_corr.Phi(); //in Radians
				V3_prot_corr.SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
				phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				//Proton kinematic variables
				double p_theta = V3_prot_corr.Theta();
				double prot_mom_corr = V3_prot_corr.Mag();
				//Proton weight
				p_acc_ratio = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);
				if ( fabs(p_acc_ratio) != p_acc_ratio ) { continue; }

			}

			TLorentzVector V4_prot_el_tot = V4_prot_corr + V4_el;

			double p_perp_tot = TMath::Sqrt(V4_prot_el_tot.Px()*V4_prot_el_tot.Px() + V4_prot_el_tot.Py()*V4_prot_el_tot.Py());
			double E_tot = V4_el.E() + V4_prot_corr.E() - m_prot + bind_en[ftarget];

			//These Histograms are events with 1 electron and  1 proton and multiple pions
			//histoweight_inc is 1/Mott_cross_sec for CLAS data
			double histoweight_inc = p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;
			double histoweight_NoMott = p_acc_ratio * e_acc_ratio * wght;

			h1_E_tot->Fill(E_tot,histoweight_inc);
			h1_E_rec_1prot->Fill(E_rec,histoweight_inc);
			h1_E_tot_1prot->Fill(E_tot,histoweight_inc);
			h2_Erec_pperp->Fill(p_perp_tot,E_rec,histoweight_inc);
			h2_Etot_pperp->Fill(p_perp_tot,E_tot,histoweight_inc);

			//---------------------------------- 1p 0pi   ----------------------------------------------

			// Main Plots

			if(num_pi_phot == 0){

				double ECalReso = (E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en];
				double EQEReso = (E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en];

				SignalEvents++;

				if (p_perp_tot < pperp_max[0]) { PMiss_FirstSlice++; }
				if (p_perp_tot > pperp_max[0] && p_perp_tot < pperp_max[1]) { PMiss_SecondSlice++; }
				if (p_perp_tot > pperp_max[1]) { PMiss_ThirdSlice++; }

				if (fabs(ECalReso)*100. < 5) { 
					ECalSignalEventsWithin5Perc++; 
					if (p_perp_tot < pperp_max[0]) { ECalSignalEventsWithin5Perc_FirstSlice++; }
					if (p_perp_tot > pperp_max[0] && p_perp_tot < pperp_max[1]) { ECalSignalEventsWithin5Perc_SecondSlice++; }
					if (p_perp_tot > pperp_max[1]) { ECalSignalEventsWithin5Perc_ThirdSlice++; }
				}

				if (fabs(EQEReso)*100. < 5) { 
					EQESignalEventsWithin5Perc++; 
					if (p_perp_tot < pperp_max[0]) { EQESignalEventsWithin5Perc_FirstSlice++; }
					if (p_perp_tot > pperp_max[0] && p_perp_tot < pperp_max[1]) { EQESignalEventsWithin5Perc_SecondSlice++; }
					if (p_perp_tot > pperp_max[1]) { EQESignalEventsWithin5Perc_ThirdSlice++; }
				}

				if (Interaction == 1) { QESignalEvents++; }
				else if (Interaction == 2) { MECSignalEvents++; }
				else if (Interaction == 3) { RESSignalEvents++; }
				else if (Interaction == 4) { DISSignalEvents++; }
				else { OtherSignalEvents++; }

				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;

				h2_Erec_pperp_newcut2->Fill(p_perp_tot,E_rec,histoweight);
				h2_Etot_pperp->Fill(p_perp_tot,E_tot,histoweight);
				h1_E_rec_cut2_new->Fill(E_rec,histoweight);
				h1_E_tot_cut2->Fill(E_tot,histoweight);
				h1_E_tot_cut2_fracfeed->Fill(ECalReso,histoweight);
				h1_E_rec_cut2_new_fracfeed->Fill(EQEReso,histoweight);
				h2_pperp_W->Fill(W_var,p_perp_tot,histoweight);
				h1_theta0->Fill((V4_beam.Vect()).Angle(V4_prot_el_tot.Vect()) *TMath::RadToDeg(),histoweight);
				h2_Ecal_Eqe->Fill(E_rec,E_tot,histoweight);
				h1_Ecal->Fill(E_tot,histoweight);
				h1_Ecal_Reso->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],histoweight);
				h2_Ecal_Etrue->Fill(E_tot,Ev,histoweight);
				h2_Etrue_Ecal->Fill(Ev,E_tot,histoweight);
				h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,histoweight);
				h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,histoweight);

				h1_xbjk_weight->Fill(x_bjk,histoweight);
				h1_Q2_weight->Fill(reco_Q2,histoweight);
				h1_Wvar_weight->Fill(W_var,histoweight);
				h1_nu_weight->Fill(nu,histoweight);
				h1_el_mom_corr->Fill(V4_el.Rho(),histoweight);
				h1_prot_mom->Fill(V3_prot_corr.Mag(),histoweight);
				h1_MissMomentum->Fill(p_perp_tot,histoweight);

				// -----------------------------------------------------------------------------------------------

				// Unweighted plots for number of events

				h1_MissMomentum_NoWeight->Fill(p_perp_tot,histoweight_NoMott);

				h1_ECal_Slice0_NoWeight->Fill(E_tot,histoweight_NoMott);
				h1_EQE_Slice0_NoWeight->Fill(E_rec,histoweight_NoMott);

				if (p_perp_tot < pperp_max[0]) { h1_ECal_Slice1_NoWeight->Fill(E_tot,histoweight_NoMott); h1_EQE_Slice1_NoWeight->Fill(E_rec,histoweight_NoMott); }
				if (p_perp_tot > pperp_max[0] && p_perp_tot < pperp_max[1]) { h1_ECal_Slice2_NoWeight->Fill(E_tot,histoweight_NoMott); h1_EQE_Slice2_NoWeight->Fill(E_rec,histoweight_NoMott); }
				if (p_perp_tot > pperp_max[1]) { h1_ECal_Slice3_NoWeight->Fill(E_tot,histoweight_NoMott); h1_EQE_Slice3_NoWeight->Fill(E_rec,histoweight_NoMott); }

				// -----------------------------------------------------------------------------------------------
				// Reconstruct xB, W, Q2 using Ecal instead of Etrue

				CalKineVars = CalculateCalKineVars(E_tot,V4_el);
				LocalWeight = histoweight;

				h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
				h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
				h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
				h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

				h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
				if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

				// Fill plots based on underlying interactions

				ECal_BreakDown[0]->Fill(E_tot,LocalWeight);
				EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
				Pmiss_BreakDown[0]->Fill(p_perp_tot,LocalWeight);
				Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
				Nu_BreakDown[0]->Fill(nu,LocalWeight);
				Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

				if (choice == 1) {
					ECal_BreakDown[Interaction]->Fill(E_tot,LocalWeight);
					EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[Interaction]->Fill(p_perp_tot,LocalWeight);
					Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
					Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
				}

				// -----------------------------------------------------------------------------------------------

				for(int i = 0; i < n_slice; i++) {

					if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){

						h1_Etot_Npi0[i]->Fill(E_tot,histoweight);
						h1_Erec_Npi0[i]->Fill(E_rec,histoweight);

					}
				}

				if (p_perp_tot < 0.2){
					
					h1_E_rec_cut005_newcut3->Fill(E_rec,histoweight);
					h2_Erec_pperp_cut3->Fill(p_perp_tot,E_rec,histoweight);
					h2_Etot_pperp->Fill(p_perp_tot,E_tot,histoweight);
				}

			} //num pi=0

			//---------------------------------- 1p 1pi   ----------------------------------------------

			if(num_pi_phot == 1){

				double N_piphot_det;
				double N_piphot_undet;
				TVector3 V3_pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 0) { //CLAS data

					pion_acc_ratio = 1; //Acceptance is 1 for CLAS data
					V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]], pyf[ind_pi_phot[0]], pzf[ind_pi_phot[0]]);
				}

				if (choice == 1) { //GENIE data

					pion_acc_ratio = 1; //Reset to 0 just to be sure
					V3_pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);
					double phi_pion = V3_pi_corr.Phi(); //in Radians
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_pi_corr.Theta();
					double pion_mom_corr = V3_pi_corr.Mag();

					if (charge_pi[0] == 1) { //acceptance for pi plus
					
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == -1) { //acceptance for pi minus. using electron acceptance map
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi[0] == 0) { //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 1 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

				}

				rotation->prot1_pi1_rot_func(V3_prot_uncorr,V3_pi_corr, charge_pi[0], &N_piphot_det,&N_piphot_undet);

				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = pion_acc_ratio * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; 
				//1proton, 1 Pion, 1 electron acceptance, GENIE weight and Mott

				if(N_piphot_det!=0){

					h1_E_rec_undetfactor->Fill(E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
					h1_E_tot_undetfactor->Fill(E_tot,(N_piphot_undet/N_piphot_det)*histoweight);
					h2_Erec_pperp_1p1pi->Fill(p_perp_tot,E_rec,(N_piphot_undet/N_piphot_det)*histoweight);
					h2_Etot_pperp->Fill(p_perp_tot,E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_E_rec_undetfactor_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],(N_piphot_undet/N_piphot_det)*histoweight);
					h1_E_tot_undetfactor_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],(N_piphot_undet/N_piphot_det)*histoweight);
					h2_pperp_W->Fill(W_var,p_perp_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-(N_piphot_undet/N_piphot_det)*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_Ecal->Fill(E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_Ecal_Reso->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-(N_piphot_undet/N_piphot_det)*histoweight);
					h2_Ecal_Etrue->Fill(E_tot,Ev,-(N_piphot_undet/N_piphot_det)*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-(N_piphot_undet/N_piphot_det)*histoweight);

					h1_xbjk_weight->Fill(x_bjk,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_Q2_weight->Fill(reco_Q2,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_Wvar_weight->Fill(W_var,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_nu_weight->Fill(nu,-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_prot_mom->Fill(V3_prot_corr.Mag(),-(N_piphot_undet/N_piphot_det)*histoweight);
					h1_MissMomentum->Fill(p_perp_tot,-(N_piphot_undet/N_piphot_det)*histoweight);

					// -----------------------------------------------------------------------------------------------
					// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot,V4_el);
					LocalWeight = -(N_piphot_undet/N_piphot_det)*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_tot,LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_perp_tot,LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

					if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_tot,LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_perp_tot,LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0; i < n_slice; i++) {

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

			if(num_pi_phot == 2) {

				const int N_2pi=2;
				TVector3 V3_2pi_corr[N_2pi],V3_2pi_rot[N_2pi],V3_p_rot;
				double P_1p0pi=0;
				double P_1p1pi[N_2pi]={1};

				double pion_acc_ratio[N_2pi] = {0};

				for (int i = 0; i < num_pi_phot; i++) {

					if (choice == 0) { //CLAS data
						V3_2pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
						pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
					}

					if (choice == 1) { //GENIE data

						V3_2pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						double phi_pion = V3_2pi_corr[i].Phi(); //in Radians
						V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_2pi_corr[i].Theta();
						double pion_mom_corr = V3_2pi_corr[i].Mag();

						if (charge_pi[i] == 1) { //acceptance for pi plus
							 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
							  if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == -1) {	//acceptance for pi minus. using electron acceptance map
							 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
							  if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
							 pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 1 Proton 2 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }
					}

				} //end loop over num_pi_phot

				rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi_corr,charge_pi,&P_1p0pi,P_1p1pi);
				//weight_pions is 1 for CLAS data
				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; 
				//1proton, 2 Pion, 1 electron acceptance, GENIE weight and Mott

				//---------------------------------- 1p 2pi->1p1pi   ----------------------------------------------

				for(int z = 0; z < N_2pi; z++){  //to consider 2 diff. 1pi states

					h1_E_tot_1p2pi->Fill(E_tot,P_1p1pi[z]*histoweight);
					h1_E_rec_1p2pi->Fill(E_rec,P_1p1pi[z]*histoweight);
					h2_Erec_pperp_1p2pi_1p1pi->Fill(p_perp_tot,E_rec,P_1p1pi[z]*histoweight);
					h2_Etot_pperp->Fill(p_perp_tot,E_tot,P_1p1pi[z]*histoweight);
					h1_E_tot_1p2pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
					h1_E_rec_1p2pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p1pi[z]*histoweight);
					h2_pperp_W->Fill(W_var,p_perp_tot,P_1p1pi[z]*histoweight);
					h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p1pi[z]*histoweight);
					h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p1pi[z]*histoweight);
					h1_Ecal->Fill(E_tot,P_1p1pi[z]*histoweight);
					h1_Ecal_Reso->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
					h2_Ecal_Etrue->Fill(E_tot,Ev,P_1p1pi[z]*histoweight);
					h2_Etrue_Ecal->Fill(Ev,E_tot,P_1p1pi[z]*histoweight);
					h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p1pi[z]*histoweight);
					h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p1pi[z]*histoweight);

					h1_xbjk_weight->Fill(x_bjk,P_1p1pi[z]*histoweight);
					h1_Q2_weight->Fill(reco_Q2,P_1p1pi[z]*histoweight);
					h1_Wvar_weight->Fill(W_var,P_1p1pi[z]*histoweight);
					h1_nu_weight->Fill(nu,P_1p1pi[z]*histoweight);
					h1_el_mom_corr->Fill(V4_el.Rho(),P_1p1pi[z]*histoweight);
					h1_prot_mom->Fill(V3_prot_corr.Mag(),P_1p1pi[z]*histoweight);
					h1_MissMomentum->Fill(p_perp_tot,P_1p1pi[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot,V4_el);
					LocalWeight = P_1p1pi[z]*histoweight;

					h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					ECal_BreakDown[0]->Fill(E_tot,LocalWeight);
					EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[0]->Fill(p_perp_tot,LocalWeight);
					Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[0]->Fill(nu,LocalWeight);
					Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

					if (choice == 1) {
						ECal_BreakDown[Interaction]->Fill(E_tot,LocalWeight);
						EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
						Pmiss_BreakDown[Interaction]->Fill(p_perp_tot,LocalWeight);
						Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
						Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
						Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

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
				h2_Etot_pperp->Fill(p_perp_tot,E_tot,-P_1p0pi*histoweight);
				h1_E_tot_1p2pi_1p0pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p0pi*histoweight);
				h1_E_rec_1p2pi_1p0pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p0pi*histoweight);
				h2_pperp_W->Fill(W_var,p_perp_tot,-P_1p0pi*histoweight);
				h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),-P_1p0pi*histoweight);
				h2_Ecal_Eqe->Fill(E_rec,E_tot,-P_1p0pi*histoweight);
				h1_Ecal->Fill(E_tot,-P_1p0pi*histoweight);
				h1_Ecal_Reso->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-P_1p0pi*histoweight);
				h2_Ecal_Etrue->Fill(E_tot,Ev,-P_1p0pi*histoweight);
				h2_Etrue_Ecal->Fill(Ev,E_tot,-P_1p0pi*histoweight);
				h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,-P_1p0pi*histoweight);
				h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,-P_1p0pi*histoweight);

				h1_xbjk_weight->Fill(x_bjk,-P_1p0pi*histoweight);
				h1_Q2_weight->Fill(reco_Q2,-P_1p0pi*histoweight);
				h1_Wvar_weight->Fill(W_var,-P_1p0pi*histoweight);
				h1_nu_weight->Fill(nu,-P_1p0pi*histoweight);
				h1_el_mom_corr->Fill(V4_el.Rho(),-P_1p0pi*histoweight);
				h1_prot_mom->Fill(V3_prot_corr.Mag(),-P_1p0pi*histoweight);
				h1_MissMomentum->Fill(p_perp_tot,-P_1p0pi*histoweight);

				// -----------------------------------------------------------------------------------------------
				// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

				CalKineVars = CalculateCalKineVars(E_tot,V4_el);
				LocalWeight = -P_1p0pi*histoweight;

				h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
				h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
				h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
				h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

				h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
				if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

				// Fill plots based on underlying interactions

				ECal_BreakDown[0]->Fill(E_tot,LocalWeight);
				EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
				Pmiss_BreakDown[0]->Fill(p_perp_tot,LocalWeight);
				Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
				Nu_BreakDown[0]->Fill(nu,LocalWeight);
				Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

				if (choice == 1) {
					ECal_BreakDown[Interaction]->Fill(E_tot,LocalWeight);
					EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[Interaction]->Fill(p_perp_tot,LocalWeight);
					Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
					Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
				}

				// -----------------------------------------------------------------------------------------------

				for(int i = 0; i < n_slice; i++){

					if (p_perp_tot<pperp_max[i] && p_perp_tot>pperp_min[i]){

						h1_Etot_bkgd_1p2pi_1p0pi[i]->Fill(E_tot,P_1p0pi*histoweight);
						h1_Erec_bkgd_1p2pi_1p0pi[i]->Fill(E_rec,P_1p0pi*histoweight);
					}
				}

			}//1p 2pi statetment ends

			//---------------------------------- 1p 3pi   ----------------------------------------------

			if(num_pi_phot == 3){

				const int N_3pi=3;
				TVector3 V3_3pi_corr[N_3pi],V3_3pi_rot[N_3pi],V3_p_rot;
				double P_1p3pi = 0;
				double pion_acc_ratio[N_3pi] = {1};

				for (int i = 0; i < num_pi_phot; i++) {

					if (choice == 0) { //CLAS data
						V3_3pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
						pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
					}

					if (choice == 1) { //GENIE data
						pion_acc_ratio[i] = 0; //Reset to 0 just to be sure
						V3_3pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Pp[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Pp[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						double phi_pion = V3_3pi_corr[i].Phi(); //in Radians
						V3_3pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_3pi_corr[i].Theta();
						double pion_mom_corr = V3_3pi_corr[i].Mag();

						if (charge_pi[i] == 1) { //acceptance for pi plus
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == -1) {	//acceptance for pi minus. using electron acceptance map
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi[i] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
							pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 3 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }
					}

				} //end loop over num_pi_phot

				rotation->prot1_pi3_rot_func(V3_prot_uncorr, V3_3pi_corr, charge_pi, &P_1p3pi);
		 		//weight_pions is 1 for CLAS data
				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2];
				//histoweight is 1/Mott_cross_sec for CLAS data
				double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec; 
				//1proton, 3 Pion, 1 electron acceptance, GENIE weight and Mott

				//---------------------------------- 1p 3pi->1p 0pi  total ?? F.H. 08/13/19 check logic here compared to 1p 2pi case ----------------------------

				h1_E_tot_1p3pi->Fill(E_tot,P_1p3pi*histoweight);
				h1_E_rec_1p3pi->Fill(E_rec,P_1p3pi*histoweight);
				h2_Erec_pperp_1p3pi->Fill(p_perp_tot,E_rec,P_1p3pi*histoweight);
				h2_Etot_pperp->Fill(p_perp_tot,E_tot,P_1p3pi*histoweight);
				h1_E_tot_1p3pi_fracfeed->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi*histoweight);
				h1_E_rec_1p3pi_fracfeed->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p3pi*histoweight);
				h2_pperp_W->Fill(W_var,p_perp_tot,P_1p3pi*histoweight);
				h1_theta0->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p3pi*histoweight);
				h2_Ecal_Eqe->Fill(E_rec,E_tot,P_1p3pi*histoweight);
				h1_Ecal->Fill(E_tot,P_1p3pi*histoweight);
				h1_Ecal_Reso->Fill((E_tot-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi*histoweight);
				h2_Ecal_Etrue->Fill(E_tot,Ev,P_1p3pi*histoweight);
				h2_Etrue_Ecal->Fill(Ev,E_tot,P_1p3pi*histoweight);
				h2_EqeEcalratio_Eqe->Fill(E_rec,E_rec/E_tot,P_1p3pi*histoweight);
				h2_EqeEcaldiff_Eqe->Fill(E_rec,E_rec-E_tot,P_1p3pi*histoweight);

				h1_xbjk_weight->Fill(x_bjk,P_1p3pi*histoweight);
				h1_Q2_weight->Fill(reco_Q2,P_1p3pi*histoweight);
				h1_Wvar_weight->Fill(W_var,P_1p3pi*histoweight);
				h1_nu_weight->Fill(nu,P_1p3pi*histoweight);
				h1_el_mom_corr->Fill(V4_el.Rho(),P_1p3pi*histoweight);
				h1_prot_mom->Fill(V3_prot_corr.Mag(),P_1p3pi*histoweight);
				h1_MissMomentum->Fill(p_perp_tot,P_1p3pi*histoweight);

				// -----------------------------------------------------------------------------------------------
				// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

				CalKineVars = CalculateCalKineVars(E_tot,V4_el);
				LocalWeight = P_1p3pi*histoweight;

				h1_nuCal_weight->Fill(CalKineVars.at(0),LocalWeight);
				h1_Q2Cal_weight->Fill(CalKineVars.at(1),LocalWeight);
				h1_xbjkCal_weight->Fill(CalKineVars.at(2),LocalWeight);
				h1_WvarCal_weight->Fill(CalKineVars.at(3),LocalWeight);

				h2_Q2_nu_weight->Fill(nu,reco_Q2,LocalWeight);
				if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector->Fill(nu,reco_Q2,LocalWeight); }

				// Fill plots based on underlying interactions

				ECal_BreakDown[0]->Fill(E_tot,LocalWeight);
				EQE_BreakDown[0]->Fill(E_rec,LocalWeight);
				Pmiss_BreakDown[0]->Fill(p_perp_tot,LocalWeight);
				Q2_BreakDown[0]->Fill(reco_Q2,LocalWeight);
				Nu_BreakDown[0]->Fill(nu,LocalWeight);
				Pe_BreakDown[0]->Fill(V4_el.Rho(),LocalWeight);

				if (choice == 1) {
					ECal_BreakDown[Interaction]->Fill(E_tot,LocalWeight);
					EQE_BreakDown[Interaction]->Fill(E_rec,LocalWeight);
					Pmiss_BreakDown[Interaction]->Fill(p_perp_tot,LocalWeight);
					Q2_BreakDown[Interaction]->Fill(reco_Q2,LocalWeight);
					Nu_BreakDown[Interaction]->Fill(nu,LocalWeight);
					Pe_BreakDown[Interaction]->Fill(V4_el.Rho(),LocalWeight);
				}

				// -----------------------------------------------------------------------------------------------

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

	for(int i = 0; i <= n_slice-1; i++) {

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

//	TH1F *h_Erec_subtruct_piplpimi_factor =(TH1F*)  h1_E_rec_cut2_new->Clone("eRecoEnergy_slice_0");
	TH1F *h_Erec_subtruct_piplpimi_factor =(TH1F*)	h1_E_rec_cut2_new->Clone("h_Erec_subtruct_piplpimi_factor");
	h_Erec_subtruct_piplpimi_factor->Add(h1_E_rec_undetfactor,-1);

//	TH1F *h_Etot_subtruct_piplpimi_factor =(TH1F*)  h1_E_tot_cut2->Clone("epRecoEnergy_slice_0");
	TH1F *h_Etot_subtruct_piplpimi_factor=(TH1F*)  h1_E_tot_cut2->Clone("h_Etot_subtruct_piplpimi_factor");
	h_Etot_subtruct_piplpimi_factor->Add(h1_E_tot_undetfactor,-1);

	TH2F *h2_Erec_pperp_1p1pisub=(TH2F*) h2_Erec_pperp_newcut2->Clone("h2_Erec_pperp_1p1pisub");
	h2_Erec_pperp_1p1pisub->Add(h2_Erec_pperp_1p1pi,-1);

	TH1F *h_Erec_subtruct_piplpimi_factor_fracfeed =(TH1F*)  h1_E_rec_cut2_new_fracfeed->Clone("h_Erec_subtruct_piplpimi_factor_fracfeed");
	h_Erec_subtruct_piplpimi_factor_fracfeed->Add(h1_E_rec_undetfactor_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_factor_fracfeed=(TH1F*)  h1_E_tot_cut2_fracfeed->Clone("h_Etot_subtruct_piplpimi_factor_fracfeed");
	h_Etot_subtruct_piplpimi_factor_fracfeed->Add(h1_E_tot_undetfactor_fracfeed,-1);

	//-----------------------------------undetected 2 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_prot=(TH1F*)	h_Erec_subtruct_piplpimi_factor->Clone("h_Erec_subtruct_piplpimi_prot");
	h_Erec_subtruct_piplpimi_prot->Add(h1_E_rec_p_bkgd,-1);

	TH1F *h_Etot_subtruct_piplpimi_prot=(TH1F*)	h_Etot_subtruct_piplpimi_factor->Clone("h_Etot_subtruct_piplpimi_prot");
	h_Etot_subtruct_piplpimi_prot->Add(h1_E_tot_p_bkgd,-1);

	TH2F *h2_Erec_pperp_2psub=(TH2F*) h2_Erec_pperp_1p1pisub->Clone("h2_Erec_pperp_2psub");
	h2_Erec_pperp_2psub->Add(h2_Erec_pperp_2p,-1);

	TH1F *h_Erec_subtruct_piplpimi_prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_factor_fracfeed->Clone("h_Erec_subtruct_piplpimi_prot_fracfeed");
	h_Erec_subtruct_piplpimi_prot_fracfeed->Add(h1_E_rec_p_bkgd_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_factor_fracfeed->Clone("h_Etot_subtruct_piplpimi_prot_fracfeed");
	h_Etot_subtruct_piplpimi_prot_fracfeed->Add(h1_E_tot_p_bkgd_fracfeed,-1);

	 //-----------------------------------undetected 3 to 2 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_32prot=(TH1F*)	h_Erec_subtruct_piplpimi_prot->Clone("h_Erec_subtruct_piplpimi_32prot");
	h_Erec_subtruct_piplpimi_32prot->Add(h1_E_rec_3pto2p);

	TH1F *h_Etot_subtruct_piplpimi_32prot=(TH1F*)	h_Etot_subtruct_piplpimi_prot->Clone("h_Etot_subtruct_piplpimi_32prot");
	h_Etot_subtruct_piplpimi_32prot->Add(h1_E_tot_3pto2p);

	TH2F *h2_Erec_pperp_32psub=(TH2F*) h2_Erec_pperp_2psub->Clone("h2_Erec_pperp_32psub");
	h2_Erec_pperp_32psub->Add(h2_Erec_pperp_321p);

	TH1F *h_Erec_subtruct_piplpimi_32prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_32prot_fracfeed");
	h_Erec_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_rec_3pto2p_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_32prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_32prot_fracfeed");
	h_Etot_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_tot_3pto2p_fracfeed);

	 //-----------------------------------undetected 3 to 1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_31prot=(TH1F*)	h_Erec_subtruct_piplpimi_32prot->Clone("h_Erec_subtruct_piplpimi_31prot");
	h_Erec_subtruct_piplpimi_31prot->Add(h1_E_rec_3pto1p,-1);

	TH1F *h_Etot_subtruct_piplpimi_31prot=(TH1F*)	h_Etot_subtruct_piplpimi_32prot->Clone("h_Etot_subtruct_piplpimi_31prot");
	h_Etot_subtruct_piplpimi_31prot->Add(h1_E_tot_3pto1p,-1);

	TH2F *h2_Erec_pperp_31psub=(TH2F*) h2_Erec_pperp_32psub->Clone("h2_Erec_pperp_31psub");
	h2_Erec_pperp_31psub->Add(h2_Erec_pperp_31p,-1);

	TH1F *h_Erec_subtruct_piplpimi_31prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_32prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_31prot_fracfeed");
	h_Erec_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_rec_3pto1p_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_31prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_32prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_31prot_fracfeed");
	h_Etot_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_tot_3pto1p_fracfeed,-1);

	 //-----------------------------------undetected 4 to 3->2->1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_43prot=(TH1F*)	h_Erec_subtruct_piplpimi_31prot->Clone("h_Erec_subtruct_piplpimi_43prot");
	h_Erec_subtruct_piplpimi_43prot->Add(h1_E_rec_4pto3p,-1);

	TH1F *h_Etot_subtruct_piplpimi_43prot=(TH1F*)	h_Etot_subtruct_piplpimi_31prot->Clone("h_Etot_subtruct_piplpimi_43prot");
	h_Etot_subtruct_piplpimi_43prot->Add(h1_E_tot_4pto3p,-1);

	TH2F *h2_Erec_pperp_43psub=(TH2F*) h2_Erec_pperp_31psub->Clone("h2_Erec_pperp_43psub");
	h2_Erec_pperp_43psub->Add(h2_Erec_pperp_4321p,-1);

	TH1F *h_Erec_subtruct_piplpimi_43prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_31prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_43prot_fracfeed");
	h_Erec_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_rec_4pto3p_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_43prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_31prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_43prot_fracfeed");
	h_Etot_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_tot_4pto3p_fracfeed,-1);

	 //-----------------------------------undetected 4 to 3->1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_431prot=(TH1F*)	h_Erec_subtruct_piplpimi_43prot->Clone("h_Erec_subtruct_piplpimi_431prot");
	h_Erec_subtruct_piplpimi_431prot->Add(h1_E_rec_43pto1p);

	TH1F *h_Etot_subtruct_piplpimi_431prot=(TH1F*)	h_Etot_subtruct_piplpimi_43prot->Clone("h_Etot_subtruct_piplpimi_431prot");
	h_Etot_subtruct_piplpimi_431prot->Add(h1_E_tot_43pto1p);

	TH2F *h2_Erec_pperp_431psub=(TH2F*) h2_Erec_pperp_43psub->Clone("h2_Erec_pperp_431psub");
	h2_Erec_pperp_431psub->Add(h2_Erec_pperp_431p);

	TH1F *h_Erec_subtruct_piplpimi_431prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_43prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_431prot_fracfeed");
	h_Erec_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_rec_43pto1p_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_431prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_43prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_431prot_fracfeed");
	h_Etot_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_tot_43pto1p_fracfeed);

	//-----------------------------------undetected 4 to 2 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_42prot=(TH1F*)	h_Erec_subtruct_piplpimi_431prot->Clone("h_Erec_subtruct_piplpimi_42prot");
	h_Erec_subtruct_piplpimi_42prot->Add(h1_E_rec_4pto2p);

	TH1F *h_Etot_subtruct_piplpimi_42prot=(TH1F*) h_Etot_subtruct_piplpimi_431prot->Clone("h_Etot_subtruct_piplpimi_42prot");
	h_Etot_subtruct_piplpimi_42prot->Add(h1_E_tot_4pto2p);

	TH2F *h2_Erec_pperp_42psub=(TH2F*) h2_Erec_pperp_431psub->Clone("h2_Erec_pperp_42psub");
	h2_Erec_pperp_42psub->Add(h2_Erec_pperp_421p);

	TH1F *h_Erec_subtruct_piplpimi_42prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_431prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_42prot_fracfeed");
	h_Erec_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_rec_4pto2p_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_42prot_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_431prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_42prot_fracfeed");
	h_Etot_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_tot_4pto2p_fracfeed);

	 //-----------------------------------undetected 4 to 1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_41prot=(TH1F*)	h_Erec_subtruct_piplpimi_42prot->Clone("h_Erec_subtruct_piplpimi_41prot");
	h_Erec_subtruct_piplpimi_41prot->Add(h1_E_rec_4pto1p,-1);

	TH1F *h_Etot_subtruct_piplpimi_41prot=(TH1F*)	h_Etot_subtruct_piplpimi_42prot->Clone("h_Etot_subtruct_piplpimi_41prot");
	h_Etot_subtruct_piplpimi_41prot->Add(h1_E_tot_4pto1p,-1);

	TH2F *h2_Erec_pperp_41psub=(TH2F*) h2_Erec_pperp_42psub->Clone("h2_Erec_pperp_41psub");
	h2_Erec_pperp_41psub->Add(h2_Erec_pperp_41p,-1);

	TH1F *h_Erec_subtruct_piplpimi_41prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_42prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_41prot_fracfeed");
	h_Erec_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_rec_4pto1p_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_41prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_42prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_41prot_fracfeed");
	h_Etot_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_tot_4pto1p_fracfeed,-1);

	//------------------------------------undetected 1p 2pi ->1 p1pi ------ --------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_1p2pi=(TH1F*)	h_Erec_subtruct_piplpimi_41prot->Clone("h_Erec_subtruct_piplpimi_1p2pi");
	h_Erec_subtruct_piplpimi_1p2pi->Add(h1_E_rec_1p2pi);

	TH1F *h_Etot_subtruct_piplpimi_1p2pi=(TH1F*)	h_Etot_subtruct_piplpimi_41prot->Clone("h_Etot_subtruct_piplpimi_1p2pi");
	h_Etot_subtruct_piplpimi_1p2pi->Add(h1_E_tot_1p2pi);

	TH2F *h2_Erec_pperp_sub_1p2pi_1p1pi=(TH2F*) h2_Erec_pperp_41psub->Clone("h2_Erec_pperp_sub_1p2pi_1p1pi");
	h2_Erec_pperp_sub_1p2pi_1p1pi->Add(h2_Erec_pperp_1p2pi_1p1pi);

	TH1F *h_Erec_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_41prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_fracfeed");
	h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Add(h1_E_rec_1p2pi_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_1p2pi_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_41prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_fracfeed");
	h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Add(h1_E_tot_1p2pi_fracfeed);

	//------------------------------------undetected 1p 2pi-> 1p 0pi  ------ --------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*)	h_Erec_subtruct_piplpimi_1p2pi->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi");
	h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Add(h1_E_rec_1p2pi_1p0pi,-1);

	TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi");
	h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Add(h1_E_tot_1p2pi_1p0pi,-1);

	TH2F *h2_Erec_pperp_sub_1p2pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p1pi->Clone("h2_Erec_pperp_sub_1p2pi_1p0pi");
	h2_Erec_pperp_sub_1p2pi_1p0pi->Add(h2_Erec_pperp_1p2pi_1p0pi,-1);

	TH1F *h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
	h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(h1_E_rec_1p2pi_1p0pi_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed");
	h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Add(h1_E_tot_1p2pi_1p0pi_fracfeed,-1);

	//------------------------------------undetected 1p 3pi-> 1p 0pi  ------ --------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_1p3pi=(TH1F*)	h_Erec_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Erec_subtruct_piplpimi_1p3pi");
	h_Erec_subtruct_piplpimi_1p3pi->Add(h1_E_rec_1p3pi);

	TH1F *h_Etot_subtruct_piplpimi_1p3pi=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi->Clone("h_Etot_subtruct_piplpimi_1p3pi");
	h_Etot_subtruct_piplpimi_1p3pi->Add(h1_E_tot_1p3pi);

	TH2F *h2_Erec_pperp_sub_1p3pi=(TH2F*) h2_Erec_pperp_sub_1p2pi_1p0pi->Clone("h2_Erec_pperp_sub_1p3pi");
	h2_Erec_pperp_sub_1p3pi->Add(h2_Erec_pperp_1p3pi);

	TH1F *h_Erec_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_1p3pi_fracfeed");
	h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Add(h1_E_rec_1p3pi_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_1p3pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p2pi_1p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_1p3pi_fracfeed");
	h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Add(h1_E_tot_1p3pi_fracfeed);


	//------------------------------------undetected 2p 2pi -> 1p 0pi  ------ --------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_2p2pi=(TH1F*)	h_Erec_subtruct_piplpimi_1p3pi->Clone("h_Erec_subtruct_piplpimi_2p2pi");
	h_Erec_subtruct_piplpimi_2p2pi->Add(h1_E_rec_2p2pi);

	TH1F *h_Etot_subtruct_piplpimi_2p2pi=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi->Clone("h_Etot_subtruct_piplpimi_2p2pi");
	h_Etot_subtruct_piplpimi_2p2pi->Add(h1_E_tot_2p2pi);

	TH2F *h2_Erec_pperp_sub_2p2pi=(TH2F*) h2_Erec_pperp_sub_1p3pi->Clone("h2_Erec_pperp_sub_2p2pi");
	h2_Erec_pperp_sub_2p2pi->Add(h2_Erec_pperp_2p2pi);

	TH1F *h_Erec_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p2pi_fracfeed");
	h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Add(h1_E_rec_2p2pi_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_2p2pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_1p3pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p2pi_fracfeed");
	h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Add(h1_E_tot_2p2pi_fracfeed);


	//------------------------------------undetected 3p 1pi -> 1p 0pi  ------ --------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_3p1pi=(TH1F*)	h_Erec_subtruct_piplpimi_2p2pi->Clone("h_Erec_subtruct_piplpimi_3p1pi");
	h_Erec_subtruct_piplpimi_3p1pi->Add(h1_E_rec_3p1pi);

	TH1F *h_Etot_subtruct_piplpimi_3p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi->Clone("h_Etot_subtruct_piplpimi_3p1pi");
	h_Etot_subtruct_piplpimi_3p1pi->Add(h1_E_tot_3p1pi);

	TH2F *h2_Erec_pperp_sub_3p1pi=(TH2F*) h2_Erec_pperp_sub_2p2pi->Clone("h2_Erec_pperp_sub_3p1pi");
	h2_Erec_pperp_sub_3p1pi->Add(h2_Erec_pperp_3p1pi);

	TH1F *h_Erec_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_3p1pi_fracfeed");
	h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Add(h1_E_rec_3p1pi_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_3p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p2pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_3p1pi_fracfeed");
	h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Add(h1_E_tot_3p1pi_fracfeed);

	//------------------------------------undetected 2p 1pi -> 2p 0pi  --------------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*)	h_Erec_subtruct_piplpimi_3p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi");
	h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Add(h1_E_rec_2p1pi_2p0pi);

	TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi");
	h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Add(h1_E_tot_2p1pi_2p0pi);

	TH2F *h2_Erec_pperp_sub_2p1pi_2p0pi=(TH2F*) h2_Erec_pperp_sub_3p1pi->Clone("h2_Erec_pperp_sub_2p1pi_2p0pi");
	h2_Erec_pperp_sub_2p1pi_2p0pi->Add(h2_Erec_pperp_2p1pi_2p0pi);

	TH1F *h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
	h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(h1_E_rec_2p1pi_2p0pi_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_3p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed");
	h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Add(h1_E_tot_2p1pi_2p0pi_fracfeed);

//------------------------------------undetected 2p 1pi -> 1p 1pi  ------ --------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*)	h_Erec_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi");
	h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Add(h1_E_rec_2p1pi_1p1pi);

	TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi");
	h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Add(h1_E_tot_2p1pi_1p1pi);

	TH2F *h2_Erec_pperp_sub_2p1pi_1p1pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_2p0pi->Clone("h2_Erec_pperp_sub_2p1pi_1p1pi");
	h2_Erec_pperp_sub_2p1pi_1p1pi->Add(h2_Erec_pperp_2p1pi_1p1pi);

	TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
	h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_rec_2p1pi_1p1pi_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_2p0pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed");
	h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Add(h1_E_tot_2p1pi_1p1pi_fracfeed);

	//------------------------------------undetected 2p 1pi ->1p 0pi  ------ --------------------------------------

//	TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*)	h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi");
	TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*)	h_Erec_subtruct_piplpimi_2p1pi_1p1pi->Clone("eRecoEnergy_slice_0");
	h_Erec_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_rec_2p1pi_1p0pi,-1);

//	TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi");
	TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi->Clone("epRecoEnergy_slice_0");
	h_Etot_subtruct_piplpimi_2p1pi_1p0pi->Add(h1_E_tot_2p1pi_1p0pi,-1);

	TH2F *h2_Erec_pperp_sub_2p1pi_1p0pi=(TH2F*) h2_Erec_pperp_sub_2p1pi_1p1pi->Clone("h2_Erec_pperp_sub_2p1pi_1p0pi");
	h2_Erec_pperp_sub_2p1pi_1p0pi->Add(h2_Erec_pperp_2p1pi_1p0pi,-1);

	TH1F *h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
	h_Erec_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(h1_E_rec_2p1pi_1p0pi_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_2p1pi_1p1pi_fracfeed->Clone("h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed");
	h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed->Add(h1_E_tot_2p1pi_1p0pi_fracfeed,-1);


	 //-----------------------------------looking only at e-, 1pi, undetected pion subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_noprot = (TH1F*)  h1_E_rec_0pi->Clone("h_Erec_subtruct_piplpimi_noprot");
	h_Erec_subtruct_piplpimi_noprot->Add(h1_E_rec_1pi_weight,-1);

	TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed = (TH1F*)  h1_E_rec_0pi_frac_feed->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed");
	h_Erec_subtruct_piplpimi_noprot_frac_feed->Add(h1_E_rec_1pi_weight_frac_feed,-1);
	 //-----------------------------------looking only at e-, 2pi undetected pion subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_noprot_2pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot->Clone("h_Erec_subtruct_piplpimi_noprot_2pi");
	h_Erec_subtruct_piplpimi_noprot_2pi->Add(h1_E_rec_2pi_weight);

	TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed2pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_frac_feed->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed2pi");
	h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Add(h1_E_rec_2pi_weight_frac_feed);

	 //-----------------------------------looking only at e-, 3pi, undetected pion subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_noprot_3pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_2pi->Clone("h_Erec_subtruct_piplpimi_noprot_3pi");
	h_Erec_subtruct_piplpimi_noprot_3pi->Add(h1_E_rec_3pi_weight);

	TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed3pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed3pi");
	h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Add(h1_E_rec_3pi_weight_frac_feed);

	 //-----------------------------------looking only at e-, 4pi, undetected pion subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_noprot_4pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_3pi->Clone("h_Erec_subtruct_piplpimi_noprot_4pi");
	h_Erec_subtruct_piplpimi_noprot_4pi->Add(h1_E_rec_4pi_weight);

	TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed4pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed4pi");
	h_Erec_subtruct_piplpimi_noprot_frac_feed4pi->Add(h1_E_rec_4pi_weight_frac_feed);

	gDirectory->Write("hist_Files", TObject::kOverwrite);
	// skim_tree->AutoSave();

	// --------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Initial # Events = " << fChain->GetEntries() << std::endl;
	std::cout << std::endl << "1e1p0pi Signal # Events = " << SignalEvents << std::endl;
	std::cout << std::endl << "Passing Rate = " << int(double(SignalEvents) / double(fChain->GetEntries())*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "PMiss Fraction 1st Slice = " << int(double(PMiss_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "PMiss Fraction 2nd Slice = " << int(double(PMiss_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "PMiss Fraction 3rd Slice = " << int(double(PMiss_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Events With ECal Within 5\% of ETrue = " << ECalSignalEventsWithin5Perc << std::endl;
	std::cout << std::endl << "ECal 5% Fraction = " << int(double(ECalSignalEventsWithin5Perc) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 1st Slice = " << int(double(ECalSignalEventsWithin5Perc_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 2nd Slice = " << int(double(ECalSignalEventsWithin5Perc_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 3rd Slice = " << int(double(ECalSignalEventsWithin5Perc_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Events With EQE Within 5\% of ETrue = " << EQESignalEventsWithin5Perc << std::endl;
	std::cout << std::endl << "EQE 5% Fraction = " << int(double(EQESignalEventsWithin5Perc) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 1st Slice = " << int(double(EQESignalEventsWithin5Perc_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 2nd Slice = " << int(double(EQESignalEventsWithin5Perc_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 3rd Slice = " << int(double(EQESignalEventsWithin5Perc_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	if (choice == 1) {

		std::cout << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "MEC Fractional Contribution = " << int(double(MECSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "RES Fractional Contribution = " << int(double(RESSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "DIS Fractional Contribution = " << int(double(DISSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;

	}

}

//End Loop function

// -------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------

double systematics::acceptance_c(double p, double cost, double phi, int particle_id,TFile* file_acceptance) {

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
