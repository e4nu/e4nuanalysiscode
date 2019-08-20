#define treeProducer_simulation_cxx
#include "treeProducer_simulation.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

#include "acceptance_c.cpp"

using namespace std;

void treeProducer_simulation::Loop() {

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	TH1D::SetDefaultSumw2();
	gRandom = new TRandom3();

	//___________________________________________________________________________________________________________________________________________________________________________________________

	// Counters

	int ECalCounter = 0, EQECounter = 0;
	int QEEvents = 0, MECEvents = 0, RESEvents = 0, DISEvents = 0;
	//___________________________________________________________________________________________________________________________________________________________________________________________

	// Setting variables

	double Ebeam = -99., fTorusCurrent = -99.;

	// Torus Current & Beam Energy

	if (E == "1_161") { Ebeam = 1.161, fTorusCurrent = 2250.; }
	if (E == "2_261") { Ebeam = 2.261, fTorusCurrent = 2250.; }
	if (E == "4_461") { Ebeam = 4.461, fTorusCurrent = 2250.; }

	// Binding Energy

	double BE = -99.;

	if (Target == "3He") { BE = He3_bind_en - D2_bind_en; }
	if (Target == "4He") { BE = He4_bind_en - H3_bind_en; }
	if (Target == "12C") { BE = C12_bind_en - B_bind_en; }
	if (Target == "56Fe") { BE = Fe_bind_en - Mn_bind_en; }

	double Mp = 0.938; // GeV
	double Mn = 0.939; // GeV
	double Ma = 6 * Mp + 6 * Mn - 0.09216; // GeV
	double MaStar = 10.26;
	
	//___________________________________________________________________________________________________________________________________________________________________________________________

	TVector3 qVector(-99.,-99.,-99.);
	TVector3* Beam = new TVector3(0,0,Ebeam);
	TVector3* PiMinus = new TVector3(-99.,-99.,-99.);
	TVector3* PiPlus = new TVector3(-99.,-99.,-99.);
	TVector3* finalStateNucleon = new TVector3(-99,-99,-99);
	TVector3* FirstFinalStateProton = new TVector3(-99,-99,-99);
	TVector3* SecondFinalStateProton = new TVector3(-99,-99,-99);
	TVector3* ElectronOut = new TVector3(-99,-99,-99); 
	TVector3* ElectronIn = new TVector3(-99,-99,-99); 
	TVector3* ElectronInReco = new TVector3(-99,-99,-99);
	TVector3* FirstElectronInReco = new TVector3(-99,-99,-99);
	TVector3* SecondElectronInReco = new TVector3(-99,-99,-99);
	TLorentzVector* ElectronOutV4 = new TLorentzVector(-99,-99,-99,-99); 
	TLorentzVector* BeamV4 = new TLorentzVector(0,0,Ebeam,Ebeam);

	//___________________________________________________________________________________________________________________________________________________________________________________________

	TString WhichMap = "e2a_maps";
	TFile* file_acceptance = TFile::Open("maps/"+WhichMap+"/"+WhichMap+"_12C_E_"+E+".root");
	TFile* file_acceptance_p = TFile::Open("maps/"+WhichMap+"/"+WhichMap+"_12C_E_"+E+"_p.root");
	TString path = file_name+"_em.root";
	TFile* file = new TFile(path,"recreate");
	std::cout << std::endl << "File " << file_name << "_em.root will be created" << std::endl << std::endl; 
	const double ProtonMass = .9383, MuonMass = .106, NeutronMass = 0.9396, ElectronMass = 0.000511; // [GeV]
	double Weight = 1., fine_struc_const = 0.007297;
	int FirstProtonMissMomentumBin = 99, SecondProtonMissMomentumBin = 99;
	int MissMomentumBin = 99;
	const int n_slice = 3;
	const double pperp_min[n_slice] = {0.,0.2,0.4};
	const double pperp_max[n_slice] = {0.2,0.4,10.};
	
	// __________________________________________________________________________________________________________________________________________________________________________________________

	// Ranges & Binning

	int NBinsEgamma = 25; double MinEgamma = 0., MaxEgamma = 0.25; TString TitleEgamma = ";E_{#gamma} (GeV);";
	int NBinsQ2 = 25; double MinQ2 = 0.4, MaxQ2 = 2.; TString TitleQ2 = ";Q^{2} (GeV^{2}/c^{2});";
	int NBinsxB = 25; double MinxB = 0.7, MaxxB = 1.3; TString TitlexB = ";x_{B};";
	int NBinsnu = 25; double Minnu = 0.2, Maxnu = 1.2; TString Titlenu = ";Energy Transfer (GeV);";
	int NBinsW = 20; double MinW = 0.7, MaxW = 1.2; TString TitleW = ";W (GeV/c^{2});";
	int NBinsPmiss = 25; double MinPmiss = 0., MaxPmiss = 1.; TString TitlePmiss = ";P^{miss}_{#perp} (GeV/c);";
	int NBinsPionMulti = 4; double MinPionMulti = -0.5, MaxPionMulti = 3.5; TString TitlePionMulti = ";Pion Multiplicity;";
	int NBinsReso = 100; double MinReso = -50., MaxReso = 10.; TString TitleReso = ";#frac{E^{cal} - E^{beam}}{E^{beam}} (%);";
	int NBinsResoQE = 100; double MinResoQE = -100., MaxResoQE = 100.; TString TitleResoQE = ";#frac{E^{QE} - E^{beam}}{E^{beam}} (%);";
	int NBinsDeltaPhiT = 18; double MinDeltaPhiT = 0., MaxDeltaPhiT = 80.; TString TitleDeltaPhiT = ";#delta#phi_{T} (degrees);";
	int NBinsDeltaAlphaT = 18; double MinDeltaAlphaT = 0., MaxDeltaAlphaT = 180.; TString TitleDeltaAlphaT = ";#delta#alpha_{T} (degrees);";

	TString TitleQ2Vsnu = ";Energy Transfer (GeV);Q^{2} (GeV^{2}/c^{2});";
	TString TitleQ2VsW = ";W (GeV/c^{2});Q^{2} (GeV^{2}/c^{2});";

	// Electron

	int NBinsElectronEnergyInclusive = 20; double MinElectronEnergyInclusive = 0.4, MaxElectronEnergyInclusive = 2.2; TString TitleElectronEnergyInclusive = ";E_{e'} (GeV);";
	int NBinsElectronEnergy = 20; double MinElectronEnergy = 1., MaxElectronEnergy = 2.1; TString TitleElectronEnergy = ";E_{e'} (GeV);";
	int NBinsElectronPhi = 45; double MinElectronPhi = 0., MaxElectronPhi = 360.; TString TitleElectronPhi = ";#phi_{e'} (degrees);";
	int NBinsElectronTheta = 15; double MinElectronTheta = 15., MaxElectronTheta = 53.; TString TitleElectronTheta = ";#theta_{e'} (degrees);";
	int NBinsElectronCosTheta = 15; double MinElectronCosTheta = 0.6, MaxElectronCosTheta = 1.; TString TitleElectronCosTheta = ";cos(#theta_{e'});";
	int NBinsElectronMom = 40; double MinElectronMom = 1.7, MaxElectronMom = 4.; TString TitleElectronMom = ";P_{e'} (GeV/c);";

	TString TitleElectronPhiVsTheta = ";#phi_{e'} (degrees);#theta_{e'} (degrees);";

	int NBinsThetaAngleRecoNuBeam = 50; double MinThetaAngleRecoNuBeam = 0., MaxThetaAngleRecoNuBeam = 70.; TString TitleThetaAngleRecoNuBeam = ";#theta_{reco #nu,beam} (degrees);";
	int NBinsCosThetaAngleRecoNuBeam = 50; double MinCosThetaAngleRecoNuBeam = 0.7, MaxCosThetaAngleRecoNuBeam = 1.; 
	TString TitleCosThetaAngleRecoNuBeam = ";cos(#theta_{reco #nu,beam}) (degrees);";
	
	// Proton

	int NBinsEp = 20; double MinEp = 0.9, MaxEp = 2.1; TString TitleEp = ";E_{p} (GeV);";
	int NBinsProtonPhi = 45; double MinProtonPhi = 0., MaxProtonPhi = 360.; TString TitleProtonPhi = ";#phi_{p} (degrees);";
	int NBinsProtonTheta = 30; double MinProtonTheta = 10., MaxProtonTheta = 120.; TString TitleProtonTheta = ";#theta_{p} (degrees);";
	int NBinsProtonCosTheta = 15; double MinProtonCosTheta = -1., MaxProtonCosTheta = 1.; TString TitleProtonCosTheta = ";cos(#theta_{p});";

	TString TitleProtonEnergyVsMissMomentum = ";P^{miss}_{#perp} (GeV/c);E_{p} (GeV);";
	TString TitleQ2VsMissMomentum = ";P^{miss}_{#perp} (GeV/c);Q^{2} (GeV^{2}/c^{2});";

	if (E == "2_261" && xBCut == "NoxBCut") {

		MinQ2 = 0.4, MaxQ2 = 2.;
		MinxB = 0., MaxxB = 1.8;
		Minnu = 0.1, Maxnu = 1.8;
		MinW = 0.6, MaxW = 2.;
		MinPmiss = 0., MaxPmiss = 1.;
		MinReso = -50., MaxReso = 10.;

		// Electron

		MinElectronEnergyInclusive = 0.4, MaxElectronEnergyInclusive = 2.2;
		MinElectronEnergy = 0.6, MaxElectronEnergy = 2.2;
		MinElectronPhi = 0., MaxElectronPhi = 360.;
		MinElectronTheta = 15., MaxElectronTheta = 53.;
		MinElectronMom = 1.7, MaxElectronMom = 4.;
	
		// Proton

		MinEp = 0.9, MaxEp = 2.2;
		MinProtonPhi = 0., MaxProtonPhi = 360.;
		MinProtonTheta = 10., MaxProtonTheta = 120.;

	}

	if (E == "4_461") {

		MinQ2 = 0.8, MaxQ2 = 5.5;
		MinxB = 0.7, MaxxB = 1.3;
		Minnu = 0.5, Maxnu =3.2;
		MinW = 0.5, MaxW = 1.5;
		MinPmiss = 0., MaxPmiss = 1.;
		MinReso = -50., MaxReso = 10.;

		// Electron

		MinElectronEnergy = 1.4, MaxElectronEnergy = 4;
		MinElectronPhi = 0., MaxElectronPhi = 360.;
		MinElectronTheta = 18., MaxElectronTheta = 53.;
		MinElectronMom = 1.7, MaxElectronMom = 4.;
	
		// Proton

		MinEp = 0.9, MaxEp = 3.7;
		MinProtonPhi = 0., MaxProtonPhi = 360.;
		MinProtonTheta = 10., MaxProtonTheta = 120.;

	}

	// ____________________________________________________________________________________________________________________________________________________________________________________________

	// Plots for 0 pions

	TH1D* EgammaPlot = new TH1D("EgammaPlot",TitleEgamma,NBinsEgamma,MinEgamma,MaxEgamma);

//	TH1D* Q2Plot = new TH1D("Q2",TitleQ2,NBinsQ2,MinQ2,MaxQ2);
	TH1D* Q2Plot = new TH1D("Q2",TitleQ2,400,0,6);

	TH1D* xBPlot = new TH1D("xB",TitlexB,NBinsxB,MinxB,MaxxB);
	TH1D* nuPlot = new TH1D("nu",Titlenu,NBinsnu,Minnu,Maxnu);

int NBinsnuMariana = 200; double MinnuMariana = 0., MaxnuMariana = 3.5;

TH1D* nuPlot_Q2_0_2 = new TH1D("nu_Q2_0_2",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_0_3 = new TH1D("nu_Q2_0_3",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_0_4 = new TH1D("nu_Q2_0_4",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_0_5 = new TH1D("nu_Q2_0_5",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_0_6 = new TH1D("nu_Q2_0_6",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_1 = new TH1D("nu_Q2_1",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_1_4 = new TH1D("nu_Q2_1_4",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_1_5 = new TH1D("nu_Q2_1_5",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_2_5 = new TH1D("nu_Q2_2_5",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);
TH1D* nuPlot_Q2_3_5 = new TH1D("nu_Q2_3_5",Titlenu,NBinsnuMariana,MinnuMariana,MaxnuMariana);

	TH1D* WPlot = new TH1D("W",TitleW,NBinsW,MinW,MaxW);

//	TH1D* MissMomentum = new TH1D("MissMomentum",TitlePmiss,NBinsPmiss,MinPmiss,MaxPmiss);
	TH1D* MissMomentum = new TH1D("MissMomentum",TitlePmiss,400,0.,1.);

//	TH1D* PionMultiPlot = new TH1D("PionMultiPlot",TitlePionMulti,NBinsPionMulti,MinPionMulti,MaxPionMulti);
	TH1D* PionMultiPlot = new TH1D("PionMultiPlot",TitlePionMulti,10,0,5);

	TH1D* PionMultiQEPlot = new TH1D("PionMultiQEPlot",TitlePionMulti,NBinsPionMulti,MinPionMulti,MaxPionMulti);

//	TH2D* Q2Vsnu = new TH2D("Q2Vsnu",TitleQ2Vsnu,NBinsnu,Minnu,Maxnu,NBinsQ2,MinQ2,MaxQ2);
	TH2D* Q2Vsnu = new TH2D("Q2Vsnu",TitleQ2Vsnu,200,0,3.5,200,0,5);
	TH2D* Q2Vsnu_FirstSector = new TH2D("Q2Vsnu_FirstSector",TitleQ2Vsnu,20,0,3.5,20,0,5);

	TH2D* Q2VsW = new TH2D("Q2VsW",TitleQ2VsW,NBinsW,MinW,MaxW,NBinsQ2,MinQ2,MaxQ2);

	// Electron

	// Inclusive Analysis

	TH1D* ElectronPhiInclusive = new TH1D("ElectronPhiInclusive",TitleElectronPhi,NBinsElectronPhi,MinElectronPhi,MaxElectronPhi);
	TH1D* ElectronThetaInclusive = new TH1D("ElectronThetaInclusive",TitleElectronTheta,NBinsElectronTheta,MinElectronTheta,MaxElectronTheta+20);
	TH1D* ElectronCosThetaInclusive = new TH1D("ElectronCosThetaInclusive",TitleElectronCosTheta,NBinsElectronCosTheta,MinElectronCosTheta,MaxElectronCosTheta);

	TH1D* Q2Incl = new TH1D("Q2Inclusive",TitleQ2,NBinsQ2,MinQ2,MaxQ2);

	double MaxnuInclusive = 1.8;
	TH1D* nuIncl = new TH1D("EePrimeInclusive",Titlenu,NBinsnu,Minnu,MaxnuInclusive);

	// Semi-Inclusive Analysis

	TH1D* ElectronEnergy = new TH1D("EePrime",TitleElectronEnergy,NBinsElectronEnergy,MinElectronEnergy,MaxElectronEnergy);
	TH1D* ElectronPhi = new TH1D("phi_ElectronOut",TitleElectronPhi,NBinsElectronPhi,MinElectronPhi,MaxElectronPhi);
	TH1D* ElectronTheta = new TH1D("theta_ElectronOut",TitleElectronTheta,NBinsElectronTheta,MinElectronTheta,MaxElectronTheta);
	TH1D* ElectronCosTheta = new TH1D("Costheta_ElectronOut",TitleElectronCosTheta,NBinsElectronCosTheta,MinElectronCosTheta,MaxElectronCosTheta);
	TH2D* ElectronPhiTheta = new TH2D("ElectronPhiTheta",TitleElectronPhiVsTheta,NBinsElectronPhi,MinElectronPhi,MaxElectronPhi,NBinsElectronTheta,MinElectronTheta,MaxElectronTheta);

	TH1D* ThetaAngleRecoNuBeamPlot = new TH1D("ThetaAngleRecoNuBeamPlot",TitleThetaAngleRecoNuBeam,NBinsThetaAngleRecoNuBeam,MinThetaAngleRecoNuBeam,MaxThetaAngleRecoNuBeam);
	TH1D* CosThetaAngleRecoNuBeamPlot = new TH1D("CosThetaAngleRecoNuBeamPlot",TitleCosThetaAngleRecoNuBeam,NBinsCosThetaAngleRecoNuBeam,MinCosThetaAngleRecoNuBeam,MaxCosThetaAngleRecoNuBeam);


	// Proton

	TH1D* ProtonEnergy = new TH1D("Ep",TitleEp,NBinsEp,MinEp,MaxEp);
	TH1D* ProtonPhi = new TH1D("phi_finalStateNucleonPlot",TitleProtonPhi,NBinsProtonPhi,MinProtonPhi,MaxProtonPhi);
	TH1D* ProtonTheta = new TH1D("theta_finalStateNucleon",TitleProtonTheta,NBinsProtonTheta,MinProtonTheta,MaxProtonTheta);
	TH1D* ProtonCosTheta = new TH1D("costheta_finalStateNucleon",TitleProtonCosTheta,NBinsProtonCosTheta,MinProtonCosTheta,MaxProtonCosTheta);
	TH1D* EcalReso = new TH1D("EcalReso",TitleReso,NBinsReso,MinReso,MaxReso);
	TH1D* EQEReso = new TH1D("EQEReso",TitleResoQE,NBinsResoQE,MinResoQE,MaxResoQE);

	TH2D* EpVsMissMomentum = new TH2D("EpVsMissMomentum",TitleProtonEnergyVsMissMomentum,NBinsPmiss,MinPmiss,MaxPmiss,NBinsEp,MinEp,MaxEp);
	TH2D* Q2VsMissMomentum = new TH2D("Q2VsMissMomentum",TitleQ2VsMissMomentum,NBinsPmiss,MinPmiss,MaxPmiss,NBinsQ2,MinQ2,MaxQ2);

	TH1D* DeltaPhiTPlot = new TH1D("DeltaPhiTPlot",TitleDeltaPhiT,NBinsDeltaPhiT,MinDeltaPhiT,MaxDeltaPhiT);
	TH1D* DeltaAlphaTPlot = new TH1D("DeltaAlphaTPlot",TitleDeltaAlphaT,NBinsDeltaAlphaT,MinDeltaAlphaT,MaxDeltaAlphaT);

	// ____________________________________________________________________________________________________________________________________________________________________________________

	int NBinsCalorimetricEnergy = 400; double MinCalorimetricEnergy = 0.,MaxCalorimetricEnergy = 6.; TString TitleCalorimetricEnergy = ";E^{cal} (GeV)";
	int NBinsLeptonicEnergy = 80; double MinLeptonicEnergy = 0.,MaxLeptonicEnergy = 6.; TString TitleLeptonicEnergy = ";E^{QE} (GeV)";

/*	TH1D* epRecoEnergy_slice_0 = new TH1D("epRecoEnergy_slice_0",TitleCalorimetricEnergy,NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);
	TH1D* epRecoEnergy_slice_1 = new TH1D("epRecoEnergy_slice_1",TitleCalorimetricEnergy,NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);
	TH1D* epRecoEnergy_slice_2 = new TH1D("epRecoEnergy_slice_2",TitleCalorimetricEnergy,NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);
	TH1D* epRecoEnergy_slice_3 = new TH1D("epRecoEnergy_slice_3",TitleCalorimetricEnergy,NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);

	TH1D* eRecoEnergy_slice_0 = new TH1D("eRecoEnergy_slice_0",TitleLeptonicEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy);
	TH1D* eRecoEnergy_slice_1 = new TH1D("eRecoEnergy_slice_1",TitleLeptonicEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy);
	TH1D* eRecoEnergy_slice_2 = new TH1D("eRecoEnergy_slice_2",TitleLeptonicEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy);
	TH1D* eRecoEnergy_slice_3 = new TH1D("eRecoEnergy_slice_3",TitleLeptonicEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy);*/

	TH2D* ECalvsEQE = new TH2D("ECalvsEQE",TitleLeptonicEnergy+TitleCalorimetricEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy,
					       NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);
	TH2D* QEECalvsEQE = new TH2D("QEECalvsEQE",TitleLeptonicEnergy+TitleCalorimetricEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy,
					       NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);
	TH2D* MECECalvsEQE = new TH2D("MECECalvsEQE",TitleLeptonicEnergy+TitleCalorimetricEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy,
					       NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);
	TH2D* RESECalvsEQE = new TH2D("RESECalvsEQE",TitleLeptonicEnergy+TitleCalorimetricEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy,
					       NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);
	TH2D* DISECalvsEQE = new TH2D("DISECalvsEQE",TitleLeptonicEnergy+TitleCalorimetricEnergy,NBinsLeptonicEnergy,MinLeptonicEnergy,MaxLeptonicEnergy,
					       NBinsCalorimetricEnergy,MinCalorimetricEnergy,MaxCalorimetricEnergy);

	// ____________________________________________________________________________________________________________________________________________________________________________________

	int n_bins;
	double *x_values = 0;

	if (E == "1_161") {

		n_bins = 38;
		x_values = new double[n_bins+1];
		for (int i=0;i<=17;i++) { x_values[i]=0.4+i*0.04; }
		for (int i=0;i<=20;i++) { x_values[i+18]=1.08+(i+1)*0.02; }
	}

	if (E == "2_261") {

		n_bins = 54;
		x_values = new double[n_bins+1];
		for (int i=0;i<=23;i++) { x_values[i]=i*0.09; }
		for (int i=0;i<=30;i++) { x_values[i+24]=2.07+(i+1)*0.03; }
	}

	if (E == "4_461") {

		n_bins = 38;
		x_values = new double[n_bins+1];
		for (int i=0;i<=21;i++) { x_values[i]=i*0.2; }
		for (int i=0;i<=16;i++) { x_values[i+22]=4.2+(i+1)*0.05; }
	}

	//_______________________________________________________________________________________________________________________________________________________________________________________

	TH1D* epRecoEnergy_slice_0 = new TH1D("epRecoEnergy_slice_0",TitleCalorimetricEnergy,n_bins,x_values);
	TH1D* epRecoEnergy_slice_1 = new TH1D("epRecoEnergy_slice_1",TitleCalorimetricEnergy,n_bins,x_values);
	TH1D* epRecoEnergy_slice_2 = new TH1D("epRecoEnergy_slice_2",TitleCalorimetricEnergy,n_bins,x_values);
	TH1D* epRecoEnergy_slice_3 = new TH1D("epRecoEnergy_slice_3",TitleCalorimetricEnergy,n_bins,x_values);

	TH1D* eRecoEnergy_slice_0 = new TH1D("eRecoEnergy_slice_0",TitleLeptonicEnergy,n_bins,x_values);
	TH1D* eRecoEnergy_slice_1 = new TH1D("eRecoEnergy_slice_1",TitleLeptonicEnergy,n_bins,x_values);
	TH1D* eRecoEnergy_slice_2 = new TH1D("eRecoEnergy_slice_2",TitleLeptonicEnergy,n_bins,x_values);
	TH1D* eRecoEnergy_slice_3 = new TH1D("eRecoEnergy_slice_3",TitleLeptonicEnergy,n_bins,x_values);

	//_______________________________________________________________________________________________________________________________________________________________________________________

	int countEvents = 0; Long64_t nbytes = 0, nb = 0;
	int EventsWithGammas = 0, EventsWithNeutralPions = 0, EventsWith1Gamma = 0, EventsWith2Gammas = 0 ;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		float cphil = 0, cphir = 0;

		//____________________________________________________________________________________________________________________________________________________________________________________

		int ProtonTagging = 0, ChargedPionPlusTagging = 0, ChargedPionMinusTagging = 0, GammaTagging = 0, RadGammaTagging = 0, NeutralPionTagging = 0;
		int PiPlusFinalStateIndex = 99, PiMinusFinalStateIndex = 99;
		vector <int> ProtonID; ProtonID.clear();
		vector <int> GammaID; GammaID.clear();

		for (int i = 0; i < nf; i++) {

			if (pdgf[i] == 2212 && pf[i] > 0.3) {

				ProtonTagging ++;
				ProtonID.push_back(i);

			}

			if (pdgf[i] == 211 && pf[i] > 0.15)  {

				ChargedPionPlusTagging ++;

			}

			if (pdgf[i] == -211 && pf[i] > 0.15)  {

				ChargedPionMinusTagging ++;

			}

			if (pdgf[i] == 22 && pf[i] < 0.3) { GammaTagging ++; GammaID.push_back(i); }
			if (pdgf[i] == 22 && pf[i] > 0.3) { RadGammaTagging ++; }
			if (pdgf[i] == 111) { NeutralPionTagging ++; }

		}

		//_____________________________________________________________________________________________________________________________________________________________________________________

		// Define the 3-vectors and the necessary variables for the particles involved

		// Outcoming e' (electron after scattering)  
    
		ElectronOut->SetXYZ(pxl,pyl,pzl);

		//double TestReso = 0.02;
		//double reso_e = TestReso; // smearing for the electron
		double reso_e = 0.005; // smearing for the electron
		double SmearedPe = gRandom->Gaus(pl,reso_e*pl);
		double SmearedEe = sqrt( SmearedPe*SmearedPe + ElectronMass * ElectronMass );
		ElectronOut->SetXYZ(SmearedPe/pl * pxl,SmearedPe/pl * pyl,SmearedPe/pl * pzl);

		ElectronOutV4->SetPxPyPzE(ElectronOut->X(),ElectronOut->Y(),ElectronOut->Z(),SmearedEe);
		if (!ElectronFiducialCut(ElectronOut) ) continue; // Electron theta & phi fiducial cuts
		double pElectronOut = ElectronOut->Mag(); 
		double theta_ElectronOut = ElectronOut->Theta();
		double phi_ElectronOut = ElectronOut->Phi()+TMath::Pi();

		// q vector  
    
		qVector.SetXYZ(Beam->X()-ElectronOut->X(),Beam->Y()-ElectronOut->Y(),Beam->Z()-ElectronOut->Z());
		double nu = Ev - El;
		TLorentzVector V4_q(qVector,nu);

//		double Q2reco = Q2;
		double Q2reco = -V4_q.Mag2();
		double xB = Q2reco / (2*ProtonMass*nu);

		// Leptonic Energy Reconstruction

		double ERecoOnlyePrime = (2*ProtonMass*BE + 2*ProtonMass*El - pow(ElectronMass,2)) / 2 / (ProtonMass - El + pElectronOut*cos(theta_ElectronOut)); 

		// Define the weight

	        // Electron

		double e_acc_ratio = acceptance_c(pElectronOut, cos(theta_ElectronOut), phi_ElectronOut, 11,file_acceptance);
		if (phi_ElectronOut < TMath::Pi()) {phi_ElectronOut += TMath::Pi();}
		else {phi_ElectronOut -= TMath::Pi();}

		// Normalize by the Mott xsection for the difference in the masses of the propagators

		double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(theta_ElectronOut)+1))/(2*pow(El,2)*pow((1-cos(theta_ElectronOut)),2.));
		if (CFGM == "_neutrino") { Mott_cross_sec = 1.; };

		double WeightIncl = wght*e_acc_ratio / Mott_cross_sec;

		//_____________________________________________________________________________________________________________________________________________________________________________________

		ElectronThetaInclusive->Fill(theta_ElectronOut*180/TMath::Pi(),WeightIncl);
		ElectronCosThetaInclusive->Fill( cos(theta_ElectronOut) ,WeightIncl);
		ElectronPhiInclusive->Fill(phi_ElectronOut*180/TMath::Pi(),WeightIncl);

		Q2Incl->Fill(Q2reco,WeightIncl);
		nuIncl->Fill(nu,WeightIncl);

		double PionMulti = ChargedPionPlusTagging + ChargedPionMinusTagging;
		// Pion Multiplicity
		if (Q2reco > 0.5) {

			if (xBCut == "NoxBCut") {
				PionMultiPlot->Fill(PionMulti,WeightIncl);
				if (ProtonTagging == 1) { 

					finalStateNucleon->SetXYZ(pxf[ProtonID[0]],pyf[ProtonID[0]],pzf[ProtonID[0]]);
					if (ProtonFiducialCut(finalStateNucleon, &cphil, &cphir) ) { // Proton theta & phi fiducial cuts
						double E_finalStateNucleon = Ef[ProtonID[0]]; 
						double pfinalStateNucleon  = finalStateNucleon->Mag();
						double phi_finalStateNucleon = finalStateNucleon->Phi() + TMath::Pi(); 
						double theta_finalStateNucleon = finalStateNucleon->Theta();
						double p_acc_ratio = acceptance_c(pfinalStateNucleon, cos(theta_finalStateNucleon), phi_finalStateNucleon, 2212,file_acceptance_p);

						PionMultiQEPlot->Fill(PionMulti,e_acc_ratio * p_acc_ratio / Mott_cross_sec); } 
					}
			}

			if (xBCut == "xBCut") {

				if (fabs(xB-1.) < 0.2) { 
					PionMultiPlot->Fill(PionMulti,WeightIncl);
					if (ProtonTagging == 1) { 

					finalStateNucleon->SetXYZ(pxf[ProtonID[0]],pyf[ProtonID[0]],pzf[ProtonID[0]]);
					if (!ProtonFiducialCut(finalStateNucleon, &cphil, &cphir) ) {  // Proton theta & phi fiducial cuts
						double E_finalStateNucleon = Ef[ProtonID[0]]; 
						double pfinalStateNucleon  = finalStateNucleon->Mag();
						double phi_finalStateNucleon = finalStateNucleon->Phi() + TMath::Pi(); 
						double theta_finalStateNucleon = finalStateNucleon->Theta();
						double p_acc_ratio = acceptance_c(pfinalStateNucleon, cos(theta_finalStateNucleon), phi_finalStateNucleon, 2212,file_acceptance_p);

						PionMultiQEPlot->Fill(PionMulti,e_acc_ratio * p_acc_ratio / Mott_cross_sec); } 

					}
				}
			}
		}
		//_____________________________________________________________________________________________________________________________________________________________________________________

		// Apply the Selection Cuts

		if (W > 2) continue;
		if (Ebeam == 4.461) { if (Q2reco < 0.8) continue; } 
		if (Ebeam == 2.261) { if (Q2reco < 0.4) continue; }
		if (Ebeam == 1.161) { if (Q2reco < 0.1) continue; }
		if (ProtonTagging != 1) {continue;}
//		if (GammaTagging != 0) { continue;} 
		if (RadGammaTagging != 0) { continue; }
		if (xBCut == "xBCut") { if (fabs(xB-1.) > 0.2) { continue;} }

		//_____________________________________________________________________________________________________________________________________________________________________________________

		//if (nfp != 1) {continue;}
		//if (nfpip + nfpim + nfpi0 != 0) {continue;}

		if (ChargedPionPlusTagging + ChargedPionMinusTagging != 0) {continue;}

		//____________________________________________________________________________________________________________________________________________________________________________________

		// 1 proton analysis

		// Hit nucleon in final state

		finalStateNucleon->SetXYZ(pxf[ProtonID[0]],pyf[ProtonID[0]],pzf[ProtonID[0]]);

		//double reso_p = TestReso; // smearing for the proton
		double reso_p = 0.01; // smearing for the proton
		//double reso_p = 0.; // smearing for the proton
		double SmearedPp = gRandom->Gaus(pf[ProtonID[0]],reso_p*pf[ProtonID[0]]);
		double SmearedEp = sqrt( SmearedPp*SmearedPp + ProtonMass * ProtonMass );
		finalStateNucleon->SetXYZ(SmearedPp/pf[ProtonID[0]] * pxf[ProtonID[0]],SmearedPp/pf[ProtonID[0]] * pyf[ProtonID[0]],SmearedPp/pf[ProtonID[0]] * pzf[ProtonID[0]]);

		if (!ProtonFiducialCut(finalStateNucleon, &cphil, &cphir) ) { continue; } // Proton theta & phi fiducial cuts
//		double E_finalStateNucleon = Ef[ProtonID[0]]; 
		double E_finalStateNucleon = SmearedEp; 
		double pfinalStateNucleon  = finalStateNucleon->Mag();
		double phi_finalStateNucleon = finalStateNucleon->Phi() + TMath::Pi(); 
		double theta_finalStateNucleon = finalStateNucleon->Theta();

		// Reconstructed Incoming Electron 

		ElectronInReco->SetXYZ(ElectronOut->X() + finalStateNucleon->X(), ElectronOut->Y() + finalStateNucleon->Y(), ElectronOut->Z() + finalStateNucleon->Z());
		double miss_momentum = sqrt(pow(ElectronInReco->X(),2)+pow(ElectronInReco->Y(),2)); 

		// Calorimetric Energy Reconstruction

		double Eecal = El + E_finalStateNucleon - ProtonMass + BE;

		// Transverse variables

		TVector3 ProtonT(pxf[ProtonID[0]],pyf[ProtonID[0]],0);
		double ProtonTMag = ProtonT.Mag();
		TVector3 MinusElectronT(-ElectronOut->X(),-ElectronOut->Y(),0);
		double MinusElectronTMag = MinusElectronT.Mag();
		TVector3 MissMomentumT(ElectronInReco->X(),ElectronInReco->Y(),0);
		double MissMomentumTMag = MissMomentumT.Mag();

		double DeltaPhiT = TMath::ACos(ProtonT*MinusElectronT/ProtonTMag/MinusElectronTMag)*180/TMath::Pi();
		double DeltaAlphaT = TMath::ACos(MissMomentumT*MinusElectronT/MissMomentumTMag/MinusElectronTMag)*180/TMath::Pi();

		// Reco neutrino angle with respect to the beam 

		double ThetaAngleRecoNuBeam = (*Beam).Angle(*ElectronInReco);
		double ThetaAngleRecoNuBeam_Deg = ThetaAngleRecoNuBeam * 180. / TMath::Pi();	
		double 	CosThetaAngleRecoNuBeam = cos(ThetaAngleRecoNuBeam);


		//_________________________________________________________________________________________________________________________________________________________________________________

		// Define the proton weight

        	// Proton

		double p_acc_ratio = acceptance_c(pfinalStateNucleon, cos(theta_finalStateNucleon), phi_finalStateNucleon, 2212,file_acceptance_p);
		if ( FSIModel == "Data" ) { p_acc_ratio = 1.; }
		if (phi_finalStateNucleon < TMath::Pi()) {phi_finalStateNucleon += TMath::Pi();}
		else {phi_finalStateNucleon -= TMath::Pi();}

		Weight = e_acc_ratio * p_acc_ratio / Mott_cross_sec;

        	if (fabs(Weight) != Weight) continue;

		Weight = wght*Weight;

		//_____________________________________________________________________________________________________________________________________________________________________________________

		if(GammaTagging!=0) {EventsWithGammas++;} 
		if(GammaTagging==1) {EventsWith1Gamma++;} 
		if(GammaTagging==2) {EventsWith2Gammas++;} 
		if(NeutralPionTagging!=0) {EventsWithNeutralPions++;}

		//_____________________________________________________________________________________________________________________________________________________________________________________

		if (qel==1) { QEEvents++; }
		if (mec==1) { MECEvents++; }
		if (res==1) { RESEvents++; }
		if (dis==1) { DISEvents++; }

		//_____________________________________________________________________________________________________________________________________________________________________________________

		countEvents ++;	// Increase the number of the events that pass the cuts by one   

		//_____________________________________________________________________________________________________________________________________________________________________________________

		// Define the missing momentum bin

		if (miss_momentum < 0.2) MissMomentumBin = 1;
		if (miss_momentum > 0.2 && miss_momentum < 0.4) MissMomentumBin = 2;
		if (miss_momentum > 0.4) MissMomentumBin = 3;

		//_____________________________________________________________________________________________________________________________________________________________________________

		// Demand that we have 1 proton and 0 charged pions
		// This is the place where we fill our main plots

		int NGammas = (int)GammaID.size();
		for (int WhichPhoton = 0; WhichPhoton < NGammas; WhichPhoton++ ) {

			EgammaPlot->Fill(Ef[GammaID[WhichPhoton]],Weight);
			
		}

		Q2Plot->Fill(Q2reco,Weight);
		nuPlot->Fill(nu,Weight);

if ( abs(Q2reco - 0.2) < 0.05 ) { nuPlot_Q2_0_2->Fill(nu,Weight); }
if ( abs(Q2reco - 0.3) < 0.05 ) { nuPlot_Q2_0_3->Fill(nu,Weight); }
if ( abs(Q2reco - 0.4) < 0.05 ) { nuPlot_Q2_0_4->Fill(nu,Weight); }
if ( abs(Q2reco - 0.5) < 0.05 ) { nuPlot_Q2_0_5->Fill(nu,Weight); }
if ( abs(Q2reco - 0.6) < 0.05 ) { nuPlot_Q2_0_6->Fill(nu,Weight); }
if ( abs(Q2reco - 1.0) < 0.05 ) { nuPlot_Q2_1->Fill(nu,Weight); }
if ( abs(Q2reco - 1.4) < 0.05 ) { nuPlot_Q2_1_4->Fill(nu,Weight); }
if ( abs(Q2reco - 1.5) < 0.05 ) { nuPlot_Q2_1_5->Fill(nu,Weight); }
if ( abs(Q2reco - 2.5) < 0.05 ) { nuPlot_Q2_2_5->Fill(nu,Weight); }
if ( abs(Q2reco - 3.5) < 0.05 ) { nuPlot_Q2_3_5->Fill(nu,Weight); }

		WPlot->Fill(W,Weight);
		xBPlot->Fill(xB,Weight);
		MissMomentum->Fill(miss_momentum,Weight);

		// 1D Electron Plots

		ElectronTheta->Fill(theta_ElectronOut*180/TMath::Pi(),Weight);
		ElectronCosTheta->Fill( cos(theta_ElectronOut) ,Weight);
		ElectronPhi->Fill(phi_ElectronOut*180/TMath::Pi(),Weight);
		ElectronEnergy->Fill(El,Weight);

		ElectronPhiTheta->Fill(phi_ElectronOut*180/TMath::Pi(),theta_ElectronOut*180/TMath::Pi(),Weight);

		ThetaAngleRecoNuBeamPlot->Fill(ThetaAngleRecoNuBeam_Deg,Weight);
		CosThetaAngleRecoNuBeamPlot->Fill(CosThetaAngleRecoNuBeam,Weight);

		// 1D Proton Plots		

		ProtonCosTheta->Fill( cos(theta_finalStateNucleon) ,Weight);	
		ProtonTheta->Fill(theta_finalStateNucleon*180/TMath::Pi(),Weight);
		ProtonPhi->Fill(phi_finalStateNucleon*180/TMath::Pi(),Weight);
		ProtonEnergy->Fill(E_finalStateNucleon,Weight);
		EcalReso->Fill((Eecal - Ebeam) / Ebeam * 100.,Weight);
		EQEReso->Fill((ERecoOnlyePrime - Ebeam) / Ebeam * 100.,Weight);

		// 2D Plots

		Q2VsW->Fill(W,Q2reco,Weight);
		Q2Vsnu->Fill(nu,Q2reco,Weight);
if ( phi_ElectronOut*180/TMath::Pi() > 0 && phi_ElectronOut*180/TMath::Pi() < 60) { Q2Vsnu_FirstSector->Fill(nu,Q2reco,Weight); }
		EpVsMissMomentum->Fill(miss_momentum,E_finalStateNucleon,Weight);
		Q2VsMissMomentum->Fill(miss_momentum,Q2reco,Weight);

		ECalvsEQE->Fill(ERecoOnlyePrime,Eecal,Weight);
		if (qel==1) { QEECalvsEQE->Fill(ERecoOnlyePrime,Eecal,Weight); }
		if (mec==1) { MECECalvsEQE->Fill(ERecoOnlyePrime,Eecal,Weight); }
		if (res==1) { RESECalvsEQE->Fill(ERecoOnlyePrime,Eecal,Weight); }
		if (dis==1) { DISECalvsEQE->Fill(ERecoOnlyePrime,Eecal,Weight); }

		DeltaPhiTPlot->Fill(DeltaPhiT,Weight);
		DeltaAlphaTPlot->Fill(DeltaAlphaT,Weight);

		if ( fabs(Eecal - Ebeam) / Ebeam *100. < 5 ) { ECalCounter ++; }
		if ( fabs(ERecoOnlyePrime - Ebeam) / Ebeam *100. < 5 ) { EQECounter ++; }

		// Reconstructed Energy Plots

		epRecoEnergy_slice_0->Fill(Eecal,Weight); eRecoEnergy_slice_0->Fill(ERecoOnlyePrime,Weight);
		if (MissMomentumBin == 1) { epRecoEnergy_slice_1->Fill(Eecal,Weight); eRecoEnergy_slice_1->Fill(ERecoOnlyePrime,Weight);}
		if (MissMomentumBin == 2) { epRecoEnergy_slice_2->Fill(Eecal,Weight); eRecoEnergy_slice_2->Fill(ERecoOnlyePrime,Weight);}
		if (MissMomentumBin == 3) { epRecoEnergy_slice_3->Fill(Eecal,Weight); eRecoEnergy_slice_3->Fill(ERecoOnlyePrime,Weight);}

	} // end of the loop over events

	// __________________________________________________________________________________________________________________________________________________________________________________________

	// Print this message after the loop over the events

	std::cout << std::endl;
	std::cout << "File " << path <<" created " << std::endl; 
	std::cout << std::endl;
	std::cout << "Efficiency = " << double(countEvents)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << std::endl << "Out of " << countEvents << " selected events, " <<  EventsWithGammas << " of them contain gammas (1 gamma = " << EventsWith1Gamma << " events, 2 gammas = " 
	<< EventsWith2Gammas << ") and " << EventsWithNeutralPions << "  of them contain neutral pions" << std::endl << std::endl;

	cout << " ECalCounter = " << ECalCounter << "   ECalCounter / countEvents (%) = " << double(ECalCounter) / double(countEvents) * 100. << std::endl;
	cout << " EQECounter = " << EQECounter << "   EQECounter / countEvents (%) = " << double(EQECounter) / double(countEvents) * 100. << std::endl << std::endl;

	std::cout << "QEEvents = " << QEEvents << "  (" << double(QEEvents)/double(countEvents)*100. << "%)" << std::endl;
	std::cout << "MECEvents = " << MECEvents << "  (" << double(MECEvents)/double(countEvents)*100. << "%)" << std::endl;
	std::cout << "RESEvents = " << RESEvents << "  (" << double(RESEvents)/double(countEvents)*100. << "%)" << std::endl;
	std::cout << "DISEvents = " << DISEvents << "  (" << double(DISEvents)/double(countEvents)*100. << "%)" << std::endl << std::endl;

	std::cout << "Total events = " << countEvents << std::endl << std::endl;

	file->Write(); file->Close(); 

} // End of the program
