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
#include "TCutG.h"
#include "TNtuple.h"
#include <TROOT.h>
#include "TSystem.h"
#include "TH3D.h"

#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>

#include "acceptance_c.cpp"

#include "Constants.h"
#include "./Secondary_Code/ToString.cpp"

using namespace std;
using namespace Constants;

void treeProducer_simulation::Loop() {

	TH1D::SetDefaultSumw2();

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	// ------------------------------------------------------------------------------------------------------------------------------------------

//	TFile* acc_e = new TFile("./AcceptanceMap/AcceptanceMaps_e.root" );
//	TNtuple* n_e = (TNtuple*)(acc_e->Get("ntuple"));
//	TFile* acc_p = new TFile("./AcceptanceMap/AcceptanceMaps_p.root" );
//	TNtuple* n_p = (TNtuple*)(acc_p->Get("ntuple"));

	TFile* acc_e = new TFile("./AcceptanceMap/AcceptanceMap_e_TH3D.root" );
	TH3D* h3_e = (TH3D*)(acc_e->Get("h3"));
	TFile* acc_p = new TFile("./AcceptanceMap/AcceptanceMap_p_TH3D.root" );
	TH3D* h3_p = (TH3D*)(acc_p->Get("h3"));

        // Loading Graphical Cuts

        TFile * graph_cuts1 = new TFile("Rey/newest_cuts.root");
        TCutG * cut_pL_pR1  = (TCutG*) graph_cuts1 -> Get("cut_pL_pR" );

	TFile* file = new TFile(file_name+"_em.root","recreate");
	std::cout << std::endl << "File " << file_name << ".root will be created" << std::endl << std::endl; 
	double Weight = 1.;
	double Ebeam = 4.325;
	double Pbeam = 4.325;

	// -------------------------------------------------------------------------------------------------------------------------

	// General Plots

	TH1D* Q2Plot = new TH1D("Q2",Constants::TitleQ2,Constants::NBinsQ2,Constants::MinQ2,Constants::MaxQ2);
	TH1D* xBPlot = new TH1D("xB",Constants::TitlexB,Constants::NBinsxB,Constants::MinxB,Constants::MaxxB);
	TH1D* nuPlot = new TH1D("nu",Constants::Titlenu,Constants::NBinsnu,Constants::Minnu,Constants::Maxnu);
	TH1D* WPlot = new TH1D("W",Constants::TitleW,Constants::NBinsW,Constants::MinW,Constants::MaxW);
	TH1D* MissMomentum = new TH1D("MissMomentum",Constants::TitlePmiss,Constants::NBinsPmiss,Constants::MinPmiss,Constants::MaxPmiss);

	TH1D* MissMomentumX = new TH1D("MissMomentumX",Constants::TitlePmissX,Constants::NBinsPmissX,Constants::MinPmissX,Constants::MaxPmissX);
	TH1D* MissMomentumY = new TH1D("MissMomentumY",Constants::TitlePmissY,Constants::NBinsPmissY,Constants::MinPmissY,Constants::MaxPmissY);
	TH1D* MissMomentumZ = new TH1D("MissMomentumZ",Constants::TitlePmissZ,Constants::NBinsPmissZ,Constants::MinPmissZ,Constants::MaxPmissZ);
	TH1D* TotalMissMomentum = new TH1D("TotalMissMomentum",Constants::TitleTotalPmiss,Constants::NBinsTotalPmiss,Constants::MinTotalPmiss,Constants::MaxTotalPmiss);

//	TH2D* Q2Vsnu = new TH2D("Q2Vsnu",Constants::TitleQ2Vsnu,Constants::NBinsnu,Constants::Minnu,Constants::Maxnu,Constants::NBinsQ2,Constants::MinQ2,Constants::MaxQ2);
//	TH2D* Q2VsW = new TH2D("Q2VsW",Constants::TitleQ2VsW,Constants::NBinsW,Constants::MinW,Constants::MaxW,Constants::NBinsQ2,Constants::MinQ2,Constants::MaxQ2);

	// Electron

	// Semi-Inclusive Analysis

	TH1D* ElectronEnergy = new TH1D("ElectronEnergy",Constants::TitleElectronEnergy,Constants::NBinsElectronEnergy,Constants::MinElectronEnergy,Constants::MaxElectronEnergy);
//	TH1D* ElectronPhi = new TH1D("phi_ElectronOut",Constants::TitleElectronPhi,Constants::NBinsElectronPhi,Constants::MinElectronPhi,Constants::MaxElectronPhi);
	TH1D* ElectronCosTheta = new TH1D("Costheta_ElectronOut",Constants::TitleElectronCosTheta,Constants::NBinsElectronCosTheta,Constants::MinElectronCosTheta,Constants::MaxElectronCosTheta);

	TH1D* ElectronMomX = new TH1D("ElectronMomX",TitleElectronMomX,NBinsElectronMomX,MinElectronMomX,MaxElectronMomX); 
	TH1D* ElectronMomY = new TH1D("ElectronMomY",TitleElectronMomY,NBinsElectronMomY,MinElectronMomY,MaxElectronMomY); 
	TH1D* ElectronMomZ = new TH1D("ElectronMomZ",TitleElectronMomZ,NBinsElectronMomZ,MinElectronMomZ,MaxElectronMomZ); 

	// Proton

	TH1D* ProtonEnergy = new TH1D("Ep",Constants::TitleEp,Constants::NBinsEp,Constants::MinEp,Constants::MaxEp);
//	TH1D* ProtonPhi = new TH1D("phi_finalStateNucleonPlot",Constants::TitleProtonPhi,Constants::NBinsProtonPhi,Constants::MinProtonPhi,Constants::MaxProtonPhi);
	TH1D* ProtonCosTheta = new TH1D("costheta_finalStateNucleon",Constants::TitleProtonCosTheta,Constants::NBinsProtonCosTheta,Constants::MinProtonCosTheta,Constants::MaxProtonCosTheta);

	TH1D* ProtonMomX = new TH1D("ProtonMomX",TitleProtonMomX,NBinsProtonMomX,MinProtonMomX,MaxProtonMomX); 
	TH1D* ProtonMomY = new TH1D("ProtonMomY",TitleProtonMomY,NBinsProtonMomY,MinProtonMomY,MaxProtonMomY); 
	TH1D* ProtonMomZ = new TH1D("ProtonMomZ",TitleProtonMomZ,NBinsProtonMomZ,MinProtonMomZ,MaxProtonMomZ); 

//	TH2D* EpVsMissMomentum = new TH2D("EpVsMissMomentum",
//		Constants::TitleProtonEnergyVsMissMomentum,Constants::NBinsPmiss,Constants::MinPmiss,Constants::MaxPmiss,Constants::NBinsEp,Constants::MinEp,Constants::MaxEp);
//	TH2D* Q2VsMissMomentum = new TH2D("Q2VsMissMomentum",
//		Constants::TitleQ2VsMissMomentum,Constants::NBinsPmiss,Constants::MinPmiss,Constants::MaxPmiss,Constants::NBinsQ2,Constants::MinQ2,Constants::MaxQ2);

//	TH2D* ElectronCosThetaVsMissMomentum = new TH2D("ElectronCosThetaVsMissMomentum",
//		Constants::TitleElectronCosThetaVsMissMomentum,
//		Constants::NBinsPmiss,Constants::MinPmiss,Constants::MaxPmiss,Constants::NBinsElectronCosTheta,Constants::MinElectronCosTheta,Constants::MaxElectronCosTheta);

	// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Reconstructed Energy


	TH1D* ECalRecoPlot = new TH1D("ECalRecoPlot",Constants::TitleCalorimetricEnergy,Constants::NBinsCalorimetricEnergy,Constants::MinCalorimetricEnergy,Constants::MaxCalorimetricEnergy);
	TH1D* ECalRecoPlotMott = new TH1D("ECalRecoPlotMottQ2",Constants::TitleCalorimetricEnergy,Constants::NBinsCalorimetricEnergy,Constants::MinCalorimetricEnergy,Constants::MaxCalorimetricEnergy);
	TH1D* ECalRecoPlotMottQ2 = new TH1D("ECalRecoPlotMott",Constants::TitleCalorimetricEnergy,Constants::NBinsCalorimetricEnergy,Constants::MinCalorimetricEnergy,Constants::MaxCalorimetricEnergy);
	TH1D* EQERecoPlot = new TH1D("EQERecoPlot",Constants::TitleLeptonicEnergy,Constants::NBinsLeptonicEnergy,Constants::MinLeptonicEnergy,Constants::MaxLeptonicEnergy);

	// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Resolution Studies

//	TH1D* ECalResoPlot = new TH1D("ECalResoPlot",Constants::TitleECalReso,Constants::NBinsECalReso,Constants::MinECalReso,Constants::MaxECalReso);
//	TH1D* EQEResoPlot = new TH1D("EQEResoPlot",Constants::TitleEQEReso,Constants::NBinsEQEReso,Constants::MinEQEReso,Constants::MaxEQEReso);

	// ---------------------------------------------------------------------------------------------------------------------------------------------------------

 	TH1D* nPhotons = new TH1D("nPhotons",";n_{#gamma};",10,0,10.);
        TH1D* EISR = new TH1D("EISR",";E_{#gamma ISR};",100,0,2.);
        TH1D* EFSR = new TH1D("EFSR",";E_{#gamma FSR};",100,0,2.);
	Q2Plot->Sumw2();
	xBPlot->Sumw2();
	nuPlot->Sumw2();
	WPlot->Sumw2();
	MissMomentum->Sumw2();
	MissMomentumX->Sumw2();
	MissMomentumY->Sumw2();
	MissMomentumZ->Sumw2();
	TotalMissMomentum->Sumw2();
	ElectronEnergy->Sumw2();
	ElectronCosTheta->Sumw2();
	ElectronMomX->Sumw2();
	ElectronMomY->Sumw2();
	ElectronMomZ->Sumw2();
	ProtonEnergy->Sumw2();
	ProtonCosTheta->Sumw2();
	ProtonMomX->Sumw2();
	ProtonMomY->Sumw2();
	ProtonMomZ->Sumw2();
	ECalRecoPlot->Sumw2();
	EQERecoPlot->Sumw2();
	ElectronMomX->Sumw2();
	ElectronMomY->Sumw2();
	ElectronMomZ->Sumw2();
	ProtonMomX->Sumw2();
	ProtonMomY->Sumw2();
	ProtonMomZ->Sumw2();
	int counter = 0;
	int countRES = 0;
	int countDIS = 0;
	int countQEL = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//	for (Long64_t jentry=0; jentry<1;jentry++) {
//	for (Long64_t jentry=0; jentry<2002;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);	nbytes += nb;

		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// -------------------------------------------------------------------------------------------------------------------------------------

		int ProtonTagging = 0, ChargedPionPlusTagging = 0, ChargedPionMinusTagging = 0, GammaTagging = 0, NeutralPionTagging = 0;
		int PiPlusFinalStateIndex = 99, PiMinusFinalStateIndex = 99;
		int nphotons =0; double eISR = 0.; double eFSR = 0.;

		vector <int> ProtonID;
		ProtonID.clear();

		double mott = (pow(1./137,2.)*(cthl+1))/(2*pow(El,2.)*pow((1-cthl),2.));
	        double mottQ2 = 1./Q2;

		for (int i = 0; i < nf; i++) {

//			if (pdgf[i] == 2212 && pf[i] > 0.3) {
			if (pdgf[i] == 2212) {

				ProtonTagging ++;
				ProtonID.push_back(i);

			}

	                if (pdgf[i] == 22) {
                           nphotons++;
                           if ( (pxf[i] < 1E-10) && (pyf[i] < 1E-10) )     eISR = Ef[i];
                           if ( fabs(cthf[i] - cthl) < 1E-10) eFSR = Ef[i];
                        }
//			if (pdgf[i] == 211 && pf[i] > 0.08)  {

//				ChargedPionPlusTagging ++;

//			}

//			if (pdgf[i] == -211 && pf[i] > 0.17)  {

//				ChargedPionMinusTagging ++;

//			}

//			if (pdgf[i] == 22 || pdgf[i] == 111) { GammaTagging ++; }

		}

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		if (ProtonTagging != 1) {continue;}
		Weight = wght;

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Beam

		TLorentzVector Beam4Vector(0,0,Ebeam,Ebeam);
		double EBeam = Beam4Vector.E();

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Electron

		double pl_smeared = gRandom->Gaus(pl,ResoP*pl);
		double thetal_smeared = gRandom->Gaus(TMath::ACos(cthl),ResoTheta);
		double phil_smeared = gRandom->Gaus(TMath::ATan2(pyl,pxl),ResoPhi);
		double EePrime = sqrt( TMath::Power(pl_smeared,2.) + TMath::Power(ElectronMass,2.) );

		TVector3 FSElectron3Vector(-99.,-99.,-99.);
		FSElectron3Vector.SetMagThetaPhi(pl_smeared,thetal_smeared,phil_smeared);

		TLorentzVector FSElectron4Vector(FSElectron3Vector,EePrime);
		double FSElectronMag = FSElectron4Vector.Rho();
		double FSElectronTheta = FSElectron4Vector.Theta();
		double FSElectronTheta_Deg = FSElectronTheta * 180. / TMath::Pi();
		double FSElectronCosTheta = cos(FSElectronTheta);
		double FSElectronPhi = FSElectron4Vector.Phi();
		double FSElectronPhi_Deg = FSElectronPhi * 180. / TMath::Pi();
		if (FSElectronPhi_Deg < 0 ) { FSElectronPhi_Deg += 360.;}
		if (FSElectronPhi_Deg > 360. ) { FSElectronPhi_Deg -= 360.;}

		// Electron Transverse Components

		TVector3 FSElectron3VectorTransverse(FSElectron4Vector.X(),FSElectron4Vector.Y(),0);
		double FSElectron3VectorTransverseMag = FSElectron3VectorTransverse.Mag();

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Proton

		double pp_smeared = gRandom->Gaus(pf[ProtonID[0]],ResoP*pf[ProtonID[0]]);
		double thetap_smeared = gRandom->Gaus(TMath::ACos(cthf[ProtonID[0]]),ResoTheta);
		double phip_smeared = gRandom->Gaus(TMath::ATan2(pyf[ProtonID[0]],pxf[ProtonID[0]]),ResoPhi);
		double EpPrime = sqrt( TMath::Power(pp_smeared,2.) + TMath::Power(ProtonMass,2.) );

		TVector3 FSProton3Vector(-99.,-99.,-99.);
		FSProton3Vector.SetMagThetaPhi(pp_smeared,thetap_smeared,phip_smeared);
		TLorentzVector FSProton4Vector(FSProton3Vector,EpPrime);
		double FSProtonMag = FSProton4Vector.Rho();
		double FSProtonTheta = FSProton4Vector.Theta();
		double FSProtonTheta_Deg = FSProtonTheta * 180. / TMath::Pi();
		double FSProtonCosTheta = cos(FSProtonTheta);
		double FSProtonPhi = FSProton4Vector.Phi();
		double FSProtonPhi_Deg = FSProtonPhi * 180. / TMath::Pi();
		if (FSProtonPhi_Deg < 0 ) { FSProtonPhi_Deg += 360.;}
		if (FSProtonPhi_Deg > 360. ) { FSProtonPhi_Deg -= 360.;}

		// Proton Transverse Components

		TVector3 FSProton3VectorTransverse(FSProton4Vector.X(),FSProton4Vector.Y(),0);
		double FSProton3VectorTransverseMag = FSProton3VectorTransverse.Mag();

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Reconstucted Incoming Electron

		TLorentzVector ISRecoElectron4Vector = FSElectron4Vector + FSProton4Vector;
		TLorentzVector q4Vector = Beam4Vector - FSElectron4Vector;

		TVector3 MissingMomentum3Vector(ISRecoElectron4Vector.X(),ISRecoElectron4Vector.Y(),0);
		double MissMomValue = MissingMomentum3Vector.Mag();

		TVector3 TotalMissingMomentum3Vector = Beam4Vector.Vect() - ISRecoElectron4Vector.Vect();
		double TotalMissMomValue = TotalMissingMomentum3Vector.Mag();

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		double Q2Value = -q4Vector.Mag2();
		double NuValue = Beam4Vector.E() - FSElectron4Vector.E();
		double xBValue = Q2Value / (2*Constants::ProtonMass*NuValue);
		double WValue = sqrt( ProtonMass * ProtonMass + 2 * ProtonMass * NuValue - Q2Value);

		// Energy Reconstruction

		double EQE = (2 * ProtonMass * FSElectron4Vector.E() - pow(ElectronMass,2) ) / 2 / ( ProtonMass - FSElectron4Vector.E() + FSElectronMag * FSElectronCosTheta); 
		double ECal = FSElectron4Vector.E() + FSProton4Vector.E() - ProtonMass + BE;
	        //if (ECal > 4.32) std::cout<<"El "<<FSElectron4Vector.E()<<" Ep "<<FSProton4Vector.E()<<" Ecal "<<ECal<<" EQE "<<EQE<<std::endl;
		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		double delta_e = (FSElectronMag - ElectronCV) / ElectronCV * 100; // %
		double delta_p = (FSProtonMag - ProtonCV) / ProtonCV * 100; // %

		int IDeTrial = (delta_e - DeltaEMin)/DeltaEStep;
		int IDpTrial = (delta_p - DeltaPMin)/DeltaPStep;

double YPTarE = FSElectronTheta - ElectronThetaCV;
double YPTarP = FSProtonTheta - ProtonThetaCV;

int YPTarETrial = (YPTarE - YPTarEMin)/YPTarEStep;
int YPTarPTrial = (YPTarP - YPTarPMin)/YPTarPStep;

//int XPTarETrial = YPTarENBins / 2;
//int XPTarPTrial = YPTarPNBins / 2;

double ZVertexTrial = gRandom->Uniform(VertexZMin,VertexZMax);

int ZvertETrial = (ZVertexTrial - ZvertEMin)/ZvertEStep;
int ZvertPTrial = (ZVertexTrial - ZvertPMin)/ZvertPStep;

double AccWeight_e = h3_e->GetBinContent(IDeTrial+1,YPTarETrial+1,ZvertETrial+1);
double AccWeight_p = h3_p->GetBinContent(IDpTrial+1,YPTarPTrial+1,ZvertPTrial+1);

if (Weight > 1.  && Weight * (AccWeight_e * AccWeight_p) > 0) std::cout<<"Weight "<<Weight<<" AccWeight_e "<<AccWeight_e<<" AccWeight_p "<<AccWeight_p<<std::endl;

Weight = Weight * (AccWeight_e * AccWeight_p);
if (Weight >0 ) std::cout<<"Final Weight "<<Weight<<std::endl;
//n_e->Draw("Acc>>h_e(1,0,100)","ID == "+ToString(IDeTrial),"goff"); 
//n_e->Draw("Acc>>h_e(1,0,100)","ID == "+ToString(IDeTrial) + "  && IY == " + ToString(YPTarETrial) + 
//			  			  "  && IX == " + ToString(XPTarETrial) + "   && Zvert == " + ToString(ZvertETrial),"goff"); 
//n_p->Draw("Acc>>h_p(1,0,100)","ID == "+ToString(IDpTrial) + "  && IY == " + ToString(YPTarPTrial) + 
//						  "  && IX == " + ToString(XPTarPTrial) + "   && Zvert == " + ToString(ZvertPTrial),"goff"); 

//TH1 *myh_e = (TH1*)gDirectory->Get("h_e"); 
//double accweight_e = myh_e->GetMean();

//TH1 *myh_p = (TH1*)gDirectory->Get("h_p"); 
//double accweight_p = myh_p->GetMean();

//if (accweight_e != 0 && accweight_p != 0) { Weight = Weight * (accweight_e * accweight_p); }
//else {Weight = 0;}
//cout << "accweight_e * accweight_p = " << accweight_e * accweight_p << endl;
//cout << "Weight = " << Weight << endl;

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Selection cuts

// ADD THEM !!!
		//if (!cut_pL_pR1->IsInside(FSProtonMag,FSElectronMag)) { continue; }
		counter ++;

		if ( fabs(delta_e) > DeltaP_CV) { continue; }
		if ( fabs(delta_p) > DeltaP_CV) { continue; }

		if ( fabs(YPTarE) > YPTarAcceptance ) { continue; }
		if ( fabs(YPTarP) > YPTarAcceptance ) { continue; }

		if (Weight == 0) continue;
		//if ( fabs(FSElectronPhi) > XPTarAcceptance ) { continue; }
		//if ( fabs(FSProtonPhi) > XPTarAcceptance ) { continue; }

		// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Making sure that all the plots are filled by exactly the same events

		if (1
	/*			    (Q2Value > MinQ2 && Q2Value < MaxQ2) &&
				    (xBValue > MinxB && xBValue < MaxxB) &&
				    (WValue > MinW && WValue < MaxW) && 
				    (NuValue > Minnu && NuValue < Maxnu) && 

				    (MissingMomentum3Vector.X() > MinPmissX && MissingMomentum3Vector.X() < MaxPmissX) &&
				    (MissingMomentum3Vector.Y() > MinPmissY && MissingMomentum3Vector.Y() < MaxPmissY) &&
				    (TotalMissingMomentum3Vector.Z() > MinPmissZ && TotalMissingMomentum3Vector.Z() < MaxPmissZ) &&
				    (TotalMissMomValue > MinTotalPmiss && TotalMissMomValue < MaxTotalPmiss) &&
				    (MissMomValue > MinPmiss && MissMomValue < MaxPmiss) &&

				    (EePrime > MinElectronEnergy && EePrime < MaxElectronEnergy) &&
				    (FSElectronCosTheta > MinElectronCosTheta && FSElectronCosTheta < MaxElectronCosTheta) &&

				    (FSElectron4Vector.X() > MinElectronMomX && FSElectron4Vector.X() < MaxElectronMomX) &&
				    (FSElectron4Vector.Y() > MinElectronMomY && FSElectron4Vector.Y() < MaxElectronMomY) &&
				    (FSElectron4Vector.Z() > MinElectronMomZ && FSElectron4Vector.Z() < MaxElectronMomZ) &&

				    (EpPrime > MinEp && EpPrime < MaxEp) &&
				    (FSProtonCosTheta > MinProtonCosTheta && FSProtonCosTheta < MaxProtonCosTheta) &&

				    (FSProton4Vector.X() > MinProtonMomX && FSProton4Vector.X() < MaxProtonMomX) &&
				    (FSProton4Vector.Y() > MinProtonMomY && FSProton4Vector.Y() < MaxProtonMomY) &&
				    (FSProton4Vector.Z() > MinProtonMomZ && FSProton4Vector.Z() < MaxProtonMomZ) &&

				    (EQE > MinLeptonicEnergy && EQE < MaxLeptonicEnergy) &&
				    (ECal > MinCalorimetricEnergy && ECal < MaxCalorimetricEnergy) 
	*/	) {

			// Exclusive Analysis

			// General Plots

			Q2Plot->Fill(Q2Value,Weight);
			xBPlot->Fill(xBValue,Weight);
			nuPlot->Fill(NuValue,Weight);
			WPlot->Fill(WValue,Weight);
			MissMomentum->Fill(MissMomValue,Weight);
			TotalMissMomentum->Fill(TotalMissMomValue,Weight);
			MissMomentumX->Fill(TotalMissingMomentum3Vector.X(),Weight);
			MissMomentumY->Fill(TotalMissingMomentum3Vector.Y(),Weight);
			MissMomentumZ->Fill(TotalMissingMomentum3Vector.Z(),Weight);

	//		Q2Vsnu->Fill(NuValue,Q2Value,Weight);
	//		Q2VsW->Fill(WValue,Q2Value,Weight);

			// -----------------------------------------------------------------------------------------------------------------------------------------------------

			// Electron

			ElectronCosTheta->Fill(FSElectronCosTheta,Weight);
//			ElectronPhi->Fill(FSElectronPhi_Deg,Weight);
			ElectronEnergy->Fill(FSElectron4Vector.E(),Weight);

			ElectronMomX->Fill(FSElectron4Vector.X(),Weight);
			ElectronMomY->Fill(FSElectron4Vector.Y(),Weight);
			ElectronMomZ->Fill(FSElectron4Vector.Z(),Weight);

			// -----------------------------------------------------------------------------------------------------------------------------------------------------

			// Proton

			ProtonCosTheta->Fill(FSProtonCosTheta,Weight);	
//			ProtonPhi->Fill(FSProtonPhi_Deg,Weight);
			ProtonEnergy->Fill(FSProton4Vector.E(),Weight);

	//		EpVsMissMomentum->Fill(MissMomValue,FSProton4Vector.E(),Weight);
	//		Q2VsMissMomentum->Fill(MissMomValue,Q2Value,Weight);
	//		ElectronCosThetaVsMissMomentum->Fill(MissMomValue,FSElectronCosTheta,Weight);

			ProtonMomX->Fill(FSProton4Vector.X(),Weight);
			ProtonMomY->Fill(FSProton4Vector.Y(),Weight);
			ProtonMomZ->Fill(FSProton4Vector.Z(),Weight);

			// -----------------------------------------------------------------------------------------------------------------------------------------------------

			// Reconstructed Energy

			ECalRecoPlot->Fill(ECal,Weight);
			ECalRecoPlotMott->Fill(ECal,Weight/mott);
			ECalRecoPlotMottQ2->Fill(ECal,Weight/mottQ2);
			EQERecoPlot->Fill(EQE,Weight);

			// Resolution Studies

//			ECalResoPlot->Fill(ECalReso,Weight);
//			EQEResoPlot->Fill(EQEReso,Weight);

			nPhotons->Fill(nphotons,Weight);
                        EISR->Fill(eISR,Weight);
                        EFSR->Fill(eFSR,Weight);

			if (res) countRES++; 
			if (qel) countQEL++; 
			if (dis) countDIS++; 

		} // End of the "same events filling" 		

	} // End of the loop over the events

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Print this message after the loop over the events

	std::cout << std::endl;
	std::cout << "counter = "<< counter << std::endl; 
	std::cout << "File " << file_name << ".root created " << std::endl; 
	std::cout << std::endl;
	std::cout << "countQEL "<<countQEL<<" countRES "<<countRES<<" countDIS "<<countDIS<<std::endl;
	file->Write(); file->Close(); 

} // End of the program
