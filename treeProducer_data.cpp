#define treeProducer_data_cxx
#include "treeProducer_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <iomanip>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

//#include "acceptance_c.cpp"

using namespace std;

void treeProducer_data::Loop() {

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	TH1D::SetDefaultSumw2();

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------

	TString WhichMap = "e2a_maps";
	//TFile* file_acceptance = TFile::Open("maps/"+WhichMap+"/"+WhichMap+"_12C_E_"+E+".root");
	TFile* file = new TFile(file_name+"_em.root","recreate");
	std::cout << std::endl << "File " << file_name << "_em.root will be created" << std::endl << std::endl; 

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------

	// Constants

	TLorentzVector Beam(0.,0.,Ebeam,Ebeam);
	double ProtonMass = 0.938, ElectronMass = 0.000511; // GeV
	double fine_struc_const = 0.007297;

	// -----------------------------------------------------------------------------------------------------------------------------------------------------------

	// Ranges & Binning


	// -----------------------------------------------------------------------------------------------------------------------------------------------------------

	int countEvents = 0; Long64_t nbytes = 0, nb = 0;

//	for (Long64_t jentry=0; jentry<nentries;jentry++) {
	for (Long64_t jentry=0; jentry<10000;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Initializing counters and creating vectors to store id of particles in the array

		int ProtonTagging = 0, ChargedPionPlusTagging = 0, ChargedPionMinusTagging = 0, GammaTagging = 0;
		vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
		ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Loop over te final state particles an counting of the particles above threshold

		for (int i = 0; i < nf; i++) {

			if (pdgf[i] == 2212 && pf[i] > 0.3) {

				ProtonTagging ++;
				ProtonID.push_back(i);

			}

			if ( pdgf[i] == 211 && pf[i] > 0.15)  {

				ChargedPionPlusTagging ++;
				PiPlusID.push_back(i);

			}

			if (pdgf[i] == -211 && pf[i] > 0.15)  {

				ChargedPionMinusTagging ++;
				PiMinusID.push_back(i);

			}

			if (pdgf[i] == 22 && pf[i] > 0.3) { 

				GammaTagging ++; 
				PhotonID.push_back(i);

			}

		}

		// Total Number of Charged Pions above threshold

		int ChargedPionTagging = ChargedPionPlusTagging+ChargedPionMinusTagging;

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Outgoing e'
    
		TLorentzVector ElectronOut(pxl,pyl,pzl,El);
		double pElectronOut = ElectronOut.Rho(); 
		double theta_ElectronOut = ElectronOut.Theta();
		double phi_ElectronOut = ElectronOut.Phi()+TMath::Pi();

		// q vector  
    
		TLorentzVector qVector = Beam-ElectronOut;

		double nu = Ev - El;
		double xB = Q2 / (2*ProtonMass*nu);

		// Leptonic Energy Reconstruction

//		double ERecoOnlyePrime = (2*ProtonMass*BE + 2*ProtonMass*El - pow(ElectronMass,2)) / 2 / (ProtonMass - El + pElectronOut*cos(theta_ElectronOut)); 

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Normalize by the Mott xsection for the difference in the masses of the propagators

		double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(theta_ElectronOut)+1))/(2*pow(El,2)*pow((1-cos(theta_ElectronOut)),2.));

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Define the weight for GENIE events using the acceptance maps

	        // Electron Weight

		double e_acc_ratio = 1.;
//		if (FSIModel != "Data") { e_acc_ratio = acceptance_c(pElectronOut, cos(theta_ElectronOut), phi_ElectronOut, 11,file_acceptance); }
		
		// Safety check

		if ( fabs(e_acc_ratio) != e_acc_ratio) { continue; }

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		// Apply the Selection Cuts

		if (W > 2) continue;
		if (Ebeam == 1.161) { if (Q2 < 0.1) continue; } 
		if (Ebeam == 2.261) { if (Q2 < 0.4) continue; } 
		if (Ebeam == 4.461) { if (Q2 < 1.) continue; } 
		if (xBCut == "xBCut") { if (fabs(xB-1.) > 0.2) continue; }

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		//Rotations

		// Plots

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------


		// Counter of events that are used to fill our plots

		countEvents ++;
	
	} // end of the loop over events

	// __________________________________________________________________________________________________________________________________________________________________________________________

	// Print this message after the loop over the events

	std::cout << std::endl;
	std::cout << "File " << file_name << "_em.root created " << std::endl; 
	std::cout << std::endl;
	std::cout << "Efficiency = " << double(countEvents)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;
	file->Write(); file->Close(); 

} // End of the program
