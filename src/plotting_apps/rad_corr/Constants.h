#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

namespace Constants {

	// Global Constants (Masses & Binding Energy)

	const double ProtonMass = .9383, MuonMass = .106, NeutronMass = 0.9396, ElectronMass = 0.000511; // [GeV]
	double BE = 0.0004; // GeV

	// Ranges & Binning

	int NBins = 30;
//	int NBins = 40;

	int NBinsQ2 = NBins; double MinQ2 = 1.4, MaxQ2 = 1.65; TString TitleQ2 = ";Q^{2} [GeV^{2}/c^{2};";
	int NBinsxB = NBins; double MinxB = 0.98, MaxxB = 1.02; TString TitlexB = ";x_{B};";
	int NBinsnu = NBins; double Minnu = 0.77, Maxnu = 0.86; TString Titlenu = ";Energy Transfer [GeV];";
	int NBinsW = NBins; double MinW = 0.92, MaxW = 0.96; TString TitleW = ";W [GeV/c^{2}];";
	int NBinsPmiss = NBins; double MinPmiss = 0., MaxPmiss = 0.035; TString TitlePmiss = ";P^{miss}_{#perp} [GeV/c];";

	int NBinsPmissX = NBins; double MinPmissX = -0.015, MaxPmissX = 0.015; TString TitlePmissX = ";P^{miss}_{x} [GeV/c];";
	int NBinsPmissY = NBins; double MinPmissY = -0.035, MaxPmissY = 0.035; TString TitlePmissY = ";P^{miss}_{y} [GeV/c];";
	int NBinsPmissZ = NBins; double MinPmissZ = -0.015, MaxPmissZ = 0.015; TString TitlePmissZ = ";P^{miss}_{z} [GeV/c];";

	int NBinsTotalPmiss = NBins; double MinTotalPmiss = 0, MaxTotalPmiss = 0.03; TString TitleTotalPmiss = ";P^{miss} [GeV/c];";
//	int NBinsReso = NBins; double MinReso = -50., MaxReso = 10.; TString TitleReso = ";#frac{E^{cal} - E^{beam}}{E^{beam}} (%);";

//	TString TitleQ2Vsnu = ";Energy Transfer (GeV);Q^{2} (GeV^{2}/c^{2});";
//	TString TitleQ2VsW = ";W (GeV/c^{2});Q^{2} (GeV^{2}/c^{2});";

	// Electron

	int NBinsElectronEnergy = NBins; double MinElectronEnergy = 3.46, MaxElectronEnergy = 3.56; TString TitleElectronEnergy = ";E_{e'} [GeV];";
//	int NBinsElectronPhi = 45; double MinElectronPhi = 0., MaxElectronPhi = 360.; TString TitleElectronPhi = ";#phi_{e'} (degrees);";
	int NBinsElectronCosTheta = NBins; double MinElectronCosTheta = 0.946, MaxElectronCosTheta = 0.954; TString TitleElectronCosTheta = ";cos(#theta_{e'});";

	int NBinsElectronMomX = NBins; double MinElectronMomX = 1.03, MaxElectronMomX = 1.13; TString TitleElectronMomX = ";P_{e',x} [GeV/c];";
	int NBinsElectronMomY = NBins; double MinElectronMomY = -0.06, MaxElectronMomY = 0.06; TString TitleElectronMomY = ";P_{e',y} [GeV/c];";
	int NBinsElectronMomZ = NBins; double MinElectronMomZ = 3.25, MaxElectronMomZ = 3.45; TString TitleElectronMomZ = ";P_{e',z} [GeV/c];";

//	TString TitleElectronPhiVsTheta = ";#phi_{e'} (degrees);#theta_{e'} (degrees);";
	
	// Proton

	int NBinsEp = NBins; double MinEp = 1.71, MaxEp = 1.81; TString TitleEp = ";E_{p} (GeV);";
//	int NBinsProtonPhi = 45; double MinProtonPhi = 0., MaxProtonPhi = 360.; TString TitleProtonPhi = ";#phi_{p} (degrees);";
	int NBinsProtonCosTheta = NBins; double MinProtonCosTheta = 0.654, MaxProtonCosTheta = 0.684; TString TitleProtonCosTheta = ";cos(#theta_{p});";

	int NBinsProtonMomX = NBins; double MinProtonMomX = -1.13, MaxProtonMomX = -1.03; TString TitleProtonMomX = ";P_{p,x} [GeV/c];";
	int NBinsProtonMomY = NBins; double MinProtonMomY = -0.06, MaxProtonMomY = 0.06; TString TitleProtonMomY = ";P_{p,y} [GeV/c];";
	int NBinsProtonMomZ = NBins; double MinProtonMomZ = 0.93, MaxProtonMomZ = 1.04; TString TitleProtonMomZ = ";P_{p,z} [GeV/c];";

//	TString TitleProtonEnergyVsMissMomentum = ";P^{miss}_{#perp} (GeV/c);E_{p} (GeV);";
//	TString TitleQ2VsMissMomentum = ";P^{miss}_{#perp} (GeV/c);Q^{2} (GeV^{2}/c^{2});";
//	TString TitleElectronCosThetaVsMissMomentum = ";P^{miss}_{#perp} (GeV/c);cos(#theta_{e'});";

	// Vertex

//	int NBinsVertexLocation = 20; 

//	double MinVertexXLocation = -4, MaxVertexXLocation = 0; TString TitleVertexXVsMissMomentum = ";P^{miss}_{#perp} (GeV/c);Vertex X Position (mm);";
//	double MinVertexYLocation = -2, MaxVertexYLocation = 3; TString TitleVertexYVsMissMomentum = ";P^{miss}_{#perp} (GeV/c);Vertex Y Position (mm);";
//	double MinVertexZLocation = -11, MaxVertexZLocation = 14; TString TitleVertexZVsMissMomentum = ";P^{miss}_{#perp} (GeV/c);Vertex Z Position (cm);";


	// Reconstructed Energy

	int NBinsCalorimetricEnergy = NBins; double MinCalorimetricEnergy = 4.3,MaxCalorimetricEnergy = 4.334; TString TitleCalorimetricEnergy = ";E_{cal} [GeV];";
	int NBinsLeptonicEnergy = NBins; double MinLeptonicEnergy = 4.3,MaxLeptonicEnergy = 4.345; TString TitleLeptonicEnergy = ";E_{QE} [GeV];";

	// Resolution Studies

//	int NBinsECalReso = 100; double MinECalReso = -6.,MaxECalReso = 6.; TString TitleECalReso = ";E^{Cal} Resolution (%);";
//	int NBinsEQEReso = 100; double MinEQEReso = -6.,MaxEQEReso = 6.; TString TitleEQEReso = ";E^{QE} Resolution(%);";

	// Spectrometer

	double ElectronCV = 3.54334, ElectronCVReso = 0.045; 
	double ProtonCV = 1.4805, ProtonCVReso = 0.045;

	double DeltaEMin = -4.;// %
	double DeltaEMax = 4.;// %
	double DeltaENBins = 8.;
	double DeltaEStep = fabs(DeltaEMax - DeltaEMin) / DeltaENBins;

	double DeltaPMin = -4.;// %
	double DeltaPMax = 4.;// %
	double DeltaPNBins = 8.;
	double DeltaPStep = fabs(DeltaEMax - DeltaEMin) / DeltaENBins;

	double YPTarEMin = -0.03;// %
	double YPTarEMax = 0.03;// %
	double YPTarENBins = 6.;
	double YPTarEStep = fabs(YPTarEMax - YPTarEMin) / YPTarENBins;

	double YPTarPMin = -0.03;// %
	double YPTarPMax = 0.03;// %
	double YPTarPNBins = 6.;
	double YPTarPStep = fabs(YPTarPMax - YPTarPMin) / YPTarPNBins;

	double XPTarEMin = -0.05;// %
	double XPTarEMax = 0.05;// %
	double XPTarENBins = 10.;
	double XPTarEStep = fabs(XPTarEMax - XPTarEMin) / XPTarENBins;

	double XPTarPMin = -0.05;// %
	double XPTarPMax = 0.05;// %
	double XPTarPNBins = 10.;
	double XPTarPStep = fabs(XPTarPMax - XPTarPMin) / XPTarPNBins;

	double ZvertEMin = -0.09;// m
	double ZvertEMax = 0.09;// m
	double ZvertENBins = 18.;
	double ZvertEStep = fabs(ZvertEMax - ZvertEMin) / ZvertENBins;

	double ZvertPMin = -0.09;// m
	double ZvertPMax = 0.09;// m
	double ZvertPNBins = 18.;
	double ZvertPStep = fabs(ZvertPMax - ZvertPMin) / ZvertPNBins;

	double ElectronThetaCV = 0.310698;
	double ProtonThetaCV = 0.852082;

//	double AngularAcceptance = 0.22*TMath::Pi() / 180.;

//	double reso = 0.0005;
//	double reso = 0.0001;

//	double EFrac = 1. / 2.;
//	double PFrac = sqrt(3.) / 2.;

//	double EFrac = 1.;
//	double PFrac = 1.;

//	double XReso = 3.09*TMath::Power(10,-3);
//	double YReso = 6.68*TMath::Power(10,-3);
//	double ZReso = 2.7*TMath::Power(10,-3);

//double coeff = 0.5;
//	double XReso = coeff*TMath::Power(10,-3);
//	double YReso = 1.*TMath::Power(10,-3);
//	double ZReso = 1.*TMath::Power(10,-3);

//	double XResoE = XReso * EFrac;	
//	double YResoE = YReso * EFrac;
//	double ZResoE = ZReso * EFrac;

//	double XResoP = XReso * PFrac;		
//	double YResoP = YReso * PFrac;		
//	double ZResoP = ZReso * PFrac;

//	double ResoP = 2.65 * TMath::Power(10,-4.); // Pmiss_z // 3.		
//	double ResoTheta =  1.7 * 0.5 * TMath::Power(10,-3.); // Pmiss_x // 1.7 
//	double ResoPhi = 4. * 1. * TMath::Power(10,-3.); // Pmiss_y // 4.

	//adi
	//double ResoP = 1.5 * TMath::Power(10,-4.); // Pmiss_z // 3.		
	//double ResoTheta = 0.8 * TMath::Power(10,-3.); //1.7* 0.2  * TMath::Power(10,-3.); // Pmiss_x // 1.7 // 0.5
	//double ResoPhi = 3.75 * TMath::Power(10,-3.);; //4. * 1. * TMath::Power(10,-3.); // Pmiss_y // 4.

	//adi
	double ResoP = 2. * TMath::Power(10,-4.); // Pmiss_z // 3.		
	double ResoTheta = 0.6 * TMath::Power(10,-3.); //1.7* 0.2  * TMath::Power(10,-3.); // Pmiss_x // 1.7 // 0.5
	double ResoPhi = 5.25 * TMath::Power(10,-3.);; //4. * 1. * TMath::Power(10,-3.); // Pmiss_y // 4.

	//double YPTarAcceptance = 0.03;
	double YPTarAcceptance = 0.22; // Based on a conversation with Rey or 25 26
	//double XPTarAcceptance = 0.02;  // 0.055
	double XPTarAcceptance = 0.055;  // 0.055

	double VertexXMin = -0.001; // m
	double VertexXMax = 0.003; // m

	double VertexYMin = -0.0007; // m
	double VertexYMax = 0.0017; // m

	double VertexZMin = -0.05; // m
	double VertexZMax = 0.05; // m

	//double DeltaP_CV = 2.225;  // 4%
	double DeltaP_CV = 1.75;  // 4%
	
}
#endif
