#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <TGaxis.h>

#include <iostream>
#include <vector>

using namespace std;

// ----------------------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------------------

void PrintAccMapsWeights() {

	gStyle->SetOptStat(0);	

//	TString Energy = "1_161";
//	TString Energy = "2_261";
	TString Energy = "4_461";

//	TString Nucleus = "4He";
//	TString Nucleus = "12C";
	TString Nucleus = "56Fe";

	TString PathToFiles = "../../myFiles/"+ Energy + "/SuSav2_RadCorr_LFGM/NoxBCut/";
	TString FileName = PathToFiles + Nucleus+"_" + Energy + "_SuSav2_RadCorr_LFGM_Plots_FSI_em.root";
	TFile* file_acceptance = TFile::Open(FileName);

	// -------------------------------------------------------------------------------

	TH1D* electrons = (TH1D*)file_acceptance->Get("h1_Electron_AccMapWeights");
	TCanvas* canElectrons = new TCanvas("canElectrons","canElectrons",205,34,1024,768);
	electrons->GetXaxis()->SetRangeUser(-0.05,1.05);
	electrons->Draw("hist");

	// -------------------------------------------------------------------------------

	TH1D* protons = (TH1D*)file_acceptance->Get("h1_Proton_AccMapWeights");
	TCanvas* canprotons = new TCanvas("canprotons","canprotons",205,34,1024,768);
	protons->GetXaxis()->SetRangeUser(-0.05,1.05);
	protons->Draw("hist");

	// -------------------------------------------------------------------------------

	TH1D* piplus = (TH1D*)file_acceptance->Get("h1_PiPlus_AccMapWeights");
	TCanvas* canpiplus = new TCanvas("canpiplus","canpiplus",205,34,1024,768);
	piplus->GetXaxis()->SetRangeUser(-0.05,1.05);
	piplus->Draw("hist");

	// -------------------------------------------------------------------------------

	TH1D* piminus = (TH1D*)file_acceptance->Get("h1_PiMinus_AccMapWeights");
	TCanvas* canpiminus = new TCanvas("canpiminus","canpiminus",205,34,1024,768);
	piminus->GetXaxis()->SetRangeUser(-0.05,1.05);
	piminus->Draw("hist");

	// ------------------------------------------------------------------------



} // End of the program
