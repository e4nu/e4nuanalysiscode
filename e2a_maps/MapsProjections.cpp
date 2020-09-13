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

void MapsProjections() {

//	TString Energy = "1_161";
//	TString Energy = "2_261";
	TString Energy = "4_461";


//	TFile* file_acceptance = TFile::Open("/home/afroditi/Downloads/e2a_solid_2261_2250_e.root");

	TFile* file_acceptance = TFile::Open("e2a_maps_12C_E_"+Energy+".root");
//	TFile* file_acceptance = TFile::Open("e2a_maps_12C_E_"+Energy+"_p.root");
//	TFile* file_acceptance = TFile::Open("e2a_maps_12C_E_"+Energy+"_pip.root");

	TH3D* reco = (TH3D*)file_acceptance->Get("Accepted Particles");
	TH3D* gen = (TH3D*)file_acceptance->Get("Generated Particles");	

//	TString Option = "xy";
	TString Option = "yz";
//	TString Option = "xz";

//	TCanvas* canReco = new TCanvas("canReco","canReco",205,34,1024,768);

	TProfile2D* recoProf2D = reco->Project3DProfile(Option);
//	gStyle->SetOptStat(0);	

//	recoProf2D->Draw("coltz");

	// ------------------------------------------------------------------------

//	TCanvas* canGen = new TCanvas("canGen","canGen",205,34,1024,768);

	TProfile2D* genProf2D = gen->Project3DProfile(Option);
//	gStyle->SetOptStat(0);	

//	genProf2D->Draw("coltz");

	// ------------------------------------------------------------------------

	TCanvas* canRatio = new TCanvas("canRatio","canRatio",205,34,1024,768);

	TProfile2D* recoProf2DClone = (TProfile2D*)(recoProf2D->Clone());
	recoProf2DClone->Divide(genProf2D);
	gStyle->SetOptStat(0);	

	recoProf2DClone->Draw("coltz");


} // End of the program
