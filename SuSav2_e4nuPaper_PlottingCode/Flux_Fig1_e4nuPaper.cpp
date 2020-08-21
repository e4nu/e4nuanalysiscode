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

#include  "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"

// ----------------------------------------------------------------------------------------------------------------

void Flux_Fig1_e4nuPaper() {

	TFile* file_Data = TFile::Open("../../myFiles/extracted_flux_pts.csv.root","readonly");
	TTree* t = (TTree*)(file_Data->Get("t"));
	
	// ----------------------------------------------------------------------------------------------------------------
	
	int TextFont = 132;
	double TextSize = 0.06;	
	
	gStyle->SetOptStat(0);	
	
	double min = 0.5;
	double max = 6;
	int bins = 50;	
	
	// ----------------------------------------------------------------------------------------------------------------
	
	TCanvas* FluxCanvas = new TCanvas("FluxCanvas","FluxCanvas",205,34,1024,768);
	FluxCanvas->SetTitle();
	FluxCanvas->SetBottomMargin(0.13);
	FluxCanvas->SetLeftMargin(0.11);	
	
	int nPoints = t->Draw("Enu_GeV:CLAS_fit:mix_fit","","goff"); 
	TGraph *CLAS_Graph = new TGraph(nPoints,t->GetV1(),t->GetV2());

	CLAS_Graph->GetXaxis()->SetRangeUser(min,max);
	CLAS_Graph->GetXaxis()->SetNdivisions(6);
	CLAS_Graph->GetXaxis()->SetTitle("E_{#nu} [GeV]");
	CLAS_Graph->GetXaxis()->SetTitleFont(TextFont);
	CLAS_Graph->GetXaxis()->SetLabelFont(TextFont);
	CLAS_Graph->GetXaxis()->SetTitleSize(TextSize);
	CLAS_Graph->GetXaxis()->SetLabelSize(TextSize);	
	CLAS_Graph->GetXaxis()->CenterTitle();
	CLAS_Graph->GetXaxis()->SetTickSize(0.02);	
	
	CLAS_Graph->GetYaxis()->SetRangeUser(0,10);
	CLAS_Graph->GetYaxis()->SetNdivisions(6);
	CLAS_Graph->GetYaxis()->SetTitle("Extracted #Phi_{e} [a.u]");
	CLAS_Graph->GetYaxis()->SetTitleFont(TextFont);
	CLAS_Graph->GetYaxis()->SetLabelFont(TextFont);
	CLAS_Graph->GetYaxis()->SetTitleSize(TextSize);
	CLAS_Graph->GetYaxis()->SetLabelSize(TextSize);	
	CLAS_Graph->GetYaxis()->CenterTitle();
	CLAS_Graph->GetYaxis()->SetTickSize(0.02);		
	CLAS_Graph->GetYaxis()->SetTitleOffset(0.7);		
	
	CLAS_Graph->SetTitle();
	CLAS_Graph->SetLineWidth(3);
	CLAS_Graph->SetLineColor(kGreen+2);	
	CLAS_Graph->Draw("ac");
	
	// ------------------------------------------------------------------------------------------------------------
	
	TGraph *mix_Graph = new TGraph(nPoints,t->GetV1(),t->GetV3());
	
	mix_Graph->SetLineWidth(3);
	mix_Graph->SetLineColor(kBlue);
	mix_Graph->SetLineStyle(5);		
	mix_Graph->Draw("c");		
	
	// ------------------------------------------------------------------------------------------------------------
	
	TLegend* leg = new TLegend(0.47,0.75,0.67,0.89);
	
	leg->SetBorderSize(0);
	leg->SetTextFont(TextFont);
	leg->SetTextSize(TextSize);
	leg->AddEntry(mix_Graph,"Extracted with GENIE","l");
	leg->AddEntry(CLAS_Graph,"Extracted with CLAS","l");	
	leg->Draw();
	
} // End of the program
