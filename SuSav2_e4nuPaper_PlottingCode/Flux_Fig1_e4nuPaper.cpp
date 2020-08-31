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
	double TextSize = 0.10;	
	int LineWidth = 8;
	
	gStyle->SetOptStat(0);	
	
	double min = 0.5;
	double max = 6;
	int bins = 50;	
	
	// ----------------------------------------------------------------------------------------------------------------
	
	TCanvas* ClasFluxCanvas = new TCanvas("ClasFluxCanvas","ClasFluxCanvas",205,34,1024,768);
	ClasFluxCanvas->SetTitle();
	ClasFluxCanvas->SetBottomMargin(0.2);
	ClasFluxCanvas->SetLeftMargin(0.15);	
	
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
	CLAS_Graph->GetXaxis()->SetTitleOffset(0.95);
	
	CLAS_Graph->GetYaxis()->SetRangeUser(0,10);
	CLAS_Graph->GetYaxis()->SetNdivisions(6);
//	CLAS_Graph->GetYaxis()->SetTitle("Extracted #Phi_{e} [a.u]");
	CLAS_Graph->GetYaxis()->SetTitle("Incident #nu_{e} Flux");
	CLAS_Graph->GetYaxis()->SetTitleFont(TextFont);
	CLAS_Graph->GetYaxis()->SetLabelFont(TextFont);
	CLAS_Graph->GetYaxis()->SetTitleSize(TextSize);
	CLAS_Graph->GetYaxis()->SetLabelSize(TextSize);	
	CLAS_Graph->GetYaxis()->CenterTitle();
	CLAS_Graph->GetYaxis()->SetTickSize(0.02);		
	CLAS_Graph->GetYaxis()->SetTitleOffset(0.7);		
	
	CLAS_Graph->SetTitle();
	CLAS_Graph->SetLineWidth(LineWidth);
	CLAS_Graph->SetLineColor(kGreen+2);	
	CLAS_Graph->Draw("ac");

	ClasFluxCanvas->SaveAs("../../myPlots/pdf/IncidentFlux.pdf");

//	TLatex* ClasLatex = new TLatex(0.65,0.8,"Incident");
//	ClasLatex->SetTextFont(TextFont);
//	ClasLatex->SetTextSize(TextSize);
//	ClasLatex->SetNDC();
//	ClasLatex->Draw();
	
	// ------------------------------------------------------------------------------------------------------------
	
	TCanvas* GENIEFluxCanvas = new TCanvas("GENIEFluxCanvas","GENIEFluxCanvas",205,34,1024,768);
	GENIEFluxCanvas->SetTitle();
	GENIEFluxCanvas->SetBottomMargin(0.2);
	GENIEFluxCanvas->SetLeftMargin(0.15);	
	
	TGraph *mix_Graph = new TGraph(nPoints,t->GetV1(),t->GetV3());
	
	mix_Graph->GetXaxis()->SetRangeUser(min,max);
	mix_Graph->GetXaxis()->SetNdivisions(6);
	mix_Graph->GetXaxis()->SetTitle("E_{#nu} [GeV]");
	mix_Graph->GetXaxis()->SetTitleFont(TextFont);
	mix_Graph->GetXaxis()->SetLabelFont(TextFont);
	mix_Graph->GetXaxis()->SetTitleSize(TextSize);
	mix_Graph->GetXaxis()->SetLabelSize(TextSize);	
	mix_Graph->GetXaxis()->CenterTitle();
	mix_Graph->GetXaxis()->SetTickSize(0.02);
	mix_Graph->GetXaxis()->SetTitleOffset(0.95);	
	
	mix_Graph->GetYaxis()->SetRangeUser(0,10);
	mix_Graph->GetYaxis()->SetNdivisions(6);
//	mix_Graph->GetYaxis()->SetTitle("Extracted #Phi_{e} [a.u]");
	mix_Graph->GetYaxis()->SetTitle("Inferred #nu_{e} Flux");
	mix_Graph->GetYaxis()->SetTitleFont(TextFont);
	mix_Graph->GetYaxis()->SetLabelFont(TextFont);
	mix_Graph->GetYaxis()->SetTitleSize(TextSize);
	mix_Graph->GetYaxis()->SetLabelSize(TextSize);	
	mix_Graph->GetYaxis()->CenterTitle();
	mix_Graph->GetYaxis()->SetTickSize(0.02);		
	mix_Graph->GetYaxis()->SetTitleOffset(0.7);		
	
	mix_Graph->SetTitle();	
	
	mix_Graph->SetLineWidth(LineWidth);
	mix_Graph->SetLineColor(kBlue);
	mix_Graph->Draw("ac");	

	TH1D* CLAS_Graph_Clone = (TH1D*)(CLAS_Graph->Clone());
	CLAS_Graph_Clone->SetLineWidth(2);
	CLAS_Graph_Clone->Draw("same c");

	GENIEFluxCanvas->SaveAs("../../myPlots/pdf/InferredFlux.pdf");

//	TLatex* GenieLatex = new TLatex(0.55,0.8,"Reconstructed");
//	GenieLatex->SetTextFont(TextFont);
//	GenieLatex->SetTextSize(TextSize);
//	GenieLatex->SetNDC();
//	GenieLatex->Draw();	
	
	// ------------------------------------------------------------------------------------------------------------
	
//	TLegend* leg = new TLegend(0.47,0.75,0.67,0.89);
//	
//	leg->SetBorderSize(0);
//	leg->SetTextFont(TextFont);
//	leg->SetTextSize(TextSize);
//	leg->AddEntry(mix_Graph,"Extracted with GENIE","l");
//	leg->AddEntry(CLAS_Graph,"Extracted with CLAS","l");	
//	leg->Draw();
	
} // End of the program
