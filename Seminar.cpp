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

// Accounting for the fact that the bin width is not constant

void ReweightPlots(TH1D* h) {

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 
					
		double content = h->GetBinContent(i);
		double error = h->GetBinError(i);
		double width = h->GetBinWidth(i);
		double newcontent = content / width;
		double newerror = error / width;				
		h->SetBinContent(i,newcontent);
		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

void ApplySystUnc(TH1D* h, double systunc) {

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 
					
		double error = h->GetBinError(i);
		double newerror = error * (1. + systunc);
		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

void Seminar() {

	// ------------------------------------------------------------------------

	SetOffsetAndSize();
	TGaxis::SetMaxDigits(4);

	int Ndivisions = 6;
	int LineWidth = 3;
	int FontStyle = 132;
	double TextSize = 0.08;

	// From Mariana's analysis note

	double SystUnc1GeV = 0.02; // 2% syst uncertainty at 1.161 GeV
	double SystUnc2GeV = 0.021; // 2.1% syst uncertainty at 2.261 GeV
	double SystUnc4GeV = 0.047; // 4.7% syst uncertainty at 4.461 GeV

	// ------------------------------------------------------------------------

//	TString nucleus = "4He"; TString LabelsOfSamples = "^{4}He"; TString JustNucleus = "He";
	TString nucleus = "12C"; TString LabelsOfSamples = "^{12}C"; TString JustNucleus = "C";
//	TString nucleus = "56Fe"; TString LabelsOfSamples = "^{56}Fe"; TString JustNucleus = "Fe";

	TString E = "1_161"; double DoubleE = 1.159;
//	TString E = "2_261"; double DoubleE = 2.257;
//	TString E = "4_461"; double DoubleE = 4.453;

	TString xBCut = "NoxBCut";

	TString FSIModel = "Data_Final"; TString FSILabel = "Data"; TString DirNames = "Data";

	// ------------------------------------------------------------------------

	TCanvas* PlotCanvas = new TCanvas(nucleus+"_"+E+"_"+xBCut,nucleus+"_"+E+"_"+xBCut,205,34,1024,768);
	PlotCanvas->SetTopMargin(0.11);
	PlotCanvas->SetBottomMargin(0.19);
	PlotCanvas->SetLeftMargin(0.15);

	// ------------------------------------------------------------------------

	TLegend* leg = new TLegend(0.17,0.48,0.37,0.87);
	leg->SetBorderSize(0);
	leg->SetTextFont(FontStyle);
	leg->SetTextSize(TextSize);

	// ------------------------------------------------------------------------

	TString ECalName = "epRecoEnergy_slice_0"; TString ECalLabel = "(e,e'p)_{1p0#pi} E^{cal}";
	//TString InclusiveEQEName = "h_Erec_subtruct_piplpimi_noprot_3pi"; TString InclusiveEQELabel = "(e,e')_{0#pi} E^{QE}";
	TString EQEName = "eRecoEnergy_slice_0"; TString EQELabel = "(e,e'p)_{1p0#pi} E^{QE}";

	// ------------------------------------------------------------------------

	TString PathToFiles = "../myFiles/"+ E + "/"+FSIModel+"/"+xBCut+"/";
	TString FileName = PathToFiles+nucleus+"_"+E+"_"+FSIModel+"_Plots_FSI_em.root";
	TFile* FileSample = TFile::Open(FileName);

	// ------------------------------------------------------------------------

	TH1D* ECalPlot = (TH1D*)( FileSample->Get(ECalName));
	ReweightPlots(ECalPlot);

	// ------------------------------------------------------------------------

	// X-axis label

	ECalPlot->GetXaxis()->SetLabelFont(FontStyle);
	ECalPlot->GetXaxis()->SetTitleFont(FontStyle);
	ECalPlot->GetXaxis()->SetLabelSize(TextSize);
	ECalPlot->GetXaxis()->SetTitleSize(TextSize);
	ECalPlot->GetXaxis()->SetTitleOffset(1.05);
	ECalPlot->GetXaxis()->SetNdivisions(Ndivisions);
	ECalPlot->GetXaxis()->SetTitle("E_{reconstructed} [GeV]");
	ECalPlot->GetXaxis()->SetLabelOffset(0.014);
	if (E == "4_461") { ECalPlot->GetXaxis()->SetRangeUser(1.,10.); }		

	TLegendEntry* lECal = leg->AddEntry(ECalPlot,ECalLabel, "");
	lECal->SetTextColor(kBlue);

	TString TitleLabel = ToString(DoubleE)+ " GeV " + LabelsOfSamples;
	ECalPlot->SetTitle(TitleLabel);
	gStyle->SetTitleTextColor(kBlack);
//	gStyle->SetTitleTextColor(kOrange+7);	
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetTitleSize(TextSize,"t");

	// ------------------------------------------------------------------------

	// Y-axis label

	ECalPlot->GetYaxis()->SetLabelFont(FontStyle);
	ECalPlot->GetYaxis()->SetTitleFont(FontStyle);
	ECalPlot->GetYaxis()->SetLabelSize(TextSize);
	ECalPlot->GetYaxis()->SetTitleSize(TextSize);
	ECalPlot->GetYaxis()->SetTitleOffset(0.95);
	ECalPlot->GetYaxis()->SetNdivisions(Ndivisions);
	ECalPlot->GetYaxis()->SetTitle("Weighted Events/GeV");
	ECalPlot->GetYaxis()->SetLabelOffset(0.014);	
	ECalPlot->GetYaxis()->SetRangeUser(0.,1.05*ECalPlot->GetMaximum());		

	// ------------------------------------------------------------------------

	// Calorimetric Reconstruction

	ECalPlot->SetLineColor(kBlue);
	ECalPlot->SetLineWidth(3);
	CenterAxisTitle(ECalPlot);
	ECalPlot->Draw("hist same");

	// ------------------------------------------------------------------------

	// Inclusive QE Reconstruction
/*
	TH1D* InclusiveEQEPlot = (TH1D*)( FileSample->Get(InclusiveEQEName));
	ReweightPlots(InclusiveEQEPlot);

	InclusiveEQEPlot->SetLineColor(kRed);
	InclusiveEQEPlot->SetLineWidth(3);
	CenterAxisTitle(InclusiveEQEPlot);
	InclusiveEQEPlot->Draw("hist same");

	TLegendEntry* lInclusiveQE = leg->AddEntry(InclusiveEQEPlot,InclusiveEQELabel, "");
	lInclusiveQE->SetTextColor(kRed);
*/
	// ------------------------------------------------------------------------

	// Exclusive QE Reconstruction

	TH1D* EQEPlot = (TH1D*)( FileSample->Get(EQEName));
	ReweightPlots(EQEPlot);

	EQEPlot->SetLineColor(kGreen-2);
	EQEPlot->SetLineWidth(3);
	CenterAxisTitle(EQEPlot);
	EQEPlot->Draw("hist same");

	TLegendEntry* lQE = leg->AddEntry(EQEPlot,EQELabel, "");
	lQE->SetTextColor(kGreen-2);

	// ------------------------------------------------------------------------

	//leg->Draw();

	// ------------------------------------------------------------------------
	
	TLatex latexIntChannel;
	latexIntChannel.SetTextFont(FontStyle);
	latexIntChannel.SetTextSize(TextSize);
	latexIntChannel.SetTextColor(kBlack);	
	latexIntChannel.DrawLatexNDC(0.2,0.8,"(e,e'p)_{1p0#pi}");	
	
	// ------------------------------------------------------------------------
	
	TLatex latexECal;
	latexECal.SetTextFont(FontStyle);
	latexECal.SetTextSize(TextSize);
	latexECal.SetTextColor(kBlue);
	double XCoordCal = 0.6;
	double YCoordCal = 0.7;
	if (E == "1_161" && nucleus == "12C") { XCoordCal = 0.55; }						
	latexECal.DrawLatexNDC(XCoordCal,YCoordCal,"E_{Cal}");
	
	// ------------------------------------------------------------------------
	
	TLatex latexEQE;
	latexEQE.SetTextFont(FontStyle);
	latexEQE.SetTextSize(TextSize);
	latexEQE.SetTextColor(kGreen-2);
	double XCoordQE = 0.55;
	double YCoordQE = 0.35;	
	if (E == "1_161" && nucleus == "12C") { XCoordQE = 0.45; }
	if (E == "2_261" && nucleus == "56Fe") { YCoordQE = 0.4; }
	if (E == "4_461" && nucleus == "56Fe") { YCoordQE = 0.47; }	
	latexEQE.DrawLatexNDC(XCoordQE,YCoordQE,"E_{QE}");		

	// ------------------------------------------------------------------------

	TLine* line = new TLine(DoubleE,0.,DoubleE,ECalPlot->GetMaximum());
	line->SetLineColor(kBlack);
	line->SetLineWidth(3);
	line->Draw();

	// ------------------------------------------------------------------------

	PlotCanvas->SaveAs("SeminarPlots/Ereco_"+nucleus+"_"+E+".pdf");
}

// ----------------------------------------------------------------------------------------------------------------
