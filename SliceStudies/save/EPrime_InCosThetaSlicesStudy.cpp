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

int TextFont = 132;
double TextSize = 0.07;

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

void PrettyPlot(TH1D* h){

		for (int i = 0; i < 2; i ++) { h->Rebin(); }
		h->Scale(1./h->Integral());
		ReweightPlots(h);
		h->SetLineWidth(3);

		h->GetXaxis()->CenterTitle();		
		h->GetXaxis()->SetRangeUser(0.5,2.5);
		h->GetXaxis()->SetTitle("E_{e'} [GeV]");
		h->GetXaxis()->SetTitleSize(TextSize);
		h->GetXaxis()->SetLabelSize(TextSize);
		h->GetXaxis()->SetTitleFont(TextFont);		
		h->GetXaxis()->SetLabelFont(TextFont);
		h->GetXaxis()->SetNdivisions(5);

		h->GetYaxis()->CenterTitle();		
		h->GetYaxis()->SetTitle("Weighted Events/GeV");
		h->GetYaxis()->SetTitleSize(TextSize);
		h->GetYaxis()->SetLabelSize(TextSize);
		h->GetYaxis()->SetTitleFont(TextFont);		
		h->GetYaxis()->SetLabelFont(TextFont);
		h->GetYaxis()->SetNdivisions(5);


}

// ----------------------------------------------------------------------------------------------------------------

void EPrime_InCosThetaSlicesStudy() {

	TFile* file_Data = TFile::Open("../../myFiles/2_261/Data_Final/NoxBCut/12C_2_261_Data_Final_Plots_FSI_em.root","readonly");
	TFile* file_SuSav2 = TFile::Open("../../myFiles/2_261/SuSav2_RadCorr_LFGM/NoxBCut/12C_2_261_SuSav2_RadCorr_LFGM_Plots_FSI_em.root","readonly");	
	
	// ----------------------------------------------------------------------------------------------------------------
	
	std::vector<int> Colors{kYellow,kGreen,kViolet,kBlack,kRed+1,kBlue,kMagenta,410,kYellow+1,kOrange+7};
	
	gStyle->SetOptStat(0);

	// ----------------------------------------------------------------------------------------------------------------	

	// SuSav2 EePrime in CosThetaPrime Slices 		

	double MinCosThetaEPrime = 0, MaxCosThetaEPrime = 1.; int CosThetaEPrimeSlices = 20;
	double CosThetaEPrimeStep = (MaxCosThetaEPrime - MinCosThetaEPrime) / CosThetaEPrimeSlices;

	//int Counter = 0;	
	
//	for (int WhichCosThetaSlice = 11 ; WhichCosThetaSlice < CosThetaEPrimeSlices; WhichCosThetaSlice++ ) {
	for (int WhichCosThetaSlice = 11 ; WhichCosThetaSlice < 19; WhichCosThetaSlice++ ) {

		TCanvas* SuSav2Canvas_EePrime_InThetaSlices = new TCanvas("SuSav2Canvas_ECal_EePrimeAndThetaSlices_CosThetaBin_"+ToString(WhichCosThetaSlice),\
								               "SuSav2Canvas_ECal_EePrimeAndThetaSlices_CosThetaBin_"+ToString(WhichCosThetaSlice),\
										205,34,1024,768);
		
		SuSav2Canvas_EePrime_InThetaSlices->SetBottomMargin(0.15);
		SuSav2Canvas_EePrime_InThetaSlices->SetLeftMargin(0.15);


		TH1D* h1_EePrime_SuSav2 = (TH1D*)file_SuSav2->Get(Form("h1_EePrime_InCosThetaE_%d_To_%d_Slices",WhichCosThetaSlice,WhichCosThetaSlice+1));

		TH1D* h1_EePrime_Data = (TH1D*)file_Data->Get(Form("h1_EePrime_InCosThetaE_%d_To_%d_Slices",WhichCosThetaSlice,WhichCosThetaSlice+1));
		

		PrettyPlot(h1_EePrime_Data);
		PrettyPlot(h1_EePrime_SuSav2);		

		h1_EePrime_Data->SetLineColor(kBlack);
		h1_EePrime_Data->GetYaxis()->SetRangeUser(0,1.25*TMath::Max(h1_EePrime_Data->GetMaximum(),h1_EePrime_SuSav2->GetMaximum()));							
		h1_EePrime_Data->Draw("hist same");
		
		h1_EePrime_SuSav2->SetLineColor(kBlue);
		h1_EePrime_SuSav2->Draw("hist same");

		TLegend* legSuSav2_EePrime_InCosThetaSlices = new TLegend(0.27,0.79,0.45,0.89);
		legSuSav2_EePrime_InCosThetaSlices->SetBorderSize(0);
		legSuSav2_EePrime_InCosThetaSlices->SetTextFont(TextFont);
		legSuSav2_EePrime_InCosThetaSlices->SetTextSize(TextSize);	
		
		legSuSav2_EePrime_InCosThetaSlices->AddEntry(h1_EePrime_SuSav2,\
						ToString(WhichCosThetaSlice*CosThetaEPrimeStep)+ " < cos(#theta_{e'}) < " + ToString( (WhichCosThetaSlice+1)*CosThetaEPrimeStep)
						,"");

		legSuSav2_EePrime_InCosThetaSlices->Draw();

	}
	

	// ----------------------------------------------------------------------------------------------------------------	


} // End of the program
