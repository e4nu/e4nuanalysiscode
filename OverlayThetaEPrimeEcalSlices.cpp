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

		h->Scale(1./h->Integral());
		ReweightPlots(h);
		h->SetLineWidth(3);

		h->GetXaxis()->CenterTitle();		
		h->GetXaxis()->SetRangeUser(0.5,2.5);
		h->GetXaxis()->SetTitle("E^{Cal} [GeV]");
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

void OverlayThetaEPrimeEcalSlices() {

	TFile* file_Data = TFile::Open("../myFiles/2_261/Data_Final/NoxBCut/12C_2_261_Data_Final_Plots_FSI_em.root","readonly");
	
	TFile* file_SuSav2 = TFile::Open("../myFiles/2_261/SuSav2_RadCorr_LFGM/NoxBCut/12C_2_261_SuSav2_RadCorr_LFGM_Plots_FSI_em.root","readonly");	
	
	// ----------------------------------------------------------------------------------------------------------------
	
//	int TextFont = 132;
//	double TextSize = 0.07;	
	
	std::vector<int> Colors{kYellow,kGreen,kViolet,kBlack,kRed+1,kBlue,kMagenta,410,kYellow+1,kOrange+7};
	
	gStyle->SetOptStat(0);	
	
	// ----------------------------------------------------------------------------------------------------------------
	
	// Data ECal in Theta Slices 
	
	TCanvas* DataCanvas_ECal_ThetaSlices = new TCanvas("DataCanvas_ECal_ThetaSlices","DataCanvas_ECal_ThetaSlices",205,34,1024,768);
	
	DataCanvas_ECal_ThetaSlices->SetBottomMargin(0.15);
	DataCanvas_ECal_ThetaSlices->SetLeftMargin(0.15);	
	
	double MinTheta = 10, MaxTheta = 60; int ThetaSlices = 10;
	int ThetaStep = (MaxTheta - MinTheta) / ThetaSlices;	
	
	TLegend* legData_ThetaSlices = new TLegend(0.2,0.35,0.4,0.85);
	legData_ThetaSlices->SetBorderSize(0);
	legData_ThetaSlices->SetTextFont(TextFont);
	legData_ThetaSlices->SetTextSize(TextSize);	
	
	
	
//	for (int WhichThetaSlice = 0 ; WhichThetaSlice < ThetaSlices; WhichThetaSlice++ ) {
	for (int WhichThetaSlice = 5 ; WhichThetaSlice < 6; WhichThetaSlice++ ) {


		TH1D* h1_ECal_InThetaSlices = (TH1D*)file_Data->Get(Form("h1_ECal_InTheta_%d_To_%d_Slices",int(MinTheta+WhichThetaSlice*ThetaStep),int(MinTheta+(WhichThetaSlice+1)*ThetaStep)));
		
		
		PrettyPlot(h1_ECal_InThetaSlices);
		h1_ECal_InThetaSlices->SetLineColor(Colors[WhichThetaSlice]);
		h1_ECal_InThetaSlices->GetYaxis()->SetRangeUser(0,10.);											
		h1_ECal_InThetaSlices->Draw("same");
		
		legData_ThetaSlices->AddEntry(h1_ECal_InThetaSlices,ToString(MinTheta+WhichThetaSlice*ThetaStep)+ "^{o} < #theta_{e'} < " + ToString(MinTheta+(WhichThetaSlice+1)*ThetaStep)+"^{o}","l");	
		
	}
	
	legData_ThetaSlices->Draw();	
	
	// ----------------------------------------------------------------------------------------------------------------
	
	// SuSav2 ECal in Theta Slices 
	
	TCanvas* SuSav2Canvas_ECal_ThetaSlices = new TCanvas("SuSav2Canvas_ECal_ThetaSlices","SuSav2Canvas_ECal_ThetaSlices",205,34,1024,768);
	
	SuSav2Canvas_ECal_ThetaSlices->SetBottomMargin(0.15);
	SuSav2Canvas_ECal_ThetaSlices->SetLeftMargin(0.15);	
	
	TLegend* legSuSav2_ThetaSlices = new TLegend(0.2,0.35,0.4,0.85);
	legSuSav2_ThetaSlices->SetBorderSize(0);
	legSuSav2_ThetaSlices->SetTextFont(TextFont);
	legSuSav2_ThetaSlices->SetTextSize(TextSize);	
	
	
	
//	for (int WhichThetaSlice = 0 ; WhichThetaSlice < ThetaSlices; WhichThetaSlice++ ) {
	for (int WhichThetaSlice = 3 ; WhichThetaSlice < 10; WhichThetaSlice++ ) {


		TH1D* h1_ECal_InThetaSlices = (TH1D*)file_SuSav2->Get(Form("h1_ECal_InTheta_%d_To_%d_Slices",int(MinTheta+WhichThetaSlice*ThetaStep),int(MinTheta+(WhichThetaSlice+1)*ThetaStep)));
		
		
		PrettyPlot(h1_ECal_InThetaSlices);
		h1_ECal_InThetaSlices->SetLineColor(Colors[WhichThetaSlice]);
		h1_ECal_InThetaSlices->GetYaxis()->SetRangeUser(0,10.);											
		h1_ECal_InThetaSlices->Draw("same");
		
		legSuSav2_ThetaSlices->AddEntry(h1_ECal_InThetaSlices,ToString(MinTheta+WhichThetaSlice*ThetaStep)+ "^{o} < #theta_{e'} < " + ToString(MinTheta+(WhichThetaSlice+1)*ThetaStep)+"^{o}","l");	
		
	}
	
	legSuSav2_ThetaSlices->Draw();
	
	// ----------------------------------------------------------------------------------------------------------------
	
	// Data ECal in EePrime Slices 
	
	TCanvas* DataCanvas_ECal_EePrimeSlices = new TCanvas("DataCanvas_ECal_EePrimeSlices","DataCanvas_ECal_EePrimeSlices",205,34,1024,768);
	
	DataCanvas_ECal_EePrimeSlices->SetBottomMargin(0.15);
	DataCanvas_ECal_EePrimeSlices->SetLeftMargin(0.15);	
	
	double MinEePrime = 0.5, MaxEePrime = 2.5; int EePrimeSlices = 10;
	double EePrimeStep = (MaxEePrime - MinEePrime) / EePrimeSlices;	
	
	TLegend* legData_EePrimeSlices = new TLegend(0.4,0.35,0.6,0.85);
	legData_EePrimeSlices->SetBorderSize(0);
	legData_EePrimeSlices->SetTextFont(TextFont);
	legData_EePrimeSlices->SetTextSize(TextSize);		
	
//	for (int WhichEePrimeSlice = 0 ; WhichEePrimeSlice < EePrimeSlices; WhichEePrimeSlice++ ) {
	for (int WhichEePrimeSlice = 2 ; WhichEePrimeSlice < 10; WhichEePrimeSlice++ ) {


		TH1D* h1_ECal_InEePrimeSlices = (TH1D*)file_Data->Get(Form("h1_ECal_InEePrime_%d_To_%d_Slices",int((MinEePrime+WhichEePrimeSlice*EePrimeStep)*1000),int((MinEePrime+(WhichEePrimeSlice+1)*EePrimeStep)*1000)));
		
		
		PrettyPlot(h1_ECal_InEePrimeSlices);
		h1_ECal_InEePrimeSlices->SetLineColor(Colors[WhichEePrimeSlice]);
		h1_ECal_InEePrimeSlices->GetYaxis()->SetRangeUser(0,10.);									
		h1_ECal_InEePrimeSlices->Draw("same");
		
		legData_EePrimeSlices->AddEntry(h1_ECal_InEePrimeSlices,ToString(MinEePrime+WhichEePrimeSlice*EePrimeStep)+ " < E_{e'} < " + ToString(MinEePrime+(WhichEePrimeSlice+1)*EePrimeStep),"l");	
		
	}
	
	legData_EePrimeSlices->Draw();
	
	// ----------------------------------------------------------------------------------------------------------------
	
	// SuSav2 ECal in EePrime Slices 
	
	TCanvas* SuSav2Canvas_ECal_EePrimeSlices = new TCanvas("SuSav2Canvas_ECal_EePrimeSlices","SuSav2Canvas_ECal_EePrimeSlices",205,34,1024,768);
	
	SuSav2Canvas_ECal_EePrimeSlices->SetBottomMargin(0.15);
	SuSav2Canvas_ECal_EePrimeSlices->SetLeftMargin(0.15);		
	
	TLegend* legSuSav2_EePrimeSlices = new TLegend(0.4,0.35,0.6,0.85);
	legSuSav2_EePrimeSlices->SetBorderSize(0);
	legSuSav2_EePrimeSlices->SetTextFont(TextFont);
	legSuSav2_EePrimeSlices->SetTextSize(TextSize);		
	
//	for (int WhichEePrimeSlice = 0 ; WhichEePrimeSlice < EePrimeSlices; WhichEePrimeSlice++ ) {
	for (int WhichEePrimeSlice = 2 ; WhichEePrimeSlice < 10; WhichEePrimeSlice++ ) {


		TH1D* h1_ECal_InEePrimeSlices = (TH1D*)file_SuSav2->Get(Form("h1_ECal_InEePrime_%d_To_%d_Slices",int((MinEePrime+WhichEePrimeSlice*EePrimeStep)*1000),int((MinEePrime+(WhichEePrimeSlice+1)*EePrimeStep)*1000)));
		
		
		PrettyPlot(h1_ECal_InEePrimeSlices);
		h1_ECal_InEePrimeSlices->SetLineColor(Colors[WhichEePrimeSlice]);
		h1_ECal_InEePrimeSlices->GetYaxis()->SetRangeUser(0,10.);									
		h1_ECal_InEePrimeSlices->Draw("same");
		
		legSuSav2_EePrimeSlices->AddEntry(h1_ECal_InEePrimeSlices,ToString(MinEePrime+WhichEePrimeSlice*EePrimeStep)+ " < E_{e'} < " + ToString(MinEePrime+(WhichEePrimeSlice+1)*EePrimeStep),"l");	
		
	}
	
	legSuSav2_EePrimeSlices->Draw();					
	
	// ----------------------------------------------------------------------------------------------------------------	

} // End of the program
