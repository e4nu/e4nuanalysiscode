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
	
	std::vector<int> Colors{kYellow,kGreen,kViolet,kBlack,kRed+1,kBlue,kMagenta,410,kYellow+1,kOrange+7};
	
	gStyle->SetOptStat(0);	
	
	double GlobalMax = 15.;
	
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
	for (int WhichThetaSlice = 3 ; WhichThetaSlice < 10; WhichThetaSlice++ ) {


		TH1D* h1_ECal_InThetaSlices = (TH1D*)file_Data->Get(Form("h1_ECal_InTheta_%d_To_%d_Slices",int(MinTheta+WhichThetaSlice*ThetaStep),int(MinTheta+(WhichThetaSlice+1)*ThetaStep)));
		
		
		PrettyPlot(h1_ECal_InThetaSlices);
		h1_ECal_InThetaSlices->SetLineColor(Colors[WhichThetaSlice]);
		h1_ECal_InThetaSlices->GetYaxis()->SetRangeUser(0,GlobalMax);										
		h1_ECal_InThetaSlices->Draw("hist same");
		
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
		h1_ECal_InThetaSlices->GetYaxis()->SetRangeUser(0,GlobalMax);										
		h1_ECal_InThetaSlices->Draw("hist same");
		
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
		h1_ECal_InEePrimeSlices->GetYaxis()->SetRangeUser(0,GlobalMax);									
		h1_ECal_InEePrimeSlices->Draw("hist same");
		
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
		h1_ECal_InEePrimeSlices->GetYaxis()->SetRangeUser(0,GlobalMax);									
		h1_ECal_InEePrimeSlices->Draw("hist same");
		
		legSuSav2_EePrimeSlices->AddEntry(h1_ECal_InEePrimeSlices,ToString(MinEePrime+WhichEePrimeSlice*EePrimeStep)+ " < E_{e'} < " + ToString(MinEePrime+(WhichEePrimeSlice+1)*EePrimeStep),"l");	
		
	}
	
	legSuSav2_EePrimeSlices->Draw();
	
	// ----------------------------------------------------------------------------------------------------------------
	
	// Data ECal in Q2 Slices 
	
	TCanvas* DataCanvas_ECal_Q2Slices = new TCanvas("DataCanvas_ECal_Q2Slices","DataCanvas_ECal_Q2Slices",205,34,1024,768);
	
	DataCanvas_ECal_Q2Slices->SetBottomMargin(0.15);
	DataCanvas_ECal_Q2Slices->SetLeftMargin(0.15);	
	
	double MinQ2 = 0., MaxQ2 = 2.; int Q2Slices = 10;
	double Q2Step = (MaxQ2 - MinQ2) / Q2Slices;	
	
	TLegend* legData_Q2Slices = new TLegend(0.4,0.3,0.6,0.85);
	legData_Q2Slices->SetBorderSize(0);
	legData_Q2Slices->SetTextFont(TextFont);
	legData_Q2Slices->SetTextSize(TextSize);		
	
//	for (int WhichQ2Slice = 0 ; WhichQ2Slice < Q2Slices; WhichQ2Slice++ ) {
	for (int WhichQ2Slice = 2 ; WhichQ2Slice < 10; WhichQ2Slice++ ) {


		TH1D* h1_ECal_InQ2Slices = (TH1D*)file_Data->Get(Form("h1_ECal_InQ2_%d_To_%d_Slices",int((MinQ2+WhichQ2Slice*Q2Step)*1000),int((MinQ2+(WhichQ2Slice+1)*Q2Step)*1000)));
		
		
		PrettyPlot(h1_ECal_InQ2Slices);
		h1_ECal_InQ2Slices->SetLineColor(Colors[WhichQ2Slice]);
		h1_ECal_InQ2Slices->GetYaxis()->SetRangeUser(0,10);									
		h1_ECal_InQ2Slices->Draw("hist same");
		
		legData_Q2Slices->AddEntry(h1_ECal_InQ2Slices,ToString(MinQ2+WhichQ2Slice*Q2Step)+ " < Q^{2} < " + ToString(MinQ2+(WhichQ2Slice+1)*Q2Step),"l");	
		
	}
	
	legData_Q2Slices->Draw();	
	
	// ----------------------------------------------------------------------------------------------------------------
	
	// SuSav2 ECal in Q2 Slices 
	
	TCanvas* SuSav2Canvas_ECal_Q2Slices = new TCanvas("SuSav2Canvas_ECal_Q2Slices","SuSav2Canvas_ECal_Q2Slices",205,34,1024,768);
	
	SuSav2Canvas_ECal_Q2Slices->SetBottomMargin(0.15);
	SuSav2Canvas_ECal_Q2Slices->SetLeftMargin(0.15);	
	
	
	TLegend* legSuSav2_Q2Slices = new TLegend(0.4,0.3,0.6,0.85);
	legSuSav2_Q2Slices->SetBorderSize(0);
	legSuSav2_Q2Slices->SetTextFont(TextFont);
	legSuSav2_Q2Slices->SetTextSize(TextSize);		
	
//	for (int WhichQ2Slice = 0 ; WhichQ2Slice < Q2Slices; WhichQ2Slice++ ) {
	for (int WhichQ2Slice = 2 ; WhichQ2Slice < 10; WhichQ2Slice++ ) {


		TH1D* h1_ECal_InQ2Slices = (TH1D*)file_SuSav2->Get(Form("h1_ECal_InQ2_%d_To_%d_Slices",int((MinQ2+WhichQ2Slice*Q2Step)*1000),int((MinQ2+(WhichQ2Slice+1)*Q2Step)*1000)));
		
		
		PrettyPlot(h1_ECal_InQ2Slices);
		h1_ECal_InQ2Slices->SetLineColor(Colors[WhichQ2Slice]);
		h1_ECal_InQ2Slices->GetYaxis()->SetRangeUser(0,10);									
		h1_ECal_InQ2Slices->Draw("hist same");
		
		legSuSav2_Q2Slices->AddEntry(h1_ECal_InQ2Slices,ToString(MinQ2+WhichQ2Slice*Q2Step)+ " < Q^{2} < " + ToString(MinQ2+(WhichQ2Slice+1)*Q2Step),"l");	
		
	}
	
	legSuSav2_Q2Slices->Draw();						
	
	// ----------------------------------------------------------------------------------------------------------------	

	// Data ECal in Q2 Slices 
	
	TCanvas* DataCanvas_ECal_EePrimeAndThetaSlices = new TCanvas("DataCanvas_ECal_EePrimeAndThetaSlices","DataCanvas_ECal_EePrimeAndThetaSlices",205,34,1024,768);
	
	DataCanvas_ECal_EePrimeAndThetaSlices->SetBottomMargin(0.15);
	DataCanvas_ECal_EePrimeAndThetaSlices->SetLeftMargin(0.15);	
	
	double MinEePrime2D = 0.5, MaxEePrime2D = 2.5; int EePrimeSlices2D = 3;
	double MinTheta2D = 10, MaxTheta2D = 60; int ThetaSlices2D = 3;	
	double EePrimeStep2D = (MaxEePrime2D - MinEePrime2D) / EePrimeSlices2D;
	int ThetaStep2D = (MaxTheta2D - MinTheta2D) / ThetaSlices2D;	
	
	TLegend* legData_EePrimeAndThetaSlices = new TLegend(0.17,0.49,0.35,0.89);
	legData_EePrimeAndThetaSlices->SetBorderSize(0);
	legData_EePrimeAndThetaSlices->SetTextFont(TextFont);
	legData_EePrimeAndThetaSlices->SetTextSize(TextSize);	

	int Counter = 0;	
	
	for (int WhichEePrimeSlice2D = 0 ; WhichEePrimeSlice2D < EePrimeSlices2D; WhichEePrimeSlice2D++ ) {

		int MinEePrimeSlice2D = (MinEePrime2D+WhichEePrimeSlice2D*EePrimeStep2D) * 1000;
		int MaxEePrimeSlice2D = (MinEePrime2D+(WhichEePrimeSlice2D+1)*EePrimeStep2D)*1000;

		for (int WhichThetaSlice2D = 1 ; WhichThetaSlice2D < ThetaSlices2D; WhichThetaSlice2D++ ) {	
		
			int MinThetaSlice2D = MinTheta2D+WhichThetaSlice2D*ThetaStep2D;
			int MaxThetaSlice2D = MinTheta2D+(WhichThetaSlice2D+1)*ThetaStep2D;


			TH1D* h1_ECal_InEePrimeAndThetaSlices = (TH1D*)file_Data->Get(Form("h1_ECal_InEePrime_%d_To_%d_InTheta_%d_To_%d_Slices",\
								MinEePrimeSlice2D,MaxEePrimeSlice2D,MinThetaSlice2D,MaxThetaSlice2D));
			
			
			PrettyPlot(h1_ECal_InEePrimeAndThetaSlices);
			h1_ECal_InEePrimeAndThetaSlices->SetLineColor(Colors[Counter+1]);
			h1_ECal_InEePrimeAndThetaSlices->GetYaxis()->SetRangeUser(0,18);									
			h1_ECal_InEePrimeAndThetaSlices->Draw("hist same");
			
			legData_EePrimeAndThetaSlices->AddEntry(h1_ECal_InEePrimeAndThetaSlices,\
								ToString(MinEePrimeSlice2D)+ " < E_{e'} < " + ToString(MaxEePrimeSlice2D) + " & " +\
								ToString(MinThetaSlice2D)+ " < #theta_{e'} < " + ToString(MaxThetaSlice2D)\
								,"l");	

			Counter++;
		
		}

	}
	
	legData_EePrimeAndThetaSlices->Draw();

	// ----------------------------------------------------------------------------------------------------------------	

	// SuSav2 ECal in Q2 Slices 
	
	TCanvas* SuSav2Canvas_ECal_EePrimeAndThetaSlices = new TCanvas("SuSav2Canvas_ECal_EePrimeAndThetaSlices","SuSav2Canvas_ECal_EePrimeAndThetaSlices",205,34,1024,768);
	
	SuSav2Canvas_ECal_EePrimeAndThetaSlices->SetBottomMargin(0.15);
	SuSav2Canvas_ECal_EePrimeAndThetaSlices->SetLeftMargin(0.15);		
	
	TLegend* legSuSav2_EePrimeAndThetaSlices = new TLegend(0.17,0.49,0.35,0.89);
	legSuSav2_EePrimeAndThetaSlices->SetBorderSize(0);
	legSuSav2_EePrimeAndThetaSlices->SetTextFont(TextFont);
	legSuSav2_EePrimeAndThetaSlices->SetTextSize(TextSize);	

	Counter = 0;	
	
	for (int WhichEePrimeSlice2D = 0 ; WhichEePrimeSlice2D < EePrimeSlices2D; WhichEePrimeSlice2D++ ) {

		int MinEePrimeSlice2D = (MinEePrime2D+WhichEePrimeSlice2D*EePrimeStep2D) * 1000;
		int MaxEePrimeSlice2D = (MinEePrime2D+(WhichEePrimeSlice2D+1)*EePrimeStep2D)*1000;

		for (int WhichThetaSlice2D = 1 ; WhichThetaSlice2D < ThetaSlices2D; WhichThetaSlice2D++ ) {	
		
			int MinThetaSlice2D = MinTheta2D+WhichThetaSlice2D*ThetaStep2D;
			int MaxThetaSlice2D = MinTheta2D+(WhichThetaSlice2D+1)*ThetaStep2D;


			TH1D* h1_ECal_InEePrimeAndThetaSlices = (TH1D*)file_SuSav2->Get(Form("h1_ECal_InEePrime_%d_To_%d_InTheta_%d_To_%d_Slices",\
								MinEePrimeSlice2D,MaxEePrimeSlice2D,MinThetaSlice2D,MaxThetaSlice2D));
			
			
			PrettyPlot(h1_ECal_InEePrimeAndThetaSlices);
			h1_ECal_InEePrimeAndThetaSlices->SetLineColor(Colors[Counter+1]);
			h1_ECal_InEePrimeAndThetaSlices->GetYaxis()->SetRangeUser(0,18);									
			h1_ECal_InEePrimeAndThetaSlices->Draw("hist same");
			
			legSuSav2_EePrimeAndThetaSlices->AddEntry(h1_ECal_InEePrimeAndThetaSlices,\
								ToString(MinEePrimeSlice2D)+ " < E_{e'} < " + ToString(MaxEePrimeSlice2D) + " & " +\
								ToString(MinThetaSlice2D)+ " < #theta_{e'} < " + ToString(MaxThetaSlice2D)\
								,"l");	

			Counter++;
		
		}

	}
	
	legSuSav2_EePrimeAndThetaSlices->Draw();

	// ----------------------------------------------------------------------------------------------------------------	


} // End of the program
