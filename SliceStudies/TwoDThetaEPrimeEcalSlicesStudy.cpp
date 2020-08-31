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

void TwoDThetaEPrimeEcalSlicesStudy() {

	TFile* file_Data = TFile::Open("../../myFiles/2_261/Data_Final/NoxBCut/12C_2_261_Data_Final_Plots_FSI_em.root","readonly");
	TFile* file_SuSav2 = TFile::Open("../../myFiles/2_261/SuSav2_RadCorr_LFGM/NoxBCut/12C_2_261_SuSav2_RadCorr_LFGM_Plots_FSI_em.root","readonly");	
	
	// ----------------------------------------------------------------------------------------------------------------
	
	std::vector<int> Colors{kYellow,kGreen,kViolet,kBlack,kRed+1,kBlue,kMagenta,410,kYellow+1,kOrange+7};
	
	gStyle->SetOptStat(0);	
	
	double GlobalMax = 15.;

	// ----------------------------------------------------------------------------------------------------------------	

	// SuSav2 ECal in EePrime and ThetaPrime Slices 		

	double MinEePrime2D = 0.5, MaxEePrime2D = 2.5; int EePrimeSlices2D = 3;
	double MinTheta2D = 10, MaxTheta2D = 60; int ThetaSlices2D = 3;	
	double EePrimeStep2D = (MaxEePrime2D - MinEePrime2D) / EePrimeSlices2D;
	int ThetaStep2D = (MaxTheta2D - MinTheta2D) / ThetaSlices2D;

	//int Counter = 0;	
	
	for (int WhichEePrimeSlice2D = 0 ; WhichEePrimeSlice2D < EePrimeSlices2D; WhichEePrimeSlice2D++ ) {

		int MinEePrimeSlice2D = (MinEePrime2D+WhichEePrimeSlice2D*EePrimeStep2D) * 1000;
		int MaxEePrimeSlice2D = (MinEePrime2D+(WhichEePrimeSlice2D+1)*EePrimeStep2D)*1000;

		for (int WhichThetaSlice2D = 0 ; WhichThetaSlice2D < ThetaSlices2D; WhichThetaSlice2D++ ) {	

			if (WhichEePrimeSlice2D == 0 && WhichThetaSlice2D == 0) { continue; }
		
			int MinThetaSlice2D = MinTheta2D+WhichThetaSlice2D*ThetaStep2D;
			int MaxThetaSlice2D = MinTheta2D+(WhichThetaSlice2D+1)*ThetaStep2D;

			TCanvas* SuSav2Canvas_ECal_EePrimeAndThetaSlices = new TCanvas("SuSav2Canvas_ECal_EePrimeAndThetaSlices_MinTheta_"\
											+ToString(MinThetaSlice2D)+"_MinEePrime_"+ToString(MinEePrimeSlice2D),\
										       "SuSav2Canvas_ECal_EePrimeAndThetaSlices_MinTheta_"\
											+ToString(MinThetaSlice2D)+"_MinEePrime_"+ToString(MinEePrimeSlice2D),\
											205,34,1024,768);
			
			SuSav2Canvas_ECal_EePrimeAndThetaSlices->SetBottomMargin(0.15);
			SuSav2Canvas_ECal_EePrimeAndThetaSlices->SetLeftMargin(0.15);


			TH1D* h1_ECal_SuSav2 = (TH1D*)file_SuSav2->Get(Form("h1_ECal_InEePrime_%d_To_%d_InTheta_%d_To_%d_Slices",\
								MinEePrimeSlice2D,MaxEePrimeSlice2D,MinThetaSlice2D,MaxThetaSlice2D));

			TH1D* h1_ECal_Data = (TH1D*)file_Data->Get(Form("h1_ECal_InEePrime_%d_To_%d_InTheta_%d_To_%d_Slices",\
								MinEePrimeSlice2D,MaxEePrimeSlice2D,MinThetaSlice2D,MaxThetaSlice2D));
			

			PrettyPlot(h1_ECal_Data);
			h1_ECal_Data->SetLineColor(kBlack);
			h1_ECal_Data->GetYaxis()->SetRangeUser(0,1.2*h1_ECal_Data->GetMaximum());									
			h1_ECal_Data->Draw("hist same");
			
			PrettyPlot(h1_ECal_SuSav2);
			h1_ECal_SuSav2->SetLineColor(kBlue);
			h1_ECal_SuSav2->Draw("hist same");

			TLegend* legSuSav2_EePrimeAndThetaSlices = new TLegend(0.17,0.79,0.35,0.89);
			legSuSav2_EePrimeAndThetaSlices->SetBorderSize(0);
			legSuSav2_EePrimeAndThetaSlices->SetTextFont(TextFont);
			legSuSav2_EePrimeAndThetaSlices->SetTextSize(TextSize);	
			
			legSuSav2_EePrimeAndThetaSlices->AddEntry(h1_ECal_SuSav2,\
								ToString(MinEePrimeSlice2D)+ " < E_{e'} < " + ToString(MaxEePrimeSlice2D) + " & " +\
								ToString(MinThetaSlice2D)+ " < #theta_{e'} < " + ToString(MaxThetaSlice2D)\
								,"");

			legSuSav2_EePrimeAndThetaSlices->Draw();	

			//Counter++;
		
		}

	}
	

	// ----------------------------------------------------------------------------------------------------------------	


} // End of the program
