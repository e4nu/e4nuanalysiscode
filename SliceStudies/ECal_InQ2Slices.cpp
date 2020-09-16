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
		h->GetXaxis()->SetTitle("");
		h->GetXaxis()->SetTitleSize(TextSize);
		h->GetXaxis()->SetLabelSize(0.);
		h->GetXaxis()->SetTitleFont(TextFont);		
		h->GetXaxis()->SetLabelFont(TextFont);
		h->GetXaxis()->SetNdivisions(4);

		h->GetYaxis()->CenterTitle();		
		h->GetYaxis()->SetTitle("");
		h->GetYaxis()->SetTitleSize(TextSize);
		h->GetYaxis()->SetLabelSize(2*TextSize);
		h->GetYaxis()->SetTitleFont(TextFont);		
		h->GetYaxis()->SetLabelFont(TextFont);
		h->GetYaxis()->SetNdivisions(5);


}

// ----------------------------------------------------------------------------------------------------------------

void ECal_InQ2Slices() {

	// ----------------------------------------------------------------------------------------------------------------
	
//	std::vector<int> Colors{kYellow,kGreen,kViolet,kBlack,kRed+1,kBlue,kMagenta,410,kYellow+1,kOrange+7};
	
	gStyle->SetOptStat(0);	
	
	double GlobalMax = 12.;

	int FontStyle = 132;
	double TextSize = 0.08;

	// ---------------------------------------------------------------------------------------------------------------

	std::vector<TString> NucleusString; NucleusString.clear();
	std::vector<TString> NucleusLatex; NucleusLatex.clear();

	NucleusString.push_back("4He"); NucleusLatex.push_back("^{4}He");
	NucleusString.push_back("12C"); NucleusLatex.push_back("^{12}C");
	NucleusString.push_back("56Fe"); NucleusLatex.push_back("^{56}Fe");

	const int NNuclei = NucleusString.size();

	// Loop over the nuclei

	for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus++) {

		// ---------------------------------------------------------------------------------------------------------------

		std::vector<TString> EnergyString; EnergyString.clear();
		std::vector<double> EnergyDouble; EnergyDouble.clear();

		EnergyString.push_back("1_161"); EnergyDouble.push_back(1.161);
		EnergyString.push_back("2_261"); EnergyDouble.push_back(2.261);
		EnergyString.push_back("4_461"); EnergyDouble.push_back(4.461);

		const int NEnergies = EnergyString.size();

		// ----------------------------------------------------------------------------------------------------------------	

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy++) {

			// ----------------------------------------------------------------------------------------------------------------	

			// Data set at 1.1 GeV only on 12C

			if (EnergyDouble[WhichEnergy] == 1.161 && NucleusString[WhichNucleus] != "12C") { continue; }

			// ----------------------------------------------------------------------------------------------------------------	

			// SuSav2 ECal in EePrime and ThetaPrime Slices 		

			double minQ2 = 0.4, maxQ2 = 1.75; int Q2Slices = 9;

			if(EnergyDouble[WhichEnergy]>1. && EnergyDouble[WhichEnergy]<2.) { minQ2 = 0.1; maxQ2 = 0.82; }
			if(EnergyDouble[WhichEnergy]>4. && EnergyDouble[WhichEnergy]<5.) { minQ2 = 0.8; maxQ2 = 4.4; }

			double Q2Step = (maxQ2 - minQ2) / Q2Slices;

			// ---------------------------------------------------------------------------------------------------------------

			// Dimensions of TPads

			double Xmin = 0.07, Xmax = 1.;
			double Ymax = 1., Ymin = 0.08;
			int PadOfHowMany = 3;
			double Xstep = (Xmax - Xmin) / double(int(Q2Slices/PadOfHowMany));
			double Ystep = ( Ymax - Ymin  ) / double(int(Q2Slices/PadOfHowMany));

			// ---------------------------------------------------------------------------------------------------------------

			TFile* file_Data = TFile::Open("../../myFiles/"+EnergyString[WhichEnergy]+"/Data_Final/NoxBCut/"+\
							NucleusString[WhichNucleus]+"_"+EnergyString[WhichEnergy]+"_Data_Final_Plots_FSI_em.root","readonly");
			TFile* file_SuSav2 = TFile::Open("../../myFiles/"+EnergyString[WhichEnergy]+"/SuSav2_RadCorr_LFGM/NoxBCut/"+\
							 NucleusString[WhichNucleus]+"_"+EnergyString[WhichEnergy]+"_SuSav2_RadCorr_LFGM_Plots_FSI_em.root","readonly");	
			TFile* file_G2018 = TFile::Open("../../myFiles/"+EnergyString[WhichEnergy]+"/hA2018_Final_RadCorr_LFGM/NoxBCut/"+\
							 NucleusString[WhichNucleus]+"_"+EnergyString[WhichEnergy]+"_hA2018_Final_RadCorr_LFGM_Plots_FSI_em.root","readonly");	

			TString CanvasName = NucleusString[WhichNucleus]+"_"+EnergyString[WhichEnergy]+"_SuSav2Canvas_ECal_Q2Slices";

			TCanvas* SuSav2Canvas_ECal_Q2Slices = new TCanvas(CanvasName,CanvasName,205,34,1300,1300);

			for (int WhichQ2Slice = 0 ; WhichQ2Slice < Q2Slices; WhichQ2Slice++ ) {	

				int MinQ2Slice = (minQ2+WhichQ2Slice*Q2Step)*1000.;
				int MaxQ2Slice = (minQ2+(WhichQ2Slice+1)*Q2Step)*1000.;

				TPad* pad = new TPad(); 

				TString TPadName = CanvasName+"_Q2Slice_"+ToString(WhichQ2Slice);
	
				SuSav2Canvas_ECal_Q2Slices->cd();
				if (WhichQ2Slice == 0) { pad = new TPad(TPadName,TPadName,0.25,0.75,0.5,1, 21); }
				if (WhichQ2Slice == 1) { pad = new TPad(TPadName,TPadName,0.5,0.75,0.75,1, 21); }
				if (WhichQ2Slice == 2) { pad = new TPad(TPadName,TPadName,0.75,0.75,1.,1, 21); }

				if (WhichQ2Slice == 3) { pad = new TPad(TPadName,TPadName,0.25,0.5,0.5,0.75,21); }
				if (WhichQ2Slice == 4) { pad = new TPad(TPadName,TPadName,0.5,0.5,0.75,0.75,21); }
				if (WhichQ2Slice == 5) { pad = new TPad(TPadName,TPadName,0.75,0.5,1.,0.75,21); }

				if (WhichQ2Slice == 6) { pad = new TPad(TPadName,TPadName,0.25,0.25,0.5,0.5,21); }
				if (WhichQ2Slice == 7) { pad = new TPad(TPadName,TPadName,0.5,0.25,0.75,0.5,21); }
				if (WhichQ2Slice == 8) { pad = new TPad(TPadName,TPadName,0.75,0.25,1.,0.5,21); }

				pad->SetFillColor(kWhite); 
				SuSav2Canvas_ECal_Q2Slices->cd();
				pad->Draw();
				pad->cd();

				pad->SetBottomMargin(0.0);
				if (WhichQ2Slice > 5) { pad->SetBottomMargin(0.14); }
				pad->SetTopMargin(0.0);
				if (WhichQ2Slice < 3) { pad->SetTopMargin(0.01); }
				pad->SetLeftMargin(0.);
				if (WhichQ2Slice == 0 || WhichQ2Slice == 3 || WhichQ2Slice == 6 ) { pad->SetLeftMargin(0.09); }
				pad->SetRightMargin(0.0);
				if (WhichQ2Slice == 2 || WhichQ2Slice == 5 || WhichQ2Slice == 8 ) { pad->SetRightMargin(0.03); }
				pad->SetFrameBorderSize(10);

				TH1D* h1_ECal_G2018 = (TH1D*)file_G2018->Get(Form("h1_ECal_InQ2_%d_To_%d_Slices",MinQ2Slice,MaxQ2Slice));

				TH1D* h1_ECal_SuSav2 = (TH1D*)file_SuSav2->Get(Form("h1_ECal_InQ2_%d_To_%d_Slices",MinQ2Slice,MaxQ2Slice));

				TH1D* h1_ECal_Data = (TH1D*)file_Data->Get(Form("h1_ECal_InQ2_%d_To_%d_Slices",MinQ2Slice,MaxQ2Slice));
					

				double LowRange = -1., HighRange = -1.;

				if (EnergyDouble[WhichEnergy] == 1.161) { LowRange = 0.55; HighRange = 1.29; }
				if (EnergyDouble[WhichEnergy] == 2.261) { LowRange = 0.67; HighRange = 2.4; }
				if (EnergyDouble[WhichEnergy] == 4.461) { LowRange = 1.5; HighRange = 4.6; }

				h1_ECal_Data->GetXaxis()->SetRangeUser(LowRange,HighRange);

				PrettyPlot(h1_ECal_Data);
				h1_ECal_Data->SetLineColor(kBlack);
				h1_ECal_Data->SetMarkerColor(kBlack);
				h1_ECal_Data->SetMarkerStyle(20);
				h1_ECal_Data->SetMarkerSize(1.2);
				h1_ECal_Data->GetYaxis()->SetRangeUser(0.01,GlobalMax);	

				if (WhichQ2Slice > 6) {	h1_ECal_Data->GetXaxis()->SetLabelSize(2*TextSize); }

				if (EnergyString[WhichEnergy] == "2_261") { h1_ECal_Data->GetXaxis()->SetRangeUser(0.5,2.5); }
										
				h1_ECal_Data->Draw("e same");

				PrettyPlot(h1_ECal_SuSav2);
				h1_ECal_SuSav2->SetLineColor(kBlack);
				h1_ECal_SuSav2->SetLineStyle(kSolid);
				h1_ECal_SuSav2->Draw("hist C same");

				PrettyPlot(h1_ECal_G2018);
				h1_ECal_G2018->SetLineColor(kBlack);
				h1_ECal_G2018->SetLineStyle(kDashed);
				h1_ECal_G2018->Draw("hist C same");

				// -------------------------------------------------------------------------------------------------------------

				TLegend* legSuSav2_Slice = new TLegend(0.11,0.62,0.15,0.92);
				legSuSav2_Slice->SetBorderSize(0);
				legSuSav2_Slice->SetTextFont(TextFont);
				legSuSav2_Slice->SetTextSize(1.2*TextSize);	

				legSuSav2_Slice->AddEntry(h1_ECal_Data,ToString(MinQ2Slice)+ " < Q^{2} < " + ToString(MaxQ2Slice)+" MeV^{2}/c^{2}","");				

				legSuSav2_Slice->Draw();

				// -------------------------------------------------------------------------------------------------------------

				if (WhichQ2Slice == 4) {

					TLegend* legSuSav2_EePrimeAndCosThetaSlices = new TLegend(0.1,0.2,0.35,0.6);
					legSuSav2_EePrimeAndCosThetaSlices->SetBorderSize(0);
					legSuSav2_EePrimeAndCosThetaSlices->SetTextFont(TextFont);
					legSuSav2_EePrimeAndCosThetaSlices->SetTextSize(1.5*TextSize);	

					legSuSav2_EePrimeAndCosThetaSlices->AddEntry(h1_ECal_Data,"Data","lep");				
					legSuSav2_EePrimeAndCosThetaSlices->AddEntry(h1_ECal_SuSav2,"SuSav2","l");
					legSuSav2_EePrimeAndCosThetaSlices->AddEntry(h1_ECal_G2018,"G2018","l");

					legSuSav2_EePrimeAndCosThetaSlices->Draw();	

					TLatex latexnucleus;
					latexnucleus.SetTextFont(FontStyle);
					latexnucleus.SetTextSize(1.5*TextSize);
					latexnucleus.DrawLatexNDC(0.55,0.35,NucleusLatex[WhichNucleus]);

				}

				TLine* line = new TLine(EnergyDouble[WhichEnergy],0.01,EnergyDouble[WhichEnergy],GlobalMax);
				line->SetLineWidth(3);
				line->SetLineStyle(5);
				line->Draw();

			} // End of the loop over the Ee' slices

			// ----------------------------------------------------------------------------------------------------------------	

			SuSav2Canvas_ECal_Q2Slices->cd();
			TPad* padWeightedEventsPerGeV = new TPad("padWeightedEventsPerGeV","padWeightedEventsPerGeV",0.2,0.3,0.23,0.8,21); 
			padWeightedEventsPerGeV->SetFillColor(kWhite); 
			padWeightedEventsPerGeV->Draw();
			padWeightedEventsPerGeV->cd();

			TLatex latexYTitlepadWeightedEventsPerGeV;
			latexYTitlepadWeightedEventsPerGeV.SetTextFont(FontStyle);
			latexYTitlepadWeightedEventsPerGeV.SetTextSize(13*TextSize);
			latexYTitlepadWeightedEventsPerGeV.SetTextColor(kBlack);
			latexYTitlepadWeightedEventsPerGeV.SetTextAngle(90);
			latexYTitlepadWeightedEventsPerGeV.DrawLatexNDC(0.8,0.03,"Weighted Events / GeV");

			// ----------------------------------------------------------------------------------------------------------------	

			SuSav2Canvas_ECal_Q2Slices->cd();
			TPad* padEcal = new TPad("padEcal","padEcal",0.45,0.15,0.8,0.22,21); 
			padEcal->SetFillColor(kWhite); 
			padEcal->Draw();
			padEcal->cd();

			TLatex latexpadEcal;
			latexpadEcal.SetTextFont(FontStyle);
			latexpadEcal.SetTextSize(9*TextSize);
			latexpadEcal.SetTextColor(kBlack);
			latexpadEcal.DrawLatexNDC(0.1,0.4,"(e,e'p)_{1p0#pi} E^{cal} [GeV]");
			
			// ----------------------------------------------------------------------------------------------------------------	

//			SuSav2Canvas_ECal_Q2Slices->SaveAs("./myPlots/"+CanvasName+".pdf");

			// ----------------------------------------------------------------------------------------------------------------

		} // End of the loop over the energies

		// ----------------------------------------------------------------------------------------------------------------

	} // End of the loop over the nuclei	

} // End of the program
