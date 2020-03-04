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

void ApplySystUnc(TH1D* h, double systunc) {

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 
					
		double error = h->GetBinError(i);
		double newerror = error * (1. + systunc);
		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

void OverlayMultiplicities_FigExtData7() {

	// ------------------------------------------------------------------------

	SetOffsetAndSize();
	TGaxis::SetMaxDigits(3);

	int Ndivisions = 3;
	int LineWidth = 3;
	int FontStyle = 132;
	double TextSize = 0.08;
	
	TString version = "v3_0_6/";

	int NECalRebin = 1;

	// From Mariana's analysis note

	double SystUnc1GeV = 0.02; // 2% syst uncertainty at 1.161 GeV
	double SystUnc2GeV = 0.021; // 2.1% syst uncertainty at 2.261 GeV
	double SystUnc4GeV = 0.047; // 4.7% syst uncertainty at 4.461 GeV

	// ------------------------------------------------------------------------

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> JustNucleus; std::vector<TString> LabelsOfSamples; 
	std::vector<TString> E; std::vector<double> DoubleE;
	std::vector<TString> LabelE; std::vector<TString> FSIModel; std::vector<TString> DirNames;  std::vector<int> BreakDownColors;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots; std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames; std::vector<TH1D*> BreakDownPlots;
	std::vector<int> Colors;
	std::vector<int> Style;

//	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He"); JustNucleus.push_back("He");
	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C"); JustNucleus.push_back("C");
//	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");  JustNucleus.push_back("Fe");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); DoubleE.push_back(1.161);
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); DoubleE.push_back(2.261);	
//	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");  DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");
 
	Colors.push_back(kBlack); Colors.push_back(kBlack); Colors.push_back(kBlue); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);

	Style.push_back(1); Style.push_back(1); Style.push_back(1); Style.push_back(1);

	BreakDownColors.push_back(kBlue); BreakDownColors.push_back(kCyan); BreakDownColors.push_back(kGreen); BreakDownColors.push_back(kMagenta);

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");

	NameOfPlots.push_back("h1_Npi"); LabelOfPlots.push_back("Multiplicities"); OutputPlotNames.push_back("PionMultiPlot");
	NameOfPlots.push_back("h1_Nprot"); LabelOfPlots.push_back("Multiplicities"); OutputPlotNames.push_back("Nproton");

	std::vector<TH1D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString WhatModelsAreIncluded = "";
	for (int LoopOverFSIModels = 0 ; LoopOverFSIModels < NFSIModels ; LoopOverFSIModels ++) { WhatModelsAreIncluded += "_"+DirNames[LoopOverFSIModels]; };

	TString RecoCalorimetry = "(e,e'p)";
	TString FSI = "FSI";

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+"Multiplicities"+"_"+xBCut[WhichxBCut],
									 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+"Multiplicities"+"_"+xBCut[WhichxBCut],
									 205,34,1024,768);

				PlotCanvas->SetLeftMargin(0.12);
				PlotCanvas->SetRightMargin(0.04);
				PlotCanvas->SetBottomMargin(0.15);
				PlotCanvas->SetTopMargin(0.0.01);

				PlotCanvas->cd();

//				TLegend* leg = new TLegend(0.2,0.35,0.5,0.55);
//				leg->SetNColumns(2);

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
						CenterAxisTitle(Plots[WhichFSIModel]);

						// --------------------------------------------------------------------------------------

						// X-axis label

						Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleOffset(0.9);

						// X-axis Title

						Plots[WhichFSIModel]->GetXaxis()->SetTitle(LabelOfPlots[WhichPlot]); 

						// X-axis Title/Tick Size

						Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize); 
						Plots[WhichFSIModel]->GetYaxis()->SetTickSize(0.02);

						// --------------------------------------------------------------------------------------
						// Y-axis label

						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);
						Plots[WhichFSIModel]->GetYaxis()->SetTitle("Weighted Events");

						// --------------------------------------------------------------------------------------

						Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.75); 

						Plots[WhichFSIModel]->SetLineWidth(LineWidth);

						// --------------------------------------------------------------------------------------

						// Scaling Factor

						double ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral();
						Plots[WhichFSIModel]->Scale(ScalingFactor);

						// ----------------------------------------------------------------------------------

						// Apply Systematic Uncertainties on Data Points

						double SystUnc = 0;
						if ( DoubleE[WhichEnergy] == 1.161 ) { SystUnc = SystUnc1GeV; }
						if ( DoubleE[WhichEnergy] == 2.261 ) { SystUnc = SystUnc2GeV; }
						if ( DoubleE[WhichEnergy] == 4.461 ) { SystUnc = SystUnc4GeV; }

						if (FSILabel[WhichFSIModel] == "Data") { ApplySystUnc(Plots[WhichFSIModel], SystUnc); }

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						double localmax = Plots[WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
						double height = 1.05;
						if ( xBCut[WhichxBCut] == "xBCut" ) { height = 1.1; }
						Plots[0]->GetYaxis()->SetRangeUser(0.,height*max);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						//Plots[0]->GetXaxis()->SetTitle(XLabel);

						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(5);
						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);

						// --------------------------------------------------------------------------------------------------

						// Multiplicity plots

						if (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" || NameOfPlots[WhichPlot] == "h1_Npi") {

							Plots[WhichFSIModel]->GetYaxis()->SetLabelOffset(-0.004);
							Plots[WhichFSIModel]->Rebin();
							Plots[0]->GetYaxis()->SetRangeUser(0.5*min,2.*max); PlotCanvas->SetLogy();
							if (FSILabel[WhichFSIModel] == "Data") { 
						 
								Plots[WhichFSIModel]->SetMarkerSize(3.); 
								if (NameOfPlots[WhichPlot] == "h1_Nprot") { 
									Plots[WhichFSIModel]->SetMarkerColor(kBlack); Plots[WhichFSIModel]->SetMarkerStyle(20); }
								else { 
									Plots[WhichFSIModel]->SetLineColor(kBlue); 
									Plots[WhichFSIModel]->SetMarkerColor(kBlue); 
									Plots[WhichFSIModel]->SetMarkerStyle(24);
								}
								gStyle->SetErrorX(0); 
								Plots[WhichFSIModel]->GetYaxis()->SetRangeUser(1E1,1E8);
								Plots[WhichFSIModel]->Draw("e same"); 
							}
							else { 
								if (NameOfPlots[WhichPlot] == "h1_Nprot") { Plots[WhichFSIModel]->SetLineColor(kBlack); }
								else { Plots[WhichFSIModel]->SetLineColor(kBlue); Plots[WhichFSIModel]->SetLineStyle(7); }
								Plots[WhichFSIModel]->Draw("hist same"); 
								gStyle->SetErrorX(0); 
								Plots[0]->Draw("e same"); 
							}
						} 

						// ----------------------------------------------------------------------------------------------------

//						if ( FSILabel[WhichFSIModel] == "Data") { leg->AddEntry(Plots[WhichFSIModel]," /","p"); }
//						else { 
//							if (NameOfPlots[WhichPlot] == "h1_Nprot") { leg->AddEntry(Plots[WhichFSIModel]," protons","l"); }
//							else { leg->AddEntry(Plots[WhichFSIModel]," #pi^{#pm}","l"); }
//						}

		                                // --------------------------------------------------------------------------------------------------

					} // End of the loop over the FSI Models 


				} // End of the loop over the plots

				// -----------------------------------------------------------------------------------------------------------------------------------------

//				leg->SetBorderSize(0);
//				leg->SetTextFont(FontStyle);
//				leg->SetTextSize(TextSize);
//				leg->Draw();

				TLatex latexData;
				latexData.SetTextFont(FontStyle);
				latexData.SetTextSize(TextSize);
				latexData.DrawLatexNDC(0.65,0.8,"Protons");

				TLatex latexGenie;
				latexGenie.SetTextFont(FontStyle);
				latexGenie.SetTextColor(kBlue);
				latexGenie.SetTextSize(TextSize);
				latexGenie.DrawLatexNDC(0.5,0.6,"#pi^{#pm}");

				// -----------------------------------------------------------------------------------------------------------------------------------------

				TString ext = "";
				if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

				PlotCanvas->SaveAs("../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]+"/"+ext+nucleus[WhichNucleus]+"_" 
					+E[WhichEnergy]+"_" +"Multiplicities"+WhatModelsAreIncluded+".pdf");

				//delete PlotCanvas;

				// -----------------------------------------------------------------------------------------------------------------------------------------

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
