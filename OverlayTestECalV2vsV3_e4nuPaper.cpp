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

void OverlayTestECalV2vsV3_e4nuPaper() {

	// ------------------------------------------------------------------------

	SetOffsetAndSize();
	TGaxis::SetMaxDigits(5);
//	TGaxis::SetExponentOffset(-0.1, 1., "y");

	int Ndivisions = 5;
	int LineWidth = 1;
	int FontStyle = 132;
	double TextSize = 0.08;
	
	TString version = "v3_0_6/";

	int NECalRebin = 1;

	// From Mariana's analysis note

	double SystUnc1GeV = 0.02; // 2% syst uncertainty at 1.161 GeV
	double SystUnc2GeV = 0.021; // 2.1% syst uncertainty at 2.261 GeV
	double SystUnc4GeV = 0.047; // 4.7% syst uncertainty at 4.461 GeV

	// ------------------------------------------------------------------------

	// Larry/Axel's suggestion for scaling the last 2 bins by EnhaceTail
//	double EnhaceTail = 1./4.;
//	double EnhaceTail = 1./2.;
	double EnhaceTail = 1./3.;
//	double EnhaceTail = 1.;

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> JustNucleus; std::vector<TString> LabelsOfSamples; 
	std::vector<TString> E; std::vector<double> DoubleE;
	std::vector<TString> LabelE; std::vector<TString> FSIModel; std::vector<TString> DirNames;  std::vector<int> BreakDownColors;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots; std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames; std::vector<TH1D*> BreakDownPlots;
	std::vector<int> Colors;
	std::vector<int> Style;

	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C"); JustNucleus.push_back("C");
	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");  JustNucleus.push_back("Fe");

	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); DoubleE.push_back(1.161);
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); DoubleE.push_back(2.261);	
	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");  DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");
 
	Colors.push_back(kBlack); Colors.push_back(kBlack); Colors.push_back(kBlue); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);

	Style.push_back(1); Style.push_back(1); Style.push_back(1); Style.push_back(1);

	BreakDownColors.push_back(kBlue); BreakDownColors.push_back(kCyan); BreakDownColors.push_back(kGreen); BreakDownColors.push_back(kMagenta);

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("v3");  DirNames.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2015_Final_NoRadCorr_RFGM"); FSILabel.push_back("v2");  DirNames.push_back("hA2015_Truth_NoRadCorr");
//	FSIModel.push_back("hA2018_Final_NoRadCorr"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");

	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> Plots_Clones;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString WhatModelsAreIncluded = "";
	for (int LoopOverFSIModels = 0 ; LoopOverFSIModels < NFSIModels ; LoopOverFSIModels ++) { WhatModelsAreIncluded += "_"+DirNames[LoopOverFSIModels]; };

	TString RecoCalorimetry = "(e,e'p)";
	TString FSI = "FSI";

	std::vector<TString> GenieFSILabel; GenieFSILabel.clear();
	GenieFSILabel.push_back("QE"); GenieFSILabel.push_back("MEC"); GenieFSILabel.push_back("RES"); GenieFSILabel.push_back("DIS");

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		TCanvas* PlotCanvas = new TCanvas(xBCut[WhichxBCut],xBCut[WhichxBCut],205,34,2000,768);
//			205,34,1024,768);
//			205,34,768,768);

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			// Loop over the energies

			for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

				// Loop over the nuclei

				for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

					// ---------------------------------------------------------------------------------------------------------------------------

					// Dimensions of TPads

					double Xmin = 0.07, Xmax = 1.;
					double Ymax = 1., Ymin = 0.1;
					double Xstep = (Xmax - Xmin) / 3.;
					double Ystep = ( Ymax - Ymin  ) / 2.;
					double XMinPad = Xmin + WhichEnergy * Xstep, XMaxPad = Xmin + ( WhichEnergy + 1) * Xstep;
					double YMinPad = Ymax - ( WhichNucleus + 1) * Ystep, YMaxPad = Ymax - WhichNucleus * Ystep;
					double space = 0.07;

					TPad* pad = new TPad(); 

					if (nucleus[WhichNucleus] == "12C") { pad = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPad,YMinPad,XMaxPad,YMaxPad, 21); }
					else { pad = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPad,YMinPad+space,XMaxPad,YMaxPad+space, 21); }

					pad->SetFillColor(kWhite); 
					PlotCanvas->cd();
					pad->Draw();
					pad->cd();

					pad->SetBottomMargin(0.15);
					pad->SetTopMargin(0.);
					pad->SetLeftMargin(0.017);
					pad->SetRightMargin(0.005);

					// ---------------------------------------------------------------------------------------------------------------------------

					// No data on 56Fe @ 1.161 GeV

					if ( nucleus[WhichNucleus] == "56Fe" && DoubleE[WhichEnergy] == 1.161 ) { delete pad; continue; }

					// ---------------------------------------------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* legGenieBlackLine = new TLegend(0.05,0.05,0.95,0.95);
					legGenieBlackLine->SetNColumns(1);

					TLegend* legGenieBreak = new TLegend(0.13,0.25,0.55,0.45);
					legGenieBreak->SetNColumns(2);

					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);
if (!FileSample) { continue; }
						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

if (!Plots[WhichFSIModel]) { continue; }

						Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
						CenterAxisTitle(Plots[WhichFSIModel]);
						Plots[WhichFSIModel]->SetLineWidth(LineWidth);

						// --------------------------------------------------------------------------------------

						// X-axis

						Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(2*TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(0.);
						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);

						// --------------------------------------------------------------------------------------

						// Y-axis label

						Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(0.);
						Plots[WhichFSIModel]->GetYaxis()->SetTitle("");
						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(0.);

						// --------------------------------------------------------------------------------------

						// Scaling Factor

						double ScalingFactor = 1. / Plots[WhichFSIModel]->Integral();  // area normalized
						Plots[WhichFSIModel]->Scale(ScalingFactor);

						//-----------------------------------------------------------------------------------------------

						//Larry's suggestion because ECal has a sharp peak and a low tail 
						//Thus we multiply the peak by EnhaceTail

						if ( DoubleE[WhichEnergy] == 2.261 ) {

							double LowE = 0.95*DoubleE[WhichEnergy];
							int LowEBin = Plots[WhichFSIModel]->FindBin(LowE);
							int HighEBin = Plots[WhichFSIModel]->GetNbinsX();

							for (int i = LowEBin + 1; i <= HighEBin; i++) { 
					
								double content = Plots[WhichFSIModel]->GetBinContent(i);
								double error = Plots[WhichFSIModel]->GetBinError(i);
								double newcontent = EnhaceTail * content;
								double newerror = EnhaceTail * error;				
								Plots[WhichFSIModel]->SetBinContent(i,newcontent);
								Plots[WhichFSIModel]->SetBinError(i,newerror);

							}

						}

						// -----------------------------------------------------------------------------------

						// Accounting for the fact that the bin width might not be constant

						ReweightPlots(Plots[WhichFSIModel]);

						// --------------------------------------------------------------------------------------

						// Rebining & ranges

						if (DoubleE[WhichEnergy] == 1.161) { Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.5,1.3); }
						if (DoubleE[WhichEnergy] == 2.261) { Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.6,2.5); }
						if (DoubleE[WhichEnergy] == 4.461) { Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(1.,4.8); }

						// ----------------------------------------------------------------------------------

						// Apply Systematic Uncertainties on Data Points

						double SystUnc = 0;
						if ( DoubleE[WhichEnergy] == 1.161 ) { SystUnc = SystUnc1GeV; }
						if ( DoubleE[WhichEnergy] == 2.261 ) { SystUnc = SystUnc2GeV; }
						if ( DoubleE[WhichEnergy] == 4.461 ) { SystUnc = SystUnc4GeV; }

						if (FSILabel[WhichFSIModel] == "Data") { ApplySystUnc(Plots[WhichFSIModel], SystUnc); }

						// ----------------------------------------------------------------------------------

						// Genie Break Down

						if (FSILabel[WhichFSIModel] == "Data") legGenieBlackLine->AddEntry(Plots[0],"Data", "lep");
						else { legGenieBlackLine->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l");  }
						//legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 

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
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						// -------------------------------------------------------------------------------------------------

						TLine* line = new TLine(0.95*DoubleE[WhichEnergy],0.,0.95*DoubleE[WhichEnergy],height*max);
						line->SetLineColor(kBlack); 
						line->SetLineWidth(LineWidth);
						if ( FSILabel[WhichFSIModel] == "Genie") { line->Draw(); }

						// --------------------------------------------------------------------------------------------------

						if (FSILabel[WhichFSIModel] == "Data") { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[WhichFSIModel]->Draw("e same"); 

						} else { 

//							Plots[WhichFSIModel]->Draw("hist same"); // draw them as histos
							Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines
							Plots[0]->Draw("e same"); 

						}

					} // End of the loop over the FSI Models 

//					if ( nucleus[WhichNucleus] == "12C" && DoubleE[WhichEnergy] == 1.161) {  

//						legGenieBreak->SetBorderSize(0); 
//						legGenieBreak->SetTextSize(1.5*TextSize);
//						legGenieBreak->SetTextFont(FontStyle);  
//						legGenieBreak->Draw();

//					}

					// -----------------------------------------------------------------------------------------------------------------------------------------

				} // End of the loop over the energies

			} // End of the loop over the nuclei

		} // End of the loop over the plots

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexCarbon;
		latexCarbon.SetTextFont(FontStyle);
		latexCarbon.SetTextSize(TextSize);
		PlotCanvas->cd();
		latexCarbon.DrawLatexNDC(0.005,0.77,"^{12}C");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexIron;
		latexIron.SetTextFont(FontStyle);
		latexIron.SetTextSize(TextSize);
		PlotCanvas->cd();
		latexIron.DrawLatexNDC(0.31,0.37,"^{56}Fe");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* pad1GeV = new TPad("pad1GeV","pad1GeV",0.09,0.85,0.26,0.95,21); 
		pad1GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad1GeV->Draw();
		pad1GeV->cd();

		TLatex latex1GeV;
		latex1GeV.SetTextFont(FontStyle);
		latex1GeV.SetTextSize(10*TextSize);
		latex1GeV.DrawLatexNDC(0.1,0.4,"1.161 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* pad2GeV = new TPad("pad2GeV","pad2GeV",0.4,0.85,0.6,0.95,21); 
		pad2GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad2GeV->Draw();
		pad2GeV->cd();

		TLatex latex2GeV;
		latex2GeV.SetTextFont(FontStyle);
		latex2GeV.SetTextSize(10*TextSize);
		latex2GeV.DrawLatexNDC(0.3,0.4,"2.261 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* pad4GeV = new TPad("pad4GeV","pad4GeV",0.7,0.85,0.9,0.95,21); 
		pad4GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad4GeV->Draw();
		pad4GeV->cd();

		TLatex latex4GeV;
		latex4GeV.SetTextFont(FontStyle);
		latex4GeV.SetTextSize(10*TextSize);
		latex4GeV.DrawLatexNDC(0.3,0.4,"4.461 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* padPmiss = new TPad("padPmiss","padPmiss",0.59,0.03,0.79,0.13,21); 
		padPmiss->SetFillColor(kWhite); 
		PlotCanvas->cd();
		padPmiss->Draw();
		padPmiss->cd();

		TLatex latexPmiss;
		latexPmiss.SetTextFont(FontStyle);
		latexPmiss.SetTextSize(10*TextSize);
		latexPmiss.DrawLatexNDC(0.2,0.3,"E^{Cal} [GeV]");

		// -------------------------------------------------------------------------------------------------
	
		TPad* padx13C12 = new TPad("padx13C12","padx13C12",0.65,0.7,0.688,0.8,21); 
		padx13C12->SetFillColor(kWhite); 
		PlotCanvas->cd();
		padx13C12->Draw();
		padx13C12->cd();

		TLatex latexx13C12;
		latexx13C12.SetTextFont(FontStyle);
		latexx13C12.SetTextSize(8*TextSize);
		latexx13C12.DrawLatexNDC(0.,0.4,"x1/3");

		// -------------------------------------------------------------------------------------------------
	
		TPad* padx13Fe56 = new TPad("padx13Fe56","padx13Fe56",0.65,0.3,0.688,0.4,21); 
		padx13Fe56->SetFillColor(kWhite); 
		PlotCanvas->cd();
		padx13Fe56->Draw();
		padx13Fe56->cd();

		TLatex latexx13F56;
		latexx13F56.SetTextFont(FontStyle);
		latexx13F56.SetTextSize(8*TextSize);
		latexx13F56.DrawLatexNDC(0.,0.4,"x1/3");

		// -------------------------------------------------------------------------------------------------
	
		TPad* padLeg = new TPad("padLeg","padLeg",0.1,0.2,0.3,0.5,21); 
		padLeg->SetFillColor(kWhite); 
		PlotCanvas->cd();
		padLeg->Draw();
		padLeg->cd();

		legGenieBlackLine->SetBorderSize(0); 
		legGenieBlackLine->SetNColumns(1); 
		legGenieBlackLine->SetTextSize(3.*TextSize); 
		legGenieBlackLine->SetTextFont(FontStyle); 
		legGenieBlackLine->Draw();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString ext = "";
		if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

		PlotCanvas->SaveAs("../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+ext+"TestGenieV2vsV3"+WhatModelsAreIncluded+".pdf");

		//delete PlotCanvas;

	} // End of the loop over the xB kinematics

} // End of the program
