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

void OverlayPmiss_FigExtData10() {

	// ------------------------------------------------------------------------

	SetOffsetAndSize();
	TGaxis::SetMaxDigits(5);

	int Ndivisions = 5;
	int LineWidth = 1;
	int FontStyle = 132;
	double TextSize = 0.08;
	
	TString version = "v3_0_6/";

	// From Mariana's analysis note

	double SystUnc1GeV = 0.02; // 2% syst uncertainty at 1.161 GeV
	double SystUnc2GeV = 0.021; // 2.1% syst uncertainty at 2.261 GeV
	double SystUnc4GeV = 0.047; // 4.7% syst uncertainty at 4.461 GeV

	// ------------------------------------------------------------------------

	// Larry/Axel's suggestion for scaling the last 2 bins by EnhaceTail

	double EnhaceTail = 1.;

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

//	BreakDownColors.push_back(kBlue); BreakDownColors.push_back(kCyan); BreakDownColors.push_back(kGreen); BreakDownColors.push_back(kMagenta);
	BreakDownColors.push_back(kBlue); BreakDownColors.push_back(429); BreakDownColors.push_back(410); BreakDownColors.push_back(610);

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");
//	FSIModel.push_back("hA2018_Final_NoRadCorr"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");

	NameOfPlots.push_back("MissMomentum"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} P_{miss}^{#perp} [GeV/c]"); OutputPlotNames.push_back("MissMomentum");

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

		TCanvas* PlotCanvas = new TCanvas(xBCut[WhichxBCut],xBCut[WhichxBCut],205,34,1600,900);
//			205,34,1024,768);
//			205,34,768,768);

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			// Loop over the energies

			for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

double MaxHeight = 22; // In order to use y-axis ticks with common scale, constraint range between (0,MaxHeight)

				// Loop over the nuclei

				for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

if (nucleus[WhichNucleus] == "56Fe") { MaxHeight = MaxHeight * 0.55; }

					// ---------------------------------------------------------------------------------------------------------------------------

					// Dimensions of TPads

					double Xmin = 0.14, Xmax = 1.;
					double Ymax = 1., Ymin = 0.05;
					double Xstep = (Xmax - Xmin) / 3.;
					double Ystep = ( Ymax - Ymin  ) / 2.;
					double XMinPad = Xmin + WhichEnergy * Xstep, XMaxPad = Xmin + ( WhichEnergy + 1) * Xstep;
if (DoubleE[WhichEnergy] == 1.161 ) { XMinPad = XMinPad - 0.05; }
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
					pad->SetTopMargin(0.0);
					if (nucleus[WhichNucleus] == "12C") { pad->SetTopMargin(0.01); }
					pad->SetLeftMargin(0.);
if (DoubleE[WhichEnergy] == 1.161 ) { pad->SetLeftMargin(0.08); }
					pad->SetRightMargin(0.0);
					if (DoubleE[WhichEnergy] == 4.461 ) { pad->SetRightMargin(0.01); }
					pad->SetFrameBorderSize(10);

					// ---------------------------------------------------------------------------------------------------------------------------

					// No data on 56Fe @ 1.161 GeV

					if ( nucleus[WhichNucleus] == "56Fe" && DoubleE[WhichEnergy] == 1.161 ) { delete pad; continue; }

					// ---------------------------------------------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* legGenieBlackLine = new TLegend(0.1,0.5,0.54,1.);
					legGenieBlackLine->SetNColumns(1);
					legGenieBlackLine->SetTextFont(FontStyle); 

					TLegend* legGenieBreak = new TLegend(0.1,0.,1.,0.5);
					legGenieBreak->SetNColumns(2);
					legGenieBreak->SetTextFont(FontStyle);

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
						Plots[WhichFSIModel]->SetLineWidth(LineWidth);

						// --------------------------------------------------------------------------------------

						// X-axis

						Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(1.2*TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(0.);
						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);

						// --------------------------------------------------------------------------------------

						// Y-axis label

						Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(0.);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelOffset(0.013);
//if (DoubleE[WhichEnergy] == 1.161) { Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(0.07); }	
Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(1.2*TextSize);


						Plots[WhichFSIModel]->GetYaxis()->SetTitle("");

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

//						if (DoubleE[WhichEnergy] == 1.161) { Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.6,1.23); }
//						if (DoubleE[WhichEnergy] == 2.261) { Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.7,2.4); }
//						if (DoubleE[WhichEnergy] == 4.461) { Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(1.3,4.6); }

for (int i = 0; i < 2; i++) { Plots[WhichFSIModel]->Rebin(); }
Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.03,0.95);
Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(3);
						// ----------------------------------------------------------------------------------

						// Apply Systematic Uncertainties on Data Points

						double SystUnc = 0;
						if ( DoubleE[WhichEnergy] == 1.161 ) { SystUnc = SystUnc1GeV; }
						if ( DoubleE[WhichEnergy] == 2.261 ) { SystUnc = SystUnc2GeV; }
						if ( DoubleE[WhichEnergy] == 4.461 ) { SystUnc = SystUnc4GeV; }

						if (FSILabel[WhichFSIModel] == "Data") { ApplySystUnc(Plots[WhichFSIModel], SystUnc); }

						// ----------------------------------------------------------------------------------

						// Genie Break Down

						if (
							FSILabel[WhichFSIModel] == "Genie"
						) {

							legGenieBlackLine->AddEntry(Plots[0],"Data", "lep"); 
							legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 

							BreakDownPlots.clear();

							for (int j = 1; j < 5; j++) {

								BreakDownPlots.push_back( (TH1D*)( FileSample->Get("Pmiss_Int_"+ToString(j)) ) );
								ReweightPlots(BreakDownPlots[j-1]);

								for (int i = 0; i < 2; i++) { BreakDownPlots[j-1]->Rebin(); }

								//-----------------------------------------------------------------------------------------------

								BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j-1]);

//								int GenieNBins = Plots[WhichFSIModel]->GetNbinsX();
//								int GenieMin = Plots[WhichFSIModel]->GetXaxis()->GetXmin();
//								int GenieMax = Plots[WhichFSIModel]->GetXaxis()->GetXmax();
//								BreakDownPlots[j-1]->SetBins(GenieNBins,GenieMin,GenieMax);
								
								BreakDownPlots[j-1]->SetLineWidth(LineWidth);
								BreakDownPlots[j-1]->SetLineStyle(Style[j-1]);
								BreakDownPlots[j-1]->Scale(ScalingFactor);

								//-----------------------------------------------------------------------------------------------

								//Larry's suggestion because ECal has a sharp peak and a low tail 
								//Thus we multiply the peak by EnhaceTail

								if ( DoubleE[WhichEnergy] == 2.261 ) {

									double LowE = 0.95*DoubleE[WhichEnergy];
									int LowEBin = Plots[WhichFSIModel]->FindBin(LowE);
									int HighEBin = Plots[WhichFSIModel]->GetNbinsX();

									for (int i = LowEBin+1; i <= HighEBin; i++) { 
					
										double content = BreakDownPlots[j-1]->GetBinContent(i);
										double error = BreakDownPlots[j-1]->GetBinError(i);
										double newcontent = EnhaceTail * content;
										double newerror = EnhaceTail * error;				
										BreakDownPlots[j-1]->SetBinContent(i,newcontent);
										BreakDownPlots[j-1]->SetBinError(i,newerror);

									}

								}

								//-----------------------------------------------------------------------------------------------

								TLegendEntry* l1Break = legGenieBreak->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								l1Break->SetTextColor(BreakDownColors[j-1]);

//								BreakDownPlots[j-1]->Draw("hist same");
								BreakDownPlots[j-1]->Draw("C hist same");

							} // end of the look over the GENIE break down

						}

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

//						double localmax = Plots[WhichFSIModel]->GetMaximum();
//						if (localmax > max) { max = localmax; }
//						double height = 1.15;
//						if ( xBCut[WhichxBCut] == "xBCut" ) { height = 1.1; }
//						if ( DoubleE[WhichEnergy] == 2.261 ) { height = 1.4; }
//						Plots[0]->GetYaxis()->SetRangeUser(0.,height*max);


Plots[0]->GetYaxis()->SetRangeUser(0.,MaxHeight);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						// -------------------------------------------------------------------------------------------------

//						TLine* line = new TLine(0.95*DoubleE[WhichEnergy],0.,0.95*DoubleE[WhichEnergy],height*max);
TLine* line = new TLine(0.95*DoubleE[WhichEnergy],0.,0.95*DoubleE[WhichEnergy],MaxHeight);
						line->SetLineColor(kBlack); 
						line->SetLineWidth(LineWidth);
						if ( FSILabel[WhichFSIModel] == "Genie" && DoubleE[WhichEnergy] == 2.261) { line->Draw(); }

						// --------------------------------------------------------------------------------------------------

						if (FSILabel[WhichFSIModel] == "Data") { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(1.5); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[WhichFSIModel]->Draw("e same"); 

						} else { 

							Plots[WhichFSIModel]->Draw("C hist same");  // "C hist same" draw them as lines // "hist same" draw them as histos
							Plots[0]->Draw("e same"); 

						}

					} // End of the loop over the FSI Models 

					// -----------------------------------------------------------------------------------------------------------------------------------------

				} // End of the loop over the energies

			} // End of the loop over the nuclei

		} // End of the loop over the plots

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexCarbon;
		latexCarbon.SetTextFont(FontStyle);
		latexCarbon.SetTextSize(0.8*TextSize);
		PlotCanvas->cd();
		latexCarbon.DrawLatexNDC(0.005,0.77,"^{12}C");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TLatex latexIron;
		latexIron.SetTextFont(FontStyle);
		latexIron.SetTextSize(0.8*TextSize);
		PlotCanvas->cd();
		latexIron.DrawLatexNDC(0.3,0.385,"^{56}Fe");

		// -----------------------------------------------------------------------------------------------------------------------------------------

//		TPad* pad1GeV = new TPad("pad1GeV","pad1GeV",0.125,0.89,0.315,0.99,21);
		TPad* pad1GeV = new TPad("pad1GeV","pad1GeV",0.225,0.89,0.415,0.99,21);
		pad1GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad1GeV->Draw();
		pad1GeV->cd();

		TLatex latex1GeV;
		latex1GeV.SetTextFont(FontStyle);
		latex1GeV.SetTextSize(5*TextSize);
		latex1GeV.DrawLatexNDC(0.04,0.45,"1.159 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

//		TPad* pad2GeV = new TPad("pad2GeV","pad2GeV",0.45,0.89,0.6,0.99,21); 
		TPad* pad2GeV = new TPad("pad2GeV","pad2GeV",0.53,0.89,0.68,0.99,21); 
		pad2GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad2GeV->Draw();
		pad2GeV->cd();

		TLatex latex2GeV;
		latex2GeV.SetTextFont(FontStyle);
		latex2GeV.SetTextSize(5*TextSize);
		latex2GeV.DrawLatexNDC(0.06,0.45,"2.257 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

//		TPad* pad4GeV = new TPad("pad4GeV","pad4GeV",0.73,0.89,0.9,0.99,21); 
		TPad* pad4GeV = new TPad("pad4GeV","pad4GeV",0.8,0.89,0.97,0.99,21); 
		pad4GeV->SetFillColor(kWhite); 
		PlotCanvas->cd();
		pad4GeV->Draw();
		pad4GeV->cd();

		TLatex latex4GeV;
		latex4GeV.SetTextFont(FontStyle);
		latex4GeV.SetTextSize(5*TextSize);
		latex4GeV.DrawLatexNDC(0.11,0.45,"4.453 GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TPad* padPmiss = new TPad("padPmiss","padPmiss",0.405,0.01,0.765,0.14,21); 
		padPmiss->SetFillColor(kWhite); 
		PlotCanvas->cd();
		padPmiss->Draw();
		padPmiss->cd();

		TLatex latexPmiss;
		latexPmiss.SetTextFont(FontStyle);
		latexPmiss.SetTextSize(5*TextSize);
		latexPmiss.DrawLatexNDC(0.1,0.5,"(e,e'p)_{1p0#pi} P_{miss}^{#perp} [GeV/c]");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the legend

		PlotCanvas->cd();
		TPad* padLegend = new TPad("padLegend","padLegend",0.07,0.19,0.3,0.53, 21); 
		padLegend->SetFillColor(kWhite); 
		padLegend->Draw();
		padLegend->cd();

		legGenieBlackLine->SetTextSize(2.*TextSize); 
		legGenieBlackLine->SetBorderSize(0); 
		legGenieBlackLine->Draw();

		legGenieBreak->SetTextSize(2.*TextSize); 
		legGenieBreak->SetBorderSize(0); 
		legGenieBreak->Draw();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis units carbon

		PlotCanvas->cd();
		TPad* padTitle = new TPad("padTitle","padTitle",0.065,0.58,0.09,1., 21); 
		padTitle->SetFillColor(kWhite); 
		padTitle->Draw();
		padTitle->cd();

		TLatex latexYTitle;
		latexYTitle.SetTextFont(FontStyle);
		latexYTitle.SetTextSize(11*TextSize);
		latexYTitle.SetTextColor(kBlack);
		latexYTitle.SetTextAngle(90);
		latexYTitle.DrawLatexNDC(0.8,0.03,"Weighted Events / GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis units iron

		PlotCanvas->cd();
		TPad* padTitleFe = new TPad("padTitleFe","padTitleFe",0.37,0.18,0.4,0.55,21); 
		padTitleFe->SetFillColor(kWhite); 
		padTitleFe->Draw();
		padTitleFe->cd();

		TLatex latexYTitleFe;
		latexYTitleFe.SetTextFont(FontStyle);
		latexYTitleFe.SetTextSize(9*TextSize);
		latexYTitleFe.SetTextColor(kBlack);
		latexYTitleFe.SetTextAngle(90);
		latexYTitleFe.DrawLatexNDC(0.8,0.03,"Weighted Events / GeV");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis 0. point

		PlotCanvas->cd();
		TPad* padTitleZero = new TPad("padTitleZero","padTitleZero",0.41,0.17,0.426,0.22,21); 
		padTitleZero->SetFillColor(kWhite); 
		padTitleZero->Draw();
		padTitleZero->cd();

		TLatex latexYTitleZero;
		latexYTitleZero.SetTextFont(FontStyle);
		latexYTitleZero.SetTextSize(20*TextSize);
		latexYTitleZero.SetTextColor(kBlack);
		latexYTitleZero.DrawLatexNDC(0.,0.1,"0");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis 1. point

		PlotCanvas->cd();
		TPad* padTitleOne = new TPad("padTitleOne","padTitleOne",0.41,0.34,0.426,0.39,21); 
		padTitleOne->SetFillColor(kWhite); 
		padTitleOne->Draw();
		padTitleOne->cd();

		TLatex latexYTitleOne;
		latexYTitleOne.SetTextFont(FontStyle);
		latexYTitleOne.SetTextSize(20*TextSize);
		latexYTitleOne.SetTextColor(kBlack);
		latexYTitleOne.DrawLatexNDC(0.,0.1,"5");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the Y-axis 2. point

		PlotCanvas->cd();
		TPad* padTitleTwo = new TPad("padTitleTwo","padTitleTwo",0.40,0.505,0.425,0.555,21); 
		padTitleTwo->SetFillColor(kWhite); 
		padTitleTwo->Draw();
		padTitleTwo->cd();

		TLatex latexYTitleTwo;
		latexYTitleTwo.SetTextFont(FontStyle);
		latexYTitleTwo.SetTextSize(13*TextSize);
		latexYTitleTwo.SetTextColor(kBlack);
		latexYTitleTwo.DrawLatexNDC(0.,0.1,"10");

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the lower X-axis to cover half zeros

		PlotCanvas->cd();
		TPad* padWhitePadOne = new TPad("padWhitePadOne","padWhitePadOne",0.425,0.15,0.445,0.185,21); 
		padWhitePadOne->SetFillColor(kWhite); 
		padWhitePadOne->Draw();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the lower X-axis to cover half zeros

		PlotCanvas->cd();
		TPad* padWhitePadTwo = new TPad("padWhitePadTwo","padWhitePadTwo",0.7,0.15,0.72,0.185,21); 
		padWhitePadTwo->SetFillColor(kWhite); 
		padWhitePadTwo->Draw();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		// Extra pad for the lower X-axis to cover half zeros

		PlotCanvas->cd();
		TPad* padWhitePadThree = new TPad("padWhitePadThree","padWhitePadThree",0.111,0.551,0.131,0.586,21); 
		padWhitePadThree->SetFillColor(kWhite); 
		padWhitePadThree->Draw();

		// -----------------------------------------------------------------------------------------------------------------------------------------

		TString ext = "";
		if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

		PlotCanvas->SaveAs("../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+ext+"FigExtData10"+WhatModelsAreIncluded+".pdf");

		//delete PlotCanvas;

	} // End of the loop over the xB kinematics

} // End of the program
