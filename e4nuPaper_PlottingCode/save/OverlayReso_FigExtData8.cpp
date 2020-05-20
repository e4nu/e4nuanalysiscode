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

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"

void OverlayReso_FigExtData8() {


	SetOffsetAndSize();

	int Ndivisions = 3;
	int LineWidth = 5;
	double TextSize = 0.07;
	
	TString version = "v3_0_6/";

	int FontStyle = 132;
	TGaxis::SetMaxDigits(3);
	TGaxis::SetExponentOffset(-0.1, 1., "y");
	gStyle->SetLineStyleString(11,"80 60");
	
	TString version = "v3_0_6/";

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> LabelsOfSamples; std::vector<TString> E;
	std::vector<TString> LabelE; std::vector<TString> JustE; std::vector<TString> FSIModel; std::vector<TString> OutputPlotNames;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots;  std::vector<TString> XLabels;

//	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C");
	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); JustE.push_back("1.159 GeV");
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); JustE.push_back("2.257 GeV");
	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV"); JustE.push_back("4.453 GeV");
 
	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie");

	xBCut.push_back("NoxBCut");
	NameOfPlots.push_back("h_Etot_subtruct_piplpimi_factor_fracfeed"); OutputPlotNames.push_back("EcalReso");
	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_factor_fracfeed"); OutputPlotNames.push_back("EQEReso"); 

	XLabels.push_back("E^{cal} Feeddown"); 
	XLabels.push_back("E^{QE} Feeddown"); 

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	const int NEnergies = E.size();
	const int NFSIModels = FSIModel.size();
	const int NPlots = NameOfPlots.size();

	TString RecoCalorimetry = "(e,e'p)";
	TString FSI = "FSI";

	TH1D* Plots[NEnergies][NFSIModels];

	// Larry's suggestion following Barak's paper for color skim

	// 12C
//	int Colors[NEnergies][NFSIModels] = {{kRed,kRed}{kGreen-3,kGreen-3}{kBlue,kBlue}};
//	int LineStyle[NEnergies] = {2,11,1};
//	int MarkerStyle[NEnergies] = {22,21,20};

	// 56Fe
	int Colors[NEnergies][NFSIModels] = {{kGreen-3,kGreen-3}{kBlue,kBlue}};
	int LineStyle[NEnergies] = {11,1};
	int MarkerStyle[NEnergies] = {21,20};

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the nuclei

		for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

			// Loop over the plots

			for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

				TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
								 nucleus[WhichNucleus]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
								 205,34,1024,768);

				PlotCanvas->SetBottomMargin(0.16);
				TLegend* leg = new TLegend(0.12,0.6,0.4,0.8);
				leg->SetNColumns(2);
				double max = -99.;

				// Loop over the energies

				for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots[WhichEnergy][WhichFSIModel] =  (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) );
						Plots[WhichEnergy][WhichFSIModel]->SetLineColor(Colors[WhichEnergy][WhichFSIModel]);
						Plots[WhichEnergy][WhichFSIModel]->SetLineWidth(LineWidth-WhichEnergy);
						CenterAxisTitle(Plots[WhichEnergy][WhichFSIModel]);

						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetTitleOffset(0.95);
						Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetNdivisions(8);
if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed") { 
	Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetRangeUser(-0.85,0.2); 
	Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetNdivisions(9); 
} else { Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->SetRangeUser(-0.81,0.07); }

						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTickSize(0.02);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);
						Plots[WhichEnergy][WhichFSIModel]->GetYaxis()->SetTitle("Weighted Events / GeV");

						if (FSIModel[WhichFSIModel] == "Data_Final") { leg->AddEntry(Plots[WhichEnergy][WhichFSIModel],"/", "ep"); }
						else { leg->AddEntry(Plots[WhichEnergy][WhichFSIModel]," "+JustE[WhichEnergy], "l"); }

						double ScalingFactor = 1. / Plots[WhichEnergy][WhichFSIModel]->Integral();

						Plots[WhichEnergy][WhichFSIModel]->Scale(ScalingFactor);

						double NBins = Plots[WhichEnergy][WhichFSIModel]->GetNbinsX(); 
				
						for (int i = 1; i <= NBins; i++) { 
					
							double content = Plots[WhichEnergy][WhichFSIModel]->GetBinContent(i);
							double error = Plots[WhichEnergy][WhichFSIModel]->GetBinError(i);
							double width = Plots[WhichEnergy][WhichFSIModel]->GetBinWidth(i);
							double newcontent = content / width;
							double newerror = error / width;				
							Plots[WhichEnergy][WhichFSIModel]->SetBinContent(i,newcontent);
							Plots[WhichEnergy][WhichFSIModel]->SetBinError(i,newerror);

						}

						double localmax = Plots[WhichEnergy][WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
						Plots[0][0]->GetYaxis()->SetRangeUser(0,1.2*max);
						PlotCanvas->Update();

						TString XLabel = Plots[WhichEnergy][WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[WhichEnergy][0]->GetXaxis()->SetTitle(XLabels[WhichPlot]);

						if (FSIModel[WhichFSIModel] == "Data_Final") { 

							Plots[WhichEnergy][WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichEnergy][WhichFSIModel]->SetMarkerColor(Colors[WhichEnergy][WhichFSIModel]);
							gStyle->SetErrorX(0);
							Plots[WhichEnergy][WhichFSIModel]->SetMarkerStyle(MarkerStyle[WhichEnergy]); 
							Plots[WhichEnergy][WhichFSIModel]->Draw("e same"); 

						}
						else { 

							Plots[WhichEnergy][WhichFSIModel]->SetLineStyle(LineStyle[WhichEnergy]);
							Plots[WhichEnergy][WhichFSIModel]->Draw("C hist same"); 
						}

					} // End of the loop over the FSI Models 


				} // End of the loop over the energies

				leg->SetBorderSize(0);
				leg->SetTextFont(FontStyle);
				leg->SetTextSize(TextSize);
				leg->Draw();

				TLatex latex;
				latex.SetTextFont(FontStyle);
				latex.SetTextSize(TextSize);
				latex.DrawLatexNDC(0.8,0.83,LabelsOfSamples[WhichNucleus]);

				TLatex latexDG;
				latexDG.SetTextFont(FontStyle);
				latexDG.SetTextSize(TextSize);
				latexDG.DrawLatexNDC(0.13,0.83,"Data/Genie");

				PlotCanvas->SaveAs("../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/FeedDown_"+nucleus[WhichNucleus]+"_" +OutputPlotNames[WhichPlot]+".pdf");

			} // End of the loop over the plots

		} // End of the loop over the nuclei

	} // End of the loop over the xB kinematics

} // End of the program
