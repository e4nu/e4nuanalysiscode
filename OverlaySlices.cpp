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

#include <iostream>
#include <vector>

using namespace std;

#include  "./Secondary_Code/CenterAxisTitle.cpp"
#include "./Secondary_Code/SetOffsetAndSize.cpp"
#include "./Secondary_Code/ToString.cpp"

void OverlaySlices() {

	int TextFont = 132;

	SetOffsetAndSize();

	std::vector<TString> xBCut; /*xBCut.push_back("xBCut");*/ xBCut.push_back("NoxBCut");

	std::vector<TString> nucleus; nucleus.push_back("12C"); /*nucleus.push_back("56Fe"); nucleus.push_back("4He");*/
	std::vector<TString> LabelsOfSamples; LabelsOfSamples.push_back("^{12}C"); /*LabelsOfSamples.push_back("^{56}Fe"); LabelsOfSamples.push_back("^{4}He");*/

	TString MainPlot = "RecoEnergy_slice_"; 
	std::vector<TString> NameOfPlots; NameOfPlots.push_back("ep"); /*NameOfPlots.push_back("e");*/

	std::vector<TString> E; /*E.push_back("1_161");*/ E.push_back("2_261");/* E.push_back("4_461");*/
	std::vector<TString> DoubleE; /*DoubleE.push_back("1_161");*/ DoubleE.push_back("2.261"); /*DoubleE.push_back("4_461");*/
	std::vector<double> doubleE; /*doubleE.push_back(1_161);*/ doubleE.push_back(2.261); /*doubleE.push_back(4.461);*/
	std::vector<TString> LabelE; /*LabelE.push_back(" @ E = 1.161 GeV");*/ LabelE.push_back(" @ E = 2.261 GeV"); /*LabelE.push_back(" @ E = 4.461 GeV");*/

	std::vector<double> EMin; /*EMin.push_back(0.5);*/ EMin.push_back(1.); /*EMin.push_back(1.5);*/
	std::vector<double> EMax; /*EMin.push_back(0.5);*/ EMax.push_back(2.6); /*EMin.push_back(5.3);*/

	std::vector<TString> Legends; 
	Legends.push_back("P^{miss}_{#perp}  < 0.2 GeV/c"); Legends.push_back("0.2 GeV/c < P^{miss}_{#perp}  < 0.4 GeV/c");  Legends.push_back("P^{miss}_{#perp}  > 0.4 GeV/c");

	std::vector<int> Colors; Colors.push_back(4); Colors.push_back(2); Colors.push_back(416);

	std::vector<TString> FSIModel; FSIModel.push_back("Data"); FSIModel.push_back("hA2018_RadCorr"); /*FSIModel.push_back("hA2018_NoRadCorr");*/

	std::vector<TString> Cuts; 
	Cuts.push_back("Q^{2} #geq 0.5 GeV^{2}/c^{2}, |x_{B} - 1| < 0.2, W < 2 GeV/c^{2}"); Cuts.push_back("Q^{2} #geq 0.5 GeV^{2}/c^{2}, W < 2 GeV/c^{2}");

	std::vector<double> CutsX; CutsX.push_back(0.12); CutsX.push_back(0.23);
	std::vector<double> CutsY; CutsY.push_back(0.83); CutsY.push_back(0.83);

	std::vector<TCanvas*> PlotCanvas;

	std::vector<TH1D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NSlices = Legends.size();
	int NPlots = NameOfPlots.size();

	TString RecoCalorimetry = "(e,e'p)";

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// Loop over the FSI Models

				for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {


					TString PathToFiles = "myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";

					TFile* FileSample = TFile::Open(PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root");

					// Loop over the plots

					for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {


						PlotCanvas.push_back(new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
										 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
										 205,34,1024,768));

						Plots.clear();

						TLegend* leg = new TLegend(0.12,0.58,0.32,0.88);

double max = -99.;

						// Loop over the slices

						for (int WhichSlice = 0; WhichSlice < NSlices; WhichSlice ++) {

							Plots.push_back( (TH1D*)(FileSample->Get(NameOfPlots[WhichPlot]+MainPlot+ToString(WhichSlice+1))) );

							Plots[WhichSlice]->GetXaxis()->SetLabelFont(TextFont);
							Plots[WhichSlice]->GetXaxis()->SetNdivisions(5);

							Plots[WhichSlice]->GetXaxis()->SetTitleFont(TextFont);
							Plots[WhichSlice]->GetYaxis()->SetLabelFont(TextFont);
							Plots[WhichSlice]->GetYaxis()->SetTitleFont(TextFont);
							Plots[WhichSlice]->GetYaxis()->SetLabelSize(0);
							Plots[WhichSlice]->GetYaxis()->SetNdivisions(10);
							Plots[WhichSlice]->GetYaxis()->SetTickSize(0.);


							CenterAxisTitle(Plots[WhichSlice]);
							Plots[WhichSlice]->GetXaxis()->SetRangeUser(EMin[WhichEnergy],EMax[WhichEnergy]);
							Plots[WhichSlice]->SetLineColor(Colors[WhichSlice]);
							leg->AddEntry(Plots[WhichSlice],Legends[WhichSlice]);
double localmax = Plots[WhichSlice]->GetMaximum();
if (max < localmax) { max = localmax; }

//							Plots[WhichSlice]->GetYaxis()->SetRangeUser(0,1.3*Plots[WhichSlice]->GetMaximum());
Plots[0]->GetYaxis()->SetRangeUser(0,1.1*max);
							Plots[WhichSlice]->SetLineWidth(3);
							Plots[WhichSlice]->Draw("hist same");

						} // End of the loop over the slices

						leg->SetBorderSize(0);
						leg->SetTextFont(TextFont);
						leg->SetTextSize(0.05);			
						leg->Draw();

						TLatex latex;
						latex.SetTextFont(TextFont);
						latex.SetTextSize(0.07);
/*						latex.DrawLatexNDC(0.22,0.93,LabelsOfSamples[WhichNucleus]+RecoCalorimetry+LabelE[WhichEnergy]);*/

						TLatex latexCuts;
						latexCuts.SetTextFont(TextFont);
						latexCuts.SetTextSize(0.06);
/*						latexCuts.DrawLatexNDC(CutsX[WhichxBCut],CutsY[WhichxBCut],Cuts[WhichxBCut]);*/

						TLatex latexFSIModel;
						latexFSIModel.SetTextFont(TextFont);
						latexFSIModel.SetTextSize(0.1);
						latexFSIModel.SetTextColorAlpha(kRed+2, 0.35);
/*						latexFSIModel.DrawLatexNDC(0.15,0.45,"#bf{"+FSIModel[WhichFSIModel]+"}");*/

						TLine* line = new TLine(doubleE[WhichEnergy],0.,doubleE[WhichEnergy],Plots[Plots.size()-3]->GetMaximum());
						line->SetLineColor(kMagenta);
						line->SetLineStyle(5);
						line->SetLineWidth(5);
						line->Draw();

						PlotCanvas[PlotCanvas.size()-1]->SaveAs("myPlots/eps/"+nucleus[WhichNucleus]+"_" +E[WhichEnergy]+"_" +NameOfPlots[WhichPlot]+"_"+FSIModel[WhichFSIModel]+".eps"); 

						PlotCanvas[PlotCanvas.size()-1]->SaveAs("myPlots/pdf/"+nucleus[WhichNucleus]+"_" +E[WhichEnergy]+"_" +NameOfPlots[WhichPlot]+"_"+FSIModel[WhichFSIModel]+".pdf"); 

					} // End of the loop over the plots			

				} // End of the loop over the FSI Models 

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
