#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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

void Create2DPlots() {

	int TextFont = 132;
	TString version = "v3_0_6/";

	SetOffsetAndSize();

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> nucleusLabel; std::vector<TString> E; std::vector<TString> LabelE;

//	xBCut.push_back("xBCut");
	xBCut.push_back("NoxBCut");

//	nucleus.push_back("3He"); nucleusLabel.push_back("^{3}He");
//	nucleus.push_back("4He"); nucleusLabel.push_back("^{4}He");
	nucleus.push_back("12C"); nucleusLabel.push_back("^{12}C"); 
//	nucleus.push_back("56Fe"); nucleusLabel.push_back("^{56}Fe"); 

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV");
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV");
//	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");


	std::vector<int> Colors; Colors.push_back(kBlack); Colors.push_back(kBlue); Colors.push_back(kRed);

	std::vector<TString> FSIModel; 
//	FSIModel.push_back("Data"); 
//	FSIModel.push_back("hA2018_Truth_RadCorr"); 
	FSIModel.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2015_Truth_NoRadCorr");

	std::vector<TString> Cuts; 
	Cuts.push_back("Q^{2} #geq 0.5 GeV^{2}/c^{2}, |x_{B} - 1| < 0.2, W < 2 GeV/c^{2}"); Cuts.push_back("Q^{2} #geq 0.5 GeV^{2}/c^{2}, W < 2 GeV/c^{2}");

	std::vector<double> CutsX; CutsX.push_back(0.14); CutsX.push_back(0.2);
	std::vector<double> CutsY; CutsY.push_back(0.83); CutsY.push_back(0.83);

	std::vector<TString> NameOfPlots;
	NameOfPlots.push_back("Q2Vsnu"); 
//	NameOfPlots.push_back("Q2Vsnu_FirstSector"); 
//	NameOfPlots.push_back("ECalvsEQE");
//	NameOfPlots.push_back("QEECalvsEQE");
//	NameOfPlots.push_back("MECECalvsEQE");
//	NameOfPlots.push_back("RESECalvsEQE");
//	NameOfPlots.push_back("DISECalvsEQE");
//	NameOfPlots.push_back("EpVsMissMomentum"); 
//	NameOfPlots.push_back("Q2VsMissMomentum"); 
//	NameOfPlots.push_back("ElectronPhiTheta");

	TString RecoCalorimetry = "(e,e'p)";

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); gStyle->SetTitleFont(TextFont,"t"); SetOffsetAndSize();

	TF1 *f1, *f2; f1 = new TF1("f1","1.5*x",0.,5); f2 = new TF1("f2","2.25*x",0.,5);
	f1->SetLineWidth(10); f2->SetLineWidth(10);
	f1->SetLineColor(6); f2->SetLineColor(6); 

	TLatex *lat1 = new TLatex(); lat1->SetTextColor(6); lat1->SetNDC(kTRUE);
	TLatex *lat2 = new TLatex(); lat2->SetTextColor(6); lat2->SetNDC(kTRUE);

	std::vector<TCanvas*> PlotCanvas;

	std::vector<TH2D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

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

					TLegend* leg = new TLegend(0.15,0.45,0.35,0.75);

					Plots.clear();

					// Loop over the plots

					for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

						PlotCanvas.push_back(new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
										 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
										 205,34,1024,768));

						Plots.push_back( (TH2D*)(FileSample->Get(NameOfPlots[WhichPlot])) );
						
						Plots[WhichPlot]->SetTitle(nucleusLabel[WhichNucleus]+RecoCalorimetry+LabelE[WhichEnergy]+" ("+FSIModel[WhichFSIModel]+")");
						Plots[WhichPlot]->GetXaxis()->SetLabelFont(TextFont);
						Plots[WhichPlot]->GetXaxis()->SetTitleFont(TextFont);
						Plots[WhichPlot]->GetYaxis()->SetLabelFont(TextFont);
						Plots[WhichPlot]->GetYaxis()->SetTitleFont(TextFont);
						Plots[WhichPlot]->GetZaxis()->SetLabelFont(TextFont);
						Plots[WhichPlot]->GetZaxis()->SetTitleFont(TextFont);
						Plots[WhichPlot]->GetZaxis()->SetLabelSize(0.03);
						Plots[WhichPlot]->Scale(1./Plots[WhichPlot]->GetMaximum());
						Plots[WhichPlot]->GetZaxis()->SetRangeUser(0,Plots[WhichPlot]->GetMaximum());
//Plots[WhichPlot]->GetXaxis()->SetRangeUser(0.5,2.5);
//Plots[WhichPlot]->GetYaxis()->SetRangeUser(0.5,2.5);
						CenterAxisTitle(Plots[WhichPlot]);
//						Plots[WhichPlot]->GetZaxis()->SetLabelSize(0);
//PlotCanvas[PlotCanvas.size()-1]->SetLogz();

if (NameOfPlots[WhichPlot] == "Q2Vsnu" || NameOfPlots[WhichPlot] == "Q2Vsnu_FirstSector") {

	if (E[WhichEnergy] == "1_161") { Plots[WhichPlot]->GetXaxis()->SetRangeUser(0,1.); Plots[WhichPlot]->GetYaxis()->SetRangeUser(0,1.);}
	if (E[WhichEnergy] == "2_261") { Plots[WhichPlot]->GetXaxis()->SetRangeUser(0.1,2.); Plots[WhichPlot]->GetYaxis()->SetRangeUser(0.1,2.); }
//	if (E[WhichEnergy] == "2_261") { Plots[WhichPlot]->GetXaxis()->SetRangeUser(0,2.5); Plots[WhichPlot]->GetYaxis()->SetRangeUser(0,2.5);}

}
						Plots[WhichPlot]->Draw("colz");

						if (NameOfPlots[WhichPlot] == "Q2Vsnu") { 

							f1->Draw("same"); f2->Draw("same"); 
							lat1->DrawLatex(0.2,0.8,"x_{B} = 1.2"); lat2->DrawLatex(0.6,0.4,"x_{B} = 0.8"); 
						}

						PlotCanvas[PlotCanvas.size()-1]->SaveAs("myPlots/eps/"+xBCut[WhichxBCut]+"/"+version
											+nucleus[WhichNucleus]+"_" +E[WhichEnergy]+"_" +NameOfPlots[WhichPlot]+"_"+FSIModel[WhichFSIModel]+".eps"); 

						PlotCanvas[PlotCanvas.size()-1]->SaveAs("myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version
											+nucleus[WhichNucleus]+"_" +E[WhichEnergy]+"_" +NameOfPlots[WhichPlot]+"_"+FSIModel[WhichFSIModel]+".pdf");

					} // End of the loop over the plots			

				} // End of the loop over the FSI Models 

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
