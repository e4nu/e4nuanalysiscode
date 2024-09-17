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

//#include  "./Secondary_Code/CenterAxisTitle.cpp"
//#include "./Secondary_Code/SetOffsetAndSize.cpp"
//#include "./Secondary_Code/ToString.cpp"

void ExtDataFig2() {


	//SetOffsetAndSize();
	TGaxis::SetMaxDigits(3);

	int Ndivisions = 10;
	int LineWidth = 3;
	int FontStyle = 132;
	double TextSize = 0.05;

	std::vector<TString> nucleus; nucleus.push_back("1H"); 
	std::vector<TString> LabelsOfSamples; LabelsOfSamples.push_back("^{1}H"); 

	std::vector<TString> E; E.push_back("4_325");

	std::vector<int> Colors; Colors.push_back(kBlack); Colors.push_back(kBlack); Colors.push_back(kBlue); Colors.push_back(kOrange); Colors.push_back(60);
	std::vector<int> LineStyle; LineStyle.push_back(1); LineStyle.push_back(1); LineStyle.push_back(1); LineStyle.push_back(2); LineStyle.push_back(8);

	std::vector<TString> FSIModel; 
        FSIModel.push_back("Data");
	FSIModel.push_back("MonoChromatic");
	FSIModel.push_back("NotRadiatedV2");
	FSIModel.push_back("New"); 

	std::vector<TString> FSILabel; 
        FSILabel.push_back("Data"); 
	FSILabel.push_back("MonoChromatic");
	FSILabel.push_back("GENIE+emMCRadCorr");
	FSILabel.push_back("GENIE+radiative corrections");

        //FSILabel.push_back("GENIE + radiative correction"); 
        //FSILabel.push_back("GENIE default"); 

	std::vector<TString> NameOfPlots;

	NameOfPlots.push_back("ECalRecoPlot");

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> Plots_Clones;

	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString WhatModelsAreIncluded = "";
	for (int LoopOverFSIModels = 0 ; LoopOverFSIModels < NFSIModels ; LoopOverFSIModels ++) { WhatModelsAreIncluded += "_"+FSIModel[LoopOverFSIModels]; };

	TString RecoCalorimetry = "(e,e'p)";
	TString FSI = "FSI";


	for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {
	
      	    for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {


		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {
		    TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot],
						      nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot],
						      205,0,600,460); 
		    Plots.clear();
                        TPad* pad1 = new TPad();
                        pad1->SetFillColor(kWhite); pad1->Draw();
                        //pad1->SetTopMargin(0.1);
                        pad1->SetBottomMargin(0.19);
                        pad1->SetLeftMargin(0.15);
                        pad1->SetRightMargin(0.04);
                        pad1->cd();
			

		    TLegend* leg = new TLegend(0.2,0.66,0.46,0.88);

		    double max = -99.;

		    // Loop over the FSI Models

		    for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {


			TString PathToFiles = "/exp/genie/app/jtena/e4nuanalysiscode/data/"; // default
			if( FSIModel[WhichFSIModel] != "Data" ) PathToFiles = "";
			TFile* FileSample = TFile::Open(PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root");

			Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );
			Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
			Plots[WhichFSIModel]->SetLineStyle(LineStyle[WhichFSIModel]);
			Plots[WhichFSIModel]->SetLineWidth(LineWidth);

			Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);
			Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
			Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
			Plots[WhichFSIModel]->GetXaxis()->SetTitle("^{1}H(e,e'p) E_{cal} [GeV]");
			//Plots[WhichFSIModel]->GetXaxis()->SetTitle("E_{e'} + E_{p} [GeV]");
			Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
			Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize+0.015);
			Plots[WhichFSIModel]->GetXaxis()->SetTitleOffset(1.05);

			Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);
			Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
			Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);
			Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
			Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize+0.015);
			Plots[WhichFSIModel]->GetYaxis()->SetTitle("# events");
			Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.9);

			if (FSIModel[WhichFSIModel] == "Data") { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "lep");}
			else { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l"); }

			double ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral();
			std::cout<<"Check scaling data integral "<<Plots[0]->Integral()<<" plot integral "<<Plots[WhichFSIModel]->Integral()<<std::endl;
			Plots[WhichFSIModel]->Scale(ScalingFactor);

			double localmax = Plots[WhichFSIModel]->GetMaximum();
			if (localmax > max) { max = localmax; }
			Plots[0]->GetYaxis()->SetRangeUser(0,1.1*max);
			//if (NameOfPlots[WhichPlot] == "ECalRecoPlot") { gPad->SetLogy(); } 
			if (FSIModel[WhichFSIModel] != "Data") { Plots[WhichFSIModel]->Draw("hist e same"); }
			else { 
				gStyle->SetErrorX(0);
				Plots[WhichFSIModel]->SetMarkerStyle(20); 
				Plots[WhichFSIModel]->SetMarkerSize(1.2); 
				Plots[WhichFSIModel]->Draw("e same"); 
			     }

   		    }   // End of the loop over the FSI Models 

		leg->SetBorderSize(0);
		leg->SetTextFont(FontStyle);
		leg->SetTextSize(0.05);			
		leg->Draw();

		if (NameOfPlots[WhichPlot] == "EQERecoPlot" || NameOfPlots[WhichPlot] == "ECalRecoPlot") {

			TLine* line = new TLine(4.3256,0,4.3256,1.1*max);
			line->SetLineColor(kBlue);
			line->SetLineStyle(2);
			line->SetLineWidth(LineWidth);
			//line->Draw();
		}

		TLatex latex;
		latex.SetTextFont(FontStyle);
		latex.SetTextSize(0.07);

 
		PlotCanvas->SaveAs(nucleus[WhichNucleus]+"_" +E[WhichEnergy]+"_" +NameOfPlots[WhichPlot]+WhatModelsAreIncluded+"hist.pdf"); 
		//delete PlotCanvas;

	} // End of the loop over the plots

    } // End of the loop over the nuclei

} // End of the loop over the energies


} // End of the program
