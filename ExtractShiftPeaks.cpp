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

int LocateBin(TH1D* h, double max, double peak) {

	int NBins = h->GetXaxis()->GetNbins();
	
	for (int i = 1; i <= NBins; i++) {

		if (h->GetBinContent(i) == max && fabs(h->GetBinCenter(i)-peak)/peak < 0.05) { return i; }

	}

	return -1;
}

// ----------------------------------------------------------------------------------------------------------------

void ExtractShiftPeaks() {

	// Range between which we want to know the fraction of reconstructed events

	double range = 0.05;

	// ---------------------------------------------------------------------------------------------------------

	std::vector<TString> xBCut; 
	std::vector<TString> nucleus;
	std::vector<TString> E; 
	std::vector<double> DoubleE;
	std::vector<TString> FSIModel; 
	std::vector<TString> DirNames;
	std::vector<TString> FSILabel; 
	std::vector<TString> NameOfPlots; 
	std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames;

//	nucleus.push_back("4He");
//	nucleus.push_back("12C");
	nucleus.push_back("56Fe");

//	E.push_back("1_161"); DoubleE.push_back(1.161);
	E.push_back("2_261"); DoubleE.push_back(2.261);	
//	E.push_back("4_461"); DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_RadCorr");

	NameOfPlots.push_back("h1_Ecal"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_EQE"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");

	std::vector<TH1D*> Plots;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString WhatModelsAreIncluded = "";
	for (int LoopOverFSIModels = 0 ; LoopOverFSIModels < NFSIModels ; LoopOverFSIModels ++) { WhatModelsAreIncluded += "_"+DirNames[LoopOverFSIModels]; };

	// Loop over the xB kinematics

	for (int WhichxBCut = 0; WhichxBCut < NxBCuts; WhichxBCut ++) {

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);
						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// ---------------------------------------------------------------------------------------------------------

						// Energy Reconstruction Percentages Calculation

						if (WhichFSIModel != 0) {

							int rebin = 0;

							Plots[0]->GetXaxis()->SetRangeUser(0.9*DoubleE[WhichEnergy],1.1*DoubleE[WhichEnergy]);
//							Plots[0]->Rebin(rebin);
							double DataMax = Plots[0]->GetMaximum();
							int DataMaxBin = LocateBin(Plots[0],DataMax,DoubleE[WhichEnergy]);
							double DataMaxBinCenter = Plots[0]->GetBinCenter(DataMaxBin);

std::cout << "DataMaxBinCenter = " << DataMaxBinCenter << std::endl;

							Plots[1]->GetXaxis()->SetRangeUser(0.85*DoubleE[WhichEnergy],1.15*DoubleE[WhichEnergy]);
//							Plots[1]->Rebin(rebin);
							double GenieMax = Plots[1]->GetMaximum();
							int GenieMaxBin = LocateBin(Plots[1],GenieMax,DoubleE[WhichEnergy]);
							double GenieMaxBinCenter = Plots[1]->GetBinCenter(GenieMaxBin);

std::cout << "GenieMaxBinCenter = " << GenieMaxBinCenter << std::endl;

							double diff = (DataMaxBinCenter - GenieMaxBinCenter) * 1000.; // MeV


							std::cout << nucleus[WhichNucleus] << "  " << E[WhichEnergy] << "  " << NameOfPlots[WhichPlot] << ", offset = " << diff << " MeV" << std::endl;

						}

					} // End of the loop over the FSI Models 

					// -----------------------------------------------------------------------------------------------------------------------------------------


				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
