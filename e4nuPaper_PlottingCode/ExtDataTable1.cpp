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

void ExtDataTable1() {

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
//	E.push_back("2_261"); DoubleE.push_back(2.261);	
	E.push_back("4_461"); DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_RadCorr");

	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
	NameOfPlots.push_back("eRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E_{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");

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

					TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 205,34,1024,768);

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);
						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						// ---------------------------------------------------------------------------------------------------------

						// Energy Reconstruction Percentages Calculation

						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos) {

							// Get percentages of events within X % of peak energy (5% for ECal/ 10% for EQE)

							// QE Energy Reconstruction

							if (string(OutputPlotNames[WhichPlot]).find("eRecoEnergy_slice") != std::string::npos) {

								double LowE = (1.-range)*DoubleE[WhichEnergy], HighE = (1.+range)*DoubleE[WhichEnergy];
								int LowEBin = Plots[WhichFSIModel]->FindBin(LowE);
								int HighEBin = Plots[WhichFSIModel]->FindBin(HighE);

								int PercLowPmiss = Plots[WhichFSIModel]->Integral(1,LowEBin) / Plots[WhichFSIModel]->Integral() * 100.;
								int PercMidPmiss = Plots[WhichFSIModel]->Integral(LowEBin+1,HighEBin) / Plots[WhichFSIModel]->Integral() * 100.;
								int PercHighPmiss = 100 - PercLowPmiss - PercMidPmiss;

							}

							// Calorimetric Energy Reconstruction

							else {

								double LowE = (1.-range)*DoubleE[WhichEnergy];
								int LowEBin = Plots[WhichFSIModel]->FindBin(LowE);

								int PercLowPmiss = Plots[WhichFSIModel]->Integral(1,LowEBin) / Plots[WhichFSIModel]->Integral() * 100.;
								int PercHighPmiss = 100 - PercLowPmiss;

							}

							// ---------------------------------------------------------------------------------------------------------

							double MaxContent = Plots[WhichFSIModel]->GetMaximum();
							int BinWithMaxContent = -99;

							int Nbins = Plots[WhichFSIModel]->GetXaxis()->GetNbins();

							for (int i = 0; i < Nbins; i++) {
								if (Plots[WhichFSIModel]->GetBinContent(i+1) ==  MaxContent) { BinWithMaxContent = i+1; }
							}

							TAxis *xaxis = Plots[WhichFSIModel]->GetXaxis();
							double binCenter = xaxis->GetBinCenter(BinWithMaxContent);

							// Percentages with respect to the peak energy

							double MinE = (1.-range)*binCenter;
							double MaxE = (1.+range)*binCenter;

							int MinBin = Plots[WhichFSIModel]->FindBin(MinE);
							int MaxBin = Plots[WhichFSIModel]->FindBin(MaxE);
							int percentage = Plots[WhichFSIModel]->Integral(MinBin,MaxBin) / Plots[WhichFSIModel]->Integral() * 100.;

							cout << endl << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": fraction within " << int(range*100.) 
							     << "% = " << percentage << endl;

						}

					} // End of the loop over the FSI Models 

					// -----------------------------------------------------------------------------------------------------------------------------------------

					delete PlotCanvas;

					cout << "---------------------------------------------------------------" << endl;

					// -----------------------------------------------------------------------------------------------------------------------------------------


				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
