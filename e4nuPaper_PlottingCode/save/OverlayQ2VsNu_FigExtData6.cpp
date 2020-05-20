#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "/home/afroditi/Dropbox/PhD/Secondary_Code/CenterAxisTitle.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/SetOffsetAndSize.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"

void OverlayQ2VsNu_FigExtData6() {

	int Ndivisions = 4;
	int FontStyle = 132;
	double TextSize = 0.08;
	
	TString version = "v3_0_6/";

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(TextSize,"t"); gStyle->SetTitleFont(FontStyle,"t"); SetOffsetAndSize();

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> LabelsOfSamples; std::vector<TString> E;  std::vector<TString> JustNucleus;
	std::vector<TString> LabelE; std::vector<TString> FSIModel; std::vector<TString> DirNames;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots; std::vector<TString> XLabelOfPlots; std::vector<TString> YLabelOfPlots;  std::vector<TString> OutputPlotNames;

//	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He");  JustNucleus.push_back("He");
	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C"); JustNucleus.push_back("C");
//	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");  JustNucleus.push_back("Fe");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV");
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.257 GeV");
//	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("GENIE");  DirNames.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("GENIE");  DirNames.push_back("hA2018_Truth_NoRadCorr");

//	NameOfPlots.push_back("h2_Q2_nu_weight"); 
	NameOfPlots.push_back("h2_Q2_nu_weight_FirstSector"); 
	XLabelOfPlots.push_back("Energy Transfer [GeV]"); 
	YLabelOfPlots.push_back("Q^{2} [GeV^{2}/c^{2}]"); 
	OutputPlotNames.push_back("Q2VsNu2D_FirstSector");

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

	TString FSI = "FSI";

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
										 205,34,2024,768);

					// ---------------------------------------------------------------------------------------------------------------------------

					// Dimensions of TPads

					double XMinPadOne = 0., XMaxPadOne = 0.5, YMinPadOne = 0., YMaxPadOne = 1.;
					double XMinPadTwo = 0.5, XMaxPadTwo = 1., YMinPadTwo = 0., YMaxPadTwo = YMaxPadOne;
					double XMinPadThree = 0.5, XMaxPadThree = 0.53, YMinPadThree = 0.07, YMaxPadThree = 0.17;
					double XMinPadFour = 0.2, XMaxPadFour = 0.8, YMinPadFour = 0.91, YMaxPadFour = 1.;

					// ---------------------------------------------------------------------------------------------------------------------------

					TPad* pad1 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadOne,YMinPadOne,XMaxPadOne,YMaxPadOne, 21); 
					pad1->SetFillColor(kWhite); pad1->Draw();
					TPad* pad2 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadTwo,YMinPadTwo,XMaxPadTwo,YMaxPadTwo,22); 
					pad2->SetFillColor(kWhite); pad2->Draw(); 
					TPad* pad3 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadThree,YMinPadThree,XMaxPadThree,YMaxPadThree,23); 
					pad3->SetFillColor(kWhite); pad3->Draw();
					TPad* pad4 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadFour,YMinPadFour,XMaxPadFour,YMaxPadFour,24); 
					pad4->SetFillColor(kWhite); pad4->Draw();
					pad1->SetBottomMargin(0.18);
					pad2->SetBottomMargin(0.18);

					pad4->cd();
					TLatex *title = new TLatex(); 
					title->SetNDC();
					title->SetTextFont(FontStyle); 
					title->SetTextColor(kBlack); 
					title->SetTextSize(0.8);
					TString myTitle = LabelsOfSamples[WhichNucleus] + " " +LabelE[WhichEnergy];
					title->DrawLatex(0.25,0.3,myTitle);

					// ---------------------------------------------------------------------------------------------------------------------------

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						if (FSILabel[WhichFSIModel] == "Data") 
							{ pad1->cd(); gStyle->SetTitleSize(TextSize,"t"); pad1->SetRightMargin(0.); pad1->SetLeftMargin(0.15); }
						else { pad2->cd(); pad2->SetLeftMargin(0.0); pad2->SetRightMargin(0.15); }

						TString PathToFiles = "../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						TH2D* Plots =  (TH2D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) ;

						CenterAxisTitle(Plots);
						Plots->SetTitleSize(TextSize,"t");
						if (FSILabel[WhichFSIModel] == "Data") { gStyle->SetTitleX(.54); }
						else { gStyle->SetTitleX(.47); }

						Plots->GetXaxis()->SetLabelFont(FontStyle);
						Plots->GetXaxis()->SetTitleFont(FontStyle);
						Plots->GetXaxis()->SetLabelSize(TextSize);
						Plots->GetXaxis()->SetTitleSize(TextSize);
						Plots->GetXaxis()->SetTitleOffset(1.);
						Plots->GetXaxis()->SetTitle(XLabelOfPlots[WhichPlot]);


						Plots->GetYaxis()->SetLabelFont(FontStyle);
						Plots->GetYaxis()->SetTitleFont(FontStyle);
						Plots->GetYaxis()->SetLabelSize(TextSize);
						Plots->GetYaxis()->SetTitleSize(TextSize);
						Plots->GetYaxis()->SetTitleOffset(0.8);
						Plots->GetYaxis()->SetTitle(YLabelOfPlots[WhichPlot]);

						// --------------------------------------------------------------------------------------------------------------------------

						// Rebinning & Ranges

						for (int i = 0; i < 0; i++) { Plots->Rebin2D(); }

						Plots->GetZaxis()->SetRangeUser(1.,Plots->GetMaximum());

						double XMin =-99.,XMax =-99.;
						double YMin =-99.,YMax =-99.;
						XMin = 0.; XMax = 1.8; Plots->GetXaxis()->SetRangeUser(XMin,XMax); 
						YMin = 0.; YMax = 2.;	Plots->GetYaxis()->SetRangeUser(YMin,YMax);

						// -----------------------------------------------------------------------------------------------------------------------------

//						Plots->Draw("coltz");
						Plots->Draw("colt");
						PlotCanvas->Update();

						// ----------------------------------------------------------------------------------------------------------------

						// TLines & TLatex

						TLatex *sample = new TLatex(); 
						sample->SetTextFont(FontStyle); 
						sample->SetTextColor(kBlack); 
						sample->SetTextSize(TextSize);
						if (FSILabel[WhichFSIModel] == "Data") { sample->DrawTextNDC(0.2,0.82,FSILabel[WhichFSIModel]); }
						else { sample->DrawTextNDC(0.05,0.82,FSILabel[WhichFSIModel]); } 

						TF1 *f1; f1 = new TF1("f1","1.876*x",0.,1.8);
						f1->SetLineWidth(2);
						f1->SetLineColor(kBlack);

						TLatex *lat1 = new TLatex(); lat1->SetTextColor(1);

						f1->Draw("same");
						lat1->SetTextFont(132);
						lat1->SetTextSize(TextSize);
						lat1->DrawLatex(1.22,1.85,"x_{B} = 1"); 

						if (FSILabel[WhichFSIModel] == "Genie" ) { Plots->GetYaxis()->SetTitle(); Plots->GetYaxis()->SetLabelSize(0.); }

						Plots->GetXaxis()->SetNdivisions(Ndivisions);
						Plots->GetYaxis()->SetNdivisions(Ndivisions);
						
						// --------------------------------------------------------------------------------------------------

//						TPaletteAxis *palette = (TPaletteAxis*)Plots->GetListOfFunctions()->FindObject("palette");
//						palette->SetX1NDC(0.87);
//						palette->SetX2NDC(0.9);
//						PlotCanvas->Modified();

					} // End of the loop over the FSI Models 

					PlotCanvas->SaveAs("../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]
							+"/FigExtData6_"+nucleus[WhichNucleus]+"_" 
							+E[WhichEnergy]+"_" +OutputPlotNames[WhichPlot]+".pdf");

					//delete PlotCanvas;

				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
