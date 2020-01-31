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

#include  "./Secondary_Code/CenterAxisTitle.cpp"
#include "./Secondary_Code/SetOffsetAndSize.cpp"
#include "./Secondary_Code/ToString.cpp"

void Create2DPlots() {

	int Ndivisions = 4;
	int FontStyle = 132;
	double TextSize = 0.08;
	
	TString version = "v3_0_6/";

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(TextSize,"t"); gStyle->SetTitleFont(FontStyle,"t"); SetOffsetAndSize();

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> LabelsOfSamples; std::vector<TString> E;  std::vector<TString> JustNucleus;
	std::vector<TString> LabelE; std::vector<TString> FSIModel; std::vector<TString> DirNames;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots; std::vector<TString> XLabelOfPlots; std::vector<TString> YLabelOfPlots;  std::vector<TString> OutputPlotNames;

////	nucleus.push_back("3He"); LabelsOfSamples.push_back("^{3}He");
//	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He");  JustNucleus.push_back("He");
//	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C"); JustNucleus.push_back("C");
	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");  JustNucleus.push_back("Fe");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV");
//	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV");
	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("hA2018_Final_NoRadCorr"); FSILabel.push_back("GENIE");  DirNames.push_back("hA2018_Truth_NoRadCorr");
//	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("GENIE");  DirNames.push_back("hA2018_Truth_NoRadCorr");

//	FSIModel.push_back("hA2018_Truth_NoRadCorr"); FSILabel.push_back("GENIE (Truth)");  DirNames.push_back("hA2018_Truth_NoRadCorr");
//	FSIModel.push_back("hN2018_Final_NoRadCorr"); FSILabel.push_back("GENIE hN2018");  DirNames.push_back("hN2018_Truth_NoRadCorr");

//	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Adi"); FSILabel.push_back("NoRad");  DirNames.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM_Adi"); FSILabel.push_back("Rad");  DirNames.push_back("hA2018_Truth_RadCorr");


//	NameOfPlots.push_back("h2_Ecal_Eqe"); XLabelOfPlots.push_back("E^{QE} (GeV)"); YLabelOfPlots.push_back("E^{cal} (GeV)"); OutputPlotNames.push_back("ECalVsEQE2D");
//	NameOfPlots.push_back("h2_Q2_nu_weight"); XLabelOfPlots.push_back("Energy Transfer (GeV)"); YLabelOfPlots.push_back("Q^{2} (GeV^{2}/c^{2})"); OutputPlotNames.push_back("Q2VsNu2D");

//	NameOfPlots.push_back("h2_Q2_nu_weight_FirstSector"); XLabelOfPlots.push_back("Energy Transfer [GeV]"); YLabelOfPlots.push_back("Q^{2} [GeV^{2}/c^{2}]"); OutputPlotNames.push_back("Q2VsNu2D_FirstSector");

//	NameOfPlots.push_back("EePrimeVsEgamma"); XLabelOfPlots.push_back("E_{#gamma} [GeV]"); YLabelOfPlots.push_back("E_{e'} [GeV]"); OutputPlotNames.push_back("EePrimeVsEgamma");

	NameOfPlots.push_back("RadCosThetaGammaEgamma"); XLabelOfPlots.push_back("cos(#theta_{#gamma})"); YLabelOfPlots.push_back("E_{#gamma} [GeV]"); OutputPlotNames.push_back("CosThetaGammaEgamma");

	NameOfPlots.push_back("RadCosDeltaThetaGammaEgamma"); XLabelOfPlots.push_back("cos(#Delta#theta_{#gamma,e'})"); YLabelOfPlots.push_back("E_{#gamma} [GeV]"); OutputPlotNames.push_back("CosDeltaThetaGammaEgamma");

//	NameOfPlots.push_back("h2_Etot_pperp"); XLabelOfPlots.push_back("P_{miss}^{#perp} [GeV/c]"); YLabelOfPlots.push_back("E^{cal} (GeV)"); OutputPlotNames.push_back("ECalVsPmiss2D");

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
//										 205,34,1024,768);
//										 205,34,2048,1656);
										 205,34,2024,768);

					// ---------------------------------------------------------------------------------------------------------------------------

					// Dimensions of TPads (pad2 will be deleted at the very end for the Ereco plots)

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
//					pad2->SetTopMargin(0.21);
//					pad1->cd();

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
//						else { pad2->cd(); pad2->SetLeftMargin(0.014); pad2->SetRightMargin(0.15); }
						else { pad2->cd(); pad2->SetLeftMargin(0.0); pad2->SetRightMargin(0.15); }

//						TCanvas* PlotCanvas = new TCanvas(FSIModel[WhichFSIModel]+"_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
//										 FSIModel[WhichFSIModel]+"_"+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
//										 205,34,1024,768);

						TString PathToFiles = "../../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						TH2D* Plots =  (TH2D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) ;

						CenterAxisTitle(Plots);
//						Plots->SetTitle(FSILabel[WhichFSIModel]+" " + LabelsOfSamples[WhichNucleus] + " " +LabelE[WhichEnergy]);
						Plots->SetTitleSize(TextSize,"t");
						if (FSILabel[WhichFSIModel] == "Data") { gStyle->SetTitleX(.54); }
						else { gStyle->SetTitleX(.47); }

						Plots->GetXaxis()->SetLabelFont(FontStyle);
						Plots->GetXaxis()->SetTitleFont(FontStyle);
						Plots->GetXaxis()->SetLabelSize(TextSize);
						Plots->GetXaxis()->SetTitleSize(TextSize);
						Plots->GetXaxis()->SetTitleOffset(1.);
						if ( NameOfPlots[WhichPlot] == "h2_Ecal_Eqe" ) 
							{ Plots->GetXaxis()->SetTitle(JustNucleus[WhichNucleus]+"(e,e'p)_{1p0#pi }"+XLabelOfPlots[WhichPlot]); }
						else { Plots->GetXaxis()->SetTitle(XLabelOfPlots[WhichPlot]); }


						Plots->GetYaxis()->SetLabelFont(FontStyle);
						Plots->GetYaxis()->SetTitleFont(FontStyle);
						Plots->GetYaxis()->SetLabelSize(TextSize);
						Plots->GetYaxis()->SetTitleSize(TextSize);
						Plots->GetYaxis()->SetTitleOffset(0.8);
						if ( NameOfPlots[WhichPlot] == "h2_Ecal_Eqe" || NameOfPlots[WhichPlot] == "h2_Etot_pperp" ) 
							{ Plots->GetYaxis()->SetTitle(JustNucleus[WhichNucleus]+"(e,e'p)_{1p0#pi} "+YLabelOfPlots[WhichPlot]); }
						else { Plots->GetYaxis()->SetTitle(YLabelOfPlots[WhichPlot]); }

						// --------------------------------------------------------------------------------------------------------------------------

						// Rebinning & Ranges

						Plots->GetZaxis()->SetRangeUser(1.,Plots->GetMaximum());
						double XMin =-99.,XMax =-99.;
						double YMin =-99.,YMax =-99.;

						if ( NameOfPlots[WhichPlot] == "h2_Ecal_Eqe" ) {
//						PlotCanvas->SetLogz();
						pad1->SetLogz();
						pad2->SetLogz();

						if (E[WhichEnergy] == "1_161") { 
							for (int i = 0; i < 1; i++) { Plots->Rebin2D(); }
							XMin = 0.4; XMax = 1.8; Plots->GetXaxis()->SetRangeUser(XMin,XMax); 
							Plots->GetYaxis()->SetRangeUser(XMin,XMax); 
						}

						if (E[WhichEnergy] == "2_261") { 
							for (int i = 0; i < 2; i++) { Plots->Rebin2D(); } 
							XMin = 0.5; XMax = 3.; Plots->GetXaxis()->SetRangeUser(XMin,XMax); 
							Plots->GetYaxis()->SetRangeUser(XMin,XMax); 
						}
					
						if (E[WhichEnergy] == "4_461") { 
							for (int i = 0; i < 4; i++) { Plots->Rebin2D(); }
							XMin = 1.5; XMax = 5.; Plots->GetXaxis()->SetRangeUser(XMin,XMax); 
							Plots->GetYaxis()->SetRangeUser(XMin,XMax); }
						}


						if ( string(NameOfPlots[WhichPlot]).find("h2_Q2_nu_weight") != std::string::npos ) {

							if (E[WhichEnergy] == "1_161") { 
								for (int i = 0; i < 1; i++) { Plots->Rebin2D(); } /*PlotCanvas->SetLogz();*/ XMin = 0.; XMax = 0.8; 
								Plots->GetXaxis()->SetRangeUser(XMin,XMax); Plots->GetYaxis()->SetRangeUser(XMin,XMax); 
							}

							if (E[WhichEnergy] == "2_261") { 
								for (int i = 0; i < 2; i++) { Plots->Rebin2D(); } XMin = 0.; XMax = 2.; 
								Plots->GetXaxis()->SetRangeUser(XMin,XMax); Plots->GetYaxis()->SetRangeUser(XMin,XMax); 
							}

							if (E[WhichEnergy] == "4_461") { 
								for (int i = 0; i < 3; i++) { Plots->Rebin2D(); }
								XMin = 0.; XMax = 6.; Plots->GetXaxis()->SetRangeUser(XMin,XMax); 
								Plots->GetYaxis()->SetRangeUser(XMin,XMax); }

						}

						if ( NameOfPlots[WhichPlot] == "h2_Q2_nu_weight_FirstSector" ) {

							if (E[WhichEnergy] == "2_261") { 
								for (int i = 0; i < 0; i++) { Plots->Rebin2D(); } 
								XMin = 0.; XMax = 1.8; Plots->GetXaxis()->SetRangeUser(XMin,XMax); 
								YMin = 0.; YMax = 2.;	Plots->GetYaxis()->SetRangeUser(YMin,YMax); 
							}

						}

						if ( NameOfPlots[WhichPlot] == "h2_Etot_pperp" ) {

//							pad1->SetLogz();
//							pad2->SetLogz();

							if (E[WhichEnergy] == "2_261") { 
								for (int i = 0; i < 2; i++) { Plots->Rebin2D(); } 
								XMin = 0.; XMax = 1.; Plots->GetXaxis()->SetRangeUser(XMin,XMax); 
								YMin = 0.5; YMax = 2.8;	Plots->GetYaxis()->SetRangeUser(YMin,YMax); 
							}

						}

						// -----------------------------------------------------------------------------------------------------------------------------

						double ScalingFactor = TMath::Power(10.,6.) / Plots->GetMaximum();
						Plots->Scale(ScalingFactor);
//						Plots->GetZaxis()->SetRangeUser(1.,Plots->GetMaximum());
						if ( NameOfPlots[WhichPlot] == "h2_Ecal_Eqe" ) { Plots->GetZaxis()->SetRangeUser(Plots->GetMaximum()*1E-2,Plots->GetMaximum()); }
						else { Plots->GetZaxis()->SetRangeUser(1.,Plots->GetMaximum()); }
						Plots->GetZaxis()->SetLabelSize(TextSize);
						Plots->GetZaxis()->SetLabelFont(FontStyle);
						Plots->GetZaxis()->SetTitle("Weighted Events");
						Plots->GetZaxis()->CenterTitle();
						Plots->GetZaxis()->SetTitleFont(FontStyle);
						Plots->GetZaxis()->SetTitleSize(TextSize);

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

						if ( NameOfPlots[WhichPlot] == "h2_Q2_nu_weight" ) {

							TF1 *f1, *f2; f1 = new TF1("f1","1.5*x",0.,5); f2 = new TF1("f2","2.25*x",0.,5);
							f1->SetLineWidth(10); f2->SetLineWidth(10);
							f1->SetLineColor(6); f2->SetLineColor(6); 

							TLatex *lat1 = new TLatex(); lat1->SetTextColor(6); lat1->SetNDC(kTRUE);
							TLatex *lat2 = new TLatex(); lat2->SetTextColor(6); lat2->SetNDC(kTRUE);

							f1->Draw("same"); f2->Draw("same"); 
							lat1->DrawLatex(0.2,0.8,"x_{B} = 1.2"); lat2->DrawLatex(0.6,0.4,"x_{B} = 0.8"); 
						}


						if ( string(NameOfPlots[WhichPlot]).find("h2_Ecal_Eqe") != std::string::npos ) {

							TF1 *f1; f1 = new TF1("f1","x",0.,6.);
							f1->SetLineWidth(2);
							f1->SetLineColor(kBlack); 
							TLatex *lat1 = new TLatex(); lat1->SetTextColor(6); lat1->SetNDC(kTRUE);
							f1->Draw("same");

							if (FSILabel[WhichFSIModel] == "Genie" ) { Plots->GetYaxis()->SetTitle(); Plots->GetYaxis()->SetLabelSize(0.); }

							Plots->GetXaxis()->SetNdivisions(Ndivisions);
							Plots->GetYaxis()->SetNdivisions(Ndivisions);
						}

						if ( NameOfPlots[WhichPlot] == "h2_Q2_nu_weight_FirstSector" ) {

							TF1 *f1; f1 = new TF1("f1","1.876*x",0.,1.8);
							f1->SetLineWidth(2);
							f1->SetLineColor(kBlack);

							TLatex *lat1 = new TLatex(); lat1->SetTextColor(1);

							f1->Draw("same");
							lat1->SetTextFont(132);
							lat1->SetTextSize(TextSize);
							lat1->DrawLatex(1.22,1.85,"x_{B} = 1"); 

//							Plots->GetXaxis()->SetTitle();
							if (FSILabel[WhichFSIModel] == "Genie" ) { Plots->GetYaxis()->SetTitle(); Plots->GetYaxis()->SetLabelSize(0.); }

							Plots->GetXaxis()->SetNdivisions(Ndivisions);
							Plots->GetYaxis()->SetNdivisions(Ndivisions);

						}

						if ( NameOfPlots[WhichPlot] == "h2_Etot_pperp" ) {

//							Plots->GetXaxis()->SetTitle();
							if (FSILabel[WhichFSIModel] == "Genie" ) { Plots->GetYaxis()->SetTitle(); Plots->GetYaxis()->SetLabelSize(0.); }

							Plots->GetXaxis()->SetNdivisions(5);
							Plots->GetYaxis()->SetNdivisions(Ndivisions);

						}

						// --------------------------------------------------------------------------------------------------

//						TPaletteAxis *palette = (TPaletteAxis*)Plots->GetListOfFunctions()->FindObject("palette");
//						palette->SetX1NDC(0.87);
//						palette->SetX2NDC(0.9);
//						PlotCanvas->Modified();

						//delete PlotCanvas;

					} // End of the loop over the FSI Models 

//					PlotCanvas->SaveAs("../../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]+"/"+nucleus[WhichNucleus]+"_" 
//							+E[WhichEnergy]+"_" +OutputPlotNames[WhichPlot]+".pdf");

						//delete PlotCanvas;

				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
