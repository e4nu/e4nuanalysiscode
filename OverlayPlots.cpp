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

void OverlayPlots() {

	// ------------------------------------------------------------------------

	SetOffsetAndSize();
	TGaxis::SetMaxDigits(3);
//	TGaxis::SetExponentOffset(-0.1, 1., "y");

	int Ndivisions = 4;
	int LineWidth = 3;
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

//	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He"); JustNucleus.push_back("He");
//	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C"); JustNucleus.push_back("C");
	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");  JustNucleus.push_back("Fe");

//	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); DoubleE.push_back(1.161);
	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); DoubleE.push_back(2.261);	
//	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");  DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");
 
//	Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kBlue); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);
	Colors.push_back(kBlack); Colors.push_back(kBlue); Colors.push_back(kBlack); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);

//	Style.push_back(9); Style.push_back(3); Style.push_back(7); Style.push_back(5);
//	Style.push_back(9); Style.push_back(9); Style.push_back(9); Style.push_back(9); // fancy dashed lines 
	Style.push_back(1); Style.push_back(1); Style.push_back(1); Style.push_back(1);

	BreakDownColors.push_back(kBlue); BreakDownColors.push_back(kCyan); BreakDownColors.push_back(kGreen); BreakDownColors.push_back(kMagenta);

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_RadCorr");

//	FSIModel.push_back("hA2018_Final_NoRadCorr"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");
//	FSIModel.push_back("hA2018_Truth_NoRadCorr"); FSILabel.push_back("Genie (Truth)");  DirNames.push_back("hA2018_Truth_NoRadCorr");
//	FSIModel.push_back("hN2018_Final_NoRadCorr"); FSILabel.push_back("Genie hN2018");  DirNames.push_back("hN2018_Truth_NoRadCorr");

//	NameOfPlots.push_back("MissMomentum"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} P_{miss}^{#perp} [GeV/c]"); OutputPlotNames.push_back("MissMomentum");
//	NameOfPlots.push_back("MissMomentum_NoWeight"); LabelOfPlots.push_back("P_{miss}^{#perp} [GeV/c]"); OutputPlotNames.push_back("MissMomentum_NoWeight");
	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
//	NameOfPlots.push_back("eRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_1"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]");  OutputPlotNames.push_back("epRecoEnergy_slice_1");
//	NameOfPlots.push_back("h1_Erec_p_bkgd_slice_sub2p1pi_2p_1"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_1");
//	NameOfPlots.push_back("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_2"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]");  OutputPlotNames.push_back("epRecoEnergy_slice_2");
//	NameOfPlots.push_back("h1_Erec_p_bkgd_slice_sub2p1pi_2p_2"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_2");
//	NameOfPlots.push_back("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_3"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]");  OutputPlotNames.push_back("epRecoEnergy_slice_3");
//	NameOfPlots.push_back("h1_Erec_p_bkgd_slice_sub2p1pi_2p_3"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_3");
//	NameOfPlots.push_back("h_Etot_subtruct_piplpimi_factor_fracfeed"); LabelOfPlots.push_back("E^{cal} Feeddown");  OutputPlotNames.push_back("EcalReso");
//	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_factor_fracfeed"); LabelOfPlots.push_back("E^{QE} Feeddown"); OutputPlotNames.push_back("EQEReso");
//	NameOfPlots.push_back("h1_el_mom_corr");  LabelOfPlots.push_back("P_{e} (GeV / c)"); OutputPlotNames.push_back("Pmu");
//	NameOfPlots.push_back("h1_prot_mom"); LabelOfPlots.push_back("P_{p} (GeV / c)"); OutputPlotNames.push_back("Pp");
//	NameOfPlots.push_back("h1_theta0"); LabelOfPlots.push_back("#theta_{0} (beam-reco angle) (degrees)"); OutputPlotNames.push_back("theta0");
//	NameOfPlots.push_back("h1_Npi"); LabelOfPlots.push_back("#pi^{+/-} Multiplicity"); OutputPlotNames.push_back("PionMultiPlot");
//	NameOfPlots.push_back("h1_Nprot"); LabelOfPlots.push_back("Proton Multiplicity"); OutputPlotNames.push_back("Nproton");
//	NameOfPlots.push_back("h1_Nphot"); LabelOfPlots.push_back("N_{#gamma}"); OutputPlotNames.push_back("Ngamma");
//	NameOfPlots.push_back("h1_Npiphot_norad"); LabelOfPlots.push_back("N_{#gamma,#pi}"); OutputPlotNames.push_back("Ngammapi");
//	NameOfPlots.push_back("h1_Q2_weight"); LabelOfPlots.push_back("Q^{2} [GeV^{2} / c^{2}]"); OutputPlotNames.push_back("Q2");
//	NameOfPlots.push_back("h1_xbjk_weight"); LabelOfPlots.push_back("x_{B}"); OutputPlotNames.push_back("xB");
//	NameOfPlots.push_back("h1_nu_weight"); LabelOfPlots.push_back("Energy Transfer [GeV]"); OutputPlotNames.push_back("nu");
//	NameOfPlots.push_back("h1_Wvar_weight"); LabelOfPlots.push_back("W [GeV / c^{2}]"); OutputPlotNames.push_back("W");
//	NameOfPlots.push_back("h1_Q2Cal_weight"); LabelOfPlots.push_back("Reconstructed Q^{2} [GeV^{2} / c^{2}]"); OutputPlotNames.push_back("Q2Cal");
//	NameOfPlots.push_back("h1_xbjkCal_weight"); LabelOfPlots.push_back("Reconstructed x_{B}"); OutputPlotNames.push_back("xBCal");
//	NameOfPlots.push_back("h1_nuCal_weight"); LabelOfPlots.push_back("Reconstructed Energy Transfer [GeV]"); OutputPlotNames.push_back("nuCal");
//	NameOfPlots.push_back("h1_WvarCal_weight"); LabelOfPlots.push_back("Reconstructed W [GeV / c^{2})]"); OutputPlotNames.push_back("WCal");
//	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_3pi"); LabelOfPlots.push_back("(e,e')_{0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("InclusiveeRecoEnergy_slice_0");
//	NameOfPlots.push_back("CosDeltaThetaElectronPhotonAboveThreshold"); LabelOfPlots.push_back("cos(#Delta#theta_{e',#gamma})");  OutputPlotNames.push_back("CosDeltaThetaElectronPhotonAboveThreshold");
//	NameOfPlots.push_back("CosDeltaPhiElectronPhotonAboveThreshold"); LabelOfPlots.push_back("cos(#Delta#phi_{e',#gamma})");  OutputPlotNames.push_back("CosDeltaPhiElectronPhotonAboveThreshold");
//	NameOfPlots.push_back("h1_E_tot_cut2"); LabelOfPlots.push_back("(e,e')_{0#pi} E^{Cal} Before Subtraction [GeV]");  OutputPlotNames.push_back("h1_E_tot_cut2");

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

		// Loop over the energies

		for (int WhichEnergy = 0; WhichEnergy < NEnergies; WhichEnergy ++) {

			// Loop over the nuclei

			for (int WhichNucleus = 0; WhichNucleus < NNuclei; WhichNucleus ++) {

				// Loop over the plots

				for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

					TCanvas* PlotCanvas = new TCanvas(nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+NameOfPlots[WhichPlot]+"_"+xBCut[WhichxBCut],
									 205,34,1024,768);
//									 205,34,768,768);

					// ---------------------------------------------------------------------------------------------------------------------------

					// Dimensions of TPads (pad2 will be deleted at the very end for the Ereco plots)

					double XMinPadOne = -99., XMaxPadOne = -99., YMinPadOne = -99., YMaxPadOne = -99.;
					double XMinPadTwo = -99., XMaxPadTwo = -99., YMinPadTwo = -99., YMaxPadTwo = -99.;

					if (
						string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos || 
						string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos ||
						NameOfPlots[WhichPlot] == "h1_Q2_weight"||
						NameOfPlots[WhichPlot] == "h1_nu_weight" ||
						NameOfPlots[WhichPlot] == "h1_Nphot" || 
						NameOfPlots[WhichPlot] == "h1_Nprot" ||
						NameOfPlots[WhichPlot] == "h1_Npi" ||
						NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
						NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed"
					   ) {

						XMinPadOne = 0., XMaxPadOne = 1., YMinPadOne = 0.0, YMaxPadOne = 1.;
						XMinPadTwo = 0., XMaxPadTwo = 1., YMinPadTwo = 0.0, YMaxPadTwo = 0.01;
					} else {
						XMinPadOne = 0., XMaxPadOne = 1., YMinPadOne = 0.2, YMaxPadOne = 1.;
						XMinPadTwo = 0., XMaxPadTwo = 1., YMinPadTwo = 0., YMaxPadTwo = 0.2;
					}

					// ----------------------------------------------------------------------------------------

					TPad* pad1 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadOne,YMinPadOne,XMaxPadOne,YMaxPadOne, 21); 
					pad1->SetFillColor(kWhite); pad1->Draw();
					TPad* pad2 = new TPad(NameOfPlots[WhichPlot],NameOfPlots[WhichPlot],XMinPadTwo,YMinPadTwo,XMaxPadTwo,YMaxPadTwo,22); 
					pad2->SetFillColor(kWhite); pad2->Draw(); 
					pad1->SetTopMargin(0.1);
					pad1->SetBottomMargin(0.19);
					if (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" /*||
							 NameOfPlots[WhichPlot] == "h1_Npi"*/) { pad1->SetLeftMargin(0.12); }
					else if (
						OutputPlotNames[WhichPlot]=="epRecoEnergy_slice_1" || 
						OutputPlotNames[WhichPlot]=="epRecoEnergy_slice_2" || 
						OutputPlotNames[WhichPlot]=="epRecoEnergy_slice_3" || 
						OutputPlotNames[WhichPlot]=="MissMomentum"
						)
						{ pad1->SetLeftMargin(0.11); }
					else { pad1->SetLeftMargin(0.1); }
					pad1->SetRightMargin(0.04);

					if (NameOfPlots[WhichPlot] == "h1_Nprot") {
						pad1->SetRightMargin(0.0);
					}

					if (NameOfPlots[WhichPlot] == "h1_Npi") {
						pad1->SetLeftMargin(0.0);
						pad1->SetRightMargin(0.12);
					}

					pad2->SetTopMargin(0.21);
					pad1->cd();

					// ---------------------------------------------------------------------------------------

					Plots.clear();

////					TLegend* leg = new TLegend(0.17,0.92,0.9,0.99);
//					TLegend* leg = new TLegend(0.3,0.92,0.9,0.99);
//					leg->SetNColumns(3);
////					leg->SetNColumns(2);

					double LegXmin = 0.7, LegYmin = 0.52, YSpread = 0.35;
					if (xBCut[WhichxBCut] == "xBCut") { LegXmin = 0.6; }
					if ( OutputPlotNames[WhichPlot]=="InclusiveeRecoEnergy_slice_0" ) { LegXmin = 0.14; LegYmin = 0.45; }
					if ( OutputPlotNames[WhichPlot]=="epRecoEnergy_slice_0" ) { LegXmin = 0.15; LegYmin = 0.5; YSpread = 0.35; }
					if ( OutputPlotNames[WhichPlot]=="MissMomentum" ) { LegXmin = 0.6; LegYmin = 0.5; YSpread = 0.35; }

					TLegend* legGenie = new TLegend(LegXmin,LegYmin,LegXmin+0.15,LegYmin+YSpread);
					legGenie->SetNColumns(1);

					TLegend* legGenieBlackLine = new TLegend(LegXmin,0.68,LegXmin+0.15,0.82);
					legGenieBlackLine->SetNColumns(1);

					TLegend* legGenieBreak = new TLegend(LegXmin,0.55,0.4,0.68);
					legGenieBreak->SetNColumns(2);

					double max = -99.;
					double min = 1E12;

					// Loop over the FSI Models

					for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

						TString PathToFiles = "../myFiles/"+ E[WhichEnergy] + "/"+FSIModel[WhichFSIModel]+"/"+xBCut[WhichxBCut]+"/";
						TString FileName = PathToFiles+nucleus[WhichNucleus]+"_"+E[WhichEnergy]+"_"+FSIModel[WhichFSIModel]+"_Plots_FSI_em.root";
						TFile* FileSample = TFile::Open(FileName);

						Plots.push_back( (TH1D*)( FileSample->Get(NameOfPlots[WhichPlot]) ) );

						Plots[WhichFSIModel]->SetLineColor(Colors[WhichFSIModel]);
						CenterAxisTitle(Plots[WhichFSIModel]);

						// --------------------------------------------------------------------------------------

						// X-axis label

						Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleOffset(1.05);

						// X-axis Title

						if (
							//string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos || 
							NameOfPlots[WhichPlot] == "h1_Q2_weight" || NameOfPlots[WhichPlot] == "h1_nu_weight" || 
							NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" ||
							NameOfPlots[WhichPlot] == "h1_Npi" ||
							NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
							NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed"
						) 
							{ Plots[WhichFSIModel]->GetXaxis()->SetTitle(LabelOfPlots[WhichPlot]); }
						else { Plots[WhichFSIModel]->GetXaxis()->SetTitle(JustNucleus[WhichNucleus]+LabelOfPlots[WhichPlot]); }

						// Y-axis Title/Tick Size

						if (	
							OutputPlotNames[WhichPlot] =="MissMomentum" ||
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_1" ||
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_0" ||
							NameOfPlots[WhichPlot] == "h1_Nphot" || 
							NameOfPlots[WhichPlot] == "h1_Nprot" ||
							NameOfPlots[WhichPlot] == "h1_Npi" ||
							NameOfPlots[WhichPlot] == "h1_nu_weight"
						)
						 { 
							Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize); 
							Plots[WhichFSIModel]->GetYaxis()->SetTickSize(0.02);
						}
						else if (
							OutputPlotNames[WhichPlot] =="InclusiveeRecoEnergy_slice_0" ||
							NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
							NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed"
						) { 
							Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize); 
							Plots[WhichFSIModel]->GetYaxis()->SetTickSize(0.);
						}
						else { 
							Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(0.);
							Plots[WhichFSIModel]->GetYaxis()->SetTickSize(0.);
						}

						// --------------------------------------------------------------------------------------

						// Y-axis label

						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(0.);
						if (
							NameOfPlots[WhichPlot] == "h1_Nphot" || 
							NameOfPlots[WhichPlot] == "h1_Nprot" || 
							/*NameOfPlots[WhichPlot] == "h1_Npi" ||*/
							OutputPlotNames[WhichPlot] =="InclusiveeRecoEnergy_slice_0" ||
							OutputPlotNames[WhichPlot] =="MissMomentum" ||
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_0" ||
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_1" ||
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_2" ||
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_3" ||
							OutputPlotNames[WhichPlot] =="nu"
						) { Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize); }

						if (string(OutputPlotNames[WhichPlot]).find("_NoWeight") != std::string::npos) 
							{ Plots[WhichFSIModel]->GetYaxis()->SetTitle("Unweighted Events / GeV"); }
						else if ( (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" /*|| 
							NameOfPlots[WhichPlot] == "h1_Npi"*/) ) 
							{ Plots[WhichFSIModel]->GetYaxis()->SetTitle("Weighted Events"); }
						else if (OutputPlotNames[WhichPlot] =="InclusiveeRecoEnergy_slice_0" || 
							 OutputPlotNames[WhichPlot] =="MissMomentum" ||
							 OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_1" ||
							   OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_0" ||
							 OutputPlotNames[WhichPlot] == "nu"
						) 
							{ Plots[WhichFSIModel]->GetYaxis()->SetTitle("Weighted Events / GeV");  }
						else if (
							NameOfPlots[WhichPlot] == "h1_Nphot" || 
							NameOfPlots[WhichPlot] == "h1_Nprot" ||
							/*NameOfPlots[WhichPlot] == "h1_Npi" ||*/
							NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
							NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed"
						) { Plots[WhichFSIModel]->GetYaxis()->SetTitle("Weighted Events"); }


						// --------------------------------------------------------------------------------------

						Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						//Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(0);
						if (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" ||
							 NameOfPlots[WhichPlot] == "h1_Npi") { Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.75); }
//						else { Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.3); }
						else if (
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_1" ||
							OutputPlotNames[WhichPlot] =="epRecoEnergy_slice_2" ||
							OutputPlotNames[WhichPlot] =="MissMomentum"
						)
							{ Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.65); }
						else { Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.6); }
						Plots[WhichFSIModel]->SetLineWidth(LineWidth);

//						if (FSILabel[WhichFSIModel] == "Data" 
//							&& !(NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" || NameOfPlots[WhichPlot] == "h1_Npi")) 
//							{ leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "lep");}
//						else { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l"); }
//						leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l");

						// --------------------------------------------------------------------------------------

						// Transverse Missing Momentum Percentages Calculation

						TLatex* LowPercPmiss = new TLatex();
						TLatex* MidPercPmiss = new TLatex();
						TLatex* HighPercPmiss = new TLatex();
						double LowPmiss = 0.2, MidPmiss = 0.4, HighPmiss = 1.;
//						double LowPmiss = 0.3, MidPmiss = 10., HighPmiss = 11.;
						TString LowPercPmissString = "", MidPercPmissString = "", HighPercPmissString = "";

						if (string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos && xBCut[WhichxBCut] == "NoxBCut") { 

							int LowPmissBin = Plots[WhichFSIModel]->FindBin(LowPmiss);
							int MidPmissBin = Plots[WhichFSIModel]->FindBin(MidPmiss);
							int HighPmissBin = Plots[WhichFSIModel]->FindBin(HighPmiss);

							int PercLowPmiss = Plots[WhichFSIModel]->Integral(1,LowPmissBin) / Plots[WhichFSIModel]->Integral() * 100.;
							int PercMidPmiss = Plots[WhichFSIModel]->Integral(LowPmissBin+1,MidPmissBin) / Plots[WhichFSIModel]->Integral() * 100.;
							int PercHighPmiss = 100. - PercLowPmiss - PercMidPmiss;
							PercMidPmiss = 100. - PercLowPmiss;

							LowPercPmissString = ToString(PercLowPmiss)+"%";
							LowPercPmiss->SetTextFont(FontStyle);
							LowPercPmiss->SetTextColor(Colors[WhichFSIModel]);
							LowPercPmiss->SetTextSize(TextSize);

							MidPercPmissString = ToString(PercMidPmiss)+"%";
							MidPercPmiss->SetTextFont(FontStyle);
							MidPercPmiss->SetTextColor(Colors[WhichFSIModel]);
							MidPercPmiss->SetTextSize(TextSize);

							HighPercPmissString = ToString(PercHighPmiss)+"%";
							HighPercPmiss->SetTextFont(FontStyle);
							HighPercPmiss->SetTextColor(Colors[WhichFSIModel]);
							HighPercPmiss->SetTextSize(TextSize);

						}

//						// ---------------------------------------------------------------------------------------------------------

						double range = 0.05;
						if (OutputPlotNames[WhichPlot] =="InclusiveeRecoEnergy_slice_0" /*|| 
						string(OutputPlotNames[WhichPlot]).find("eRecoEnergy_slice") != std::string::npos*/) { range = 0.1; }

						// Energy Reconstruction Percentages Calculation

						TLatex* LowPercEReco = new TLatex();
						TLatex* MidPercEReco = new TLatex();
						TLatex* HighPercEReco = new TLatex();
						TString LowPercERecoString = "", MidPercERecoString = "", HighPercERecoString = "";

						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos) {

							// QE Energy Reconstruction

							if (string(OutputPlotNames[WhichPlot]).find("eRecoEnergy_slice") != std::string::npos) {

								double LowE = (1.-range)*DoubleE[WhichEnergy], HighE = (1.+range)*DoubleE[WhichEnergy];
								int LowEBin = Plots[WhichFSIModel]->FindBin(LowE);
								int HighEBin = Plots[WhichFSIModel]->FindBin(HighE);

								int PercLowPmiss = Plots[WhichFSIModel]->Integral(1,LowEBin) / Plots[WhichFSIModel]->Integral() * 100.;
								int PercMidPmiss = Plots[WhichFSIModel]->Integral(LowEBin+1,HighEBin) / Plots[WhichFSIModel]->Integral() * 100.;
								int PercHighPmiss = 100 - PercLowPmiss - PercMidPmiss;

								LowPercERecoString = ToString(PercLowPmiss)+"%";
								LowPercEReco->SetTextFont(FontStyle);
								LowPercEReco->SetTextColor(Colors[WhichFSIModel]);
								LowPercEReco->SetTextSize(TextSize);

								MidPercERecoString = ToString(PercMidPmiss)+"%";
								MidPercEReco->SetTextFont(FontStyle);
								MidPercEReco->SetTextColor(Colors[WhichFSIModel]);
								MidPercEReco->SetTextSize(TextSize);

								HighPercERecoString = ToString(PercHighPmiss)+"%";
								HighPercEReco->SetTextFont(FontStyle);
								HighPercEReco->SetTextColor(Colors[WhichFSIModel]);
								HighPercEReco->SetTextSize(TextSize);

							}

							// Calorimetric Energy Reconstruction Percentages Calculation

							else {

								double LowE = (1.-range)*DoubleE[WhichEnergy];
								int LowEBin = Plots[WhichFSIModel]->FindBin(LowE);

								int PercLowPmiss = Plots[WhichFSIModel]->Integral(1,LowEBin) / Plots[WhichFSIModel]->Integral() * 100.;
								int PercHighPmiss = 100 - PercLowPmiss;

								LowPercERecoString = ToString(PercLowPmiss)+"%";
								LowPercEReco->SetTextFont(FontStyle);
								LowPercEReco->SetTextColor(Colors[WhichFSIModel]);
								LowPercEReco->SetTextSize(TextSize);

								HighPercERecoString = ToString(PercHighPmiss)+"%";
								HighPercEReco->SetTextFont(FontStyle);
								HighPercEReco->SetTextColor(Colors[WhichFSIModel]);
								HighPercEReco->SetTextSize(TextSize);

								// -----------------------------------------------------------------------------------------------

							}

						}

						// -------------------------------------------------------------------------------------

						// Get percentages of events within X % of peak energy (5% for ECal/ 10% for EQE)

						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos) {

							double MaxContent = Plots[WhichFSIModel]->GetMaximum();
							int BinWithMaxContent = -99;

							int Nbins = Plots[WhichFSIModel]->GetXaxis()->GetNbins();

							for (int i = 0; i < Nbins; i++) {
								if (Plots[WhichFSIModel]->GetBinContent(i+1) ==  MaxContent) { BinWithMaxContent = i+1; }
							}

							TAxis *xaxis = Plots[WhichFSIModel]->GetXaxis();
							double binCenter = xaxis->GetBinCenter(BinWithMaxContent);

							// With respect to the peak energy
							double MinE = (1.-range)*binCenter;
							double MaxE = (1.+range)*binCenter;

							// With respect to the true beam energy
//							double MinE = (1.-range)*DoubleE[WhichEnergy];
//							double MaxE = (1.+range)*DoubleE[WhichEnergy];

							int MinBin = Plots[WhichFSIModel]->FindBin(MinE);
							int MaxBin = Plots[WhichFSIModel]->FindBin(MaxE);
							int percentage = Plots[WhichFSIModel]->Integral(MinBin,MaxBin) / Plots[WhichFSIModel]->Integral() * 100.;

							cout << endl;	
							cout << FSILabel[WhichFSIModel] << "  " << LabelOfPlots[WhichPlot] << ": fraction within " << int(range*100.) << "% = " << percentage << endl;
							cout << endl;

						}

						//-----------------------------------------------------------------------------------------------

						//Larry's suggestion because ECal has a sharp peak and a low tail 
						//Thus we multiply the peak by EnhaceTail

						if ( (OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_0" || OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_1") 
							&& DoubleE[WhichEnergy] == 2.261 && nucleus[WhichNucleus] == "12C" ) {

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

						// --------------------------------------------------------------------------------------

						// Scaling Factor

//						double ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral(); // default

						double ScalingFactor = 1.;

//						double ScalingFactor = 1. / Plots[WhichFSIModel]->Integral(); // area normalized

						if (
							NameOfPlots[WhichPlot] == "h1_Nphot" || 
							NameOfPlots[WhichPlot] == "h1_Nprot" ||
							NameOfPlots[WhichPlot] == "h1_Npi") { ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral();}

						else { ScalingFactor = 1. / Plots[WhichFSIModel]->Integral(); }  // area normalized

//						double ScalingFactor = 1. / Plots[WhichFSIModel]->GetMaximum(); // peak at 1

//						double ScalingFactor = Plots[0]->GetEntries() / Plots[WhichFSIModel]->GetEntries();
//						double ScalingFactor = Plots[0]->GetMaximum() / Plots[WhichFSIModel]->GetMaximum();
//						double ScalingFactor = 18E8 / Plots[WhichFSIModel]->GetMaximum();
//						double ScalingFactor = 1.;
						Plots[WhichFSIModel]->Scale(ScalingFactor);

						// -----------------------------------------------------------------------------------

						// Accounting for the fact that the bin width might not be constant

						if ( !(
							NameOfPlots[WhichPlot] == "h1_Nphot" || 
							NameOfPlots[WhichPlot] == "h1_Nprot" ||
							NameOfPlots[WhichPlot] == "h1_Npi") ) { ReweightPlots(Plots[WhichFSIModel]); }

						// --------------------------------------------------------------------------------------

						// Rebining & ranges

//						if (string(OutputPlotNames[WhichPlot]).find("epRecoEnergy_slice") != std::string::npos) 
//							{ for (int i = 0; i < NECalRebin; i++) { Plots[WhichFSIModel]->Rebin(); } }

						if (string(OutputPlotNames[WhichPlot]).find("epRecoEnergy_slice") != std::string::npos && nucleus[WhichNucleus] == "12C" 
							&& DoubleE[WhichEnergy] == 2.261) 
							{ Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.6,2.5); }

						if (string(OutputPlotNames[WhichPlot]).find("epRecoEnergy_slice") != std::string::npos && nucleus[WhichNucleus] == "56Fe" 
							&& DoubleE[WhichEnergy] == 4.461) 
							{ Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(1.,5.); }

						if (string(OutputPlotNames[WhichPlot]).find("eRecoEnergy_slice") != std::string::npos) 
							{ for (int i = 0; i < 0; i++) { Plots[WhichFSIModel]->Rebin(); } }

//						if (NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed") 
//							{ for (int i = 0; i < 1; i++) { Plots[WhichFSIModel]->Rebin(); } }

						if (NameOfPlots[WhichPlot] == "h1_el_mom_corr") { for (int i = 0; i < 1; i++) 
							{ Plots[WhichFSIModel]->Rebin(); } }

						if (NameOfPlots[WhichPlot] == "h1_prot_mom") 
							{ for (int i = 0; i < 2; i++) { Plots[WhichFSIModel]->Rebin(); } }

						if (NameOfPlots[WhichPlot] == "h1_Wvar_weight" || NameOfPlots[WhichPlot] == "h1_WvarCal_weight") 
							{ for (int i = 0; i < 3; i++) { Plots[WhichFSIModel]->Rebin(); } }

						if (NameOfPlots[WhichPlot] == "h1_xbjk_weight" || NameOfPlots[WhichPlot] == "h1_xbjkCal_weight") 
							{ for (int i = 0; i < 2; i++) { Plots[WhichFSIModel]->Rebin(); } Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.,2.);}

						if (NameOfPlots[WhichPlot] == "h1_Q2_weight" || NameOfPlots[WhichPlot] == "h1_Q2Cal_weight") {
							for (int i = 0; i < 2; i++) { Plots[WhichFSIModel]->Rebin();} 
							Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.,2.);
							if (DoubleE[WhichEnergy] == 4.461) { 
								for (int i = 0; i < 2; i++) { Plots[WhichFSIModel]->Rebin(); } 
								Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(1.,4.);  }
						}

						if ( (NameOfPlots[WhichPlot] == "h1_Q2_weight" || NameOfPlots[WhichPlot] == "h1_Q2Cal_weight") && xBCut[WhichxBCut] == "xBCut") 
							{ Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.4,2.5); }

						if (NameOfPlots[WhichPlot] == "h1_theta0") 
							{ Plots[WhichFSIModel]->Rebin(); Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.,60.); }

						if (NameOfPlots[WhichPlot] == "h1_nu_weight" || NameOfPlots[WhichPlot] == "h1_nuCal_weight") 
							{ for (int i = 0; i < 3; i++) {Plots[WhichFSIModel]->Rebin();} Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.,1.7); }

						if ( (NameOfPlots[WhichPlot] == "h1_nu_weight" || NameOfPlots[WhichPlot] == "h1_nuCal_weight") && xBCut[WhichxBCut] == "xBCut") 
							{ Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.2,1.2); }

						if (string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos) { for (int i = 0; i < 2; i++) 
							{ Plots[WhichFSIModel]->Rebin();} }

						if (NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed") 
							{ for (int i = 0; i < 0; i++) {Plots[WhichFSIModel]->Rebin();} Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(-0.7,0.1); }

						if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed" && xBCut[WhichxBCut] == "NoxBCut") 
							{ for (int i = 0; i < 0; i++) {Plots[WhichFSIModel]->Rebin();} Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(-0.8,0.3); }

						if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed" && xBCut[WhichxBCut] == "xBCut") 
							{ for (int i = 0; i < 0; i++) {Plots[WhichFSIModel]->Rebin();} Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(-0.2,0.2); }

//						// --------------------------------------------------------------------------------------

//						// Scaling Factor

////						double ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral(); // default

//						double ScalingFactor = 1. / Plots[WhichFSIModel]->Integral(); // area normalized

////						double ScalingFactor = 1. / Plots[WhichFSIModel]->GetMaximum(); // peak at 1

////						double ScalingFactor = Plots[0]->GetEntries() / Plots[WhichFSIModel]->GetEntries();
////						double ScalingFactor = Plots[0]->GetMaximum() / Plots[WhichFSIModel]->GetMaximum();
////						double ScalingFactor = 18E8 / Plots[WhichFSIModel]->GetMaximum();
////						double ScalingFactor = 1.;
//						Plots[WhichFSIModel]->Scale(ScalingFactor);

						// ----------------------------------------------------------------------------------

						// Apply Systematic Uncertainties on Data Points

						double SystUnc = 0;
						if ( DoubleE[WhichEnergy] == 1.161 ) { SystUnc = SystUnc1GeV; }
						if ( DoubleE[WhichEnergy] == 2.261 ) { SystUnc = SystUnc2GeV; }
						if ( DoubleE[WhichEnergy] == 4.461 ) { SystUnc = SystUnc4GeV; }

						if (FSILabel[WhichFSIModel] == "Data") { ApplySystUnc(Plots[WhichFSIModel], SystUnc); }

						// ----------------------------------------------------------------------------------

						// Genie Break Down
/*
						if (
							FSILabel[WhichFSIModel] == "Genie" && 
							//FSILabel[WhichFSIModel] == "Rad" &&
							( (NameOfPlots[WhichPlot] == "MissMomentum" )  || 
							NameOfPlots[WhichPlot] == "h1_Q2_weight"|| 
							NameOfPlots[WhichPlot] == "h1_nu_weight" ||
							NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_3pi" ||
							NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0" ||
							NameOfPlots[WhichPlot] == "h1_el_mom_corr"
							) 

						) {

							//if ( xBCut[WhichxBCut] == "NoxBCut") { 
								//if (Plots[WhichPlot] == "MissMomentum") {
								legGenie->AddEntry(Plots[0],"Data", "lep"); 
								legGenieBlackLine->AddEntry(Plots[0],"Data", "lep"); 
								//}
								legGenie->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 
								legGenieBlackLine->AddEntry(Plots[WhichFSIModel],"GENIE (Total)", "l"); 
							//}
							//else { legGenie->AddEntry(Plots[WhichFSIModel],"GENIE", "l"); }

							BreakDownPlots.clear();

							for (int j = 1; j < 5; j++) {

								if (NameOfPlots[WhichPlot] == "MissMomentum") 
									{ BreakDownPlots.push_back( (TH1D*)( FileSample->Get("Pmiss_Int_"+ToString(j)) ) ); }
								if (NameOfPlots[WhichPlot] == "h1_Q2_weight") 
									{ BreakDownPlots.push_back( (TH1D*)( FileSample->Get("Q2_Int_"+ToString(j)) ) ); }
								if (NameOfPlots[WhichPlot] == "h1_nu_weight") 
									{ BreakDownPlots.push_back( (TH1D*)( FileSample->Get("Nu_Int_"+ToString(j)) ) ); }
								if (NameOfPlots[WhichPlot] == "h1_el_mom_corr") 
									{ BreakDownPlots.push_back( (TH1D*)( FileSample->Get("Pe_Int_"+ToString(j)) ) ); }
								if (NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_noprot_3pi") 
								{ 
									BreakDownPlots.push_back( (TH1D*)( FileSample->Get("InclusiveEQE_Int_"+ToString(j)) ) ); 
									if (xBCut[WhichxBCut] == "xBCut" && DoubleE[WhichEnergy] == 1.161) {
										BreakDownPlots[j-1]->GetXaxis()->SetRangeUser(1.05,1.35);
										Plots[0]->GetXaxis()->SetRangeUser(1.05,1.3);
									}

									if (xBCut[WhichxBCut] == "xBCut" && DoubleE[WhichEnergy] == 4.461) {
										BreakDownPlots[j-1]->GetXaxis()->SetRangeUser(3.5,5.5);
										Plots[0]->GetXaxis()->SetRangeUser(3.5,5.5);
									}
								}
								if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0") 
								{ 
									BreakDownPlots.push_back( (TH1D*)( FileSample->Get("ECal_Int_"+ToString(j)) ) ); 
									if (xBCut[WhichxBCut] == "xBCut" && nucleus[WhichNucleus] == "12C") {
										BreakDownPlots[j-1]->GetXaxis()->SetRangeUser(1.5,2.4);
										Plots[0]->GetXaxis()->SetRangeUser(1.5,2.4);
									}

									if (xBCut[WhichxBCut] == "xBCut" && nucleus[WhichNucleus] == "56Fe") {
										BreakDownPlots[j-1]->GetXaxis()->SetRangeUser(1.7,4.8);
										Plots[0]->GetXaxis()->SetRangeUser(1.7,4.8);
									}
								}

								ReweightPlots(BreakDownPlots[j-1]);
								if ( 
									NameOfPlots[WhichPlot] != "h_Erec_subtruct_piplpimi_noprot_3pi" 
									&& NameOfPlots[WhichPlot] != "epRecoEnergy_slice_0"
									&& NameOfPlots[WhichPlot] != "h1_el_mom_corr"
								) 
									{ for (int i = 0; i < 2; i++) { BreakDownPlots[j-1]->Rebin(); } }
								if (NameOfPlots[WhichPlot] == "h1_nu_weight") { BreakDownPlots[j-1]->Rebin(); }
								//if (NameOfPlots[WhichPlot] == "epRecoEnergy_slice_0") { BreakDownPlots[j-1]->Rebin(); }
								if (NameOfPlots[WhichPlot] == "h1_el_mom_corr") { BreakDownPlots[j-1]->Rebin(); }

								//-----------------------------------------------------------------------------------------------

								//Larry's suggestion because ECal has a sharp peak and a low tail 
								//Thus we multiply the peak by EnhaceTail

								if ( (OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_0") 
									&& DoubleE[WhichEnergy] == 2.261 && nucleus[WhichNucleus] == "12C" ) {

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

								BreakDownPlots[j-1]->SetLineColor(BreakDownColors[j-1]);

								if (    NameOfPlots[WhichPlot] != "h_Erec_subtruct_piplpimi_noprot_3pi" 
									&& NameOfPlots[WhichPlot] != "epRecoEnergy_slice_0") ) {

									int GenieNBins = Plots[WhichFSIModel]->GetNbinsX();
									int GenieMin = Plots[WhichFSIModel]->GetXaxis()->GetXmin();
									int GenieMax = Plots[WhichFSIModel]->GetXaxis()->GetXmax();
									BreakDownPlots[j-1]->SetBins(GenieNBins,GenieMin,GenieMax);
								}

								BreakDownPlots[j-1]->SetLineWidth(3);
								BreakDownPlots[j-1]->SetLineStyle(Style[j-1]);
								BreakDownPlots[j-1]->Scale(ScalingFactor);
								TLegendEntry* l1 = legGenie->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								l1->SetTextColor(BreakDownColors[j-1]);

								TLegendEntry* l1Break = legGenieBreak->AddEntry(BreakDownPlots[j-1],GenieFSILabel[j-1], "l");
								l1Break->SetTextColor(BreakDownColors[j-1]);

//								BreakDownPlots[j-1]->Draw("hist same");
								BreakDownPlots[j-1]->Draw("C hist same");

							} // end of the look over the GENIE break down

						}
*/
						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						double localmax = Plots[WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
//						double height = 1.4;
						double height = 1.05;
//						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice_3") != std::string::npos) { height = 1.5; }
						if ( xBCut[WhichxBCut] == "xBCut" ) { height = 1.1; }
						Plots[0]->GetYaxis()->SetRangeUser(0.,height*max);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						if (
							string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos ||
							NameOfPlots[WhichPlot] == "h1_Nphot" || 
							NameOfPlots[WhichPlot] == "h1_Nprot" ||
							NameOfPlots[WhichPlot] == "h1_Npi" ||
							NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
							NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed"
						)
						{ 
							Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(5); 
							//if (FSILabel[WhichFSIModel] == "Genie") { Plots[WhichFSIModel]->SetLineColor(kRed);} 
						}
						else { Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions); }
						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);

						// --------------------------------------------------------------------------------------------------

						// Multiplicity plots

						if (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" || NameOfPlots[WhichPlot] == "h1_Npi") {

							Plots[WhichFSIModel]->GetYaxis()->SetLabelOffset(-0.004);
							Plots[WhichFSIModel]->Rebin();
							Plots[0]->GetYaxis()->SetRangeUser(0.5*min,2.*max); pad1->SetLogy();
							if (FSILabel[WhichFSIModel] == "Data") { 
						 
								Plots[WhichFSIModel]->SetMarkerSize(3.); 
								if (NameOfPlots[WhichPlot] == "h1_Nprot") { 
									Plots[WhichFSIModel]->SetMarkerColor(kBlack); Plots[WhichFSIModel]->SetMarkerStyle(20); }
								else { 
									Plots[WhichFSIModel]->SetMarkerColor(kBlue); 
									Plots[WhichFSIModel]->SetMarkerStyle(24);
								}
								gStyle->SetErrorX(0); 
								Plots[WhichFSIModel]->GetYaxis()->SetRangeUser(1E1,1E8);
								Plots[WhichFSIModel]->Draw("e same"); 
							}
							else { 
//								if (NameOfPlots[WhichPlot] == "h1_Nprot") { Plots[WhichFSIModel]->SetLineColor(kBlack); }
//								else { Plots[WhichFSIModel]->SetLineColor(kBlue); Plots[WhichFSIModel]->SetLineStyle(7); }
								Plots[WhichFSIModel]->Draw("hist same"); 
								gStyle->SetErrorX(0); 
								Plots[0]->Draw("e same"); 
							}
						} 

						// All the other plots

						else { 
						//	Plots[WhichFSIModel]->Draw("C hist same"); 
							if (FSILabel[WhichFSIModel] == "Data") { 

								Plots[WhichFSIModel]->SetMarkerStyle(20); 
								Plots[WhichFSIModel]->SetMarkerSize(2.); 
								Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

								if ( OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_1" 
									|| OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_2"
									|| OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_3"  ) { Plots[WhichFSIModel]->SetMarkerSize(3.); }
								gStyle->SetErrorX(0); // Removing the horizontal errors
								Plots[WhichFSIModel]->Draw("e same"); 
							} else { 
//								Plots[WhichFSIModel]->Draw("hist same"); // draw them as histos
								Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines
								Plots[0]->Draw("e same"); 
							}
						}

						// ---------------------------------------------------------------------------------------------------

						if (
							NameOfPlots[WhichPlot] == "MissMomentum" && 
							!(nucleus[WhichNucleus] == "12C" && DoubleE[WhichEnergy] == 2.261) && 
							FSILabel[WhichFSIModel] == "Data" 
						) {

							TLatex* myEnergy = new TLatex();
							myEnergy->SetTextFont(FontStyle);
							myEnergy->SetTextColor(kBlack);
							myEnergy->SetTextSize(TextSize);
							TString myEnergyString = ToString(DoubleE[WhichEnergy])+" GeV";
							pad1->cd();
							myEnergy->DrawLatexNDC(0.65,0.8,myEnergyString);

						}

						// -----------------------------------------------------------------------------------

						double SpaceBetweenPerc = 0.08;

						// Transverse Missing Momentum Percentages

						if (
							string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos && 
							xBCut[WhichxBCut] == "NoxBCut" &&
							nucleus[WhichNucleus] == "12C" &&
							DoubleE[WhichEnergy] == 2.261
						) { 

							pad1->cd();
//							double LowPmiss = 0.2, MidPmiss = 0.4, HighPmiss = 1.;

							TLine* line1 = new TLine(LowPmiss,0.,LowPmiss,height*max);
							line1->SetLineColor(kBlack); line1->SetLineWidth(LineWidth);
							line1->Draw(); 

							TLine* line2 = new TLine(MidPmiss,0.,MidPmiss,height*max);
							line2->SetLineColor(kBlack); line2->SetLineWidth(LineWidth);
							line2->Draw(); 

//							LowPercPmiss->DrawTextNDC(0.15,0.83-WhichFSIModel*SpaceBetweenPerc,LowPercPmissString);
//							MidPercPmiss->DrawTextNDC(0.31,0.83-WhichFSIModel*SpaceBetweenPerc,MidPercPmissString);
//							HighPercPmiss->DrawTextNDC(0.5,0.83-WhichFSIModel*SpaceBetweenPerc,HighPercPmissString);

// Add the 2 lines below if we want the % in the Pmiss plot
//							LowPercPmiss->DrawTextNDC(0.25,0.83-WhichFSIModel*SpaceBetweenPerc,LowPercPmissString);
//							MidPercPmiss->DrawTextNDC(0.38,0.83-WhichFSIModel*SpaceBetweenPerc,MidPercPmissString);
//							HighPercPmiss->DrawTextNDC(0.5,0.83-WhichFSIModel*SpaceBetweenPerc,HighPercPmissString);

							TLatex* myNucleus = new TLatex();
							myNucleus->SetTextFont(FontStyle);
							myNucleus->SetTextColor(kBlack);
							myNucleus->SetTextSize(TextSize);
							myNucleus->DrawLatexNDC(0.65,0.83,JustNucleus[WhichNucleus]+"(e,e'p)_{1p0#pi}");
//							myNucleus->DrawLatexNDC(0.7,0.58,JustNucleus[WhichNucleus]+"(e,e'p)_{1p0#pi}");
							

						}

						// -----------------------------------------------------------------------------------------------------------------------------------

						// Energy Reconstruction

						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos) {

							pad1->cd();

							// QE Energy Reconstruction

							if (string(OutputPlotNames[WhichPlot]).find("eRecoEnergy_slice") != std::string::npos) {

////								LowPercEReco->DrawTextNDC(0.3,0.83-WhichFSIModel*SpaceBetweenPerc,LowPercERecoString);
//								LowPercEReco->DrawTextNDC(0.4,0.83-WhichFSIModel*SpaceBetweenPerc,LowPercERecoString);

//								MidPercEReco->DrawTextNDC(0.61,0.83-WhichFSIModel*SpaceBetweenPerc,MidPercERecoString);
//								HighPercEReco->DrawTextNDC(0.77,0.83-WhichFSIModel*SpaceBetweenPerc,HighPercERecoString);

								TLine* line1 = new TLine(LowE,0.,LowE,height*max);
								line1->SetLineColor(kBlack); line1->SetLineWidth(LineWidth);
								//line1->Draw(); 

								TLine* line2 = new TLine(HighE,0.,HighE,height*max);
								line2->SetLineColor(kBlack); line2->SetLineWidth(LineWidth);
								//line2->Draw(); 

							}

							// Calorimetric Energy Reconstruction

							else {

//								if ( OutputPlotNames[WhichPlot] != "epRecoEnergy_slice_3") {
//									if (DoubleE[WhichEnergy] == 4.461) { LowPercEReco->DrawTextNDC(0.37,0.83-WhichFSIModel*SpaceBetweenPerc,LowPercERecoString); }
//									else { LowPercEReco->DrawTextNDC(0.32,0.83-WhichFSIModel*SpaceBetweenPerc,LowPercERecoString); }
//								}

//								if ( OutputPlotNames[WhichPlot] != "epRecoEnergy_slice_3") {
//									if (DoubleE[WhichEnergy] == 4.461) { HighPercEReco->DrawTextNDC(0.79,0.83-WhichFSIModel*SpaceBetweenPerc,HighPercERecoString); }
//									else { HighPercEReco->DrawTextNDC(0.77,0.83-WhichFSIModel*SpaceBetweenPerc,HighPercERecoString); }
//								}

								TLine* line1 = new TLine(LowE,0.,LowE,height*max);
								line1->SetLineColor(kBlack); line1->SetLineWidth(LineWidth);
								if ( (OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_0" || 
								      OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_1") 
								   && nucleus[WhichNucleus] == "12C") { line1->Draw(); }

								// -----------------------------------------------------------------------------------------------

							}

						}

		                                // ---------------------------------------------------------------------------------------------------

		                                // Place of Data / Genie labelling

						TLatex latexData;
						latexData.SetTextFont(FontStyle);
						latexData.SetTextSize(0.1);
						latexData.SetTextColor(Colors[WhichFSIModel]);

						TLatex latexxB;
						latexxB.SetTextFont(FontStyle);
						latexxB.SetTextSize(0.1);
						latexxB.SetTextColor(Colors[WhichFSIModel]);

						TLatex latexScale;
						latexScale.SetTextFont(FontStyle);
						latexScale.SetTextSize(0.1);
						latexScale.SetTextColor(Colors[WhichFSIModel]);

						if ( (OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_0" || OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_1") 
							&& DoubleE[WhichEnergy] == 2.261 && FSILabel[WhichFSIModel] == "Data" )
							{ latexData.SetTextColor(kBlack); }	

//						if (string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos && xBCut[WhichxBCut] == "NoxBCut" 
//							&& nucleus[WhichNucleus] == "12C" && FSILabel[WhichFSIModel] == "Data") 
//							{ latexData.SetTextSize(TextSize);
//							latexData.DrawLatexNDC(0.55+0.2*WhichFSIModel,0.33-0.22*WhichFSIModel,FSILabel[WhichFSIModel]); }
//						if (string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos && xBCut[WhichxBCut] == "xBCut" && 
//							FSILabel[WhichFSIModel] == "Data") 
//							{ latexData.SetTextSize(TextSize);latexData.DrawLatexNDC(0.29+0.2*WhichFSIModel,0.49-0.22*WhichFSIModel,
//							FSILabel[WhichFSIModel]); }

//						if ( OutputPlotNames[WhichPlot] == "InclusiveeRecoEnergy_slice_0" && FSILabel[WhichFSIModel] == "Data" )
////							{ latexData.DrawLatexNDC(0.3-0.05*WhichFSIModel,0.22+0.3*WhichFSIModel,FSILabel[WhichFSIModel]); }
////							{ latexData.DrawLatexNDC(0.75-0.05*WhichFSIModel,0.43+0.3*WhichFSIModel,FSILabel[WhichFSIModel]); }
//							{ latexData.DrawLatexNDC(0.82-0.05*WhichFSIModel,0.43+0.3*WhichFSIModel,FSILabel[WhichFSIModel]); }

//						if ( OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_1") 
//							{ latexData.DrawLatexNDC(0.51-0.2*WhichFSIModel,0.3-0.05*WhichFSIModel,FSILabel[WhichFSIModel]); }
//						if ( OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_2") 
//							{ latexData.DrawLatexNDC(0.45-0.2*WhichFSIModel,0.4-0.05*WhichFSIModel,FSILabel[WhichFSIModel]); }
//						if ( OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_3") 
//							{ latexData.DrawLatexNDC(0.45-0.3*WhichFSIModel,0.68-0.05*WhichFSIModel,FSILabel[WhichFSIModel]); }

						if ( (OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_0")
							&& DoubleE[WhichEnergy] == 2.261 && nucleus[WhichNucleus] == "12C" ) {
 
							latexScale.SetTextSize(TextSize); 
							latexScale.DrawLatexNDC(0.86,0.47,"x1/3"); 
//							if (FSILabel[WhichFSIModel] == "Data") {
//								latexData.DrawLatexNDC(0.4,0.43+0.3*WhichFSIModel,FSILabel[WhichFSIModel]);
//							}
						}

						if ( (OutputPlotNames[WhichPlot] == "epRecoEnergy_slice_1")
							&& DoubleE[WhichEnergy] == 2.261 && nucleus[WhichNucleus] == "12C" ) 
							{ latexData.SetTextSize(TextSize); latexData.DrawLatexNDC(0.86,0.5,"x1/3"); }

//						if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos && DoubleE[WhichEnergy] == 2.261 
//							&& nucleus[WhichNucleus] == "56Fe" ) 
//							{ latexData.DrawLatexNDC(0.42-0.2*WhichFSIModel,0.42-0.1*WhichFSIModel,FSILabel[WhichFSIModel]); }

//						else if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos) 
//							{ latexData.DrawLatexNDC(0.42-0.2*WhichFSIModel,0.52-0.1*WhichFSIModel,FSILabel[WhichFSIModel]); }
//						else if (string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos) 
//							{ latexData.DrawLatexNDC(0.51-0.2*WhichFSIModel,0.32-0.05*WhichFSIModel,FSILabel[WhichFSIModel]); }

						if (
							/*NameOfPlots[WhichPlot] == "h1_Nphot" ||*/ 
							NameOfPlots[WhichPlot] == "h1_Nprot"
							/*|| NameOfPlots[WhichPlot] == "h1_Npi"*/
						)
							{ 
								if (FSILabel[WhichFSIModel] == "Genie") { latexData.SetTextColor(kBlack); }
								latexData.DrawLatexNDC(0.56-0.08*WhichFSIModel,0.4+0.35*WhichFSIModel,FSILabel[WhichFSIModel]); 
							}

						if (
							( (NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
							   NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed") && xBCut[WhichxBCut] == "NoxBCut")
						)
							{ 
								if (FSILabel[WhichFSIModel] == "Genie") { latexData.SetTextColor(kRed); latexxB.DrawLatexNDC(0.78,0.8,"All x_{B}");}
								latexData.DrawLatexNDC(0.47-0.2*WhichFSIModel,0.36+0.22*WhichFSIModel,FSILabel[WhichFSIModel]); 
							}

						if (
							( (NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
							   NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed") && xBCut[WhichxBCut] == "xBCut")
						)
							{ 
								if (FSILabel[WhichFSIModel] == "Genie") { latexxB.DrawLatexNDC(0.78,0.8,"x_{B} ~ 1");}
							}

		                                // --------------------------------------------------------------------------------------------------

					} // End of the loop over the FSI Models 

//					leg->SetBorderSize(0);
//					leg->SetTextFont(FontStyle);
//					leg->SetTextSize(TextSize);
////					leg->Draw();

					legGenie->SetBorderSize(0);
					legGenie->SetTextFont(FontStyle);

					legGenieBlackLine->SetBorderSize(0);
					legGenieBlackLine->SetTextFont(FontStyle);

					legGenieBreak->SetBorderSize(0);
					legGenieBreak->SetTextFont(FontStyle);

					if (xBCut[WhichxBCut] == "NoxBCut") { legGenie->SetTextSize(TextSize); }
					if (
						xBCut[WhichxBCut] == "xBCut" || 
						OutputPlotNames[WhichPlot]=="InclusiveeRecoEnergy_slice_0") 
					{ legGenie->SetTextSize(TextSize); }

//					if ( NameOfPlots[WhichPlot] == "MissMomentum" ) { legGenie->Draw(); }
					if ( OutputPlotNames[WhichPlot]=="InclusiveeRecoEnergy_slice_0"  && nucleus[WhichNucleus] == "12C" && DoubleE[WhichEnergy] == 1.161) { 

						legGenieBlackLine->SetNColumns(1); 
						legGenieBlackLine->SetTextSize(TextSize-0.03); 
						legGenieBlackLine->Draw(); 

						legGenieBreak->SetTextSize(TextSize-0.03); 
						legGenieBreak->Draw();

						TLatex* myNucleus = new TLatex();
						myNucleus->SetTextFont(FontStyle);
						myNucleus->SetTextColor(kBlack);
						myNucleus->SetTextSize(TextSize-0.02);
						myNucleus->DrawLatexNDC(0.14,0.85,JustNucleus[WhichNucleus]+"(e,e')_{0#pi}");

						TLatex* myEbeam = new TLatex();
						myEbeam->SetTextFont(FontStyle);
						myEbeam->SetTextColor(kAzure+4);
						myEbeam->SetTextSize(TextSize-0.02);
						myEbeam->DrawLatexNDC(0.67,0.34,"E_{beam}");

						TLatex* myArrow = new TLatex();
						myArrow->SetTextFont(FontStyle);
						myArrow->SetTextColor(kAzure+4);
						myArrow->SetTextSize(1.2*TextSize);
						myArrow->DrawLatex(1.161,0.1,"#Downarrow");

					}

					if ( 
						(NameOfPlots[WhichPlot] == "MissMomentum" && nucleus[WhichNucleus] == "12C" && DoubleE[WhichEnergy] == 2.261) 
						|| (OutputPlotNames[WhichPlot]=="epRecoEnergy_slice_0" && nucleus[WhichNucleus] == "12C" && DoubleE[WhichEnergy] == 2.261) 
					) 
						{ legGenie->SetTextSize(TextSize-0.04);legGenie->Draw(); }

					TLatex latex;
					latex.SetTextFont(FontStyle);
					latex.SetTextSize(TextSize);
//					if (NameOfPlots[WhichPlot] == "MissMomentum") { latex.DrawLatexNDC(0.5,0.65,LabelsOfSamples[WhichNucleus]+LabelE[WhichEnergy]); }
//					else if (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" || NameOfPlots[WhichPlot] == "h1_Npi") 
//						{ latex.DrawLatexNDC(0.5,0.7,LabelsOfSamples[WhichNucleus]+LabelE[WhichEnergy]); }
//					else { latex.DrawLatexNDC(0.14,0.65,LabelsOfSamples[WhichNucleus]+LabelE[WhichEnergy]); }

					pad2->cd(); pad2->SetGrid();
					Plots_Clones.clear();
		
					for (int WhichFSIModel = 0 ; WhichFSIModel < NFSIModels ; WhichFSIModel ++ ) {
		
						Plots_Clones.push_back( (TH1D*)(Plots[WhichFSIModel]->Clone())) ;
						Plots_Clones[WhichFSIModel]->GetXaxis()->SetTitleSize(0.0);
						Plots_Clones[WhichFSIModel]->GetXaxis()->SetLabelSize(0.0);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetRangeUser(0,2);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetLabelSize(0.13);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleSize(0.22);

						// Residual

						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.2);
//						Plots_Clones[WhichFSIModel]->Add(Plots[0],-1);
//						Plots_Clones[WhichFSIModel]->Divide(Plots[0]);
//						Plots_Clones[WhichFSIModel]->Scale(1./Plots[0]->Integral());
//						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitle("#frac{MC-Data}{Integral}");			
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetLabelSize(0.21);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleSize(0.21);
						Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.21);
//						Plots_Clones[WhichFSIModel]->GetYaxis()->SetRangeUser(-0.03,0.03);

						Plots_Clones[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);

//						Plots_Clones[WhichFSIModel]->Draw("hist same");

						if (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" || NameOfPlots[WhichPlot] == "h1_Npi") { 
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetLabelOffset(0.005);
							Plots_Clones[WhichFSIModel]->Divide(Plots[0]);

							Plots_Clones[WhichFSIModel]->GetYaxis()->SetNdivisions(3);
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitle("#frac{Genie}{Data}");	
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetRangeUser(0.,2.);
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetLabelSize(0.3);
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleSize(0.3);
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.13);

							Plots_Clones[WhichFSIModel]->Draw("hist same");
						} else { 
							Plots_Clones[WhichFSIModel]->Add(Plots[0],-1);
							Plots_Clones[WhichFSIModel]->Scale(1./Plots[0]->Integral());
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetTitle("#frac{Genie-Data}{Integral}");	
							Plots_Clones[WhichFSIModel]->GetYaxis()->SetRangeUser(-0.03,0.03);		
//							Plots_Clones[WhichFSIModel]->Draw("C hist same"); 
							Plots_Clones[WhichFSIModel]->Draw("hist same"); 
						}
											}

					if (
						string(OutputPlotNames[WhichPlot]).find("RecoEnergy_slice") != std::string::npos || 
						string(NameOfPlots[WhichPlot]).find("MissMomentum") != std::string::npos ||
                                                NameOfPlots[WhichPlot] == "h1_Q2_weight" ||
						NameOfPlots[WhichPlot] == "h1_nu_weight" ||
						NameOfPlots[WhichPlot] == "h1_Nphot" || 
						NameOfPlots[WhichPlot] == "h1_Nprot" ||
						NameOfPlots[WhichPlot] == "h1_Npi" ||
						NameOfPlots[WhichPlot] == "h_Etot_subtruct_piplpimi_factor_fracfeed" ||
						NameOfPlots[WhichPlot] == "h_Erec_subtruct_piplpimi_factor_fracfeed"
					   ) { delete pad2; }

					// -----------------------------------------------------------------------------------------------------------------------------------------

					TString ext = "";
					if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

//					PlotCanvas->SaveAs("../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]+"/"+ext+nucleus[WhichNucleus]+"_" 
//						+E[WhichEnergy]+"_" +OutputPlotNames[WhichPlot]+WhatModelsAreIncluded+".pdf");

					//delete PlotCanvas;

					// -----------------------------------------------------------------------------------------------------------------------------------------


				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
