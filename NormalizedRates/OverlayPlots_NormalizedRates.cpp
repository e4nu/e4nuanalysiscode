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

void OverlayPlots_NormalizedRates() {

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
	double EnhaceTail = 1./3.;

	std::vector<TString> xBCut; std::vector<TString> nucleus; std::vector<TString> JustNucleus; std::vector<TString> LabelsOfSamples; 
	std::vector<TString> E; std::vector<double> DoubleE;
	std::vector<TString> LabelE; std::vector<TString> FSIModel; std::vector<TString> DirNames;  std::vector<int> BreakDownColors;
	std::vector<TString> FSILabel; std::vector<TString> NameOfPlots; std::vector<TString> LabelOfPlots;  
	std::vector<TString> OutputPlotNames; std::vector<TH1D*> BreakDownPlots;
	std::vector<int> Colors;
	std::vector<int> Style;

//	nucleus.push_back("4He"); LabelsOfSamples.push_back("^{4}He"); JustNucleus.push_back("He");
	nucleus.push_back("12C"); LabelsOfSamples.push_back("^{12}C"); JustNucleus.push_back("C");
//	nucleus.push_back("56Fe"); LabelsOfSamples.push_back("^{56}Fe");  JustNucleus.push_back("Fe");

	E.push_back("1_161"); LabelE.push_back(" @ E = 1.161 GeV"); DoubleE.push_back(1.161);
//	E.push_back("2_261"); LabelE.push_back(" @ E = 2.261 GeV"); DoubleE.push_back(2.261);	
//	E.push_back("4_461"); LabelE.push_back(" @ E = 4.461 GeV");  DoubleE.push_back(4.461);

	xBCut.push_back("NoxBCut");
//	xBCut.push_back("xBCut");
 
//	Colors.push_back(kBlack); Colors.push_back(kRed); Colors.push_back(kBlue); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);
	Colors.push_back(kBlack); Colors.push_back(kBlack); Colors.push_back(kBlack); Colors.push_back(kMagenta); Colors.push_back(kGreen); Colors.push_back(kOrange + 7);

//	Style.push_back(9); Style.push_back(3); Style.push_back(7); Style.push_back(5);
//	Style.push_back(9); Style.push_back(9); Style.push_back(9); Style.push_back(9); // fancy dashed lines 
	Style.push_back(1); Style.push_back(kDashed); Style.push_back(1); Style.push_back(1);

	BreakDownColors.push_back(kBlue); BreakDownColors.push_back(kCyan); BreakDownColors.push_back(kGreen); BreakDownColors.push_back(kMagenta);

	FSIModel.push_back("Data_Final"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");
	FSIModel.push_back("hA2018_Final_RadCorr_LFGM"); FSILabel.push_back("G2018");  DirNames.push_back("hA2018_Truth_RadCorr");
//	FSIModel.push_back("SuSav2_NoRadCorr_LFGM"); FSILabel.push_back("SuSav2 NoRad");  DirNames.push_back("hA2018_Truth_RadCorr");
	FSIModel.push_back("SuSav2_RadCorr_LFGM"); FSILabel.push_back("SuSav2");  DirNames.push_back("hA2018_Truth_RadCorr");
//	FSIModel.push_back("SuSav2_02_11a_NoRadCorr_LFGM"); FSILabel.push_back("SuSav2");  DirNames.push_back("hA2018_Truth_RadCorr");

//	FSIModel.push_back("Data_Final_NoChargedPions"); FSILabel.push_back("Data"); DirNames.push_back("Data");
//	FSIModel.push_back("hA2018_Final_RadCorr_LFGM_NoChargedPions"); FSILabel.push_back("G2018");  DirNames.push_back("hA2018_Truth_RadCorr");
//	FSIModel.push_back("SuSav2_RadCorr_LFGM_NoChargedPions"); FSILabel.push_back("SuSav2");  DirNames.push_back("hA2018_Truth_RadCorr");

//	FSIModel.push_back("hA2018_Final_NoRadCorr"); FSILabel.push_back("Genie");  DirNames.push_back("hA2018_Truth_NoRadCorr");
//	FSIModel.push_back("hA2018_Truth_NoRadCorr"); FSILabel.push_back("Genie (Truth)");  DirNames.push_back("hA2018_Truth_NoRadCorr");
//	FSIModel.push_back("hN2018_Final_NoRadCorr_LFGM"); FSILabel.push_back("Genie hN2018");  DirNames.push_back("hN2018_Truth_NoRadCorr");
//	FSIModel.push_back("hN2018_Final_RadCorr_LFGM"); FSILabel.push_back("Genie hN2018");  DirNames.push_back("hN2018_Truth_NoRadCorr");

//	NameOfPlots.push_back("MissMomentum"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} P_{T} [GeV/c]"); OutputPlotNames.push_back("MissMomentum");
//	NameOfPlots.push_back("epRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]"); OutputPlotNames.push_back("epRecoEnergy_slice_0");
	NameOfPlots.push_back("eRecoEnergy_slice_0"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_1"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]");  OutputPlotNames.push_back("epRecoEnergy_slice_1");
//	NameOfPlots.push_back("h1_Erec_p_bkgd_slice_sub2p1pi_2p_1"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_1");
//	NameOfPlots.push_back("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_2"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]");  OutputPlotNames.push_back("epRecoEnergy_slice_2");
//	NameOfPlots.push_back("h1_Erec_p_bkgd_slice_sub2p1pi_2p_2"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_2");
//	NameOfPlots.push_back("h1_Etot_p_bkgd_slice_sub2p1pi_1p0pi_3"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{cal} [GeV]");  OutputPlotNames.push_back("epRecoEnergy_slice_3");
//	NameOfPlots.push_back("h1_Erec_p_bkgd_slice_sub2p1pi_2p_3"); LabelOfPlots.push_back("(e,e'p)_{1p0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("eRecoEnergy_slice_3");
//	NameOfPlots.push_back("h_Etot_subtruct_piplpimi_factor_fracfeed"); LabelOfPlots.push_back("E^{cal} Feeddown");  OutputPlotNames.push_back("EcalReso");
//	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_factor_fracfeed"); LabelOfPlots.push_back("E^{QE} Feeddown"); OutputPlotNames.push_back("EQEReso");
	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_3pi"); LabelOfPlots.push_back("(e,e')_{0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("InclusiveeRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_EQE_FullyInclusive"); LabelOfPlots.push_back("(e,e') E^{QE} [GeV]");  OutputPlotNames.push_back("FullyInclusiveeRecoEnergy_slice_0");
	NameOfPlots.push_back("h1_EQE_FullyInclusive_IrregBins"); LabelOfPlots.push_back("(e,e') E^{QE} [GeV]");  OutputPlotNames.push_back("FullyInclusiveeRecoEnergy_slice_0_IrregBins");
	NameOfPlots.push_back("h1_EQE_FullyInclusive_IrregBins_NoPions"); LabelOfPlots.push_back("(e,e') E^{QE} [GeV]");  OutputPlotNames.push_back("FullyInclusiveeRecoEnergy_slice_0_IrregBins_NoPions");
//	NameOfPlots.push_back("h1_E_rec"); LabelOfPlots.push_back("(e,e') E^{QE} [GeV]");  OutputPlotNames.push_back("FullyInclusiveeRecoEnergy_slice_0");
//	NameOfPlots.push_back("h1_E_tot_cut2"); LabelOfPlots.push_back("(e,e')_{0#pi} E^{Cal} Before Subtraction [GeV]");  OutputPlotNames.push_back("h1_E_tot_cut2");

	std::vector<TH1D*> Plots;
	std::vector<TH1D*> Plots_Clones;

	int NxBCuts = xBCut.size();
	int NNuclei = nucleus.size();
	int NEnergies = E.size();
	int NFSIModels = FSIModel.size();
	int NPlots = NameOfPlots.size();

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

					PlotCanvas->SetBottomMargin(0.2);
					PlotCanvas->SetLeftMargin(0.15);

					// ---------------------------------------------------------------------------------------

					Plots.clear();

					TLegend* leg = new TLegend(0.2,0.7,0.5,0.89);
					leg->SetNColumns(1);

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

						// --------------------------------------------------------------------------------------

						// X-axis label

						Plots[WhichFSIModel]->GetXaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetXaxis()->SetLabelSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleSize(TextSize);
						Plots[WhichFSIModel]->GetXaxis()->SetTitleOffset(1.05);

						// X-axis Title

						Plots[WhichFSIModel]->GetXaxis()->SetTitle(JustNucleus[WhichNucleus]+LabelOfPlots[WhichPlot]);

						// Y-axis Title/Tick Size

						Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize); 
						Plots[WhichFSIModel]->GetYaxis()->SetTickSize(0.02);

						// --------------------------------------------------------------------------------------

						// Y-axis label

						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);
						Plots[WhichFSIModel]->GetYaxis()->SetTitle("Weighted Events / GeV");

						// --------------------------------------------------------------------------------------

						Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.65); 

						Plots[WhichFSIModel]->SetLineWidth(LineWidth);

						if (FSILabel[WhichFSIModel] == "Data") { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "lep");}
						else { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l"); }

						// --------------------------------------------------------------------------------------

						// Scaling Factor

//						double ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral(); // default

						double ScalingFactor = 1.;

//						double ScalingFactor = 1. / Plots[WhichFSIModel]->Integral(); // area normalized

//						if (
//							NameOfPlots[WhichPlot] == "h1_Nphot" || 
//							NameOfPlots[WhichPlot] == "h1_Nprot" ||
//							NameOfPlots[WhichPlot] == "h1_Npi") { ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral();}

//						else { ScalingFactor = 1. / Plots[WhichFSIModel]->Integral(); }  // area normalized

//						double ScalingFactor = 1. / Plots[WhichFSIModel]->GetMaximum(); // peak at 1

						ScalingFactor = Plots[0]->GetEntries() / Plots[WhichFSIModel]->GetEntries();
//						double ScalingFactor = Plots[0]->GetMaximum() / Plots[WhichFSIModel]->GetMaximum();
//						double ScalingFactor = 18E8 / Plots[WhichFSIModel]->GetMaximum();
//						double ScalingFactor = 1.;

//						ScalingFactor = Plots[0]->Integral() / Plots[WhichFSIModel]->Integral();
						Plots[WhichFSIModel]->Scale(ScalingFactor);

						// -----------------------------------------------------------------------------------

						// Accounting for the fact that the bin width might not be constant

						ReweightPlots(Plots[WhichFSIModel]);

						// --------------------------------------------------------------------------------------

						// Rebining & ranges

//						if (string(OutputPlotNames[WhichPlot]).find("epRecoEnergy_slice") != std::string::npos) 
//							{ for (int i = 0; i < NECalRebin; i++) { Plots[WhichFSIModel]->Rebin(); } }

						if ( NameOfPlots[WhichPlot] == "h1_EQE_FullyInclusive" && DoubleE[WhichEnergy] == 1.161) 
							{ for (int i = 0; i < 4; i++) { Plots[WhichFSIModel]->Rebin(); Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.4,1.7); } }

						if ( NameOfPlots[WhichPlot] == "h1_EQE_FullyInclusive" && DoubleE[WhichEnergy] == 2.261) 
							{ for (int i = 0; i < 4; i++) { Plots[WhichFSIModel]->Rebin(); Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.4,3.); } }

						if ( NameOfPlots[WhichPlot] == "h1_EQE_FullyInclusive" && DoubleE[WhichEnergy] == 4.461) 
							{ for (int i = 0; i < 5; i++) { Plots[WhichFSIModel]->Rebin(); Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(1.5,6.); } }

						if (string(OutputPlotNames[WhichPlot]).find("epRecoEnergy_slice") != std::string::npos && nucleus[WhichNucleus] == "12C" 
							&& DoubleE[WhichEnergy] == 2.261) 
							{ Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.6,2.5); }

						if (string(OutputPlotNames[WhichPlot]).find("epRecoEnergy_slice") != std::string::npos && nucleus[WhichNucleus] == "56Fe" 
							&& DoubleE[WhichEnergy] == 4.461) 
							{ Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(1.,5.); }

						if (string(OutputPlotNames[WhichPlot]).find("eRecoEnergy_slice") != std::string::npos) 
							{ for (int i = 0; i < 0; i++) { Plots[WhichFSIModel]->Rebin(); } }

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

						// ---------------------------------------------------------------------------------------------------

						// Max, min, title & # divisions

						double localmax = Plots[WhichFSIModel]->GetMaximum();
						if (localmax > max) { max = localmax; }
						double height = 1.05;
						if ( xBCut[WhichxBCut] == "xBCut" ) { height = 1.1; }
						Plots[0]->GetYaxis()->SetRangeUser(0.,height*max);

						double localmin = Plots[WhichFSIModel]->GetBinContent(Plots[WhichFSIModel]->FindBin(4)); // multiplicity 4 is the highest one in data
						if (localmin < min && localmin != 0) { min = localmin; }

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);
						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);

						// --------------------------------------------------------------------------------------------------

						// Multiplicity plots

						if (NameOfPlots[WhichPlot] == "h1_Nphot" || NameOfPlots[WhichPlot] == "h1_Nprot" || NameOfPlots[WhichPlot] == "h1_Npi") {

							Plots[WhichFSIModel]->GetYaxis()->SetLabelOffset(-0.004);
							Plots[WhichFSIModel]->Rebin();
							Plots[0]->GetYaxis()->SetRangeUser(0.5*min,2.*max);
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
								Plots[WhichFSIModel]->SetLineStyle(Style[WhichFSIModel]); 
								Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines
								Plots[0]->Draw("e same"); 
							}
						}

		                                // ---------------------------------------------------------------------------------------------------
		                                // --------------------------------------------------------------------------------------------------

					} // End of the loop over the FSI Models 

					leg->SetBorderSize(0);
					leg->SetTextFont(FontStyle);
					leg->SetTextSize(TextSize);
					leg->Draw(); // Just data + e.g. susav2

					// -----------------------------------------------------------------------------------------------------------------------------------------

					TString ext = "";
					if ( xBCut[WhichxBCut] == "xBCut" ) { ext = "xB_"; } 

//					PlotCanvas->SaveAs("../myPlots/pdf/"+xBCut[WhichxBCut]+"/"+version+nucleus[WhichNucleus]+"/"+E[WhichEnergy]+"/"+ext+nucleus[WhichNucleus]+"_" 
//						+E[WhichEnergy]+"_" +OutputPlotNames[WhichPlot]+".pdf");

					//delete PlotCanvas;

					// -----------------------------------------------------------------------------------------------------------------------------------------


				} // End of the loop over the plots

			} // End of the loop over the nuclei

		} // End of the loop over the energies

	} // End of the loop over the xB kinematics

} // End of the program
