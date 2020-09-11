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
//#include "/home/afroditi/Dropbox/PhD/Secondary_Code/ToString.cpp"
#include "/home/afroditi/Dropbox/PhD/Secondary_Code/myFunctions.cpp"

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

	NameOfPlots.push_back("h_Erec_subtruct_piplpimi_noprot_3pi"); LabelOfPlots.push_back("(e,e')_{0#pi} E^{QE} [GeV]");  OutputPlotNames.push_back("InclusiveeRecoEnergy_slice_0");

//	NameOfPlots.push_back("h1_EQE_FullyInclusive"); LabelOfPlots.push_back("(e,e') E^{QE} [GeV]");  OutputPlotNames.push_back("FullyInclusiveeRecoEnergy_slice_0");

	NameOfPlots.push_back("h1_EQE_FullyInclusive_IrregBins"); LabelOfPlots.push_back("(e,e') E^{QE} [GeV]");  OutputPlotNames.push_back("FullyInclusiveeRecoEnergy_slice_0_IrregBins");

//	NameOfPlots.push_back("h1_EQE_FullyInclusive_IrregBins_NoPions"); LabelOfPlots.push_back("(e,e') E^{QE} [GeV]");  OutputPlotNames.push_back("FullyInclusiveeRecoEnergy_slice_0_IrregBins_NoPions");
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

	// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// For now, numbers specific to carbon

	std::map<TString,double> IntegratedCharge; // mC

	IntegratedCharge["1_161"]=0.19;
	IntegratedCharge["2_261"]=1.79;
	IntegratedCharge["4_461"]=2.14;

	std::map<TString,double> TargetLength; // cm

	TargetLength["1_161"]=0.1; 
	TargetLength["2_261"]=0.1;
	TargetLength["4_461"]=0.1;

	// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

					double DataIntegral = 1.;

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

						Plots[WhichFSIModel]->GetYaxis()->SetTitleSize(TextSize-0.01); 
						Plots[WhichFSIModel]->GetYaxis()->SetTickSize(0.02);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelSize(TextSize);
						Plots[WhichFSIModel]->GetYaxis()->SetTitle("Weighted Events / (GeV mC cm)");

						// --------------------------------------------------------------------------------------

						Plots[WhichFSIModel]->GetYaxis()->SetTitleFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetLabelFont(FontStyle);
						Plots[WhichFSIModel]->GetYaxis()->SetTitleOffset(0.85); 

						Plots[WhichFSIModel]->SetLineWidth(LineWidth);

						if (FSILabel[WhichFSIModel] == "Data") { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "lep");}
						else { leg->AddEntry(Plots[WhichFSIModel],FSILabel[WhichFSIModel], "l"); }

						// --------------------------------------------------------------------------------------

						// Scaling Factor
						// Scale to data integral after dividing by integrated charge & length

						if (FSILabel[WhichFSIModel] == "Data") { 
							Plots[WhichFSIModel]->Scale(1. / (IntegratedCharge[E[WhichEnergy]] * TargetLength[E[WhichEnergy]]) );
							DataIntegral =  Plots[0]->Integral();
						}

						double ScalingFactor = DataIntegral / Plots[WhichFSIModel]->Integral();
						Plots[WhichFSIModel]->Scale(ScalingFactor);

						// -----------------------------------------------------------------------------------

						// Accounting for the fact that the bin width might not be constant

						ReweightPlots(Plots[WhichFSIModel]);

						// --------------------------------------------------------------------------------------

						// Rebining & ranges


//						if (DoubleE[WhichEnergy] == 1.161)  
//							{ for (int i = 0; i < 4; i++) { Plots[WhichFSIModel]->Rebin(); Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.4,1.7); } }

//						if ( DoubleE[WhichEnergy] == 2.261) 
//							{ for (int i = 0; i < 4; i++) { Plots[WhichFSIModel]->Rebin(); Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(0.4,3.); } }

//						if ( DoubleE[WhichEnergy] == 4.461) 
//							{ for (int i = 0; i < 4; i++) { Plots[WhichFSIModel]->Rebin(); Plots[WhichFSIModel]->GetXaxis()->SetRangeUser(1.5,6.); } }

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

						TString XLabel = Plots[WhichFSIModel]->GetXaxis()->GetTitle();
						Plots[0]->GetXaxis()->SetTitle(XLabel);

						Plots[WhichFSIModel]->GetXaxis()->SetNdivisions(Ndivisions);
						Plots[WhichFSIModel]->GetYaxis()->SetNdivisions(Ndivisions);

						// --------------------------------------------------------------------------------------------------

						if (FSILabel[WhichFSIModel] == "Data") { 

							Plots[WhichFSIModel]->SetMarkerStyle(20); 
							Plots[WhichFSIModel]->SetMarkerSize(2.); 
							Plots[WhichFSIModel]->SetMarkerColor(kBlack); 

							gStyle->SetErrorX(0); // Removing the horizontal errors
							Plots[WhichFSIModel]->Draw("e same"); 

						} else { 

//							Plots[WhichFSIModel]->Draw("hist same"); // draw them as histos
							Plots[WhichFSIModel]->SetLineStyle(Style[WhichFSIModel]); 
							Plots[WhichFSIModel]->Draw("C hist same");  // draw them as lines
							Plots[0]->Draw("e same"); 
						
							double acc = 3.;

							if (FSILabel[WhichFSIModel] == "G2018") {
							
								TLatex* latex = new TLatex();
								latex->SetTextFont(FontStyle);
								latex->SetTextSize(TextSize);
								latex->DrawLatexNDC(0.17,0.6,"SF G2018 = "+ToString(round(ScalingFactor,acc)));

							}

							if (FSILabel[WhichFSIModel] == "SuSav2") {
							
								TLatex* latex = new TLatex();
								latex->SetTextFont(FontStyle);
								latex->SetTextSize(TextSize);
								latex->DrawLatexNDC(0.17,0.5,"SF SuSav2 = "+ToString(round(ScalingFactor,acc)));

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
