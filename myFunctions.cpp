#include "TMath.h"
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGaxis.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

#include "AfroConstants.h"

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------

double round(double var,int acc = 0) 
{ 
    double value = (int)(var * TMath::Power(10.,acc) + .5); 
    return (double)value / TMath::Power(10.,acc); 
} 

// -------------------------------------------------------------------------------------------------------------------------------------

double IntegratedXSec(TH1D* h, int MinBin = -1, int MaxBin = -1) {

	int NBinsX = h->GetXaxis()->GetNbins();

	double IntegratedXSec = 0;

	double LocalMinBin = 0;
	if (MinBin != -1) { LocalMinBin = MinBin; }

	double LocalMaxBin = NBinsX;
	if (MaxBin != -1) { LocalMaxBin = MaxBin; }

	for (int WhichXBin = LocalMinBin; WhichXBin < LocalMaxBin; WhichXBin++) {

		double BinWidth = h->GetBinWidth(WhichXBin+1);
		double BinEntry = h->GetBinContent(WhichXBin+1);

		IntegratedXSec += BinEntry * BinWidth;

	}

	return IntegratedXSec;

}

// -------------------------------------------------------------------------------------------------------------------------------------

void GlobalSettings() {

	TGaxis::SetMaxDigits(5);

	gStyle->SetTitleSize(TextSize,"t"); 
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetOptStat(0);	

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();

}

// -------------------------------------------------------------------------------------------------------------------------------------

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

void ApplyAcceptanceCorrUnc(TH1D* h) {

	double SystUnc = AcceptanceCorrUnc;

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 

		double error = h->GetBinError(i);
		double content = h->GetBinContent(i);
		double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(SystUnc*content,2.));
		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

void ApplyOverallNormUnc(TH1D* h) {

	double SystUnc = OverallNormUnc;

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 

		double error = h->GetBinError(i);
		double content = h->GetBinContent(i);
		double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(SystUnc*content,2.));
		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

void ApplySystUnc(TH1D* h, TString Energy) {

	double SystUnc = 0;
	if ( Energy == "1_161" ) { SystUnc = SystUnc1GeV; }
	if ( Energy == "2_261" ) { SystUnc = SystUnc2GeV; }
	if ( Energy == "4_461" ) { SystUnc = SystUnc4GeV; }

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 
					
//		double error = h->GetBinError(i);
//		double newerror = error * (1. + systunc);

		double error = h->GetBinError(i);
		double content = h->GetBinContent(i);
		double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(SystUnc*content,2.));
		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

void ApplySectorSystUnc(TH1D* h, TString Energy) {

	double SystUnc = 0;
	if ( Energy == "1_161" ) { SystUnc = SectorSystUnc1GeV; }
	if ( Energy == "2_261" ) { SystUnc = SectorSystUnc2GeV; }
	if ( Energy == "4_461" ) { SystUnc = SectorSystUnc4GeV; }

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 
					
//		double error = h->GetBinError(i);
//		double newerror = error * (1. + systunc);

		double error = h->GetBinError(i);
		double content = h->GetBinContent(i);
		double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(SystUnc*content,2.));
		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

std::vector<double> GetUncertaintyBand(std::vector<TH1D*> h) {

	int NPlots = h.size();

	int NBins = h[0]->GetNbinsX();

	std::vector<double> UncVector; UncVector.clear();

	for (int i = 1; i <= NBins; i++) {  

		std::vector<double> BinVector; BinVector.clear();

		for (int WhichPlot = 0; WhichPlot < NPlots; WhichPlot ++) {

			double content = h[WhichPlot]->GetBinContent(i);
			BinVector.push_back(content);

		} // End of the loop over the plots

		auto min = std::min_element(std::begin(BinVector), std::end(BinVector));

		auto max = std::max_element(std::begin(BinVector), std::end(BinVector));

		double spread = *max - *min;

		UncVector.push_back(0.5*TMath::Abs(spread)); // take half the spread as an uncertainty

	} // End of the loop over the bins

	return UncVector;

}

// ----------------------------------------------------------------------------------------------------------------

void ApplySectorSystUnc(TH1D* h, std::vector<double> sectorsystunc) {

	double NBins = h->GetNbinsX(); 
				
	for (int i = 1; i <= NBins; i++) { 
					
//		double error = h->GetBinError(i);
//		double newerror = error * (1. + systunc);

		double error = h->GetBinError(i);
		double sectorerorr = sectorsystunc.at(i-1);
		double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(sectorerorr,2.));

		h->SetBinError(i,newerror);

	}

}

// ----------------------------------------------------------------------------------------------------------------

TH1D* VectorToHistSystUnc(TH1D* h, std::vector<double> sectorsystunc, TString name) {

	double NBins = h->GetXaxis()->GetNbins(); 

	TH1D* ClonePlot = (TH1D*)(h->Clone("Unc_"+name));
				
	for (int i = 1; i <= NBins; i++) { 

		double unc = sectorsystunc[i-1];

		ClonePlot->SetBinContent(i,unc);
		ClonePlot->SetBinError(i,0);

	}

	return ClonePlot;

}

// -------------------------------------------------------------------------------------------------------------------------------------

TString ToStringInt(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

// -------------------------------------------------------------------------------------------------------------------------------------

TString ToStringDouble(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

// -------------------------------------------------------------------------------------------------------------------------------------

void PrettyDoubleXSecPlot(TH1D* h) {

	// ----------------------------------------------------------------------------------------------------------------

	h->SetLineWidth(LineWidth);

	// ----------------------------------------------------------------------------------------------------------------

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);

	// ----------------------------------------------------------------------------------------------------------------

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(TextSize); 
	//h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTitle(DoubleXSecTitle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(1.05);
	h->GetYaxis()->SetNdivisions(Ndivisions);

	return;	

}

// -------------------------------------------------------------------------------------------------------------------------------------

void PrettyGraph(TGraph* h) {

	// ----------------------------------------------------------------------------------------------------------------

	h->SetLineWidth(LineWidth);

	// ----------------------------------------------------------------------------------------------------------------

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);

	// ----------------------------------------------------------------------------------------------------------------

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(TextSize); 
	//h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTitle(DoubleXSecTitle);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(1.05);
	h->GetYaxis()->SetNdivisions(Ndivisions);

	return;	

}

// -------------------------------------------------------------------------------------------------------------------------------------

void AbsoluteXSecScaling(TH1D* h, TString Sample, TString Nucleus, TString E) {  

	double SF = 1.;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	// Data sets

	if (Sample == "Data") { 

		SF = 1. / ( IntegratedCharge[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus] ) * ConversionFactorCm2ToMicroBarn;

	}

	else if (Sample == "Data_FilterRuns") { 

		SF = 1. / (IntegratedCharge_FilterRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn;

	}

	else if (Sample == "Data_NewFilterRuns") { 

		SF = 1. / (IntegratedCharge_NewFilterRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn;

	}

	else if (Sample == "GoodRunList_Data") { 

		SF = 1. / (IntegratedCharge_GoodRunList_AllRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn;

	}

	else if (Sample == "LowCurrent_GoodRunList_Data") { 

		SF = 1. / (IntegratedCharge_GoodRunList_LowCurrentRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn;


	}

	else if (Sample == "HighCurrent_GoodRunList_Data") { 

		SF = 1. / (IntegratedCharge_GoodRunList_HighCurrentRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn;

	}

	else if (Sample == "Pinned Data" || Sample == "Pinned Data No Rotations") { 

		SF = 1. / (IntegratedCharge_PinnedFiles[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn;

	}

	else if (Sample == "Mikhail Data") { 

		SF = 1. / (IntegratedCharge_MikhailFiles[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn;

	}

	// -------------------------------------------------------------------------------------------------------------------------------------------------

	// Simulation sets

	else if (Sample == "SuSav2" || Sample == "SuSav2_NoAccMaps" 
		/*|| Sample == "SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc_NoThetaCut" || Sample == "hA2018_Final_RadCorr_LFGM_Truth_WithoutFidAcc_NoThetaCut"*/) { 

				SF = (SuSav2GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (SuSav2NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "SuSav2 NoRad") { 

				SF = (SuSav2GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (NoRadSuSav2NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "SuSav2 Master NoRad") { 

				SF = (SuSav2GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (MasterNoRadSuSav2NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "SuSav2 Master Rad") { 

				SF = (SuSav2GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (MasterRadSuSav2NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "SuSav2 Rad Schwinger") { 

				SF = (SuSav2GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (RadSchwingerSuSav2NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "SuSav2 Rad Updated Schwinger") { 

				SF = (SuSav2GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (RadUpdatedSchwingerSuSav2NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "G2018 Rad Updated Schwinger") { 

				SF = (G2018GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (RadUpdatedSchwingerG2018NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "G2018 NoRad") { 

				SF = (G2018GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (NoRadG2018NumberEvents[std::make_pair(Nucleus, E)] ) ) ;

	}

	else if (Sample == "G2018") { 

		SF = ( G2018GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (G2018NumberEvents[std::make_pair(Nucleus, E)] ) );

	}

	else if (Sample == "G2018 Master Rad") { 

		SF = ( G2018GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (MasterRadG2018NumberEvents[std::make_pair(Nucleus, E)] ) );

	}

	else if (Sample == "G2018 Master NoRad") { 

		SF = ( G2018GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (MasterNoRadG2018NumberEvents[std::make_pair(Nucleus, E)] ) );

	}

	else if (Sample == "G2018 QE Only") { 

		SF = ( QEG2018GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (QEMasterRadG2018NumberEvents[std::make_pair(Nucleus, E)] ) );

	}

//	else if (Sample == "G18_02c NoRad") { 

//		SF = ( G18_02cGenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
//					ConversionFactorCm2ToMicroBarn / (NoRadG18_02cNumberEvents[std::make_pair(Nucleus, E)] ) );

//	}

//	else if (Sample == "G18_02d NoRad") { 

//		SF = ( G18_02dGenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
//					ConversionFactorCm2ToMicroBarn / (NoRadG18_02dNumberEvents[std::make_pair(Nucleus, E)] ) );

//	}

	else {

		std::cout << "Craaaaaaaaaaaaaaap !!!!!!!!! What is the SF in AbsoluteXSecScaling for " << h->GetName() << " in " << Sample << "???????????????" << std::endl;

	}		

	h->Scale(SF);

}

// -------------------------------------------------------------------------------------------------------------------------------------

void AbsoluteXSec2DScaling(TH2D* h, TString Sample, TString Nucleus, TString E) {  

	double SF = 1.;

	// ----------------------------------------------------------------------------------------------------------------------------------------

	// Data sets

	if (Sample == "Data") { 

		SF = 1. / ( IntegratedCharge[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus] ) * ConversionFactorCm2ToMicroBarn / dOmega;

	}

	if (Sample == "Data_FilterRuns") { 

		SF = 1. / (IntegratedCharge_FilterRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn / dOmega;

	}

	if (Sample == "Data_NewFilterRuns") { 

		SF = 1. / (IntegratedCharge_NewFilterRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn / dOmega;

	}

	if (Sample == "GoodRunList_Data") { 

		SF = 1. / (IntegratedCharge_GoodRunList_AllRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn / dOmega;

	}

	if (Sample == "LowCurrent_GoodRunList_Data") { 

		SF = 1. / (IntegratedCharge_GoodRunList_LowCurrentRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn / dOmega;


	}

	if (Sample == "HighCurrent_GoodRunList_Data") { 

		SF = 1. / (IntegratedCharge_GoodRunList_HighCurrentRuns[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn / dOmega;

	}

	if (Sample == "Pinned Data" || Sample == "Pinned Data No Rotations") { 

		SF = 1. / (IntegratedCharge_PinnedFiles[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn / dOmega;

	}

	if (Sample == "Mikhail Data") { 

		SF = 1. / (IntegratedCharge_MikhailFiles[std::make_pair(Nucleus, E)] *\
						    TargetLength[std::make_pair(Nucleus, E)] *\
						    TargetDensity[std::make_pair(Nucleus, E)] *\
						    OverallUnitConversionFactor / MassNumber[Nucleus]) * ConversionFactorCm2ToMicroBarn / dOmega;

	}

	// -------------------------------------------------------------------------------------------------------------------------------------------------

	// Simulation sets

	if (Sample == "SuSav2" || Sample == "SuSav2_NoAccMaps" ) { 

				SF = (SuSav2GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (SuSav2NumberEvents[std::make_pair(Nucleus, E)] *\
					dOmega) ) ;

	}

	if (Sample == "G2018") { 

		SF = ( G2018GenieXSec[std::make_pair(Nucleus, E)] * TMath::Power(10.,-38.) *\
					ConversionFactorCm2ToMicroBarn / (G2018NumberEvents[std::make_pair(Nucleus, E)] *\
					dOmega) );

	}	

	h->Scale(SF);

}

// -------------------------------------------------------------------------------------------------------------------------------------

void ApplyRebinningTProfile(TProfile* h, TString Energy, TString PlotVar) {

	// -----------------------------------------------------------------------------------------------------------------------------

	if (string(PlotVar).find("Omega") != std::string::npos) {

		if (Energy == "1_161") { for (int i = 0; i < 5; i++) { h->Rebin(); } }
		if (Energy == "2_261") { for (int i = 0; i < 5; i++) { h->Rebin(); } }
		if (Energy == "4.461") { for (int i = 0; i < 6; i++) { h->Rebin(); } }

	} else if (string(PlotVar).find("EcalReso") != std::string::npos || string(PlotVar).find("ECalReso") != std::string::npos || string(PlotVar).find("h_Etot_subtruct_piplpimi_factor_fracfeed") != std::string::npos ) {

	} else if (
		string(PlotVar).find("T2KEQEReso") != std::string::npos
	) {

		for (int i = 0; i < 2; i++) { h->Rebin();} 

	} else if (
		string(PlotVar).find("EQEReso") != std::string::npos || 
		string(PlotVar).find("h_Erec_subtruct_piplpimi_factor_fracfeed") != std::string::npos ||
		string(PlotVar).find("h_Erec_subtruct_piplpimi_noprot_frac_feed") != std::string::npos
	) {

	} else if (string(PlotVar).find("EQE") != std::string::npos || string(PlotVar).find("eReco") != std::string::npos || string(PlotVar).find("Erec") != std::string::npos) {

	} else if (string(PlotVar).find("cal") != std::string::npos || string(PlotVar).find("Cal") != std::string::npos || string(PlotVar).find("epReco") != std::string::npos || string(PlotVar).find("Etot") != std::string::npos || string(PlotVar).find("E_tot") != std::string::npos) {

	} else if (string(PlotVar).find("PT") != std::string::npos || string(PlotVar).find("MissMomentum") != std::string::npos) {

		for (int i = 0; i < 2; i++) { h->Rebin();} 

	} else if (string(PlotVar).find("DeltaAlphaT") != std::string::npos ) {

		for (int i = 0; i < 1; i++) { h->Rebin();} 

	} else if (string(PlotVar).find("DeltaPhiT") != std::string::npos ) {

		for (int i = 0; i < 1; i++) { h->Rebin();} 

	} else if (string(PlotVar).find("Wvar") != std::string::npos || string(PlotVar).find("W_") != std::string::npos ) {

		for (int i = 0; i < 1; i++) { h->Rebin();} 

	} else if (string(PlotVar).find("T2KEQEReso") != std::string::npos || string(PlotVar).find("W_") != std::string::npos ) {

		for (int i = 0; i < 1; i++) { h->Rebin();} 

	}

	else { cout << "Aaaaaaaaaaaah ! How do I rebin this plot in ApplyRebinningTProfile ?" << endl; }

	return;	

}

// -------------------------------------------------------------------------------------------------------------------------------------

void ApplyRebinning(TH1D* h, TString Energy, TString PlotVar) {

	// -----------------------------------------------------------------------------------------------------------------------------

	if (string(PlotVar).find("Omega_FullyInclusive") != std::string::npos) {

		if (Energy == "1_161") { for (int i = 0; i < 5; i++) { h->Rebin(); } }

	} else if (string(PlotVar).find("Omega") != std::string::npos) {

		if (Energy == "1_161") { for (int i = 0; i < 5; i++) { h->Rebin(); } }
		if (Energy == "2_261") { for (int i = 0; i < 5; i++) { h->Rebin(); } }
		if (Energy == "4.461") { for (int i = 0; i < 6; i++) { h->Rebin(); } }

	} else if (string(PlotVar).find("EcalReso") != std::string::npos || string(PlotVar).find("ECalReso") != std::string::npos || string(PlotVar).find("h_Etot_subtruct_piplpimi_factor_fracfeed") != std::string::npos ) {

	} else if (
		string(PlotVar).find("T2KEQEReso") != std::string::npos
	) {

		for (int i = 0; i < 2; i++) { h->Rebin();}

	} else if (
		string(PlotVar).find("EQEReso") != std::string::npos || 
		string(PlotVar).find("h_Erec_subtruct_piplpimi_factor_fracfeed") != std::string::npos ||
		string(PlotVar).find("h_Erec_subtruct_piplpimi_noprot_frac_feed") != std::string::npos

	) {

	} else if (string(PlotVar).find("EQE") != std::string::npos || string(PlotVar).find("eReco") != std::string::npos || string(PlotVar).find("Erec") != std::string::npos) {

	} else if (string(PlotVar).find("cal") != std::string::npos || string(PlotVar).find("Cal") != std::string::npos || string(PlotVar).find("epReco") != std::string::npos || string(PlotVar).find("Etot") != std::string::npos || string(PlotVar).find("E_tot") != std::string::npos) {

	} else if (string(PlotVar).find("PT") != std::string::npos || string(PlotVar).find("MissMomentum") != std::string::npos) {

		for (int i = 0; i < 2; i++) { h->Rebin();} 

	} else if (string(PlotVar).find("DeltaAlphaT") != std::string::npos ) {

		for (int i = 0; i < 1; i++) { h->Rebin();} 

	} else if (string(PlotVar).find("DeltaPhiT") != std::string::npos ) {

		for (int i = 0; i < 1; i++) { h->Rebin();} 

	} else if (string(PlotVar).find("Wvar") != std::string::npos || string(PlotVar).find("W_") != std::string::npos ) {

		for (int i = 0; i < 2; i++) { h->Rebin();} 

	}

	else { cout << "Aaaaaaaaaaaah ! How do I rebin this plot ?" << endl; }

	return;	

}

// -------------------------------------------------------------------------------------------------------------------------------------

void ApplyRange(TH1D* h, TString Energy, TString PlotVar) {

	// -----------------------------------------------------------------------------------------------------------------------------

	if (string(PlotVar).find("Omega") != std::string::npos) {

		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.,0.7); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(0.,1.5); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(0.5,3.); }

	} else if (
		string(PlotVar).find("EcalReso") != std::string::npos || string(PlotVar).find("ECalReso") != std::string::npos || 
		string(PlotVar).find("h_Etot_subtruct_piplpimi_factor_fracfeed") != std::string::npos ||
		string(PlotVar).find("h1_Ecal_Reso") != std::string::npos ||
		string(PlotVar).find("h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed") != std::string::npos

	) {

		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(-0.7,0.06); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(-0.7,0.06); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(-0.7,0.03); }

	} else if (
		string(PlotVar).find("T2KEQEReso") != std::string::npos
		) {

		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(-0.75,0.39); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(-0.69,0.21); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(-0.75,0.21); }


	} else if (
		string(PlotVar).find("EQEReso") != std::string::npos || 
		string(PlotVar).find("h_Erec_subtruct_piplpimi_factor_fracfeed") != std::string::npos ||
		string(PlotVar).find("h_Erec_subtruct_piplpimi_factor_fracfeed") != std::string::npos ||
		string(PlotVar).find("h_Erec_subtruct_piplpimi_noprot_frac_feed") != std::string::npos
		) {

		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(-0.75,0.21); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(-0.69,0.21); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(-0.75,0.21); }

	} else if (string(PlotVar).find("EQE") != std::string::npos || string(PlotVar).find("eReco") != std::string::npos || string(PlotVar).find("Erec") != std::string::npos) {

		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.47,1.4); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(0.7,2.6); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(1.9,5.2); }

	} else if (string(PlotVar).find("Etot") != std::string::npos || string(PlotVar).find("Cal") != std::string::npos || string(PlotVar).find("cal") != std::string::npos || string(PlotVar).find("epReco") != std::string::npos || string(PlotVar).find("E_tot") != std::string::npos || string(PlotVar).find("h1_Ecal_SuperFine") != std::string::npos) {

//		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.57,1.23); } // default, but now in the Ecal 6-pannel, we need to expand the range for smaller box 
		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.5,1.23); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(0.67,2.4); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(1.5,4.6); }

	} else if (string(PlotVar).find("PT") != std::string::npos || string(PlotVar).find("MissMomentum") != std::string::npos) {

	} else if (string(PlotVar).find("DeltaAlphaT") != std::string::npos ) {

	} else if (string(PlotVar).find("DeltaPhiT") != std::string::npos ) {

	} else if (string(PlotVar).find("Wvar") != std::string::npos || string(PlotVar).find("W_") != std::string::npos ) {

		h->GetXaxis()->SetRangeUser(0.6,1.5);

	}


	else { cout << "Aaaaaaaaaaaah ! How do I set the range for this plot ?" << endl; }

	return;	

}

// -----------------------------------------------------------------------------------------------------------------------------

void ApplyRange(TGraph* h, TString Energy, TString PlotVar) {

	// -----------------------------------------------------------------------------------------------------------------------------

	if (string(PlotVar).find("Omega") != std::string::npos) {

		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.,0.7); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(0.,1.5); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(0.5,3.); }

	} else if (string(PlotVar).find("EcalReso") != std::string::npos || string(PlotVar).find("ECalReso") != std::string::npos || string(PlotVar).find("h_Etot_subtruct_piplpimi_factor_fracfeed") != std::string::npos ) {

		h->GetXaxis()->SetRangeUser(-0.81,0.07);

	} else if (
		string(PlotVar).find("EQEReso") != std::string::npos || 
		string(PlotVar).find("h_Erec_subtruct_piplpimi_factor_fracfeed") != std::string::npos ||
		string(PlotVar).find("h_Erec_subtruct_piplpimi_noprot_frac_feed") != std::string::npos
		) {

		h->GetXaxis()->SetRangeUser(-0.85,0.2);

	} else if (string(PlotVar).find("EQE") != std::string::npos || string(PlotVar).find("eReco") != std::string::npos || string(PlotVar).find("Erec") != std::string::npos) {

		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.47,1.4); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(0.7,2.6); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(2.,5.); }

	} else if (string(PlotVar).find("Etot") != std::string::npos || string(PlotVar).find("Cal") != std::string::npos || string(PlotVar).find("cal") != std::string::npos || string(PlotVar).find("epReco") != std::string::npos || string(PlotVar).find("E_tot") != std::string::npos || string(PlotVar).find("h1_Ecal_SuperFine") != std::string::npos) {

//		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.57,1.23); } // default, but now in the Ecal 6-pannel, we need to expand the range for smaller box 
		if (Energy == "1_161") { h->GetXaxis()->SetRangeUser(0.5,1.23); }
		if (Energy == "2_261") { h->GetXaxis()->SetRangeUser(0.67,2.4); }
		if (Energy == "4_461") { h->GetXaxis()->SetRangeUser(1.5,4.6); }

	} else if (string(PlotVar).find("PT") != std::string::npos || string(PlotVar).find("MissMomentum") != std::string::npos) {

	} else if (string(PlotVar).find("DeltaAlphaT") != std::string::npos ) {

	} else if (string(PlotVar).find("DeltaPhiT") != std::string::npos ) {

	} else if (string(PlotVar).find("Wvar") != std::string::npos || string(PlotVar).find("W_") != std::string::npos ) {

		h->GetXaxis()->SetRangeUser(0.6,1.5);

	}


	else { cout << "Aaaaaaaaaaaah ! How do I set the range for this plot ?" << endl; }

	return;	

}
// -------------------------------------------------------------------------------------------------------------------------------------

void UniversalE4vFunction(TH1D* h, TString DataSetLabel, TString nucleus, TString E, TString name) {

	// Scale to obtain absolute double differential cross sections 
	// Use charge, density and length for data samples
	// Use total number of events in genie sample and relevant genie cross sections for simulation

	AbsoluteXSecScaling(h,DataSetLabel,nucleus,E);

	// Rebin if necessary

	ApplyRebinning(h,E,name);

	// Division by the bin width

	ReweightPlots(h);

	// Use relevant ranges
			
	ApplyRange(h,E,name);

	// if data sample: 

	// 	apply systematics due to rotations et al

	if (string(DataSetLabel).find("Data") != std::string::npos) { ApplySystUnc(h, E); }

	// 	apply acceptance systematics using sector-by -sector uncertainties

	if (string(DataSetLabel).find("Data") != std::string::npos) { ApplySectorSystUnc(h, E); }

	//	3% overall normalization uncertainty 
	//	2% electron efficiency
	// 	not included in plots, will be mentioned in text 

	//if (string(DataSetLabel).find("Data") != std::string::npos) { ApplyOverallNormUnc(h); }

	//	acceptance constant correction uncertainty
	//if (string(DataSetLabel).find("Data") != std::string::npos) { ApplyAcceptanceCorrUnc(h); }
	
}

// -------------------------------------------------------------------------------------------------------------------------------------

TH1D* AcceptanceCorrection(TH1D* h, TString ScaleToDataSet, TString nucleus, TString E, TString name, TString xBCut, TString Extension = "") {

	TH1D::SetDefaultSumw2();

	std::vector<TH1D*> Plots; Plots.clear();
	std::vector<TString> FSIModel; FSIModel.clear();

	// --------------------------------------------------------------------------------------	

	TString AlternativeModel = "hA2018_Final";
	if (ScaleToDataSet == "hA2018_Final") { AlternativeModel = "SuSav2"; }

	// Unfolding using SuSav2
	// keep in mind that the Rad G2018 sample is questionable
	
	FSIModel.push_back("SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger"+Extension); // 0: SuSav2 Rad for radiation correction
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc"+Extension); // 1: Reco 1p0pi SuSav2 NoRad plot for average & for radiation correction
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc_Offset"+Extension); // 2: Reco 1p0pi G2018 NoRad Offset plot for average
	FSIModel.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc"+Extension); // 3: True 1p0pi SuSav2 NoRad plot for average
	FSIModel.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Offset"+Extension); // 4: True 1p0pi G2018 NoRad plot Offset for average

	if (
		string(name).find("T2KEQEReso") != std::string::npos ||
		name == "h_Erec_subtruct_piplpimi_noprot_3pi" || 
		name == "h_Erec_subtruct_piplpimi_noprot_frac_feed" ||
		name == "h_Erec_subtruct_piplpimi_noprot_frac_feed3pi"
	) {

		FSIModel[0] = "SuSav2_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger"+Extension;
		FSIModel[1] = "SuSav2_NoRadCorr_LFGM_Truth0pi_WithFidAcc"+Extension;
		FSIModel[2] = "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithFidAcc"+Extension;
		FSIModel[3] = "SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc"+Extension;
		FSIModel[4] = "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc"+Extension;

	}

	// --------------------------------------------------------------------------------------	

	int NFSIModels = FSIModel.size();

	for (int WhichFSIModel = 0; WhichFSIModel < NFSIModels; WhichFSIModel ++) {

		// --------------------------------------------------------------------------------------

		TString PathToFiles = GlobalPathToFiles + E + "/" + FSIModel[WhichFSIModel] + "/" + xBCut + "/";
		TString FileName = PathToFiles + nucleus +"_" + E + "_" + FSIModel[WhichFSIModel] + "_Plots_FSI_em.root";
		TFile* FileSample = TFile::Open(FileName);

		Plots.push_back( (TH1D*)( FileSample->Get(name) ) );

		UniversalE4vFunction(Plots[WhichFSIModel],FSIModelsToLabels[FSIModel[WhichFSIModel]],nucleus,E,name);

		// --------------------------------------------------------------------------------------

	}

	// --------------------------------------------------------------------------------------	

	// Clone initial histo before you start setting the bin contents

	TH1D* OverallClone = (TH1D*)h->Clone();	

	// --------------------------------------------------------------------------------------	

	// Acceptance Correction

	// Two simulation models: get their average & use it as a correction factor
	// 1: Reco 1p0pi SuSav2 NoRad plot for average
	// 2: Reco 1p0pi G2018 Offset NoRad plot for average
	// 3: True 1p0pi SuSav2 NoRad plot for average
	// 4: True 1p0pi G2018 Offset NoRad plot for average

	TH1D* CorrectionSuSav2 = (TH1D*)Plots[3]->Clone();	
	CorrectionSuSav2->Divide(Plots[1]);
	TH1D* CorrectionG2018 = (TH1D*)Plots[4]->Clone();	
	CorrectionG2018->Divide(Plots[2]);

	TH1D* Average = (TH1D*)(CorrectionSuSav2->Clone());
	Average->Add(CorrectionG2018);
	Average->Scale(0.5);

	// Radiation Correction	// Use SuSav2

	TH1D* RadCorrection = (TH1D*)Plots[1]->Clone();
	RadCorrection->Divide(Plots[0]);

	// --------------------------------------------------------------------------------------	

	int NBins = OverallClone->GetXaxis()->GetNbins();

	double AccCorrTolerance = 30;

	for (int WhichBin = 0; WhichBin < NBins; WhichBin++) {

		double AccCorr = 0.;
		double RadCorr = 0.;

		double NewBinContent = 0.;
		double NewBinError = 0.;		

		//if (Plots[0]->GetBinContent(WhichBin + 1) > 0) { 

			AccCorr = Average->GetBinContent(WhichBin + 1);

			// Sanity checks for acceptance corrections 
			if (AccCorr < 0 || AccCorr > AccCorrTolerance) { 

				double CorrectionSuSav2Bin = CorrectionSuSav2->GetBinContent(WhichBin + 1); 
				double CorrectionG2018Bin = CorrectionG2018->GetBinContent(WhichBin + 1); 
				
				if (CorrectionSuSav2Bin > 0 && CorrectionSuSav2Bin < AccCorrTolerance) { AccCorr = CorrectionSuSav2Bin; } 
				else if (CorrectionG2018Bin > 0 && CorrectionG2018Bin < AccCorrTolerance) { AccCorr = CorrectionG2018Bin; }
				else { AccCorr = 0.; } 

			}

			RadCorr = RadCorrection->GetBinContent(WhichBin + 1);

			NewBinContent = h->GetBinContent(WhichBin + 1) * AccCorr * RadCorr;
			NewBinError = h->GetBinError(WhichBin + 1) * AccCorr * RadCorr;

//cout << "AccCorr = " << AccCorr << endl;
//cout << "RadCorr = " << RadCorr << endl;
//cout << "h->GetBinContent(WhichBin + 1) = " << h->GetBinContent(WhichBin + 1) << "   NewBinContent = " << NewBinContent << endl;
//cout << "h->GetBinError(WhichBin + 1) = " << h->GetBinError(WhichBin + 1) << "   NewBinError = " << NewBinError << endl << endl;

		//}

		OverallClone->SetBinContent(WhichBin + 1, NewBinContent);
		OverallClone->SetBinError(WhichBin + 1, NewBinError);

	}

	// --------------------------------------------------------------------------------------	

//	// To quickly plot the acceptance correction

//	TString CanvasName = "AccCorrCanvas";
//	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
//	Average->GetXaxis()->SetRangeUser(0.1,0.7);
//	Average->GetYaxis()->SetTitle("Detector Acceptance Correction");
//	Average->Draw("e hist");

//	// To quickly plot the radiation correction

//	TString RadCanvasName = "RadCorrCanvas";
//	TCanvas* RadPlotCanvas = new TCanvas(RadCanvasName,RadCanvasName,205,34,1024,768);
//	RadCorrection->GetXaxis()->SetRangeUser(0.1,0.7);
//	RadCorrection->GetYaxis()->SetTitle("Radiative Correction");
//	RadCorrection->Draw("e hist");

	// --------------------------------------------------------------------------------------

	// Obtain acceptance correction uncertainy using non radiative samples

	// First, grab the files with smearing at truth level (_Truth_WithoutFidAcc_Offset), otherwise infinities at the edges 
	// Second, grab the G2018 true 1p0pi files with offset (_Truth_WithFidAcc_Offset), otherwise infinities at the double ratio

	std::vector<TH1D*> PlotsOffset; PlotsOffset.clear();
	std::vector<TString> FSIModelOffset; FSIModelOffset.clear();

	FSIModelOffset.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc"+Extension); // main reco plots for unfolding uncertainty with smearing
	FSIModelOffset.push_back("SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing"+Extension); // main plots for unfolding uncertainty with smearing
	FSIModelOffset.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc_Offset"+Extension); // alternative model plots for acceptance correction uncertainty with smearing & offset 
	FSIModelOffset.push_back("hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing_Offset"+Extension); // alternative model plots for acceptance correction uncertainty with smearing & offset

	if (
		string(name).find("T2KEQEReso") != std::string::npos ||
		name == "h_Erec_subtruct_piplpimi_noprot_3pi" || 
		name == "h_Erec_subtruct_piplpimi_noprot_frac_feed" ||
		name == "h_Erec_subtruct_piplpimi_noprot_frac_feed3pi"

	) {

		FSIModelOffset[0] = "SuSav2_NoRadCorr_LFGM_Truth0pi_WithFidAcc";
		FSIModelOffset[1] = "SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc_Smearing";
		FSIModelOffset[2] = "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithFidAcc";
		FSIModelOffset[3] = "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc_Smearing";

	}

	// --------------------------------------------------------------------------------------	

	int NFSIModelsOffset = FSIModelOffset.size();

	for (int WhichFSIModel = 0; WhichFSIModel < NFSIModelsOffset; WhichFSIModel ++) {

		// --------------------------------------------------------------------------------------

		TString PathToFilesOffset = GlobalPathToFiles + E + "/" + FSIModelOffset[WhichFSIModel] + "/" + xBCut + "/";
		TString FileNameOffset = PathToFilesOffset + nucleus +"_" + E + "_" + FSIModelOffset[WhichFSIModel] + "_Plots_FSI_em.root";
		TFile* FileSampleOffset = TFile::Open(FileNameOffset);

		PlotsOffset.push_back( (TH1D*)( FileSampleOffset->Get(name) ) );

		UniversalE4vFunction(PlotsOffset[WhichFSIModel],FSIModelsToLabels[FSIModelOffset[WhichFSIModel]],nucleus,E,name);

		// --------------------------------------------------------------------------------------

	}

	// --------------------------------------------------------------------------------------
		

	// Test mode for acceptance correction uncertainty

	TH1D* NominalModelRatio = (TH1D*)PlotsOffset[1]->Clone();
	NominalModelRatio->Divide(PlotsOffset[0]); // true / reco unfolding for nominal model

	TH1D* AlternativeModelRatio = (TH1D*)PlotsOffset[3]->Clone();
	AlternativeModelRatio->Divide(PlotsOffset[2]); // true / reco unfolding for alternative model

	TH1D* average = (TH1D*)(NominalModelRatio->Clone());
	average->Add(AlternativeModelRatio);
	average->Scale(0.5);

	TH1D* Spread = (TH1D*)(NominalModelRatio->Clone());
	Spread->Add(AlternativeModelRatio,-1);
	Spread->Divide(average);
	Spread->Scale(1./TMath::Sqrt(12.));

	int NBinsSpread = Spread->GetXaxis()->GetNbins();

	// make sure that the bin entries are positive

	for (int WhichBin = 0; WhichBin < NBinsSpread; WhichBin++) {

		double BinContent = Spread->GetBinContent(WhichBin+1);
		if (BinContent < 0) { Spread->SetBinContent(WhichBin+1,-BinContent); }

	}

	// --------------------------------------------------------------------------------------	

	// Special case for Ecal
	// To avoid the infinities around the peak
	// We merge the (4) bins around the peak

	double DoubleE = -99., reso = 0.;
	if (E == "1_161") { DoubleE = 1.161; reso = 0.07; }
	if (E == "2_261") { DoubleE = 2.261; reso = 0.08; }
	if (E == "4_461") { DoubleE = 4.461; reso = 0.06; }

	double sum = 0; int nbins = 0;

	if (string(name).find("epRecoEnergy_slice") != std::string::npos || name == "h1_Ecal_Reso" || name == "h_Etot_subtruct_piplpimi_2p1pi_1p0pi_fracfeed" ) {

		// Loop over the bins and take the average of the bins around the peak

		for (int WhichBin = 1; WhichBin <= NBinsSpread; WhichBin++) {

			double BinCenter = Spread->GetBinCenter(WhichBin);
			double BinContent = TMath::Abs(Spread->GetBinContent(WhichBin));

			if (BinCenter > (1-reso) * DoubleE && BinCenter < (1+reso) * DoubleE ) {

				sum += BinContent; nbins++;

			}

		}

		sum = sum / double(nbins);

		// -----------------------------------------------------------------------------

		for (int WhichBin = 1; WhichBin <= NBinsSpread; WhichBin++) {

			double BinCenter = Spread->GetBinCenter(WhichBin);
			double BinContent = Spread->GetBinContent(WhichBin);

			if (BinCenter > (1-reso) * DoubleE && BinCenter < (1+reso) * DoubleE ) {

				Spread->SetBinContent(WhichBin,sum);

			}

		}

	}

	// now back to the main plot
	// add the acceptance correction errors in quadrature 
	// along with the other uncertainties (they have already been added)
	// Keep in mind the special case for Ecal, where we merge things around the peak to avoid infinities

	TH1D* AccCorrUncFracPlot = (TH1D*)(OverallClone->Clone("AccCorrUncFracPlot"));

	for (int WhichBin = 0; WhichBin < NBinsSpread; WhichBin++) {

		double SpreadBinContent = Spread->GetBinContent(WhichBin+1); // 0.XY so literally percentage error
		double XSecBinError = OverallClone->GetBinError(WhichBin+1); // pre existing bin error
		double XSecBinEntry = OverallClone->GetBinContent(WhichBin+1); // bin CV
		double AccCorrError = SpreadBinContent * XSecBinEntry;

		AccCorrUncFracPlot->SetBinContent(WhichBin+1,SpreadBinContent);
		AccCorrUncFracPlot->SetBinError(WhichBin+1,0.0000001);

		double NewXSecBinError = TMath::Sqrt( TMath::Power(XSecBinError,2.) + TMath::Power(AccCorrError,2.) );  
		OverallClone->SetBinError(WhichBin+1,NewXSecBinError); // final bin error


	double BinCenter = OverallClone->GetBinCenter(WhichBin+1);
	//cout << "Bin center = " << BinCenter << " XSecBinError = " << XSecBinError << "  SpreadBinContent = " << SpreadBinContent << "  AccCorrError = " << AccCorrError << endl;
	//cout << "Bin center = " << BinCenter << " XSecBinError = " << XSecBinError << "  NewXSecBinError = " << NewXSecBinError << endl;

	}

	// --------------------------------------------------------------------------------------	

	// To quickly plot the acceptance correction unc factor

//	TString AccCorrUncCanvasName = "AccCorrUncCanvas";
//	TCanvas* AccCorrUncPlotCanvas = new TCanvas(AccCorrUncCanvasName,AccCorrUncCanvasName,205,34,1024,768);
//	AccCorrUncFracPlot->GetXaxis()->SetRangeUser(0.1,0.7);
//	AccCorrUncFracPlot->GetYaxis()->SetTitle("Detector Acceptance Correction Fractional Uncertainty");
//	AccCorrUncFracPlot->Draw("e hist");

	// --------------------------------------------------------------------------------------	

	return OverallClone;

}

// -------------------------------------------------------------------------------------------------------------------------------------

void UniversalE4v2DFunction(TH2D* h, TString DataSetLabel, TString nucleus, TString E, TString name) {

	// Scale to obtain absolute double differential cross sections 
	// Use charge, density and length for data samples
	// Use total number of events in genie sample and relevant genie cross sections for simulation

	AbsoluteXSec2DScaling(h,DataSetLabel,nucleus,E);

	// Division by the bin width

//	ReweightPlots(h);

//	// Rebin is necessary

//	ApplyRebinning(h,E,name);

//	// Use relevant ranges
//			
//	ApplyRange(h,E,name);

//	// if data sample: 
//	//                 apply systematics due to rotations et al

//	if (string(DataSetLabel).find("Data") != std::string::npos) { ApplySystUnc(h, E); }

//	//                 apply acceptance systematics using sector-by -sector uncertainties

//	if (string(DataSetLabel).find("Data") != std::string::npos) { ApplySectorSystUnc(h, E); }
}


// -------------------------------------------------------------------------------------------------------------------------------------

double computeMean(std::vector<double> numbers,std::vector<double> weights = {}) {

	if(numbers.empty()) return 0;

	int NBins = (int)(numbers.size());

	if (numbers.size() != weights.size()) { 

		weights.resize( NBins );

		for (int bin = 0; bin < NBins; bin ++) {

			weights[bin] = 1.;

		}

	}

	double total = 0;
	double sumweights = 0;

	for (int number = 0; number < NBins; number ++) {

		total += numbers[number] * weights[number];
		sumweights += weights[number];
	}

	if (sumweights == 0) { return 0.; }

	double average = total / sumweights;

	return average;
}

// -------------------------------------------------------------------------------------------------------------------------------------

double computeStd(double mean, std::vector<double> numbers,std::vector<double> weights = {}) {

	if(numbers.empty()) return 0;

	int NBins = (int)(numbers.size());

	if (numbers.size() != weights.size()) { 

		weights.resize( NBins );

		for (int bin = 0; bin < NBins; bin ++) {

			weights[bin] = 1.;

		}

	}

	double DiffToMean = 0;
	double sumweights = 0;
	int M = 0; // Non zero weights

	for (int number = 0; number < NBins; number ++) {

		DiffToMean += weights[number]*(numbers[number] - mean)*(numbers[number] - mean);
		sumweights += weights[number];
		if (weights[number] != 0) { M++; }

	}

	if (sumweights == 0) { return 0.; }

	return sqrt(DiffToMean / ( (double)(M-1)/(double)(M) * sumweights) );
}

// -------------------------------------------------------------------------------------------------------------------------------------

double SumSqDiffInBin(std::vector<double> numbers, double mean = 0) {

	if(numbers.empty()) return 0;

	int NBins = (int)(numbers.size());

	double sum = 0;

	for (int number = 0; number < NBins; number ++) {

		sum += TMath::Power(numbers[number] - mean,2.);

	}

	return sum;

}

// -------------------------------------------------------------------------------------------------------------------------------------

double Chi2(TH1D* h1,TH1D* h2, int LowBin = -1, int HighBin = -1) {

	int NBinsX = h1->GetXaxis()->GetNbins();

	double chi2 = 0;
	
	if (LowBin == -1) { LowBin = 0; }
	if (HighBin == -1) { HighBin = NBinsX; }	

	for (int WhichXBin = LowBin; WhichXBin < HighBin; WhichXBin++) {

		double h1Entry = h1->GetBinContent(WhichXBin+1);
		double h1Error = h1->GetBinError(WhichXBin+1);
		double h2Entry = h2->GetBinContent(WhichXBin+1);
		double h2Error = h2->GetBinError(WhichXBin+1);

		double num = TMath::Power(h1Entry - h2Entry,2.);
		double den = TMath::Power(h1Error,2.) + TMath::Power(h2Error,2.);
		if (den != 0) { chi2 += (num / den); }

	}

	return chi2;

}

// -------------------------------------------------------------------------------------------------------------------------------------

TPad* CreateTPad(int WhichEnergy, int WhichNucleus, double Energy, TString nucleus, TString name, TCanvas* PlotCanvas) {

	double XMinPad = Xmin + WhichEnergy * Xstep, XMaxPad = Xmin + ( WhichEnergy + 1) * Xstep;
	if (Energy == 1.161 ) { XMinPad = XMinPad - 0.05; }
	double YMinPad = Ymax - ( WhichNucleus + 1) * Ystep, YMaxPad = Ymax - WhichNucleus * Ystep;

	TPad* pad = new TPad(); 

	if (nucleus == "12C") 
	{ pad = new TPad(name,name,XMinPad,YMinPad,XMaxPad,YMaxPad, 21); }
	else { 
		if (Energy == 2.261) 
			{ pad = new TPad(name,name,XMinPad-0.031,YMinPad+space,XMaxPad,YMaxPad+space, 21); }
		else { pad = new TPad(name,name,XMinPad,YMinPad+space,XMaxPad,YMaxPad+space, 21); }
	}

	pad->SetFillColor(kWhite); 
	PlotCanvas->cd();
	pad->Draw();
	pad->cd();

	pad->SetBottomMargin(0.15);

	pad->SetTopMargin(0.0);
	if (nucleus == "12C") { pad->SetTopMargin(0.01); }

	pad->SetRightMargin(0.0);
	if (Energy == 4.461 ) { pad->SetRightMargin(0.01); }

	pad->SetLeftMargin(0.);
	if (Energy == 1.161 ) { pad->SetLeftMargin(0.135); }
	if (Energy == 2.261 && nucleus == "56Fe") { pad->SetLeftMargin(0.1); }	

	return pad;

}

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetFrameWidthNDC()
  {
    return 1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetUxmax()
  {
    gPad->Update();
    return gPad->GetUxmax();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetUxmin()
  {
    gPad->Update();
    return gPad->GetUxmin();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetUymin()
  {
    gPad->Update();
    return gPad->GetUymin();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetUymax()
  {
    gPad->Update();
    return gPad->GetUymax();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetFrameWidthAxis()
  {
    return GetUxmax() - GetUxmin();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  AxisToUser(const double x, const bool logx)
  {
    return logx ? std::pow(10.0, x) : x;
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  AxisToX(const double ax)
  {
    return AxisToUser(ax, gPad->GetLogx());
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetRatioWidthNDCAxis()
  {
    return GetFrameWidthNDC() / GetFrameWidthAxis();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  NDCtoX(const double ndcx)
  {
    const double ax =
        (ndcx - gPad->GetLeftMargin()) / GetRatioWidthNDCAxis() + GetUxmin();
    return AxisToX(ax);
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetFrameHeightNDC()
  {
    return 1.0 - gPad->GetTopMargin() - gPad->GetBottomMargin();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetFrameHeightAxis()
  {
    return GetUymax() - GetUymin();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  GetRatioHeightNDCAxis()
  {
    return GetFrameHeightNDC() / GetFrameHeightAxis();
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  AxisToY(const double ay)
  {
    return AxisToUser(ay, gPad->GetLogy());
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  NDCtoY(const double ndcy)
  {
    const double ay =
        (ndcy - gPad->GetBottomMargin()) / GetRatioHeightNDCAxis() + GetUymin();
    return AxisToY(ay);
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  UserToAxis(const double x, const bool logx)
  {
    return logx ? std::log10(x) : x;
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  XtoAxis(const double x)
  {
    return UserToAxis(x, gPad->GetLogx());
  }

// -------------------------------------------------------------------------------------------------------------------------------------

  double
  XtoNDC(const double x)
  {
    const double ax = XtoAxis(x);
    return GetRatioWidthNDCAxis() * (ax - GetUxmin()) + gPad->GetLeftMargin();
  }

// -------------------------------------------------------------------------------------------------------------------------------------



