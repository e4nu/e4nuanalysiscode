#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <TF1.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

	// ------------------------------------------------------------------------------------

	// Constants

	// ------------------------------------------------------------------------------------

	// Dimensions of TPads

	double Xmin = 0.14, Xmax = 1.;
	double Ymax = 1., Ymin = 0.05;
	double Xstep = (Xmax - Xmin) / 3.;
	double Ystep = ( Ymax - Ymin  ) / 2.;
	double space = 0.07;

	// ------------------------------------------------------------------------------------

	const int Ndivisions = 6;
	const int LineWidth = 3;
	const int FontStyle = 132;
	const double TextSize = 0.07;
	const int MarkerStyle = 20;
	const int MarkerSize = 2.;
	const double acc = 3.;

	// ------------------------------------------------------------------------------------

	const TString GlobalPathToFiles = "/home/afroditi/Dropbox/PhD/myCode/30th_Refactorization/myFiles/";		
	const TString version = "v3_0_6/";	
/*	const TString DoubleXSecTitle = "#frac{d^{2}#sigma}{d#Omega dE} [#frac{#mub}{sr GeV nucleus}]";*/
	const TString DoubleXSecTitle = "Normalized Yield";
	const TString DoubleAccCorrXSecTitle = "Cross Section";

	// Conversion factors for extraction of absolute cross sections

	const double ConversionFactorCm2ToMicroBarn = TMath::Power(10.,30.);

	// Uncertainties from Mariana's analysis note

	const double SystUnc1GeV = 0.02; // 2% syst uncertainty at 1.161 GeV
	const double SystUnc2GeV = 0.021; // 2.1% syst uncertainty at 2.261 GeV
	const double SystUnc4GeV = 0.047; // 4.7% syst uncertainty at 4.461 GeV

	// Sector Uncertainties from Afro's study

	const double SectorSystUnc1GeV = 0.06; // 6% syst uncertainty at 1.161 GeV
	const double SectorSystUnc2GeV = 0.06; // 6% syst uncertainty at 2.261 GeV
	const double SectorSystUnc4GeV = 0.06; // 6% syst uncertainty at 4.461 GeV

	// 20 % overall normalization uncertainty 

//	const double OverallNormUnc = 0.2;
	const double OverallNormUnc = 0.;

	// 20 % acceptance correction uncertainty 

//	const double AcceptanceCorrUnc = 0.2;
	const double AcceptanceCorrUnc = 0.;

	// Larry/Axel's suggestion for scaling down the last 2 bins by EnhaceTail

	const double EnhaceTail = 1./1.;

	// Clas dOmega 

//	const double dOmega = 0.02; // sr
	const double dOmega = 0.015; // sr

	// 1e -> 1.6x10^-19 C
	// 1C -> 6.25x10^18 e
	// 1mC -> 6.25x10^15 e
	// Thus the numbers above have to be multiplied by this number to make sure that we refer to electrons and not charge

	double ConversionFactorChargeToElectrons = 6.25*TMath::Power(10.,15.);

	// Avogadro constant: 6x10^23
	// number of atoms in 12 grams of the isotope 12C
	// 1 gr -> 6x10^23 / 12 = 5x10^22 atoms
//	double ConversionFactorGramToAtoms = 5*TMath::Power(10.,22);
	double AvogadroNumber = 6.02*TMath::Power(10.,23);
	double OverallUnitConversionFactor = ConversionFactorChargeToElectrons * AvogadroNumber;

	// ----------------------------------------------------------------------------------------------------------------------------------------------

	// Maps 

	// Plot names to plot labels

	static std::map<TString,TString> PlotNamesToLabels =
	{
		{ "DeltaAlphaT", "(e,e'p)_{1p0#pi} #delta#alpha_{T} [deg]" },
		{ "DeltaPhiT", "(e,e'p)_{1p0#pi} #delta#phi_{T} [deg]" },
		{ "ECalReso", "E^{cal} Feeddown" },
		{ "EQEReso", "E^{QE} Feeddown" },
		{ "ECal", "(e,e'p)_{1p0#pi} E^{Cal} [GeV]" },
		{ "DeltaPT", "(e,e'p)_{1p0#pi} P_{T} [GeV/c]" }
	};

	// FSI Models to Labels

	static std::map<TString,TString> FSIModelsToLabels =
	{
		{ "MikhailCook_Data", "Mikhail Data" },
		{ "Pinned_Data_Final", "Pinned Data" },
		{ "Pinned_Data_Final_XSec", "Pinned Data" },
		{ "Pinned_Data_NewFiducials_SixSectors", "Pinned Data" },
		{ "Pinned_Data_Final_SixSectors", "Pinned Data" },

		// ------------------------------------------------------------------------------------------------------

		{ "SuSav2_RadCorr_LFGM_Schwinger", "SuSav2 Rad Schwinger" },

		{ "SuSav2_RadCorr_LFGM_UpdatedSchwinger", "SuSav2 Rad Updated Schwinger" },
		{ "SuSav2_RadCorr_LFGM_UpdatedSchwinger_Offset", "SuSav2 Rad Updated Schwinger" },
		{ "SuSav2_RadCorr_LFGM_UpdatedSchwinger_Test", "SuSav2 Rad Updated Schwinger" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger", "SuSav2 Rad Updated Schwinger" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger_XSec", "SuSav2 Rad Updated Schwinger" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger_Offset", "SuSav2 Rad Updated Schwinger" },
		{ "SuSav2_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger", "SuSav2 Rad Updated Schwinger" },
		{ "SuSav2_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger_Offset", "SuSav2 Rad Updated Schwinger" },

		{ "SuSav2_RadCorr_LFGM", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_XSec", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithFidAcc", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithFidAcc_XSec", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc_Offset", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc_XSec", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc_NoThetaCut", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth0pi_WithFidAcc", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth0pi_WithoutFidAcc", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_SixSectors", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Filtered", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithFidAcc_Filtered", "SuSav2" },
		{ "SuSav2_RadCorr_LFGM_Truth_WithoutFidAcc_Filtered", "SuSav2" },

		{ "SuSav2_NoRadCorr_LFGM", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_Test", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_Offset", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing_Offset", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing_XSec", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc_Smearing", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc_Offset", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc_Test", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth0pi_WithFidAcc_Test", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_XSec", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithoutFidAcc_XSec", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc_XSec", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth_WithFidAcc_XSec", "SuSav2 NoRad" },
		{ "SuSav2_NoRadCorr_LFGM_Truth0pi_WithFidAcc", "SuSav2 NoRad" },

		// ------------------------------------------------------------------------------------------------------

		{ "hA2018_Final_RadCorr_LFGM", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Offset", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Truth_WithFidAcc", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Truth_WithFidAcc_Offset", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Truth_WithoutFidAcc", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Truth_WithoutFidAcc_NoThetaCut", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Truth_WithoutFidAcc_Offset", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Truth0pi_WithFidAcc", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_Truth0pi_WithoutFidAcc", "G2018" },
		{ "hA2018_Final_RadCorr_LFGM_SixSectors", "G2018" },

		{ "hA2018_Final_RadCorr_LFGM_UpdatedSchwinger", "G2018 Rad Updated Schwinger" },
		{ "hA2018_Final_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger", "G2018 Rad Updated Schwinger" },
		{ "hA2018_Final_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger", "G2018 Rad Updated Schwinger" },
		{ "hA2018_Final_RadCorr_LFGM_UpdatedSchwinger_Offset", "G2018 Rad Updated Schwinger" },
		{ "hA2018_Final_RadCorr_LFGM_Truth_WithFidAcc_UpdatedSchwinger_Offset", "G2018 Rad Updated Schwinger" },
		{ "hA2018_Final_RadCorr_LFGM_Truth0pi_WithFidAcc_UpdatedSchwinger_Offset", "G2018 Rad Updated Schwinger" },

		{ "hA2018_Final_NoRadCorr_LFGM", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Offset", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Playground", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Offset", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Offset_XSec", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing_Offset_XSec", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_Smearing_Offset", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc_Smearing", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithoutFidAcc_Offset", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_XSec", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc_XSec", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc_Offset_XSec", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithoutFidAcc_XSec", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth_WithFidAcc_Offset", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithFidAcc", "G2018 NoRad" },
		{ "hA2018_Final_NoRadCorr_LFGM_Truth0pi_WithFidAcc_Offset", "G2018 NoRad" },

		// ------------------------------------------------------------------------------------------------------

		{ "SuSav2_MasterNoRad", "SuSav2 Master NoRad" },
		{ "SuSav2_MasterNoRad_Truth_WithFidAcc", "SuSav2 Master NoRad" },
		{ "SuSav2_MasterNoRad_Truth0pi_WithFidAcc", "SuSav2 Master NoRad" },
		{ "SuSav2_MasterNoRad_Truth_WithoutFidAcc", "SuSav2 Master NoRad" },
		{ "SuSav2_MasterNoRad_Truth0pi_WithoutFidAcc", "SuSav2 Master NoRad" },

		{ "SuSav2_MasterRad", "SuSav2 Master Rad" },
		{ "SuSav2_MasterRad_Truth_WithFidAcc", "SuSav2 Master Rad" },
		{ "SuSav2_MasterRad_Truth0pi_WithFidAcc", "SuSav2 Master Rad" },
		{ "SuSav2_MasterRad_Truth_WithoutFidAcc", "SuSav2 Master Rad" },
		{ "SuSav2_MasterRad_Truth0pi_WithoutFidAcc", "SuSav2 Master Rad" },

		// ------------------------------------------------------------------------------------------------------

		{ "G2018_MasterNoRad", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth_WithFidAcc", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth0pi_WithFidAcc", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth_WithoutFidAcc", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth0pi_WithoutFidAcc", "G2018 Master NoRad" },

		{ "G2018_MasterRad", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth_WithFidAcc", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth0pi_WithFidAcc", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth_WithoutFidAcc", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth0pi_WithoutFidAcc", "G2018 Master Rad" },

		{ "G2018_MasterRad_QEOnly", "G2018 QE Only" },

		// ------------------------------------------------------------------------------------------------------

		{ "G2018_MasterNoRad_Offset", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth_WithFidAcc_Offset", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth0pi_WithFidAcc_Offset", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth_WithoutFidAcc_Offset", "G2018 Master NoRad" },
		{ "G2018_MasterNoRad_Truth0pi_WithoutFidAcc_Offset", "G2018 Master NoRad" },

		{ "G2018_MasterRad_Offset", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth_WithFidAcc_Offset", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth0pi_WithFidAcc_Offset", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth_WithoutFidAcc_Offset", "G2018 Master Rad" },
		{ "G2018_MasterRad_Truth0pi_WithoutFidAcc_Offset", "G2018 Master Rad" },

		// ------------------------------------------------------------------------------------------------------

		{ "GTEST18_02c_NoRadCorr_LFGM_Truth_WithFidAcc", "GTEST18_02c NoRad" },

		// ------------------------------------------------------------------------------------------------------

		{ "GTEST18_02d_NoRadCorr_LFGM_Truth_WithFidAcc", "GTEST18_02d NoRad" },

		// ------------------------------------------------------------------------------------------------------



	};

	// Mass Numbers

	static std::map<TString,double> MassNumber =
	{
		{ "1H", 1 },
		{ "4He", 4 },
		{ "12C", 12 },
		{ "CH2", 14 },
		{ "56Fe", 56 }
	};

	// mC // Regular files // Data_Final

	static std::map<std::pair<TString,TString>,double> IntegratedCharge =
	{
		{ std::make_pair("4He", "2_261"), 1.08 },
		{ std::make_pair("4He", "4_461"), 0.87 },
//		{ std::make_pair("12C", "1_161"), 0.19 },
		{ std::make_pair("12C", "1_161"), 0.079 },
		{ std::make_pair("12C", "2_261"), 1.79 },
		{ std::make_pair("12C", "4_461"), 2.14 },
		{ std::make_pair("56Fe", "2_261"), 0.22 },
		{ std::make_pair("56Fe", "4_461"), 0.29 }
	};

	// mC // Good run list all runs // Mikhail cook pass3

	static std::map<std::pair<TString,TString>,double> IntegratedCharge_MikhailFiles =
	{
		{ std::make_pair("4He", "2_261"), 0. },
		{ std::make_pair("4He", "4_461"), 0. },
		{ std::make_pair("12C", "1_161"), 0. },
		{ std::make_pair("12C", "2_261"), 2.47238 },
		{ std::make_pair("12C", "4_461"), 2.06258 },
		{ std::make_pair("CH2", "1_161"), 0. },
		{ std::make_pair("CH2", "2_261"), 0.31532 },
		{ std::make_pair("CH2", "4_461"), 0.17789 },
		{ std::make_pair("56Fe", "2_261"), 0.20595 },
		{ std::make_pair("56Fe", "4_461"), 0.24576 }
	};

	// mC // Good run list all runs // pinned files by Stuart F

	static std::map<std::pair<TString,TString>,double> IntegratedCharge_PinnedFiles =
	{
		{ std::make_pair("4He", "2_261"), 1.16584 },
		{ std::make_pair("4He", "4_461"), 0.97884 },
		{ std::make_pair("12C", "1_161"), 0.079 },
//		{ std::make_pair("CH2", "1_161"), 0.070707652 },

//		{ std::make_pair("CH2", "1_161"), 0.0794 },
		{ std::make_pair("CH2", "1_161"), 0.068 }, // L.W. Dec 17 2020

		{ std::make_pair("12C", "2_261"), 2.83649 },
//		{ std::make_pair("12C", "2_261"), 0.007609864 },
		{ std::make_pair("12C", "4_461"), 2.31146 },
		{ std::make_pair("56Fe", "2_261"), 0.217238 },
		{ std::make_pair("56Fe", "4_461"), 0.308581 }
	};

	// mC // Filtered runs

	static std::map<std::pair<TString,TString>,double> IntegratedCharge_FilterRuns =
	{
		{ std::make_pair("4He", "2_261"), 0.024727582 },
		{ std::make_pair("4He", "4_461"), 0. },
		{ std::make_pair("12C", "1_161"), 0.05392387 },
		{ std::make_pair("12C", "2_261"), 0.060684561 },
		{ std::make_pair("12C", "4_461"), 0.099557913 },
		{ std::make_pair("56Fe", "2_261"), 0. },
		{ std::make_pair("56Fe", "4_461"), 0.089507691 }
	};

	// mC // New Filtered runs

	static std::map<std::pair<TString,TString>,double> IntegratedCharge_NewFilterRuns =
	{
		{ std::make_pair("4He", "2_261"), 0. },
		{ std::make_pair("4He", "4_461"), 0. },
		{ std::make_pair("12C", "1_161"), 0.070707652 },
		{ std::make_pair("12C", "2_261"), 0. },
		{ std::make_pair("12C", "4_461"), 0. },
		{ std::make_pair("56Fe", "2_261"), 0. },
		{ std::make_pair("56Fe", "4_461"), 0. }
	};

	// mC // Good run list all runs
	static std::map<std::pair<TString,TString>,double> IntegratedCharge_GoodRunList_AllRuns =
	{
		{ std::make_pair("4He", "2_261"), 1.16584 },
		{ std::make_pair("4He", "4_461"), 0.97884 },
		{ std::make_pair("12C", "1_161"), 0.18432 },
//		{ std::make_pair("12C", "2_261"), 2.8682265 },
		{ std::make_pair("12C", "2_261"), 0.007609864 },
		{ std::make_pair("12C", "4_461"), 2.31146 },
		{ std::make_pair("56Fe", "2_261"), 0.217238 },
		{ std::make_pair("56Fe", "4_461"), 0.308581 }
	};

	// mC // Good run list low current runs

	static std::map<std::pair<TString,TString>,double> IntegratedCharge_GoodRunList_LowCurrentRuns =
	{
		{ std::make_pair("4He", "2_261"), 0. },
		{ std::make_pair("4He", "4_461"), 0. },
		{ std::make_pair("12C", "1_161"), 0.079 },
		{ std::make_pair("12C", "2_261"), 0. },
		{ std::make_pair("12C", "4_461"), 0. },
		{ std::make_pair("56Fe", "2_261"), 0. },
		{ std::make_pair("56Fe", "4_461"), 0. }
	};

	// mC // Good run list high current runs

	static std::map<std::pair<TString,TString>,double> IntegratedCharge_GoodRunList_HighCurrentRuns =
	{
		{ std::make_pair("4He", "2_261"), 0. },
		{ std::make_pair("4He", "4_461"), 0. },
		{ std::make_pair("12C", "1_161"), 0.105 },
		{ std::make_pair("12C", "2_261"), 0. },
		{ std::make_pair("12C", "4_461"), 0. },
		{ std::make_pair("56Fe", "2_261"), 0. },
		{ std::make_pair("56Fe", "4_461"), 0. }
	};

	// cm

	static std::map<std::pair<TString,TString>,double> TargetLength =
	{
		{ std::make_pair("4He", "2_261"), 4.3 }, // 3.72-4.99 // taking the average
		{ std::make_pair("4He", "4_461"), 4.3 }, // 3.72-4.99 // taking the average
		{ std::make_pair("CH2", "1_161"), 0.07 },
		{ std::make_pair("12C", "1_161"), 0.1 },
		{ std::make_pair("12C", "2_261"), 0.1 },
		{ std::make_pair("12C", "4_461"), 0.1 },
		{ std::make_pair("56Fe", "2_261"), 0.015 },
		{ std::make_pair("56Fe", "4_461"), 0.015 }
	};

	// g/cm^3

	static std::map<std::pair<TString,TString>,double> TargetDensity =
	{
		{ std::make_pair("4He", "2_261"), 0.125 },
		{ std::make_pair("4He", "4_461"), 0.125 },
		{ std::make_pair("CH2", "1_161"), 1.392 },
		{ std::make_pair("12C", "1_161"), 1.786 },
		{ std::make_pair("12C", "2_261"), 1.786 },
		{ std::make_pair("12C", "4_461"), 1.786 },
		{ std::make_pair("56Fe", "2_261"), 7.872 },
		{ std::make_pair("56Fe", "4_461"), 7.872 }
	};

	// SuSav2 GENIE spline xsec // 10^{-38} cm^2

	static std::map<std::pair<TString,TString>,double> SuSav2GenieXSec =
	{
		{ std::make_pair("4He", "2_261"), 8.30934e+07 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"), 3.01721e+07 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"), 1.28967e+09 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"), 2.1024e+08 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"), 8.36795e+07 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 9.66272e+08 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 3.84607e+08 } // Q2 > 0.8
	};

	// Rad+Schwinger SuSav2 GENIE number events 

	static std::map<std::pair<TString,TString>,double> RadSchwingerSuSav2NumberEvents =
	{
		{ std::make_pair("4He", "2_261"), 0 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"), 0 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"), 29300000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"), 27700000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"), 26500000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 0 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.8
	};

	// Rad+UpdatedSchwinger SuSav2 GENIE number events 

	static std::map<std::pair<TString,TString>,double> RadUpdatedSchwingerSuSav2NumberEvents =
	{

		{ std::make_pair("4He", "2_261"), 42900000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"), 44700000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"), 33500000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"), 34100000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"), 37700000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 34500000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 49900000 } // Q2 > 0.8

/*		{ std::make_pair("4He", "2_261"), 81100000 }, // Q2 > 0.4*/
/*		{ std::make_pair("4He", "4_461"), 92300000 }, // Q2 > 0.8*/
/*		{ std::make_pair("12C", "1_161"), 99300000 }, // Q2 > 0.1*/
/*		{ std::make_pair("12C", "2_261"), 99800000 }, // Q2 > 0.4*/
/*		{ std::make_pair("12C", "4_461"), 106400000 }, // Q2 > 0.8*/
/*		{ std::make_pair("56Fe", "2_261"), 83700000 }, // Q2 > 0.4*/
/*		{ std::make_pair("56Fe", "4_461"), 99800000 } // Q2 > 0.8*/
	};

	// No Rad SuSav2 GENIE number events 

	static std::map<std::pair<TString,TString>,double> NoRadSuSav2NumberEvents =
	{

		// Produced using adi's rad branch
/*		{ std::make_pair("4He", "2_261"), 0 }, // Q2 > 0.4*/
/*		{ std::make_pair("4He", "4_461"), 0 }, // Q2 > 0.8*/
/*		{ std::make_pair("12C", "1_161"), 49400000 }, // Q2 > 0.1*/
/*		{ std::make_pair("12C", "2_261"), 0 }, // Q2 > 0.4*/
/*		{ std::make_pair("12C", "4_461"), 0 }, // Q2 > 0.8*/
/*		{ std::make_pair("56Fe", "2_261"), 0 }, // Q2 > 0.4*/
/*		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.8*/

		// master branch production
		{ std::make_pair("4He", "2_261"),  20000000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  17700000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  19800000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  174600000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  164300000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 167000000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 190600000 } // Q2 > 0.8
	};

	// Rad SuSav2 GENIE number events 

	static std::map<std::pair<TString,TString>,double> SuSav2NumberEvents =
	{
		{ std::make_pair("4He", "2_261"),  19800000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  20000000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  39700000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  201300000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  179600000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 152800000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 130800000 } // Q2 > 0.8
	};

	// G2018 GENIE spline xsec // 10^{-38} cm^2

	static std::map<std::pair<TString,TString>,double> G2018GenieXSec =
	{
		{ std::make_pair("1H", "1_161"),  1.4515324e+08 }, // Q2 > 0.1
		//(double) 1.1038574e+09 // Q2 > 0.02
		{ std::make_pair("1H", "2_261"),  20943873. }, // Q2 > 0.4
		{ std::make_pair("1H", "4_461"),  8521094.7 }, // Q2 > 0.8
		{ std::make_pair("4He", "2_261"), 6.55943e+07 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"), 2.73355e+07 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  1.10931e+09 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"), 1.96812e+08 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"), 8.20065e+07 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"),9.02272e+08  }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 3.8295048e+08 } // Q2 > 0.8
	};

	// QE Only G2018 GENIE spline xsec // 10^{-38} cm^2

	static std::map<std::pair<TString,TString>,double> QEG2018GenieXSec =
	{
		{ std::make_pair("1H", "2_261"),  0. }, // Q2 > 0.4
		{ std::make_pair("1H", "4_461"),  0. }, // Q2 > 0.8
		{ std::make_pair("4He", "2_261"), 0. }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"), 0. }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  6.0386491e+08 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"), 0. }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"), 0. }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"),0.  }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 0. } // Q2 > 0.8
	};

	// Rad G2018 GENIE number events 

	static std::map<std::pair<TString,TString>,double> G2018NumberEvents =
	{
		{ std::make_pair("1H", "1_161"),   3000000 }, // Q2 > 0.1
		{ std::make_pair("4He", "2_261"),  20000000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  9000000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  50000000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  227100000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  150100000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 50000000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 150100000 } // Q2 > 0.8

	};

	// Rad+UpdatedSchwinger G2018 GENIE number events 

	static std::map<std::pair<TString,TString>,double> RadUpdatedSchwingerG2018NumberEvents =
	{

		{ std::make_pair("4He", "2_261"), 37500000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"), 36700000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"), 44200000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"), 44700000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"), 41400000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 36700000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 40100000 } // Q2 > 0.8

/*		{ std::make_pair("4He", "2_261"), 98700000 }, // Q2 > 0.4*/
/*		{ std::make_pair("4He", "4_461"), 99400000 }, // Q2 > 0.8*/
/*		{ std::make_pair("12C", "1_161"), 57300000 }, // Q2 > 0.1*/
/*		{ std::make_pair("12C", "2_261"), 47400000 }, // Q2 > 0.4*/
/*		{ std::make_pair("12C", "4_461"), 99900000 }, // Q2 > 0.8*/
/*		{ std::make_pair("56Fe", "2_261"), 99600000 }, // Q2 > 0.4*/
/*		{ std::make_pair("56Fe", "4_461"), 99300000 } // Q2 > 0.8*/
	};

	// No Rad G2018 GENIE number events 

	static std::map<std::pair<TString,TString>,double> NoRadG2018NumberEvents =
	{

/*		{ std::make_pair("4He", "2_261"), 0 }, // Q2 > 0.4*/
/*		{ std::make_pair("4He", "4_461"), 0 }, // Q2 > 0.8*/
/*		{ std::make_pair("12C", "1_161"), 0 }, // Q2 > 0.1*/
/*		{ std::make_pair("12C", "2_261"), 0 }, // Q2 > 0.4*/
/*		{ std::make_pair("12C", "4_461"), 0 }, // Q2 > 0.8*/
/*		{ std::make_pair("56Fe", "2_261"), 0 }, // Q2 > 0.4*/
/*		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.8*/

//		{ std::make_pair("1H", "1_161"), 0 }, // Q2 > 0.1
		{ std::make_pair("4He", "2_261"),  69000000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  70400000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  43900000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  53200000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  67900000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 69600000 }, // Q2 > 0.4
//		{ std::make_pair("56Fe", "4_461"), 49749540 } // Q2 > 0.8 // Feb 18
		{ std::make_pair("56Fe", "4_461"), 81100000 } // Q2 > 0.8 // Default
//		{ std::make_pair("56Fe", "4_461"), 46500000 } // Q2 > 0.8

	};

	//  ------------------------------------------------------------------------------
	//  ------------------------------------------------------------------------------

	// Master No Rad G2018 number events 

	static std::map<std::pair<TString,TString>,double> MasterNoRadG2018NumberEvents =
	{

		{ std::make_pair("4He", "2_261"),  48100000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  48100000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  49700000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  49700000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  49500000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 48600000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 48300000 } // Q2 > 0.8

	};

	// Master Rad G2018 number events 

	static std::map<std::pair<TString,TString>,double> MasterRadG2018NumberEvents =
	{

		{ std::make_pair("4He", "2_261"),  48600000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  48100000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  47600000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  48800000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  47500000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 47900000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 47700000 } // Q2 > 0.8

	};

	// QE Only Master Rad G2018 number events 

	static std::map<std::pair<TString,TString>,double> QEMasterRadG2018NumberEvents =
	{

		{ std::make_pair("4He", "2_261"),  0 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  0 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  200000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  0 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  0 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 0 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.8

	};

	// Master No Rad SuSav2 number events 

	static std::map<std::pair<TString,TString>,double> MasterNoRadSuSav2NumberEvents =
	{

		{ std::make_pair("4He", "2_261"),  47200000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  47400000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  48000000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  47800000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  47500000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 0 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 47900000 } // Q2 > 0.8

	};

	// Master Rad SuSav2 number events 

	static std::map<std::pair<TString,TString>,double> MasterRadSuSav2NumberEvents =
	{

		{ std::make_pair("4He", "2_261"),  47100000 }, // Q2 > 0.4
		{ std::make_pair("4He", "4_461"),  47200000 }, // Q2 > 0.8
		{ std::make_pair("12C", "1_161"),  47800000 }, // Q2 > 0.1
		{ std::make_pair("12C", "2_261"),  47800000 }, // Q2 > 0.4
		{ std::make_pair("12C", "4_461"),  47400000 }, // Q2 > 0.8
		{ std::make_pair("56Fe", "2_261"), 45700000 }, // Q2 > 0.4
		{ std::make_pair("56Fe", "4_461"), 46000000 } // Q2 > 0.8

	};

	//  ------------------------------------------------------------------------------

	// No Rad GTEST18_02c GENIE number events 

	static std::map<std::pair<TString,TString>,double> NoRadGTEST18_02cNumberEvents =
	{

		// master branch production
		{ std::make_pair("4He", "2_261"),  0 }, // Q2 > 0.02
		{ std::make_pair("4He", "4_461"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "1_161"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "2_261"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "4_461"),  0 }, // Q2 > 0.02
		{ std::make_pair("56Fe", "2_261"), 0 }, // Q2 > 0.02
		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.02
	};

	// GTEST18_02c GENIE spline xsec // 10^{-38} cm^2

	static std::map<std::pair<TString,TString>,double> GTEST18_02cGenieXSec =
	{
		{ std::make_pair("1H", "1_161"),  0 }, // Q2 > 0.02
		{ std::make_pair("1H", "2_261"),  0 }, // Q2 > 0.02
		{ std::make_pair("1H", "4_461"),  0 }, // Q2 > 0.02
		{ std::make_pair("4He", "2_261"), 0 }, // Q2 > 0.02
		{ std::make_pair("4He", "4_461"), 0 }, // Q2 > 0.02
		{ std::make_pair("12C", "1_161"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "2_261"), 0 }, // Q2 > 0.02
		{ std::make_pair("12C", "4_461"), 0 }, // Q2 > 0.02
		{ std::make_pair("56Fe", "2_261"),0  }, // Q2 > 0.02
		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.02
	};

	// No Rad GTEST18_02d GENIE number events 

	static std::map<std::pair<TString,TString>,double> NoRadGTEST18_02dNumberEvents =
	{

		// master branch production
		{ std::make_pair("4He", "2_261"),  0 }, // Q2 > 0.02
		{ std::make_pair("4He", "4_461"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "1_161"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "2_261"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "4_461"),  0 }, // Q2 > 0.02
		{ std::make_pair("56Fe", "2_261"), 0 }, // Q2 > 0.02
		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.02
	};

	// GTEST18_02d GENIE spline xsec // 10^{-38} cm^2

	static std::map<std::pair<TString,TString>,double> GTEST18_02dGenieXSec =
	{
		{ std::make_pair("1H", "1_161"),  0 }, // Q2 > 0.02
		{ std::make_pair("1H", "2_261"),  0 }, // Q2 > 0.02
		{ std::make_pair("1H", "4_461"),  0 }, // Q2 > 0.02
		{ std::make_pair("4He", "2_261"), 0 }, // Q2 > 0.02
		{ std::make_pair("4He", "4_461"), 0 }, // Q2 > 0.02
		{ std::make_pair("12C", "1_161"),  0 }, // Q2 > 0.02
		{ std::make_pair("12C", "2_261"), 0 }, // Q2 > 0.02
		{ std::make_pair("12C", "4_461"), 0 }, // Q2 > 0.02
		{ std::make_pair("56Fe", "2_261"),0  }, // Q2 > 0.02
		{ std::make_pair("56Fe", "4_461"), 0 } // Q2 > 0.02
	};


	//  ------------------------------------------------------------------------------
	//  ------------------------------------------------------------------------------

//	const std::vector<int> BreakDownColors{kBlue,429,410,610}; // QE, MEC, RES, DIS
	const std::vector<int> BreakDownColors{kBlue+1,kRed-3,kGreen+1,kOrange+1};

	const std::vector<int> SectorColors{kBlack,610,410,kRed+1,kGreen+3,kBlue};
	const std::vector<int> Style{1,1,kDashed,1,1};
	const std::vector<TString> GenieFSILabel{"QE","MEC","RES","DIS"};
	const std::vector<int> DataSetColors{1,1,1,1,1};

	//  ------------------------------------------------------------------------------

	// XSec labels

	TString XSecEcalLabel = "#frac{d#sigma}{dE_{cal}} #left[#frac{#mub}{GeV}#right]";
	TString XSecEQELabel = "#frac{d#sigma}{dE_{QE}} #left[#frac{#mub}{GeV}#right]";

	TString ResoXSecEcalLabel = "#frac{d#sigma}{dE_{cal}^{Feed}} #left[#mub#right]";
	TString ResoXSecEQELabel = "#frac{d#sigma}{dE_{QE}^{Feed}} #left[#mub#right]";

	//  ------------------------------------------------------------------------------


#endif
