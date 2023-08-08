#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TEfficiency.h"
#include <iomanip>
#include <filesystem>

std::string GetAxisLabel( std::string observable, unsigned int id_axis ){
	std::string x_axis, y_axis ;
	if( observable == "ECal") { x_axis = "E_{Cal} [GeV]"; y_axis  = "d#sigma/dE_{Cal} #left[#mub GeV^{-1}#right]"; }
	else if ( observable == "pfl_theta") { x_axis = "#theta_{e'} [deg]"; y_axis  = "d#sigma/d#theta_{e'} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "pfl_phi") { x_axis = "#phi_{e'} [deg]"; y_axis  = "d#sigma/d#phi_{e'} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "pfl") { x_axis = "p_{e'} [GeV/c]"; y_axis  = "d#sigma/dp_{e'} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "proton_mom") { x_axis = "p_{p} [GeV/c]"; y_axis  = "d#sigma/dp_{p} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "proton_theta") { x_axis = "#theta_{p} [deg]"; y_axis  = "d#sigma/d#theta_{p} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "proton_phi") { x_axis = "E_{Cal} [GeV]"; y_axis  = "d#sigma/dE_{Cal} #left[#mub GeV^{-1}#right]"; }
	else if ( observable == "pim_mom") { x_axis = "p_{#pi^{-}} [GeV/c]"; y_axis  = "d#sigma/dp_{#pi^{-}} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "pim_theta") { x_axis = "#theta_{#pi^{-}} [deg]"; y_axis  = "d#sigma/d#theta_{#pi^{-}} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "pip_mom") { x_axis = "p_{#pi^{+}} [GeV/c]"; y_axis  = "d#sigma/dp_{#pi^{+}} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "pip_theta") { x_axis = "#theta_{#pi^{+}} [deg]"; y_axis  = "d#sigma/d#theta_{#pi^{+}} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "RecoW") { x_axis = "W [GeV]"; y_axis  = "d#sigma/dW #left[#mub GeV^{-1}#right#right]"; }
	else if ( observable == "RecoQELEnu") { x_axis = "E^{QE} [GeV]"; y_axis  = "d#sigma/dE^{QE} #left[#mub GeV^{-1}#right#right]"; }
	else if ( observable == "RecoXBJK") { x_axis = "x_{BJK} [GeV]"; y_axis  = "d#sigma/dx_{BJK} #left[#mub GeV^{-1}#right]"; }
	else if ( observable == "RecoQ2") { x_axis = "Q^{2} [GeV]"; y_axis  = "d#sigma/dQ^{2}} #left[#mub GeV^{-1}#right]"; }
	else if ( observable == "Recoq3") { x_axis = "q_{3} [GeV]"; y_axis  = "d#sigma/dq_{3} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "DeltaPT") { x_axis = "#deltap_{T} [GeV]"; y_axis  = "d#sigma/d#deltap_{T} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "HadDeltaPT") { x_axis = "#deltap_{T}^{had} [GeV]"; y_axis  = "d#sigma/d#deltap_{T}^{had} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "DeltaPhiT") { x_axis = "#delta#phi_{T} [deg]"; y_axis  = "d#sigma/d#delta#phi_{T} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "HadDeltaPhiT") { x_axis = "#delta#phi_{T}^{had} [deg]"; y_axis  = "d#sigma/d#delta#phi_{T}^{had} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "AlphaT") { x_axis = "#alpha_{T} [deg]"; y_axis  = "d#sigma/d#alpha_{T} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "HadAlphaT") { x_axis = "#alpha_{T}^{had} [deg]"; y_axis  = "d#sigma/d#alpha_{T}^{had} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "RecoEnergyTransfer") { x_axis = "#omega [GeV]"; y_axis  = "d#sigma/d#omega #left[#mub GeV^{-1}#right]"; }
  else if ( observable == "HadSystemMass") { x_axis = "M_{X}[GeV]"; y_axis = "d#sigma/dM_{X} #left[#mub GeV^{-1}#right]"; }
	if( id_axis ==0 ) return x_axis ;
	return y_axis ;
}


std::vector<double> GetUniformBinning( unsigned int nbins, double min, double max){
  std::vector<double> binning ;
  double step = (max-min)/nbins;
  for( unsigned int i = 0 ; i < nbins + 1 ; ++i ){
    binning.push_back( min + i * step ) ;
  }
  return binning ;
}

std::vector<double> GetBinning( std::string observable, double EBeam ){
	std::vector<double> binning ;
	if( observable == "ECal") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.8, 1.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 1, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 1.5, 5 );
  }	else if ( observable == "pfl_theta") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 20, 45 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 20, 45 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 15, 45 );
  }	else if ( observable == "pfl_phi") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 180 );
  } else if ( observable == "pfl") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.35, 0.9 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0.5, 1.7 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 1, 3.8 );
  } else if ( observable == "proton_mom") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.2, 1.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0.2, 2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 4 );
  }	else if ( observable == "proton_theta") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 110 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 15, 120 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 80 );
  } else if ( observable == "proton_phi") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 180 );
  } else if ( observable == "pim_mom") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.1, 0.7 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 10, 0, 1.5 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 1.6 );
  }	else if ( observable == "pim_theta") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 15, 0, 140 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 30, 140 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 30, 150 );
  }	else if ( observable == "pip_mom") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.1, 0.6);
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 1.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
  }	else if ( observable == "pip_theta") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 15, 0, 150 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 180 );
  }	else if ( observable == "RecoW") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 16, 1, 1.5 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 18, 1, 1.9 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 1.2, 2.5 );
  }	else if ( observable == "RecoXBJK") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 0.9 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 0.9 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 1 );
  } else if ( observable == "RecoQ2") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.15, 0.6 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0.3, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 15, 0.9, EBeam+0.2 );
	} else if ( observable == "RecoQELEnu") {
	    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0., 1.3 );
	    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0.3, EBeam+0.2 );
	    else if( EBeam == 4.461 ) binning = GetUniformBinning( 15, 0.9, EBeam+0.2 );
  }	else if ( observable == "Recoq3") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, EBeam+0.2 );
  }	else if ( observable == "DeltaPT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, EBeam+0.2 );
  }	else if ( observable == "HadDeltaPT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 1);
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 1 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 1 );
  }	else if ( observable == "DeltaPhiT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0,180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 100 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 180 );
  } else if ( observable == "HadDeltaPhiT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 100 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 100 );
  }	else if ( observable == "AlphaT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 180 );
  } else if ( observable == "HadAlphaT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 15, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 15, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 15, 0, 180 );
  }	else if ( observable == "RecoEnergyTransfer") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.3, 0.8 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0.5, 2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 4 );
  } else if ( observable == "HadSystemMass"){
		if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 1, 1.6 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 1, 2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 2.7 );
	}

	return binning ;
}

std::vector<double> GetAdditionalBinning( std::string second_observable, double EBeam ) {
	// In some cases, we might want to add additional plots.
	// In particular, we might want to break the plot into additonal plots as a function of a second observable
	// This function returns the binning for this second observable.
	// The "binning" corresponds to the ranges of interest
	std::vector<double> binning ;
	std::vector<double> original_binning = GetBinning( second_observable, EBeam ) ;
	if( second_observable == "ECal" ) {
		binning.push_back(original_binning[0]);
		binning.push_back(EBeam*(1-0.05));
		binning.push_back(original_binning[original_binning.size()-1]);
	}	else if ( second_observable == "HadDeltaPT" || second_observable == "DeltaPT" ){
		binning.push_back(original_binning[0]);
		binning.push_back(0.2);
		//binning.push_back(0.4);
		binning.push_back(original_binning[original_binning.size()-1]);
	} else if ( second_observable == "HadAlphaT" || second_observable == "AlphaT" ){
		binning.push_back(original_binning[0]);
		binning.push_back(45);
		binning.push_back(original_binning[original_binning.size()-1]);
	}
	return binning;
}

std::string GetAlternativeObs( std::string observable ){
	if( observable == "ECal" ) return "HadDeltaPT" ;
	else if( observable == "HadDeltaPT" || observable == "DeltaPT") return "ECal" ;
	else if( observable == "HadAlphaT"  || observable == "AlphaT" ) return "ECal" ;
	return "";
}

std::string GetObsName( std::string observable ){
	if( observable == "ECal" ) return "E_{Cal}" ;
	else if( observable == "HadDeltaPT" ) return "#deltap^{had}_{T}" ;
	else if( observable == "HadAlphaT"  ) return "#alpha^{had}_{T}" ;
	return "";
}

std::string GetUnit( std::string observable ){
	if( observable == "ECal" ) return "[GeV]" ;
	else if( observable == "HadDeltaPT" ) return "[GeV/c]" ;
	else if( observable == "HadAlphaT"  ) return "[deg]" ;
	return "";
}

double GetMaximum( std::vector<TH1D*> predictions){
	double max = 0;
	for( unsigned int i = 0 ; i < predictions.size();++i){
		if ( max < predictions[i] -> GetMaximum() ) max = predictions[i] -> GetMaximum();
	}
	return max*(1 + 0.12);
}

void StandardFormat( TH1D * prediction, std::string title, int color, int style, std::string observable, double y_max = 0 ) {
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetPaperSize(20,26);
  gStyle->SetTitleFont(132,"pad");
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineStyleString(2,"[12 12]");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  prediction -> SetLineColor(color);
  prediction -> SetLineStyle(style);
	prediction -> SetMarkerStyle(style);
  prediction -> SetMarkerColor(color);
  prediction -> SetLineWidth(2);

  prediction -> SetTitle(title.c_str());
	//prediction -> SetTitleFont(13);
	prediction -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
	prediction -> GetYaxis()->SetTitle(GetAxisLabel(observable,1).c_str());
  prediction -> GetXaxis()->CenterTitle();
  prediction -> GetYaxis()->CenterTitle();

  if( y_max == 0 ) y_max = (prediction -> GetMaximum()) * ( 1+0.2 );
	int FontStyle = 132;
  prediction->GetXaxis()->SetTitleOffset(0.8);
	prediction->GetXaxis()->SetLabelSize(0.04);
  prediction->GetXaxis()->SetTitleSize(0.06);
  prediction->GetXaxis()->SetNdivisions(6);
	prediction->GetXaxis()->SetLabelFont(FontStyle);
	prediction->GetXaxis()->SetTitleFont(FontStyle);

  prediction->GetYaxis()->SetNdivisions(6);
  prediction->GetYaxis()->SetTitleOffset(0.8);
  prediction->GetYaxis()->SetLabelSize(0.04);
  prediction->GetYaxis()->SetTitleSize(0.06);
	prediction->GetYaxis()->SetLabelFont(43);
	prediction->GetYaxis()->SetLabelFont(FontStyle);
	prediction->GetYaxis()->SetTitleFont(FontStyle);
  prediction->GetYaxis()->SetRangeUser(0,y_max);
  prediction->GetYaxis()->SetMaxDigits(1) ;
	prediction->SetTitleFont(FontStyle);

  return;
}

std::string compute_acceptance(std::vector<std::string> mc_files, std::string observable, std::string title,
                               std::string input_MC_location, std::string output_location,  std::string output_file_name,
                               int id_sector = -1 /*all*/ ) {

	// Define trees
  std::vector<TFile*> files_mcrecoacc, files_mctrueacc;
  std::vector<TTree*> trees_mcrecoacc, trees_mctrueacc ;
	// Define Hists
  std::vector<TH1D*>  hists_recoacc, hists_trueacc, hists_recoacc_0, hists_trueacc_0, hists_recoacc_1, hists_trueacc_1,
                      hists_recoacc_2, hists_trueacc_2, hists_recoacc_3, hists_trueacc_3,
                      hists_recoacc_4, hists_trueacc_4, hists_recoacc_5, hists_trueacc_5 ;
	std::vector<std::vector<TH1D*>> hists_recoacc_slices, hists_trueacc_slices;
  std::vector<TTree*> trees;
  std::vector<TH1D*>  hists, ratios, ratios_0, ratios_1, ratios_2, ratios_3, ratios_4, ratios_5 ;
	std::vector<std::vector<TH1D*>> ratios_slices ;
  std::vector<double> binning ;
	// Get energy from tree to define range
	double BeamE ;

  for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
    files_mcrecoacc.push_back(new TFile((input_MC_location+mc_files[i]+"_truereco.root").c_str(),"ROOT"));
    files_mctrueacc.push_back(new TFile((input_MC_location+mc_files[i]+"_true.root").c_str(),"ROOT"));
    if( !files_mcrecoacc[i] ) { std::cout << "ERROR: the "<< mc_files[i] << "_truereco.root does not exist." <<std::endl; return "";}
    if( !files_mctrueacc[i] ) { std::cout << "ERROR: the "<< mc_files[i] << "_true.root  does not exist." <<std::endl; return "";}
    trees_mcrecoacc.push_back( (TTree*)files_mcrecoacc[i]->Get("MCCLAS6Tree"));
    trees_mctrueacc.push_back( (TTree*)files_mctrueacc[i]->Get("MCCLAS6Tree"));
    if( !trees_mctrueacc[i] || !trees_mcrecoacc[i] ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return "";}

		trees_mctrueacc[0]->SetBranchAddress("BeamE",&BeamE);
		trees_mctrueacc[0]->GetEntry(0);
		binning = GetBinning(observable,BeamE);

    hists_recoacc.push_back( new TH1D( ("Reco MC ACC Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc.push_back( new TH1D( ("True MC ACC Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_0.push_back( new TH1D( ("Reco MC ACC Sector  0 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_0.push_back( new TH1D( ("True MC ACC Sector  0 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_1.push_back( new TH1D( ("Reco MC ACC Sector  1 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_1.push_back( new TH1D( ("True MC ACC Sector  1 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_2.push_back( new TH1D( ("Reco MC ACC Sector  2 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_2.push_back( new TH1D( ("True MC ACC Sector  2 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_3.push_back( new TH1D( ("Reco MC ACC Sector  3 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_3.push_back( new TH1D( ("True MC ACC Sector  3 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_4.push_back( new TH1D( ("Reco MC ACC Sector  4 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_4.push_back( new TH1D( ("True MC ACC Sector  4 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_5.push_back( new TH1D( ("Reco MC ACC Sector  5 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_5.push_back( new TH1D( ("True MC ACC Sector  5 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;

		std::vector<double> addbinning = GetAdditionalBinning( GetAlternativeObs(observable), BeamE );
		if ( addbinning.size() > 0 ) {
			// Adding additional histograms for slices calculation
			std::vector<TH1D*> temp_reco_slices, temp_true_slices ;
			for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
				std::string name = "MC Acceptance for " ;
				if ( k == 0 ) name += GetAlternativeObs(observable) + " < " + std::to_string(addbinning[k+1]) ;
				else if ( k == addbinning.size()-2 ) name += GetAlternativeObs(observable) + " > " + std::to_string(addbinning[k]) ;
				else name += std::to_string(addbinning[k]) + " < " +  GetAlternativeObs(observable) + " < " + std::to_string(addbinning[k]) ;
				name += ". Model "+std::to_string(i);
				temp_reco_slices.push_back( new TH1D( ("Reco " + name).c_str(), "", binning.size()-1, &binning[0] ) ) ;
				temp_true_slices.push_back( new TH1D( ("True " + name).c_str(), "", binning.size()-1, &binning[0] ) ) ;
			}
			hists_recoacc_slices.push_back(temp_reco_slices);
			hists_trueacc_slices.push_back(temp_true_slices);
		}

    unsigned int initial_size_trees = trees.size();
    unsigned int initial_size_hists = hists.size();
    trees.push_back(trees_mcrecoacc[i]);
    trees.push_back(trees_mctrueacc[i]);
    hists.push_back(hists_recoacc[i]);   //0
    hists.push_back(hists_trueacc[i]);   //1
    hists.push_back(hists_recoacc_0[i]); //2
    hists.push_back(hists_trueacc_0[i]); //3
    hists.push_back(hists_recoacc_1[i]); //4
    hists.push_back(hists_trueacc_1[i]); //5
    hists.push_back(hists_recoacc_2[i]); //6
    hists.push_back(hists_trueacc_2[i]); //7
    hists.push_back(hists_recoacc_3[i]); //8
    hists.push_back(hists_trueacc_3[i]); //9
    hists.push_back(hists_recoacc_4[i]); //10
    hists.push_back(hists_trueacc_4[i]); //11
    hists.push_back(hists_recoacc_5[i]); //12
    hists.push_back(hists_trueacc_5[i]); //13

    // Set condition for new hists
		for( unsigned int id = initial_size_hists ; id < hists.size(); id++ ){
			hists[id] -> Sumw2() ;
		}

    // OBSERVABLE DEFINITION:
    double TotWeight ;
    double ECal,Recoq3,RecoW;
    double pfl,pfl_theta,pfl_phi;
    double proton_mom,proton_phi,proton_theta;
    double pim_mom,pim_theta,pim_phi;
    double pip_mom,pip_theta,pip_phi;
    double HadAlphaT, HadDeltaPT, HadDeltaPhiT ;
    double AlphaT, DeltaPT, DeltaPhiT ;
    double RecoXBJK, RecoEnergyTransfer, RecoQ2, HadSystemMass, RecoQELEnu ;
    long NEntries ;
    bool IsBkg ;
    int ElectronSector ;
    for ( unsigned int j = initial_size_trees ; j < trees.size() ; ++j ){
      NEntries = trees[j] -> GetEntries() ;
      trees[j] -> SetBranchAddress("TotWeight",&TotWeight);
      trees[j] -> SetBranchAddress("IsBkg",&IsBkg);
      trees[j] -> SetBranchAddress("ECal",&ECal);
      trees[j] -> SetBranchAddress("pfl_theta",&pfl_theta);
      trees[j] -> SetBranchAddress("pfl_phi",&pfl_phi);
      trees[j] -> SetBranchAddress("pfl",&pfl);
      trees[j] -> SetBranchAddress("proton_mom",&proton_mom);
      trees[j] -> SetBranchAddress("proton_theta",&proton_theta);
      trees[j] -> SetBranchAddress("proton_phi",&proton_phi);
      trees[j] -> SetBranchAddress("pim_mom",&pim_mom);
      trees[j] -> SetBranchAddress("pim_theta",&pim_theta);
      trees[j] -> SetBranchAddress("pim_phi",&pim_phi);
      trees[j] -> SetBranchAddress("pip_mom",&pip_mom);
      trees[j] -> SetBranchAddress("pip_theta",&pip_theta);
      trees[j] -> SetBranchAddress("pip_phi",&pip_phi);
      trees[j] -> SetBranchAddress("RecoW",&RecoW);
			trees[j] -> SetBranchAddress("RecoQELEnu",&RecoQELEnu);
      trees[j] -> SetBranchAddress("Recoq3",&Recoq3);
      trees[j] -> SetBranchAddress("RecoXBJK",&RecoXBJK);
      trees[j] -> SetBranchAddress("RecoQ2",&RecoQ2);
      trees[j] -> SetBranchAddress("RecoEnergyTransfer",&RecoEnergyTransfer);
      trees[j] -> SetBranchAddress("AlphaT",&AlphaT);
      trees[j] -> SetBranchAddress("HadAlphaT",&HadAlphaT);
      trees[j] -> SetBranchAddress("DeltaPT",&DeltaPT);
      trees[j] -> SetBranchAddress("HadDeltaPT",&HadDeltaPT);
      trees[j] -> SetBranchAddress("DeltaPhiT",&DeltaPhiT);
      trees[j] -> SetBranchAddress("HadDeltaPhiT",&HadDeltaPhiT);
      trees[j] -> SetBranchAddress("ElectronSector",&ElectronSector);
      trees[j] -> SetBranchAddress("HadSystemMass", &HadSystemMass);
      for( int k = 0 ; k < NEntries; ++k ) {
        trees[j]->GetEntry(k) ;
        double content = 0 ;
        double w = TotWeight ;

        if( observable == "ECal") content = ECal ;
        else if ( observable == "pfl") content = pfl ;
        else if ( observable == "pfl_theta") content = pfl_theta ;
        else if ( observable == "pfl_phi") content = pfl_phi ;
        else if ( observable == "proton_mom") content = proton_mom ;
        else if ( observable == "proton_theta") content = proton_theta  ;
        else if ( observable == "proton_phi") content = proton_phi  ;
        else if ( observable == "pim_mom") content = pim_mom ;
        else if ( observable == "pim_theta") content = pim_theta  ;
        else if ( observable == "pim_phi") content = pim_phi  ;
        else if ( observable == "pip_mom") content = pip_mom ;
        else if ( observable == "pip_theta") content = pip_theta  ;
        else if ( observable == "pip_phi") content = pip_phi  ;
        else if ( observable == "RecoW") content = RecoW ;
        else if ( observable == "Recoq3") content = Recoq3 ;
				else if ( observable == "RecoQELEnu") content = RecoQELEnu ;
        else if ( observable == "RecoXBJK") content = RecoXBJK ;
        else if ( observable == "RecoQ2") content = RecoQ2 ;
        else if ( observable == "RecoEnergyTransfer") content = RecoEnergyTransfer ;
        else if ( observable == "AlphaT") content = AlphaT ;
        else if ( observable == "HadAlphaT") content = HadAlphaT ;
        else if ( observable == "DeltaPT") content = DeltaPT ;
        else if ( observable == "HadDeltaPT") content = HadDeltaPT ;
        else if ( observable == "DeltaPhiT") content = DeltaPhiT ;
        else if ( observable == "HadDeltaPhiT") content = HadDeltaPhiT ;
        else if ( observable == "HadSystemMass") content = HadSystemMass ;
        // Fill the per Sector  histogram
        hists[2*(ElectronSector+1)+(j-initial_size_trees)+initial_size_hists] -> Fill( content, w ) ;
        hists[2*(ElectronSector+1)+(j-initial_size_trees)+initial_size_hists] -> SetLineWidth(3);

        if( id_sector > 0 ) {
          // Compute only for Sector  of interest
          if( id_sector != ElectronSector ) continue ;
        }

				hists[j+initial_size_hists-initial_size_trees] -> Fill( content, w ) ;
        hists[j+initial_size_hists-initial_size_trees] -> SetLineWidth(3);

				std::string alt_obs = GetAlternativeObs(observable) ;
				double content_2 ;
				if ( alt_obs == "ECal" ) content_2 = ECal ;
				else if ( alt_obs == "HadAlphaT" ) content_2 = HadAlphaT ;
				else if ( alt_obs == "HadDeltaPT" ) content_2 = HadDeltaPT ;

        // Fill sliced histogram
				if( addbinning.size() != 0 ) {
					for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
						if( content_2 > addbinning[l] && content_2 < addbinning[l+1] ){
							if( j == initial_size_trees ) hists_recoacc_slices[i][l] -> Fill( content, w ) ;
							else if( j == initial_size_trees + 1 ) hists_trueacc_slices[i][l] -> Fill( content, w ) ;
						}
					}
      	}
			}
    }

    ratios.push_back( (TH1D*)hists_trueacc[i]->Clone() ) ;
    ratios[i] -> Divide( hists_recoacc[i] );
    ratios[i] -> SetName(("Acceptance_model_"+std::to_string(i)).c_str());
    StandardFormat( ratios[i], title, kBlack+i+1, 2+i, observable ) ;
    ratios[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios[i] -> GetYaxis()->SetTitle("Acceptance correction");

    ratios_0.push_back( (TH1D*)hists_trueacc_0[i]->Clone() ) ;
    ratios_0[i] -> Divide( hists_recoacc_0[i] );
    ratios_0[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_0").c_str());
    StandardFormat( ratios_0[i], title, kOrange+1+i, 2+i, observable ) ;
    ratios_0[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_0[i] -> GetYaxis()->SetTitle("Acceptance correction e-Sector  0");

    ratios_1.push_back( (TH1D*)hists_trueacc_1[i]->Clone() ) ;
    ratios_1[i] -> Divide( hists_recoacc_1[i] );
    ratios_1[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_1").c_str());
    StandardFormat( ratios_1[i], title, kPink+4-i, 2+i, observable ) ;
    ratios_1[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_1[i] -> GetYaxis()->SetTitle("Acceptance correction e-Sector 1");

    ratios_2.push_back( (TH1D*)hists_trueacc_2[i]->Clone() ) ;
    ratios_2[i] -> Divide( hists_recoacc_2[i] );
    ratios_2[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_2").c_str());
    StandardFormat( ratios_2[i], title, kViolet+5-i, 2+i, observable ) ;
    ratios_2[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_2[i] -> GetYaxis()->SetTitle("Acceptance correction e-Sector 2");

    ratios_3.push_back( (TH1D*)hists_trueacc_3[i]->Clone() ) ;
    ratios_3[i] -> Divide( hists_recoacc_3[i] );
    ratios_3[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_3").c_str());
    StandardFormat( ratios_3[i], title, kAzure-5+i, 2+i, observable ) ;
    ratios_3[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_3[i] -> GetYaxis()->SetTitle("Acceptance correction e-Sector 3");

    ratios_4.push_back( (TH1D*)hists_trueacc_4[i]->Clone() ) ;
    ratios_4[i] -> Divide( hists_recoacc_4[i] );
    ratios_4[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_4").c_str());
    StandardFormat( ratios_4[i], title, kTeal-7-i, 2+i, observable ) ;
    ratios_4[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_4[i] -> GetYaxis()->SetTitle("Acceptance correction e-Sector 4");

    ratios_5.push_back( (TH1D*)hists_trueacc_5[i]->Clone() ) ;
    ratios_5[i] -> Divide( hists_recoacc_5[i] );
    ratios_5[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_5").c_str());
    StandardFormat( ratios_5[i], title, kGreen-3-i, 2+i, observable ) ;
    ratios_5[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_5[i] -> GetYaxis()->SetTitle("Acceptance correction e-Sector 5");

		std::vector<TH1D*> temp_ratios_slices ;
		if( hists_trueacc_slices.size() != 0 ) {
			for( unsigned int l = 0 ; l < hists_trueacc_slices[i].size() ; ++l ){
					temp_ratios_slices.push_back( (TH1D*)hists_trueacc_slices[i][l]->Clone() );
					temp_ratios_slices[l] -> Divide( hists_recoacc_slices[i][l] );
					StandardFormat( temp_ratios_slices[l], title, kGreen-3-i, 2+i, observable ) ;
					std::string name = "Acceptance for slice " + std::to_string(l) ;
		    	temp_ratios_slices[l] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
		    	temp_ratios_slices[l] -> GetYaxis()->SetTitle(name.c_str());
			}
		}
		ratios_slices.push_back(temp_ratios_slices);
  }

  TH1D* ratio = (TH1D*)ratios[0]->Clone();
  ratio -> SetName("Acceptance");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio->Add(ratios[i]);
  }
  ratio -> Scale( 1./mc_files.size() );
  StandardFormat( ratio, title, kBlack, 1, observable ) ;
  ratio -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio -> GetYaxis()->SetTitle("Acceptance correction");

  TH1D* ratio_0 = (TH1D*)ratios_0[0]->Clone();
  ratio_0 -> SetName("Acceptance_0");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_0->Add(ratios_0[i]);
  }
  ratio_0 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_0, title, kOrange+1, 1, observable ) ;
  ratio_0 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_0 -> GetYaxis()->SetTitle("Acceptance correction e-Sector 0");

  TH1D* ratio_1 = (TH1D*)ratios_1[0]->Clone();
  ratio_1 -> SetName("Acceptance_1");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_1->Add(ratios_1[i]);
  }
  ratio_1 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_1, title, kPink+4, 1, observable ) ;
  ratio_1 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_1 -> GetYaxis()->SetTitle("Acceptance correction e-Sector 1");

  TH1D* ratio_2 = (TH1D*)ratios_2[0]->Clone();
  ratio_2 -> SetName("Acceptance_2");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_2->Add(ratios_2[i]);
  }
  ratio_2 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_2, title, kViolet+5, 1, observable ) ;
  ratio_2 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_2 -> GetYaxis()->SetTitle("Acceptance correction e-Sector 2");

  TH1D* ratio_3 = (TH1D*)ratios_3[0]->Clone();
  ratio_3 -> SetName("Acceptance_3");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_3->Add(ratios_3[i]);
  }
  ratio_3 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_3, title, kAzure-5, 1, observable ) ;
  ratio_3 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_3 -> GetYaxis()->SetTitle("Acceptance correction e-Sector 3");

  TH1D* ratio_4 = (TH1D*)ratios_4[0]->Clone();
  ratio_4 -> SetName("Acceptance_4");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_4->Add(ratios_4[i]);
  }
  ratio_4 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_4, title, kTeal-7, 1, observable ) ;
  ratio_4 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_4 -> GetYaxis()->SetTitle("Acceptance correction e-Sector 4");

  TH1D* ratio_5 = (TH1D*)ratios_5[0]->Clone();
  ratio_5 -> SetName("Acceptance_5");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_5->Add(ratios_5[i]);
  }
  ratio_5 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_5, title, kGreen-3, 1, observable ) ;
  ratio_5 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_5 -> GetYaxis()->SetTitle("Acceptance correction e-Sector 5");

	std::vector<double> addbinning = GetAdditionalBinning( GetAlternativeObs(observable), BeamE ) ;
	if ( addbinning.size() > 0 ) {
		// Adding additional histograms for slices calculation
		std::vector<TH1D*> temp_reco_slices, temp_true_slices ;
		for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
			std::string name = "MC Acceptance for " ;
			if ( k == 0 ) name += GetAlternativeObs(observable) + " < " + std::to_string(addbinning[k+1]) ;
			else if ( k == addbinning.size()-2 ) name += GetAlternativeObs(observable) + " > " + std::to_string(addbinning[k]) ;
			else name += std::to_string(addbinning[k]) + " < " +  GetAlternativeObs(observable) + " < " + std::to_string(addbinning[k]) ;
		}
	}

	std::vector<TH1D*> ratio_slices ;
	if( ratios_slices.size() != 0 ) {
		for( unsigned l = 0 ; l < ratios_slices[0].size(); ++l ){
			std::string name = "MC Acceptance for " ;
			if ( l == 0 ) name += GetAlternativeObs(observable) + " < " + std::to_string(addbinning[l+1]) ;
			else if ( l == addbinning.size()-2 ) name += GetAlternativeObs(observable) + " > " + std::to_string(addbinning[l]) ;
			else name += std::to_string(addbinning[l]) + " < " +  GetAlternativeObs(observable) + " < " + std::to_string(addbinning[l]) ;

			TH1D * temp_slice_ratio = (TH1D*)ratios_slices[0][l]->Clone();
			temp_slice_ratio -> SetName(("Acceptance_Slice_"+std::to_string(l)).c_str());
			for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
			 	temp_slice_ratio->Add(ratios_slices[i][l]);
			}
			temp_slice_ratio -> Scale( 1./mc_files.size() );
			StandardFormat( temp_slice_ratio, title, kGreen-3, 1, observable ) ;
			temp_slice_ratio -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
			temp_slice_ratio -> GetYaxis()->SetTitle("Acceptance");
			temp_slice_ratio ->SetTitle((name+". Slice "+std::to_string(l)).c_str());
			ratio_slices.push_back(temp_slice_ratio);
		}
  }
  std::string output_name = output_file_name+"_acceptance_correction_"+observable ;
  if( id_sector > 0 ) output_name += "_sector_"+std::to_string(id_sector) ;
  std::string acc_file = "/AcceptanceFiles/"+output_name ;

	std::filesystem::path acceptance_path{(output_location+"/AcceptanceFiles").c_str()};
	if( ! std::filesystem::exists(acceptance_path) ) std::filesystem::create_directory(acceptance_path);

  TFile outputFile ((output_location+acc_file+".root").c_str(),"RECREATE");

  TCanvas * c_1 = new TCanvas("c_1","c_1",200,10,700,500);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

	// Store total contribution (averaged)
  ratio->Write();
  ratio_0->Write();
  ratio_1->Write();
  ratio_2->Write();
  ratio_3->Write();
  ratio_4->Write();
  ratio_5->Write();

	// Store per model
	for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
    ratios_0[i] -> Write();
    ratios_1[i] -> Write();
    ratios_2[i] -> Write();
    ratios_3[i] -> Write();
    ratios_4[i] -> Write();
    ratios_5[i] -> Write();
  }

	for ( unsigned int i = 0 ; i < ratio_slices.size(); ++i ){
		ratio_slices[i] -> Write() ;
	}

	// Plot it
	ratio->Draw("hist err");
  for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
    ratios[i]->Draw("hist err same");
    ratios[i]->Write();
  }
  ratio->Draw("hist err same");
  //teff->Draw("AP");

  c_1->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_total.root").c_str());
  c_1->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_total.pdf").c_str());
  delete c_1 ;

  // Draw total xsec per sectors
  TCanvas * c_sector_2 = new TCanvas("c_sector_2","c_sector_2",200,10,700,500);
  c_sector_2->cd();
  TPad *pad_sector = new TPad("pad1","",0,0,1,1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->Divide(3,2);

  TPad *pad_sector_0 = (TPad*)pad_sector->cd(1);
  pad_sector_0 -> cd();
  pad_sector_0 -> SetBottomMargin(0.15);
  pad_sector_0 -> SetLeftMargin(0.15);
  ratio_0 -> GetYaxis()->SetTitleOffset(1.2);
  ratio_0 -> Draw("hist err");
  for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
     ratios_0[i] -> Draw("hist err same");
	}

  TPad *pad_sector_1 = (TPad*)pad_sector->cd(2);
  pad_sector_1 -> cd();
  pad_sector_1 -> SetBottomMargin(0.15);
  pad_sector_1 -> SetLeftMargin(0.15);
  ratio_1 -> GetYaxis()->SetTitleOffset(1.2);
  ratio_1 -> Draw("hist err");
	for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
		 ratios_1[i] -> Draw("hist err same");
	}

  TPad *pad_sector_2 = (TPad*)pad_sector->cd(3);
  pad_sector_2 -> cd();
  pad_sector_2 -> SetBottomMargin(0.15);
  pad_sector_2 -> SetLeftMargin(0.15);
  ratio_2 -> GetYaxis()->SetTitleOffset(1.2);
  ratio_2 -> Draw("hist");
	for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
		 ratios_2[i] -> Draw("hist same");
	}

  TPad *pad_sector_3 = (TPad*)pad_sector->cd(4);
  pad_sector_3 -> cd();
  pad_sector_3 -> SetBottomMargin(0.15);
  pad_sector_3 -> SetLeftMargin(0.15);
  ratio_3 -> GetYaxis()->SetTitleOffset(1.2);
  ratio_3 -> Draw("hist err");
	for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
		 ratios_3[i] -> Draw("hist err same");
	}

  TPad *pad_sector_4 = (TPad*)pad_sector->cd(5);
  pad_sector_4 -> cd();
  pad_sector_4 -> SetBottomMargin(0.15);
  pad_sector_4 -> SetLeftMargin(0.15);
  ratio_4 -> GetYaxis()->SetTitleOffset(1.2);
  ratio_4 -> Draw("hist err");
	for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
		 ratios_4[i] -> Draw("hist err same");
	}

  TPad *pad_sector_5 = (TPad*)pad_sector->cd(6);
  pad_sector_5 -> cd();
  pad_sector_5 -> SetBottomMargin(0.15);
  pad_sector_5 -> SetLeftMargin(0.15);
  ratio_5 -> GetYaxis()->SetTitleOffset(1.2);
  ratio_5 -> Draw("hist err");
	for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
		 ratios_5[i] -> Draw("hist err same");
	}

  c_sector_2->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_persector.root").c_str());
  c_sector_2->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_persector.pdf").c_str());
  delete c_sector_2 ;

	for (size_t i = 0; i < files_mcrecoacc.size(); i++) {
		delete files_mcrecoacc[i] ;
		delete files_mctrueacc[i] ;
	}
	outputFile.Close();
  return acc_file ;
}

void NormalizeHist( TH1D * h, double normalization_factor ){
    // Data normalization
    h -> Scale( normalization_factor );
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

void CorrectData(TH1D* h, TH1D* acc) {

  double NBins = h->GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double content = h->GetBinContent(i);
    double error = h->GetBinError(i);
    if( h->GetBinContent(i) != 0 ) h->SetBinContent(i, content * acc->GetBinContent(i));
  }
}

// Input paramters:
// MC_file_name : true MC file, without detector effects, after e4nu analysis
// data_file_name: data file, after e4nu Analysis
// acceptance_file_name: acceptance file obtained with compute_acceptance.C
// target target_pdg
// beam energy
// Number of events in original MC file (pre-analysis)
// id_secotr: ID of the Sector  of interest

void Plot1DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
                std::string acceptance_file_name, std::string observable,
                std::string title, std::string data_name, std::vector<std::string> model,
                std::string input_MC_location, std::string input_data_location, std::string output_location,
								std::string output_file_name, int id_sector = -1 /*all*/ ) {

  TCanvas * c1 = new TCanvas("c1","c1",200,10,700,500);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  std::vector<TFile*> files_true_MC ;
  for( unsigned int id = 0 ; id < MC_files_name.size(); ++id ){
      files_true_MC.push_back(new TFile((input_MC_location+MC_files_name[id]+"_true.root").c_str(),"ROOT"));
      if( !files_true_MC[id] ) { std::cout << "ERROR: the "<< input_MC_location<<MC_files_name[id]<<"_true.root does not exist." << std::endl; return ;}
  }
  TFile * file_data = new TFile((input_data_location+data_file_name+".root").c_str(),"READ");
  TFile * file_acceptance = new TFile((output_location+acceptance_file_name+".root").c_str(),"READ");
  if( !file_data ) { std::cout << "ERROR: the "<< input_data_location << data_file_name << ".root does not exist." <<std::endl; return ;}
  if( !file_acceptance ) { std::cout << "ERROR: the "<< output_location << acceptance_file_name << ".root does not exist." <<std::endl; return ;}

  TH1D * h_acceptance = (TH1D*)file_acceptance->Get("Acceptance");
  TH1D * h_acceptance_0 = (TH1D*)file_acceptance->Get("Acceptance_0");
  TH1D * h_acceptance_1 = (TH1D*)file_acceptance->Get("Acceptance_1");
  TH1D * h_acceptance_2 = (TH1D*)file_acceptance->Get("Acceptance_2");
  TH1D * h_acceptance_3 = (TH1D*)file_acceptance->Get("Acceptance_3");
  TH1D * h_acceptance_4 = (TH1D*)file_acceptance->Get("Acceptance_4");
  TH1D * h_acceptance_5 = (TH1D*)file_acceptance->Get("Acceptance_5");

	// Get Tree for main model
  TTree * tree_true = (TTree*)files_true_MC[0]->Get("MCCLAS6Tree");
	// Get configured energy, used for plotting
	double BeamE ;
	tree_true->SetBranchAddress("BeamE",&BeamE);
	tree_true->GetEntry(0);

	// Get Acceptance for slices
	std::vector<TH1D*> h_acc_slices ;
	std::vector<double> addbinning = GetAdditionalBinning( GetAlternativeObs(observable), BeamE ) ;
	if ( addbinning.size() > 0 ) {
		for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
			h_acc_slices.push_back( (TH1D*)file_acceptance->Get(("Acceptance_Slice_"+std::to_string(k)).c_str() ) ) ;
			if( !h_acc_slices[k] ) { std::cout << "ERROR: Slice acceptance empty"<<std::endl;return;}
		}
	}

  // For submodels only total prediction is plotted
  std::vector<TTree*> tree_submodels ;
  for( unsigned int id = 1 ; id < MC_files_name.size() ; ++id ){
    tree_submodels.push_back( (TTree*)files_true_MC[id]->Get("MCCLAS6Tree") );
    if( !tree_submodels[id - 1] ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return ;}
  }

  TTree * tree_data = (TTree*)file_data->Get("CLAS6Tree");

  if( !h_acceptance ) { std::cout << "ERROR: Acceptance is not defined"<<std::endl; return ; }
  if( !tree_true || !tree_data ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return ;}

  // Create histogram for total and total xsec per sector
  TH1D * hist_true = (TH1D*) h_acceptance ->Clone();
  hist_true -> SetName( "MC True") ;
  hist_true -> Reset("ICE");
  TH1D * hist_true_0 = (TH1D*) h_acceptance_0 ->Clone();
  hist_true_0 -> SetName( "MC True Sector  0") ;
  hist_true_0 -> Reset("ICE");
  TH1D * hist_true_1 = (TH1D*) h_acceptance_1 ->Clone();
  hist_true_1 -> SetName( "MC True Sector  1") ;
  hist_true_1 -> Reset("ICE");
  TH1D * hist_true_2 = (TH1D*) h_acceptance_2 ->Clone();
  hist_true_2 -> SetName( "MC True Sector  2") ;
  hist_true_2 -> Reset("ICE");
  TH1D * hist_true_3 = (TH1D*) h_acceptance_3 ->Clone();
  hist_true_3 -> SetName( "MC True Sector  3") ;
  hist_true_3 -> Reset("ICE");
  TH1D * hist_true_4 = (TH1D*) h_acceptance_4 ->Clone();
  hist_true_4 -> SetName( "MC True Sector  4") ;
  hist_true_4 -> Reset("ICE");
  TH1D * hist_true_5 = (TH1D*) h_acceptance_5 ->Clone();
  hist_true_5 -> SetName( "MC True Sector  5") ;
  hist_true_5 -> Reset("ICE");

  // Breakdown histograms for total (all sectors only):
  TH1D * hist_true_QEL = (TH1D*) h_acceptance ->Clone();
  hist_true_QEL -> SetName( "MC True QEL") ;
  hist_true_QEL -> Reset("ICE");
  TH1D * hist_true_RES = (TH1D*) h_acceptance ->Clone();
  hist_true_RES -> SetName( "MC True QEL") ;
  hist_true_RES -> Reset("ICE");
  TH1D * hist_true_SIS = (TH1D*) h_acceptance ->Clone();
  hist_true_SIS -> SetName( "MC True QEL") ;
  hist_true_SIS -> Reset("ICE");
  TH1D * hist_true_MEC = (TH1D*) h_acceptance ->Clone();
  hist_true_MEC -> SetName( "MC True QEL") ;
  hist_true_MEC -> Reset("ICE");
  TH1D * hist_true_DIS = (TH1D*) h_acceptance ->Clone();
  hist_true_DIS -> SetName( "MC True QEL") ;
  hist_true_DIS -> Reset("ICE");

  // Same per model - only total prediction
  std::vector<TH1D*> hists_true_submodel;
  for( unsigned int id = 1 ; id < MC_files_name.size(); ++id ){
    hists_true_submodel.push_back( (TH1D*) h_acceptance ->Clone() );
    hists_true_submodel[id - 1] -> SetName( ("MC True Model "+std::to_string(id)).c_str()) ;
    hists_true_submodel[id - 1] -> Reset("ICE");
    hists_true_submodel[id - 1] -> SetLineWidth(3);
  }

  // Create hist for each slice on true
	std::vector<TH1D*> h_total_slices, h_QEL_slices, h_RES_slices, h_DIS_slices, h_MEC_slices, h_SIS_slices ;
	if ( addbinning.size() > 0 ) {
		for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
			h_total_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
			h_total_slices[k] -> SetName( ("MC True Slice "+std::to_string(k)).c_str() ) ;
		  h_total_slices[k] -> Reset("ICE");
			h_QEL_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
			h_QEL_slices[k] -> SetName( ("MC True QEL Slice "+std::to_string(k)).c_str() ) ;
		  h_QEL_slices[k] -> Reset("ICE");
			h_RES_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
			h_RES_slices[k] -> SetName( ("MC True RES Slice "+std::to_string(k)).c_str() ) ;
		  h_RES_slices[k] -> Reset("ICE");
			h_SIS_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
			h_SIS_slices[k] -> SetName( ("MC True SIS Slice "+std::to_string(k)).c_str() ) ;
		  h_SIS_slices[k] -> Reset("ICE");
			h_MEC_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
			h_MEC_slices[k] -> SetName( ("MC True MEC Slice "+std::to_string(k)).c_str() ) ;
		  h_MEC_slices[k] -> Reset("ICE");
			h_DIS_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
			h_DIS_slices[k] -> SetName( ("MC True DIS Slice "+std::to_string(k)).c_str() ) ;
		  h_DIS_slices[k] -> Reset("ICE");
		}
	}


  // total and per sector
  TH1D * hist_data = (TH1D*) h_acceptance ->Clone();
  hist_data -> SetName( "Data") ;
  hist_data -> Reset("ICE");
  TH1D * hist_data_0 = (TH1D*) h_acceptance_0 ->Clone();
  hist_data_0 -> SetName( "Data Sector  0") ;
  hist_data_0 -> Reset("ICE");
  TH1D * hist_data_1 = (TH1D*) h_acceptance_1 ->Clone();
  hist_data_1 -> SetName( "Data Sector  1") ;
  hist_data_1 -> Reset("ICE");
  TH1D * hist_data_2 = (TH1D*) h_acceptance_2 ->Clone();
  hist_data_2 -> Reset("ICE");
	hist_data_2 -> SetName( "Data Sector  2") ;
  TH1D * hist_data_3 = (TH1D*) h_acceptance_3 ->Clone();
  hist_data_3 -> SetName( "Data Sector  3") ;
  hist_data_3 -> Reset("ICE");
  TH1D * hist_data_4 = (TH1D*) h_acceptance_4 ->Clone();
  hist_data_4 -> SetName( "Data Sector  4") ;
  hist_data_4 -> Reset("ICE");
  TH1D * hist_data_5 = (TH1D*) h_acceptance_5 ->Clone();
  hist_data_5 -> SetName( "Data Sector  5") ;
  hist_data_5 -> Reset("ICE");

	// Create hist for each slice on data
	std::vector<TH1D*> h_data_slices  ;
	if ( addbinning.size() > 0 ) {
		for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
			h_data_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
			h_data_slices[k] -> SetName( ("Data Slice "+std::to_string(k)).c_str() ) ;
		  h_data_slices[k] -> Reset("ICE");
		}
	}

  std::vector<TTree*> trees = { tree_true, tree_data };
  std::vector<TH1D*> hists = { hist_true, hist_data, hist_true_0, hist_data_0,
                               hist_true_1, hist_data_1, hist_true_2, hist_data_2,
                               hist_true_3, hist_data_3, hist_true_4, hist_data_4,
                               hist_true_5, hist_data_5};

  unsigned int size_primary_trees = trees.size();
  unsigned int size_primary_hists = hists.size();
  // Adding total predictions for alternative models
  for( unsigned int id = 1 ; id < MC_files_name.size(); ++id ){
      trees.push_back( tree_submodels[id - 1] );
      hists.push_back( hists_true_submodel[id - 1] );
  }

  // OBSERVABLE DEFINITION:
  double TotWeight ;
  double ECal,Recoq3,RecoW;
  double pfl,pfl_theta,pfl_phi;
  double proton_mom,proton_phi,proton_theta;
  double pim_mom,pim_theta,pim_phi;
  double pip_mom,pip_theta,pip_phi;
  double HadAlphaT, HadDeltaPT, HadDeltaPhiT ;
  double AlphaT, DeltaPT, DeltaPhiT ;
  double RecoXBJK, RecoEnergyTransfer, RecoQ2, HadSystemMass, RecoQELEnu ;
  long NEntries ;
  bool IsBkg ;
  int ElectronSector ;
  bool QEL, RES, DIS, MEC;
  double MCNormalization, DataNormalization ;

  for ( unsigned int i = 0 ; i < trees.size() ; ++i ){

    NEntries = trees[i] -> GetEntries() ;
    trees[i] -> SetBranchAddress("TotWeight",&TotWeight);
    trees[i] -> SetBranchAddress("IsBkg",&IsBkg);
    trees[i] -> SetBranchAddress("ECal",&ECal);
    trees[i] -> SetBranchAddress("pfl_theta",&pfl_theta);
    trees[i] -> SetBranchAddress("pfl_phi",&pfl_phi);
    trees[i] -> SetBranchAddress("pfl",&pfl);
    trees[i] -> SetBranchAddress("proton_mom",&proton_mom);
    trees[i] -> SetBranchAddress("proton_theta",&proton_theta);
    trees[i] -> SetBranchAddress("proton_phi",&proton_phi);
    trees[i] -> SetBranchAddress("pim_mom",&pim_mom);
    trees[i] -> SetBranchAddress("pim_theta",&pim_theta);
    trees[i] -> SetBranchAddress("pim_phi",&pim_phi);
    trees[i] -> SetBranchAddress("pip_mom",&pip_mom);
    trees[i] -> SetBranchAddress("pip_theta",&pip_theta);
    trees[i] -> SetBranchAddress("pip_phi",&pip_phi);
    trees[i] -> SetBranchAddress("RecoW",&RecoW);
    trees[i] -> SetBranchAddress("Recoq3",&Recoq3);
		trees[i] -> SetBranchAddress("RecoQELEnu",&RecoQELEnu);
    trees[i] -> SetBranchAddress("RecoXBJK",&RecoXBJK);
    trees[i] -> SetBranchAddress("RecoQ2",&RecoQ2);
    trees[i] -> SetBranchAddress("RecoEnergyTransfer",&RecoEnergyTransfer);
    trees[i] -> SetBranchAddress("AlphaT",&AlphaT);
    trees[i] -> SetBranchAddress("HadAlphaT",&HadAlphaT);
    trees[i] -> SetBranchAddress("DeltaPT",&DeltaPT);
    trees[i] -> SetBranchAddress("HadDeltaPT",&HadDeltaPT);
    trees[i] -> SetBranchAddress("DeltaPhiT",&DeltaPhiT);
    trees[i] -> SetBranchAddress("HadDeltaPhiT",&HadDeltaPhiT);
    trees[i] -> SetBranchAddress("ElectronSector",&ElectronSector);
		trees[i] -> SetBranchAddress("HadSystemMass", &HadSystemMass);
		// Only fill true info for the first model:
    if( i == 0 ) trees[i] -> SetBranchAddress("QEL",&QEL);
    if( i == 0 ) trees[i] -> SetBranchAddress("RES",&RES);
    if( i == 0 ) trees[i] -> SetBranchAddress("MEC",&MEC);
    if( i == 0 ) trees[i] -> SetBranchAddress("DIS",&DIS);

    // Only second tree corresponds to data
    if( i != 1 ) trees[i] -> SetBranchAddress("MCNormalization", &MCNormalization );
    else trees[i] -> SetBranchAddress("DataNormalization",&DataNormalization );

    for( int j = 0 ; j < NEntries ; ++j ) {
      trees[i]->GetEntry(j) ;
      double content = 0 ;
      double w = TotWeight ;

      if( observable == "ECal") content = ECal ;
      else if ( observable == "pfl") content = pfl ;
      else if ( observable == "pfl_theta") content = pfl_theta ;
      else if ( observable == "pfl_phi") content = pfl_phi ;
      else if ( observable == "proton_mom") content = proton_mom ;
      else if ( observable == "proton_theta") content = proton_theta ;
      else if ( observable == "proton_phi") content = proton_phi ;
      else if ( observable == "pim_mom") content = pim_mom ;
      else if ( observable == "pim_theta") content = pim_theta ;
      else if ( observable == "pim_phi") content = pim_phi ;
      else if ( observable == "pip_mom") content = pip_mom ;
      else if ( observable == "pip_theta") content = pip_theta ;
      else if ( observable == "pip_phi") content = pip_phi ;
      else if ( observable == "RecoW") content = RecoW ;
      else if ( observable == "Recoq3") content = Recoq3 ;
			else if ( observable == "RecoQELEnu") content = RecoQELEnu ;
      else if ( observable == "RecoXBJK") content = RecoXBJK ;
      else if ( observable == "RecoQ2") content = RecoQ2 ;
      else if ( observable == "RecoEnergyTransfer") content = RecoEnergyTransfer ;
      else if ( observable == "AlphaT") content = AlphaT ;
      else if ( observable == "HadAlphaT") content = HadAlphaT ;
      else if ( observable == "DeltaPT") content = DeltaPT ;
      else if ( observable == "HadDeltaPT") content = HadDeltaPT ;
      else if ( observable == "DeltaPhiT") content = DeltaPhiT ;
      else if ( observable == "HadDeltaPhiT") content = HadDeltaPhiT ;
			else if ( observable == "HadSystemMass") content = HadSystemMass ;
      unsigned int id_hist = i ;
      // Fill the per Sector  histogram. Only for primary model
      if( i < size_primary_trees ) hists[size_primary_trees*(ElectronSector+1)+i] -> Fill( content, w ) ;
      if( i > size_primary_trees - 1 ) id_hist = size_primary_hists + ( i - size_primary_trees );

      if( id_sector > 0 ) {
      	// Compute only for Sector  of interest
	      if( id_sector != ElectronSector ) continue ;
      }
      hists[id_hist] -> Fill( content, w ) ;
      hists[id_hist] -> SetLineWidth(3);

      if( i == 0 ){
        if( QEL ) hist_true_QEL -> Fill( content, w ) ;
        if( RES ) {
          hist_true_RES -> Fill( content, w ) ;
        }
        if( DIS ) {
          if( RecoW < 1.7 ) hist_true_SIS -> Fill( content, w ) ;
          else hist_true_DIS -> Fill( content, w ) ;
        }
        if( MEC ) hist_true_MEC -> Fill( content, w ) ;
      }

			// Fill slices
			if( addbinning.size() != 0 ) {
				std::string alt_obs = GetAlternativeObs(observable) ;
				double content_2 ;
				if ( alt_obs == "ECal" ) content_2 = ECal ;
				else if ( alt_obs == "HadAlphaT" ) content_2 = HadAlphaT ;
				else if ( alt_obs == "HadDeltaPT" ) content_2 = HadDeltaPT ;

				for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
					if( content_2 > addbinning[l] && content_2 < addbinning[l+1] ){
						if( i == 0 /* MC */ ) {
							h_total_slices[l] -> Fill( content, w ) ;
							// Fill also breakdown for slice
							if( QEL ) h_QEL_slices[l] -> Fill( content, w ) ;
			        if( RES ) {
			          h_RES_slices[l] -> Fill( content, w ) ;
			        }
			        if( DIS ) {
			          if( RecoW < 1.7 ) h_SIS_slices[l] -> Fill( content, w ) ;
			          else h_DIS_slices[l] -> Fill( content, w ) ;
			        }
			        if( MEC ) h_MEC_slices[l] -> Fill( content, w ) ;
						} else if( i == 1 /* Data */ ) h_data_slices[l] -> Fill( content, w ) ;
					}
				}
			}
    }
  }

  // Normalize data
  NormalizeHist(hist_data, DataNormalization );
  NormalizeHist(hist_data_0, DataNormalization );
  NormalizeHist(hist_data_1, DataNormalization );
  NormalizeHist(hist_data_2, DataNormalization );
  NormalizeHist(hist_data_3, DataNormalization );
  NormalizeHist(hist_data_4, DataNormalization );
  NormalizeHist(hist_data_5, DataNormalization );

  // Store uncorrected data
  TH1D * hist_data_uncorr = (TH1D*) hist_data ->Clone();
  hist_data_uncorr -> SetName( "Uncorrected Data") ;
  TH1D * hist_data_uncorr_0 = (TH1D*) hist_data_0 ->Clone();
  hist_data_uncorr_0 -> SetName( "Uncorrected Data Sector  0") ;
  TH1D * hist_data_uncorr_1 = (TH1D*) hist_data_1 ->Clone();
  hist_data_uncorr_1 -> SetName( "Uncorrected Data Sector  1") ;
  TH1D * hist_data_uncorr_2 = (TH1D*) hist_data_2 ->Clone();
  hist_data_uncorr_2 -> SetName( "Uncorrected Data Sector  2") ;
  TH1D * hist_data_uncorr_3 = (TH1D*) hist_data_3 ->Clone();
  hist_data_uncorr_3 -> SetName( "Uncorrected Data Sector  3") ;
  TH1D * hist_data_uncorr_4 = (TH1D*) hist_data_4 ->Clone();
  hist_data_uncorr_4 -> SetName( "Uncorrected Data Sector  4") ;
  TH1D * hist_data_uncorr_5 = (TH1D*) hist_data_5 ->Clone();
  hist_data_uncorr_5 -> SetName( "Uncorrected Data Sector  5") ;

  // Correct data for detector acceptance :
  CorrectData(hist_data, h_acceptance);
  CorrectData(hist_data_0, h_acceptance_0);
  CorrectData(hist_data_1, h_acceptance_1);
  CorrectData(hist_data_2, h_acceptance_2);
  CorrectData(hist_data_3, h_acceptance_3);
  CorrectData(hist_data_4, h_acceptance_4);
  CorrectData(hist_data_5, h_acceptance_5);

	// Normalize data from slices
	if( addbinning.size() != 0 ) {
		for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
			NormalizeHist(h_data_slices[l], DataNormalization );
			CorrectData(h_data_slices[l], h_acc_slices[l] );
		}
	}

  // Normalize MC
  NormalizeHist(hist_true, MCNormalization);
  NormalizeHist(hist_true_0, MCNormalization);
  NormalizeHist(hist_true_1, MCNormalization);
  NormalizeHist(hist_true_2, MCNormalization);
  NormalizeHist(hist_true_3, MCNormalization);
  NormalizeHist(hist_true_4, MCNormalization);
  NormalizeHist(hist_true_5, MCNormalization);
  NormalizeHist(hist_true_QEL, MCNormalization);
  NormalizeHist(hist_true_RES, MCNormalization);
  NormalizeHist(hist_true_SIS, MCNormalization);
  NormalizeHist(hist_true_MEC, MCNormalization);
  NormalizeHist(hist_true_DIS, MCNormalization);

  for( unsigned int id = 0 ; id < hists_true_submodel.size() ; ++id ){
    NormalizeHist( hists_true_submodel[id], MCNormalization);
		StandardFormat( hists_true_submodel[id], title, kBlack, 2+id, observable ) ;
  }

	// Normalize true from slices
	if( addbinning.size() != 0 ) {
		for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
			NormalizeHist(h_total_slices[l], MCNormalization );
			NormalizeHist(h_QEL_slices[l], MCNormalization );
			NormalizeHist(h_RES_slices[l], MCNormalization );
			NormalizeHist(h_SIS_slices[l], MCNormalization );
			NormalizeHist(h_DIS_slices[l], MCNormalization );
			NormalizeHist(h_MEC_slices[l], MCNormalization );

			std::vector<TH1D*> all_slices{h_total_slices[l],h_data_slices[l]};
			double y_max_total = GetMaximum( all_slices );

      // Add Slice information in title
			std::string title_subname = title ;
			std::string alt_obs = GetAlternativeObs(observable);
			if ( l == 0 ) {
        std::ostringstream o1 ;
				o1 << std::fixed<< std::setprecision(1) << addbinning[l+1] ;
				title_subname += " " + GetObsName(alt_obs) + "<" + o1.str() +" "+GetUnit(alt_obs) ;
			} else if ( l == addbinning.size()-2 ) {
				std::ostringstream o1 ;
				o1 << std::fixed<< std::setprecision(1) << addbinning[l] ;
				title_subname += " " + GetObsName(alt_obs) + ">" + o1.str() +" "+GetUnit(alt_obs) ;
			} else {
				std::ostringstream o1, o2 ;
				o1 << std::fixed<< std::setprecision(1) << addbinning[l] ;
				o2 << std::fixed<< std::setprecision(1) << addbinning[l+1] ;
				title_subname += " " + o1.str() + "<"+ GetObsName(alt_obs) + "<" + o2.str()+" "+GetUnit(alt_obs) ;
      }

			StandardFormat( h_total_slices[l], title_subname, kBlack, 1, observable, y_max_total ) ;
			StandardFormat( h_QEL_slices[l], title_subname, kBlue-3, 1, observable, y_max_total ) ;
			StandardFormat( h_RES_slices[l], title_subname, kGreen+2, 1, observable, y_max_total ) ;
			StandardFormat( h_SIS_slices[l], title_subname, kOrange, 1, observable, y_max_total ) ;
			StandardFormat( h_MEC_slices[l], title_subname, kMagenta-3, 1, observable, y_max_total ) ;
			StandardFormat( h_DIS_slices[l], title_subname, kCyan+1, 1, observable, y_max_total ) ;
			StandardFormat( h_data_slices[l], title_subname, kBlack, 8, observable, y_max_total ) ;

		}
	}

  // Find absolute y max
	std::vector<TH1D*> temp_check = {hist_true,hist_data};
	for( unsigned int id = 0 ; id < hists_true_submodel.size() ; ++id ){
		temp_check.push_back(hists_true_submodel[id]);
	}
	double y_max_total = GetMaximum( temp_check );
  // Format plots
  StandardFormat( hist_data, title, kBlack, 8, observable, y_max_total ) ;
  StandardFormat( hist_data_0, title+" Sector  0", kOrange+1, 8, observable ) ;
  StandardFormat( hist_data_1, title+" Sector  1", kPink+4, 8, observable ) ;
  StandardFormat( hist_data_2, title+" Sector  2", kViolet+5, 8, observable ) ;
  StandardFormat( hist_data_3, title+" Sector  3", kAzure-5, 8, observable ) ;
  StandardFormat( hist_data_4, title+" Sector  4", kTeal-7, 8, observable ) ;
  StandardFormat( hist_data_5, title+" Sector  5", kGreen-3, 8, observable ) ;
  hist_data -> SetLineWidth(0);

  StandardFormat( hist_data_uncorr, title, kRed, 8, observable, y_max_total ) ;
  StandardFormat( hist_data_uncorr_0, title+" Sector  0", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_1, title+" Sector  1", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_2, title+" Sector  2", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_3, title+" Sector  3", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_4, title+" Sector  4", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_5, title+" Sector  5", kRed, 8, observable ) ;

  StandardFormat( hist_true, title, kBlack, 1, observable, y_max_total ) ;
  StandardFormat( hist_true_0, title+" Sector  0", kBlack, 1, observable ) ;
  StandardFormat( hist_true_1, title+" Sector  1", kBlack, 1, observable ) ;
  StandardFormat( hist_true_2, title+" Sector  2", kBlack, 1, observable ) ;
  StandardFormat( hist_true_3, title+" Sector  3", kBlack, 1, observable ) ;
  StandardFormat( hist_true_4, title+" Sector  4", kBlack, 1, observable ) ;
  StandardFormat( hist_true_5, title+" Sector  5", kBlack, 1, observable ) ;

	StandardFormat( hist_true_QEL, title, kBlue-3, 1, observable ) ;
	StandardFormat( hist_true_RES, title, kGreen+2, 1, observable ) ;
	StandardFormat( hist_true_SIS, title, kOrange, 1, observable ) ;
  StandardFormat( hist_true_MEC, title, kMagenta-3, 1, observable ) ;
	StandardFormat( hist_true_DIS, title, kCyan+1, 1, observable ) ;

  // Draw total xsec (all sectors):
  hist_true -> Draw("hist");
  hist_true_QEL -> Draw("hist same");
  hist_true_RES -> Draw("hist same");
  hist_true_SIS -> Draw("hist same");
  hist_true_MEC -> Draw("hist same");
  hist_true_DIS -> Draw("hist same");
  for( unsigned int id = 0 ; id < hists_true_submodel.size() ; ++id ){
    hists_true_submodel[id] -> SetLineWidth(3);
    hists_true_submodel[id] -> Draw("hist same");
  }
  hist_data -> Draw(" err same ");

  std::string output_name = output_file_name+"_dxsec_d"+observable ;
  if( id_sector > 0 ) output_name += "_sector_"+std::to_string(id_sector) ;
	std::filesystem::path totalxsec_path{(output_location+"/TotalXSec/").c_str()};
	if( ! std::filesystem::exists(totalxsec_path) ) std::filesystem::create_directory(totalxsec_path);
  c1->SaveAs((output_location+"/TotalXSec/"+output_name+".root").c_str());
  c1->SaveAs((output_location+"/TotalXSec/"+output_name+".pdf").c_str());
  delete c1 ;

	TCanvas* c_slices = new TCanvas("c_slices","c_slices",200,10,700,500);
	c_slices->cd();
  TPad *pad_slices = new TPad("pad1","",0,0,1,1);
  pad_slices->Draw();
  pad_slices->cd();
  pad_slices->SetBottomMargin(0.15);
  pad_slices->SetLeftMargin(0.15);

	// Normalize true from slices
	if( addbinning.size() != 0 ) {
		pad_slices->Divide(addbinning.size()-1,0);
		for( unsigned int l = 0 ; l < addbinning.size()-1 ; ++l ){
			TPad *pad_slice_i = (TPad*)pad_slices->cd(1+l);
		  pad_slice_i -> cd();
		  pad_slice_i -> SetBottomMargin(0.15);
		  pad_slice_i -> SetLeftMargin(0.15);
		  h_total_slices[l] -> GetYaxis()->SetTitleOffset(1.2);

			h_total_slices[l]->Draw("hist");
			h_QEL_slices[l]->Draw("hist same ");
			h_RES_slices[l]->Draw("hist same ");
			h_SIS_slices[l]->Draw("hist same ");
			h_MEC_slices[l]->Draw("hist same ");
			h_DIS_slices[l]->Draw("hist same ");
			h_data_slices[l]->Draw("err same ");

		}
	}
	output_name = output_file_name+"_dxsec_d"+observable+"_"+GetAlternativeObs(observable)+"Slices" ;
  c_slices->SaveAs((output_location+"/TotalXSec/"+output_name+".root").c_str());
  c_slices->SaveAs((output_location+"/TotalXSec/"+output_name+".pdf").c_str());
	delete c_slices;

  // Draw total xsec per sectors
  TCanvas * c_sector = new TCanvas("c_sector","c_sector",200,10,700,500);
  c_sector->cd();
  TPad *pad_sector = new TPad("pad1","",0,0,1,1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->Divide(3,2);

  TPad *pad_sector_0 = (TPad*)pad_sector->cd(1);
  pad_sector_0 -> cd();
  pad_sector_0 -> SetBottomMargin(0.15);
  pad_sector_0 -> SetLeftMargin(0.15);
  hist_true_0 -> GetYaxis()->SetTitleOffset(1.2);
  hist_data_0 -> SetMarkerSize(0.7);
  hist_true_0 -> Draw("hist");
  hist_data_0 -> Draw(" err same ");
  //hist_data_uncorr_0 -> Draw(" err same ");

  TPad *pad_sector_1 = (TPad*)pad_sector->cd(2);
  pad_sector_1 -> cd();
  pad_sector_1 -> SetBottomMargin(0.15);
  pad_sector_1 -> SetLeftMargin(0.15);
  hist_true_1 -> GetYaxis()->SetTitleOffset(1.2);
  hist_data_1 -> SetMarkerSize(0.7);
  hist_true_1 -> Draw("hist");
  hist_data_1 -> Draw(" err same ");
  //hist_data_uncorr_1 -> Draw(" err same ");

  TPad *pad_sector_2 = (TPad*)pad_sector->cd(3);
  pad_sector_2 -> cd();
  pad_sector_2 -> SetBottomMargin(0.15);
  pad_sector_2 -> SetLeftMargin(0.15);
  hist_true_2 -> GetYaxis()->SetTitleOffset(1.2);
  hist_data_2 -> SetMarkerSize(0.7);
  hist_true_2 -> Draw("hist");
  hist_data_2 -> Draw(" err same ");
  //hist_data_uncorr_2 -> Draw(" err same ");

  TPad *pad_sector_3 = (TPad*)pad_sector->cd(4);
  pad_sector_3 -> cd();
  pad_sector_3 -> SetBottomMargin(0.15);
  pad_sector_3 -> SetLeftMargin(0.15);
  hist_true_3 -> GetYaxis()->SetTitleOffset(1.2);
  hist_data_3 -> SetMarkerSize(0.7);
  hist_true_3 -> Draw("hist");
  hist_data_3 -> Draw(" err same ");
  //hist_data_uncorr_3 -> Draw(" err same ");

  TPad *pad_sector_4 = (TPad*)pad_sector->cd(5);
  pad_sector_4 -> cd();
  pad_sector_4 -> SetBottomMargin(0.15);
  pad_sector_4 -> SetLeftMargin(0.15);
  hist_true_4 -> GetYaxis()->SetTitleOffset(1.2);
  hist_data_4 -> SetMarkerSize(0.7);
  hist_true_4 -> Draw("hist");
  hist_data_4 -> Draw(" err same ");
  //hist_data_uncorr_4 -> Draw(" err same ");

  TPad *pad_sector_5 = (TPad*)pad_sector->cd(6);
  pad_sector_5 -> cd();
  pad_sector_5 -> SetBottomMargin(0.15);
  pad_sector_5 -> SetLeftMargin(0.15);
  hist_true_5 -> GetYaxis()->SetTitleOffset(1.2);
  hist_data_5 -> SetMarkerSize(0.7);
  hist_true_5 -> Draw("hist");
  hist_data_5 -> Draw(" err same ");
  //hist_data_uncorr_5 -> Draw(" err same ");

  output_name = MC_files_name[0]+"_dxsec_d"+observable+"_persector" ;
	std::filesystem::path xsecpersector_path{(output_location+"XSecPerSector/").c_str()};
	if( ! std::filesystem::exists(xsecpersector_path) ) std::filesystem::create_directory(xsecpersector_path);
  c_sector->SaveAs((output_location+"XSecPerSector/"+output_name+".root").c_str());
  c_sector->SaveAs((output_location+"XSecPerSector/"+output_name+".pdf").c_str());
  delete c_sector ;

  // Store legend in separate file
  TCanvas * c_leg = new TCanvas("c_leg","c_leg");
  c_leg->cd();
  double LegXmin = 0.1, LegYmin = 0.65, YSpread = 0.25;
	TLegend* leg = new TLegend(LegXmin,LegYmin,LegXmin+0.9,LegYmin+YSpread);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.08);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  leg->AddEntry(hist_true,("GENIE "+model[0]).c_str(),"l");
  leg->AddEntry(hist_true_QEL,(model[0]+" EMQEL").c_str(),"l");
  leg->AddEntry(hist_true_RES,(model[0]+" EMRES").c_str(),"l");
  leg->AddEntry(hist_true_SIS,(model[0]+" EMSIS").c_str(),"l");
  leg->AddEntry(hist_true_MEC,(model[0]+" EMMEC").c_str(),"l");
  leg->AddEntry(hist_true_DIS,(model[0]+" EMDIS").c_str(),"l");
  for( unsigned int id = 1 ; id < model.size() ; ++id ){
     leg->AddEntry(hists_true_submodel[id-1],("GENIE "+model[id]).c_str(),"l");
  }
  leg->AddEntry(hist_data, data_name.c_str(), "lp");
  leg->Draw();
  output_name = MC_files_name[0] ;
  c_leg->SaveAs((output_location+output_name+"_legend.root").c_str());
  c_leg->SaveAs((output_location+output_name+"_legend.pdf").c_str());

  delete c_leg;
}

void Plot1DXSec(){

	// Add without the .root
  std::string mc_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/mc_files/";
  std::string data_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/data_files/";
  std::string output_location ;
	std::filesystem::path output_path ;
  std::vector<std::string> mc_files;
  std::vector<std::string> model_names ;

  std::string file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_4461MeV_clas6data";
  std::string title = "e^{12}C 1p1#pi^{-} at 4.416 GeV";
  std::string data_name = "CLAS6 data";

	//std::vector<std::string> observables = {"RecoW","pfl","pfl_theta","pim_mom","pim_theta","pip_mom","pip_theta","proton_mom","proton_theta",
  //                                        "HadAlphaT","HadDeltaPT","HadDeltaPhiT","ECal","RecoQ2","RecoEnergyTransfer","RecoXBJK","HadSystemMass"};
	std::vector<std::string> observables = {"HadAlphaT", "HadDeltaPT", "ECal"};

  // To be defined in loop
	std::vector<double> binning;
	std::string acceptance_file;
	std::string output_name;

  for( unsigned int i = 0 ; i < observables.size(); ++i ){

		output_name = "e4nuanalysis_1p1pimanalysis_Averaged_Q2_01_e_on_1000060120_1161MeV_NoRad" ;
		output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1p1pim_1GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
		mc_files = {"e4nuanalysis_1p1pimanalysis_GEM21_11a_Q2_01_e_on_1000060120_1161MeV_NoRad",
		            "e4nuanalysis_1p1pimanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad",
							  "e4nuanalysis_1p1pimanalysis_G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad"};
		file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_1161MeV_clas6data";
		title = "e^{12}C 1p1#pi^{-} at 1.116 GeV";
		model_names = { "GEM21_11a","G18_10a", "G18_10b" } ;
		acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
		model_names.push_back("No FSI");
		mc_files.push_back("e4nuanalysis_1p1pimanalysis_G18_10a_NoFSI_Q2_01_e_on_1000060120_1161MeV_NoRad");
		Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;

/*
		output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1p1pim_2GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
    mc_files = {"e4nuanalysis_1p1pimanalysis_GEM21_11a_Q2_04_e_on_1000060120_2261MeV_NoRad",
		            "e4nuanalysis_1p1pimanalysis_G18_10a_Q2_04_e_on_1000060120_2261MeV_NoRad",
		            "e4nuanalysis_1p1pimanalysis_G18_10b_Q2_04_e_on_1000060120_2261MeV_NoRad"};
    file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_2261MeV_clas6data";
    title = "e^{12}C 1p1#pi^{-} at 2.216 GeV";
		output_name = "e4nuanalysis_1p1pimanalysis_Averaged_Q2_04_e_on_1000060120_2261MeV_NoRad" ;
    model_names = { "GEM21_11a","G18_10a", "G18_10b" } ;
	  acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
		model_names.push_back("No FSI");
		mc_files.push_back("e4nuanalysis_1p1pimanalysis_G18_10a_NoFSI_Q2_04_e_on_1000060120_2261MeV_NoRad");
    Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;

		output_name = "e4nuanalysis_1p1pimanalysis_Averaged_Q2_08_e_on_1000060120_4461MeV_NoRad" ;
    output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1p1pim_4GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
		mc_files = {"e4nuanalysis_1p1pimanalysis_G18_10a_Q2_08_e_on_1000060120_4461MeV_NoRad","e4nuanalysis_1p1pimanalysis_G18_10b_Q2_08_e_on_1000060120_4461MeV_NoRad"};
    file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_4461MeV_clas6data";
    title = "e^{12}C 1p1#pi^{-} at 4.416 GeV";
		model_names = { "G18_10a", "G18_10b" } ;
		acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
		Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;

    // Pi+ Analysis
		output_name = "e4nuanalysis_1p1pipanalysis_Averaged_Q2_01_e_on_1000060120_1161MeV_NoRad" ;
		output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1p1pip_1GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
		mc_files = {"e4nuanalysis_1p1pipanalysis_GEM21_11a_Q2_01_e_on_1000060120_1161MeV_NoRad",
	              "e4nuanalysis_1p1pipanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad",
							  "e4nuanalysis_1p1pipanalysis_G18_10b_Q2_01_e_on_1000060120_1161MeV_NoRad"};
    file_data = "e4nuanalysis_1p1pipanalysis_e_on_1000060120_1161MeV_clas6data";
    title = "e^{12}C 1p1#pi^{+} at 1.116 GeV";
    model_names = { "GEM21_11a", "G18_10a", "G18_10b" } ;
		acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
		model_names.push_back("No FSI");
		mc_files.push_back("e4nuanalysis_1p1pipanalysis_G18_10a_NoFSI_Q2_01_e_on_1000060120_1161MeV_NoRad");
		Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;

		output_name = "e4nuanalysis_1p1pipanalysis_Averaged_Q2_04_e_on_1000060120_2261MeV_NoRad" ;
		output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1p1pip_2GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
		mc_files = {"e4nuanalysis_1p1pipanalysis_GEM21_11a_Q2_04_e_on_1000060120_2261MeV_NoRad",
		            "e4nuanalysis_1p1pipanalysis_G18_10a_Q2_04_e_on_1000060120_2261MeV_NoRad",
							  "e4nuanalysis_1p1pipanalysis_G18_10b_Q2_04_e_on_1000060120_2261MeV_NoRad"};
    file_data = "e4nuanalysis_1p1pipanalysis_e_on_1000060120_2261MeV_clas6data";
    title = "e^{12}C 1p1#pi^{+} at 2.216 GeV";
    model_names = { "GEM21_11a", "G18_10a", "G18_10b" } ;
		acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
		model_names.push_back("No FSI");
		mc_files.push_back("e4nuanalysis_1p1pipanalysis_G18_10a_NoFSI_Q2_04_e_on_1000060120_2261MeV_NoRad");
		Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;

		output_name = "e4nuanalysis_1p1pipanalysis_Averaged_Q2_08_e_on_1000060120_4461MeV_NoRad" ;
		output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1p1pip_4GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
		mc_files = {"e4nuanalysis_1p1pipanalysis_G18_10a_Q2_08_e_on_1000060120_4461MeV_NoRad",
								"e4nuanalysis_1p1pipanalysis_G18_10b_Q2_08_e_on_1000060120_4461MeV_NoRad"};
		file_data = "e4nuanalysis_1p1pipanalysis_e_on_1000060120_4461MeV_clas6data";
		title = "e^{12}C 1p1#pi^{+} at 4.416 GeV";
		model_names = { "G18_10a", "G18_10b" } ;
		acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
		Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;

		output_name = "e4nuanalysis_1pimanalysis_Averaged_Q2_04_e_on_1000060120_2261MeV_NoRad" ;
		output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1pim_2GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
		mc_files = {"e4nuanalysis_1pimanalysis_GEM21_11a_Q2_04_e_on_1000060120_2261MeV_NoRad",
								"e4nuanalysis_1pimanalysis_G18_10a_Q2_04_e_on_1000060120_2261MeV_NoRad",
							  "e4nuanalysis_1pimanalysis_G18_10b_Q2_04_e_on_1000060120_2261MeV_NoRad"};
		file_data = "e4nuanalysis_1pimanalysis_e_on_1000060120_2261MeV_clas6data";
		title = "e^{12}C 1#pi^{-} at 2.216 GeV";
		model_names = { "GEM21_11a", "G18_10a", "G18_10b" } ;
		acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
		model_names.push_back("No FSI");
		mc_files.push_back("e4nuanalysis_1pimanalysis_G18_10a_NoFSI_Q2_04_e_on_1000060120_2261MeV_NoRad");
		Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;

    output_name = "e4nuanalysis_1pipanalysis_Averaged_Q2_04_e_on_1000060120_2261MeV_NoRad" ;
    output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1pip_2GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
    mc_files = {"e4nuanalysis_1pipanalysis_GEM21_11a_Q2_04_e_on_1000060120_2261MeV_NoRad",
		    				"e4nuanalysis_1pipanalysis_G18_10a_Q2_04_e_on_1000060120_2261MeV_NoRad",
		 	   		  	"e4nuanalysis_1pipanalysis_G18_10b_Q2_04_e_on_1000060120_2261MeV_NoRad"};
    file_data = "e4nuanalysis_1pimanalysis_e_on_1000060120_2261MeV_clas6data";
    title = "e^{12}C 1#pi^{+} at 2.216 GeV";
    model_names = { "GEM21_11a", "G18_10a", "G18_10b" } ;
    acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
    model_names.push_back("No FSI");
    mc_files.push_back("e4nuanalysis_1pipanalysis_G18_10a_NoFSI_Q2_04_e_on_1000060120_2261MeV_NoRad");
    Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;


    output_name = "e4nuanalysis_1pipanalysis_Averaged_Q2_01_e_on_1000060120_1161MeV_NoRad" ;
    output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files_1pim_1GeV/";
		output_path = output_location.c_str();
		if( ! std::filesystem::exists(output_path) ) std::filesystem::create_directory(output_path);
    mc_files = {"e4nuanalysis_1pimanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad"};
    file_data = "e4nuanalysis_1pimanalysis_e_on_1000060120_1161MeV_clas6data";
    title = "e^{12}C 1#pi^{-} at 1.116 GeV";
    model_names = { "G18_10a" } ;
    acceptance_file = compute_acceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
    Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, -1 ) ;
   */
  }
}
