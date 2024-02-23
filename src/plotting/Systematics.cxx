#include "plotting/Systematics.h"
#include "plotting/PlottingUtils.h"
#include "plotting/XSecUtils.h"
#include "TLegend.h"
#include <iomanip>
#include "TMath.h"
#include <filesystem>
#include <sstream>
#include <iostream>
#include <string>

using namespace e4nu ;
using namespace e4nu::systematics ;

void systematics::ComputeHistSyst( std::vector<std::string> input_files, std::vector<std::string> tags, std::string observable, bool is_data,
				   std::string input_location, std::string output_location, std::string analysis_id ){
  if( input_files.size() ==0 ) return ;

  std::string name_fullpath = output_location+"/"+input_files[0]+"_"+observable+"_systematics.root";
  std::unique_ptr<TFile> myFile( TFile::Open(name_fullpath.c_str(), "RECREATE") );

  std::vector<TFile*> ifiles ;
  for( unsigned int id = 0 ; id < input_files.size(); ++id ){
    ifiles.push_back(new TFile((input_location+"/"+input_files[id]+".root").c_str(),"ROOT"));
    if( !ifiles[id] ) { std::cout << "ERROR: the "<< input_location<<input_files[id]<<".root does not exist." << std::endl; return ;}
  }

  // Get Tree for main model
  std::string treename = "MCCLAS6Tree";
  if( is_data ) treename = "CLAS6Tree";

  std::vector<TH1D*> hists;
  for( unsigned int i = 0 ; i < input_files.size(); ++i ) {
    double BeamE ;
    TTree * tree = (TTree*)ifiles[i]->Get(treename.c_str());
    if( !tree ) return ;
    tree ->SetBranchAddress("BeamE",&BeamE);
    tree->GetEntry(0);
    std::vector<double> binning = plotting::GetBinning(observable,BeamE,analysis_id);
    if( binning.size()==0 ) return ;

    // Create histogram for total and total xsec per sector
    TH1D * hist = new TH1D( tags[i].c_str(), tags[i].c_str(), binning.size()-1, &binning[0] ) ;
    if( !hist ) return ;

    // OBSERVABLE DEFINITION:
    double TotWeight ;
    double ECal,Recoq3,RecoW;
    double pfl,pfl_theta,pfl_phi;
    double proton_mom,proton_phi,proton_theta;
    double pim_mom,pim_theta,pim_phi;
    double pip_mom,pip_theta,pip_phi;
    double HadAlphaT, HadDeltaPT, HadDeltaPTx, HadDeltaPTy, HadDeltaPhiT ;
    double AlphaT, DeltaPT, DeltaPhiT ;
    double RecoXBJK, RecoEnergyTransfer, RecoQ2, HadSystemMass, RecoQELEnu ;
    double MissingEnergy, MissingAngle, MissingMomentum ;
    double InferedNucleonMom ;
    double HadronsAngle,Angleqvshad ;
    double AdlerAngleThetaP, AdlerAnglePhiP, AdlerAngleThetaPi, AdlerAnglePhiPi ;
    long NEntries ;
    bool IsBkg ;
    double MCNormalization, DataNormalization ;
    std::vector<double> mc_norm ;

    if( !tree ) continue ;
    NEntries = tree -> GetEntries() ;

    tree -> SetBranchAddress("TotWeight",&TotWeight);
    tree -> SetBranchAddress("IsBkg",&IsBkg);
    tree -> SetBranchAddress("ECal",&ECal);
    tree -> SetBranchAddress("pfl_theta",&pfl_theta);
    tree -> SetBranchAddress("pfl_phi",&pfl_phi);
    tree -> SetBranchAddress("pfl",&pfl);
    tree -> SetBranchAddress("proton_mom",&proton_mom);
    tree -> SetBranchAddress("proton_theta",&proton_theta);
    tree -> SetBranchAddress("proton_phi",&proton_phi);
    tree -> SetBranchAddress("pim_mom",&pim_mom);
    tree -> SetBranchAddress("pim_theta",&pim_theta);
    tree -> SetBranchAddress("pim_phi",&pim_phi);
    tree -> SetBranchAddress("pip_mom",&pip_mom);
    tree -> SetBranchAddress("pip_theta",&pip_theta);
    tree -> SetBranchAddress("pip_phi",&pip_phi);
    tree -> SetBranchAddress("RecoW",&RecoW);
    tree -> SetBranchAddress("Recoq3",&Recoq3);
    tree -> SetBranchAddress("RecoQELEnu",&RecoQELEnu);
    tree -> SetBranchAddress("RecoXBJK",&RecoXBJK);
    tree -> SetBranchAddress("RecoQ2",&RecoQ2);
    tree -> SetBranchAddress("RecoEnergyTransfer",&RecoEnergyTransfer);
    tree -> SetBranchAddress("AlphaT",&AlphaT);
    tree -> SetBranchAddress("HadAlphaT",&HadAlphaT);
    tree -> SetBranchAddress("DeltaPT",&DeltaPT);
    tree -> SetBranchAddress("HadDeltaPT",&HadDeltaPT);
    tree -> SetBranchAddress("HadDeltaPTx",&HadDeltaPTx);
    tree -> SetBranchAddress("HadDeltaPTy",&HadDeltaPTy);
    tree -> SetBranchAddress("DeltaPhiT",&DeltaPhiT);
    tree -> SetBranchAddress("HadDeltaPhiT",&HadDeltaPhiT);
    tree -> SetBranchAddress("HadSystemMass", &HadSystemMass);
    tree -> SetBranchAddress("MissingEnergy", &MissingEnergy);
    tree -> SetBranchAddress("MissingAngle", &MissingAngle);
    tree -> SetBranchAddress("MissingMomentum", &MissingMomentum);
    tree -> SetBranchAddress("InferedNucleonMom", &InferedNucleonMom);
    tree -> SetBranchAddress("HadronsAngle", &HadronsAngle);
    tree -> SetBranchAddress("AdlerAngleThetaP", &AdlerAngleThetaP);
    tree -> SetBranchAddress("AdlerAnglePhiP", &AdlerAnglePhiP);
    tree -> SetBranchAddress("AdlerAngleThetaPi", &AdlerAngleThetaPi);
    tree -> SetBranchAddress("AdlerAnglePhiPi", &AdlerAnglePhiPi);
    tree -> SetBranchAddress("Angleqvshad",&Angleqvshad);

    // Only second tree corresponds to data
    if( !is_data ) {
      tree -> SetBranchAddress("MCNormalization", &MCNormalization );
    } else {
      tree -> SetBranchAddress("DataNormalization",&DataNormalization );
    }

    double norm =0;
    for( int j = 0 ; j < NEntries ; ++j ) {
      tree->GetEntry(j) ;
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
      else if ( observable == "HadDeltaPTx") content = HadDeltaPTx ;
      else if ( observable == "HadDeltaPTy") content = HadDeltaPTy ;
      else if ( observable == "DeltaPhiT") content = DeltaPhiT ;
      else if ( observable == "HadDeltaPhiT") content = HadDeltaPhiT ;
      else if ( observable == "HadSystemMass") content = HadSystemMass ;
      else if ( observable == "MissingEnergy") content = MissingEnergy ;
      else if ( observable == "MissingEnergy") content = MissingEnergy ;
      else if ( observable == "MissingAngle") content = MissingAngle ;
      else if ( observable == "MissingMomentum") content = MissingMomentum ;
      else if ( observable == "InferedNucleonMom") content = InferedNucleonMom ;
      else if ( observable == "HadronsAngle") content = HadronsAngle ;
      else if ( observable == "AdlerAngleThetaP") content = AdlerAngleThetaP ;
      else if ( observable == "AdlerAnglePhiP") content = AdlerAnglePhiP ;
      else if ( observable == "AdlerAngleThetaPi") content = AdlerAngleThetaPi ;
      else if ( observable == "AdlerAnglePhiPi") content = AdlerAnglePhiPi ;
      else if ( observable == "Angleqvshad") content = Angleqvshad ;

      norm = DataNormalization ;
      if( !is_data ) norm = MCNormalization ;

      hist -> Fill( content, w ) ;
      hist -> SetLineWidth(3);
    }// entries loop

    plotting::NormalizeHist(hist,norm);
    myFile->WriteObject(hist,tags[i].c_str());
    hists.push_back(hist);

  } // input files

  if( hists.size() > 1 ) {
    for( unsigned int i = 1 ; i < hists.size() ; ++i ) {
      TH1D * diff = (TH1D*) hists[0]->Clone();
      diff -> Add( hists[i], -1);
      diff -> SetName( ("diff_def_vs_"+std::to_string(i)).c_str()) ;
      myFile->WriteObject(diff, ("diff_def_vs_"+std::to_string(i)).c_str() );
    }
  }

}


void systematics::AddSystematic( TH1D & hist, const double rel_error, const std::string name ) {
  double NBins = hist.GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double error = hist.GetBinError(i);
    double content = hist.GetBinContent(i);
    double newerror = TMath::Sqrt( TMath::Power(error,2.) + TMath::Power(rel_error*content,2.));
    hist.SetBinError(i,newerror);
  }
}

void systematics::AddSystematic( TH1D & hist, const TH1D & hist_w_error ) {
  double NBins = hist.GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double stat_error = hist.GetBinError(i);
		double syst_error = hist_w_error.GetBinError(i);
    double newerror = TMath::Sqrt( TMath::Power(stat_error,2.) + TMath::Power(syst_error,2.));
    hist.SetBinError(i,newerror);
  }
}

void systematics::SectorVariationError( TH1D & hist, const std::vector<TH1D*> h_per_sector ) {

	double diff = 0 ;
	double raw_sector_variance = 0;

	for( unsigned j = 0 ; j < hist.GetNbinsX() ; ++j ){
		double var_j = h_per_sector[0]->GetBinError(j);
		double content_j = h_per_sector[0]->GetBinContent(j);
		double var_total_stat = var_j;
		double mean_content = 0 ;
		for( unsigned int i = 1 ; i < h_per_sector.size() ; ++i ){
			var_total_stat += h_per_sector[i]->GetBinError(j);
			content_j = h_per_sector[i]->GetBinContent(j);
			mean_content += content_j/var_j;
		}
		mean_content /= (1./var_total_stat);
		raw_sector_variance = pow( content_j - mean_content, 2 )/(h_per_sector.size()-1);
		double corrected_sector_variance = raw_sector_variance - pow(var_total_stat,2) ; // need to remove double counting of stat error
		double hist_error = TMath::Sqrt( pow(hist.GetBinError(j),2) + pow(corrected_sector_variance,2));
		hist.SetBinError(j,hist_error);
	}

}
