#include "plotting/XSecUtils.h"
#include "plotting/Systematics.h"
#include "TLegend.h"
#include <iomanip>
#include <filesystem>
#include <sstream>
#include <iostream>
#include <string>

using namespace e4nu ;
using namespace e4nu::plotting ;

void plotting::Plot1DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
			  std::string acceptance_file_name, std::string observable,
			  std::string title, std::string data_name, std::vector<std::string> model,
			  std::string input_MC_location, std::string input_data_location, std::string output_location,
			  std::string output_file_name, bool plot_data, std::map<string,double> systematic_map,
			  std::string analysis_id, bool store_root ){

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
  TFile * file_data = nullptr ;
  if( plot_data ) file_data = new TFile((input_data_location+data_file_name+".root").c_str(),"READ");

  TFile * file_acceptance = new TFile((output_location+acceptance_file_name+".root").c_str(),"READ");
  if( !file_data && plot_data ) { std::cout << "ERROR: the "<< input_data_location << data_file_name << ".root does not exist." <<std::endl; return ;}
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
  std::vector<double> addbinning = GetAdditionalBinning( GetAlternativeObs(observable), BeamE, analysis_id ) ;
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

  TTree * tree_data = nullptr ;
  if( plot_data ) tree_data = (TTree*)file_data->Get("CLAS6Tree");

  if( !h_acceptance ) { std::cout << "ERROR: Acceptance is not defined"<<std::endl; return ; }
  if( !tree_true || (!tree_data && plot_data ) ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return ;}

  // Create histogram for total and total xsec per sector
  TH1D * hist_true = (TH1D*) h_acceptance ->Clone();
  hist_true -> SetName( "MC_True") ;
  hist_true -> Reset();

  TH1D * hist_true_0 = (TH1D*) h_acceptance_0 ->Clone();
  hist_true_0 -> SetName( "MC_True_Sector_0") ;
  hist_true_0 -> Reset();
  TH1D * hist_true_1 = (TH1D*) h_acceptance_1 ->Clone();
  hist_true_1 -> SetName( "MC_True_Sector_1") ;
  hist_true_1 -> Reset();
  TH1D * hist_true_2 = (TH1D*) h_acceptance_2 ->Clone();
  hist_true_2 -> SetName( "MC_True_Sector_2") ;
  hist_true_2 -> Reset();
  TH1D * hist_true_3 = (TH1D*) h_acceptance_3 ->Clone();
  hist_true_3 -> SetName( "MC_True_Sector_3") ;
  hist_true_3 -> Reset();
  TH1D * hist_true_4 = (TH1D*) h_acceptance_4 ->Clone();
  hist_true_4 -> SetName( "MC_True_Sector_4") ;
  hist_true_4 -> Reset();
  TH1D * hist_true_5 = (TH1D*) h_acceptance_5 ->Clone();
  hist_true_5 -> SetName( "MC_True_Sector_5") ;
  hist_true_5 -> Reset();

  // Breakdown histograms for total (all sectors only):
  TH1D * hist_true_QEL = (TH1D*) h_acceptance ->Clone();
  hist_true_QEL -> SetName( "MC_True_QEL") ;
  hist_true_QEL -> Reset();
  TH1D * hist_true_RES_Delta = (TH1D*) h_acceptance ->Clone();
  hist_true_RES_Delta -> SetName( "MC_True_RES_Delta") ;
  hist_true_RES_Delta -> Reset();
  TH1D * hist_true_RES = (TH1D*) h_acceptance ->Clone();
  hist_true_RES -> SetName( "MC_True_RES") ;
  hist_true_RES -> Reset();
  TH1D * hist_true_SIS = (TH1D*) h_acceptance ->Clone();
  hist_true_SIS -> SetName( "MC_True_SIS") ;
  hist_true_SIS -> Reset();
  TH1D * hist_true_MEC = (TH1D*) h_acceptance ->Clone();
  hist_true_MEC -> SetName( "MC_True_MEC") ;
  hist_true_MEC -> Reset();
  TH1D * hist_true_DIS = (TH1D*) h_acceptance ->Clone();
  hist_true_DIS -> SetName( "MC_True_DIS") ;
  hist_true_DIS -> Reset();

  // Same per model - only total prediction
  std::vector<TH1D*> hists_true_submodel;
  for( unsigned int id = 1 ; id < MC_files_name.size(); ++id ){
    hists_true_submodel.push_back( (TH1D*) h_acceptance ->Clone() );
    hists_true_submodel[id - 1] -> SetName( ("MC_True_Model_"+std::to_string(id)).c_str()) ;
    hists_true_submodel[id - 1] -> Reset();
    hists_true_submodel[id - 1] -> SetLineWidth(3);
  }

  // Create hist for each slice on true
  std::vector<TH1D*> h_total_slices, h_QEL_slices, h_RES_Delta_slices, h_RES_slices, h_DIS_slices, h_MEC_slices, h_SIS_slices ;
  if ( addbinning.size() > 0 ) {
    for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
      h_total_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_total_slices[k] -> SetName( ("MC_True_Slice"+std::to_string(k)).c_str() ) ;
      h_total_slices[k] -> Reset();
      h_QEL_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_QEL_slices[k] -> SetName( ("MC_True_QEL_Slice"+std::to_string(k)).c_str() ) ;
      h_QEL_slices[k] -> Reset();
      h_RES_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_RES_slices[k] -> SetName( ("MC_True_RES_Slice"+std::to_string(k)).c_str() ) ;
      h_RES_slices[k] -> Reset();
      h_RES_Delta_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_RES_Delta_slices[k] -> SetName( ("MC_True_RES_Delta_Slice"+std::to_string(k)).c_str() ) ;
      h_RES_Delta_slices[k] -> Reset();
      h_SIS_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_SIS_slices[k] -> SetName( ("MC_True_SIS_Slice"+std::to_string(k)).c_str() ) ;
      h_SIS_slices[k] -> Reset();
      h_MEC_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_MEC_slices[k] -> SetName( ("MC_True_MEC_Slice"+std::to_string(k)).c_str() ) ;
      h_MEC_slices[k] -> Reset();
      h_DIS_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_DIS_slices[k] -> SetName( ("MC_True_DIS_Slice"+std::to_string(k)).c_str() ) ;
      h_DIS_slices[k] -> Reset();
    }
  }


  // total and per sector
  TH1D * hist_data=nullptr, * hist_data_0=nullptr, * hist_data_1=nullptr, * hist_data_2=nullptr, * hist_data_3=nullptr, * hist_data_4=nullptr, * hist_data_5 =nullptr;
  if( plot_data ) {
    hist_data = (TH1D*) h_acceptance ->Clone();
    hist_data -> SetName( "Data") ;
    hist_data -> Reset();

    hist_data_0 = (TH1D*) h_acceptance_0 ->Clone();
    hist_data_0 -> SetName( "Data_Sector_0") ;
    hist_data_0 -> Reset();
    hist_data_1 = (TH1D*) h_acceptance_1 ->Clone();
    hist_data_1 -> SetName( "Data_Sector_1") ;
    hist_data_1 -> Reset();
    hist_data_2 = (TH1D*) h_acceptance_2 ->Clone();
    hist_data_2 -> Reset();
    hist_data_2 -> SetName( "Data_Sector_2") ;
    hist_data_3 = (TH1D*) h_acceptance_3 ->Clone();
    hist_data_3 -> SetName( "Data_Sector_3") ;
    hist_data_3 -> Reset();
    hist_data_4 = (TH1D*) h_acceptance_4 ->Clone();
    hist_data_4 -> SetName( "Data_Sector_4") ;
    hist_data_4 -> Reset();
    hist_data_5 = (TH1D*) h_acceptance_5 ->Clone();
    hist_data_5 -> SetName( "Data_Sector_5") ;
    hist_data_5 -> Reset();
  }
  // Create hist for each slice on data
  std::vector<TH1D*> h_data_slices  ;
  if ( addbinning.size() > 0 && plot_data ) {
    for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
      h_data_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_data_slices[k] -> SetName( ("Data_Slice_"+std::to_string(k)).c_str() ) ;
      h_data_slices[k] -> Reset();
    }
  }

  std::vector<TTree*> trees = { tree_true } ;
  if( plot_data ) trees.push_back( tree_data );

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

  // If data is plot, the position of its trees and histograms is 1.
  // Otherwise we set it to a big number so it is ignored
  unsigned int id_data = 9999;
  if( plot_data ) id_data = 1 ;

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
  int ElectronSector ;
  bool QEL, RES, DIS, MEC;
  double MCNormalization, DataNormalization ;
  std::vector<double> mc_norm ;
  int resid;

  for ( unsigned int i = 0 ; i < trees.size() ; ++i ){
    NEntries = trees[i] -> GetEntries() ;
    if( !trees[i] ) continue ;
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
    trees[i] -> SetBranchAddress("HadDeltaPTx",&HadDeltaPTx);
    trees[i] -> SetBranchAddress("HadDeltaPTy",&HadDeltaPTy);
    trees[i] -> SetBranchAddress("DeltaPhiT",&DeltaPhiT);
    trees[i] -> SetBranchAddress("HadDeltaPhiT",&HadDeltaPhiT);
    trees[i] -> SetBranchAddress("ElectronSector",&ElectronSector);
    trees[i] -> SetBranchAddress("HadSystemMass", &HadSystemMass);
    trees[i] -> SetBranchAddress("MissingEnergy", &MissingEnergy);
    trees[i] -> SetBranchAddress("MissingAngle", &MissingAngle);
    trees[i] -> SetBranchAddress("MissingMomentum", &MissingMomentum);
    trees[i] -> SetBranchAddress("InferedNucleonMom", &InferedNucleonMom);
    trees[i] -> SetBranchAddress("HadronsAngle", &HadronsAngle);
    trees[i] -> SetBranchAddress("AdlerAngleThetaP", &AdlerAngleThetaP);
    trees[i] -> SetBranchAddress("AdlerAnglePhiP", &AdlerAnglePhiP);
    trees[i] -> SetBranchAddress("AdlerAngleThetaPi", &AdlerAngleThetaPi);
    trees[i] -> SetBranchAddress("AdlerAnglePhiPi", &AdlerAnglePhiPi);
    trees[i] -> SetBranchAddress("Angleqvshad",&Angleqvshad);

    // Only fill true info for the first model:
    if( i == 0 ) trees[i] -> SetBranchAddress("QEL",&QEL);
    if( i == 0 ) trees[i] -> SetBranchAddress("RES",&RES);
    if( i == 0 ) trees[i] -> SetBranchAddress("MEC",&MEC);
    if( i == 0 ) trees[i] -> SetBranchAddress("DIS",&DIS);
    if( i == 0 ) trees[i] -> SetBranchAddress("resid", &resid);

    // Only second tree corresponds to data
    if( i != id_data ) {
      trees[i] -> SetBranchAddress("MCNormalization", &MCNormalization );
    } else {
      trees[i] -> SetBranchAddress("DataNormalization",&DataNormalization );
    }

    for( int j = 0 ; j < NEntries ; ++j ) {
      trees[i]->GetEntry(j) ;
      double content = 0 ;
      double w = TotWeight ;
      if( i != id_data && j == 0 ) mc_norm.push_back(MCNormalization);
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

      unsigned int id_hist = i ;
      // Fill the per Sector  histogram. Only for primary model
      if( hists[size_primary_trees*(ElectronSector+1)+i] ) {
	if( i < size_primary_trees ) {
	  hists[size_primary_trees*(ElectronSector+1)+i] -> Fill( content, w ) ;
	}
      }
      if( i > size_primary_trees - 1 ) id_hist = size_primary_hists + ( i - size_primary_trees );

      if( hists[id_hist] ) {
	hists[id_hist] -> Fill( content, w ) ;
	hists[id_hist] -> SetLineWidth(3);
      }

      if( i == 0 ){
        if( QEL ) hist_true_QEL -> Fill( content, w ) ;
        if( RES ) {
	  if( resid == 0 ) hist_true_RES_Delta -> Fill( content, w ) ;
          else hist_true_RES -> Fill( content, w ) ;
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
	double content_2 = 0 ;
	if ( alt_obs == "ECal" ) content_2 = ECal ;
	else if ( alt_obs == "HadAlphaT" ) content_2 = HadAlphaT ;
	else if ( alt_obs == "HadDeltaPT" ) content_2 = HadDeltaPT ;

	for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
	  if( content_2 > addbinning[l] && content_2 < addbinning[l+1] ){
	    if( i == 0 ) { //MC
	      h_total_slices[l] -> Fill( content, w ) ;
	      // Fill also breakdown for slice
	      if( QEL ) h_QEL_slices[l] -> Fill( content, w ) ;
	      if( RES ) {
		if( resid == 0 ) h_RES_Delta_slices[l] -> Fill( content, w ) ;
		else h_RES_slices[l] -> Fill( content, w ) ;
	      }
	      if( DIS ) {
		if( RecoW < 1.7 ) h_SIS_slices[l] -> Fill( content, w ) ;
		else h_DIS_slices[l] -> Fill( content, w ) ;
	      }
	      if( MEC ) h_MEC_slices[l] -> Fill( content, w ) ;
	    } else if( i == 1 && plot_data ) h_data_slices[l] -> Fill( content, w ) ;
	  }
	}
      }
    }
  }


  // Store uncorrected data
  TH1D * hist_data_eventrate=nullptr, * hist_data_eventrate_0=nullptr, * hist_data_eventrate_1=nullptr, * hist_data_eventrate_2=nullptr, * hist_data_eventrate_3=nullptr, * hist_data_eventrate_4=nullptr, * hist_data_eventrate_5=nullptr ;
  //Store before normalising
  if( plot_data ) {
    hist_data_eventrate = (TH1D*) hist_data ->Clone();
    hist_data_eventrate -> SetName( "Data") ;
    hist_data_eventrate_0 = (TH1D*) hist_data_0 ->Clone();
    hist_data_eventrate_0 -> SetName( "Data_Sector_0") ;
    hist_data_eventrate_1 = (TH1D*) hist_data_1 ->Clone();
    hist_data_eventrate_1 -> SetName( "Data_Sector_1") ;
    hist_data_eventrate_2 = (TH1D*) hist_data_2 ->Clone();
    hist_data_eventrate_2 -> SetName( "Data_Sector_2") ;
    hist_data_eventrate_3 = (TH1D*) hist_data_3 ->Clone();
    hist_data_eventrate_3 -> SetName( "Data_Sector_3") ;
    hist_data_eventrate_4 = (TH1D*) hist_data_4 ->Clone();
    hist_data_eventrate_4 -> SetName( "Data_Sector_4") ;
    hist_data_eventrate_5 = (TH1D*) hist_data_5 ->Clone();
    hist_data_eventrate_5 -> SetName( "Data_Sector_5") ;

    // Normaize by bin width
    NormalizeHist(hist_data_eventrate, 1 );
    NormalizeHist(hist_data_eventrate_0, 1 );
    NormalizeHist(hist_data_eventrate_1, 1 );
    NormalizeHist(hist_data_eventrate_2, 1 );
    NormalizeHist(hist_data_eventrate_3, 1 );
    NormalizeHist(hist_data_eventrate_4, 1 );
    NormalizeHist(hist_data_eventrate_5, 1 );
  }


  // Store uncorrected data
  TH1D * hist_data_correventrate=nullptr, * hist_data_correventrate_0=nullptr, * hist_data_correventrate_1=nullptr, * hist_data_correventrate_2=nullptr, * hist_data_correventrate_3=nullptr, * hist_data_correventrate_4=nullptr, * hist_data_correventrate_5=nullptr ;
  if( plot_data ) {
    // Correct data for detector acceptance :
    CorrectData(hist_data, h_acceptance);
    CorrectData(hist_data_0, h_acceptance_0);
    CorrectData(hist_data_1, h_acceptance_1);
    CorrectData(hist_data_2, h_acceptance_2);
    CorrectData(hist_data_3, h_acceptance_3);
    CorrectData(hist_data_4, h_acceptance_4);
    CorrectData(hist_data_5, h_acceptance_5);

		//Adding Acceptance correction systematics from model dependence
		systematics::AddSystematic( *hist_data, *h_acceptance );
		systematics::AddSystematic( *hist_data_0, *h_acceptance );
		systematics::AddSystematic( *hist_data_1, *h_acceptance );
		systematics::AddSystematic( *hist_data_2, *h_acceptance );
		systematics::AddSystematic( *hist_data_3, *h_acceptance );
		systematics::AddSystematic( *hist_data_4, *h_acceptance );
		systematics::AddSystematic( *hist_data_5, *h_acceptance );

    hist_data_correventrate = (TH1D*) hist_data ->Clone();
    hist_data_correventrate -> SetName( "Corrected_Event_Rate_Data") ;
    hist_data_correventrate_0 = (TH1D*) hist_data_0 ->Clone();
    hist_data_correventrate_0 -> SetName( "Corrected_Event_Rate_Data_Sector_0") ;
    hist_data_correventrate_1 = (TH1D*) hist_data_1 ->Clone();
    hist_data_correventrate_1 -> SetName( "Corrected_Event_Rate_Data_Sector_1") ;
    hist_data_correventrate_2 = (TH1D*) hist_data_2 ->Clone();
    hist_data_correventrate_2 -> SetName( "Corrected_Event_Rate_Data_Sector_2") ;
    hist_data_correventrate_3 = (TH1D*) hist_data_3 ->Clone();
    hist_data_correventrate_3 -> SetName( "Corrected_Event_Rate_Data_Sector_3") ;
    hist_data_correventrate_4 = (TH1D*) hist_data_4 ->Clone();
    hist_data_correventrate_4 -> SetName( "Corrected_Event_Rate_Data_Sector_4") ;
    hist_data_correventrate_5 = (TH1D*) hist_data_5 ->Clone();
    hist_data_correventrate_5 -> SetName( "Corrected_Event_Rate_Data_Sector_5") ;

    // Normaize by bin width
    NormalizeHist(hist_data_correventrate, 1 );
    NormalizeHist(hist_data_correventrate_0, 1 );
    NormalizeHist(hist_data_correventrate_1, 1 );
    NormalizeHist(hist_data_correventrate_2, 1 );
    NormalizeHist(hist_data_correventrate_3, 1 );
    NormalizeHist(hist_data_correventrate_4, 1 );
    NormalizeHist(hist_data_correventrate_5, 1 );

    // Normalize data from slices
    if( addbinning.size() != 0 ) {
      for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
	NormalizeHist(h_data_slices[l], DataNormalization );
	CorrectData(h_data_slices[l], h_acc_slices[l] );

	for( auto it = systematic_map.begin() ; it != systematic_map.end() ; ++it ) {
	  systematics::AddSystematic(*h_data_slices[l], it->second, it->first) ;
	}
      }
    }
  }

  NormalizeHist(hist_true, mc_norm[0]);
  NormalizeHist(hist_true_0, mc_norm[0]);
  NormalizeHist(hist_true_1, mc_norm[0]);
  NormalizeHist(hist_true_2, mc_norm[0]);
  NormalizeHist(hist_true_3, mc_norm[0]);
  NormalizeHist(hist_true_4, mc_norm[0]);
  NormalizeHist(hist_true_5, mc_norm[0]);
  NormalizeHist(hist_true_QEL, mc_norm[0]);
  NormalizeHist(hist_true_RES, mc_norm[0]);
  NormalizeHist(hist_true_RES_Delta, mc_norm[0]);
  NormalizeHist(hist_true_SIS, mc_norm[0]);
  NormalizeHist(hist_true_MEC, mc_norm[0]);
  NormalizeHist(hist_true_DIS, mc_norm[0]);

  // add systematics
  if( plot_data ) {
		NormalizeHist(hist_data, DataNormalization );
    NormalizeHist(hist_data_0, DataNormalization );
    NormalizeHist(hist_data_1, DataNormalization );
    NormalizeHist(hist_data_2, DataNormalization );
    NormalizeHist(hist_data_3, DataNormalization );
    NormalizeHist(hist_data_4, DataNormalization );
    NormalizeHist(hist_data_5, DataNormalization );

		// Add sector variation ERROR
		systematics::SectorVariationError( *hist_data, {hist_data_0,hist_data_1,hist_data_2,hist_data_3,hist_data_4,hist_data_5}) ;

    //adding systematics from systematic map. Relative systematic added to all bins
    for( auto it = systematic_map.begin() ; it != systematic_map.end() ; ++it ) {
      std::cout << " Adding " << it->second*100 << " % systematic on " << it->first << std::endl;
      systematics::AddSystematic(*hist_data, it->second, it->first) ;
      systematics::AddSystematic(*hist_data_0, it->second, it->first) ;
      systematics::AddSystematic(*hist_data_1, it->second, it->first) ;
      systematics::AddSystematic(*hist_data_2, it->second, it->first) ;
      systematics::AddSystematic(*hist_data_3, it->second, it->first) ;
      systematics::AddSystematic(*hist_data_4, it->second, it->first) ;
      systematics::AddSystematic(*hist_data_5, it->second, it->first) ;
    }
  }

  for( unsigned int id = 0 ; id < hists_true_submodel.size() ; ++id ){
    NormalizeHist( hists_true_submodel[id], mc_norm[id+1]);
    StandardFormat( hists_true_submodel[id], title, kBlack, 2+id, observable ) ;
  }

  // Normalize true from slices
  if( addbinning.size() != 0 ) {
    for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
      NormalizeHist(h_total_slices[l], mc_norm[0] );
      NormalizeHist(h_QEL_slices[l], mc_norm[0] );
      NormalizeHist(h_RES_Delta_slices[l], mc_norm[0] );
      NormalizeHist(h_RES_slices[l], mc_norm[0] );
      NormalizeHist(h_SIS_slices[l], mc_norm[0] );
      NormalizeHist(h_DIS_slices[l], mc_norm[0] );
      NormalizeHist(h_MEC_slices[l], mc_norm[0] );
    }
  }

  // Store histograms for plotting
  std::vector<TH1D*> mc_hists = {hist_true};
  for( unsigned int id = 0 ; id < hists_true_submodel.size() ; ++id ){
    mc_hists.push_back(hists_true_submodel[id]);
  }

  std::vector<TH1D*> breakdown = {hist_true_QEL,hist_true_RES_Delta,hist_true_RES,hist_true_SIS,hist_true_MEC,hist_true_DIS};
  std::vector<TH1D*> mc_per_sector = {hist_true_0,hist_true_1,hist_true_2,hist_true_3,hist_true_4,hist_true_5};
  std::vector<TH1D*> data_per_sector = {hist_data_0,hist_data_1,hist_data_2,hist_data_3,hist_data_4,hist_data_5};
  std::vector<std::vector<TH1D*>> all_slices = {h_total_slices,h_QEL_slices,h_RES_Delta_slices,h_RES_slices,h_SIS_slices,h_MEC_slices,h_DIS_slices,h_data_slices};

  // Plot Total, XSector, Legend
  plotting::PlotEventRate( hist_data_eventrate, observable, title, data_name, input_data_location, output_location,
			   output_file_name+"_eventrate", analysis_id, store_root ) ;

  plotting::PlotEventRate( hist_data_correventrate, observable, title, data_name, input_data_location, output_location,
			   output_file_name+"_correventrate", analysis_id, store_root ) ;

  plotting::PlotXsecDataTotal( hist_data, observable, title, data_name, input_data_location, output_location, output_file_name,
			       systematic_map, analysis_id, store_root);

  plotting::PlotTotalXSec( mc_hists, breakdown, hist_data, observable, title, data_name, model, input_MC_location,
			   input_data_location, output_location, output_file_name, systematic_map, analysis_id, store_root);

  plotting::PlotPerSector( mc_per_sector, data_per_sector, observable, title, data_name, model,input_MC_location,
			   input_data_location, output_location, output_file_name, systematic_map, analysis_id, store_root ) ;

  plotting::PlotLegend( mc_hists, breakdown, hist_data, observable, data_name, model, output_location, output_file_name, store_root );

  //  plotting::PlotSlices(  all_slices, addbinning, observable, title, data_name, model, input_MC_location, input_data_location,
  //			 output_location, output_file_name,  analysis_id, store_root ) ;
}

void plotting::PlotXsecDataTotal( TH1D * data, std::string observable, std::string title, std::string data_name,
				  std::string input_data_location, std::string output_location,
				  std::string output_file_name, std::map<string,double> systematic_map,
				  std::string analysis_id, bool store_root ) {
  TCanvas * c1 = new TCanvas("c1","c1",200,10,700,500);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Format plots
  if( data ) {
    StandardFormat( data, title, kBlack, 8, observable ) ;
    data -> SetLineStyle(1);
  }

  if( data ) data -> Draw(" err ");

  std::string output_name = output_file_name+"_onlydata_dxsec_d"+observable ;

  std::filesystem::path totalxsec_path{(output_location+"/TotalXSec/").c_str()};
  if( ! std::filesystem::exists(totalxsec_path) ) std::filesystem::create_directory(totalxsec_path);

  if( store_root ) {
		TFile root_file((output_location+"/TotalXSec/"+output_name+".root").c_str(),"recreate");
		data->Write();
		c1->Write();
	}
  c1->SaveAs((output_location+"/TotalXSec/"+output_name+".pdf").c_str());
  delete c1 ;

}

void plotting::PlotTotalXSec( std::vector<TH1D*> mc_hists, std::vector<TH1D*> breakdown, TH1D * data, std::string observable,
			      std::string title, std::string data_name, std::vector<std::string> model,
			      std::string input_MC_location, std::string input_data_location, std::string output_location,
			      std::string output_file_name, std::map<string,double> systematic_map,
			      std::string analysis_id, bool store_root ) {
  TCanvas * c1 = new TCanvas("c1","c1",200,10,700,500);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Find absolute y max
  std::vector<TH1D*> temp_check = {mc_hists[0]};
  if(data) temp_check.push_back(data);

  for( unsigned int id = 1 ; id < mc_hists.size() ; ++id ){
    temp_check.push_back(mc_hists[id]);
  }
  double y_max_total = GetMaximum( temp_check );
  // Format plots
  if( data ) {
    StandardFormat( data, title, kBlack, 8, observable, y_max_total ) ;
    data -> SetLineStyle(1);
  }

  StandardFormat( mc_hists[0], title, kBlack, 1, observable, y_max_total ) ;
  if( breakdown.size() == 6 ) {
    StandardFormat( breakdown[0], title, kBlue-3, 1, observable ) ;
    StandardFormat( breakdown[1], title, kRed-4, 1, observable ) ;
    StandardFormat( breakdown[2], title, kGreen+2, 1, observable ) ;
    StandardFormat( breakdown[3], title, kOrange, 1, observable ) ;
    StandardFormat( breakdown[4], title, kMagenta-3, 1, observable ) ;
    StandardFormat( breakdown[5], title, kCyan+1, 1, observable ) ;
  }

  // Draw total xsec (all sectors):
  mc_hists[0] -> Draw("hist");
  for( unsigned int i = 1; i < mc_hists.size() ; ++i ) mc_hists[i] -> Draw("hist same");
  for( unsigned int i = 0; i < breakdown.size() ; ++i ) breakdown[i] -> Draw("hist same");

  if( data ) data -> Draw(" err same ");

  if( data && observable=="ECal" && plotting::PlotZoomIn(analysis_id) == true ){
    // Add a sub-pad1
    TPad * sub_pad = new TPad("subpad","",0.2,0.2,0.85,0.85);
    sub_pad->SetFillStyle(4000);
    sub_pad->Draw();
    sub_pad->cd();
    sub_pad -> SetBottomMargin(0.15);
    sub_pad -> SetLeftMargin(0.15);

    TH1D* tmp_hist_true = (TH1D*)mc_hists[0]->Clone();
    TH1D* tmp_hist_data ;
    if( data ) tmp_hist_data = (TH1D*)data->Clone();
    tmp_hist_true->SetTitle("");

    //tmp_hist_true->GetXaxis()->SetRangeUser(0,BeamE*(1-0.1));
    if( data ) tmp_hist_true->GetYaxis()->SetRangeUser(0,tmp_hist_data->GetBinContent(tmp_hist_data->GetMaximumBin())*(1+0.25));

    tmp_hist_true -> Draw("hist");
    for( unsigned int i = 1; i < mc_hists.size() ; ++i ) mc_hists[i] -> Draw("hist same");
    for( unsigned int i = 0; i < breakdown.size() ; ++i ) breakdown[i] -> Draw("hist same");

    if( data ) data -> Draw(" err same ");
  }
  std::string output_name = output_file_name+"_dxsec_d"+observable ;

  std::filesystem::path totalxsec_path{(output_location+"/TotalXSec/").c_str()};
  if( ! std::filesystem::exists(totalxsec_path) ) std::filesystem::create_directory(totalxsec_path);

  // Print out integral for debugging
  double data_integral ;
  if( data ) data_integral= data->Integral("width") ;
  double mc_integral = mc_hists[0]->Integral("width") ;
  if( data ) std::cout << " Total integrated cross section (data) " << data_integral << std::endl;
  std::cout << " Total integrated cross section (mc) " << mc_integral << std::endl;

  if( store_root ) {
		TFile root_file((output_location+"/TotalXSec/"+output_name+".root").c_str(),"recreate");
		data->Write();
		for( unsigned int i = 0; i < mc_hists.size() ; ++i ) mc_hists[i] -> Write();
		for( unsigned int i = 0; i < breakdown.size() ; ++i ) breakdown[i] -> Write();
		c1->Write();
	}
  c1->SaveAs((output_location+"/TotalXSec/"+output_name+".pdf").c_str());
  delete c1 ;

}

void plotting::PlotEventRate( TH1D * data, std::string observable, std::string title, std::string data_name, std::string input_data_location,
			      std::string output_location, std::string output_file_name, std::string analysis_id, bool store_root ) {
  TCanvas * c1 = new TCanvas("c1","c1",200,10,700,500);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  // Format plots
  if( data ) {
    StandardFormat( data, title, kBlack, 8, observable, 0, "Counts/Bin Width" ) ;
    data -> SetLineStyle(1);
  }

  if( data ) {
    data -> Draw(" err ");
  }

  double data_integral ;
  if( data ) data_integral= data->Integral("width") ;
  if( data ) std::cout << " Total integrated event rate (data) " << data_integral << std::endl;

  std::string output_name = output_file_name+"_Nevents_"+observable ;

  std::filesystem::path totalxsec_path{(output_location+"/EventRate/").c_str()};
  if( ! std::filesystem::exists(totalxsec_path) ) std::filesystem::create_directory(totalxsec_path);

  if( store_root ) {
		TFile root_file((output_location+"/EventRate/"+output_name+".root").c_str(),"recreate");
		data->Write();
		c1->Write();
	}
  c1->SaveAs((output_location+"/EventRate/"+output_name+".pdf").c_str());
  delete c1 ;

}

void plotting::PlotLegend( std::vector<TH1D*> mc_hists, std::vector<TH1D*> breakdown, TH1D * data, std::string observable,
			   std::string data_name, std::vector<std::string> model,std::string output_location,
			   std::string output_file_name, bool store_root ) {

  // Store legend in separate file
  TCanvas * c_leg = new TCanvas("c_leg","c_leg");
  c_leg->cd();
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  double LegXmin = 0.1, LegYmin = 0.65, YSpread = 0.25;
  TLegend* leg = new TLegend(LegXmin,LegYmin,LegXmin+0.9,LegYmin+YSpread);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.08);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  std::string model_def = "MC";
  if ( model.size() > 0 ) model_def = model[0];
  if( mc_hists[0] ) {
    leg->AddEntry(mc_hists[0],("GENIE "+model_def).c_str(),"l");
    leg->AddEntry(breakdown[0],(model_def+" EMQEL").c_str(),"l");
    leg->AddEntry(breakdown[1],(model_def+" EMRES P33(1232)").c_str(),"l");
    leg->AddEntry(breakdown[2],(model_def+" EMRES Others").c_str(),"l");
    leg->AddEntry(breakdown[3],(model_def+" EMSIS").c_str(),"l");
    leg->AddEntry(breakdown[4],(model_def+" EMMEC").c_str(),"l");
    leg->AddEntry(breakdown[5],(model_def+" EMDIS").c_str(),"l");
  }

  if( mc_hists.size() > 1 ) {
    for( unsigned int id = 1 ; id < mc_hists.size() ; ++id ){
      if( !mc_hists[id] ) continue ;
      std::string model_id = "Model " + std::to_string(id) ;
      if( model.size() == mc_hists.size() ) model_id = model[id];
      leg->AddEntry(mc_hists[id],("GENIE "+model_id).c_str(),"l");
    }
  }

  std::string data_def = "CLAS6 data";
  if( data_name!="") data_def = data_name ;
  if(data)  leg->AddEntry(data, data_def.c_str(), "lp");
  leg->Draw();
  std::string output_name = output_file_name ;
  if( store_root ) c_leg->SaveAs((output_location+"/"+output_name+"_legend.root").c_str());
  c_leg->SaveAs((output_location+"/"+output_name+"_legend.pdf").c_str());

  delete c_leg;
}

void plotting::PlotPerSector( std::vector<TH1D*> mc_per_sector,std::vector<TH1D*> data_per_sector, std::string observable,
			      std::string title, std::string data_name, std::vector<std::string> model,
			      std::string input_MC_location, std::string input_data_location, std::string output_location,
			      std::string output_file_name, std::map<string,double> systematic_map,
			      std::string analysis_id, bool store_root ) {

  // Format plots
  if( data_per_sector.size()!=0 ) {
    if( data_per_sector[0] ) StandardFormat( data_per_sector[0], title+" Sector  0", kOrange+1, 8, observable ) ;
    if( data_per_sector[1] ) StandardFormat( data_per_sector[1], title+" Sector  1", kPink+4, 8, observable ) ;
    if( data_per_sector[2] ) StandardFormat( data_per_sector[2], title+" Sector  2", kViolet+5, 8, observable ) ;
    if( data_per_sector[3] ) StandardFormat( data_per_sector[3], title+" Sector  3", kAzure-5, 8, observable ) ;
    if( data_per_sector[4] ) StandardFormat( data_per_sector[4], title+" Sector  4", kTeal-7, 8, observable ) ;
    if( data_per_sector[5] ) StandardFormat( data_per_sector[5], title+" Sector  5", kGreen-3, 8, observable ) ;
  }

  if( mc_per_sector.size() != 0 ) {
    StandardFormat( mc_per_sector[0], title+" Sector  0", kBlack, 1, observable ) ;
    StandardFormat( mc_per_sector[1], title+" Sector  1", kBlack, 1, observable ) ;
    StandardFormat( mc_per_sector[2], title+" Sector  2", kBlack, 1, observable ) ;
    StandardFormat( mc_per_sector[3], title+" Sector  3", kBlack, 1, observable ) ;
    StandardFormat( mc_per_sector[4], title+" Sector  4", kBlack, 1, observable ) ;
    StandardFormat( mc_per_sector[5], title+" Sector  5", kBlack, 1, observable ) ;
  }

  // Draw total xsec per sectors
  TCanvas * c_sector = new TCanvas("c_sector","c_sector",200,10,700,500);
  c_sector->cd();
  TPad *pad_sector = new TPad("pad1","",0,0,1,1);
  pad_sector->Draw();
  pad_sector->cd();
  pad_sector->SetBottomMargin(0.15);
  pad_sector->SetLeftMargin(0.15);
  pad_sector->Divide(3,2);

  for( unsigned int i = 0 ; i < 6 ; ++i ) {
    TPad *pad_sector_i = (TPad*)pad_sector->cd(i+1);
    pad_sector_i -> cd();
    pad_sector_i -> SetBottomMargin(0.15);
    pad_sector_i -> SetLeftMargin(0.15);
    mc_per_sector[i] -> GetYaxis()->SetTitleOffset(1.2);
    mc_per_sector[i] -> Draw("hist");

    if(data_per_sector.size()!=0 && data_per_sector[i]) {
      data_per_sector[i] -> SetMarkerSize(0.7);
      data_per_sector[i] -> Draw(" err same ");
    }
  }

  std::string output_name = output_file_name+"_dxsec_d"+observable+"_persector" ;
  std::filesystem::path xsecpersector_path{(output_location+"/XSecPerSector/").c_str()};
  if( ! std::filesystem::exists(xsecpersector_path) ) std::filesystem::create_directory(xsecpersector_path);
  if( store_root ) c_sector->SaveAs((output_location+"/XSecPerSector/"+output_name+".root").c_str());
  c_sector->SaveAs((output_location+"/XSecPerSector/"+output_name+".pdf").c_str());


  for( unsigned int i = 0 ; i < 6 ; ++i ) {
    TPad *pad_sector_i = (TPad*)pad_sector->cd(i+1);
    pad_sector_i -> cd();
    pad_sector_i -> SetBottomMargin(0.15);
    pad_sector_i -> SetLeftMargin(0.15);

    if(data_per_sector.size()!=0 && data_per_sector[i]) {
      data_per_sector[i] -> SetMarkerSize(0.7);
      data_per_sector[i] -> Draw(" err ");
    }
  }

  output_name = output_file_name+"_dataonly_dxsec_d"+observable+"_persector" ;
  xsecpersector_path =(output_location+"/XSecPerSector/").c_str();
  if( ! std::filesystem::exists(xsecpersector_path) ) std::filesystem::create_directory(xsecpersector_path);

	if( store_root ) {
		TFile root_file((output_location+"/XSecPerSector/"+output_name+".root").c_str(),"recreate");
		for( unsigned int i = 0 ; i < 6 ; ++i ) {
			mc_per_sector[i] -> Write();
			if(data_per_sector.size()!=0 && data_per_sector[i]) {
				data_per_sector[i] -> Write();
			}
		}
		c_sector->Write();
	}

  c_sector->SaveAs((output_location+"/XSecPerSector/"+output_name+".pdf").c_str());
  delete c_sector ;
}

void plotting::PlotSlices( std::vector<std::vector<TH1D*>> all_slices, std::vector<double> addbinning, std::string observable,
			   std::string title, std::string data_name, std::vector<std::string> model,
			   std::string input_MC_location, std::string input_data_location, std::string output_location,
			   std::string output_file_name, std::map<string,double> systematic_map,
			   std::string analysis_id, bool store_root ) {

  for( unsigned int l = 0 ; l < addbinning.size()-1 ; l++ ){
    double y_max_total = GetMaximum( all_slices[0] );

    // Add Slice information in title
    std::string title_subname = title ;
    std::string alt_obs = GetAlternativeObs(observable);
    if ( l == 0 ) {
      std::ostringstream o1 ;
      o1 << std::fixed<< std::setprecision(1) << addbinning[l+1] ;
      title_subname += " " + plotting::GetObsName(alt_obs) + "<" + o1.str() +" "+plotting::GetUnit(alt_obs) ;
    } else if ( l == addbinning.size()-2 ) {
      std::ostringstream o1 ;
      o1 << std::fixed<< std::setprecision(1) << addbinning[l] ;
      title_subname += " " + plotting::GetObsName(alt_obs) + ">" + o1.str() +" "+plotting::GetUnit(alt_obs) ;
    } else {
      std::ostringstream o1, o2 ;
      o1 << std::fixed<< std::setprecision(1) << addbinning[l] ;
      o2 << std::fixed<< std::setprecision(1) << addbinning[l+1] ;
      title_subname += " " + o1.str() + "<"+ plotting::GetObsName(alt_obs) + "<" + o2.str()+" "+plotting::GetUnit(alt_obs) ;
    }

    if( all_slices[0][l]) StandardFormat( all_slices[0][l], title_subname, kBlack, 1, observable, y_max_total ) ;
    if( all_slices[1][l]) StandardFormat( all_slices[1][l], title_subname, kBlue-3, 1, observable, y_max_total ) ;
    if( all_slices[2][l]) StandardFormat( all_slices[2][l], title_subname, kRed-4, 1, observable, y_max_total ) ;
    if( all_slices[3][l]) StandardFormat( all_slices[3][l], title_subname, kGreen+2, 1, observable, y_max_total ) ;
    if( all_slices[4][l]) StandardFormat( all_slices[4][l], title_subname, kOrange, 1, observable, y_max_total ) ;
    if( all_slices[5][l]) StandardFormat( all_slices[5][l], title_subname, kMagenta-3, 1, observable, y_max_total ) ;
    if( all_slices[6][l]) StandardFormat( all_slices[6][l], title_subname, kCyan+1, 1, observable, y_max_total ) ;
    if( all_slices.size() == 8 && all_slices[7][l]) StandardFormat( all_slices[7][l], title_subname, kBlack, 8, observable, y_max_total ) ;
  }

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
      if( all_slices[0][l]) {
	all_slices[0][l] -> GetYaxis()->SetTitleOffset(1.2);
	all_slices[0][l]->Draw("hist");
      }

      for( unsigned i = 1 ; i < all_slices.size() -1 ; ++i ) {
	if( all_slices[i][l] ) all_slices[i][l] ->Draw("same");
      }
      if( all_slices.size() == 8 && all_slices[8][l] ) all_slices[8][l] -> Draw("err same ");
    }
  }

  std::string output_name = output_file_name+"_dxsec_d"+observable+"_"+GetAlternativeObs(observable)+"_Slices" ;
  if( store_root ) c_slices->SaveAs((output_location+"/TotalXSec/"+output_name+".root").c_str());
  c_slices->SaveAs((output_location+"/TotalXSec/"+output_name+".pdf").c_str());

  delete c_slices;
}
