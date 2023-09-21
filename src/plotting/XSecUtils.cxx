#include "plotting/XSecUtils.h"
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
			  std::string output_file_name, bool plot_data ){

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

  TTree * tree_data ;
  if( plot_data ) tree_data = (TTree*)file_data->Get("CLAS6Tree");

  if( !h_acceptance ) { std::cout << "ERROR: Acceptance is not defined"<<std::endl; return ; }
  if( !tree_true || (!tree_data && plot_data ) ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return ;}

  // Create histogram for total and total xsec per sector
  TH1D * hist_true = (TH1D*) h_acceptance ->Clone();
  hist_true -> SetName( "MC True") ;
  hist_true -> Reset();

  TH1D * hist_true_0 = (TH1D*) h_acceptance_0 ->Clone();
  hist_true_0 -> SetName( "MC True Sector  0") ;
  hist_true_0 -> Reset();
  TH1D * hist_true_1 = (TH1D*) h_acceptance_1 ->Clone();
  hist_true_1 -> SetName( "MC True Sector  1") ;
  hist_true_1 -> Reset();
  TH1D * hist_true_2 = (TH1D*) h_acceptance_2 ->Clone();
  hist_true_2 -> SetName( "MC True Sector  2") ;
  hist_true_2 -> Reset();
  TH1D * hist_true_3 = (TH1D*) h_acceptance_3 ->Clone();
  hist_true_3 -> SetName( "MC True Sector  3") ;
  hist_true_3 -> Reset();
  TH1D * hist_true_4 = (TH1D*) h_acceptance_4 ->Clone();
  hist_true_4 -> SetName( "MC True Sector  4") ;
  hist_true_4 -> Reset();
  TH1D * hist_true_5 = (TH1D*) h_acceptance_5 ->Clone();
  hist_true_5 -> SetName( "MC True Sector  5") ;
  hist_true_5 -> Reset();

  // Breakdown histograms for total (all sectors only):
  TH1D * hist_true_QEL = (TH1D*) h_acceptance ->Clone();
  hist_true_QEL -> SetName( "MC True QEL") ;
  hist_true_QEL -> Reset();
  TH1D * hist_true_RES_Delta = (TH1D*) h_acceptance ->Clone();
  hist_true_RES_Delta -> SetName( "MC True RES Delta") ;
  hist_true_RES_Delta -> Reset();
  TH1D * hist_true_RES = (TH1D*) h_acceptance ->Clone();
  hist_true_RES -> SetName( "MC True RES") ;
  hist_true_RES -> Reset();
  TH1D * hist_true_SIS = (TH1D*) h_acceptance ->Clone();
  hist_true_SIS -> SetName( "MC True SIS") ;
  hist_true_SIS -> Reset();
  TH1D * hist_true_MEC = (TH1D*) h_acceptance ->Clone();
  hist_true_MEC -> SetName( "MC True MEC") ;
  hist_true_MEC -> Reset();
  TH1D * hist_true_DIS = (TH1D*) h_acceptance ->Clone();
  hist_true_DIS -> SetName( "MC True DIS") ;
  hist_true_DIS -> Reset();

  // Same per model - only total prediction
  std::vector<TH1D*> hists_true_submodel;
  for( unsigned int id = 1 ; id < MC_files_name.size(); ++id ){
    hists_true_submodel.push_back( (TH1D*) h_acceptance ->Clone() );
    hists_true_submodel[id - 1] -> SetName( ("MC True Model "+std::to_string(id)).c_str()) ;
    hists_true_submodel[id - 1] -> Reset();
    hists_true_submodel[id - 1] -> SetLineWidth(3);
  }

  // Create hist for each slice on true
  std::vector<TH1D*> h_total_slices, h_QEL_slices, h_RES_Delta_slices, h_RES_slices, h_DIS_slices, h_MEC_slices, h_SIS_slices ;
  if ( addbinning.size() > 0 ) {
    for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
      h_total_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_total_slices[k] -> SetName( ("MC True Slice "+std::to_string(k)).c_str() ) ;
      h_total_slices[k] -> Reset();
      h_QEL_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_QEL_slices[k] -> SetName( ("MC True QEL Slice "+std::to_string(k)).c_str() ) ;
      h_QEL_slices[k] -> Reset();
      h_RES_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_RES_slices[k] -> SetName( ("MC True RES Slice "+std::to_string(k)).c_str() ) ;
      h_RES_slices[k] -> Reset();
      h_RES_Delta_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_RES_Delta_slices[k] -> SetName( ("MC True RES Delta Slice "+std::to_string(k)).c_str() ) ;
      h_RES_Delta_slices[k] -> Reset();
      h_SIS_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_SIS_slices[k] -> SetName( ("MC True SIS Slice "+std::to_string(k)).c_str() ) ;
      h_SIS_slices[k] -> Reset();
      h_MEC_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_MEC_slices[k] -> SetName( ("MC True MEC Slice "+std::to_string(k)).c_str() ) ;
      h_MEC_slices[k] -> Reset();
      h_DIS_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_DIS_slices[k] -> SetName( ("MC True DIS Slice "+std::to_string(k)).c_str() ) ;
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
    hist_data_0 -> SetName( "Data Sector  0") ;
    hist_data_0 -> Reset();
    hist_data_1 = (TH1D*) h_acceptance_1 ->Clone();
    hist_data_1 -> SetName( "Data Sector  1") ;
    hist_data_1 -> Reset();
    hist_data_2 = (TH1D*) h_acceptance_2 ->Clone();
    hist_data_2 -> Reset();
    hist_data_2 -> SetName( "Data Sector  2") ;
    hist_data_3 = (TH1D*) h_acceptance_3 ->Clone();
    hist_data_3 -> SetName( "Data Sector  3") ;
    hist_data_3 -> Reset();
    hist_data_4 = (TH1D*) h_acceptance_4 ->Clone();
    hist_data_4 -> SetName( "Data Sector  4") ;
    hist_data_4 -> Reset();
    hist_data_5 = (TH1D*) h_acceptance_5 ->Clone();
    hist_data_5 -> SetName( "Data Sector  5") ;
    hist_data_5 -> Reset();
  }
  // Create hist for each slice on data
  std::vector<TH1D*> h_data_slices  ;
  if ( addbinning.size() > 0 && plot_data ) {
    for( unsigned int k = 0 ; k < addbinning.size()-1 ; k++ ){
      h_data_slices.push_back( (TH1D*) h_acc_slices[k] ->Clone() ) ;
      h_data_slices[k] -> SetName( ("Data Slice "+std::to_string(k)).c_str() ) ;
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
  double HadronsAngle ; 
  long NEntries ;
  bool IsBkg ;
  int ElectronSector ;
  bool QEL, RES, DIS, MEC;
  double MCNormalization, DataNormalization ;
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

    // Only fill true info for the first model:
    if( i == 0 ) trees[i] -> SetBranchAddress("QEL",&QEL);
    if( i == 0 ) trees[i] -> SetBranchAddress("RES",&RES);
    if( i == 0 ) trees[i] -> SetBranchAddress("MEC",&MEC);
    if( i == 0 ) trees[i] -> SetBranchAddress("DIS",&DIS);
    if( i == 0 ) trees[i] -> SetBranchAddress("resid", &resid);

    // Only second tree corresponds to data
    if( i != 1 ) trees[i] -> SetBranchAddress("MCNormalization", &MCNormalization );
    else {
      if( plot_data ) trees[i] -> SetBranchAddress("DataNormalization",&DataNormalization );
    }

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

      unsigned int id_hist = i ;
      // Fill the per Sector  histogram. Only for primary model
      if( hists[size_primary_trees*(ElectronSector+1)+i] ) { 
	if( i < size_primary_trees ) hists[size_primary_trees*(ElectronSector+1)+i] -> Fill( content, w ) ;
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

  // Normalize data
  if( plot_data ) { 
    NormalizeHist(hist_data, DataNormalization );
    NormalizeHist(hist_data_0, DataNormalization );
    NormalizeHist(hist_data_1, DataNormalization );
    NormalizeHist(hist_data_2, DataNormalization );
    NormalizeHist(hist_data_3, DataNormalization );
    NormalizeHist(hist_data_4, DataNormalization );
    NormalizeHist(hist_data_5, DataNormalization );
  }

  // Store uncorrected data
  TH1D * hist_data_uncorr=nullptr, * hist_data_uncorr_0=nullptr, * hist_data_uncorr_1=nullptr, * hist_data_uncorr_2=nullptr, * hist_data_uncorr_3=nullptr, * hist_data_uncorr_4=nullptr, * hist_data_uncorr_5=nullptr ;
  if( plot_data ) { 
    hist_data_uncorr = (TH1D*) hist_data ->Clone();
    hist_data_uncorr -> SetName( "Uncorrected Data") ;
    hist_data_uncorr_0 = (TH1D*) hist_data_0 ->Clone();
    hist_data_uncorr_0 -> SetName( "Uncorrected Data Sector  0") ;
    hist_data_uncorr_1 = (TH1D*) hist_data_1 ->Clone();
    hist_data_uncorr_1 -> SetName( "Uncorrected Data Sector  1") ;
    hist_data_uncorr_2 = (TH1D*) hist_data_2 ->Clone();
    hist_data_uncorr_2 -> SetName( "Uncorrected Data Sector  2") ;
    hist_data_uncorr_3 = (TH1D*) hist_data_3 ->Clone();
    hist_data_uncorr_3 -> SetName( "Uncorrected Data Sector  3") ;
    hist_data_uncorr_4 = (TH1D*) hist_data_4 ->Clone();
    hist_data_uncorr_4 -> SetName( "Uncorrected Data Sector  4") ;
    hist_data_uncorr_5 = (TH1D*) hist_data_5 ->Clone();
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
  NormalizeHist(hist_true_RES_Delta, MCNormalization);
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
      NormalizeHist(h_RES_Delta_slices[l], MCNormalization );
      NormalizeHist(h_RES_slices[l], MCNormalization );
      NormalizeHist(h_SIS_slices[l], MCNormalization );
      NormalizeHist(h_DIS_slices[l], MCNormalization );
      NormalizeHist(h_MEC_slices[l], MCNormalization );

      std::vector<TH1D*> all_slices{h_total_slices[l]};
      if( plot_data ) all_slices.push_back(h_data_slices[l]);
      double y_max_total = GetMaximum( all_slices );

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

      StandardFormat( h_total_slices[l], title_subname, kBlack, 1, observable, y_max_total ) ;
      StandardFormat( h_QEL_slices[l], title_subname, kBlue-3, 1, observable, y_max_total ) ;
      StandardFormat( h_RES_Delta_slices[l], title_subname, kRed-4, 1, observable, y_max_total ) ;
      StandardFormat( h_RES_slices[l], title_subname, kGreen+2, 1, observable, y_max_total ) ;
      StandardFormat( h_SIS_slices[l], title_subname, kOrange, 1, observable, y_max_total ) ;
      StandardFormat( h_MEC_slices[l], title_subname, kMagenta-3, 1, observable, y_max_total ) ;
      StandardFormat( h_DIS_slices[l], title_subname, kCyan+1, 1, observable, y_max_total ) ;
      if( plot_data ) StandardFormat( h_data_slices[l], title_subname, kBlack, 8, observable, y_max_total ) ;

    }
  }

  // Find absolute y max
  std::vector<TH1D*> temp_check = {hist_true};
  if(plot_data) temp_check.push_back(hist_data);

  for( unsigned int id = 0 ; id < hists_true_submodel.size() ; ++id ){
    temp_check.push_back(hists_true_submodel[id]);
  }
  double y_max_total = GetMaximum( temp_check );
  // Format plots
  if( plot_data ) { 
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
  }

  StandardFormat( hist_true, title, kBlack, 1, observable, y_max_total ) ;
  StandardFormat( hist_true_0, title+" Sector  0", kBlack, 1, observable ) ;
  StandardFormat( hist_true_1, title+" Sector  1", kBlack, 1, observable ) ;
  StandardFormat( hist_true_2, title+" Sector  2", kBlack, 1, observable ) ;
  StandardFormat( hist_true_3, title+" Sector  3", kBlack, 1, observable ) ;
  StandardFormat( hist_true_4, title+" Sector  4", kBlack, 1, observable ) ;
  StandardFormat( hist_true_5, title+" Sector  5", kBlack, 1, observable ) ;

  StandardFormat( hist_true_QEL, title, kBlue-3, 1, observable ) ;
  StandardFormat( hist_true_RES_Delta, title, kRed-4, 1, observable ) ;
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
  if( plot_data ) hist_data -> Draw(" err same ");
 
  if( observable=="ECal"){
    // Add a sub-pad1
    TPad * sub_pad = new TPad("subpad","",0.2,0.2,0.85,0.85);
    sub_pad->SetFillStyle(4000);
    sub_pad->Draw();
    sub_pad->cd();
    sub_pad -> SetBottomMargin(0.15);
    sub_pad -> SetLeftMargin(0.15);

    TH1D* tmp_hist_true = (TH1D*)hist_true->Clone();
    TH1D* tmp_hist_data = (TH1D*)hist_data->Clone();
    tmp_hist_true->SetTitle("");

    //tmp_hist_true->GetXaxis()->SetRangeUser(0,BeamE*(1-0.1));
    tmp_hist_data->GetXaxis()->SetRangeUser(0,BeamE*(1-0.02));
    tmp_hist_true->GetYaxis()->SetRangeUser(0,tmp_hist_data->GetBinContent(tmp_hist_data->GetMaximumBin())*(1+0.25));

    tmp_hist_true -> Draw("hist");
    hist_true_QEL -> Draw("hist same");
    hist_true_RES -> Draw("hist same");
    hist_true_SIS -> Draw("hist same");
    hist_true_MEC -> Draw("hist same");
    hist_true_DIS -> Draw("hist same");
    for( unsigned int id = 0 ; id < hists_true_submodel.size() ; ++id ){
      hists_true_submodel[id] -> SetLineWidth(3);
      hists_true_submodel[id] -> Draw("hist same");
    }
    if( plot_data ) hist_data -> Draw(" err same ");

  }
  std::string output_name = output_file_name+"_dxsec_d"+observable ;
  
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
      h_RES_Delta_slices[l]->Draw("hist same ");
      h_RES_slices[l]->Draw("hist same ");
      h_SIS_slices[l]->Draw("hist same ");
      h_MEC_slices[l]->Draw("hist same ");
      h_DIS_slices[l]->Draw("hist same ");
      if( plot_data ) h_data_slices[l]->Draw("err same ");

    }
  }
  if( addbinning.size() != 0 ) {
    output_name = output_file_name+"_dxsec_d"+observable+"_"+GetAlternativeObs(observable)+"_Slices" ;
    c_slices->SaveAs((output_location+"/TotalXSec/"+output_name+".root").c_str());
    c_slices->SaveAs((output_location+"/TotalXSec/"+output_name+".pdf").c_str());
  }
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
  hist_true_0 -> Draw("hist");
  if(plot_data) {
    hist_data_0 -> SetMarkerSize(0.7);
    hist_data_0 -> Draw(" err same ");
    //hist_data_uncorr_0 -> Draw(" err same ");
  }

  TPad *pad_sector_1 = (TPad*)pad_sector->cd(2);
  pad_sector_1 -> cd();
  pad_sector_1 -> SetBottomMargin(0.15);
  pad_sector_1 -> SetLeftMargin(0.15);
  hist_true_1 -> GetYaxis()->SetTitleOffset(1.2);
  hist_true_1 -> Draw("hist");
  if(plot_data) {
    hist_data_1 -> SetMarkerSize(0.7);
    hist_data_1 -> Draw(" err same ");
    //hist_data_uncorr_1 -> Draw(" err same ");
  }

  TPad *pad_sector_2 = (TPad*)pad_sector->cd(3);
  pad_sector_2 -> cd();
  pad_sector_2 -> SetBottomMargin(0.15);
  pad_sector_2 -> SetLeftMargin(0.15);
  hist_true_2 -> GetYaxis()->SetTitleOffset(1.2);
  hist_true_2 -> Draw("hist");
  if(plot_data) {
    hist_data_2 -> SetMarkerSize(0.7);
    hist_data_2 -> Draw(" err same ");
    //hist_data_uncorr_2 -> Draw(" err same ");
  }

  TPad *pad_sector_3 = (TPad*)pad_sector->cd(4);
  pad_sector_3 -> cd();
  pad_sector_3 -> SetBottomMargin(0.15);
  pad_sector_3 -> SetLeftMargin(0.15);
  hist_true_3 -> GetYaxis()->SetTitleOffset(1.2);
  hist_true_3 -> Draw("hist");
  if(plot_data) {
    hist_data_3 -> SetMarkerSize(0.7);
    hist_data_3 -> Draw(" err same ");
    //hist_data_uncorr_3 -> Draw(" err same ");
  }

  TPad *pad_sector_4 = (TPad*)pad_sector->cd(5);
  pad_sector_4 -> cd();
  pad_sector_4 -> SetBottomMargin(0.15);
  pad_sector_4 -> SetLeftMargin(0.15);
  hist_true_4 -> GetYaxis()->SetTitleOffset(1.2);
  hist_true_4 -> Draw("hist");
  if(plot_data) {
    hist_data_4 -> SetMarkerSize(0.7);
    hist_data_4 -> Draw(" err same ");
    //hist_data_uncorr_4 -> Draw(" err same ");
  }

  TPad *pad_sector_5 = (TPad*)pad_sector->cd(6);
  pad_sector_5 -> cd();
  pad_sector_5 -> SetBottomMargin(0.15);
  pad_sector_5 -> SetLeftMargin(0.15);
  hist_true_5 -> GetYaxis()->SetTitleOffset(1.2);
  hist_true_5 -> Draw("hist");
  if(plot_data) {
    hist_data_5 -> SetMarkerSize(0.7);
    hist_data_5 -> Draw(" err same ");
    //hist_data_uncorr_5 -> Draw(" err same ");
  }

  output_name = MC_files_name[0]+"_dxsec_d"+observable+"_persector" ;
  std::filesystem::path xsecpersector_path{(output_location+"/XSecPerSector/").c_str()};
  if( ! std::filesystem::exists(xsecpersector_path) ) std::filesystem::create_directory(xsecpersector_path);
  c_sector->SaveAs((output_location+"/XSecPerSector/"+output_name+".root").c_str());
  c_sector->SaveAs((output_location+"/XSecPerSector/"+output_name+".pdf").c_str());
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
  leg->AddEntry(hist_true_RES_Delta,(model[0]+" EMRES P33(1232)").c_str(),"l");
  leg->AddEntry(hist_true_RES,(model[0]+" EMRES Others").c_str(),"l");
  leg->AddEntry(hist_true_SIS,(model[0]+" EMSIS").c_str(),"l");
  leg->AddEntry(hist_true_MEC,(model[0]+" EMMEC").c_str(),"l");
  leg->AddEntry(hist_true_DIS,(model[0]+" EMDIS").c_str(),"l");

  if( model.size() > 1 ) {
    for( unsigned int id = 1 ; id < model.size() -1 ; ++id ){
      leg->AddEntry(hists_true_submodel[id-1],("GENIE "+model[id]).c_str(),"l");
    }
  }

  if(plot_data)  leg->AddEntry(hist_data, data_name.c_str(), "lp");
  leg->Draw();
  output_name = MC_files_name[0] ;
  c_leg->SaveAs((output_location+"/"+output_name+"_legend.root").c_str());
  c_leg->SaveAs((output_location+"/"+output_name+"_legend.pdf").c_str());

  delete c_leg;
}
