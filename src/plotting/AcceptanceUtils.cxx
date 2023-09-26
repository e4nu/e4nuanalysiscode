#include <iomanip>
#include <filesystem>
#include "plotting/AcceptanceUtils.h"

using namespace e4nu ;
using namespace e4nu::plotting ;

std::string plotting::ComputeAcceptance(std::vector<std::string> mc_files, std::string observable, std::string title,
					std::string input_MC_location, std::string output_location,  std::string output_file_name, 
					std::string analysis_id, bool store_root ) {

  // Define trees
  std::vector<TFile*> files_mcrecoacc, files_mctrueacc;
  std::vector<TTree*> trees_mcrecoacc, trees_mctrueacc ;

  // Define Hists
  // The _# correspond to histograms for each sector
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
    binning = plotting::GetBinning(observable,BeamE,analysis_id);
    
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

    std::vector<double> addbinning = GetAdditionalBinning( GetAlternativeObs(observable), BeamE, analysis_id );
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
    double HadAlphaT, HadDeltaPT, HadDeltaPTx, HadDeltaPTy, HadDeltaPhiT ;
    double AlphaT, DeltaPT, DeltaPhiT ;
    double RecoXBJK, RecoEnergyTransfer, RecoQ2, HadSystemMass, RecoQELEnu ;
    double MissingEnergy, MissingAngle, MissingMomentum ;
    double InferedNucleonMom ;
    double HadronsAngle, Angleqvshad;
    double AdlerAngleThetaP, AdlerAnglePhiP, AdlerAngleThetaPi, AdlerAnglePhiPi ; 
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
      trees[j] -> SetBranchAddress("HadDeltaPTx",&HadDeltaPTx);
      trees[j] -> SetBranchAddress("HadDeltaPTy",&HadDeltaPTy);
      trees[j] -> SetBranchAddress("DeltaPhiT",&DeltaPhiT);
      trees[j] -> SetBranchAddress("HadDeltaPhiT",&HadDeltaPhiT);
      trees[j] -> SetBranchAddress("ElectronSector",&ElectronSector);
      trees[j] -> SetBranchAddress("HadSystemMass", &HadSystemMass);
      trees[j] -> SetBranchAddress("MissingEnergy", &MissingEnergy);
      trees[j] -> SetBranchAddress("MissingAngle", &MissingAngle);
      trees[j] -> SetBranchAddress("MissingMomentum", &MissingMomentum);
      trees[j] -> SetBranchAddress("InferedNucleonMom", &InferedNucleonMom);
      trees[j] -> SetBranchAddress("HadronsAngle",&HadronsAngle);
      trees[j] -> SetBranchAddress("AdlerAngleThetaP",&AdlerAngleThetaP);
      trees[j] -> SetBranchAddress("AdlerAnglePhiP",&AdlerAnglePhiP);
      trees[j] -> SetBranchAddress("AdlerAngleThetaPi",&AdlerAngleThetaPi);
      trees[j] -> SetBranchAddress("AdlerAnglePhiPi",&AdlerAnglePhiPi);
      trees[j] -> SetBranchAddress("Angleqvshad",&Angleqvshad);

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
	else if ( observable == "HadDeltaPTx") content = HadDeltaPTx ;
	else if ( observable == "HadDeltaPTy") content = HadDeltaPTy ;
        else if ( observable == "DeltaPhiT") content = DeltaPhiT ;
        else if ( observable == "HadDeltaPhiT") content = HadDeltaPhiT ;
        else if ( observable == "HadSystemMass") content = HadSystemMass ;
	else if ( observable == "MissingEnergy") content = MissingEnergy ;
	else if ( observable == "MissingMomentum") content = MissingMomentum ;
	else if ( observable == "MissingAngle") content = MissingAngle ;
	else if ( observable == "InferedNucleonMom") content = InferedNucleonMom ;
	else if ( observable == "HadronsAngle" ) content = HadronsAngle ; 
	else if ( observable == "AdlerAngleThetaP" ) content = AdlerAngleThetaP ; 
	else if ( observable == "AdlerAnglePhiP" ) content = AdlerAnglePhiP ; 
	else if ( observable == "AdlerAngleThetaPi" ) content = AdlerAngleThetaPi ; 
	else if ( observable == "AdlerAnglePhiPi" ) content = AdlerAnglePhiPi ; 
	else if ( observable == "Angleqvshad" ) content = Angleqvshad ; 

        // Fill the per Sector  histogram
        hists[2*(ElectronSector+1)+(j-initial_size_trees)+initial_size_hists] -> Fill( content, w ) ;
        hists[2*(ElectronSector+1)+(j-initial_size_trees)+initial_size_hists] -> SetLineWidth(3);

	hists[j+initial_size_hists-initial_size_trees] -> Fill( content, w ) ;
        hists[j+initial_size_hists-initial_size_trees] -> SetLineWidth(3);

	std::string alt_obs = GetAlternativeObs(observable) ;
	double content_2 = 0 ;
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

  std::vector<double> addbinning = GetAdditionalBinning( GetAlternativeObs(observable), BeamE, analysis_id ) ;
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

  if( store_root ) c_1->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_total.root").c_str());
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

  if( store_root ) c_sector_2->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_persector.root").c_str());
  c_sector_2->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_persector.pdf").c_str());
  delete c_sector_2 ;

  for (size_t i = 0; i < files_mcrecoacc.size(); i++) {
    delete files_mcrecoacc[i] ;
    delete files_mctrueacc[i] ;
  }
  outputFile.Close();
  return acc_file ;
}
