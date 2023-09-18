#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"

void compute_2Dacceptance( std::string file_name, std::string observable_1, unsigned int nbins_1, double min_1, double max_1,
			   std::string observable_2, unsigned int nbins_2, double min_2, double max_2, int id_sector = -1 /*all*/ ) {

  TCanvas * c1 = new TCanvas("c1","c1",800,600);

  TFile * file_mcrecoacc = new TFile((file_name+"_truereco.root").c_str(),"ROOT");
  TFile * file_mctrueacc = new TFile((file_name+"_true.root").c_str(),"READ");

  if( !file_mcrecoacc ) { std::cout << "ERROR: the "<< file_name << "_truereco.root does not exist." <<std::endl; return ;}
  if( !file_mctrueacc ) { std::cout << "ERROR: the "<< file_name << "_true.root  does not exist." <<std::endl; return ;}

  TTree * tree_mcrecoacc = (TTree*)file_mcrecoacc->Get("MCCLAS6Tree");
  TTree * tree_mctrueacc = (TTree*)file_mctrueacc->Get("MCCLAS6Tree");

  if( !tree_mcrecoacc || !tree_mctrueacc ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return ;}

  TH2D * hist_recoacc = new TH2D( "Reco MC ACC", "", nbins_1, min_1, max_1, nbins_2, min_2, max_2 );
  TH2D * hist_trueacc = new TH2D( "True MC ACC", "", nbins_1, min_1, max_1, nbins_2, min_2, max_2 );

  std::vector<TTree*> trees = { tree_mcrecoacc, tree_mctrueacc };
  std::vector<TH2D*> hists = { hist_recoacc, hist_trueacc };

  // OBSERVABLE DEFINITION:
  double TotWeight ;
  double ECal,Recoq3,RecoW;
  double pfl,pfl_theta,pfl_phi;
  double proton_mom,proton_phi,proton_theta;
  double pi_mom,pi_theta,pi_phi;
  double HadAlphaT, HadDeltaPT, HadDeltaPhiT ;
  long NEntries ;
  bool IsBkg ;
  int ElectronSector ;
  for ( unsigned int i = 0 ; i < trees.size() ; ++i ){

    NEntries = trees[i] -> GetEntries() ;
    trees[i] -> SetBranchAddress("TotWeight",&TotWeight);
    trees[i] -> SetBranchAddress("IsBkg",&IsBkg);
    trees[i] -> SetBranchAddress("pfl",&pfl);
    trees[i] -> SetBranchAddress("pfl_theta",&pfl_theta);
    trees[i] -> SetBranchAddress("pfl_phi",&pfl_phi);

    trees[i] -> SetBranchAddress("proton_mom",&proton_mom);
    trees[i] -> SetBranchAddress("proton_theta",&proton_theta);
    trees[i] -> SetBranchAddress("proton_phi",&proton_phi);

    trees[i] -> SetBranchAddress("pim_mom",&pi_mom);
    trees[i] -> SetBranchAddress("pim_theta",&pi_theta);
    trees[i] -> SetBranchAddress("pim_phi",&pi_phi);

    trees[i] -> SetBranchAddress("ECal",&ECal);
    trees[i] -> SetBranchAddress("Recoq3",&Recoq3);
    trees[i] -> SetBranchAddress("RecoW",&RecoW);

    trees[i] -> SetBranchAddress("ElectronSector",&ElectronSector);


    for( int j = 0 ; j < NEntries ; ++j ) {
      trees[i]->GetEntry(j) ;
      double content_1 = 0 ;
      double content_2 = 0 ;
      double w = TotWeight ;

      if( observable_1 == "ECal") content_1 = ECal ;
      else if ( observable_1 == "pfl_theta") content_1 = pfl_theta * 180 / TMath::Pi();
      else if ( observable_1 == "pfl_phi") content_1 = pfl_phi * 180 / TMath::Pi();
      else if ( observable_1 == "pfl") content_1 = pfl ;
      else if ( observable_1 == "proton_mom") content_1 = proton_mom ;
      else if ( observable_1 == "proton_theta") content_1 = proton_theta * 180 / TMath::Pi() ;
      else if ( observable_1 == "proton_phi") content_1 = proton_phi * 180 / TMath::Pi() ;
      else if ( observable_1 == "pim_mom") content_1 = pi_mom ;
      else if ( observable_1 == "pim_theta") content_1 = pi_theta * 180 / TMath::Pi() ;
      else if ( observable_1 == "RecoW") content_1 = RecoW ;
      else if ( observable_1 == "Recoq3") content_1 = Recoq3 ;

      if( observable_2 == "ECal") content_2 = ECal ;
      else if ( observable_2 == "pfl_theta") content_2 = pfl_theta * 180 / TMath::Pi();
      else if ( observable_2 == "pfl_phi") content_2 = pfl_phi * 180 / TMath::Pi();
      else if ( observable_2 == "pfl") content_2 = pfl ;
      else if ( observable_2 == "proton_mom") content_2 = proton_mom ;
      else if ( observable_2 == "proton_theta") content_2 = proton_theta * 180 / TMath::Pi() ;
      else if ( observable_2 == "proton_phi") content_2 = proton_phi * 180 / TMath::Pi() ;
      else if ( observable_2 == "pim_mom") content_2 = pi_mom ;
      else if ( observable_2 == "pim_theta") content_2 = pi_theta * 180 / TMath::Pi() ;
      else if ( observable_2 == "RecoW") content_2 = RecoW ;
      else if ( observable_2 == "Recoq3") content_2 = Recoq3 ;

      if( id_sector > 0 ) {
	// Compute only for sector of interest
	if( id_sector != ElectronSector ) continue ;
      }
      hists[i] -> Fill( content_1, content_2, w ) ;
      hists[i] -> SetLineWidth(3);
    }
  }

  TH2D * ratio = (TH2D*)hist_trueacc->Clone();
  ratio -> Divide(hist_recoacc);
  ratio -> GetXaxis()->SetTitle(observable_1.c_str());
  ratio -> GetYaxis()->SetTitle(observable_2.c_str());
  ratio -> GetZaxis()->SetTitle("Acceptance correction");
  ratio->Draw("COLZ");

  std::string output_name = file_name+"_2Dacceptance_correction_"+observable_1+"_vs_"+observable_2 ;
  if( id_sector > 0 ) output_name += "_sector_"+std::to_string(id_sector) ;
  c1->SaveAs((output_name+".root").c_str());
}

void compute_2Dacceptance(){

  // EDIT :
  std::string file_name = "e4nuanalysis_genie_GEM21_C12_2261MeV_1p0piwpi0_1M_Q4" ;
  std::string observable_1 = "proton_mom";
  int nbins_1 = 30 ;
  double min_1 = 0;
  double max_1 = 2.4;

  std::string observable_2 = "proton_theta";
  int nbins_2 = 30 ;
  double min_2 = 20;
  double max_2 = 120;

  compute_2Dacceptance( file_name, observable_1, nbins_1, min_1, max_1, observable_2, nbins_2, min_2, max_2 ) ;
}
