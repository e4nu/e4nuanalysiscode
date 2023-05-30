#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"

void compute_acceptance(std::string file_name, std::string observable, unsigned int nbins, double min, double max, int id_sector = -1 /*all*/ ) {

  TCanvas * c1 = new TCanvas("c1","c1",800,600);

  TFile * file_mcrecoacc = new TFile((file_name+"_truereco.root").c_str(),"ROOT");
  TFile * file_mctrueacc = new TFile((file_name+"_true.root").c_str(),"READ");

  if( !file_mcrecoacc ) { std::cout << "ERROR: the "<< file_name << "_truereco.root does not exist." <<std::endl; return ;}
  if( !file_mctrueacc ) { std::cout << "ERROR: the "<< file_name << "_true.root  does not exist." <<std::endl; return ;}

  TTree * tree_mcrecoacc = (TTree*)file_mcrecoacc->Get("MCCLAS6Tree");
  TTree * tree_mctrueacc = (TTree*)file_mctrueacc->Get("MCCLAS6Tree");

  if( !tree_mcrecoacc || !tree_mctrueacc ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return ;}

  TH1D * hist_recoacc = new TH1D( "Reco MC ACC", "", nbins, min, max );
  TH1D * hist_trueacc = new TH1D( "True MC ACC", "", nbins, min, max );

  std::vector<TTree*> trees = { tree_mcrecoacc, tree_mctrueacc };
  std::vector<TH1D*> hists = { hist_recoacc, hist_trueacc };

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
      double content = 0 ;
      double w = TotWeight ;

      if( observable == "ECal") content = ECal ;
      else if ( observable == "pfl_theta") content = pfl_theta * 180 / TMath::Pi();
      else if ( observable == "pfl_phi") content = pfl_phi * 180 / TMath::Pi();
      else if ( observable == "pfl") content = pfl ;
      else if ( observable == "proton_mom") content = proton_mom ;
      else if ( observable == "proton_theta") content = proton_theta * 180 / TMath::Pi() ;
      else if ( observable == "proton_phi") content = proton_phi * 180 / TMath::Pi() ;
      else if ( observable == "pim_mom") content = pi_mom ;
      else if ( observable == "pim_theta") content = pi_theta * 180 / TMath::Pi() ;
      else if ( observable == "RecoW") content = RecoW ;
      else if ( observable == "Recoq3") content = Recoq3 ;

      if( id_sector > 0 ) {
	// Compute only for sector of interest
	if( id_sector != ElectronSector ) continue ;
      }
      hists[i] -> Fill( content, w ) ;
      hists[i] -> SetLineWidth(3);
    }
  }

  TH1D * ratio = (TH1D*)hist_trueacc->Clone();
  ratio -> Divide(hist_recoacc);
  ratio -> GetXaxis()->SetTitle(observable.c_str());
  ratio -> GetYaxis()->SetTitle("Acceptance correction");
  ratio->Draw("hist err");

  std::string output_name = file_name+"_acceptance_correction_"+observable ;
  if( id_sector > 0 ) output_name += "_sector_"+std::to_string(id_sector) ;
  c1->SaveAs((output_name+".root").c_str());
}

void compute_acceptance(){

  // EDIT :
  std::string file_name = "e4nuanalysis_genie_GEM21_C12_2261MeV_1p0piwpi0_1M_Q4" ;
  std::string observable = "proton_mom";
  int nbins = 60 ;
  double min = 0;
  double max = 2.4;

  compute_acceptance( file_name, observable, nbins, min, max ) ;
}
