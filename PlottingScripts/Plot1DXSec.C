#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"


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

void CorrectData(TH1D* h, TH1D* acc) {

  double NBins = h->GetNbinsX();
  for (int i = 1; i <= NBins; i++) {
    double content = h->GetBinContent(i);
    double error = h->GetBinError(i);
    if( h->GetBinContent(i) != 0 ) h->SetBinContent(i, content * acc->GetBinContent(i));
  }
}

void NormalizeData( TH1D * h, unsigned int target_pdg = 1000060120, unsigned int beam_E = 1.161 ){
    // Data normalization
    const double ConversionFactorCm2ToMicroBarn = TMath::Power(10.,30.);
    double ConversionFactorChargeToElectrons = 6.25*TMath::Power(10.,15.);
    double AvogadroNumber = 6.02*TMath::Power(10.,23);
    double OverallUnitConversionFactor = ConversionFactorChargeToElectrons * AvogadroNumber;
    int MassNumber = 12 ;
    double integrated_ch = 0.079 ;
    double lenght = 0.1;
    double density = 1.786;

    if( beam_E == 2.216 && target_pdg == 1000060120) {
      integrated_ch = 1.79 ;
    }

    double data_scaling = ConversionFactorCm2ToMicroBarn / ( integrated_ch * lenght * density * OverallUnitConversionFactor / MassNumber ) ;

    h -> Scale( data_scaling );
    ReweightPlots(h);
}


void NormalizeMC( TH1D * h, unsigned int target_pdg = 1000060120, unsigned int beam_E = 1.161, double N_events = 19800000 ){
    // Data normalization
    const double ConversionFactorCm2ToMicroBarn = TMath::Power(10.,30.);
    double xsec_susa = 1.28967e+09 ;

    if( beam_E == 2.216 && target_pdg == 1000060120) {
      xsec_susa = 2.1024e+08 ;
    }

    double mc_scaling = xsec_susa * ConversionFactorCm2ToMicroBarn * TMath::Power(10.,-38.) / N_events ;

    h -> Scale( mc_scaling );
    ReweightPlots(h);
}

// Input paramters:
// MC_file_name : true MC file, without detector effects, after e4nu analysis
// data_file_name: data file, after e4nu Analysis
// acceptance_file_name: acceptance file obtained with compute_acceptance.C
// target target_pdg
// beam energy
// Number of events in original MC file (pre-analysis)
// id_secotr: ID of the sector of interest

void Plot1DXSec(std::string MC_file_name, std::string data_file_name,
                std::string acceptance_file_name, std::string observable,
                unsigned int target_pdg, unsigned int E_beam, long N_events,
                int id_sector = -1 /*all*/ ) {

  TCanvas * c1 = new TCanvas("c1","c1",800,600);

  TFile * file_true_MC = new TFile((MC_file_name).c_str(),"ROOT");
  TFile * file_data = new TFile((data_file_name).c_str(),"READ");
  TFile * file_acceptance = new TFile((acceptance_file_name).c_str(),"READ");

  if( !file_true_MC ) { std::cout << "ERROR: the "<< file_true_MC << " does not exist." <<std::endl; return ;}
  if( !file_data ) { std::cout << "ERROR: the "<< file_data << " does not exist." <<std::endl; return ;}
  if( !file_acceptance ) { std::cout << "ERROR: the "<< file_acceptance << " does not exist." <<std::endl; return ;}

  TH1D * h_acceptance = (TH1D*)file_acceptance->Get("Acceptance");
  TTree * tree_true = (TTree*)file_true_MC->Get("MCCLAS6Tree");
  TTree * tree_data = (TTree*)file_data->Get("CLAS6Tree");

  if( !h_acceptance ) { std::cout << "ERROR: Acceptance is not defined"<<std::endl; return ; }
  if( !tree_true || !tree_data ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return ;}

  TH1D * hist_true = (TH1D*) h_acceptance ->Clone();
  hist_true -> SetName( "MC True") ;
  hist_true -> Reset("ICE");
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
  TH1D * hist_data = (TH1D*) h_acceptance ->Clone();
  hist_data -> SetName( "Data") ;
  hist_data -> Reset("ICE");

  std::vector<TTree*> trees = { tree_true, tree_data };
  std::vector<TH1D*> hists = { hist_true, hist_data };

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
  bool QEL, RES, DIS, MEC;
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

    trees[i] -> SetBranchAddress("QEL",&QEL);
    trees[i] -> SetBranchAddress("RES",&RES);
    trees[i] -> SetBranchAddress("MEC",&MEC);
    trees[i] -> SetBranchAddress("DIS",&DIS);

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

      if( i == 0 ){
        if( QEL ) hist_true_QEL -> Fill( content, w ) ;
        if( RES ) {
          hist_true_RES -> Fill( content, w ) ;
          hist_true_SIS -> Fill( content, w ) ;
        }
        if( DIS ) hist_true_DIS -> Fill( content, w ) ;
        if( MEC ) hist_true_MEC -> Fill( content, w ) ;
      }
    }
  }

  CorrectData(hist_data, h_acceptance);
  NormalizeData(hist_data, target_pdg, E_beam );
  NormalizeMC(hist_true, target_pdg, E_beam, N_events );
  NormalizeMC(hist_true_QEL, target_pdg, E_beam, N_events );
  NormalizeMC(hist_true_RES, target_pdg, E_beam, N_events );
  NormalizeMC(hist_true_SIS, target_pdg, E_beam, N_events );
  NormalizeMC(hist_true_MEC, target_pdg, E_beam, N_events );
  NormalizeMC(hist_true_DIS, target_pdg, E_beam, N_events );

  hist_data->SetLineColor(kBlack);
  hist_true->SetLineColor(kBlack);
  hist_true_QEL->SetLineColor(kBlue-3);
  hist_true_RES->SetLineColor(kGreen+2);
  hist_true_SIS->SetLineColor(kOrange);
  hist_true_MEC->SetLineColor(kMagenta-3);
  hist_true_DIS->SetLineColor(kCyan+1);

  hist_true -> Draw("hist");
  hist_true_QEL -> Draw("hist same");
  hist_true_RES -> Draw("hist same");
  hist_true_SIS -> Draw("hist same");
  hist_true_MEC -> Draw("hist same");
  hist_true_DIS -> Draw("hist same");
  hist_data -> Draw(" same ");

  std::string output_name = "e4nuanalysis_"+std::to_string(target_pdg)+"_"+std::to_string(E_beam)+"GeV_dxsec_d"+observable ;
  if( id_sector > 0 ) output_name += "_sector_"+std::to_string(id_sector) ;
  c1->SaveAs((output_name+".root").c_str());
}

void Plot1DXSec(){

  // EDIT :
  std::string file_data = "e4nuanalysis_clas6data_C12_2261MeV_1p0pi.root";
  std::string file_mc = "e4nuanalysis_genie_GEM21_C12_2261MeV_1p0piwpi0_1M_Q4_wreso_true.root" ;
  std::string acceptance_file = "e4nuanalysis_genie_GEM21_C12_2261MeV_1p0piwpi0_1M_Q4_wreso_acceptance_correction_ECal.root";
  std::string observable = "ECal";

  Plot1DXSec( file_mc, file_data, acceptance_file, "ECal", 1000060120, 1.161, 19800000, -1 ) ;
}
