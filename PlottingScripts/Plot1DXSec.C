#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TTree.h"

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
	else if ( observable == "RecoXBJK") { x_axis = "x_{BJK} [GeV]"; y_axis  = "d#sigma/dx_{BJK}} #left[#mub GeV^{-1}#right]"; }
	else if ( observable == "RecoQ2") { x_axis = "Q^{2} [GeV]"; y_axis  = "d#sigma/dQ^{2}} #left[#mub GeV^{-1}#right]"; }
	else if ( observable == "Recoq3") { x_axis = "q_{3} [GeV]"; y_axis  = "d#sigma/dq_{3} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "DeltaPT") { x_axis = "#deltap_{T} [GeV]"; y_axis  = "d#sigma/d#deltap_{T} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "HadDeltaPT") { x_axis = "#deltap_{T}^{had} [GeV]"; y_axis  = "d#sigma/d#deltap_{T}^{had} #left[#mub #left(GeV/c#right)^{-1}#right]"; }
	else if ( observable == "DeltaPhiT") { x_axis = "#delta#phi_{T} [deg]"; y_axis  = "d#sigma/d#delta#phi_{T} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "HadDeltaPhiT") { x_axis = "#delta#phi_{T}^{had} [deg]"; y_axis  = "d#sigma/d#delta#phi_{T}^{had} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "AlphaT") { x_axis = "#alpha_{T} [deg]"; y_axis  = "d#sigma/d#alpha_{T} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "HadAlphaT") { x_axis = "#alpha_{T}^{had} [deg]"; y_axis  = "d#sigma/d#alpha_{T}^{had} #left[#mub deg^{-1}#right]"; }
	else if ( observable == "RecoEnergyTransfer") { x_axis = "#omega [GeV]"; y_axis  = "d#sigma/d#omega #left[#mub GeV^{-1}#right]"; }

	if( id_axis ==0 ) return x_axis ;
	return y_axis ;
}

void StandardFormat( TH1D * prediction, std::string title, int color, int style, std::string observable ) {
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
	prediction->SetTitleFont(FontStyle);

  return;
}

std::string compute_acceptance(std::vector<std::string> mc_files, std::string observable, std::string title,
                               std::vector<double> binning,
                               std::string input_MC_location, std::string output_location,
                               int id_sector = -1 /*all*/ ) {

  std::vector<TFile*> files_mcrecoacc, files_mctrueacc;
  std::vector<TTree*> trees_mcrecoacc, trees_mctrueacc ;
  std::vector<TH1D*>  hists_recoacc, hists_trueacc, hists_recoacc_0, hists_trueacc_0, hists_recoacc_1, hists_trueacc_1,
                      hists_recoacc_2, hists_trueacc_2, hists_recoacc_3, hists_trueacc_3,
                      hists_recoacc_4, hists_trueacc_4, hists_recoacc_5, hists_trueacc_5 ;
  std::vector<TTree*> trees;
  std::vector<TH1D*>  hists, ratios, ratios_0, ratios_1, ratios_2, ratios_3, ratios_4, ratios_5 ;

  for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
    files_mcrecoacc.push_back(new TFile((input_MC_location+mc_files[i]+"_truereco.root").c_str(),"ROOT"));
    files_mctrueacc.push_back(new TFile((input_MC_location+mc_files[i]+"_true.root").c_str(),"ROOT"));
    if( !files_mcrecoacc[i] ) { std::cout << "ERROR: the "<< mc_files[i] << "_truereco.root does not exist." <<std::endl; return "";}
    if( !files_mctrueacc[i] ) { std::cout << "ERROR: the "<< mc_files[i] << "_true.root  does not exist." <<std::endl; return "";}
    trees_mcrecoacc.push_back( (TTree*)files_mcrecoacc[i]->Get("MCCLAS6Tree"));
    trees_mctrueacc.push_back( (TTree*)files_mctrueacc[i]->Get("MCCLAS6Tree"));
    if( !trees_mctrueacc[i] || !trees_mcrecoacc[i] ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return "";}

    hists_recoacc.push_back( new TH1D( ("Reco MC ACC Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc.push_back( new TH1D( ("True MC ACC Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_0.push_back( new TH1D( ("Reco MC ACC Sector 0 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_0.push_back( new TH1D( ("True MC ACC Sector 0 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_1.push_back( new TH1D( ("Reco MC ACC Sector 1 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_1.push_back( new TH1D( ("True MC ACC Sector 1 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_2.push_back( new TH1D( ("Reco MC ACC Sector 2 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_2.push_back( new TH1D( ("True MC ACC Sector 2 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_3.push_back( new TH1D( ("Reco MC ACC Sector 3 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_3.push_back( new TH1D( ("True MC ACC Sector 3 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_4.push_back( new TH1D( ("Reco MC ACC Sector 4 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_4.push_back( new TH1D( ("True MC ACC Sector 4 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_recoacc_5.push_back( new TH1D( ("Reco MC ACC Sector 5 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
    hists_trueacc_5.push_back( new TH1D( ("True MC ACC Sector 5 Model "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;

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

    // OBSERVABLE DEFINITION:
    double TotWeight ;
    double ECal,Recoq3,RecoW;
    double pfl,pfl_theta,pfl_phi;
    double proton_mom,proton_phi,proton_theta;
    double pim_mom,pim_theta,pim_phi;
    double pip_mom,pip_theta,pip_phi;
    double HadAlphaT, HadDeltaPT, HadDeltaPhiT ;
    double AlphaT, DeltaPT, DeltaPhiT ;
    double RecoXBJK, RecoEnergyTransfer, RecoQ2 ;
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
        else if ( observable == "RecoXBJK") content = RecoXBJK ;
        else if ( observable == "RecoQ2") content = RecoQ2 ;
        else if ( observable == "RecoEnergyTransfer") content = RecoEnergyTransfer ;
        else if ( observable == "AlphaT") content = AlphaT ;
        else if ( observable == "HadAlphaT") content = HadAlphaT ;
        else if ( observable == "DeltaPT") content = DeltaPT ;
        else if ( observable == "HadDeltaPT") content = HadDeltaPT ;
        else if ( observable == "DeltaPhiT") content = DeltaPhiT ;
        else if ( observable == "HadDeltaPhiT") content = HadDeltaPhiT ;

        // Fill the per sector histogram
        hists[2*(ElectronSector+1)+(j-initial_size_trees)+initial_size_hists] -> Fill( content, w ) ;
        hists[2*(ElectronSector+1)+(j-initial_size_trees)+initial_size_hists] -> SetLineWidth(3);

        if( id_sector > 0 ) {
          // Compute only for sector of interest
          if( id_sector != ElectronSector ) continue ;
        }
        hists[j] -> Fill( content, w ) ;
        hists[j] -> SetLineWidth(3);
      }
    }

    ratios.push_back( (TH1D*)hists_trueacc[i]->Clone() ) ;
    ratios[i] -> Divide( hists_recoacc[i] );
    ratios[i] -> SetName(("Acceptance_model_"+std::to_string(i)).c_str());
    StandardFormat( ratios[i], title, kBlack, 2+i, observable ) ;
    ratios[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios[i] -> GetYaxis()->SetTitle("Acceptance correction");

    ratios_0.push_back( (TH1D*)hists_trueacc_0[i]->Clone() ) ;
    ratios_0[i] -> Divide( hists_recoacc_0[i] );
    ratios_0[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_0").c_str());
    StandardFormat( ratios_0[i], title, kOrange+1+i, 2+i, observable ) ;
    ratios_0[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_0[i] -> GetYaxis()->SetTitle("Acceptance correction Sector 0");

    ratios_1.push_back( (TH1D*)hists_trueacc_1[i]->Clone() ) ;
    ratios_1[i] -> Divide( hists_recoacc_1[i] );
    ratios_1[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_1").c_str());
    StandardFormat( ratios_1[i], title, kPink+4-i, 2+i, observable ) ;
    ratios_1[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_1[i] -> GetYaxis()->SetTitle("Acceptance correction Sector 1");

    ratios_2.push_back( (TH1D*)hists_trueacc_2[i]->Clone() ) ;
    ratios_2[i] -> Divide( hists_recoacc_2[i] );
    ratios_2[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_2").c_str());
    StandardFormat( ratios_2[i], title, kViolet+5-i, 2+i, observable ) ;
    ratios_2[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_2[i] -> GetYaxis()->SetTitle("Acceptance correction Sector 2");

    ratios_3.push_back( (TH1D*)hists_trueacc_3[i]->Clone() ) ;
    ratios_3[i] -> Divide( hists_recoacc_3[i] );
    ratios_3[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_3").c_str());
    StandardFormat( ratios_3[i], title, kAzure-5+i, 2+i, observable ) ;
    ratios_3[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_3[i] -> GetYaxis()->SetTitle(("Acceptance_model sector_3").c_str());

    ratios_4.push_back( (TH1D*)hists_trueacc_4[i]->Clone() ) ;
    ratios_4[i] -> Divide( hists_recoacc_4[i] );
    ratios_4[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_4").c_str());
    StandardFormat( ratios_4[i], title, kTeal-7-i, 2+i, observable ) ;
    ratios_4[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_4[i] -> GetYaxis()->SetTitle("Acceptance correction Sector 4");

    ratios_5.push_back( (TH1D*)hists_trueacc_5[i]->Clone() ) ;
    ratios_5[i] -> Divide( hists_recoacc_5[i] );
    ratios_5[i] -> SetName(("Acceptance_model_"+std::to_string(i)+"_sector_5").c_str());
    StandardFormat( ratios_5[i], title, kGreen-3-i, 2+i, observable ) ;
    ratios_5[i] -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
    ratios_5[i] -> GetYaxis()->SetTitle("Acceptance correction Sector 5");
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
  ratio_0 -> GetYaxis()->SetTitle("Acceptance correction Sector 0");

  TH1D* ratio_1 = (TH1D*)ratios_1[0]->Clone();
  ratio_1 -> SetName("Acceptance_1");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_1->Add(ratios_1[i]);
  }
  ratio_1 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_1, title, kPink+4, 1, observable ) ;
  ratio_1 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_1 -> GetYaxis()->SetTitle("Acceptance correction Sector 1");

  TH1D* ratio_2 = (TH1D*)ratios_2[0]->Clone();
  ratio_2 -> SetName("Acceptance_2");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_2->Add(ratios_2[i]);
  }
  ratio_2 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_2, title, kViolet+5, 1, observable ) ;
  ratio_2 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_2 -> GetYaxis()->SetTitle("Acceptance correction Sector 2");

  TH1D* ratio_3 = (TH1D*)ratios_3[0]->Clone();
  ratio_3 -> SetName("Acceptance_3");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_3->Add(ratios_3[i]);
  }
  ratio_3 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_3, title, kAzure-5, 1, observable ) ;
  ratio_3 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_3 -> GetYaxis()->SetTitle("Acceptance correction Sector 3");

  TH1D* ratio_4 = (TH1D*)ratios_4[0]->Clone();
  ratio_4 -> SetName("Acceptance_4");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_4->Add(ratios_4[i]);
  }
  ratio_4 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_4, title, kTeal-7, 1, observable ) ;
  ratio_4 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_4 -> GetYaxis()->SetTitle("Acceptance correction Sector 4");

  TH1D* ratio_5 = (TH1D*)ratios_5[0]->Clone();
  ratio_5 -> SetName("Acceptance_5");
  for( unsigned int i = 1 ; i < mc_files.size() ; ++i ) {
     ratio_5->Add(ratios_5[i]);
  }
  ratio_5 -> Scale( 1./mc_files.size() );
  StandardFormat( ratio_5, title, kGreen-3, 1, observable ) ;
  ratio_5 -> GetXaxis()->SetTitle(GetAxisLabel(observable,0).c_str());
  ratio_5 -> GetYaxis()->SetTitle("Acceptance correction Sector 5");

  std::string output_name = mc_files[0]+"_acceptance_correction_"+observable ;
  if( id_sector > 0 ) output_name += "_sector_"+std::to_string(id_sector) ;
  std::string acc_file = "/AcceptanceFiles/"+output_name ;
  TFile outputFile ((output_location+acc_file+".root").c_str(),"RECREATE");

  TCanvas * c_1 = new TCanvas("c_1","c_1",200,10,700,500);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  ratio->Draw("hist err");
  ratio->Write();
  ratio_0->Draw("hist err");
  ratio_0->Write();
  ratio_1->Draw("hist err");
  ratio_1->Write();
  ratio_2->Draw("hist err");
  ratio_2->Write();
  ratio_3->Draw("hist err");
  ratio_3->Write();
  ratio_4->Draw("hist err");
  ratio_4->Write();
  ratio_5->Draw("hist err");
  ratio_5->Write();
  for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
    ratios[i]->Draw("hist err");
    ratios[i]->Write();
  }

  c_1->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_total.root").c_str());
  c_1->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_total.pdf").c_str());
  delete c_1 ;

  for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) {
    ratios_0[i] -> Draw("hist err");
    ratios_0[i] -> Write();
    ratios_1[i] -> Draw("hist err");
    ratios_1[i] -> Write();
    ratios_2[i] -> Draw("hist err");
    ratios_2[i] -> Write();
    ratios_3[i] -> Draw("hist err");
    ratios_3[i] -> Write();
    ratios_4[i] -> Draw("hist err");
    ratios_4[i] -> Write();
    ratios_5[i] -> Draw("hist err");
    ratios_5[i] -> Write();
  }
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
  ratios_0[0] -> GetYaxis()->SetTitleOffset(1.2);
  ratios_0[0] -> Draw("hist");

  TPad *pad_sector_1 = (TPad*)pad_sector->cd(2);
  pad_sector_1 -> cd();
  pad_sector_1 -> SetBottomMargin(0.15);
  pad_sector_1 -> SetLeftMargin(0.15);
  ratios_1[0] -> GetYaxis()->SetTitleOffset(1.2);
  ratios_1[0] -> Draw("hist");

  TPad *pad_sector_2 = (TPad*)pad_sector->cd(3);
  pad_sector_2 -> cd();
  pad_sector_2 -> SetBottomMargin(0.15);
  pad_sector_2 -> SetLeftMargin(0.15);
  ratios_2[0] -> GetYaxis()->SetTitleOffset(1.2);
  ratios_2[0] -> Draw("hist");

  TPad *pad_sector_3 = (TPad*)pad_sector->cd(4);
  pad_sector_3 -> cd();
  pad_sector_3 -> SetBottomMargin(0.15);
  pad_sector_3 -> SetLeftMargin(0.15);
  ratios_3[0] -> GetYaxis()->SetTitleOffset(1.2);
  ratios_3[0] -> Draw("hist");

  TPad *pad_sector_4 = (TPad*)pad_sector->cd(5);
  pad_sector_4 -> cd();
  pad_sector_4 -> SetBottomMargin(0.15);
  pad_sector_4 -> SetLeftMargin(0.15);
  ratios_4[0] -> GetYaxis()->SetTitleOffset(1.2);
  ratios_4[0] -> Draw("hist");

  TPad *pad_sector_5 = (TPad*)pad_sector->cd(6);
  pad_sector_5 -> cd();
  pad_sector_5 -> SetBottomMargin(0.15);
  pad_sector_5 -> SetLeftMargin(0.15);
  ratios_5[0] -> GetYaxis()->SetTitleOffset(1.2);
  ratios_5[0] -> Draw("hist");

  c_sector_2->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_persector.root").c_str());
  c_sector_2->SaveAs((output_location+"/AcceptanceFiles/"+output_name+"_persector.pdf").c_str());
  delete c_sector_2 ;

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
// id_secotr: ID of the sector of interest

void Plot1DXSec(std::vector<std::string> MC_files_name, std::string data_file_name,
                std::string acceptance_file_name, std::string observable,
                std::string title, std::string data_name, std::vector<std::string> model,
                std::string input_MC_location, std::string input_data_location, std::string output_location,
                int id_sector = -1 /*all*/ ) {

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
  TTree * tree_true = (TTree*)files_true_MC[0]->Get("MCCLAS6Tree");

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
  hist_true_0 -> SetName( "MC True Sector 0") ;
  hist_true_0 -> Reset("ICE");
  TH1D * hist_true_1 = (TH1D*) h_acceptance_1 ->Clone();
  hist_true_1 -> SetName( "MC True Sector 1") ;
  hist_true_1 -> Reset("ICE");
  TH1D * hist_true_2 = (TH1D*) h_acceptance_2 ->Clone();
  hist_true_2 -> SetName( "MC True Sector 2") ;
  hist_true_2 -> Reset("ICE");
  TH1D * hist_true_3 = (TH1D*) h_acceptance_3 ->Clone();
  hist_true_3 -> SetName( "MC True Sector 3") ;
  hist_true_3 -> Reset("ICE");
  TH1D * hist_true_4 = (TH1D*) h_acceptance_4 ->Clone();
  hist_true_4 -> SetName( "MC True Sector 4") ;
  hist_true_4 -> Reset("ICE");
  TH1D * hist_true_5 = (TH1D*) h_acceptance_5 ->Clone();
  hist_true_5 -> SetName( "MC True Sector 5") ;
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

  // total and per sector
  TH1D * hist_data = (TH1D*) h_acceptance ->Clone();
  hist_data -> SetName( "Data") ;
  hist_data -> Reset("ICE");
  TH1D * hist_data_0 = (TH1D*) h_acceptance_0 ->Clone();
  hist_data_0 -> SetName( "Data Sector 0") ;
  hist_data_0 -> Reset("ICE");
  TH1D * hist_data_1 = (TH1D*) h_acceptance_1 ->Clone();
  hist_data_1 -> SetName( "Data Sector 1") ;
  hist_data_1 -> Reset("ICE");
  TH1D * hist_data_2 = (TH1D*) h_acceptance_2 ->Clone();
  hist_data_2 -> SetName( "Data Sector 2") ;
  hist_data_2 -> Reset("ICE");
  TH1D * hist_data_3 = (TH1D*) h_acceptance_3 ->Clone();
  hist_data_3 -> SetName( "Data Sector 3") ;
  hist_data_3 -> Reset("ICE");
  TH1D * hist_data_4 = (TH1D*) h_acceptance_4 ->Clone();
  hist_data_4 -> SetName( "Data Sector 4") ;
  hist_data_4 -> Reset("ICE");
  TH1D * hist_data_5 = (TH1D*) h_acceptance_5 ->Clone();
  hist_data_5 -> SetName( "Data Sector 5") ;
  hist_data_5 -> Reset("ICE");

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
  double RecoXBJK, RecoEnergyTransfer, RecoQ2 ;
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
    trees[i] -> SetBranchAddress("QEL",&QEL);
    trees[i] -> SetBranchAddress("RES",&RES);
    trees[i] -> SetBranchAddress("MEC",&MEC);
    trees[i] -> SetBranchAddress("DIS",&DIS);

    if( i == 0 ) trees[i] -> SetBranchAddress("MCNormalization", &MCNormalization );
    else if ( i == 1 ) trees[i] -> SetBranchAddress("DataNormalization",&DataNormalization );

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
      else if ( observable == "RecoXBJK") content = RecoXBJK ;
      else if ( observable == "RecoQ2") content = RecoQ2 ;
      else if ( observable == "RecoEnergyTransfer") content = RecoEnergyTransfer ;
      else if ( observable == "AlphaT") content = AlphaT ;
      else if ( observable == "HadAlphaT") content = HadAlphaT ;
      else if ( observable == "DeltaPT") content = DeltaPT ;
      else if ( observable == "HadDeltaPT") content = HadDeltaPT ;
      else if ( observable == "DeltaPhiT") content = DeltaPhiT ;
      else if ( observable == "HadDeltaPhiT") content = HadDeltaPhiT ;

      unsigned int id_hist = i ;
      // Fill the per sector histogram. Only for primary model
      if( i != 2 ) hists[size_primary_trees*(ElectronSector+1)+i] -> Fill( content, w ) ;
      if( i == 2 ) id_hist = size_primary_hists ;

      if( id_sector > 0 ) {
      	// Compute only for sector of interest
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
  hist_data_uncorr_0 -> SetName( "Uncorrected Data Sector 0") ;
  TH1D * hist_data_uncorr_1 = (TH1D*) hist_data_1 ->Clone();
  hist_data_uncorr_1 -> SetName( "Uncorrected Data Sector 1") ;
  TH1D * hist_data_uncorr_2 = (TH1D*) hist_data_2 ->Clone();
  hist_data_uncorr_2 -> SetName( "Uncorrected Data Sector 2") ;
  TH1D * hist_data_uncorr_3 = (TH1D*) hist_data_3 ->Clone();
  hist_data_uncorr_3 -> SetName( "Uncorrected Data Sector 3") ;
  TH1D * hist_data_uncorr_4 = (TH1D*) hist_data_4 ->Clone();
  hist_data_uncorr_4 -> SetName( "Uncorrected Data Sector 4") ;
  TH1D * hist_data_uncorr_5 = (TH1D*) hist_data_5 ->Clone();
  hist_data_uncorr_5 -> SetName( "Uncorrected Data Sector 5") ;

  // Correct data for detector acceptance :
  CorrectData(hist_data, h_acceptance);
  CorrectData(hist_data_0, h_acceptance_0);
  CorrectData(hist_data_1, h_acceptance_1);
  CorrectData(hist_data_2, h_acceptance_2);
  CorrectData(hist_data_3, h_acceptance_3);
  CorrectData(hist_data_4, h_acceptance_4);
  CorrectData(hist_data_5, h_acceptance_5);

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

  // Format plots
  StandardFormat( hist_data, title, kBlack, 8, observable ) ;
  StandardFormat( hist_data_0, title+" Sector 0", kOrange+1, 8, observable ) ;
  StandardFormat( hist_data_1, title+" Sector 1", kPink+4, 8, observable ) ;
  StandardFormat( hist_data_2, title+" Sector 2", kViolet+5, 8, observable ) ;
  StandardFormat( hist_data_3, title+" Sector 3", kAzure-5, 8, observable ) ;
  StandardFormat( hist_data_4, title+" Sector 4", kTeal-7, 8, observable ) ;
  StandardFormat( hist_data_5, title+" Sector 5", kGreen-3, 8, observable ) ;
  hist_data -> SetLineWidth(0);

  StandardFormat( hist_data_uncorr, title, kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_0, title+" Sector 0", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_1, title+" Sector 1", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_2, title+" Sector 2", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_3, title+" Sector 3", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_4, title+" Sector 4", kRed, 8, observable ) ;
  StandardFormat( hist_data_uncorr_5, title+" Sector 5", kRed, 8, observable ) ;

  StandardFormat( hist_true, title, kBlack, 1, observable ) ;
  StandardFormat( hist_true_0, title+" Sector 0", kBlack, 1, observable ) ;
  StandardFormat( hist_true_1, title+" Sector 1", kBlack, 1, observable ) ;
  StandardFormat( hist_true_2, title+" Sector 2", kBlack, 1, observable ) ;
  StandardFormat( hist_true_3, title+" Sector 3", kBlack, 1, observable ) ;
  StandardFormat( hist_true_4, title+" Sector 4", kBlack, 1, observable ) ;
  StandardFormat( hist_true_5, title+" Sector 5", kBlack, 1, observable ) ;

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

  std::string output_name = MC_files_name[0]+"_dxsec_d"+observable ;
  if( id_sector > 0 ) output_name += "_sector_"+std::to_string(id_sector) ;
  c1->SaveAs((output_location+"/TotalXSec/"+output_name+".root").c_str());
  c1->SaveAs((output_location+"/TotalXSec/"+output_name+".pdf").c_str());
  delete c1 ;

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
     leg->AddEntry(hists_true_submodel[id],("GENIE "+model[id]).c_str(),"l");
  }
  leg->AddEntry(hist_data, data_name.c_str(), "lp");
  leg->Draw();
  output_name = MC_files_name[0] ;
  c_leg->SaveAs((output_location+output_name+"_legend.root").c_str());
  c_leg->SaveAs((output_location+output_name+"_legend.pdf").c_str());

  delete c_leg;
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
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, EBeam+0.2 );
  }	else if ( observable == "pfl_theta") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 10, 20, 70 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 20, 60 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 15, 50 );
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
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 100 );
  } else if ( observable == "proton_phi") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 180 );
  } else if ( observable == "pim_mom") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0.1, 0.7 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 1.5 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 0, 1.6 );
  }	else if ( observable == "pim_theta") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 15, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 30, 150 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 30, 150 );
  }	else if ( observable == "pip_mom") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
  }	else if ( observable == "pip_theta") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 180 );
  }	else if ( observable == "RecoW") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 16, 1, 1.5 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 18, 1, 1.9 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 10, 1.2, 2.5 );
  }	else if ( observable == "RecoXBJK") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
  } else if ( observable == "RecoQ2") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
  }	else if ( observable == "Recoq3") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
  }	else if ( observable == "DeltaPT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
  }	else if ( observable == "HadDeltaPT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, EBeam+0.2 );
  }	else if ( observable == "DeltaPhiT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0,180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 180 );
  } else if ( observable == "HadDeltaPhiT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 180 );
  }	else if ( observable == "AlphaT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 180 );
  } else if ( observable == "HadAlphaT") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 15, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 15, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 15, 0, 180 );
  }	else if ( observable == "RecoEnergyTransfer") {
    if( EBeam == 1.161 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 2.261 ) binning = GetUniformBinning( 20, 0, 180 );
    else if( EBeam == 4.461 ) binning = GetUniformBinning( 20, 0, 180 );
  }

	return binning ;
}

void Plot1DXSec(){

	// Add without the .root
  std::string mc_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/mc_files/";
  std::string data_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/data_files/";
  std::string output_location = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/PionAnalysis/1p1pi/output_files/";
  std::vector<std::string> mc_files = {"e4nuanalysis_1p1pimanalysis_G18_10a_Q2_08_e_on_1000060120_4461MeV_NoRad"};
  std::vector<std::string> model_names = { "G18_10a" } ;
  if( model_names.size() != mc_files.size() ) { std::cout<< " Need same number of models and names" << std::endl; return ; }

  std::string file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_4461MeV_clas6data";
  std::string title = "e^{12}C 1p1#pi^{-} at 4.416 GeV";
  std::string data_name = "CLAS6 data";

  std::vector<std::string> observables = {"RecoW","pfl","pfl_theta","pim_mom","pim_theta","proton_mom","proton_theta","HadAlphaT","HadDeltaPT","HadDeltaPhiT"};
  for( unsigned int i = 0 ; i < observables.size(); ++i ){
    mc_files = {"e4nuanalysis_1p1pimanalysis_G18_10a_Q2_08_e_on_1000060120_4461MeV_NoRad"};
    file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_4461MeV_clas6data";
    title = "e^{12}C 1p1#pi^{-} at 4.416 GeV";
    std::vector<double> binning = GetBinning(observables[i],4.461);
	  std::string acceptance_file = compute_acceptance( mc_files, observables[i], title, binning, mc_location, output_location ) ;
    Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, -1 ) ;

    mc_files = {"e4nuanalysis_1p1pimanalysis_G18_10a_Q2_04_e_on_1000060120_2261MeV_NoRad"};
    file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_2261MeV_clas6data";
    title = "e^{12}C 1p1#pi^{-} at 2.216 GeV";

    binning = GetBinning(observables[i],2.261);
	  acceptance_file = compute_acceptance( mc_files, observables[i], title, binning, mc_location, output_location ) ;
    Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, -1 ) ;

    mc_files = {"e4nuanalysis_1p1pimanalysis_G18_10a_Q2_01_e_on_1000060120_1161MeV_NoRad"};
    file_data = "e4nuanalysis_1p1pimanalysis_e_on_1000060120_1161MeV_clas6data";
    title = "e^{12}C 1p1#pi^{-} at 1.116 GeV";

    binning = GetBinning(observables[i],1.161);
    acceptance_file = compute_acceptance( mc_files, observables[i], title, binning, mc_location, output_location ) ;
    Plot1DXSec( mc_files, file_data, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, -1 ) ;
  }
}
