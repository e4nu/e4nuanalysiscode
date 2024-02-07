// __________________________________________________________________________
/* This app is used to plot the estimated vs true background obtained from  */
/* The output root file computed with bkg debug mode ("DebugBkg true")      */
/* An example of how to run the analysis with this mode is available in     */
/* ConfFiles/example_configuration.txt                                      */
// __________________________________________________________________________
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "utils/RadiativeCorrUtils.h"
#include "plotting/PlottingUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TGaxis.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf;
using namespace e4nu::plotting;
/////////////////////////////////////////////////////////////////
// Options:                                                    //
// * input-files : comma sepparated list of data root files    //
// * output-file : file to store plots (in root, pdf... format)//
// * observable : observable used for the x axis definition    //
// * legend-list : tags for the legend                         //
// * analysis-key: i.e. 1p1pim                                 //
// * sector : default all                                      //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting histograms..." << std::endl;

  if ( argc == 0 ) {
    std::cout << "Input files not specified. Abort..." << std::endl;
    return 0;
  }

  std::vector<string> input_files, legend_list ;
  std::string output_file = "histogram_data";
  std::string observable = "ECal";
  int sector = -9999 ; // all
  std::string analysis_key = "1p1pim";
  double EBeam = 4.461 ;
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-files",argc,argv)) {
      string input = GetArg("input-files",argc,argv);
      stringstream ss(input);
      while( ss.good() )
	{
	  string substr;
	  getline( ss, substr, ',' );
	  input_files.push_back( substr );
	}
      if( input_files.size() == 0 ) return 0;
    } else { return 0 ;}

    if( ExistArg("legend-list",argc,argv)) {
      string input = GetArg("legend-list",argc,argv);
      stringstream ss(input);
      while( ss.good() )
	{
	  string substr;
	  getline( ss, substr, ',' );
	  legend_list.push_back( substr );
	}
      if( input_files.size() == 0 ) return 0;
    }

    if( ExistArg("output-file",argc,argv)) {
      output_file = GetArg("output-file",argc,argv);
    }
    if( ExistArg("observable",argc,argv)) {
      observable = GetArg("observable",argc,argv);
    }
    if( ExistArg("beam-energy",argc,argv)) {
      EBeam = stof(GetArg("beam-energy",argc,argv));
      std::cout << EBeam << std::endl;
    }

  }

  // Color palette
  std::vector<int> color_list = {kBlack,kBlue-4,kOrange+1,kGreen+3,kMagenta+2,kRed-4,kSpring+7,kTeal-6,kCyan+3,kAzure-1,kBlue-8,kViolet-4,kMagenta-10,kPink-9} ;

  TCanvas * c  = new TCanvas("","",800,800);
  c->SetFillColor(0);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);

  auto legend = new TLegend(0.2,0.65,0.55,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(132);
  legend->SetTextSize(0.08);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03);

  std::vector<TFile*> in_root_files ;
  std::vector<TTree*> in_trees;
  std::vector<TH1D*> hists;
  for( unsigned int id = 0 ; id < input_files.size(); ++id ){
    in_root_files.push_back(new TFile((input_files[id]).c_str(),"ROOT"));
    if( !in_root_files[id] ) { std::cout << "ERROR: the "<< input_files[id]<<" does not exist." << std::endl; return 0;}
  }

  for( unsigned int id = 0 ; id < input_files.size() ; ++id ){
    in_trees.push_back( (TTree*)in_root_files[id]->Get("CLAS6Tree") );
    if( !in_trees[id] ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return 0;}
  }

  std::vector<double> binning = plotting::GetBinning( observable, EBeam, analysis_key );
  for( unsigned int i = 0 ; i < input_files.size(); ++i ){
    hists.push_back( new TH1D( ("Dataset "+std::to_string(i)).c_str(), "", binning.size()-1, &binning[0] ) ) ;
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
  double HadronsAngle,Angleqvshad ;
  double AdlerAngleThetaP, AdlerAnglePhiP, AdlerAngleThetaPi, AdlerAnglePhiPi ;
  long NEntries ;
  bool IsBkg ;
  int ElectronSector ;

  for ( unsigned int i = 0 ; i < in_trees.size() ; ++i ){
    NEntries = in_trees[i] -> GetEntries() ;
    if( !in_trees[i] ) continue ;
    in_trees[i] -> SetBranchAddress("TotWeight",&TotWeight);
    in_trees[i] -> SetBranchAddress("IsBkg",&IsBkg);
    in_trees[i] -> SetBranchAddress("ECal",&ECal);
    in_trees[i] -> SetBranchAddress("pfl_theta",&pfl_theta);
    in_trees[i] -> SetBranchAddress("pfl_phi",&pfl_phi);
    in_trees[i] -> SetBranchAddress("pfl",&pfl);
    in_trees[i] -> SetBranchAddress("proton_mom",&proton_mom);
    in_trees[i] -> SetBranchAddress("proton_theta",&proton_theta);
    in_trees[i] -> SetBranchAddress("proton_phi",&proton_phi);
    in_trees[i] -> SetBranchAddress("pim_mom",&pim_mom);
    in_trees[i] -> SetBranchAddress("pim_theta",&pim_theta);
    in_trees[i] -> SetBranchAddress("pim_phi",&pim_phi);
    in_trees[i] -> SetBranchAddress("pip_mom",&pip_mom);
    in_trees[i] -> SetBranchAddress("pip_theta",&pip_theta);
    in_trees[i] -> SetBranchAddress("pip_phi",&pip_phi);
    in_trees[i] -> SetBranchAddress("RecoW",&RecoW);
    in_trees[i] -> SetBranchAddress("Recoq3",&Recoq3);
    in_trees[i] -> SetBranchAddress("RecoQELEnu",&RecoQELEnu);
    in_trees[i] -> SetBranchAddress("RecoXBJK",&RecoXBJK);
    in_trees[i] -> SetBranchAddress("RecoQ2",&RecoQ2);
    in_trees[i] -> SetBranchAddress("RecoEnergyTransfer",&RecoEnergyTransfer);
    in_trees[i] -> SetBranchAddress("AlphaT",&AlphaT);
    in_trees[i] -> SetBranchAddress("HadAlphaT",&HadAlphaT);
    in_trees[i] -> SetBranchAddress("DeltaPT",&DeltaPT);
    in_trees[i] -> SetBranchAddress("HadDeltaPT",&HadDeltaPT);
    in_trees[i] -> SetBranchAddress("HadDeltaPTx",&HadDeltaPTx);
    in_trees[i] -> SetBranchAddress("HadDeltaPTy",&HadDeltaPTy);
    in_trees[i] -> SetBranchAddress("DeltaPhiT",&DeltaPhiT);
    in_trees[i] -> SetBranchAddress("HadDeltaPhiT",&HadDeltaPhiT);
    in_trees[i] -> SetBranchAddress("ElectronSector",&ElectronSector);
    in_trees[i] -> SetBranchAddress("HadSystemMass", &HadSystemMass);
    in_trees[i] -> SetBranchAddress("MissingEnergy", &MissingEnergy);
    in_trees[i] -> SetBranchAddress("MissingAngle", &MissingAngle);
    in_trees[i] -> SetBranchAddress("MissingMomentum", &MissingMomentum);
    in_trees[i] -> SetBranchAddress("InferedNucleonMom", &InferedNucleonMom);
    in_trees[i] -> SetBranchAddress("HadronsAngle", &HadronsAngle);
    in_trees[i] -> SetBranchAddress("AdlerAngleThetaP", &AdlerAngleThetaP);
    in_trees[i] -> SetBranchAddress("AdlerAnglePhiP", &AdlerAnglePhiP);
    in_trees[i] -> SetBranchAddress("AdlerAngleThetaPi", &AdlerAngleThetaPi);
    in_trees[i] -> SetBranchAddress("AdlerAnglePhiPi", &AdlerAnglePhiPi);
    in_trees[i] -> SetBranchAddress("Angleqvshad",&Angleqvshad);

    for( int j = 0 ; j < NEntries ; ++j ) {
      in_trees[i]->GetEntry(j) ;
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
      hists[i]->Fill(content,w);
    }
  }

  double ymax = 0 ;
  for( unsigned int i = 0 ; i < input_files.size(); ++i ){
    double max = 0 ;
    for( unsigned int k = 0 ; k < hists[i]->GetNbinsX() ; ++k ) {
      if ( hists[i]->GetBinContent(k) > max ) max = hists[i]->GetBinContent(k) ;
    }
    if( max > ymax ) ymax = max;
  }
  ymax *= ( 1+0.2 ) ;

  // Normalise for bin size
  for( unsigned int i = 0 ; i < input_files.size(); ++i ){
    plotting::StandardFormat( hists[i], "", color_list[i], 1, observable, ymax, "Counts/Bin Width");
    hists[i] -> SetLineStyle(1);
    hists[i] -> SetMarkerStyle(8);
    hists[i] -> GetYaxis() -> TAxis::SetMaxDigits(3);
    plotting::NormalizeHist(hists[i], 1 );
    if( i == 0 ) hists[i] -> Draw("err");
    else hists[i] -> Draw("err same");

    if ( legend_list.size() == input_files.size()){
      legend->AddEntry(hists[i],(legend_list[i]).c_str());
    }
  }
  hists[0] -> Draw("err same");
  if ( legend_list.size() == input_files.size()) legend->Draw();

  std::string output_name = output_file+"_Nevents_"+observable ;

  c->SaveAs((output_file+".root").c_str());
  c->SaveAs((output_file+".pdf").c_str());

  return 0 ;

}
