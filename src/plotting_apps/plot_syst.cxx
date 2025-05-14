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
using namespace e4nu::conf;
using namespace e4nu::plotting;
/////////////////////////////////////////////////////////////////
// Options:                                                    //
// * input-files : comma sepparated list of data root files    //
// * output-file : file to store plots (in root, pdf... format)//
// * observable : observable used for the x axis definition    //
// * legend-list : tags for the legend                         //
// * beam-energy : energy                                      //
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
  std::string output_file = "histogram_syst";
  std::string observable = "ECal";
  int sector = -9999 ; // all
  std::string analysis_key = "1p1pim";
  double EBeam = 1 ;
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
    double energy = 0 ;
    if( ExistArg("beam-energy",argc,argv)) {
      energy = stof(GetArg("beam-energy",argc,argv));
    }

    if( energy > 0  && energy < 1.5 ) EBeam = 1.161 ;
    else if( energy >1.5 && energy < 3 ) EBeam = 2.261 ;
    else if( energy >3 && energy < 5 ) EBeam = 4.461 ;
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
  long NEntries ;

  for ( unsigned int i = 0 ; i < in_trees.size() ; ++i ){
    NEntries = in_trees[i] -> GetEntries() ;
    if( !in_trees[i] ) continue ;
    plotting::SetAnalysisBranch( in_trees[i] ) ;

    for( int j = 0 ; j < NEntries ; ++j ) {
      in_trees[i]->GetEntry(j) ;
      double content = content = GetObservable(observable);
      double w = EventWght * AccWght ;
      //if( scale_mott ) w *= MottXSecScale;
      hists[i]->Fill(content,w);
    }
  }

  // Store systematics in histogram
  std::vector<TH1D*> hist_syst ;

  for( unsigned int i = 0 ; i < hists.size() - 1; ++i ) {
    hist_syst.push_back( new TH1D( ("Systematics"+to_string(i)).c_str(), "", binning.size()-1, &binning[0] ));
    for( unsigned int j = 1 ; j < binning.size()-1 ; ++j ) {
      double err = pow(hists[i]->GetBinContent(j)-hists[i+1]->GetBinContent(j),2) ;
      if( err < 0 ) err = 0;
      if( hists[i]->GetBinContent(j) != 0 ) hist_syst[i] -> SetBinContent( j, sqrt(err)/hists[i]->GetBinContent(j)*100 );
      else hist_syst[i] -> SetBinContent( j, 0 );
    }
  }

  for( unsigned int i = 0 ; i < hists.size() -1 ; ++i ) {
    hist_syst[i]->SetLineColor(kBlue);
    hist_syst[i]->SetMarkerColor(kBlue);
    hist_syst[i]->SetMarkerStyle(8);
    hist_syst[i]->GetXaxis()->SetTitle(observable.c_str());
    hist_syst[i]->GetYaxis()->SetTitle("#sigma_{RMS}/x_{Nom}[%]");

    if( i == 0 ) hist_syst[i]->Draw();
    else hist_syst[i]->Draw("same");
  }
  std::string output_name = output_file+"_Nevents_"+observable ;

  c->SaveAs((output_file+".root").c_str());
  c->SaveAs((output_file+".pdf").c_str());

  return 0 ;

}
