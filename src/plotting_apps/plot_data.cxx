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
// * analysis-key: i.e. 1p1pim                                 //
// * sector : default all                                      //
// * beam-energy : energy                                      //
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
    if( ExistArg("analysis-key",argc,argv)) {
      analysis_key = GetArg("analysis-key",argc,argv);
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
  for ( unsigned int i = 0 ; i < in_trees.size() ; ++i ){
    NEntries = in_trees[i] -> GetEntries() ;
    if( !in_trees[i] ) continue ;
    plotting::SetAnalysisBranch( in_trees[i] ) ;

    for( int j = 0 ; j < NEntries ; ++j ) {
      in_trees[i]->GetEntry(j) ;
      double content = GetObservable(observable);
      double w = EventWght * AccWght ;
      //if( scale_mott )
      //w *= MottXSecScale;
      hists[i]->Fill(content,w);
    }
  }

  double ymax = 0 ;
  for( unsigned int i = 0 ; i < input_files.size(); ++i ){
    double max = 0 ;
    for( int k = 0 ; k < hists[i]->GetNbinsX() ; ++k ) {
      if ( hists[i]->GetBinContent(k) > max ) max = hists[i]->GetBinContent(k) ;
    }
    if( max > ymax ) ymax = max;
  }
  ymax *= ( 1+0.2 ) ;

  // Normalise for bin size
  for( unsigned int i = 0 ; i < input_files.size(); ++i ){
    plotting::StandardFormat( hists[i], "", color_list[i], 1, observable, false, ymax, "Counts/Bin Width");
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
