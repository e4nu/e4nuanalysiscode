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

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf;
using namespace e4nu::plotting;
/////////////////////////////////////////////////////////////////
// Options:                                                    //
// * input-files : comma sepparated list of root file          //
//                 the first histogram is the reference one    //
// * output-file : file to store plots (in root, pdf... format)//
// * observable : observable used for the x axis definition    //
// * legend-list : tags for the legend                         //
// * color-list : list of colors for each hist                 //
// * plot-diff: plots difference with default                  //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting histograms..." << std::endl;

  if ( argc == 0 ) {
    std::cout << "Input files not specified. Abort..." << std::endl;
    return 0;
  }

  // Color palette
  std::vector<int> color_list = {kBlack,kBlue-4,kOrange+1,kGreen+3,kMagenta+2,kRed-4,kSpring+7,kTeal-6,kCyan+3,kAzure-1,kBlue-8,kViolet-4,kMagenta-10,kPink-9} ;

  std::vector<string> input_files, legend_list ;
  std::vector<int> colors;
  std::string output_file = "histograms";
  std::string observable = "ECal";
  bool plot_diff = false;
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
    } else { 
      for( int i = 0 ; i < input_files.size() ; ++i ) legend_list.push_back( "Model "+to_string(i));
    }

    if( ExistArg("output-file",argc,argv)) {
      output_file = GetArg("output-file",argc,argv);
    }
    if( ExistArg("observable",argc,argv)) {
      observable = GetArg("observable",argc,argv);
    }
    if( ExistArg("color-list",argc,argv)) {
      string input = GetArg("color-list",argc,argv);
      stringstream ss(input);
      while( ss.good() )
	{
	  string substr;
	  getline( ss, substr, ',' );
	  colors.push_back( stoi(substr) );
	}
    }
    if( ExistArg("plot-diff",argc,argv)) plot_diff = true ;

  }

  TCanvas * c  = new TCanvas("","",800,800);
  c->SetFillColor(0);
  TPad *pad11 = new TPad("pad1","",0,0,1,1);
  pad11->Draw();
  pad11->cd();
  pad11->SetBottomMargin(0.15);
  pad11->SetLeftMargin(0.15);

  // Top canvas - distribution
  // Bottom canvas - relative difference
  TPad *pad1 = new TPad("pad1","",0,0.4,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,0.4);

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


  auto legend = new TLegend(0,0,1,1);
  legend->SetLineWidth(0);
  vector<TFile*> in_root_files ;
  vector<TH1D*> hists ; 
  vector<TH1D*> hists_diff ;
  TFile * def = new TFile( input_files[0].c_str(), "ROOTFile" );
  TH1D* hist_def = (TH1D*) def -> Get( observable.c_str() ) ; 
  if(! hist_def ) return 0;

  for( unsigned int i = 0 ; i < input_files.size() ; ++i ) { 
    in_root_files.push_back(new TFile( input_files[i].c_str(), "ROOTFile" )) ;
    if( !in_root_files[i] ) {
      std::cout << " root file " << input_files[i] << "does not exist" << std::endl;
      return 0;
    }
    std::cout <<input_files[i]<<std::endl;
    hists.push_back( (TH1D*) in_root_files[i]->Get( observable.c_str() ) ) ;
    if( !hists[i] ) return 0;
    if( plot_diff ) { 
	hists_diff.push_back( (TH1D*) hist_def->Clone() ) ;
	hists_diff[i]->Scale(-1.);
    	hists_diff[i]->Add(hists[i]);
    	hists_diff[i]->Divide(hists[i]);
    	hists_diff[i]->Scale(100.); // Relative error in %
    }

    int color = color_list[i] ;
    if( colors.size() == input_files.size() ) color = colors[i];
  
    if(plot_diff){
    	hists_diff[i] -> SetLineColor(color);
    	hists_diff[i] -> SetMarkerStyle(8);
    	hists_diff[i] -> SetMarkerColor(color);
    	hists_diff[i] -> SetLineWidth(2);
    	hists_diff[i]->GetYaxis()->SetLabelSize(0.08);
    	hists_diff[i]->GetXaxis()->SetLabelSize(0.08);
    	hists_diff[i]->GetXaxis()->SetTitleSize(0.08);
    	hists_diff[i]->GetYaxis()->SetTitleSize(0.08);
    	hists_diff[i] -> SetTitle("");
    	hists_diff[i] -> GetXaxis()->SetTitle(plotting::GetAxisLabel(observable,0).c_str());
    	hists_diff[i] -> GetYaxis()->SetTitle(plotting::GetAxisLabel(observable,1).c_str());
    	hists_diff[i] -> GetXaxis()->CenterTitle();
    	hists_diff[i] -> GetYaxis()->CenterTitle();
    	hists_diff[i] -> GetYaxis()->SetTitle("Rel.Diff");
    	hists_diff[i]->GetYaxis()->SetRangeUser(-50,50);
    }    
  }

  // Setting formatt
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
  
  if( plot_diff ) { 
    pad1->Draw();
    pad1->cd();
    pad1->SetBottomMargin(0.015);
    pad1->SetLeftMargin(0.1);
  }

 double integral_0 = hists[0]->Integral() ;

 double y_max =0; 
 for( unsigned int i = 0 ; i < input_files.size() ; ++i ) {   
    double integral_i = hists[i]->Integral();
    double scaling = 1; 
    if(integral_i!=0) scaling = round(integral_0/integral_i);
    if( scaling != 1 ) legend_list[i] += " (x"+to_string((int)scaling)+")" ;
    legend->AddEntry(hists[i], legend_list[i].c_str());  
    hists[i]->Scale(scaling);

    y_max = plotting::GetMaximum( hists ) ; 
  }

  for( unsigned int i = 0 ; i < input_files.size() ; ++i ) {
    hists[i]->GetYaxis()->SetRangeUser(0,y_max);
    hists[i]->GetYaxis()->SetMaxDigits(3) ;   
    hists[i]->SetLineStyle(1);
    int color = color_list[i];
    if( colors.size() == input_files.size() ) color = colors[i];
    hists[i]->SetLineColor(color);
    hists[i]->SetMarkerColor(color);
    hists[i]->SetTitle("");

    if( i == 0 ) hists[i] -> Draw("ERR");
    else  hists[i] -> Draw("ERR same");
  }

  // Top canvas - distribution
  // Bottom canvas - relative difference
  if( plot_diff ) { 
    c->cd();
    pad2->Draw();
    pad2->cd();
    pad2->SetBottomMargin(0.25);
    pad2->SetLeftMargin(0.1);
    pad2->SetTopMargin(0.05);
    for( unsigned int i = 0 ; i < input_files.size() ; ++i ) {   
      if( i == 0 ) hists_diff[i]->Draw("hist ERR");
      else  hists_diff[i]->Draw("hist ERR same");
    }
  }
  c->SaveAs((output_file+".root").c_str());
  c->SaveAs((output_file+".pdf").c_str());

  c->Clear();
  if( legend_list.size() != 0) legend->Draw();

  c->SaveAs((output_file+"_legend.root").c_str());
  c->SaveAs((output_file+"_legend.pdf").c_str());
  return 0 ;

}
