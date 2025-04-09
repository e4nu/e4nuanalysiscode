// __________________________________________________________________________
/* This app is used to plot the estimated vs true background obtained from  */
/* The output root file computed with bkg debug mode ("DebugBkg true")      */
/* An example of how to run the analysis with this mode is available in     */
/* ConfFiles/example_configuration.txt                                      */
// __________________________________________________________________________
#include <iostream>
#include <vector>
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
// * input-file : root file containing histograms              //
// * output-file : file to store plots (in root, pdf... format)//
// * observable : observable used for the x axis definition    //
// * bkg-string : defaulted to _TotTrueBkg                     //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting e4nu analysis background check..." << std::endl;

  if ( argc == 0 ) {
    std::cout << "Input file not specified. Abort..." << std::endl;
    return 0;
  }

  std::string input_file ;
  std::string output_file = "bakground_debug.root";
  std::string observable = "ECal";
  std::string bkg_string = "";
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-file",argc,argv)) {
      input_file = GetArg("input-file",argc,argv);
      std::cout << " Using " << input_file << " for plotting "<<std::endl;
    } else { return 0 ;}

    if( ExistArg("output-file",argc,argv)) {
      output_file = GetArg("output-file",argc,argv);
    }
    if( ExistArg("observable",argc,argv)) {
      observable = GetArg("observable",argc,argv);
    }
    if( ExistArg("bkg-string",argc,argv)) {
      bkg_string = GetArg("bkg-string",argc,argv);
    }
  }

  TFile* in_root_file = new TFile( input_file.c_str(), "ROOTFile" ) ;
  if( !in_root_file ) {
    std::cout << " root file " << input_file << "does not exist" << std::endl;
    return 0;
  }

  TCanvas * c  = new TCanvas("","",800,800);
  c->SetFillColor(0);
  // Top canvas - distribution
  // Bottom canvas - relative difference
  TPad *pad1 = new TPad("pad1","",0,0.4,1,1);
  TPad *pad2 = new TPad("pad2","",0,0,1,0.4);

  TH1D * h_tot_true = (TH1D*) in_root_file->Get( (observable+"_TotTrueBkg"+bkg_string).c_str() ) ;
  TH1D * h_tot_est = (TH1D*) in_root_file->Get( (observable+"_TotEstBkg"+bkg_string).c_str() ) ;

  TH1D * h_signal = (TH1D*) in_root_file->Get( (observable+"_OnlySignal").c_str() ) ;
  TH1D * h_total_bkg = (TH1D*) in_root_file->Get( (observable+"_TotTrueBkg").c_str() ) ;

  TH1D * h_diff_true = (TH1D*) h_tot_true->Clone();
  double fraction = h_total_bkg->GetEntries()/h_signal->GetEntries();
  std::cout << " Fraction Background Events: " << fraction << std::endl;

  // Scale by number of signal Events
  // h_tot_true->Scale(1./h_signal->GetEntries());
  // h_tot_est->Scale(1./h_signal->GetEntries());

  // Scale by Bin Width
  // for (int k = 1; k <= h_tot_true->GetNbinsX(); k++)
  // {
  //   double width = h_tot_true->GetBinWidth(k);
  //   // h_tot_true->SetBinContent(k, h_tot_true->GetBinContent(k) / width );
  //   // h_tot_true->SetBinError(k, h_tot_true->GetBinError(k) / width );
  //   // h_tot_est->SetBinContent(k, h_tot_est->GetBinContent(k) / width / fraction );
  //   // h_tot_est->SetBinError(k, h_tot_est->GetBinError(k) / width / fraction );
  // }

  for( unsigned int j = 1 ; j < h_diff_true -> GetNbinsX() ; ++j ){
    double err = pow(h_tot_true->GetBinContent(j)-h_tot_est->GetBinContent(j),2) ;

    // Substract stat error
    err -= pow(h_tot_true->GetBinError(j),2) ;
    err -= pow(h_tot_est->GetBinError(j),2) ;

    if( err < 0 ) err = 0 ;
    
    err = sqrt(err);
    if( h_tot_est->GetBinContent(j) != 0 ) {
      h_diff_true->SetBinContent(j,err/h_tot_est->GetBinContent(j)*100*fraction);
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

  h_tot_true -> SetLineColor(kRed);
  h_tot_true -> SetMarkerStyle(8);
  h_tot_true -> SetMarkerColor(kRed);
  h_tot_true -> SetLineWidth(2);
  h_tot_true -> GetYaxis()->SetLabelSize(0.06);
  h_tot_true -> GetYaxis()->SetTitleSize(0.06);
  h_tot_true -> SetTitle("Total Background");
  h_tot_true -> GetXaxis()->SetTitle(plotting::GetAxisLabel(observable,0).c_str());
  h_tot_true -> GetYaxis()->SetTitle(plotting::GetAxisLabel(observable,1).c_str());
  h_tot_true -> GetXaxis()->CenterTitle();
  h_tot_true -> GetYaxis()->CenterTitle();
  double max = plotting::GetMaximum({h_tot_true,h_tot_est});
  h_tot_true->GetYaxis()->SetRangeUser(0,max);

  h_tot_est -> SetLineColor(kBlue);
  h_tot_est -> SetMarkerStyle(8);
  h_tot_est -> SetMarkerColor(kBlue);
  h_tot_est -> SetLineWidth(2);

  h_diff_true -> SetLineColor(kBlack);
  h_diff_true -> SetMarkerStyle(8);
  h_diff_true -> SetMarkerColor(kBlack);
  h_diff_true -> SetLineWidth(2);

  h_diff_true -> SetTitle("");
  h_diff_true -> GetXaxis()->SetTitle(plotting::GetAxisLabel(observable,0).c_str());
  h_diff_true -> GetYaxis()->SetTitle("#sigma_{BKG}/Est[%]");
  h_diff_true -> GetXaxis()->CenterTitle();
  h_diff_true -> GetYaxis()->CenterTitle();
  h_diff_true -> GetYaxis()->SetNdivisions(6);
  h_diff_true -> GetYaxis()->SetRangeUser(-50,50);
  h_diff_true -> GetXaxis()->SetLabelSize(0.08);
  h_diff_true -> GetYaxis()->SetLabelSize(0.08);
  h_diff_true -> GetXaxis()->SetTitleOffset(1.2);
  h_diff_true -> GetXaxis()->SetTitleSize(0.08);
  h_diff_true -> GetYaxis()->SetTitleSize(0.08);

  auto legend = new TLegend(0.15,0.6,0.35,0.8);
  legend->AddEntry(h_tot_true, "Total true background");
  legend->AddEntry(h_tot_est, "Total estimated background");
  h_tot_true -> GetYaxis()->SetTitle("#Events");

  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.015);
  pad1->SetLeftMargin(0.1);
  h_tot_true -> Draw("hist ERR");
  h_tot_est  -> Draw("hist ERR same");
  legend->Draw();

  // Top canvas - distribution
  // Bottom canvas - relative difference
  c->cd();
  pad2->Draw();
  pad2->cd();
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.1);
  pad2->SetTopMargin(0.05);
  h_diff_true->Draw("hist");
  c->SaveAs(output_file.c_str());

  return 0 ;
}
