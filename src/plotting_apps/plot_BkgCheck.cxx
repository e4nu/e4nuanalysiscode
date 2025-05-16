// __________________________________________________________________________
/* This app is used to plot the estimated vs true background obtained from  */
/* The output root file computed with bkg debug mode ("DebugBkg true")      */
/* An example of how to run the analysis with this mode is available in     */
/* ConfFiles/example_configuration.txt                                      */
// __________________________________________________________________________
#include <iostream>
#include <vector>
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
  std::string output_file = "bakground_debug";
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

  for( unsigned int j = 1 ; j < h_diff_true -> GetNbinsX() ; ++j ){
    double err = pow(h_tot_true->GetBinContent(j)-h_tot_est->GetBinContent(j),2) ;

    // Substract stat error
    err -= pow(h_tot_true->GetBinError(j),2) ;
    err -= pow(h_tot_est->GetBinError(j),2) ;

    if( err < 0 ) err = 0 ;
    err = sqrt(err);

    if( h_tot_est->GetBinContent(j) != 0 && h_tot_est->GetBinContent(j) > h_tot_est->Integral()*0.01 ) {
      h_diff_true->SetBinContent(j, err/h_tot_est->GetBinContent(j)*100);
      h_diff_true->SetBinError(j,0);
    } else {
      h_diff_true->SetBinContent(j,0);
      h_diff_true->SetBinError(j,0);
    }
  }

  h_tot_true -> SetLineColor(kBlack);
  h_tot_true -> SetMarkerStyle(8);
  h_tot_true -> SetMarkerColor(kBlack);
  h_tot_true -> SetLineWidth(2);
  h_tot_true -> GetYaxis()->SetLabelSize(0.06);
  h_tot_true -> GetYaxis()->SetTitleSize(0.06);
  h_tot_true -> SetTitle("Total Background");
  h_tot_true -> GetXaxis()->SetTitle(plotting::GetAxisLabel(observable,0).c_str());

  std::vector<TH1D*> vector_hists = {h_tot_true,h_tot_est};
  double max = plotting::GetMaximum(vector_hists);
  h_tot_true->GetYaxis()->SetRangeUser(0,max);

  h_tot_est -> SetLineColor(kBlue);
  h_tot_est -> SetMarkerStyle(8);
  h_tot_est -> SetMarkerColor(kBlue);
  h_tot_est -> SetLineWidth(2);

  h_diff_true -> SetLineColor(kBlue);
  h_diff_true -> SetMarkerStyle(8);
  h_diff_true -> SetMarkerColor(kBlue);
  h_diff_true -> SetLineWidth(2);

  h_diff_true -> SetTitle("");
  h_diff_true -> GetXaxis()->SetTitle(plotting::GetAxisLabel(observable,0).c_str());

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
  c->SaveAs((output_file+".root").c_str());


  // Compute systematics
  // Method 1 :
  TH1D * h_method1 = (TH1D*)h_diff_true->Clone();

  for( unsigned int j = 1 ; j < h_method1 -> GetNbinsX() ; ++j ){
    double fraction_j = h_total_bkg->GetBinContent(j)/h_signal->GetBinContent(j);
    if(h_total_bkg->GetBinContent(j)<500) fraction_j=0;
    if(h_signal->GetBinContent(j)<500) fraction_j=0;

    h_method1->SetBinContent(j,h_diff_true->GetBinContent(j)*fraction_j);
    h_method1->SetBinError(j,0);
  }

  // Method 2 :
  // Compute average
  double sum_weights =0, sum_err = 0;
  // Compute average error
  for( unsigned int j = 1 ; j < h_diff_true -> GetNbinsX() ; ++j ){
    double err = pow(h_tot_true->GetBinContent(j)-h_tot_est->GetBinContent(j),2) ;

    // Substract stat error
    err -= pow(h_tot_true->GetBinError(j),2) ;
    err -= pow(h_tot_est->GetBinError(j),2) ;

    if( err < 0 ) err = 0 ;
    err = sqrt(err);

    if( h_tot_est->GetBinContent(j) != 0 && h_tot_est->GetBinContent(j) > h_tot_est->Integral()*0.01 ) {
      sum_err += err /h_tot_est->GetBinContent(j) ;
      sum_weights += 1 ;
    }
  }

  double av_err = sum_err / sum_weights ;
  std::cout << " Weighted error [%]:" << sum_err / sum_weights * 100 << std::endl;

  TH1D * h_method2 = (TH1D*)h_diff_true->Clone();
  for( unsigned int j = 1 ; j < h_method1 -> GetNbinsX() ; ++j ){
    double fraction_j = h_total_bkg->GetBinContent(j)/h_signal->GetBinContent(j);
    if(h_total_bkg->GetBinContent(j)<500) fraction_j=0;
    if(h_signal->GetBinContent(j)<500) fraction_j=0;

    h_method2->SetBinContent(j,av_err*fraction_j*100);
    h_method2->SetBinError(j,0);
  }

  TFile *ROOTOut = new TFile((output_file+"_"+observable+"_syst.root").c_str(),"RECREATE");
  h_method1->SetName(("BkgSyst_Method1_"+observable).c_str());
  h_method1->Write();
  h_method2->SetName(("BkgSyst_Method2_"+observable).c_str());
  h_method2->Write();
  ROOTOut->Close();

  return 0 ;
}
