#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <utility>

#include <string>
#include <map>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "THStack.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLine.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TGraphErrors.h"


void UseE4NuStyle () {

  TStyle * e4nuStyle= new TStyle("e4nu","e4nu plots style");

  e4nuStyle->SetFrameBorderMode(0);
  e4nuStyle->SetCanvasBorderMode(0);
  e4nuStyle->SetPadBorderMode(0);
  e4nuStyle->SetPadColor(0);
  e4nuStyle->SetCanvasColor(0);
  e4nuStyle->SetStatColor(0);
  e4nuStyle->SetFillColor(0);
  e4nuStyle->SetLegendBorderSize(1);

  e4nuStyle->SetPaperSize(20,26);
  e4nuStyle->SetPadTopMargin(0.05);
  e4nuStyle->SetPadRightMargin(0.05);
  e4nuStyle->SetPadBottomMargin(0.2);
  e4nuStyle->SetPadLeftMargin(0.2);

  e4nuStyle->SetTextFont(132);
  e4nuStyle->SetTextSize(0.08);
  e4nuStyle->SetLabelFont(132,"x");
  e4nuStyle->SetLabelFont(132,"y");
  e4nuStyle->SetLabelFont(132,"z");
  e4nuStyle->SetLabelSize(0.05,"x");
  e4nuStyle->SetTitleSize(0.06,"x");
  e4nuStyle->SetLabelSize(0.05,"y");
  e4nuStyle->SetTitleSize(0.06,"y");
  e4nuStyle->SetLabelSize(0.05,"z");
  e4nuStyle->SetTitleSize(0.06,"z");
  e4nuStyle->SetLabelFont(132,"t");
  e4nuStyle->SetTitleFont(132,"x");
  e4nuStyle->SetTitleFont(132,"y");
  e4nuStyle->SetTitleFont(132,"z");
  e4nuStyle->SetTitleFont(132,"t");
  e4nuStyle->SetTitleFillColor(0);
  e4nuStyle->SetTitleX(0.25);
  e4nuStyle->SetTitleFontSize(0.08);
  e4nuStyle->SetTitleFont(132,"pad");

  e4nuStyle->SetMarkerStyle(20);
  e4nuStyle->SetHistLineWidth(1.85);
  e4nuStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  e4nuStyle->SetErrorX(0.001);

  e4nuStyle->SetTitleBorderSize(1);
  e4nuStyle->SetTitleX(0);
  e4nuStyle->SetTitleH(.05);

  e4nuStyle->SetOptStat(0);
  e4nuStyle->SetOptFit(0);

  e4nuStyle->SetPadTickX(1);
  e4nuStyle->SetPadTickY(1);

  e4nuStyle->SetPalette(1,0);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
				   NCont);
  e4nuStyle->SetNumberContours(NCont);
  e4nuStyle->cd();

  return ;
}

void StandardFormatHist( TH1D * prediction, std::string x_label, std::string y_label, std::string title, int color, int style ) {

  UseE4NuStyle();
  prediction -> SetLineColor(color);
  prediction -> SetLineStyle(style);
  prediction -> SetLineWidth(2);

  prediction -> SetTitle(title.c_str());
  prediction -> GetXaxis() -> SetTitle(x_label.c_str());
  prediction -> GetYaxis() -> SetTitle(y_label.c_str());
  prediction -> GetXaxis()->CenterTitle();
  prediction -> GetYaxis()->CenterTitle();

  prediction->GetXaxis()->SetLabelSize(0.07);
  prediction->GetXaxis()->SetTitleOffset(1.);
  prediction->GetXaxis()->SetRangeUser(0.,1.25);
  prediction->GetXaxis()->SetTitleSize(0.07);
  prediction->GetXaxis()->SetNdivisions(6);
  prediction->GetYaxis()->SetNdivisions(5);
  prediction->GetYaxis()->SetTitleOffset(0.8);
  prediction->GetYaxis()->SetLabelSize(0.07);
  prediction->GetYaxis()->SetTitleSize(0.07);
  prediction->GetYaxis()->SetRangeUser(0., 1.);
  prediction->GetXaxis()->SetRangeUser(0., 10);
  return;
}

void e4nu_pmult(){

  std::unique_ptr<TFile> e4nu_analysis_nofiducial( TFile::Open("output_test_1p0pi_10M.root") );
  std::unique_ptr<TFile> e4nu_analysis_fiducial( TFile::Open("output_test_1p0pi_10M_fiducial.root") );

  auto tree_nofiducial = e4nu_analysis_nofiducial->Get<TTree>("MCTree");
  auto tree_fiducial = e4nu_analysis_fiducial->Get<TTree>("MCTree");

  int TrueNProtons,RecoNProtons;

  tree_nofiducial->SetBranchAddress("TrueNProtons", &TrueNProtons);
  tree_nofiducial->SetBranchAddress("RecoNProtons", &RecoNProtons);

  int TrueNProtonsWF,RecoNProtonsWF;
  tree_fiducial->SetBranchAddress("TrueNProtons", &TrueNProtonsWF);
  tree_fiducial->SetBranchAddress("RecoNProtons", &RecoNProtonsWF);

  TH1D * h_nofiducial_Truepmult = new TH1D("proton_mult_true", "proton_mult", 10, 0, 10 );
  TH1D * h_nofiducial_Recopmult = new TH1D("proton_mult_reco", "proton_mult", 10, 0, 10 );
  TH1D * h_fiducial_Recopmult = new TH1D("proton_mult_recowf", "proton_mult", 10, 0, 10 );

  for (int iEntry = 0; tree_nofiducial->LoadTree(iEntry) >= 0; ++iEntry) {
    // Load the data for the given tree entry
    tree_nofiducial->GetEntry(iEntry);
    h_nofiducial_Truepmult->Fill(TrueNProtons);
    h_nofiducial_Recopmult->Fill(RecoNProtons);
  }

  h_nofiducial_Truepmult->Scale(1./h_nofiducial_Truepmult->GetEntries());
  h_nofiducial_Recopmult->Scale(1./h_nofiducial_Recopmult->GetEntries());

  for (int iEntry = 0; tree_fiducial->LoadTree(iEntry) >= 0; ++iEntry) {
    // Load the data for the given tree entry
    tree_fiducial->GetEntry(iEntry);
    h_fiducial_Recopmult->Fill(RecoNProtonsWF);
  }
  h_fiducial_Recopmult->Scale(1./h_fiducial_Recopmult->GetEntries());

  TCanvas * c  = new TCanvas("c","c",800,600);

  h_nofiducial_Truepmult-> SetLineColor( kBlack ) ;
  h_nofiducial_Recopmult-> SetLineColor( kRed ) ;
  h_fiducial_Recopmult-> SetLineColor( kGreen ) ;

  StandardFormatHist( h_nofiducial_Truepmult, "Proton multiplicity", "Nentries/Total", "", kBlack, 1 )  ;
  StandardFormatHist( h_nofiducial_Recopmult, "Proton multiplicity", "Nentries/Total", "", kRed, 1 )  ;
  StandardFormatHist( h_fiducial_Recopmult, "Proton multiplicity", "Nentries/Total", "", kBlue, 1 )  ;

  h_nofiducial_Truepmult->Draw("HIST");
  h_nofiducial_Recopmult ->Draw("HIST SAME");
  h_fiducial_Recopmult->Draw("HIST SAME");


  TLegend* leg = new TLegend(.5,.75, 0.95,0.95,"");
  leg->SetFillStyle(0);

  leg -> AddEntry( h_nofiducial_Truepmult, "True MC proton multiplicity");
  leg -> AddEntry( h_nofiducial_Recopmult, "Reconstructed proton multiplicity before fiducial cut");
  leg -> AddEntry( h_fiducial_Recopmult, "Reconstructed proton multiplicity after fiducial cut");

  leg -> Draw() ;
  c->SaveAs("proton_mult_study.pdf");
  c->SaveAs("proton_mult_study.root");


}
