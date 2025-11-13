#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>
#include <TStyle.h>
#include <TLegend.h>
#include <vector>
#include <string>
#include <sstream>
#include "plotting/PlottingUtils.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::plotting;

int main( int argc, char* argv[] ) {

  if ( argc == 0 ) {
    std::cout << "Input files not specified. Abort..." << std::endl;
    return 0;
  }

  std::vector<string> input_files, legend_list ;
  std::string output_file = "histogram_analized_data";
  std::string correction = "";
  bool show_ratio = false, log_scale = false;
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

    if( ExistArg("correction",argc,argv) ){
      correction = GetArg("correction",argc,argv);
      show_ratio = true;
    }
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
    if( ExistArg("plot-ratio",argc,argv)) {
      show_ratio = true ;
    }
    if( ExistArg("log-scale",argc,argv)) {
      log_scale = true ;
    }
  }

  // Open the ROOT files
  TFile* file_default = TFile::Open(input_files[0].c_str());
  TFile* file_corrected = TFile::Open(input_files[1].c_str());

  if (!file_default || !file_corrected || file_default->IsZombie() || file_corrected->IsZombie()) {
    std::cerr << "Error: Could not open one of the files." << std::endl;
    return 1;
  }

  // Retrieve histograms
  TH1D* hist_default = dynamic_cast<TH1D*>(file_default->Get("Corrected_Event_Rate_Data"));
  TH1D* hist_corrected = dynamic_cast<TH1D*>(file_corrected->Get("Corrected_Event_Rate_Data"));

  if (!hist_default || !hist_corrected) {
    std::cerr << "Error: Could not retrieve one of the histograms." << std::endl;
    file_default->Close();
    file_corrected->Close();
    return 1;
  }

  hist_corrected->SetLineColor(kBlack);
  hist_corrected->SetMarkerColor(kBlack);
  hist_default->SetLineColor(kOrange + 1);
  hist_default->SetMarkerColor(kOrange + 1);

  if( show_ratio ) {
    hist_default->GetXaxis()->SetLabelSize(0.);
    hist_default->GetXaxis()->SetTitleSize(0.);
    hist_corrected->GetXaxis()->SetLabelSize(0.);
    hist_corrected->GetXaxis()->SetTitleSize(0.);
  }

  // Visual comparison
  TCanvas* canvas = new TCanvas("canvas", "Histogram Comparison", 800, 600);
  TPad* pad1 ;
  TPad* pad2 ;
  gStyle->SetOptStat(0);

  if ( show_ratio ) {
    pad1 = new TPad("pad1", "Pad 1", 0.0, 0.3, 1, 1.0) ;
    pad2 = new TPad("pad2", "Pad 2", 0.0, 0., 1.0, 0.3);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.1);
    pad1->SetTopMargin(0.1);
    pad2->SetBottomMargin(0.4); // Reduce margin to avoid overlap
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.1);
    pad1->SetBottomMargin(0.02); // Reduce margin to avoid overlap
    pad2->SetTopMargin(0.02);
    pad1->Draw();
    pad2->Draw();
    // Pad 1 (Histograms)
    pad1->cd();
  } else {
    pad1 = new TPad("pad1", "Pad 1", 0.0, 0.0, 1, 1.0) ;
    pad1->SetLeftMargin(0.25);
    pad1->SetRightMargin(0.05);
    pad1->SetBottomMargin(0.2);
    pad1->SetTopMargin(0.1); // Reduce margin to avoid overlap
    pad1->Draw();

    // Pad 1 (Histograms)
    pad1->cd();

  }
  if( log_scale ) pad1->SetLogy();

  // Adjust Y-axis range to fit both histograms
  double maxhist_default = hist_default->GetMaximum();
  double maxhist_corrected = hist_corrected->GetMaximum();
  double maxRange = std::max(maxhist_default, maxhist_corrected) * 1.1; // Add 10% padding
  double minRange = 0;

  if( log_scale ) {
    maxRange *= 2;
    minRange = 10 ;
  }
  hist_default->GetYaxis()->SetTitleOffset(0.8);
  hist_default->GetYaxis()->SetTitle("Corr. Counts/Bin Width");
  hist_default->GetYaxis()->SetRangeUser(minRange, maxRange); // Set dynamic Y-axis range
  hist_default->Draw("hist err");
  hist_corrected->Draw("hist err same");


  auto legend = new TLegend(0.2,0.65,0.75,0.85);
  if( legend_list.size() == input_files.size() ) {
    legend->SetBorderSize(0);
    legend->SetTextFont(132);
    legend->SetTextSize(0.08);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    legend->AddEntry(hist_default,legend_list[0].c_str());
    legend->AddEntry(hist_corrected,legend_list[1].c_str());
    legend->Draw();
  }

  if(show_ratio && correction == false ){
    // Pad 2 (Difference Histogram)
    pad2->cd();
    TH1D* herr1 = (TH1D*)hist_corrected->Clone();

    // Compute RMS per bin
    double max_hist = hist_corrected->GetMaximum();
    for (int i = 1; i <= hist_corrected->GetNbinsX(); ++i) {
      double x1 = hist_corrected->GetBinContent(i); // from default file
      double x2 = hist_default->GetBinContent(i);   // from shifted file
      double binErrordef = hist_corrected->GetBinError(i);      // Get error of bin i
      // Mean per bin
      double xmean = 0.5 * (x1 + x2);

      if( xmean < 0.2 * max_hist ) {
        herr1->SetBinContent(i, 0);
        continue ;
      }
      // RMS formula
      double xRMS = ((x1 - xmean) * (x1 - xmean) + (x2 - xmean) * (x2 - xmean)) / 2.0 ;// - pow(binErrordef, 2)  ;
      if( xRMS <= 0 ) xRMS = 0;
      // Optionally, express as a relative RMS (%)
      double value = (xmean != 0.0) ? (sqrt(xRMS) / xmean * 100.0) : 0.0;

      herr1->SetBinContent(i, value);
    }

    // Compute average accross Bins
    double sum = 0.0;
    double count = 0;

    for (int i = 1; i <= herr1->GetNbinsX(); ++i) {
      double val = herr1->GetBinContent(i);
      double width = herr1->GetBinWidth(i);
      if (val > 0) { // Only include non-empty bins
        sum += val * width;
        count += width ;
      }
    }
    std::cout << " average is " << sum / count << std::endl;

    herr1->SetLineColor(kOrange + 1);
    herr1->SetLineWidth(3);
    herr1->GetYaxis()->SetRangeUser(0, 15);
    herr1->GetYaxis()->SetLabelSize(0.2);
    herr1->GetXaxis()->SetLabelSize(0.2);
    herr1->GetYaxis()->SetTitleSize(0.16);
    herr1->GetXaxis()->SetTitleSize(0.18);
    herr1->GetYaxis()->SetTitleOffset(0.3);
    herr1->GetXaxis()->SetTitleOffset(1.);
    herr1->Draw("HIST");
    herr1->GetYaxis()->SetTitle("x_{RMS}/Def[%]");
  } else if ( correction != "" ) {
    pad2->cd();

    TFile* file_correction = TFile::Open(correction.c_str());
    if( !file_correction ) return 0 ;
    TH1D* correction_factor = dynamic_cast<TH1D*>(file_correction->Get("Acceptance"));
    if( !correction_factor ) return 0 ;
    correction_factor->GetYaxis()->SetLabelSize(0.2);
    correction_factor->GetXaxis()->SetLabelSize(0.2);
    correction_factor->GetYaxis()->SetTitleSize(0.16);
    correction_factor->GetXaxis()->SetTitleSize(0.18);
    correction_factor->GetYaxis()->SetTitleOffset(0.3);
    correction_factor->GetXaxis()->SetTitleOffset(1.);
    correction_factor->Draw("hist err");
    correction_factor->GetYaxis()->SetTitle("Corr.");
  }

  // Clean up
  canvas->SaveAs("comparison.pdf"); // Save the comparison plot
  canvas->SaveAs("comparison.root");
  return 0;
}
