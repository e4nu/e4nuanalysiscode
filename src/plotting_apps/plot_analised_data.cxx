#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>
#include <TStyle.h>

int main( int argc, char* argv[] ) {

  // Paths to ROOT files
  const char* file1Path = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots_new//EventRate/clas6analysis_1p1pim_4GeV_corr_event_rate_Nevents_ECal.root";
  const char* filedefPath = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots_new//EventRate/clas6analysis_1p1pim_4GeV_3DegShift_corr_event_rate_Nevents_ECal.root";

  // Open the ROOT files
  TFile* file1 = TFile::Open(file1Path);
  TFile* filedef = TFile::Open(filedefPath);

  if (!file1 || !filedef || file1->IsZombie() || filedef->IsZombie()) {
    std::cerr << "Error: Could not open one of the files." << std::endl;
    return 1;
  }

  // Retrieve histograms
  TH1* hist1 = dynamic_cast<TH1*>(file1->Get("Corrected_Event_Rate_Data"));
  TH1* histdef = dynamic_cast<TH1*>(filedef->Get("Corrected_Event_Rate_Data"));

  if (!hist1 || !histdef) {
    std::cerr << "Error: Could not retrieve one of the histograms." << std::endl;
    file1->Close();
    filedef->Close();
    return 1;
  }
  histdef->SetLineColor(kBlack);
  histdef->SetMarkerColor(kBlack);
  hist1->SetLineColor(kOrange + 1);
  hist1->SetMarkerColor(kOrange + 1);
  hist1->GetXaxis()->SetLabelSize(0.);
  hist1->GetXaxis()->SetTitleSize(0.);
  histdef->GetXaxis()->SetLabelSize(0.);
  histdef->GetXaxis()->SetTitleSize(0.);
  // Visual comparison
  TCanvas* canvas = new TCanvas("canvas", "Histogram Comparison", 800, 600);
  TPad* pad1 = new TPad("pad1", "Pad 1", 0.0, 0.3, 1, 1.0);
  TPad* pad2 = new TPad("pad2", "Pad 2", 0.0, 0., 1.0, 0.3);

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
  gStyle->SetOptStat(0);

  // Adjust Y-axis range to fit both histograms
  double maxHist1 = hist1->GetMaximum();
  double maxHistdef = histdef->GetMaximum();
  double maxRange = std::max(maxHist1, maxHistdef) * 1.1; // Add 10% padding

  hist1->GetYaxis()->SetRangeUser(0, maxRange); // Set dynamic Y-axis range
  hist1->Draw("hist err");
  histdef->Draw("hist err same");

  // Pad 2 (Difference Histogram)
  pad2->cd();
  TH1D* herr1 = (TH1D*)histdef->Clone();

  // Compute difference
  for (int i = 1; i <= histdef->GetNbinsX(); ++i) { // Loop from bin 1 to nBins (excluding underflow and overflow bins)
    double binContentdef = histdef->GetBinContent(i);  // Get content of bin i
    double binErrordef = histdef->GetBinError(i);      // Get error of bin i

    double binContent1 = hist1->GetBinContent(i);  // Get content of bin i
    double binError1 = hist1->GetBinError(i);      // Get error of bin i

    double err2_1 = pow(binContentdef - binContent1, 2)/12 - pow(binErrordef, 2) - pow(binError1, 2);
    if (err2_1 < 0) err2_1 = 0;

    double err_1 = sqrt(err2_1);
    if (binContentdef != 0)
    herr1->SetBinContent(i, err_1 / binContentdef * 100);
    else
    herr1->SetBinContent(i, 0);
  }

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
  herr1->GetYaxis()->SetTitle("#sigma_{fid}/Def[%]");

  // Clean up
  canvas->SaveAs("comparison.pdf"); // Save the comparison plot
  canvas->SaveAs("comparison.root");
  return 0;
}
