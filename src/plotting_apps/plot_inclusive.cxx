#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <TGraph.h>
#include <vector>
#include <string>

double compute_error_scale_chi2(TH1* h_data, TGraph* g) {
  int nbins = h_data->GetNbinsX();
  double chi2sum = 0.0;
  int used_bins = 0;

  for (int ib = 1; ib <= nbins; ++ib) {
    double x = h_data->GetBinCenter(ib);
    double hval = h_data->GetBinContent(ib);
    double herr = h_data->GetBinError(ib);

    if (herr <= 0) continue;             // skip bins with zero (or invalid) error
    if (h_data->GetBinContent(ib) == 0 && h_data->GetEntries() == 0) continue; // optionally skip empty hist

    double gval = g->Eval(x);            // linear interp by default
    double resid = hval - gval;
    chi2sum += (resid * resid) / (herr * herr);
    ++used_bins;
  }

  if (used_bins <= 0) {
    std::cerr << "No valid bins for chi2 calculation. Returning scale=1.\n";
    return 1.0;
  }

  double ndf = used_bins; // minus number of fitted params if any (e.g. minus 1 if you fit a normalization)
  double f = std::sqrt(chi2sum / ndf);
  std::cout << " Chi2 /ndf " << chi2sum / ndf << std::endl;
  if (f < 1.0) f = 1.0; // optional: do not shrink errors
  return f;
}

// apply the scale to histogram (multiply each bin error by f)
void add_histogram_errors(TH1* h_data, double f) {
  int nbins = h_data->GetNbinsX();
  for (int ib = 1; ib <= nbins; ++ib) {
    double olderr = h_data->GetBinError(ib);
    h_data->SetBinError(ib, sqrt( pow(olderr,2) + pow( f * h_data->GetBinContent(ib),2)) );
  }
}

double ComputeRelativeNormError(TH1* h_data, TGraph* g)
{
  double tol = 0.1;
  if (!h_data || !g) throw std::runtime_error("Null pointer to TH1D or TGraph.");

  auto chi2_over_ndf = [&](double r) {
    int nbins = h_data->GetNbinsX();
    double chi2sum = 0.0;
    int used_bins = 0;

    for (int ib = 1; ib <= nbins; ++ib) {
      double x = h_data->GetBinCenter(ib);
      double hval = h_data->GetBinContent(ib);
      double herr = h_data->GetBinError(ib);

      if (herr <= 0) continue;             // skip bins with zero (or invalid) error
      if (h_data->GetBinContent(ib) == 0 && h_data->GetEntries() == 0) continue; // optionally skip empty hist
      if( x < 0.15 ) continue ;
      double gval = g->Eval(x);            // linear interp by default
      double resid = hval - gval;
      chi2sum += (resid * resid) / (herr * herr + (r*hval)*(r*hval));
      ++used_bins;
    }
    return chi2sum / used_bins;
  };

  // Check chi2 at r=0
  double chi2_0 = chi2_over_ndf(0.0);
  std::cout << " chi2_over_ndf(r_high) " << chi2_over_ndf(0) << std::endl;
  if (chi2_0 <= 1.0) return 0.0;

  // Bracket the solution
  double r_low = 0.0;
  double r_high = 1.0;

  std::cout << " chi2_over_ndf(r_high) " << chi2_over_ndf(r_high) << std::endl;
  while (chi2_over_ndf(r_high) > 1.0) {
    r_high +=0.005;
    std::cout << " chi2_over_ndf(r_high) " << chi2_over_ndf(r_high) << std::endl;
    if (r_high > 1e8) throw std::runtime_error("Failed to bracket solution for r.");
  }

  // Bisection
  for (int iter = 0; iter < 200; iter++) {
    double r_mid = 0.5 * (r_low + r_high);
    double val = chi2_over_ndf(r_mid);

    if (std::abs(r_high - r_low) < tol)
    return r_mid;

    if (val > 1.0)
    r_low = r_mid;
    else
    r_high = r_mid;
  }

  return 0.5 * (r_low + r_high);
}

int main( int argc, char* argv[] ) {
  bool substract = false ;

  // Paths to ROOT files
  // const char* file1Path = "/Users/juliatenavidal/Desktop/Research/E4Nu/FinalPionProductionAnalysis//plots_new//TotalXSec/clas6analysis_Inclusive_1GeV_with_breakdown_dxsec_dRecoEnergyTransfer.root";
  // const char *file_fit = "/Users/juliatenavidal/Desktop/Research/E4Nu/FinalPionProductionAnalysis/e4nuanalysiscode/data/InclusiveBodekChristy_C_1GeV.csv";

  // const char* file1Path = "/Users/juliatenavidal/Desktop/Research/E4Nu/FinalPionProductionAnalysis//plots_new//TotalXSec/clas6analysis_Inclusive_2GeV_with_breakdown_dxsec_dRecoEnergyTransfer.root";
  // const char *file_fit = "/Users/juliatenavidal/Desktop/Research/E4Nu/FinalPionProductionAnalysis/e4nuanalysiscode/data/InclusiveBodekChristy_C_2GeV_29.5deg.csv";

  const char* file1Path = "/Users/juliatenavidal/Desktop/Research/E4Nu/FinalPionProductionAnalysis//plots_new//TotalXSec/clas6analysis_Inclusive_4GeV_with_breakdown_dxsec_dRecoEnergyTransfer.root";
  const char *file_fit = "/Users/juliatenavidal/Desktop/Research/E4Nu/FinalPionProductionAnalysis/e4nuanalysiscode/data/InclusiveBodekChristy_C_4GeV.csv";

  // Fill Graph
  std::ifstream inputFile ;
  inputFile.open(file_fit);
  std::vector<double> omega, xsec ;
  if (inputFile.is_open()){
    std::string line;

    while (std::getline(inputFile, line)) {
      std::stringstream ss(line);
      std::string col1, col2;

      std::getline(ss, col1, ';');
      std::getline(ss, col2, ';');

      omega.push_back(atof(col1.c_str()));
      xsec.push_back(atof(col2.c_str()));
    }
    inputFile.close();
  }
  auto g = new TGraph(xsec.size()-1,&omega[0],&xsec[0]);

  // Open the ROOT files
  TFile* file1 = TFile::Open(file1Path);

  if (!file1 || file1->IsZombie() ) {
    std::cerr << "Error: Could not open one of the files." << std::endl;
    return 1;
  }

  // Retrieve histograms
  TH1* h_data = dynamic_cast<TH1*>(file1->Get("Data"));
  TH1* h_genie = dynamic_cast<TH1*>(file1->Get("MC_True"));

  if (!h_data || !h_genie) {
    std::cerr << "Error: Could not retrieve one of the histograms." << std::endl;
    file1->Close();
    return 1;
  }

  h_data->SetLineColor(kBlack);
  h_genie->SetLineColor(kBlue);
  h_data->SetMarkerColor(kBlack);
  g->SetLineColor(kRed);

  // Clone data and add normalization uncertainty
  TH1D * h_data_wErr = (TH1D*) h_data->Clone();

  // Compute normalization error from data and fit:
  double norm_err = ComputeRelativeNormError( h_data, g );
  std::cout <<  "ComputeRelativeNormError "<< ComputeRelativeNormError(h_data,g) * 100 << " % "<< std::endl;

  add_histogram_errors(h_data_wErr,norm_err);

  // Visual comparison
  TCanvas* canvas = new TCanvas("canvas", "Histogram Comparison", 600, 600);
  gStyle->SetOptStat(0);

  // Adjust Y-axis range to fit both histograms
  double maxdata = h_data->GetMaximum();
  double maxgenie = h_genie->GetMaximum();
  double maxRange = std::max(maxdata, maxgenie) * 1.1; // Add 10% padding

  h_data->GetYaxis()->SetRangeUser(0, maxRange); // Set dynamic Y-axis range

  // draw uncertainty band
  h_data_wErr->Draw("E1");

  // draw data points on top
  h_data->SetMarkerStyle(20);
  h_data_wErr->SetMarkerSize(1.2);
  h_data->SetMarkerSize(1.2);
  h_data->SetLineColor(kBlack);
  h_data_wErr->SetLineColor(12);


  h_data->Draw("E1 SAME");

  // fit
  g->Draw("L SAME");

  //h_genie->Draw("hist err same");

  // Clean up
  canvas->SaveAs("Inclusive_comparison.pdf"); // Save the comparison plot
  canvas->SaveAs("Inclusive_comparison.root");
  return 0;
}
