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

int main( int argc, char* argv[] ) {

  // Paths to ROOT files
  // const char* file1Path = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_Inclusive_1GeV_with_breakdown_dxsec_dRecoEnergyTransfer.root";
  // const char *file_fit = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/e4nuanalysiscode/data/InclusiveBodekChristy_C_1GeV.csv";
  // const char* file1Path = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_Inclusive_2GeV_with_breakdown_dxsec_dRecoEnergyTransfer.root";
  // const char *file_fit = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/e4nuanalysiscode/data/InclusiveBodekChristy_C_2GeV_29.5deg.csv";
  const char* file1Path = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_Inclusive_4GeV_with_breakdown_dxsec_dRecoEnergyTransfer.root";
  const char *file_fit = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/e4nuanalysiscode/data/InclusiveBodekChristy_C_4GeV.csv";

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
      std::cout << atof(col1.c_str()) << " " << atof(col2.c_str()) << std::endl;
    }
    //while (inputFile.good()) cout << (char) inputFile.get();
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

  // Visual comparison
  TCanvas* canvas = new TCanvas("canvas", "Histogram Comparison", 600, 600);
  gStyle->SetOptStat(0);

  // Adjust Y-axis range to fit both histograms
  double maxdata = h_data->GetMaximum();
  double maxgenie = h_genie->GetMaximum();
  double maxRange = std::max(maxdata, maxgenie) * 1.1; // Add 10% padding

  h_data->GetYaxis()->SetRangeUser(0, maxRange); // Set dynamic Y-axis range
  h_data->Draw("hist err");
  g->Draw("same");
  h_genie->Draw("hist err same");


  // Clean up
  canvas->SaveAs("Inclusive_comparison.pdf"); // Save the comparison plot
  canvas->SaveAs("Inclusive_comparison.root");
  return 0;
}
