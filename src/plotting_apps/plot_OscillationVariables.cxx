// __________________________________________________________________________
/*   */
// __________________________________________________________________________
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "utils/RadiativeCorrUtils.h"
#include "plotting/PlottingUtils.h"
#include "plotting/XSecUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TH2D.h"
#include "TGraph2D.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf;
using namespace e4nu::plotting;
/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting histograms..." << std::endl;

  if ( argc == 0 ) {
    std::cout << "Input files not specified. Abort..." << std::endl;
    return 0;
  }

  std::string input_file, legend_list ;
  std::string output_file = "histogram_data";
  std::string observable = "ECal";
  std::string analysis_key = "1p1pim";

  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-file",argc,argv)) {
      input_file = GetArg("input-file",argc,argv);
    } else { return 0 ;}

    if( ExistArg("output-file",argc,argv)) {
      output_file = GetArg("output-file",argc,argv);
    }
  }

  TCanvas * c  = new TCanvas("","",800,800);
  c->SetFillColor(0);
  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);
  // gStyle->SetPalette(kLightTemperature);

  GetMissingEnergyGraph( (input_file).c_str() );

  // Now you can draw the graph after all entries have been processed.
  graph_oscillations->SetTitle("3D Scatter Plot;Efl;HadSystemMass;MissingEnergy");
  graph_oscillations->SetMarkerStyle(20);
  graph_oscillations->SetMarkerColor(kRed); // Change color if desired
  graph_oscillations->Draw("COLZ");

  std::string output_name = "test_GraphOscillations" ;

  c->SaveAs((output_name+".root").c_str());
  c->SaveAs((output_name+".pdf").c_str());

  // Pt < 0.2
  if( graph_oscillations_1 ) {
    // gStyle->SetPalette(kLightTemperature);
    graph_oscillations_1->SetTitle("3D Scatter Plot;Efl;HadSystemMass;MissingEnergy");
    graph_oscillations_1->SetMarkerStyle(20);
    graph_oscillations_1->SetMarkerColor(kRed); // Change color if desired
    graph_oscillations_1->Draw("COLZ");
    output_name = "test_GraphOscillations_1" ;

    c->SaveAs((output_name+".root").c_str());
    c->SaveAs((output_name+".pdf").c_str());
  }

  if( graph_oscillations_2 ) {
    graph_oscillations_2->SetTitle("3D Scatter Plot;Efl;HadSystemMass;MissingEnergy");
    graph_oscillations_2->SetMarkerStyle(20);
    graph_oscillations_2->SetMarkerColor(kRed); // Change color if desired
    graph_oscillations_2->Draw("COLZ");

    output_name = "test_GraphOscillations_2" ;

    c->SaveAs((output_name+".root").c_str());
    c->SaveAs((output_name+".pdf").c_str());
  }

  if( graph_oscillations_3 ) {
    graph_oscillations_3->SetTitle("3D Scatter Plot;Efl;HadSystemMass;MissingEnergy");
    graph_oscillations_3->SetMarkerStyle(20);
    graph_oscillations_3->SetMarkerColor(kRed); // Change color if desired
    graph_oscillations_3->Draw("COLZ");

    output_name = "test_GraphOscillations_3" ;

    c->SaveAs((output_name+".root").c_str());
    c->SaveAs((output_name+".pdf").c_str());
  }

  // TFile* true_root_file = new TFile("/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_1p1pim_1GeV_with_breakdown_dxsec_dMissingEnergy.root","ROOT") ;
  // TFile* reco_root_file = new TFile("/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_1p1pim_1GeV_with_breakdown_dxsec_dCorrMissingEnergy.root","ROOT") ;
  // if( !true_root_file || !reco_root_file ) {
  //   std::cout << " ERROR: Files do not exist."<<std::endl;
  //   return 0 ;
  // }
  // TH1D * mc_data = (TH1D*)true_root_file->Get("MC_True");
  // TH1D * mc_recodata = (TH1D*)reco_root_file->Get("MC_True");
  // TH1D * true_data = (TH1D*)true_root_file->Get("Data");
  // TH1D * reco_data = (TH1D*)reco_root_file->Get("Data");
  // reco_data->SetMarkerColor(kRed);
  // reco_data->SetLineColor(kRed);
  // mc_recodata->SetLineColor(kRed);
  // mc_data->GetXaxis()->SetLabelSize(0.05);
  // mc_data->GetXaxis()->SetTitleSize(0.08);
  // mc_data->Draw("hist");
  // mc_recodata->Draw("hist same");
  // true_data->Draw("same");
  // reco_data->Draw("same");
  // c->SaveAs("testing.root");
  //
  //
  // TFile* true_data_root_2Dfile = new TFile("/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_1p1pim_1GeV_with_breakdown_dxsec_dMissingEnergy_vs_HadDeltaPT_data.root","ROOT") ;
  // TFile* true_mc_root_2Dfile = new TFile("/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_1p1pim_1GeV_with_breakdown_dxsec_dMissingEnergy_vs_HadDeltaPT_mc.root","ROOT") ;
  // TFile* reco_data_root_2Dfile = new TFile("/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_1p1pim_1GeV_with_breakdown_dxsec_dCorrMissingEnergy_vs_HadDeltaPT_data.root","ROOT") ;
  // TFile* reco_mc_root_2Dfile = new TFile("/Users/juliatenavidal/Desktop/Postdoc/e4nu/FinalPionProductionAnalysis/plots//TotalXSec/clas6analysis_1p1pim_1GeV_with_breakdown_dxsec_dCorrMissingEnergy_vs_HadDeltaPT_mc.root","ROOT") ;
  // if( !true_data_root_2Dfile || !true_mc_root_2Dfile || !reco_data_root_2Dfile || !reco_mc_root_2Dfile ) {
  //   std::cout << " ERROR: Files do not exist."<<std::endl;
  //   return 0 ;
  // }
  // TH2D * mc_data_2D = (TH2D*)true_mc_root_2Dfile->Get("MC_True");
  // TH2D * mc_recodata_2D = (TH2D*)reco_mc_root_2Dfile->Get("MC_True");
  // TH2D * true_data_2D = (TH2D*)true_data_root_2Dfile->Get("Data");
  // TH2D * reco_data_2D = (TH2D*)reco_data_root_2Dfile->Get("Data");
  //
  // std::vector<TH2D*> mc_hists = {mc_data_2D,mc_recodata_2D};
  // std::vector<TH2D*> data_hists = {true_data_2D,reco_data_2D};
  //
  // plotting::PlotProjection( mc_hists, data_hists, "MissingEnergy", "MissingEnergy", c, true, "", 0.2, 0.4 );
  // c->SaveAs("testing2D.root");

  return 0 ;

}
