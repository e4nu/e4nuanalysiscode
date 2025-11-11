#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "plotting/PlottingUtils.h"
#include "utils/DetectorUtils.h"
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
// * input_file : comma sepparated list of data root files    //
// * output-file : file to store plots (in root, pdf... format)//
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting histograms..." << std::endl;

  if ( argc == 0 ) {
    std::cout << "Input files not specified. Abort..." << std::endl;
    return 0;
  }

  string input_file ;
  string topology = "1p1pim" ;
  std::string output_file = "hadron_sector";

  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-file",argc,argv)) {
      input_file = GetArg("input-file",argc,argv);
    } else { return 0 ;}

    if( ExistArg("output-file",argc,argv)) {
      output_file = GetArg("output-file",argc,argv);
    }
    if( ExistArg("topology",argc,argv)) {
      topology = GetArg("topology",argc,argv);
    }
  }

  // Color palette
  std::vector<int> color_list = {kBlack,kBlue-4,kOrange+1,kGreen+3,kMagenta+2,kRed-4,kSpring+7,kTeal-6,kCyan+3,kAzure-1,kBlue-8,kViolet-4,kMagenta-10,kPink-9} ;

  TCanvas * c  = new TCanvas("","",800,800);
  c->SetFillColor(0);

  TFile* in_root_file = new TFile(input_file.c_str(),"ROOT");
  if( !in_root_file ) { std::cout << "ERROR: the "<< input_file <<" does not exist." << std::endl; return 0;}

  TTree* in_tree = (TTree*)in_root_file->Get("CLAS6Tree") ;
  if( !in_tree ) { std::cout << "ERROR: the threes do not exist." <<std::endl; return 0;}

  TH1D* hist_ps = new TH1D( "", "", 6, -0.5,5.5);
  TH1D* hist_pis = new TH1D( "", "", 6, -0.5,5.5);

  // OBSERVABLE DEFINITION:
  long NEntries = in_tree -> GetEntries() ;
  plotting::SetAnalysisBranch( in_tree ) ;

  for( int j = 0 ; j < NEntries ; ++j ) {
    in_tree->GetEntry(j) ;
    double w = EventWght * AccWght ;

    // Convert proton phi to radians:
    proton_phi *= TMath::DegToRad() ;
    int proton_sector = std::abs((int)utils::GetSector(proton_phi) - (int)ElectronSector);
    hist_ps->Fill(proton_sector,w);

    double pion_phi = 0;
    if( topology == "1p1pip" ) pion_phi = pip_phi;
    else if( topology == "1p1pim" ) pion_phi = pim_phi;

    pion_phi *= TMath::DegToRad();
    int pi_sector = std::abs((int) utils::GetSector(pion_phi) - (int)ElectronSector);
    hist_pis->Fill(pi_sector,w);
  }

  // Normalise
  hist_ps->Scale(1./hist_ps->Integral());
  hist_pis->Scale(1./hist_pis->Integral());

  // Set axis names
  hist_ps->GetYaxis()->SetTitle("Normalized Event Rate");
  hist_pis->GetYaxis()->SetTitle("Normalized Event Rate");

  hist_ps->GetXaxis()->SetTitle("Hadron Sector - e-Sector");
  hist_pis->GetXaxis()->SetTitle("Hadron Sector - e-Sector");

  hist_ps->Draw("hist err");
  hist_pis->SetLineColor(kRed);
  hist_pis->Draw("hist err same");

  auto legend = new TLegend(0.2,0.65,0.55,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(132);
  legend->SetTextSize(0.08);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03);
  legend->AddEntry(hist_ps,"Proton");
  legend->AddEntry(hist_pis,"Pion");
  legend->Draw();

  std::string output_name = output_file ;

  c->SaveAs((output_file+".root").c_str());
  c->SaveAs((output_file+".pdf").c_str());

  return 0 ;

}
