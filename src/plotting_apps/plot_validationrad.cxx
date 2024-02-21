// __________________________________________________________________________
/* This app is used to plot the validation figure from the nature paper    */
/* We use the ConfFiles/mc_conf/clas6mc_1panalysis_eH_4325MeV.txt          */
/* to compute the MC prediction. The data is stored in a root file         */
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
#include "TH3D.h"
#include "TLegend.h"
#include <TRandom.h>

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf;
using namespace e4nu::plotting;
/////////////////////////////////////////////////////////////////
// Options:                                                    //
// * input-file : analised MC file                             //
/////////////////////////////////////////////////////////////////

double DeltaEMin = -4.;// %
double DeltaEMax = 4.;// %
double DeltaENBins = 8.;
double DeltaEStep = fabs(DeltaEMax - DeltaEMin) / DeltaENBins;

double DeltaPMin = -4.;// %
double DeltaPMax = 4.;// %
double DeltaPNBins = 8.;
double DeltaPStep = fabs(DeltaEMax - DeltaEMin) / DeltaENBins;

double YPTarEMin = -0.03;// %
double YPTarEMax = 0.03;// %
double YPTarENBins = 6.;
double YPTarEStep = fabs(YPTarEMax - YPTarEMin) / YPTarENBins;

double YPTarPMin = -0.03;// %
double YPTarPMax = 0.03;// %
double YPTarPNBins = 6.;
double YPTarPStep = fabs(YPTarPMax - YPTarPMin) / YPTarPNBins;

double XPTarEMin = -0.05;// %
double XPTarEMax = 0.05;// %
double XPTarENBins = 10.;
double XPTarEStep = fabs(XPTarEMax - XPTarEMin) / XPTarENBins;

double XPTarPMin = -0.05;// %
double XPTarPMax = 0.05;// %
double XPTarPNBins = 10.;
double XPTarPStep = fabs(XPTarPMax - XPTarPMin) / XPTarPNBins;

double ZvertEMin = -0.09;// m
double ZvertEMax = 0.09;// m
double ZvertENBins = 18.;
double ZvertEStep = fabs(ZvertEMax - ZvertEMin) / ZvertENBins;

double ZvertPMin = -0.09;// m
double ZvertPMax = 0.09;// m
double ZvertPNBins = 18.;
double ZvertPStep = fabs(ZvertPMax - ZvertPMin) / ZvertPNBins;

double ElectronThetaCV = 0.310698;
double ProtonThetaCV = 0.852082;

double VertexZMin = -0.05; // m
double VertexZMax = 0.05; // m

//double DeltaP_CV = 2.225;  // 4%
double DeltaP_CV = 1.75;  // 4%

int main( int argc, char* argv[] ) {
  
  std::cout << "Plot validation of radiative correction using QEL data..." << std::endl;

  if ( argc == 0 ) {
    std::cout << "Input files not specified. Abort..." << std::endl;
    return 0;
  }

  string input_file ;
  std::string output_file = "radcorr_verification.root";
  std::string observable = "ECal";

  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("input-file",argc,argv)) {
      input_file = GetArg("input-file",argc,argv);
    } else { return 0 ;}
  }

  TCanvas * c  = new TCanvas("","",800,800);
  c->SetFillColor(0);
  // Top canvas - distribution
  // Bottom canvas - relative difference
  TPad *pad1 = new TPad("pad1","",0,0.4,1,1);

  char * env = std::getenv("E4NUANALYSIS") ;
  std::string path( env ) ;
  path += "/" ;

  std::cout << " Data from : "<<path<<"data/1H_4_325_Data_Plots_FSI_em.root" << std::endl;
  TFile * data = new TFile( (path+"data/1H_4_325_Data_Plots_FSI_em.root").c_str(), "ROOTFile" ); 
  TFile * mc = new TFile( input_file.c_str(), "ROOTFile" );
  if( !data ) return 0 ; 
  if( !mc ) return 0 ; 

  TH1D* hist_data = (TH1D*) data->Get("ECalRecoPlot");
  if( !hist_data ) return 0 ;
  
  TTree * mc_tree = (TTree*) mc->Get("MCCLAS6Tree");
  TH1D* hist_mc = (TH1D*) hist_data->Clone();
  hist_mc->Reset();

  // Access acceptance correction
  string acc_path = path+"/data/AcceptanceMaps/CLAS6";
  TFile* acc_e = new TFile((acc_path+"/AcceptanceMap_e_TH3D.root").c_str() );
  TH3D* h3_e = (TH3D*)(acc_e->Get("h3"));
  TFile* acc_p = new TFile((acc_path+"/AcceptanceMap_p_TH3D.root").c_str() );
  TH3D* h3_p = (TH3D*)(acc_p->Get("h3"));

  double ECal = 0 ; 
  double TotWeight = 1 ;
  bool QEL = false ;
  long NEntries = mc_tree -> GetEntries() ;
  double pfl, pfl_phi, pfl_theta ; 
  double proton_mom, proton_phi, proton_theta ;
  mc_tree -> SetBranchAddress("TotWeight",&TotWeight);
  mc_tree -> SetBranchAddress("ECal",&ECal);
  mc_tree -> SetBranchAddress("QEL",&QEL);
  mc_tree -> SetBranchAddress("pfl",&pfl);
  mc_tree -> SetBranchAddress("pfl_theta",&pfl_theta);
  mc_tree -> SetBranchAddress("pfl_theta",&pfl_theta);
  mc_tree -> SetBranchAddress("proton_mom",&proton_mom);
  mc_tree -> SetBranchAddress("proton_theta",&proton_theta);
  mc_tree -> SetBranchAddress("proton_phi",&proton_phi);

  for( int j = 0 ; j < NEntries ; ++j ) {
    mc_tree->GetEntry(j) ;

    double FSElectronMag = pfl;
    double FSElectronTheta = pfl_theta;
    double FSElectronTheta_Deg = pfl_theta * 180. / TMath::Pi();
    double FSElectronCosTheta = cos(FSElectronTheta);
    double FSElectronPhi = pfl_phi;
    double FSElectronPhi_Deg = FSElectronPhi * 180. / TMath::Pi();
    //if (FSElectronPhi_Deg < 0 ) { FSElectronPhi_Deg += 360.;}
    // if (FSElectronPhi_Deg > 360. ) { FSElectronPhi_Deg -= 360.;}
    
    double FSProtonMag = proton_mom;
    double FSProtonTheta = proton_theta ; 
    double FSProtonTheta_Deg = FSProtonTheta * 180. / TMath::Pi();
    double FSProtonCosTheta = cos(FSProtonTheta);
    double FSProtonPhi = proton_phi;
    double FSProtonPhi_Deg = FSProtonPhi * 180. / TMath::Pi();
    //    if (FSProtonPhi_Deg < 0 ) { FSProtonPhi_Deg += 360.;}
    // if (FSProtonPhi_Deg > 360. ) { FSProtonPhi_Deg -= 360.;}


    double ElectronCV = 3.54334, ElectronCVReso = 0.045; 
    double ProtonCV = 1.4805, ProtonCVReso = 0.045;
    
    double delta_e = (FSElectronMag - ElectronCV) / ElectronCV * 100; // %
    double delta_p = (FSProtonMag - ProtonCV) / ProtonCV * 100; // %
    
    int IDeTrial = (delta_e - DeltaEMin)/DeltaEStep;
    int IDpTrial = (delta_p - DeltaPMin)/DeltaPStep;
    
    double YPTarE = FSElectronTheta - ElectronThetaCV;
    double YPTarP = FSProtonTheta - ProtonThetaCV;
  
    int YPTarETrial = (YPTarE - YPTarEMin)/YPTarEStep;
    int YPTarPTrial = (YPTarP - YPTarPMin)/YPTarPStep;
    
    //int XPTarETrial = YPTarENBins / 2;
    //int XPTarPTrial = YPTarPNBins / 2;
    
    double ZVertexTrial = gRandom->Uniform(VertexZMin,VertexZMax);
    
    int ZvertETrial = (ZVertexTrial - ZvertEMin)/ZvertEStep;
    int ZvertPTrial = (ZVertexTrial - ZvertPMin)/ZvertPStep;
    
    double AccWeight_e = h3_e->GetBinContent(IDeTrial+1,YPTarETrial+1,ZvertETrial+1);
    double AccWeight_p = h3_p->GetBinContent(IDpTrial+1,YPTarPTrial+1,ZvertPTrial+1);
        
    std::cout << " AccWeight_e " << AccWeight_e << std::endl;
    std::cout << " AccWeight_p " << AccWeight_p << std::endl;

    hist_mc -> Fill( ECal, TotWeight * AccWeight_e * AccWeight_p ) ; 
    //    std::cout << " ECal " << ECal << " w " << TotWeight << std::endl;
  }
  double data_integral = hist_data->Integral();
  double mc_integral = hist_mc->Integral();

  hist_mc->Scale(data_integral/mc_integral);
  hist_mc->SetLineColor(kBlue);
  //  hist_data->Draw("hist");
  hist_mc->Draw("hist same");


  c->SaveAs(output_file.c_str());
  return 0 ;

}
