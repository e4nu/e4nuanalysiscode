/**
 * This file contains utils specific for detector
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date October 2022                                                                                                                                                                                              
 **/

#include "utils/DetectorUtils.h"
#include "conf/AccpetanceMapsI.h"
#include "TMath.h"
#include "TH3D.h"
#include "TFile.h"

using namespace e4nu;

double utils::GetAcceptanceMapWeight( const int pdg, const TLorentzVector p4mom, const int target, const double EBeam, const std::string local_path ) {
  TFile * file_acceptance = TFile::Open( conf::GetAcceptanceFile( pdg, target, EBeam, local_path ).c_str() ) ; 
  if( !file_acceptance ) return 1. ;
 
  TH3D * acc = (TH3D*) file_acceptance -> Get("Accepted Particles");
  TH3D * gen = (TH3D*) file_acceptance -> Get("Generated Particles");
  if( !acc || !gen ) return 1. ;

  double phi = p4mom.Phi() ; 
  //map 330 till 360 to [-30:0] for the acceptance map histogram
  if(phi > (2*TMath::Pi() - TMath::Pi()/6.) ) { phi -= 2*TMath::Pi(); }
  
  double pbin_gen = gen->GetXaxis()->FindBin(p4mom.P());
  double tbin_gen = gen->GetYaxis()->FindBin(p4mom.CosTheta());
  double phibin_gen = gen->GetZaxis()->FindBin(phi*180/TMath::Pi());
  double num_gen = gen->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

  double pbin_acc = acc->GetXaxis()->FindBin(p4mom.P());
  double tbin_acc = acc->GetYaxis()->FindBin(p4mom.CosTheta());
  double phibin_acc = acc->GetZaxis()->FindBin(phi*180/TMath::Pi());
  double num_acc = acc->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

  double acc_ratio = (double)num_acc / (double)num_gen;
  //  double acc_err = (double)sqrt(acc_ratio*(1-acc_ratio)) / (double)num_gen;

  return acc_ratio;
} 

unsigned int utils::GetSector( double phi ) {
  phi *= TMath::RadToDeg() ; 
  phi += 30 ; //Add 30 degree for plotting and photon phi cut
  if ( phi < 0 ) phi += 360 ; //Add 360 so that electron phi is between 0 and 360 degree

  return (unsigned int) phi / 60 ; 
}

bool utils::IsValidSector( const double phi, const double EBeam, const bool use_all ) {
  if( use_all ) return true ; 
  return true ; 

  unsigned int sector = utils::GetSector( phi ) ; 
  if ( ( sector == 2 || sector == 4 ) && EBeam == 1.161 ) return false ; 
  else if ( ( sector == 2 || sector == 3 || sector == 4 ) && EBeam == 2.261 ) return false ; 
  return true ; 

}
