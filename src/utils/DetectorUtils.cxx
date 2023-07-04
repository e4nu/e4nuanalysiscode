/**
 * This file contains utils specific for detector
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date October 2022                                                                                                                                                                                              
 **/
#include <iostream>
#include "utils/DetectorUtils.h"
#include "conf/AccpetanceMapsI.h"
#include "TMath.h"
#include "TFile.h"

using namespace e4nu;

double utils::GetAcceptanceMapWeight( TH3D & acc, TH3D & gen, const TLorentzVector p4mom ){

  double phi = p4mom.Phi() + TMath::Pi() ; 
  if(phi > (2*TMath::Pi() - TMath::Pi()/6.) ) { phi -= 2*TMath::Pi(); }
  phi *= 180/TMath::Pi() ;

  //Redefinition of the phi angle
  //because the acceptance maps are defined between (-30,330)
  //  phi -= 30 ; 

  double pbin_gen = gen.GetXaxis()->FindBin(p4mom.P());
  double tbin_gen = gen.GetYaxis()->FindBin(p4mom.CosTheta());
  double phibin_gen = gen.GetZaxis()->FindBin(phi);
  double num_gen = gen.GetBinContent(pbin_gen, tbin_gen, phibin_gen);

  double pbin_acc = acc.GetXaxis()->FindBin(p4mom.P());
  double tbin_acc = acc.GetYaxis()->FindBin(p4mom.CosTheta());
  double phibin_acc = acc.GetZaxis()->FindBin(phi);
  double num_acc = acc.GetBinContent(pbin_acc, tbin_acc, phibin_acc);

  double acc_w = num_acc / num_gen ; 
  if( fabs(acc_w) != acc_w ) return 0 ; // Safety check
  return acc_w ; 
} 

unsigned int utils::GetSector( double phi ) {
  phi *= TMath::RadToDeg() ; 
  phi += 30 ; //Add 30 degree for plotting and photon phi cut
  if ( phi < 0 ) phi += 360 ; //Add 360 so that electron phi is between 0 and 360 degree

  return (unsigned int) phi / 60. ; 
}

bool utils::IsValidSector( const double phi, const double EBeam, const bool use_all ) {
  if( use_all ) return true ; 

  unsigned int sector = utils::GetSector( phi ) ; 
  if ( ( sector == 2 || sector == 4 ) && EBeam == 1.161 ) return false ; 
  else if ( ( sector == 2 || sector == 3 || sector == 4 ) && EBeam == 2.261 ) return false ; 
  return true ; 

}
