/**
 * This file contains utils specific for detector
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date October 2022                                                                                                                                                                                              
 **/
#include "utils/DetectorUtils.h"
#include "TMath.h"

using namespace e4nu;
using namespace e4nu::utils ; 

unsigned int GetSector( double phi ) {
  phi *= TMath::RadToDeg() ; 
  phi += 30 ; //Add 30 degree for plotting and photon phi cut
  if ( phi < 0 ) phi += 360 ; //Add 360 so that electron phi is between 0 and 360 degree

  return (unsigned int) phi / 60 ; 
}

bool IsValidSector( const double phi, const double EBeam, const bool use_all ) {
  if( use_all ) return true ; 
  return true ; 

  unsigned int sector = utils::GetSector( phi ) ; 
  if ( ( sector == 2 || sector == 4 ) && EBeam == 1.161 ) return false ; 
  else if ( ( sector == 2 || sector == 3 || sector == 4 ) && EBeam == 2.261 ) return false ; 
  return true ; 

}
