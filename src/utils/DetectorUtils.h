/**
 * This file contains utils specific for detector
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _DETECTOR_UTILS_H_
#define _DETECTOR_UTILS_H_

#include "TLorentzVector.h"
#include "TH3D.h"
#include "TFile.h"

namespace e4nu {
  namespace utils
  {
    double GetAcceptanceMapWeight( TH3D & h_acc, TH3D & h_gen, const TLorentzVector p4mom );
    unsigned int GetSector( double phi ) ;
    bool IsValidSector( const double phi, const double EBeam, const bool use_all ) ;
  }
}

#endif
