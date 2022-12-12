/**
 * This file contains utils specific for detector
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _DETECTOR_UTILS_H_
#define _DETECTOR_UTILS_H_

#include "TLorentzVector.h"

namespace e4nu {
  namespace utils
  {
    double GetAcceptanceMapWeight( const int pdg, const TLorentzVector p4mom, const int target, const double EBeam, const std::string local_path ) ; 
    unsigned int GetSector( double phi ) ;
    bool IsValidSector( const double phi, const double EBeam, const bool use_all ) ;
  }
}

#endif
