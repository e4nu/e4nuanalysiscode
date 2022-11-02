/**
 * This file contains utils specific for detector
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _DETECTOR_UTILS_H_
#define _DETECTOR_UTILS_H_

namespace e4nu {
  namespace utils
  {
    double GetSector( double phi ) ;
    bool IsValidSector( const double phi, const double EBeam, const bool use_all ) ;
  }
}

#endif
