/**
 * This file contains detector specific configurables
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _DETECTOR_I_H_
#define _DETECTOR_I_H_

namespace e4nu { 
  namespace DetectorI {
      bool GetQ2Cut( double & Q2cut, const double Ebeam, bool apply_Q2cut = true ) {
	if ( Ebeam >= 1. && Ebeam < 2. ) {
	  Q2cut = 0.1;
	  return true ;
	} else if ( Ebeam >= 2. && Ebeam < 3. ) {
	  Q2cut = 0.4 ;
	  return true ;
	} else if ( Ebeam >= 4. && Ebeam < 5. ) {
	  Q2cut= 0.8 ;
	  return true ;
	}
	return false ;
      }
    }
}

#endif _DETECTOR_I_H_
