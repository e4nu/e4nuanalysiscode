/**
 * This file contains detector/analysis specific configurables
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/
#include "conf/ParticleI.h"

#ifndef _ANALYSIS_I_H_
#define _ANALYSIS_I_H_

namespace e4nu { 
  namespace AnalysisI {

    // Phase space limits
    const double kQ2Max = 1.75 ; 
    const double kQ2Min = 0.4 ; 
    const double kMinTheta = 10 ;
    const double kMaxTheta = 60 ;
    const double kMinEePrime = 0.5 ; 
    const double kMaxEePrime = 2.5 ;
      
    //    if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.) { minQ2 = 0.1; maxQ2 = 0.82; }
    // if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) { minQ2 = 1.1; maxQ2 = 4.4; }
    const double GetMinMomentumCut( const int particle_pdg, const double EBeam ) { 
      min_p = 0 ;
      if( particle_pdg == ParticleI::kPdgElectron ) {
	if( E == 1.161 /*GeV*/ ) min_p = 0.4 ; 
	else if ( E == 0.55 /*GeV*/ ) min_p = 0.55 ;
	else if ( E == 4.461 /*GeV*/ ) min_p = 1.1 ; 
      }
      return min_p ; 
    }

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

    const double GetWCut( const double Ebeam ) {
      if ( Ebeam < 2 /*GeV*/ ) return 2. ; 
      return 0. ;
    }
  }
}

#endif _DETECTOR_I_H_
