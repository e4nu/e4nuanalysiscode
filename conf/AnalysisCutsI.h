/**
 * This file contains detector/analysis specific configurables
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _ANALYSIS_I_H_
#define _ANALYSIS_I_H_

namespace e4nu { 
  namespace conf {

    // Phase space limits
    const double kQ2Max = 1.75 ; 
    const double kQ2Min = 0.4 ; 
    const double kMinTheta = 10 ;
    const double kMaxTheta = 60 ;
    const double kMinEePrime = 0.5 ; 
    const double kMaxEePrime = 2.5 ;
      
    double GetMinMomentumCut( const int particle_pdg, const double EBeam ) ; 
    bool GetQ2Cut( double & Q2cut, const double Ebeam, const bool apply_Q2cut = true ) ;
    double GetWCut( const double Ebeam ) ;
  }
}

#endif
