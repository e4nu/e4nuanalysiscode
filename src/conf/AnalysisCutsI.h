/**
 * This file contains detector/analysis specific configurables
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _ANALYSISCUTS_I_H_
#define _ANALYSISCUTS_I_H_

namespace e4nu { 
  namespace conf {

    double GetMinMomentumCut( const int particle_pdg, const double EBeam ) ; 
    bool ValidPhiOpeningAngle( double phi /*rad*/ ) ;
    bool GoodSectorPhiSlice( double phi /*rad*/ ) ; 
    bool GetQ2Cut( double & Q2cut, const double Ebeam ) ; 
    bool GetWCut( double & WCut, const double Ebeam ) ;
  }
}

#endif
