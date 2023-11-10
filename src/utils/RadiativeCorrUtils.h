#ifndef _RADUTILS_H_
#define _RADUTILS_H_
#include "TLorentzVector.h"

namespace e4nu {
  namespace utils
  {
    // QEL
    double RadCorrQELVertex( const double Q2 ) ; 
    double RadCorrQELVacumm( const double Q2 ) ; 
    double RadCorrQELRealRad( const double Q2, const double E, const double Ep, const double theta) ; 
    double RadCorrQELRealProtonD1( const double Q2, const double E, const double Ep, const double theta) ; 
    double RadCorrQELRealProtonD20( const double Q2, const double E, const double Ep, const double theta, const double deltaE, const TLorentzVector nucleon ) ;
    double RadCorrQELRealProtonD21( const double Q2, const double E, const double Ep, const double theta) ; 

  }
}

#endif
