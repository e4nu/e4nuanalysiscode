/**                                     
 * \info This script contains general information on CLAS6 constants
 **/

#ifndef _CLAS6CONST_I_H_
#define _CLAS6CONST_I_H_

#include "conf/TargetI.h"

namespace e4nu {
  namespace conf {
    double GetIntegratedCharge( const unsigned int tgt_pdg, const double EBeam ) ; 
    double GetIntegratedChargeFilterRuns( const unsigned int tgt_pdg, const double EBeam ); 
    double GetIntegratedChargeNewFilterRuns( const unsigned int tgt_pdg, const double EBeam ); 
    double GetIntegratedChargeGoodRunListAllRuns( const unsigned int tgt_pdg, const double EBeam ); 
    double GetIntegratedChargeGoodRunListLowRuns( const unsigned int tgt_pdg, const double EBeam ); 
    double GetIntegratedChargeGoodRunListHighRuns( const unsigned int tgt_pdg, const double EBeam ); 
    double GetTargetLength( const unsigned int tgt_pdg ) ;
    double GetTargetDensity( const unsigned int tgt_pdg ) ;
  }
}
#endif
