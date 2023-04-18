/**                                                                                                                                                                                            * \info This script contains general information on ficucial cut parameters
 **/

#include <iostream>
#include "conf/CLAS6ConstantsI.h"
#include "conf/TargetI.h"

using namespace e4nu;

double conf::GetIntegratedCharge( const unsigned int tgt_pdg, const double EBeam ) {
  double ic = 0 ;// mC

  if( tgt_pdg == kPdgHe4 ) { 
    if( EBeam == 2.261 ) ic = 1.08 ;
    else if ( EBeam == 4.461 ) ic = 0.87 ;
  } else if ( tgt_pdg == kPdgC12 ) { 
    if( EBeam == 1.161 ) ic = 0.19 ;
    else if( EBeam == 2.261 ) ic = 1.79 ;
    else if ( EBeam == 4.461 ) ic = 2.14 ;
  } else if ( tgt_pdg == kPdgFe56 ) { 
    if( EBeam == 2.261 ) ic = 0.22 ;
    else if ( EBeam == 4.461 ) ic = 0.29 ;
  } 
  return ic ; 
}

double conf::GetIntegratedChargeNewFilterRuns( const unsigned int tgt_pdg, const double EBeam ){
  double ic = 0 ;// mC
  if ( tgt_pdg == kPdgC12 ) { 
    if( EBeam == 1.161 ) ic = 0.070707652 ;
  } 
  return ic ; 
}

double conf::GetIntegratedChargeGoodRunListAllRuns( const unsigned int tgt_pdg, const double EBeam ){
  double ic = 0 ;// mC
  if ( tgt_pdg == kPdgC12 ) { 
    if( EBeam == 1.161 ) ic = 0.18432 ;
  } 
  return ic ; 
}

double conf::GetIntegratedChargeGoodRunListLowRuns( const unsigned int tgt_pdg, const double EBeam ){
  double ic = 0 ;// mC
  if ( tgt_pdg == kPdgC12 ) { 
    if( EBeam == 1.161 ) ic = 0.079 ;
  } 
  return ic ; 
}

double conf::GetIntegratedChargeGoodRunListHighRuns( const unsigned int tgt_pdg, const double EBeam ){
  double ic = 0 ;// mC
  if ( tgt_pdg == kPdgC12 ) { 
    if( EBeam == 1.161 ) ic = 0.105 ;
  } 
  return ic ; 
}

double conf::GetTargetLength( const unsigned int tgt_pdg ) {
  double length = 0 ; // cm 

  if( tgt_pdg == kPdgHe4 ) {
    length = 4.3 ; 
  } else if ( tgt_pdg == kPdgC12 ) { 
    length = 0.1 ;
  } else if ( tgt_pdg == kPdgFe56 ) { 
    length = 0.015 ; 
  } 
  return length ; 
}


double conf::GetTargetDensity( const unsigned int tgt_pdg ) {
  double density = 0 ; // g/cm^2
  
  if( tgt_pdg == kPdgHe4 ) {
	density = 0.125 ; 
  } else if ( tgt_pdg == kPdgC12 ) { 
    density = 1.786 ;
  } else if ( tgt_pdg == kPdgFe56 ) { 
    density = 7.872 ; 
  } 
  return density ; 
}
