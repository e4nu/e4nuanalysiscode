/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _E4NUANALYSIS_H_
#define _E4NUANALYSIS_H_

#include <iostream>
#include "analysis/MCAnalysisI.h"
#include "analysis/CLASAnalysisI.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class E4NuAnalysis : public MCAnalysisI, public CLASAnalysisI { 
  public : 
    E4NuAnalysis(); 
    E4NuAnalysis( const std::string conf_file ) ;
    E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) ;
    ~E4NuAnalysis();

    bool LoadData( const std::string file ) ; 
    bool LoadData( const std::string file, const unsigned int nmax ) ; 
    bool Finalise( const std::string out_file ) ;

    // Main Analyse function
    bool Analyse(void) ; 

  private : 
    e4nu::EventI * GetValidEvent( const unsigned int event_id ) ;
    unsigned int GetNEvents( void ) const ;

    long int fNEventsAfterEMomCut = 0 ; 
    long int fNEventsAfterPhiCut = 0 ; 
    long int fNEventsAfterQ2Cut = 0 ; 
    long int fNEventsAfterWCut = 0 ; 
    long int fNEventsAfterPhiOpeningAngleCut = 0 ; 
    long int fNEventsAfterThetaCut = 0 ;     
    
  };
}

#endif
