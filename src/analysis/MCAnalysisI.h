/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _MCANALYSIS_I_H_
#define _MCANALYSIS_I_H_

#include <iostream>
#include <map>
#include "utils/Fiducial.h"
#include "analysis/AnalysisI.h"
#include "physics/MCEventHolder.h"
#include "physics/MCEvent.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class MCAnalysisI : virtual public AnalysisI {
  public : 
    virtual ~MCAnalysisI();

  protected :
 
    MCAnalysisI();
 
    bool LoadData( const std::string file ) ; 
    bool LoadData( const std::string file, const unsigned int nmax ) ; 
    virtual unsigned int GetNEvents( void ) const ;
    e4nu::EventI * GetValidEvent( const unsigned int event_id ) ;
    void SmearParticles( MCEvent * event ) ;
    bool Finalise( const std::string out_file ) ;

  private :

    e4nu::EventI * GetEvent( const unsigned int event_id ) ;

    MCEventHolder * fData ; 
    Fiducial * kFiducialCut ;


    // Store Statistics after cuts
    unsigned int fEventsBeforeCuts = 0 ; 
    unsigned int fNEventsAfterFiducial = 0 ; 
    unsigned int fNEventsAfterEMomCut = 0 ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
