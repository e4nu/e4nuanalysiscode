/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _MCANALYSIS_I_H_
#define _MCANALYSIS_I_H_

#include <iostream>
#include <map>
#include "TH3D.h"
#include "utils/Fiducial.h"
#include "analysis/ConfigureI.h"
#include "physics/MCEventHolder.h"
#include "physics/MCEvent.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class MCAnalysisI: virtual public ConfigureI {
  public : 
    virtual ~MCAnalysisI();

  protected :
 
    MCAnalysisI();
 
    bool LoadData(void);
    unsigned int GetNEvents( void ) const ;
    EventI * GetValidEvent( const unsigned int event_id ) ;
    void SmearParticles( MCEvent * event ) ;
    bool Finalise( ) ; 

  private :

    e4nu::EventI * GetEvent( const unsigned int event_id ) ;

    MCEventHolder * fData = nullptr ; 

    std::map<int,std::unique_ptr<TFile>> kAcceptanceMap;
    std::map<int,std::unique_ptr<TH3D>> kAccMap ; 
    std::map<int,std::unique_ptr<TH3D>> kGenMap ; 
    std::unique_ptr<Fiducial> kFiducialCut ;

    // Store Statistics after cuts
    long int fEventsBeforeCuts = 0 ; 
    long int fNEventsAfterFiducial = 0 ; 
    long int fNEventsAfterTopologyCut = 0 ; 
    long int fNBkgEvents = 0 ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
