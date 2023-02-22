/**
 * This class deals with the analysis of GENIE MC data
 * 
 * It will classify events as signal or background, provided a configuration file where the Topology is defined
 * The background is stored for the background substraction correction
 * 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _MCANALYSIS_I_H_
#define _MCANALYSIS_I_H_

#include <iostream>
#include <map>
#include "TH3D.h"
#include "utils/Fiducial.h"
#include "analysis/AnalysisI.h"
#include "physics/MCEventHolder.h"
#include "physics/MCEvent.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class MCAnalysisI: virtual public AnalysisI {
  public : 

    virtual ~MCAnalysisI();

  protected :
 
    MCAnalysisI();
 
    bool LoadData(void);
    unsigned int GetNEvents( void ) const ;
    EventI * GetValidEvent( const unsigned int event_id ) ;
    bool SubtractBackground( void ) ;
    bool Finalise(void) ; 
    bool StoreTree(MCEvent * event);

  private :

    void SmearParticles( MCEvent * event ) ;
    e4nu::EventI * GetEvent( const unsigned int event_id ) ;

    MCEventHolder * fData = nullptr ; 
    std::map<int,std::unique_ptr<TFile>> kAcceptanceMap;
    std::map<int,std::unique_ptr<TH3D>> kAccMap ; 
    std::map<int,std::unique_ptr<TH3D>> kGenMap ; 

    // Store Statistics after cuts
    long int kNEventsBeforeCuts = 0 ; 
    long int kNEventsAfterFiducial = 0 ; 
    long int kNEventsAfterTopologyCut = 0 ; 
    long int kNBkgEvents = 0 ; 

    // XSec value
    double kXSec = 0 ; 

    // Event Holder for signal and background
    std::map<int,std::vector<e4nu::MCEvent>> kAnalysedEventHolder;

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
