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

using namespace e4nu::conf ;

namespace e4nu {
  class MCCLAS6AnalysisI: virtual public AnalysisI {
  public :
    virtual ~MCCLAS6AnalysisI();

  protected :
    MCCLAS6AnalysisI();

    bool LoadData(void);
    unsigned int GetNEvents( void ) const ;
    Event * GetValidEvent( const unsigned int event_id ) ;
    bool Finalise( std::map<int,std::vector<e4nu::Event>> & event_holder ) ;
    bool StoreTree(Event event);

  private :

    void SmearParticles( Event & event ) ;
    void ApplyAcceptanceCorrection( Event & event ) ;
    Event * GetEvent( const unsigned int event_id ) ;

    MCEventHolder * fData = nullptr ;
    std::map<int,std::unique_ptr<TFile>> kAcceptanceMap;
    std::map<int,std::unique_ptr<TH3D>> kAccMap ;
    std::map<int,std::unique_ptr<TH3D>> kGenMap ;

    // XSec value
    double kXSec = 0 ;

    void Initialize(void) ;
    void Clear(void);

  };
}

#endif
