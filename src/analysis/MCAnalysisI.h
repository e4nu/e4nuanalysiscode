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
    bool Finalise( std::map<int,std::vector<e4nu::EventI*>> & event_holder ) ; 
    bool StoreTree(MCEvent * event);

    // ID for Background historams
    unsigned int id_signal, id_tottruebkg, id_totestbkg, id_acccorr;

  private :

    void SmearParticles( MCEvent * event ) ;
    void ApplyMomentumCut( MCEvent * event ) ;
    bool ApplyFiducialCut( MCEvent * event ) ; 
    void ApplyAcceptanceCorrection( MCEvent * event ) ;
    EventI * GetEvent( const unsigned int event_id ) ;
    void PlotBkgInformation( EventI * event ) ;
    
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
