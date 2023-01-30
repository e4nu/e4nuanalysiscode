/**
 * \info These parameters are configurable 
 * the default values are set here
 **/

#ifndef _ANALYSIS_I_H_
#define _ANALYSIS_I_H_

#include <vector>
#include <map>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "physics/EventI.h"
#include "utils/Fiducial.h"
#include "conf/FiducialCutI.h"
#include "analysis/BackgroundI.h"

namespace e4nu { 

  class AnalysisI : public BackgroundI {

  public: 
    // Default constructor
    AnalysisI() ;
    AnalysisI( const std::string input_file ) ;
    AnalysisI( const double EBeam, const unsigned int TargetPdg ) ;

    bool Analyse( EventI * event ) ; 
    double GetElectronMinTheta( TLorentzVector emom ) ;
    bool Finalise(void) const ; 

    virtual ~AnalysisI();
      
  protected: 

    TF1 * fElectronFit ; 

    // Cuts counters
    long int fNEventsAfterEMomCut = 0 ; 
    long int fNEventsAfterEThetaCut = 0 ; 
    long int fNEventsAfterPhiCut = 0 ; 
    long int fNEventsAfterQ2Cut = 0 ; 
    long int fNEventsAfterWCut = 0 ; 
    long int fNEventsAfterPhiOpeningAngleCut = 0 ; 
    long int fNEventsAfterThetaCut = 0 ;     

  };
}

#endif
