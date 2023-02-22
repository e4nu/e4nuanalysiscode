/**
 * This class deals with common analysis features which are present in data and MC
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
#include "analysis/BackgroundI.h"

namespace e4nu { 

  class AnalysisI : public BackgroundI {

  public: 
    // Default constructor
    AnalysisI() ;
    AnalysisI( const std::string input_file ) ;
    AnalysisI( const double EBeam, const unsigned int TargetPdg ) ;

    bool Analyse( EventI * event ) ; 
    bool Finalise(void) const ; 

    virtual ~AnalysisI();
      
  protected: 

    // Cuts counters
    long int kNEventsAfterEMomCut = 0 ; 
    long int kNEventsAfterEThetaCut = 0 ; 
    long int kNEventsAfterPhiCut = 0 ; 
    long int kNEventsAfterQ2Cut = 0 ; 
    long int kNEventsAfterWCut = 0 ; 
    long int kNEventsAfterPhiOpeningAngleCut = 0 ; 
    long int kNEventsAfterThetaCut = 0 ;     

  };
}

#endif
