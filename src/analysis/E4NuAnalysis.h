/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _E4NUANALYSIS_H_
#define _E4NUANALYSIS_H_

#include <iostream>
#include "TF1.h"
#include "analysis/MCAnalysisI.h"
#include "analysis/CLASAnalysisI.h"
#include "utils/Subtraction.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class E4NuAnalysis : public MCAnalysisI, public CLASAnalysisI { 
  public : 
    E4NuAnalysis(); 
    E4NuAnalysis( const std::string conf_file ) ;
    E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) ;
    virtual ~E4NuAnalysis();

    bool LoadData(void) ; 

    // Main Analyse function
    bool Analyse(void) ; 
    bool SubstractBackground(void) ; 
    bool Finalise(void);

  private : 
    e4nu::EventI * GetValidEvent( const unsigned int event_id ) ;
    unsigned int GetNEvents( void ) const ;

    Subtraction * fRotation;

    void Initialize(void) ; 
    
  };
}

#endif
