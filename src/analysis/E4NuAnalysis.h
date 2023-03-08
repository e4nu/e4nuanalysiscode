/**
 * This class  is the interface for any analysis
 * Summary: This class deals with a generic E4nu analysis. 
 * It will treat the data according to it's configuration as MC or CLAS data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _E4NUANALYSIS_H_
#define _E4NUANALYSIS_H_

#include <iostream>
#include "TF1.h"
#include "analysis/MCAnalysisI.h"
#include "analysis/CLASAnalysisI.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class E4NuAnalysis : public MCAnalysisI, public CLASAnalysisI { 
  public : 
    E4NuAnalysis(); 
    E4NuAnalysis( const std::string conf_file ) ;
    E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) ;

    bool LoadData(void) ; 
    bool Analyse(void) ; 
    bool SubtractBackground( void ) ;
    bool Finalise(void);

    virtual ~E4NuAnalysis();

  private : 

    e4nu::EventI * GetValidEvent( const unsigned int event_id ) ;
    unsigned int GetNEvents( void ) const ;

    void Initialize(void) ; 
    
  };
}

#endif
