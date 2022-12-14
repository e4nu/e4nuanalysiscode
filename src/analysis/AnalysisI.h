/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _ANALYSIS_I_H_
#define _ANALYSIS_I_H_

#include <iostream>
#include "conf/ConfigureI.h"
#include "physics/EventI.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class AnalysisI : public ConfigureI {
  public : 

    // Store events on root file for further analysis
    bool StoreAnalysis( const std::string out_file ) { ; }

    virtual ~AnalysisI();

  protected : 
    AnalysisI(); 
    AnalysisI( const std::string conf_file ) ;
    AnalysisI( const double EBeam, const unsigned int TargetPdg ) ; 

    bool fIsDataLoaded ;

  private :
    void Initialize(void) ;

  };
}

#endif
