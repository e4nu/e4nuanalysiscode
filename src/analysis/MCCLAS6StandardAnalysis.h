/**
 * This class deals with the analysis of GENIE MC data
 * 
 * It will classify events as signal or background, provided a configuration file where the Topology is defined
 * The background is stored for the background substraction correction
 * 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _MCANALYSISSTANDARD_I_H_
#define _MCANALYSISSTANDARD_I_H_

#include <iostream>
#include "analysis/MCCLAS6AnalysisI.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class MCCLAS6StandardAnalysis: virtual public MCCLAS6AnalysisI {
  public : 
    virtual ~MCCLAS6StandardAnalysis();

  protected : 
    MCCLAS6StandardAnalysis();
 
    EventI * GetValidEvent( const unsigned int event_id ) ;

  };
}

#endif
