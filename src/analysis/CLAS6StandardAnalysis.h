/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLAS6STANDARDANALYSIS_I_H_
#define _CLAS6STANDARDANALYSIS_I_H_

#include <iostream>
#include "analysis/CLAS6AnalysisI.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class CLAS6StandardAnalysis : virtual public CLAS6AnalysisI {
  public : 
    CLAS6StandardAnalysis(); 

  protected : 
    virtual ~CLAS6StandardAnalysis();

    Event * GetValidEvent( const unsigned int event_id ) ;

  };
}

#endif
