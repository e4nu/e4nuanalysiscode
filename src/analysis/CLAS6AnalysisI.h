/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLAS6ANALYSIS_I_H_
#define _CLAS6ANALYSIS_I_H_

#include <iostream>
#include "analysis/AnalysisI.h"
#include "physics/CLAS6EventHolder.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class CLAS6AnalysisI : virtual public AnalysisI {
  public : 
    CLAS6AnalysisI(); 

    bool LoadData(void) ; 
    e4nu::EventI * GetEvent( const unsigned int event_id ) ;
    unsigned int GetNEvents( void ) const ;

    // Load Data from root file:
    virtual ~CLAS6AnalysisI();

  private :

    CLAS6EventHolder * fData ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
