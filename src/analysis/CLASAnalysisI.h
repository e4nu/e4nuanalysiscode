/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLASANALYSIS_I_H_
#define _CLASANALYSIS_I_H_

#include <iostream>
#include "analysis/AnalysisI.h"
#include "physics/CLASEventHolder.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class CLASAnalysisI : virtual public AnalysisI {
  public : 
    CLASAnalysisI(); 

    bool LoadData(void) ; 
    e4nu::EventI * GetEvent( const unsigned int event_id ) ;
    unsigned int GetNEvents( void ) const ;

    // Load Data from root file:
    virtual ~CLASAnalysisI();

  private :

    CLASEventHolder * fData ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
