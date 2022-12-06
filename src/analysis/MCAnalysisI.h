/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _MCANALYSIS_I_H_
#define _MCANALYSIS_I_H_

#include <iostream>
#include "analysis/AnalysisI.h"
#include "physics/MCEventHolder.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class MCAnalysisI : virtual public AnalysisI {
  public : 
    MCAnalysisI(); 

    bool LoadData( const std::string file ) ; 
    bool LoadData( const std::string file, const unsigned int nmax ) ; 
    void GetEvent( const unsigned int event_id ) const ;

    // Load Data from root file:
    virtual ~MCAnalysisI();

  private :

    MCEventHolder * fData ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
