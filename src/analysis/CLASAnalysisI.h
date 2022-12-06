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

    bool LoadData( const std::string file ) ; 
    bool LoadData( const std::string file, const unsigned int nmax ) ; 
    void GetEvent( const unsigned int event_id ) const ;

    // Load Data from root file:
    virtual ~CLASAnalysisI();

  private :

    CLASEventHolder * fData ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
