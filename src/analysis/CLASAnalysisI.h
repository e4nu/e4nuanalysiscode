/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLASANALYSIS_I_H_
#define _CLASANALYSIS_I_H_

#include <iostream>
#include "analysis/ConfigureI.h"
#include "physics/CLASEventHolder.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class CLASAnalysisI : virtual public ConfigureI {
  public : 
    CLASAnalysisI(); 

    bool LoadData( const std::string file ) ; 
    bool LoadData( const std::string file, const unsigned int nmax ) ; 
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
