/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLAS6ANALYSIS_I_H_
#define _CLAS6ANALYSIS_I_H_

#include <iostream>
#include "analysis/AnalysisI.h"
#include "physics/CLAS6EventHolder.h"
#include "physics/CLAS6Event.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class CLAS6AnalysisI : virtual public AnalysisI {
  public : 
    CLAS6AnalysisI(); 

  protected : 
    virtual ~CLAS6AnalysisI();

    bool LoadData(void);
    unsigned int GetNEvents( void ) const ;
    EventI * GetValidEvent( const unsigned int event_id ) ;
    e4nu::EventI * GetEvent( const unsigned int event_id ) ;
    bool Finalise( std::map<int,std::vector<e4nu::EventI*>> & event_holder ) ; 
    bool StoreTree(CLAS6Event * event);

  private :

    CLAS6EventHolder * fData = nullptr ; 

    // Store Statistics after cuts
    long int kNEventsBeforeCuts = 0 ; 
    long int kNEventsAfterTopologyCut = 0 ; 
    long int kNBkgEvents = 0 ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
