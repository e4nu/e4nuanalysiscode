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
    bool SubtractBackground( void ) ;
    bool Finalise(void) ; 
    bool StoreTree(CLAS6Event * event);

  private :

    e4nu::EventI * GetEvent( const unsigned int event_id ) ;

    CLAS6EventHolder * fData = nullptr ; 

    // Store Statistics after cuts
    long int kNEventsBeforeCuts = 0 ; 
    long int kNEventsAfterTopologyCut = 0 ; 
    long int kNBkgEvents = 0 ; 

    // Event Holder for signal and background
    std::map<int,std::vector<e4nu::CLAS6Event>> kAnalysedEventHolder;

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
