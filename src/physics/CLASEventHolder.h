/**
 * This class contains all the information of a release
 * Implementation of EventHolderI for CLAS Data - GENIE 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLAS_EVENT_HOLDER_H_
#define _CLAS_EVENT_HOLDER_H_

#include "physics/EventHolderI.h"
#include "physics/CLASEvent.h"
#include "physics/EventI.h"


namespace e4nu {
  class CLASEventHolder : EventHolderI {
  public: 

    CLASEventHolder(); 
    CLASEventHolder( const std::string root_file, const int first_event, const int maxevents ) ; 
    CLASEventHolder( const std::vector<std::string> root_file_list ) ; 
    
    bool LoadBranch(void) ;
    bool LoadAllEvents(void) ;
    bool LoadEvent( const unsigned int event_id ) ;
    unsigned int GetNEventsChain(void) ;
    
    ~CLASEventHolder();

  private : 
    void Initialize(void) ;
    void Clear(void) ; 
    bool InitChain(void) ;

    // Members for root file


  } ;
}

#endif
