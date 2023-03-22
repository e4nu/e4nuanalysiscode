/**
 * This class contains all the information of a release
 * Implementation of EventHolderI for CLAS6 Data - GENIE 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLAS6_EVENT_HOLDER_H_
#define _CLAS6_EVENT_HOLDER_H_

#include "physics/EventHolderI.h"
#include "physics/CLAS6Event.h"
#include "physics/EventI.h"


namespace e4nu {
  class CLAS6EventHolder : EventHolderI {
  public: 

    CLAS6EventHolder(); 
    CLAS6EventHolder( const std::string root_file, const int first_event, const int maxevents ) ; 
    CLAS6EventHolder( const std::vector<std::string> root_file_list ) ; 
    
    bool LoadBranch(void) ;
    bool LoadAllEvents(void) ;
    bool LoadEvent( const unsigned int event_id ) ;
    unsigned int GetNEventsChain(void) ;
    
    ~CLAS6EventHolder();

  private : 
    void Initialize(void) ;
    void Clear(void) ; 
    bool InitChain(void) ;

    // Members for root file


  } ;
}

#endif
