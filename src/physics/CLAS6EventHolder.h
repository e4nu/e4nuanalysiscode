/**
 * This class contains all the information of a release
 * Implementation of EventHolderI for CLAS6 Data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _CLAS6_EVENT_HOLDER_H_
#define _CLAS6_EVENT_HOLDER_H_

#include "physics/EventHolderI.h"

namespace e4nu {
  class CLAS6EventHolder : public EventHolderI {
  public: 

    CLAS6EventHolder(); 
    CLAS6EventHolder( const std::string root_file, const unsigned int first_event, const unsigned int maxevents ) ; 
    CLAS6EventHolder( const std::vector<std::string> root_file_list ) ; // add first and last 
    
    bool LoadBranch(void) ;
    
    e4nu::EventI * GetEvent(const unsigned int event_id) ;

    ~CLAS6EventHolder();

  private : 
    void Initialize(void) ;
    void Clear(void) ; 
    bool InitChain(void) ;

  } ;
}

#endif
