/**
 * This class contains all the information of a release
 * It is an interface - therefore it is independent of the release type, which can be either MC or data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _EVENT_HOLDER_I_H_
#define _EVENT_HOLDER_I_H_

#include <iostream>
#include <string>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//#include "physics/MCEvent.h"
#include "physics/EventI.h"

namespace e4nu {
  class EventHolderI {
  public : 
    virtual ~EventHolderI();

    unsigned int GetNEvents(void) const { return fMaxEvents ; } 

  protected : 
    EventHolderI(); 
    EventHolderI( const std::string root_file ) ; 
    EventHolderI( const std::string root_file, const int nmaxevents ) ; 
    EventHolderI( const std::vector<std::string> root_file_list ) ; 
    
    bool LoadMembers( const std::string file ) ; // returns tree number in TChain

    virtual bool LoadBranch(void) = 0 ; 
    virtual e4nu::EventI * GetEvent(const unsigned int event_id) = 0 ;

    std::unique_ptr<TChain> fEventHolderChain ; // Can contain more than one tree
    //TChain * fEventHolderChain ; 

    // Members
    bool fIsConfigured ; 
    int fMaxEvents ; 

  private :

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
