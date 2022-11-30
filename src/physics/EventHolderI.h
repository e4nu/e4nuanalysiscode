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

namespace e4nu {
  class EventHolderI {
  public: 
    EventHolderI(); 
    EventHolderI( const std::string root_file ) ; 
    EventHolderI( const std::vector<std::string> root_file_list ) ; 
    ~EventHolderI();

  protected : 
    bool LoadMembers( const std::string file ) ; // returns tree number in TChain
    
    TChain * fEventHolderChain ; // Can contain more than one tree

    // Declaration of leaf types
    Int_t           iev;
        
    // List of branches
    TBranch        *b_iev; 
    
  private : 
    void Initialize(void) ;
    void Clear(void); 
    bool InitChain(void);

    // Members
    bool fIsConfigured ; 
  };
}

#endif
