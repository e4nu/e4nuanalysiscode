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
    
    e4nu::Event * GetEvent(const unsigned int event_id) ;

    ~CLAS6EventHolder();

  private : 
    void Initialize(void) ;
    void Clear(void) ; 
    
    /* NOTICE:
     The CLAS event format is different from the GENIE MC format
     The code GetCharge_FilterData.C applies some generic CLAS-specific cuts
     and stores the information in a genie-like format. 
    */

    // Members for root file
    Int_t iev = 0 ;
    Int_t tgt = 0 ;
    Double_t Ev = 0 ;
    Double_t pxv = 0 ;
    Double_t pyv = 0 ;
    Double_t pzv = 0 ;
    Double_t El = 0 ;
    Double_t pxl = 0 ;
    Double_t pyl  = 0 ;
    Double_t pzl = 0 ;
    Int_t  nfp = 0 ;
    Int_t  nfn = 0 ;
    Int_t  nfpip = 0 ;
    Int_t  nfpim = 0 ;
    Int_t  nfpi0 = 0 ;
    Int_t nf = 0 ;
    Int_t pdgf[120] ;
    Double_t Ef[120] ;
    Double_t pxf[120] ;
    Double_t pyf[120] ;
    Double_t pzf[120] ;
    Double_t vtxx = 0 ;
    Double_t vtxy = 0 ;
    Double_t vtxz = 0 ;
    Double_t vtxt = 0 ;

    TBranch * b_iev = nullptr ;   
    TBranch * b_tgt = nullptr ;   
    TBranch * b_Ev = nullptr ;   
    TBranch * b_pxv = nullptr ;   
    TBranch * b_pyv = nullptr ;   
    TBranch * b_pzv = nullptr ;   
    TBranch * b_El = nullptr ;   
    TBranch * b_pxl = nullptr ;   
    TBranch * b_pyl = nullptr ;   
    TBranch * b_pzl = nullptr ;   
    TBranch * b_nfp = nullptr ;   
    TBranch * b_nfn = nullptr ;   
    TBranch * b_nfpip = nullptr ;   
    TBranch * b_nfpim = nullptr ;   
    TBranch * b_nfpi0 = nullptr ;   
    TBranch * b_nf = nullptr ;   
    TBranch * b_pdgf = nullptr ;   
    TBranch * b_Ef = nullptr ;   
    TBranch * b_pxf = nullptr ;   
    TBranch * b_pyf = nullptr ;   
    TBranch * b_pzf = nullptr ;   
    TBranch * b_vtxx = nullptr ;   
    TBranch * b_vtxy = nullptr ;   
    TBranch * b_vtxz = nullptr ;   
    TBranch * b_vtxt = nullptr ;   

  } ;
}

#endif
