/**
 * This class contains all the information of a release
 * Implementation of EventHolderI for MC Data - GENIE 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _MC_EVENT_HOLDER_H_
#define _MC_EVENT_HOLDER_H_

#include "physics/EventHolderI.h"

namespace e4nu {
  class MCEventHolder : public EventHolderI {
  public: 

    MCEventHolder(); 
    MCEventHolder( const std::string root_file ) ; 
    MCEventHolder( const std::string root_file, const int maxevents ) ; 
    MCEventHolder( const std::vector<std::string> root_file_list ) ; 
    
    bool LoadBranch(void) ;
    
    e4nu::EventI * GetEvent(const unsigned int event_id) ;

    ~MCEventHolder();

  private : 
    void Initialize(void) ;
    void Clear(void) ; 
    bool InitChain(void) ;

    // Members for root file
    Int_t iev = 0 ;
    Int_t tgt = 0 ;
    Bool_t qel = false ;
    Bool_t mec = false ;
    Bool_t res = false ;
    Bool_t dis = false ;
    Bool_t em = false ;
    Bool_t cc = false ;
    Bool_t nc = false ;
    Double_t wght = 0 ;
    Double_t xs = 0 ;
    Double_t ys = 0 ;
    Double_t ts = 0 ;
    Double_t Q2s = 0 ;
    Double_t Ws = 0 ;
    Double_t x = 0 ;
    Double_t y = 0 ;
    Double_t t = 0 ;
    Double_t Q2 = 0 ;
    Double_t W = 0 ;
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
    Int_t  nfkp = 0 ;
    Int_t  nfkm = 0 ;
    Int_t  nfk0 = 0 ;
    Int_t  nfem = 0 ;
    Int_t  nfother = 0 ;
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
    TBranch * b_qel = nullptr ;   
    TBranch * b_mec = nullptr ;   
    TBranch * b_res = nullptr ;   
    TBranch * b_dis = nullptr ;   
    TBranch * b_em = nullptr ;   
    TBranch * b_cc = nullptr ;   
    TBranch * b_nc = nullptr ;   
    TBranch * b_wght = nullptr ;   
    TBranch * b_xs = nullptr ;   
    TBranch * b_ys = nullptr ;   
    TBranch * b_ts = nullptr ;   
    TBranch * b_Q2s = nullptr ;   
    TBranch * b_Ws = nullptr ;   
    TBranch * b_x = nullptr ;   
    TBranch * b_y = nullptr ;   
    TBranch * b_t = nullptr ;   
    TBranch * b_Q2 = nullptr ;   
    TBranch * b_W = nullptr ;   
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
    TBranch * b_nfkp = nullptr ;   
    TBranch * b_nfkm = nullptr ;   
    TBranch * b_nfk0 = nullptr ;   
    TBranch * b_nfem = nullptr ;   
    TBranch * b_nfother = nullptr ;   
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
