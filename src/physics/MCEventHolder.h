/**
 * This class contains all the information of a release
 * Implementation of EventHolderI for MC Data - GENIE 
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _MC_EVENT_HOLDER_H_
#define _MC_EVENT_HOLDER_H_

#include "physics/EventHolderI.h"

namespace e4nu {
  class MCEventHolder : EventHolderI {
  public: 
    MCEventHolder(); 
    MCEventHolder( const std::string root_file ) ; 
    MCEventHolder( const std::vector<std::string> root_file_list ) ; 
    
    bool LoadEvents(void) ;

    ~MCEventHolder();

  private : 
    void Initialize(void) ;
    void Clear(void); 
    bool InitChain(void);

    // Members for root file
    Int_t neu;
    Int_t fspl;
    Int_t tgt;
    Int_t Z;
    Int_t A;
    Int_t hitnuc;
    Int_t hitqrk;
    Int_t resid;
    Bool_t sea;
    Bool_t qel;
    Bool_t mec;
    Bool_t res;
    Bool_t dis;
    Bool_t coh;
    Bool_t dfr;
    Bool_t imd;
    Bool_t imdanh;
    Bool_t singlek;
    Bool_t nuel;
    Bool_t em;
    Bool_t cc;
    Bool_t nc;
    Bool_t charm;
    Int_t neut_code;
    Int_t nuance_code;
    Double_t wght;
    Double_t xs;
    Double_t ys;
    Double_t ts;
    Double_t Q2s;
    Double_t Ws;
    Double_t x;
    Double_t y;
    Double_t t;
    Double_t Q2;
    Double_t W;
    Double_t EvRF;
    Double_t Ev;
    Double_t pxv;
    Double_t pyv;
    Double_t pzv;
    Double_t En;
    Double_t pxn;
    Double_t pyn;
    Double_t pzn;
    Double_t El;
    Double_t pxl;
    Double_t pyl;
    Double_t pzl;
    Double_t pl;
    Double_t cthl;
    Int_t  nfp;
    Int_t  nfn;
    Int_t  nfpip;
    Int_t  nfpim;
    Int_t  nfpi0;
    Int_t  nfkp;
    Int_t  nfkm;
    Int_t  nfk0;
    Int_t  nfem;
    Int_t  nfother;
    Int_t  nip;
    Int_t  nin;
    Int_t  nipip;
    Int_t  nipim;
    Int_t  nipi0;
    Int_t  nikp;
    Int_t  nikm;
    Int_t  nik0;
    Int_t  niem;
    Int_t  niother;
    Int_t  ni;
    Int_t  pdgi[2];
    Int_t  resc[1];
    Double_t Ei[2];
    Double_t pxi[2];
    Double_t pyi[2];
    Double_t pzi[2];
    Int_t nf;
    Int_t pdgf[120];
    Double_t Ef[120];
    Double_t pxf[120];
    Double_t pyf[120];
    Double_t pzf[120];
    Double_t pf[120]; 
    Double_t cthf[120];
    Double_t vtxx;
    Double_t vtxy;
    Double_t vtxz;
    Double_t vtxt;
    Double_t sumKEf;
    Double_t calresp0;

    TBranch * b_iev;   
    TBranch * b_neu;   
    TBranch * b_fspl;   
    TBranch * b_tgt;   
    TBranch * b_Z;   
    TBranch * b_A;   
    TBranch * b_hitnuc;   
    TBranch * b_hitqrk;   
    TBranch * b_resid;   
    TBranch * b_sea;   
    TBranch * b_qel;   
    TBranch * b_mec;   
    TBranch * b_res;   
    TBranch * b_dis;   
    TBranch * b_coh;   
    TBranch * b_dfr;   
    TBranch * b_imd;   
    TBranch * b_imdanh;   
    TBranch * b_singlek;   
    TBranch * b_nuel;   
    TBranch * b_em;   
    TBranch * b_cc;   
    TBranch * b_nc;   
    TBranch * b_charm;   
    TBranch * b_neut_code;   
    TBranch * b_nuance_code;   
    TBranch * b_wght;   
    TBranch * b_xs;   
    TBranch * b_ys;   
    TBranch * b_ts;   
    TBranch * b_Q2s;   
    TBranch * b_Ws;   
    TBranch * b_x;   
    TBranch * b_y;   
    TBranch * b_t;   
    TBranch * b_Q2;   
    TBranch * b_W;   
    TBranch * b_EvRF;   
    TBranch * b_Ev;   
    TBranch * b_pxv;   
    TBranch * b_pyv;   
    TBranch * b_pzv;   
    TBranch * b_En;   
    TBranch * b_pxn;   
    TBranch * b_pyn;   
    TBranch * b_pzn;   
    TBranch * b_El;   
    TBranch * b_pxl;   
    TBranch * b_pyl;   
    TBranch * b_pzl;   
    TBranch * b_pl;   
    TBranch * b_cthl;   
    TBranch * b_nfp;   
    TBranch * b_nfn;   
    TBranch * b_nfpip;   
    TBranch * b_nfpim;   
    TBranch * b_nfpi0;   
    TBranch * b_nfkp;   
    TBranch * b_nfkm;   
    TBranch * b_nfk0;   
    TBranch * b_nfem;   
    TBranch * b_nfother;   
    TBranch * b_nip;   
    TBranch * b_nin;   
    TBranch * b_nipip;   
    TBranch * b_nipim;   
    TBranch * b_nipi0;   
    TBranch * b_nikp;   
    TBranch * b_nikm;   
    TBranch * b_nik0;   
    TBranch * b_niem;   
    TBranch * b_niother;   
    TBranch * b_ni;   
    TBranch * b_pdgi;   
    TBranch * b_resc;   
    TBranch * b_Ei;   
    TBranch * b_pxi;   
    TBranch * b_pyi;   
    TBranch * b_pzi;   
    TBranch * b_nf;   
    TBranch * b_pdgf;   
    TBranch * b_Ef;   
    TBranch * b_pxf;   
    TBranch * b_pyf;   
    TBranch * b_pzf;   
    TBranch * b_pf;   
    TBranch * b_cthf;   
    TBranch * b_vtxx;   
    TBranch * b_vtxy;   
    TBranch * b_vtxz;   
    TBranch * b_vtxt;   
    TBranch * b_sumKEf;   
    TBranch * b_calresp0;   

  };
}

#endif
