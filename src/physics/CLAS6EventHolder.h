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

    // root file members
    // Declaration of leaf types
    UChar_t npart;
    UInt_t  runnb;
    UInt_t  evntid;
    UChar_t evstat;
    Char_t  evntype;
    Int_t   evntclas;
    Float_t q_l;
    Float_t t_l;
    Float_t tr_time;
    Float_t rf_time;
    Int_t   l2bit;
    Int_t   l3bit;
    UChar_t helicity;
    Int_t   hlsc;
    Int_t   intt;
    Int_t   helicity_cor;
    Int_t   gpart;  
    Int_t   id[40]; 
    Int_t   stat[40];
    Int_t   dc[40];
    Int_t   cc[40];
    Int_t   sc[40];
    Int_t   ec[40];
    Int_t   lec[40];
    Int_t   st[40];
    Float_t p[40];
    Float_t m[40];
    Int_t   q[40];
    Float_t b[40];
    Float_t cx[40];
    Float_t cy[40];
    Float_t cz[40];
    Float_t vx[40];
    Float_t vy[40];
    Float_t vz[40];
    Int_t   dc_part;
    Int_t   dc_sect[40];
    Int_t   dc_trk[40];
    Int_t   dc_stat[40];
    Int_t   tb_st[40];
    Float_t dc_xsc[40];
    Float_t dc_ysc[40];
    Float_t dc_zsc[40];
    Float_t dc_cxsc[40];
    Float_t dc_cysc[40];
    Float_t dc_czsc[40];
    Float_t dc_vx[40];
    Float_t dc_vy[40];
    Float_t dc_vz[40];
    Float_t dc_vr[40];
    Float_t tl1_cx[40];
    Float_t tl1_cy[40];
    Float_t tl1_cz[40];
    Float_t tl1_x[40];
    Float_t tl1_y[40];
    Float_t tl1_z[40];
    Float_t tl1_r[40];
    Float_t dc_c2[40];
    Int_t   ec_part;
    Int_t   ec_stat[40];
    Int_t   ec_sect[40];
    Int_t   ec_whol[40];
    Int_t   ec_inst[40];
    Int_t   ec_oust[40];
    Float_t etot[40];
    Float_t ec_ei[40];
    Float_t ec_eo[40];
    Float_t ec_t[40];
    Float_t ec_r[40];
    Float_t ech_x[40];
    Float_t ech_y[40];
    Float_t ech_z[40];
    Float_t ec_m2[40];
    Float_t ec_m3[40];
    Float_t ec_m4[40];
    Float_t ec_c2[40];
    Int_t   sc_part;
    Int_t   sc_sect[40];
    Int_t   sc_hit[40];
    Int_t   sc_pd[40];
    Int_t   sc_stat[40];
    Float_t edep[40];
    Float_t sc_t[40];
    Float_t sc_r[40];
    Float_t sc_c2[40];
    Int_t   cc_part;
    Int_t   cc_sect[40];
    Int_t   cc_hit[40];
    Int_t   cc_segm[40];
    Int_t   nphe[40]; 
    Float_t cc_t[40];
    Float_t cc_r[40]; 
    Float_t cc_c2[40]; 
    Int_t   lac_part;
    Int_t   lec_sect[40]; 
    Int_t   lec_hit[40];  
    Int_t   lec_stat[40];  
    Float_t lec_etot[40];  
    Float_t lec_ein[40];  
    Float_t lec_t[40];  
    Float_t lec_r[40];  
    Float_t lec_x[40];  
    Float_t lec_y[40];  
    Float_t lec_z[40];  
    Float_t lec_c2[40];  

    // List of branches
    TBranch * b_npart;   
    TBranch * b_runnb;   
    TBranch * b_evntid;   
    TBranch * b_evstat;   
    TBranch * b_evntype;   
    TBranch * b_evntclas;   
    TBranch * b_q_l;   
    TBranch * b_t_l;   
    TBranch * b_tr_time;   
    TBranch * b_rf_time;   
    TBranch * b_l2bit;   
    TBranch * b_l3bit;   
    TBranch * b_helicity;   
    TBranch * b_hlsc;   
    TBranch * b_intt;   
    TBranch * b_helicity_cor;   
    TBranch * b_gpart;   
    TBranch * b_id;   
    TBranch * b_stat;   
    TBranch * b_dc;   
    TBranch * b_cc;   
    TBranch * b_sc;   
    TBranch * b_ec;   
    TBranch * b_lec;   
    TBranch * b_st;   
    TBranch * b_p;   
    TBranch * b_m;   
    TBranch * b_q;   
    TBranch * b_b;   
    TBranch * b_cx;   
    TBranch * b_cy;   
    TBranch * b_cz;   
    TBranch * b_vx;   
    TBranch * b_vy;   
    TBranch * b_vz;   
    TBranch * b_dc_part;   
    TBranch * b_dc_sect;   
    TBranch * b_dc_trk;   
    TBranch * b_dc_stat;   
    TBranch * b_tb_st;   
    TBranch * b_dc_xsc;   
    TBranch * b_dc_ysc;   
    TBranch * b_dc_zsc;   
    TBranch * b_dc_cxsc;   
    TBranch * b_dc_cysc;   
    TBranch * b_dc_czsc;   
    TBranch * b_dc_vx;   
    TBranch * b_dc_vy;   
    TBranch * b_dc_vz;   
    TBranch * b_dc_vr;   
    TBranch * b_tl1_cx;   
    TBranch * b_tl1_cy;   
    TBranch * b_tl1_cz;   
    TBranch * b_tl1_x;   
    TBranch * b_tl1_y;   
    TBranch * b_tl1_z;   
    TBranch * b_tl1_r;   
    TBranch * b_dc_c2;   
    TBranch * b_ec_part;   
    TBranch * b_ec_stat;   
    TBranch * b_ec_sect;   
    TBranch * b_ec_whol;   
    TBranch * b_ec_inst;   
    TBranch * b_ec_oust;   
    TBranch * b_etot;   
    TBranch * b_ec_ei;   
    TBranch * b_ec_eo;   
    TBranch * b_ec_t;   
    TBranch * b_ec_r;   
    TBranch * b_ech_x;   
    TBranch * b_ech_y;   
    TBranch * b_ech_z;   
    TBranch * b_ec_m2;   
    TBranch * b_ec_m3;   
    TBranch * b_ec_m4;   
    TBranch * b_ec_c2;   
    TBranch * b_sc_part;   
    TBranch * b_sc_sect;   
    TBranch * b_sc_hit;   
    TBranch * b_sc_pd;   
    TBranch * b_sc_stat;   
    TBranch * b_edep;   
    TBranch * b_sc_t;   
    TBranch * b_sc_r;   
    TBranch * b_sc_c2;   
    TBranch * b_cc_part;   
    TBranch * b_cc_sect;   
    TBranch * b_cc_hit;   
    TBranch * b_cc_segm;   
    TBranch * b_nphe;   
    TBranch * b_cc_t;   
    TBranch * b_cc_r;   
    TBranch * b_cc_c2;   
    TBranch * b_lac_part;   
    TBranch * b_lec_sect;   
    TBranch * b_lec_hit;   
    TBranch * b_lec_stat;   
    TBranch * b_lec_etot;   
    TBranch * b_lec_ein;   
    TBranch * b_lec_t;   
    TBranch * b_lec_r;   
    TBranch * b_lec_x;   
    TBranch * b_lec_y;   
    TBranch * b_lec_z;   
    TBranch * b_lec_c2;   

  } ;
}

#endif
