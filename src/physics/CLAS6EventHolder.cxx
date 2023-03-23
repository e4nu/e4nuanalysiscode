// _______________________________________________
/*
 * CLAS6EventHolder implementation 
 */

#include "physics/CLAS6EventHolder.h"
#include "physics/CLAS6Event.h"

using namespace e4nu ; 

CLAS6EventHolder::CLAS6EventHolder(): EventHolderI() { this -> Initialize(); } 

CLAS6EventHolder::~CLAS6EventHolder() {
  //  this -> Clear();
}

CLAS6EventHolder::CLAS6EventHolder( const std::string file, const unsigned first_event, const unsigned int nmaxevents ): EventHolderI( file, first_event , nmaxevents ) { 
  this -> Initialize() ; 
}

CLAS6EventHolder::CLAS6EventHolder( const std::vector<std::string> files ): EventHolderI( files ) { 
  this -> Initialize() ; 
}


bool CLAS6EventHolder::LoadBranch(void) {
  if( !fEventHolderChain ) return false ; 


  fEventHolderChain -> SetBranchAddress("npart", &npart, &b_npart);
  fEventHolderChain -> SetBranchAddress("runnb", &runnb, &b_runnb);
  fEventHolderChain -> SetBranchAddress("evntid", &evntid, &b_evntid);
  fEventHolderChain -> SetBranchAddress("evstat", &evstat, &b_evstat);
  fEventHolderChain -> SetBranchAddress("evntype", &evntype, &b_evntype);
  fEventHolderChain -> SetBranchAddress("evntclas", &evntclas, &b_evntclas);
  fEventHolderChain -> SetBranchAddress("q_l", &q_l, &b_q_l);
  fEventHolderChain -> SetBranchAddress("t_l", &t_l, &b_t_l);
  fEventHolderChain -> SetBranchAddress("tr_time", &tr_time, &b_tr_time);
  fEventHolderChain -> SetBranchAddress("rf_time", &rf_time, &b_rf_time);
  fEventHolderChain -> SetBranchAddress("l2bit", &l2bit, &b_l2bit);
  fEventHolderChain -> SetBranchAddress("l3bit", &l3bit, &b_l3bit);
  fEventHolderChain -> SetBranchAddress("helicity", &helicity, &b_helicity);
  fEventHolderChain -> SetBranchAddress("hlsc", &hlsc, &b_hlsc);
  fEventHolderChain -> SetBranchAddress("intt", &intt, &b_intt);
  fEventHolderChain -> SetBranchAddress("helicity_cor", &helicity_cor, &b_helicity_cor);
  fEventHolderChain -> SetBranchAddress("gpart", &gpart, &b_gpart);
  fEventHolderChain -> SetBranchAddress("id", id, &b_id);
  fEventHolderChain -> SetBranchAddress("stat", stat, &b_stat);
  fEventHolderChain -> SetBranchAddress("dc", dc, &b_dc);
  fEventHolderChain -> SetBranchAddress("cc", cc, &b_cc);
  fEventHolderChain -> SetBranchAddress("sc", sc, &b_sc);
  fEventHolderChain -> SetBranchAddress("ec", ec, &b_ec);
  fEventHolderChain -> SetBranchAddress("lec", lec, &b_lec);
  fEventHolderChain -> SetBranchAddress("st", st, &b_st);
  fEventHolderChain -> SetBranchAddress("p", p, &b_p);
  fEventHolderChain -> SetBranchAddress("m", m, &b_m);
  fEventHolderChain -> SetBranchAddress("q", q, &b_q);
  fEventHolderChain -> SetBranchAddress("b", b, &b_b);
  fEventHolderChain -> SetBranchAddress("cx", cx, &b_cx);
  fEventHolderChain -> SetBranchAddress("cy", cy, &b_cy);
  fEventHolderChain -> SetBranchAddress("cz", cz, &b_cz);
  fEventHolderChain -> SetBranchAddress("vx", vx, &b_vx);
  fEventHolderChain -> SetBranchAddress("vy", vy, &b_vy);
  fEventHolderChain -> SetBranchAddress("vz", vz, &b_vz);
  fEventHolderChain -> SetBranchAddress("dc_part", &dc_part, &b_dc_part);
  fEventHolderChain -> SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
  fEventHolderChain -> SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
  fEventHolderChain -> SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
  fEventHolderChain -> SetBranchAddress("tb_st", tb_st, &b_tb_st);
  fEventHolderChain -> SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
  fEventHolderChain -> SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
  fEventHolderChain -> SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
  fEventHolderChain -> SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
  fEventHolderChain -> SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
  fEventHolderChain -> SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
  fEventHolderChain -> SetBranchAddress("dc_vx", dc_vx, &b_dc_vx);
  fEventHolderChain -> SetBranchAddress("dc_vy", dc_vy, &b_dc_vy);
  fEventHolderChain -> SetBranchAddress("dc_vz", dc_vz, &b_dc_vz);
  fEventHolderChain -> SetBranchAddress("dc_vr", dc_vr, &b_dc_vr);
  fEventHolderChain -> SetBranchAddress("tl1_cx", tl1_cx, &b_tl1_cx);
  fEventHolderChain -> SetBranchAddress("tl1_cy", tl1_cy, &b_tl1_cy);
  fEventHolderChain -> SetBranchAddress("tl1_cz", tl1_cz, &b_tl1_cz);
  fEventHolderChain -> SetBranchAddress("tl1_x", tl1_x, &b_tl1_x);
  fEventHolderChain -> SetBranchAddress("tl1_y", tl1_y, &b_tl1_y);
  fEventHolderChain -> SetBranchAddress("tl1_z", tl1_z, &b_tl1_z);
  fEventHolderChain -> SetBranchAddress("tl1_r", tl1_r, &b_tl1_r);
  fEventHolderChain -> SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
  fEventHolderChain -> SetBranchAddress("ec_part", &ec_part, &b_ec_part);
  fEventHolderChain -> SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
  fEventHolderChain -> SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
  fEventHolderChain -> SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
  fEventHolderChain -> SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
  fEventHolderChain -> SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
  fEventHolderChain -> SetBranchAddress("etot", etot, &b_etot);
  fEventHolderChain -> SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
  fEventHolderChain -> SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
  fEventHolderChain -> SetBranchAddress("ec_t", ec_t, &b_ec_t);
  fEventHolderChain -> SetBranchAddress("ec_r", ec_r, &b_ec_r);
  fEventHolderChain -> SetBranchAddress("ech_x", ech_x, &b_ech_x);
  fEventHolderChain -> SetBranchAddress("ech_y", ech_y, &b_ech_y);
  fEventHolderChain -> SetBranchAddress("ech_z", ech_z, &b_ech_z);
  fEventHolderChain -> SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
  fEventHolderChain -> SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
  fEventHolderChain -> SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
  fEventHolderChain -> SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
  fEventHolderChain -> SetBranchAddress("sc_part", &sc_part, &b_sc_part);
  fEventHolderChain -> SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
  fEventHolderChain -> SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
  fEventHolderChain -> SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
  fEventHolderChain -> SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
  fEventHolderChain -> SetBranchAddress("edep", edep, &b_edep);
  fEventHolderChain -> SetBranchAddress("sc_t", sc_t, &b_sc_t);
  fEventHolderChain -> SetBranchAddress("sc_r", sc_r, &b_sc_r);
  fEventHolderChain -> SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
  fEventHolderChain -> SetBranchAddress("cc_part", &cc_part, &b_cc_part);
  fEventHolderChain -> SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
  fEventHolderChain -> SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
  fEventHolderChain -> SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
  fEventHolderChain -> SetBranchAddress("nphe", nphe, &b_nphe);
  fEventHolderChain -> SetBranchAddress("cc_t", cc_t, &b_cc_t);
  fEventHolderChain -> SetBranchAddress("cc_r", cc_r, &b_cc_r);
  fEventHolderChain -> SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
  fEventHolderChain -> SetBranchAddress("lac_part", &lac_part, &b_lac_part);
  fEventHolderChain -> SetBranchAddress("lec_sect", lec_sect, &b_lec_sect);
  fEventHolderChain -> SetBranchAddress("lec_hit", lec_hit, &b_lec_hit);
  fEventHolderChain -> SetBranchAddress("lec_stat", lec_stat, &b_lec_stat);
  fEventHolderChain -> SetBranchAddress("lec_etot", lec_etot, &b_lec_etot);
  fEventHolderChain -> SetBranchAddress("lec_ein", lec_ein, &b_lec_ein);
  fEventHolderChain -> SetBranchAddress("lec_t", lec_t, &b_lec_t);
  fEventHolderChain -> SetBranchAddress("lec_r", lec_r, &b_lec_r);
  fEventHolderChain -> SetBranchAddress("lec_x", lec_x, &b_lec_x);
  fEventHolderChain -> SetBranchAddress("lec_y", lec_y, &b_lec_y);
  fEventHolderChain -> SetBranchAddress("lec_z", lec_z, &b_lec_z);
  fEventHolderChain -> SetBranchAddress("lec_c2", lec_c2, &b_lec_c2);  

  return true ;
}

EventI * CLAS6EventHolder::GetEvent(const unsigned int event_id) {

  if ( event_id > (unsigned int) fMaxEvents ) return nullptr ; 

  //  static 
  CLAS6Event * event = new CLAS6Event() ; 
  fEventHolderChain -> GetEntry( event_id ) ; 
  /*
  event -> SetEventID( evntid ) ;
  event -> SetEventWeight( 1 ) ;

  //  event -> SetTargetPdg( tgt ) ; 
  event -> SetInLeptPdg( 11 ) ;
  event -> SetOutLeptPdg( 11 ) ; 
  
  event -> SetInLeptonKinematics( Ev, pxv, pyv, pzv ) ; 
  event -> SetOutLeptonKinematics( El, pxl, pyl, pzl ) ; 
  event -> SetInUnCorrLeptonKinematics( Ev, pxv, pyv, pzv ) ; 
  event -> SetOutUnCorrLeptonKinematics( El, pxl, pyl, pzl ) ; 
  
  event -> SetNProtons( nfp ) ; 
  event -> SetNNeutrons( nfn ) ; 
  event -> SetNPiP( nfpip ) ; 
  event -> SetNPiM( nfpim ) ; 
  event -> SetNPi0( nfpi0 ) ; 
  event -> SetNKP( nfkp ) ;
  event -> SetNKM( nfkm ) ; 
  event -> SetNK0( nfk0 ) ; 
  event -> SetNEM( nfem ) ; 
  event -> SetNOther( nfother ) ; 
  
  event -> SetVertex( vtxx, vtxy, vtxz, vtxt ) ; 


  // Set final state particle kinematics
  for ( unsigned int p = 0 ; p < (unsigned int) nf ; ++p ) {
    event -> SetFinalParticle( pdgf[p], Ef[p], pxf[p], pyf[p], pzf[p] ) ; 
    event -> SetFinalParticleUnCorr( pdgf[p], Ef[p], pxf[p], pyf[p], pzf[p] ) ; 
  }
*/
  return event ; 

}

void CLAS6EventHolder::Initialize() { 
  this -> LoadBranch();
}

void CLAS6EventHolder::Clear() { 

}
