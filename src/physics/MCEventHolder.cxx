// _______________________________________________
/*
 * MCEventHolder implementation 
 */

#include "physics/MCEventHolder.h"
#include "physics/MCEvent.h"

using namespace e4nu ; 

MCEventHolder::MCEventHolder(): EventHolderI() { this->Initialize(); } 

MCEventHolder::~MCEventHolder() {
  this->Clear();
}

MCEventHolder::MCEventHolder( const std::string file ): EventHolderI( file ) { 
  this->Initialize() ; 
}

MCEventHolder::MCEventHolder( const std::string file, const int nmaxevents ): EventHolderI( file, nmaxevents ) { 
  this->Initialize() ; 
}

MCEventHolder::MCEventHolder( const std::vector<std::string> files ): EventHolderI( files ) { 
  this->Initialize() ; 
}


bool MCEventHolder::LoadBranch(void) {
  if( !fEventHolderChain ) return false ; 

  fEventHolderChain->SetBranchAddress("iev", &iev, &b_iev);
  fEventHolderChain->SetBranchAddress("neu", &neu, &b_neu);
  fEventHolderChain->SetBranchAddress("fspl", &fspl, &b_fspl);
  fEventHolderChain->SetBranchAddress("tgt", &tgt, &b_tgt);
  fEventHolderChain->SetBranchAddress("Z", &Z, &b_Z);
  fEventHolderChain->SetBranchAddress("A", &A, &b_A);
  fEventHolderChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
  fEventHolderChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
  fEventHolderChain->SetBranchAddress("resid", &resid, &b_resid);
  fEventHolderChain->SetBranchAddress("sea", &sea, &b_sea);
  fEventHolderChain->SetBranchAddress("qel", &qel, &b_qel);
  fEventHolderChain->SetBranchAddress("mec", &mec, &b_mec);//go down from here
  fEventHolderChain->SetBranchAddress("res", &res, &b_res);
  fEventHolderChain->SetBranchAddress("dis", &dis, &b_dis);
  fEventHolderChain->SetBranchAddress("coh", &coh, &b_coh);
  fEventHolderChain->SetBranchAddress("dfr", &dfr, &b_dfr);
  fEventHolderChain->SetBranchAddress("imd", &imd, &b_imd);
  fEventHolderChain->SetBranchAddress("imdanh", &imdanh, &b_imdanh);
  fEventHolderChain->SetBranchAddress("singlek", &singlek, &b_singlek);
  fEventHolderChain->SetBranchAddress("nuel", &nuel, &b_nuel);
  fEventHolderChain->SetBranchAddress("em", &em, &b_em);
  fEventHolderChain->SetBranchAddress("cc", &cc, &b_cc);
  fEventHolderChain->SetBranchAddress("nc", &nc, &b_nc);
  fEventHolderChain->SetBranchAddress("charm", &charm, &b_charm);
  fEventHolderChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
  fEventHolderChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
  fEventHolderChain->SetBranchAddress("wght", &wght, &b_wght);
  fEventHolderChain->SetBranchAddress("xs", &xs, &b_xs);
  fEventHolderChain->SetBranchAddress("ys", &ys, &b_ys);
  fEventHolderChain->SetBranchAddress("ts", &ts, &b_ts);
  fEventHolderChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
  fEventHolderChain->SetBranchAddress("Ws", &Ws, &b_Ws);
  fEventHolderChain->SetBranchAddress("x", &x, &b_x);
  fEventHolderChain->SetBranchAddress("y", &y, &b_y);
  fEventHolderChain->SetBranchAddress("t", &t, &b_t);
  fEventHolderChain->SetBranchAddress("Q2", &Q2, &b_Q2);
  fEventHolderChain->SetBranchAddress("W", &W, &b_W);
  fEventHolderChain->SetBranchAddress("EvRF", &EvRF, &b_EvRF);
  fEventHolderChain->SetBranchAddress("Ev", &Ev, &b_Ev);
  fEventHolderChain->SetBranchAddress("pxv", &pxv, &b_pxv);
  fEventHolderChain->SetBranchAddress("pyv", &pyv, &b_pyv);
  fEventHolderChain->SetBranchAddress("pzv", &pzv, &b_pzv);
  fEventHolderChain->SetBranchAddress("En", &En, &b_En);
  fEventHolderChain->SetBranchAddress("pxn", &pxn, &b_pxn);
  fEventHolderChain->SetBranchAddress("pyn", &pyn, &b_pyn);
  fEventHolderChain->SetBranchAddress("pzn", &pzn, &b_pzn);
  fEventHolderChain->SetBranchAddress("El", &El, &b_El);
  fEventHolderChain->SetBranchAddress("pxl", &pxl, &b_pxl);
  fEventHolderChain->SetBranchAddress("pyl", &pyl, &b_pyl);
  fEventHolderChain->SetBranchAddress("pzl", &pzl, &b_pzl);
  fEventHolderChain->SetBranchAddress("pl", &pl, &b_pl);
  fEventHolderChain->SetBranchAddress("cthl", &cthl, &b_cthl);
  fEventHolderChain->SetBranchAddress("nfp", &nfp, &b_nfp);
  fEventHolderChain->SetBranchAddress("nfn", &nfn, &b_nfn);
  fEventHolderChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
  fEventHolderChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
  fEventHolderChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
  fEventHolderChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
  fEventHolderChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
  fEventHolderChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
  fEventHolderChain->SetBranchAddress("nfem", &nfem, &b_nfem);
  fEventHolderChain->SetBranchAddress("nfother", &nfother, &b_nfother);
  fEventHolderChain->SetBranchAddress("nip", &nip, &b_nip);
  fEventHolderChain->SetBranchAddress("nin", &nin, &b_nin);
  fEventHolderChain->SetBranchAddress("nipip", &nipip, &b_nipip);
  fEventHolderChain->SetBranchAddress("nipim", &nipim, &b_nipim);
  fEventHolderChain->SetBranchAddress("nipi0", &nipi0, &b_nipi0);
  fEventHolderChain->SetBranchAddress("nikp", &nikp, &b_nikp);
  fEventHolderChain->SetBranchAddress("nikm", &nikm, &b_nikm);
  fEventHolderChain->SetBranchAddress("nik0", &nik0, &b_nik0);
  fEventHolderChain->SetBranchAddress("niem", &niem, &b_niem);
  fEventHolderChain->SetBranchAddress("niother", &niother, &b_niother);
  fEventHolderChain->SetBranchAddress("ni", &ni, &b_ni);
  fEventHolderChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
  fEventHolderChain->SetBranchAddress("resc", resc, &b_resc);
  fEventHolderChain->SetBranchAddress("Ei", Ei, &b_Ei);
  fEventHolderChain->SetBranchAddress("pxi", pxi, &b_pxi);
  fEventHolderChain->SetBranchAddress("pyi", pyi, &b_pyi);
  fEventHolderChain->SetBranchAddress("pzi", pzi, &b_pzi);
  fEventHolderChain->SetBranchAddress("nf", &nf, &b_nf);
  fEventHolderChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
  fEventHolderChain->SetBranchAddress("Ef", Ef, &b_Ef);
  fEventHolderChain->SetBranchAddress("pxf", pxf, &b_pxf);
  fEventHolderChain->SetBranchAddress("pyf", pyf, &b_pyf);
  fEventHolderChain->SetBranchAddress("pzf", pzf, &b_pzf);
  fEventHolderChain->SetBranchAddress("pf", pf, &b_pf);
  fEventHolderChain->SetBranchAddress("cthf", cthf, &b_cthf);
  fEventHolderChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
  fEventHolderChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
  fEventHolderChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
  fEventHolderChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
  fEventHolderChain->SetBranchAddress("sumKEf", &sumKEf, &b_sumKEf);
  fEventHolderChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);
  return true ;
}

EventI * MCEventHolder::GetEvent(const unsigned int event_id) {

  if ( event_id > (unsigned int) fMaxEvents ) return nullptr ; 

  MCEvent * event = new MCEvent() ; 
  fEventHolderChain->GetEntry( event_id ) ; 
  
  event -> SetEventID( iev ) ;
  event -> SetWeight( wght ) ;
  event -> SetIsEM( em ) ;   
  event -> SetIsCC( cc ) ; 
  event -> SetIsNC( nc ) ; 
  event -> SetIsQEL( qel ) ; 
  event -> SetIsRES( res ) ; 
  event -> SetIsDIS( dis ) ; 
  event -> SetIsMEC( mec ) ; 
  event -> SetTargetPdg( tgt ) ; 
  event -> SetInLeptPdg( 11 ) ;
  event -> SetOutLeptPdg( 11 ) ; 
  
  event -> SetInLeptonKinematics( Ev, pxv, pyv, pzv ) ; 
  event -> SetOutLeptonKinematics( El, pxl, pyl, pzl ) ; 
  
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
  event -> SetHitNuclPdg( hitnuc ) ; 
  event -> SetTrueQ2s( Q2s ) ; 
  event -> SetTrueWs( Ws ) ;
  event -> SetTruexs( xs ) ; 
  event -> SetTrueys( ys ) ; 
  event -> SetTrueQ2( Q2 ) ; 
  event -> SetTrueW( W ) ;
  event -> SetTruex( x ) ; 
  event -> SetTruey( y ) ; 

  // Set final state particle kinematics
  for ( unsigned int p = 0 ; p < (unsigned int) nf ; ++p ) {
    event -> SetFinalParticle( pdgf[p], Ef[p], pxf[p], pyf[p], pzf[p] ) ; 
  }
  return event ; 
}

void MCEventHolder::Initialize() { 
  this->LoadBranch();
}

void MCEventHolder::Clear() { 

    delete b_iev ;   
    delete b_neu ;   
    delete b_fspl ;   
    delete b_tgt ;   
    delete b_Z ;   
    delete b_A ;   
    delete b_hitnuc ;   
    delete b_hitqrk ;   
    delete b_resid ;   
    delete b_sea ;   
    delete b_qel ;   
    delete b_mec ;   
    delete b_res ;   
    delete b_dis ;   
    delete b_coh ;   
    delete b_dfr ;   
    delete b_imd ;   
    delete b_imdanh ;   
    delete b_singlek ;   
    delete b_nuel ;   
    delete b_em ;   
    delete b_cc ;   
    delete b_nc ;   
    delete b_charm ;   
    delete b_neut_code ;   
    delete b_nuance_code ;   
    delete b_wght ;   
    delete b_xs ;   
    delete b_ys ;   
    delete b_ts ;   
    delete b_Q2s ;   
    delete b_Ws ;   
    delete b_x ;   
    delete b_y ;   
    delete b_t ;   
    delete b_Q2 ;   
    delete b_W ;   
    delete b_EvRF ;   
    delete b_Ev ;   
    delete b_pxv ;   
    delete b_pyv ;   
    delete b_pzv ;   
    delete b_En ;   
    delete b_pxn ;   
    delete b_pyn ;   
    delete b_pzn ;   
    delete b_El ;   
    delete b_pxl ;   
    delete b_pyl ;   
    delete b_pzl ;   
    delete b_pl ;   
    delete b_cthl ;   
    delete b_nfp ;   
    delete b_nfn ;   
    delete b_nfpip ;   
    delete b_nfpim ;   
    delete b_nfpi0 ;   
    delete b_nfkp ;   
    delete b_nfkm ;   
    delete b_nfk0 ;   
    delete b_nfem ;   
    delete b_nfother ;   
    delete b_nip ;   
    delete b_nin ;   
    delete b_nipip ;   
    delete b_nipim ;   
    delete b_nipi0 ;   
    delete b_nikp ;   
    delete b_nikm ;   
    delete b_nik0 ;   
    delete b_niem ;   
    delete b_niother ;   
    delete b_ni ;   
    delete b_pdgi ;   
    delete b_resc ;   
    delete b_Ei ;   
    delete b_pxi ;   
    delete b_pyi ;   
    delete b_pzi ;   
    delete b_nf ;   
    delete b_pdgf ;   
    delete b_Ef ;   
    delete b_pxf ;   
    delete b_pyf ;   
    delete b_pzf ;   
    delete b_pf ;   
    delete b_cthf ;   
    delete b_vtxx ;   
    delete b_vtxy ;   
    delete b_vtxz ;   
    delete b_vtxt ;   
    delete b_sumKEf ;   
    delete b_calresp0 ;   

}
