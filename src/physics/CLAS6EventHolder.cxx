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

  /* NOTICE:
     The CLAS event format is different from the GENIE MC format
     The code GetCharge_FilterData.C applies some generic CLAS-specific cuts
     and stores the information in a genie-like format. 
     This code will print a warning if the wrong format is used. 
     The available GENIE-like CLAS6 data can be found here: 
     /w/hallb-scifs17exp/clas/claseg2/apapadop/GetCharge_genie_filtered_data_e2a_ep_%s_%s_neutrino6_united4_radphot_test_100M.root
  */

  fEventHolderChain->SetBranchAddress("iev", &iev, &b_iev);
  fEventHolderChain->SetBranchAddress("tgt", &tgt, &b_tgt);
  fEventHolderChain->SetBranchAddress("Ev", &Ev, &b_Ev);
  fEventHolderChain->SetBranchAddress("pxv", &pxv, &b_pxv);
  fEventHolderChain->SetBranchAddress("pyv", &pyv, &b_pyv);
  fEventHolderChain->SetBranchAddress("pzv", &pzv, &b_pzv);
  fEventHolderChain->SetBranchAddress("El", &El, &b_El);
  fEventHolderChain->SetBranchAddress("pxl", &pxl, &b_pxl);
  fEventHolderChain->SetBranchAddress("pyl", &pyl, &b_pyl);
  fEventHolderChain->SetBranchAddress("pzl", &pzl, &b_pzl);
  fEventHolderChain->SetBranchAddress("nfp", &nfp, &b_nfp);
  fEventHolderChain->SetBranchAddress("nfn", &nfn, &b_nfn);
  fEventHolderChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
  fEventHolderChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
  fEventHolderChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
  fEventHolderChain->SetBranchAddress("nf", &nf, &b_nf);
  fEventHolderChain->SetBranchAddress("Ef", Ef, &b_Ef);
  fEventHolderChain->SetBranchAddress("pxf", pxf, &b_pxf);
  fEventHolderChain->SetBranchAddress("pyf", pyf, &b_pyf);
  fEventHolderChain->SetBranchAddress("pzf", pzf, &b_pzf);
  fEventHolderChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
  fEventHolderChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
  fEventHolderChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
  fEventHolderChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
  fEventHolderChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);

  return true ;
}

EventI * CLAS6EventHolder::GetEvent(const unsigned int event_id) {

  if ( event_id > (unsigned int) fMaxEvents ) return nullptr ; 

  //  static 
  CLAS6Event * event = new CLAS6Event() ; 
  if ( event_id > (unsigned int) fMaxEvents ) return nullptr ; 
  fEventHolderChain -> GetEntry( event_id ) ; 

  event -> SetEventID( iev ) ;
  event -> SetEventWeight( 1. ) ;
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
  event -> SetVertex( vtxx, vtxy, vtxz, vtxt ) ; 

  // Set final state particle kinematics
  for ( unsigned int p = 0 ; p < (unsigned int) nf ; ++p ) {
    event -> SetFinalParticle( pdgf[p], Ef[p], pxf[p], pyf[p], pzf[p] ) ; 
  }

  return event ; 
}

void CLAS6EventHolder::Initialize() { 
  this -> LoadBranch();
}

void CLAS6EventHolder::Clear() { 

}
