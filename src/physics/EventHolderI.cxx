// _______________________________________________
/*
 * EventHolder Interface base class
 * 
 */
#include <iostream>
#include "physics/EventHolderI.h"

using namespace e4nu ; 

EventHolderI::EventHolderI() { 
  this->Initialize() ;
  std::cout << "Initialized member" << std::endl;
}

EventHolderI::~EventHolderI() {
  this->Clear();
}

EventHolderI::EventHolderI( const std::string file ) { 
  this->Initialize() ; 
  std::cout<< "Loading "<< file << " ... \n" ;
  if( this->LoadMembers( file ) ) fIsConfigured = true ; 
}

EventHolderI::EventHolderI( const std::vector<std::string> files ) { 
  this->Initialize() ; 
  /*
  for ( unsigned int file_id = 0 ; file_id < files.size() ; ++file_id ) {
    fIsConfigured *= this->LoadMembers( files[file_id] ) ; 
  }
  */
}

bool EventHolderI::LoadMembers( const std::string file ) {
  if ( ! fEventHolderChain -> Add( file.c_str() ) ) return false ; 
  return true ; 
} 

bool EventHolderI::InitChain(void) {
  if( !fEventHolderChain ) return false ; 
  fEventHolderChain->SetBranchAddress("iev", &iev, &b_iev);
  return true ; 
}

void EventHolderI::Initialize() { 
  fEventHolderChain = new TChain("gst","e4nu_analysis") ; 
  fIsConfigured = true ; 
}

void EventHolderI::Clear() { 
  fEventHolderChain = nullptr ;
  b_iev = nullptr ; 
}
