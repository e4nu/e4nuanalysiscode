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
}

EventHolderI::~EventHolderI() {
  this->Clear();
}

EventHolderI::EventHolderI( const std::string file ) { 
  this->Initialize() ; 
  std::cout<< "Loading "<< file << " ... \n" ;
  if( this->LoadMembers( file ) ) fIsConfigured = true ; 
}

EventHolderI::EventHolderI( const std::string file, const int nmaxevents ) { 
  this->Initialize() ; 
  fMaxEvents = nmaxevents ; 
  std::cout<< "Loading "<< fMaxEvents << " from " << file << " ... \n" ;
  if( this->LoadMembers( file ) ) fIsConfigured = true ; 
}

EventHolderI::EventHolderI( const std::vector<std::string> files ) { 
  this->Initialize() ; 
  for ( unsigned int file_id = 0 ; file_id < files.size() ; ++file_id ) {
    fIsConfigured *= this->LoadMembers( files[file_id] ) ; 
  }
}

bool EventHolderI::LoadMembers( const std::string file ) {
  if ( ! fEventHolderChain -> Add( file.c_str() ) ) return false ; 
  return true ; 
} 

void EventHolderI::Initialize() { 
  fEventHolderChain = new TChain("gst","e4nu_analysis") ; 
  fIsConfigured = true ; 
  fMaxEvents = -1 ; 
}

void EventHolderI::Clear() { 
  fEventHolderChain = nullptr ;
  for ( unsigned int i = 0 ; i < fEvents.size() ; ++i ) {
    delete fEvents[i] ;
  }
  fEvents.clear() ; 
  fMaxEvents = -1 ;

}
