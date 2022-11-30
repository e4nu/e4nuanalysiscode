// _______________________________________________
/*
 * MCEventHolder implementation 
 */

#include "physics/MCEventHolder.h"

using namespace e4nu ; 

MCEventHolder::MCEventHolder(): EventHolderI() { ; } 

MCEventHolder::~MCEventHolder() {
  this->Clear();
}

MCEventHolder::MCEventHolder( const std::string file ): EventHolderI( file ) { 
  this->Initialize() ; 
  std::cout << "Initialized member MC " << std::endl;
  LoadEvents();
}

MCEventHolder::MCEventHolder( const std::vector<std::string> files ): EventHolderI( files ) { 
  this->Initialize() ; 
}


bool MCEventHolder::LoadEvents(void) {
  if( !fEventHolderChain ) return false ; 
   
  fEventHolderChain->SetBranchAddress("iev", &fiev, &fb_iev);

  unsigned int nevents = fEventHolderChain->GetEntries();
  std::cout << "Loading "<< nevents <<" events ..." <<std::endl;

  for( unsigned int i = 0 ; i < nevents ; ++i ) { 
    fEventHolderChain->GetEntry( i ) ; 
    std::cout << " iev = " << fiev << std::endl;
  }  
  return true ; 
}

void MCEventHolder::Initialize() { 
  
}

void MCEventHolder::Clear() { 
  //  b_iev = nullptr ;   
}
