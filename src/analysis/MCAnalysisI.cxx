// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/MCAnalysisI.h"
#include "utils/ParticleUtils.h"
#include "conf/ParticleI.h"

using namespace e4nu ; 

MCAnalysisI::MCAnalysisI() {
  this->Initialize() ;
}

MCAnalysisI::~MCAnalysisI() {;}

bool MCAnalysisI::LoadData( const std::string file ) {
  if( ! fIsDataLoaded ) fData = new MCEventHolder(file);
  else return fIsDataLoaded ; 

  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

bool MCAnalysisI::LoadData( const std::string file, const unsigned int nmax ) {

  if( ! fIsDataLoaded ) fData = new MCEventHolder( file, nmax ) ;
  else return fIsDataLoaded ; 

  fIsDataLoaded = true ;

  return fIsDataLoaded ; 
}

EventI * MCAnalysisI::GetEvent( const unsigned int event_id ) {
  return fData -> GetEvent(event_id) ; 
}

EventI * MCAnalysisI::GetValidEvent( const unsigned int event_id ) {

  MCEvent * event = (MCEvent*) fData -> GetEvent(event_id) ; 

  if ( event -> GetTargetPdg() != GetConfiguredTarget() ) return nullptr ; 

  TLorentzVector in_mom = event -> GetInLepton4Mom() ; 
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 

  double EBeam = GetConfiguredEBeam() ; 
  if ( in_mom.E() != EBeam ) return nullptr ; 

  if( ApplyReso() ) {
    this -> SmearParticles( event ) ; 
  }

  return event ; 
}

void MCAnalysisI::SmearParticles( MCEvent * event ) {
  
  TLorentzVector in_mom = event -> GetInLepton4Mom() ; 
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 
  std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticles4Mom() ;

  double EBeam = GetConfiguredEBeam() ; 

  utils::ApplyResolution( conf::kPdgElectron , in_mom, EBeam ) ; 
  event -> EventI::SetOutLeptonKinematics( in_mom ) ; 
  
  utils::ApplyResolution( conf::kPdgElectron , out_mom, EBeam ) ; 
  event -> EventI::SetOutLeptonKinematics( out_mom ) ; 

  for( std::map<int,std::vector<TLorentzVector>>::iterator it = part_map.begin() ; it != part_map.end() ; ++it ) {
    std::cout << it->first << std::endl;

    std::vector<TLorentzVector> vtemp ; 
    for( unsigned int i = 0 ; i < (it->second).size() ; ++i ) { 
      TLorentzVector temp = (it->second)[i] ; 
      utils::ApplyResolution( it->first, temp, EBeam ) ;
      vtemp.push_back(temp) ; 
    }
    part_map[it->first] = vtemp ; 
  }
  event -> EventI::SetFinalParticlesKinematics( part_map ) ; 
} 

unsigned int MCAnalysisI::GetNEvents( void ) const {
  return (unsigned int) fData ->GetNEvents() ; 
}
void MCAnalysisI::Initialize() { 

}
