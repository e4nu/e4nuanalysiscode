// _______________________________________________
/*
 * Analysis Interface base class
 * 
 */
#include <iostream>
#include "analysis/MCAnalysisI.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "utils/DetectorUtils.h"
#include "conf/FiducialCutI.h"
#include "conf/AccpetanceMapsI.h"

using namespace e4nu ; 

MCAnalysisI::MCAnalysisI() {
  this->Initialize() ;
}

MCAnalysisI::~MCAnalysisI() {
  delete fData;
  for( auto it = kAccMap.begin() ; it != kAccMap.end() ; ++it ) {
    delete kAccMap[it->first] ; 
    delete kGenMap[it->first] ; 
  }

  for( auto it=kAcceptanceMap.begin(); it != kAcceptanceMap.end() ; ++it ) {
    kAcceptanceMap[it->first]->Close();
    delete kAcceptanceMap[it->first];
  }
  kAcceptanceMap.clear();
  kAccMap.clear();
  kGenMap.clear();
}

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
  if( !event ) return nullptr ; 

  ++fEventsBeforeCuts ;

  TLorentzVector in_mom = event -> GetInLepton4Mom() ; 
  TLorentzVector out_mom = event -> GetOutLepton4Mom() ; 

  //  std::cout << "out_mom.P() = " << out_mom.P() << std::endl;

  // Check run is correct
  double EBeam = GetConfiguredEBeam() ; 
  if ( in_mom.E() != EBeam ) return nullptr ;  
  if ( (unsigned int) event -> GetTargetPdg() != GetConfiguredTarget() ) return nullptr ; 

  // Check weight is physical
  double wght = event->GetWeight() ; 
  if ( wght < 0 || wght > 10 ) return nullptr ; 

  // Apply Fiducial volume cuts
  if( ApplyFiducial() ) {
    if (! kFiducialCut -> EFiducialCut(EBeam, out_mom.Vect() ) ) return nullptr ; 
    ++fNEventsAfterFiducial ;
  }
    
  // Apply smaring to particles
  if( ApplyReso() ) {
    this -> SmearParticles( event ) ; 
  }

  // Get Topology Definition
  std::map<int,unsigned int> Topology = GetTopology(); 
  std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticles4Mom() ;

  // Signal event
  bool is_signal = false ; 
  double acc_wght = 1 ; 

  for( auto it = Topology.begin() ; it != Topology.end() ; ++it ) {
    if( it->first == conf::kPdgElectron ) {
      // Apply acceptance for electrons
      acc_wght *= utils::GetAcceptanceMapWeight( kAccMap[it->first], kGenMap[it->first], out_mom ) ; 
      if( fabs( acc_wght ) != acc_wght ) return nullptr ; 
      is_signal = true ; 
    } else if( part_map.count( it->first ) && part_map[it->first].size() == it->second ) {
      // If we have an exclusive event, apply acceptance weight as well for signal event
      is_signal = true ; 
      for( unsigned int i = 0 ; i < it->second ; ++i ) {
	// Do I need particle id ? 
	acc_wght *= utils::GetAcceptanceMapWeight( kAccMap[it->first], kGenMap[it->first], part_map[it->first][i] );
	if( fabs( acc_wght) != acc_wght ) return nullptr ; 
      } 
    } else {
      is_signal = false ; 
    }

    if( is_signal && ApplyAccWeights() ) wght *= acc_wght ; 
      
    // Background event
    // BACKGROUND SUBSTRACTION METHOD HERE
  }

  // Apply Mottxsec weight
  //
  if( IsElectronData() ) {
    wght *= utils::GetXSecScale(out_mom, EBeam, true ) ; 
  }

  // Set Final weight
  event->SetWeight( wght ) ; 

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

  // Initialize fiducial for this run
  kFiducialCut = std::unique_ptr<Fiducial>( new Fiducial()) ;
  double EBeam = GetConfiguredEBeam() ; 

  if( ApplyFiducial() ) {
    kFiducialCut -> InitPiMinusFit( EBeam ) ; 
    kFiducialCut -> InitEClimits(); 
    kFiducialCut -> up_lim1_ec -> Eval(60) ;
    kFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), GetConfiguredTarget() , EBeam ) ;
    kFiducialCut -> SetFiducialCutParameters( EBeam ) ;
 }

  // Initialize acceptance map histograms from file
  if( ApplyAccWeights() ) { 
    kAcceptanceMap = conf::GetAcceptanceFileMap( GetConfiguredTarget(), EBeam ) ; 

    kAccMap[conf::kPdgElectron] = (TH3D*) kAcceptanceMap[conf::kPdgElectron] -> Get("Accepted Particles") ;
    kAccMap[conf::kPdgProton] = (TH3D*) kAcceptanceMap[conf::kPdgProton] -> Get("Accepted Particles") ;
    kAccMap[conf::kPdgPiP] = (TH3D*) kAcceptanceMap[conf::kPdgPiP] -> Get("Accepted Particles") ;
    kAccMap[conf::kPdgPiM] = (TH3D*) kAcceptanceMap[conf::kPdgPiM] -> Get("Accepted Particles") ;
    
    kGenMap[conf::kPdgElectron] = (TH3D*) kAcceptanceMap[conf::kPdgElectron] -> Get("Generated Particles") ;
    kGenMap[conf::kPdgProton] = (TH3D*) kAcceptanceMap[conf::kPdgProton] -> Get("Generated Particles") ;
    kGenMap[conf::kPdgPiP] = (TH3D*) kAcceptanceMap[conf::kPdgPiP] -> Get("Generated Particles") ;
    kGenMap[conf::kPdgPiM] = (TH3D*) kAcceptanceMap[conf::kPdgPiM] -> Get("Generated Particles") ;
  }  

}

bool MCAnalysisI::Finalise( const std::string out_file ) {

  std::cout << " Total Number of Events Processed = " << fEventsBeforeCuts << std::endl;
  std::cout << " Events after fiducial cuts = " << fNEventsAfterFiducial << std::endl;

  return true ; 
}
