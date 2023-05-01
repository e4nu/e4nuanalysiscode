// _______________________________________________
/*
 * Event Interface base class
 * 
 */
#include <iostream>
#include "physics/EventI.h"
#include "conf/ParticleI.h"
#include "utils/DetectorUtils.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"

using namespace e4nu ; 

EventI::EventI() { 
  this->Initialize() ;
}

EventI::~EventI() {
  this->Clear();
}

void EventI::SetOutLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fOutLepton.SetPxPyPzE( px, py, pz, E ) ; 
  return ; 
}

void EventI::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fInLepton.SetPxPyPzE( px, py, pz, E ) ; 
  return ; 
} 

void EventI::SetFinalParticle( const int pdg, const double E, const double px, const double py, const double pz ) {
  TLorentzVector mom;
  mom.SetPxPyPzE( px, py, pz, E ) ; 
  if( fFinalParticles.find(pdg) == fFinalParticles.end() ) {
    std::vector<TLorentzVector> vct ;
    vct.push_back(mom); 
    fFinalParticles.insert( std::pair<int,std::vector<TLorentzVector>>(pdg, vct) ) ; 
  } else {
    fFinalParticles[pdg].push_back( mom ) ; 
  }
}

void EventI::SetOutUnCorrLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fOutLeptonUnCorr.SetPxPyPzE( px, py, pz, E ) ; 
  return ; 
}

void EventI::SetInUnCorrLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fInLeptonUnCorr.SetPxPyPzE( px, py, pz, E ) ; 
  return ; 
} 

void EventI::SetFinalParticleUnCorr( const int pdg, const double E, const double px, const double py, const double pz ) {
  TLorentzVector mom;
  mom.SetPxPyPzE( px, py, pz, E ) ; 
  if( fFinalParticlesUnCorr.find(pdg) == fFinalParticlesUnCorr.end() ) {
    std::vector<TLorentzVector> vct ;
    vct.push_back(mom); 
    fFinalParticlesUnCorr.insert( std::pair<int,std::vector<TLorentzVector>>(pdg, vct) ) ; 
  } else {
    fFinalParticlesUnCorr[pdg].push_back( mom ) ; 
  }
}

void EventI::StoreAnalysisRecord( unsigned int analysis_step ) {
  double weight = this->GetTotalWeight() ; 
  std::map<int,std::vector<TLorentzVector>> part_map = this->GetFinalParticles4Mom() ; 
  std::vector<int> pdg_list ; 
  for( auto it = part_map.begin() ; it != part_map.end() ; ++it ) { 
    for( unsigned int i = 0 ; i < part_map[it->first].size() ; ++i ) pdg_list.push_back( it->first ) ; 
  }
  std::pair<std::vector<int>,double> pair ( pdg_list, weight ) ; 
  fAnalysisRecord[analysis_step] = pair ; 
  return ; 
}

unsigned int EventI::GetEventMultiplicity( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) {
  // return number of charged particles in event
  unsigned int multiplicity = 0 ; 
  for( auto it = hadronic_system.begin() ; it != hadronic_system.end() ; ++it ) {
    multiplicity += std::abs( utils::GetParticleCharge( it->first ) ); 
  }
  return multiplicity ; 
}

unsigned int EventI::GetNSignalParticles( std::map<int,std::vector<TLorentzVector>> hadronic_system, const std::map<int,unsigned int> topology ) {
  unsigned int N_signal = 0 ; 
  for( auto it = hadronic_system.begin() ; it != hadronic_system.end() ; ++it ) {
    if( it->first == conf::kPdgElectron ) continue ; 
    if( topology.find(it->first) != topology.end() ) {
      N_signal += hadronic_system[it->first].size() ;
    }
  }

  return N_signal ;
}

int EventI::GetEventTotalVisibleCharge( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) {
  // return number of charged particles in event
  unsigned int charge = 0 ; 
  for( auto it = hadronic_system.begin() ; it != hadronic_system.end() ; ++it ) {
    charge += utils::GetParticleCharge( it->first ) ; 
  }
  return charge ; 
}

double EventI::GetObservable( const std::string observable ) {
  unsigned int target = fTargetPdg ; 
  double EBeam = GetInLepton4Mom().E();
  TLorentzVector ef4mom = GetOutLepton4Mom() ;
  TLorentzVector p4mom, pip4mom, pim4mom ; 
  bool event_wproton = false, event_wpip = false, event_wpim = false ; 
  
  if( fFinalParticles.find( conf::kPdgProton ) != fFinalParticles.end() ) { 
    event_wproton = true ;
    double max_mom = 0 ; 
    for ( unsigned int i = 0 ; i < fFinalParticles[conf::kPdgProton].size() ; ++i ) {
      if( fFinalParticles[conf::kPdgProton][i].P() > max_mom ) {
	max_mom = fFinalParticles[conf::kPdgProton][i].P() ; 
	p4mom = fFinalParticles[conf::kPdgProton][i] ;
      }
    }
  }

  if( fFinalParticles.find( conf::kPdgPiP ) != fFinalParticles.end() ) { 
    event_wpip = true ;
    double max_mom = 0 ; 
    for ( unsigned int i = 0 ; i < fFinalParticles[conf::kPdgPiP].size() ; ++i ) {
      if( fFinalParticles[conf::kPdgPiP][i].P() > max_mom ) {
	max_mom = fFinalParticles[conf::kPdgPiP][i].P() ; 
	pip4mom = fFinalParticles[conf::kPdgPiP][i] ;
      }
    }
  }

  if( fFinalParticles.find( conf::kPdgPiM ) != fFinalParticles.end() ) { 
    event_wpim = true ;
    double max_mom = 0 ; 
    for ( unsigned int i = 0 ; i < fFinalParticles[conf::kPdgPiM].size() ; ++i ) {
      if( fFinalParticles[conf::kPdgPiM][i].P() > max_mom ) {
	max_mom = fFinalParticles[conf::kPdgPiM][i].P() ; 
	pim4mom = fFinalParticles[conf::kPdgPiM][i] ;
      }
    }
  }

  if( observable == "ECal" ) {
    if ( event_wproton == false ) return 0 ; 
    return utils::GetECal( ef4mom.E(), GetFinalParticles4Mom(), target ) ; 
  } else if ( observable == "RecoEnu" ) {
    return utils::GetRecoEnu( ef4mom, target ) ;
  } else if ( observable == "QELRecoEnu" ) {
    return utils::GetQELRecoEnu( ef4mom, target ) ;
  } else if ( observable == "EnergyTransfer" ) {
    return utils::GetEnergyTransfer( ef4mom, EBeam ) ; 
  } else if ( observable == "RecoQ2" ) {
    return utils::GetRecoQ2( ef4mom, EBeam ) ; 
  } else if ( observable == "RecoXBJK" ) {
    return utils::GetRecoXBJK( ef4mom, EBeam ) ;
  } else if ( observable == "RecoW" ) {
    return utils::GetRecoW(ef4mom,EBeam ); 
  } else if ( observable == "DeltaPT" ) {
    if ( event_wproton == false ) return 0; 
    return utils::DeltaPT( ef4mom.Vect(), p4mom.Vect() ).Mag() ;
  } else if ( observable == "DeltaAlphaT" ) {
    if ( event_wproton == false ) return 0; 
    return utils::DeltaAlphaT( ef4mom.Vect(), p4mom.Vect() ) ; 
  } else if ( observable == "DeltaPhiT" ) {
    if ( event_wproton == false ) return 0; 
    return utils::DeltaPhiT( ef4mom.Vect(), p4mom.Vect() ) ;
  } else if ( observable == "LeadingPMom" ) {
    return p4mom.P() ; 
  } else if ( observable == "OutEMom" ) {
    return ef4mom.P() ; 
  } else if ( observable == "OutEPhi" ) {
    ef4mom = GetOutLepton4Mom() ;
    ef4mom.SetPhi( ef4mom.Phi()+TMath::Pi() ) ;
    return ef4mom.Phi()* 180 / TMath::Pi();
  } else if ( observable == "Sector" ) {
    ef4mom = GetOutLepton4Mom() ;
    ef4mom.SetPhi( ef4mom.Phi()+TMath::Pi() ) ;
    return utils::GetSector( ef4mom.Phi() ) ; 
  } else if ( observable == "Weight" ) { 
    return this->GetEventWeight() ;
  } else if ( observable == "LeadingPiPMom" ) {
    if( !event_wpip ) return 0 ; 
    return pip4mom.P() ; 
  } else if ( observable == "LeadingPiMMom" ) {
    if( !event_wpim ) return 0 ; 
    return pim4mom.P() ; 
  } else if ( observable == "LeadingPiPTheta" ) {
    if( !event_wpip ) return 0 ; 
    return pip4mom.Theta() * 180 / TMath::Pi() ; 
  } else if ( observable == "LeadingPiMTheta" ) {
    if( !event_wpim ) return 0 ; 
    return pim4mom.Theta() * 180 / TMath::Pi() ; 
  } else if ( observable == "HadSystemDeltaAlphaT" ) {
    return utils::DeltaAlphaT( ef4mom, GetFinalParticles4Mom() ) ;
  } else if ( observable == "HadSystemDeltaPhiT" ) {
    return utils::DeltaPhiT( ef4mom, GetFinalParticles4Mom() ) ;
  } else if ( observable == "HadSystemDeltaPT" ) {
    return utils::DeltaPT( ef4mom, GetFinalParticles4Mom() ).Mag() ;
  }

  std::cout << observable << " is NOT defined " << std::endl;
  return 0 ; 
}


void EventI::SetMottXSecWeight(void) { 
  // Set Mott XSec
  double reco_Q2 = utils::GetRecoQ2( this->GetOutLepton4Mom(), this->GetInLepton4Mom().E() ) ;
  fMottXSecWght = std::pow( reco_Q2, 2 ) ; 
}

TVector3 EventI::GetRecoq3() const { 
  return utils::GetRecoq3( fOutLepton, fInLepton.E() ) ; 
}

void EventI::Initialize() { 
  fFinalParticles.clear() ; 
  fFinalParticlesUnCorr.clear() ; 
  fIsMC = false ; 
  fEventID = 0 ; 
  fWeight = 0 ; 
  fEventID = 0 ; 
  fTargetPdg = 0 ; 
  fInLeptPdg = 11 ; 
  fOutLeptPdg = 11 ; 
  fNP = 0 ; 
  fNN = 0 ; 
  fNPiP = 0 ; 
  fNPiM = 0 ; 
  fNPi0 = 0 ; 
  fNKM = 0 ; 
  fNKP = 0 ; 
  fNK0 = 0 ; 
  fNEM = 0 ; 
  fNOther = 0 ;

  fInLepton.SetPxPyPzE( 0,0,0,0 ) ;
  fOutLepton.SetPxPyPzE( 0,0,0,0 ) ;
  fInLeptonUnCorr.SetPxPyPzE( 0,0,0,0 ) ;
  fOutLeptonUnCorr.SetPxPyPzE( 0,0,0,0 ) ;
  
  fAnalysisRecord.clear();
}

void EventI::Clear() { 

  fFinalParticles.clear() ; 
  fFinalParticlesUnCorr.clear() ; 
  fAnalysisRecord.clear();
}
