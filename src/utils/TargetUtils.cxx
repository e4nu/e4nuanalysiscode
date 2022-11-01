//____________________________________________________________________________
/*
 * Julia Tena Vidal jtenavidal \at tauex.tau.ac.il 
 */

#include "src/utils/TargetUtils.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"

using namespace e4nu ; 

unsigned double utils::GetECalOffset( const unsigned int target_pdg ) {
  double ECalOffset = 0 ; 
  if ( target_pdg == kPdgHe3 ) ECalOffset = ParticleI::kECalOffsetHe3; 
  else if ( target_pdg == kPdgHe4 || target_pdg == kPdgC12 ) ECalOffset = ParticleI::kECalOffsetC12 ;
  else if ( target_pdg == kPdgFe56 ) ECalOffset = ParticleI::kECalOffsetFe56 ;
  return ECalOffset ; 
}

unsigned double utils::GetBindingEnergy( const unsigned int target_pdg ) {
  double binding_energy = 0. ; 
  if ( target_pdg == kPdgHe3 ) binding_energy = TargetI::kBEHe3 - TargetI::kBED2 + utils::GetECalOffset( target_pdg ) ;
  else if ( target_pdg == kPdgHe4 ) binding_energy = TargetI::kBEHe4 - TargetI::kBEHe3 + utils::GetECalOffset( target_pdg ) ;
  else if ( target_pdg == kPdgC12 ) binding_energy = TargetI::kBEC12 - TargetI::kBEB + utils::GetECalOffset( target_pdg ) ;
  else if ( target_pdg == kPdgFe56 ) binding_energy = TargetI::kBEFe54 - TargetI::kBEMn + utils::GetECalOffset( target_pdg ) ;
  else if ( target_pdg == kPdgCH2 ) binding_energy = TargetI::kBEC12 - TargetI::kBEB ;
  return binding_energy;
}

unsigned int utils::GetTargetNProtons( const unsigned int target_pdg ) {
  
  if ( target_pdg == kPdgHe3 || target_pdg == kPdgHe4 ) return 2 ; 
  else if ( target_pdg == kPdgC12 || target_pdg == kPdgCH2 ) return 6 ; 
  else if ( target_pdg == kPdgFe56) return 26 ; 
  return 0;
}

unsigned int utils::GetTargetNNeutrons( const unsigned int target_pdg ) { 
  if ( target_pdg == kPdgHe3 ) return 1 ; 
  if ( target_pdg == kPdgHe4 ) return 2 ; 
  else if ( target_pdg == kPdgC12 || target_pdg == kPdgCH2 ) return 6 ; 
  else if ( target_pdg == kPdgFe56) return 30 ; 
  return 0;
}

unsigned double utils::GetTargetMass( const unsigned int target_pdg ) {
  double mass = 0 ; 
  int n_protons = utils::GetTargetNProtons( target_mass ) ; 
  int n_neutrons = utils::GetTargetNNeutrons( target_mass ) ; 
  double BE = utils::GetBindingEnergy( target_pdg ) ; 

  if ( target_pdg == kPdgCH2 ){ BE = utils::GetBindingEnergy( kPdgC12 ) ; } // IS THIS CORRECT ? 
 
  mass = n_protons * ParticleI::kMassProton + n_neutrons * ParticleI::kMassNeutron - BE ; 
  return mass ;
}

unsigned double utils::GetResidualTargetMass( const unsigned int target_pdg ) {
  double mass =0 ;

  if ( target_pdg == kPdgHe3 ) mass = ParticleI::kMassProton + ParticleI::kMassNeutron - TargetI::KBED2 ; 
  else if ( target_pdg == kPdgHe4 ) mass = ParticleI::kMassProton + 2*ParticleI::kMassNeutron - TargetI::KBEHe3 ;
  else if ( target_pdg == kPdgC12 ) mass = 5*ParticleI::kMassProton + 6*ParticleI::kMassNeutron - TargetI::KBEB ;
  else if ( target_pdg == kPdgFe56 ) mass = 25*ParticleI::kMassProton + 30*ParticleI::kMassNeutron - TargetI::KBEMn ;
  else if ( target_pdg == kPdgCH2 ) mass = 25*ParticleI::kMassProton + 30*ParticleI::kMassNeutron - TargetI::KBEMn ;                                                                                                                            
  return mass ;
}


