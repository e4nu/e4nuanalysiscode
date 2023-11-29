//____________________________________________________________________________
/*
 * Julia Tena Vidal jtenavidal \at tauex.tau.ac.il 
 */
#include <iostream>
#include "utils/TargetUtils.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"

using namespace e4nu ; 

double utils::GetECalOffset( const unsigned int target_pdg ) {
  double ECalOffset = 0 ; 
  if ( target_pdg == conf::kPdgHe3 ) ECalOffset = conf::kECalOffsetHe3; 
  else if ( target_pdg == conf::kPdgHe4 || target_pdg == conf::kPdgC12 ) ECalOffset = conf::kECalOffsetC12 ;
  else if ( target_pdg == conf::kPdgFe56 ) ECalOffset = conf::kECalOffsetFe56 ;
  return ECalOffset ; 
}

double utils::GetBindingEnergy( const unsigned int target_pdg ) {
  double binding_energy = 0. ; 
  if ( target_pdg == conf::kPdgH ) binding_energy = conf::kBEH;
  else if ( target_pdg == conf::kPdgHe3 ) binding_energy = conf::kBEHe3 - conf::kBED2 + utils::GetECalOffset( target_pdg ) ;
  else if ( target_pdg == conf::kPdgHe4 ) binding_energy = conf::kBEHe4 - conf::kBEHe3 + utils::GetECalOffset( target_pdg ) ;
  else if ( target_pdg == conf::kPdgC12 ) binding_energy = conf::kBEC12 - conf::kBEB + utils::GetECalOffset( target_pdg ) ;
  else if ( target_pdg == conf::kPdgFe56 ) binding_energy = conf::kBEFe56 - conf::kBEMn + utils::GetECalOffset( target_pdg ) ;
  //else if ( target_pdg == conf::kPdgCH2 ) binding_energy = conf::kBEC12 - conf::kBEB ;
  return binding_energy;
}

unsigned int utils::GetTargetNProtons( const unsigned int target_pdg ) {
  if ( target_pdg == conf::kPdgH ) return 1 ;
  else if ( target_pdg == conf::kPdgHe3 || target_pdg == conf::kPdgHe4 ) return 2 ; 
  else if ( target_pdg == conf::kPdgC12 ) return 6 ; 
  else if ( target_pdg == conf::kPdgH ) return 6 ; 
  else if ( target_pdg == conf::kPdgFe56) return 26 ; 
  return 0;
}

unsigned int utils::GetTargetNNeutrons( const unsigned int target_pdg ) { 
  if ( target_pdg == conf::kPdgHe3 ) return 1 ; 
  if ( target_pdg == conf::kPdgHe4 ) return 2 ; 
  else if ( target_pdg == conf::kPdgC12 ) return 6 ; 
  else if ( target_pdg == conf::kPdgFe56 ) return 30 ; 
  return 0;
}

double utils::GetTargetMass( const unsigned int target_pdg ) {
  double mass = 0 ; 
  int n_protons = utils::GetTargetNProtons( target_pdg ) ; 
  int n_neutrons = utils::GetTargetNNeutrons( target_pdg ) ; 
  double BE = utils::GetBindingEnergy( target_pdg ) ; 
 
  mass = n_protons * conf::kProtonMass + n_neutrons * conf::kNeutronMass - BE ; 
  return mass ;
}

double utils::GetResidualTargetMass( const unsigned int target_pdg ) {
  double mass =0 ;

  if ( target_pdg == conf::kPdgHe3 ) mass = conf::kProtonMass + conf::kNeutronMass - conf::kBED2 ; 
  else if ( target_pdg == conf::kPdgHe4 ) mass = conf::kProtonMass + 2*conf::kNeutronMass - conf::kBEHe3 ;
  else if ( target_pdg == conf::kPdgC12 ) mass = 5*conf::kProtonMass + 6*conf::kNeutronMass - conf::kBEB ;
  else if ( target_pdg == conf::kPdgFe56 ) mass = 25*conf::kProtonMass + 30*conf::kNeutronMass - conf::kBEMn ;
  return mass ;
}

unsigned int utils::GetMassNumber( const unsigned int target_pdg ){
  double mass =0 ;

  if ( target_pdg == conf::kPdgHe3 ) mass = 3 ; 
  else if ( target_pdg == conf::kPdgHe4 ) mass = 4 ; 
  else if ( target_pdg == conf::kPdgC12 ) mass = 12 ; 
  else if ( target_pdg == conf::kPdgFe56 ) mass = 56 ;
  return mass ;
}
