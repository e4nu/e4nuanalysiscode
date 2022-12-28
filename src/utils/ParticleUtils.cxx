//____________________________________________________________________________
/*
 * Julia Tena Vidal jtenavidal \at tauex.tau.ac.il 
 */

#include "utils/ParticleUtils.h"
#include "conf/ParticleI.h"
#include <TMath.h>
#include <TRandom3.h>

using namespace e4nu ; 

double utils::GetParticleResolucion( const int particle_pdg, const double Ebeam ) {
  double resolution = 0 ; 

  if ( particle_pdg == conf::kPdgProton ) resolution = conf::kProtonRes ; 
  else if ( particle_pdg == conf::kPdgElectron ) resolution = conf::kElectronRes ; 
  else if ( particle_pdg == conf::kPdgPiP || particle_pdg == conf::kPdgPiM || particle_pdg == conf::kPdgPi0 ) resolution = conf::kPionRes ; // also pi0?

  if ( Ebeam == 1.161 ) resolution *= 3; // Is it only this value or beam_E>1.1 GeV ? 

  return resolution ; 
}

double utils::GetParticleMass( const int pdg ) {
  double mass = 0 ; 
  if( pdg == conf::kPdgProton ) mass = conf::kProtonMass ; 
  else if ( pdg == conf::kPdgElectron ) mass = conf::kElectronMass ; 
  else if ( pdg == conf::kPdgPiP ) mass = conf::kPiPMass ; 
  else if ( pdg == conf::kPdgPiM ) mass = conf::kPiMMass ; 
  else if ( pdg == conf::kPdgNeutron ) mass = conf::kNeutronMass ;

  return mass ; 
}

void utils::ApplyResolution( const int pdg, TLorentzVector & mom, const double EBeam ) {
  double res = utils::GetParticleResolucion( pdg, EBeam ) ;
  double p = mom.P() ;
  double M = GetParticleMass( pdg ) ;

  gRandom = new TRandom3() ; 
  gRandom->SetSeed(10);

  double SmearedPe = gRandom->Gaus(p,res*p);
  double SmearedE = sqrt( pow( SmearedPe,2 ) + pow( M,2 ) ) ; 

  mom.SetPxPyPzE( SmearedPe/p * mom.Px(), SmearedPe/p *mom.Py(), SmearedPe/p *mom.Pz(), SmearedE ) ; 
  delete gRandom ; 
}

int utils::GetParticleCharge( const int pdg ) {
  if( pdg == conf::kPdgElectron ) return conf::kElectronCharge ; 
  else if( pdg == conf::kPdgProton ) return conf::kProtonCharge ; 
  else if ( pdg == conf::kPdgPiP ) return conf::kPiPCharge ; 
  else if ( pdg == conf::kPdgPiM ) return conf::kPiMCharge ; 
  else if ( pdg == conf::kPdgPi0 ) return conf::kPi0Charge ;
  else if ( pdg == conf::kPdgNeutron ) return conf::kNeutronCharge ; 
  else if ( pdg == conf::kPdgPhoton ) return conf::kPhotonCharge ; 
  return 0 ; 
}
