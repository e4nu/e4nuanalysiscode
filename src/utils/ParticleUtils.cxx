//____________________________________________________________________________
/*
 * Julia Tena Vidal jtenavidal \at tauex.tau.ac.il 
 */

#include "src/utils/ParticleUtils.h"
#include "conf/ParticleI.h"

using namespace e4nu ; 

unsigned double utils::GetParticleResolucion( int particle_pdg, double Ebeam, bool apply_resolution = true ) {
  double resolution = 0 ; 
  if ( !apply_resolution ) return resolution ; 

  if ( particle_pdg == kPdgProton ) resolution = ParticleI::kProtonRes ; 
  else if ( particle_pdg == kPdgElectron ) resolution = ParticleI::kElectronRes ; 
  else if ( particle_pdg == kPdgPiM || particle_pdg == kPdgPim || particle_pdg == kPdgPi0 ) resolution = ParticleI::kPionRes ; // also pi0?

  if ( Ebeam == 1.161 ) resolution *= 3; // Is it only this value or beam_E>1.1 GeV ? 
  return resolution ; 
}
