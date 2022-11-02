//____________________________________________________________________________
/*
 * Julia Tena Vidal jtenavidal \at tauex.tau.ac.il 
 */

#include "src/utils/ParticleUtils.h"
#include "conf/ParticleI.h"

using namespace e4nu ; 

double utils::GetParticleResolucion( const int particle_pdg, const double Ebeam, const bool apply_resolution ) {
  double resolution = 0 ; 
  if ( !apply_resolution ) return resolution ; 

  if ( particle_pdg == conf::kPdgProton ) resolution = conf::kProtonRes ; 
  else if ( particle_pdg == conf::kPdgElectron ) resolution = conf::kElectronRes ; 
  else if ( particle_pdg == conf::kPdgPiP || particle_pdg == conf::kPdgPiM || particle_pdg == conf::kPdgPi0 ) resolution = conf::kPionRes ; // also pi0?

  if ( Ebeam == 1.161 ) resolution *= 3; // Is it only this value or beam_E>1.1 GeV ? 
  return resolution ; 
}
