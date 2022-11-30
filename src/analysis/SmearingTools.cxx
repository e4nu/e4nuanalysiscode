//____________________________________________________________________________
/*
 * Julia Tena Vidal jtenavidal \at tauex.tau.ac.il 
 */

#include "analysis/SmearingTools.h"
#include "utils/ParticleUtils.h"

#include <iostream>
#include <cmath>
#include "TRandom3.h"

using namespace e4nu ; 
using namespace e4nu::analysis ; 

bool ApplySmearing( TLorentzVector & particle, const int pdg, const bool is_MC ) {
  if( !is_MC ) return false ;  

  double resolution ;
  if( ! utils::GetParticleResolution( resolution, pdg ) ) return false ; 

  double mass ; 
  if ( ! utils::GetParticleMass( mass, pdg ) ) return false ;
 
  gRandom = new TRandom3();
  double mom = particle.P() ;
  double smear_mom = gRandom->Gaus( mom, resolution * mom ) ; 
  double smear_E = sqrt( std::pow( smear_mom, 2 ) + std::pow( mass, 2 ) ); 

  // Apply smearing
  particle.SetPxPyPzE( smear_mom/mom * particle.Px(), smear_mom / mom * particle.Py(), smear_mom / mom * particle.Pz(), smear_E ) ; 

  particle.SetPhi( particle.Phi() + TMath::Pi() ) ; // Vec.Phi() is between (-180,180), GENIE coordinate system flipped with respect to CLAS  

  return true; 
}
