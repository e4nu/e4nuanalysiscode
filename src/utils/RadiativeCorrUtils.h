#ifndef _RADUTILS_H_
#define _RADUTILS_H_
#include "TLorentzVector.h"
#include <string>
namespace e4nu {
  namespace utils
  {
    // QEL
    double RadCorrQELVertex( const double Q2 ) ; 
    double RadCorrQELVacumm( const double Q2 ) ; 
    double RadCorrQELRealRad( const double Q2, const double E, const double Ep, const double theta) ; 
    double SIMCRadCorrQELRadInElectron( const double EBeam, const TLorentzVector electron_vertex, const int tgt, const double thickness, const double max_Ephoton );
    double SIMCRadCorrQELRadOutElectron( const TLorentzVector electron_vertex, TLorentzVector & out_electron, const double EBeam, const double Q2, const int tgt, const double thickness, const double max_Ephoton, const std::string model ) ;
    // Probability funcitons
    double SIMCBFactor( const double tgt_pdg ) ;
    double SIMCEnergyLoss( const double EBeam, const TLorentzVector particle, const int p_pdg, const double tgt_pdg, const double thickness, const double max_Ephoton ) ;
    double SimpleEnergyLoss(const double EBeam, const double tgt_pdg, const double thickness, const double max_Ephoton ) ;
    double VanderhagenELoss( const double Q2 , const double Ee ) ;

    // Compute TLorentzVector for emited photon
    TLorentzVector GetEmittedHardPhoton( TLorentzVector electron, double eloss ) ; 
  }
}

#endif
