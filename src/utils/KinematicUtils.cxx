/**                                                                                                                                                                                                                   * This file contains utils specific for particles                                                                                                                                                                    * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                   * \date October 2022                                                                                                                                                                                                 **/
#include <iostream>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "src/utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "src/utils/TargetUtils.h"

using namespace e4nu ;
using namespace e4nu::utils;


double GetQELRecoEnu( const TLorentzVector & leptonf, const unsigned int target_pdg ) {
  double Ereco = 0 ; 
  double E = leptonf.E();
  double P = leptonf.P();
  double CosTh = TMath::Cos( leptonf.Theta() ); 
  double BE = utils::GetBindingEnergy( target_pdg ) ;

  Ereco = std::pow( conf::kProtonMass,2) - pow( conf::kElectronMass,2) + 2*E*( conf::kNeutronMass - BE ) - std::pow( conf::kNeutronMass -  BE ,2) ;
  Ereco /= ( conf::kNeutronMass - BE ) - E + P * CosTh ; 
  Ereco /= 2. ;

  return Ereco ; 
} 

double GetEnergyTransfer( const TLorentzVector & leptonf, const double Ebeam ) {
  return Ebeam - leptonf.E() ; 
}

TVector3 GetRecoq3( const TLorentzVector & leptonf, const double EBeam ) {
  TLorentzVector beam ( 0,0,EBeam,EBeam) ;
  return ( leptonf - beam ).Vect() ;
}

double GetRecoQ2( const TLorentzVector & leptonf, const double EBeam ) {
  TLorentzVector beam ( 0,0,EBeam,EBeam) ;
  return -( leptonf - beam ).Mag2() ; 
}

double GetRecoXBJK( const TLorentzVector & leptonf, const double EBeam ) {
  double nu = utils::GetEnergyTransfer( leptonf, EBeam ) ; 
  double Q2 = utils::GetRecoQ2( leptonf, EBeam ) ; 
  return Q2 / ( 2 * conf::kProtonMass * nu ) ;
}

double GetRecoW( const TLorentzVector & leptonf, const double EBeam ) {
  double mp = conf::kProtonMass ; 
  double nu = utils::GetEnergyTransfer( leptonf, EBeam ) ; 
  TVector3 q3 = utils::GetRecoq3( leptonf, EBeam ) ; 
  double W2 = std::pow( mp + nu, 2 ) - q3.Mag2() ; 
  if ( W2 < 0 ) return 0 ;
  return TMath::Sqrt( W2 ) ; 
} 

double GetXSecScale( const TLorentzVector & leptonf, const double EBeam, const bool is_electron ) {
  double scale = 1 ; 
  if ( is_electron ) {
    double reco_Q2 = utils::GetRecoQ2( leptonf, EBeam ) ;
    scale /= std::pow( reco_Q2, 2 ) ; 
  }
  // Add additional scalings ?
  return scale ; 
}

TVector3 GetPT( const TVector3 p, const double EBeam ) {
  TLorentzVector beam ( 0,0,EBeam,EBeam) ;
  TVector3 beam_dir = beam.Vect().Unit();
  double vect_parallel = p.Dot(beam_dir);

  // Calculate transverse vector:
  TVector3 vect_T = p - (vect_parallel * beam_dir ) ; 
  return vect_T ; 
}

double DeltaAlphaT( const TVector3 p1 , const TVector3 p2, const double EBeam ) {
  TVector3 P1T_dir = utils::GetPT(p1,EBeam).Unit();
  TVector3 DeltaPT_dir = utils::DeltaPT(p1, p2, EBeam).Unit();

  return acos(-P1T_dir.Dot(DeltaPT_dir));
}

TVector3 DeltaPT( const TVector3 p1 , const TVector3 p2, const double EBeam ) {
  TVector3 P1_T = utils::GetPT(p1, EBeam);
  TVector3 P2_T = utils::GetPT(p2, EBeam);

  return P1_T + P2_T;
}

double DeltaPhiT( const TVector3 p1 , const TVector3 p2, const double EBeam ) {
  TVector3 P1T_dir = utils::GetPT(p1, EBeam).Unit();
  TVector3 P2T_dir = utils::GetPT(p2, EBeam).Unit();
  return acos(-P1T_dir.Dot(P2T_dir)); 
}





