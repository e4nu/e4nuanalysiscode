#ifndef _KINEMATIC_UTILS_H_
#define _KINEMATIC_UTILS_H_

#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"

namespace e4nu {
  namespace utils {
    double GetECal( const double Ef, const std::map<int,std::vector<TLorentzVector>> particle_map, const int tgt ) ;
    double GetRecoEnu( const TLorentzVector & leptonf, const unsigned int target_pdg ) ;
    double GetQELRecoEnu( const TLorentzVector & leptonf, const unsigned int target_pdg ) ;
    double GetEnergyTransfer( const TLorentzVector & leptonf, const double Ebeam ) ;
    double GetNuECal( const TLorentzVector & leptonf, const double ECal ) ;
    TVector3 GetRecoq3( const TLorentzVector & leptonf, const double EBeam ) ;
    double GetRecoQ2( const TLorentzVector & leptonf, const double EBeam ) ;
    double GetRecoXBJK( const TLorentzVector & leptonf, const double EBeam ) ;
    double GetRecoW( const TLorentzVector & leptonf, const double EBeam ) ;
    TVector3 GetPT( const TVector3 p ) ;
    double DeltaAlphaT( const TVector3 p1 , const TVector3 p2 ) ; 
    double DeltaAlphaT( const TLorentzVector out_electron , const std::map<int,std::vector<TLorentzVector>> hadrons ) ; 
    TVector3 DeltaPT( const TVector3 p1 , const TVector3 p2 ); 
    TVector3 DeltaPT(  const TLorentzVector out_electron , const std::map<int,std::vector<TLorentzVector>> hadrons ) ;
    double DeltaPhiT( const TVector3 p1 , const TVector3 p2 );
    double DeltaPhiT(  const TLorentzVector out_electron , const std::map<int,std::vector<TLorentzVector>> hadrons ) ; 
    double HadSystemMass( const std::map<int,std::vector<TLorentzVector>> hadrons ) ;
  }
}

#endif



