#ifndef _KINEMATIC_UTILS_H_
#define _KINEMATIC_UTILS_H_

#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"

namespace e4nu
{
  namespace utils
  {
    double GetECal(const double Ef, const std::map<int, std::vector<TLorentzVector>> particle_map, const int tgt);
    double GetRecoEnu(const TLorentzVector &leptonf, const unsigned int target_pdg);
    double GetQELRecoEnu(const TLorentzVector &leptonf, const unsigned int target_pdg);
    double GetEnergyTransfer(const TLorentzVector &leptonf, const double Ebeam);
    double GetNuECal(const TLorentzVector &leptonf, const double ECal);
    TVector3 GetRecoq3(const TLorentzVector &leptonf, const double EBeam);
    double GetRecoQ2(const TLorentzVector &leptonf, const double EBeam);
    double GetRecoXBJK(const TLorentzVector &leptonf, const double EBeam);
    double GetRecoW(const TLorentzVector &leptonf, const double EBeam);
    TVector3 GetPT(const TVector3 p);
    double DeltaAlphaT(const TVector3 p1, const TVector3 p2);
    double DeltaAlphaT(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons);
    TVector3 DeltaPT(const TVector3 p1, const TVector3 p2);
    TVector3 DeltaPT(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons);
    double DeltaPTx(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons);
    double DeltaPTy(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons);
    double DeltaPhiT(const TVector3 p1, const TVector3 p2);
    double DeltaPhiT(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons);
    double HadSystemMass(const std::map<int, std::vector<TLorentzVector>> hadrons);
    TLorentzVector Missing4Momenta(const double EBeam, const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons, const int tgt);
    double InferedNucleonMom(const double EBeam, const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons, const int tgt);
    double Angle(const TVector3 p1, const TVector3 p2);
    double AngleBeamDir(TVector3 vector);
    TLorentzVector TotHadron(const std::map<int, std::vector<TLorentzVector>> hadrons);
    TLorentzVector VectorInHadFrame(TLorentzVector vector, const std::map<int, std::vector<TLorentzVector>> hadrons);
    TLorentzVector FindParticle(const unsigned int particle_pdg, const std::map<int, std::vector<TLorentzVector>> hadrons);
    double GetAdlerAngleTheta(const double EBeam, const TLorentzVector leptonf, const std::map<int, std::vector<TLorentzVector>> hadrons, const unsigned int particle_pdg);
    double GetAdlerAnglePhi(const double EBeam, const TLorentzVector leptonf, const std::map<int, std::vector<TLorentzVector>> hadrons, const unsigned int particle_pdg);
    double GetRecoEvPionProduction(const TLorentzVector out_electron, const TLorentzVector out_pion);
    double GetRecoWPionProduction(const TLorentzVector out_electron, const TLorentzVector out_pion);
  }
}

#endif
