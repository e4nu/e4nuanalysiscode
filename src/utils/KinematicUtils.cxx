#include <iostream>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "utils/TargetUtils.h"
#include "utils/ParticleUtils.h"
#include "utils/TargetUtils.h"

using namespace e4nu;

double utils::GetECal(const double Ef, const std::map<int, std::vector<TLorentzVector>> particle_map, const int tgt)
{
  double ECal = Ef; // Add energy of outgoing lepton
  for (auto it = particle_map.begin(); it != particle_map.end(); ++it)
  {
    // Calculate ECal for visible particles
    for (unsigned int i = 0; i < (it->second).size(); ++i)
    {
      if (it->first == 11)
        continue;
      ECal += (it->second)[i].E(); // Add Kinetic energy of hadrons
      if (it->first == conf::kPdgProton)
        ECal += utils::GetBindingEnergy(tgt) - utils::GetParticleMass(it->first); // Correct for proton binding energy
    }
  }

  return ECal;
}

double utils::GetRecoEnu(const TLorentzVector &leptonf, const unsigned int target_pdg)
{
  return (conf::kProtonMass * utils::GetBindingEnergy(target_pdg) + conf::kProtonMass * leptonf.E()) / (conf::kProtonMass - leptonf.E() + leptonf.Rho() * TMath::Cos(leptonf.Theta()));
}

double utils::GetQELRecoEnu(const TLorentzVector &leptonf, const unsigned int target_pdg)
{
  double Ereco = 0;
  double E = leptonf.E();
  double P = leptonf.P();
  double CosTh = TMath::Cos(leptonf.Theta());
  double BE = utils::GetBindingEnergy(target_pdg);

  Ereco = std::pow(conf::kProtonMass, 2) - pow(conf::kElectronMass, 2) + 2 * E * (conf::kNeutronMass - BE) - std::pow(conf::kNeutronMass - BE, 2);
  Ereco /= (conf::kNeutronMass - BE) - E + P * CosTh;
  Ereco /= 2.;

  return Ereco;
}

double utils::GetEnergyTransfer(const TLorentzVector &leptonf, const double Ebeam)
{
  return Ebeam - leptonf.E();
}

double utils::GetNuECal(const TLorentzVector &leptonf, const double ECal)
{
  TLorentzVector nu(0, 0, ECal, ECal);
  return -(leptonf - nu).E();
}

TVector3 utils::GetRecoq3(const TLorentzVector &leptonf, const double EBeam)
{
  TLorentzVector beam(0, 0, EBeam, EBeam);
  return (beam - leptonf).Vect();
}

double utils::GetRecoQ2(const TLorentzVector &leptonf, const double EBeam)
{
  TLorentzVector beam(0, 0, EBeam, EBeam);
  return -(leptonf - beam).Mag2();
}

double utils::GetRecoXBJK(const TLorentzVector &leptonf, const double EBeam)
{
  double nu = utils::GetEnergyTransfer(leptonf, EBeam);
  double Q2 = utils::GetRecoQ2(leptonf, EBeam);
  return Q2 / (2 * conf::kProtonMass * nu);
}

double utils::GetRecoW(const TLorentzVector &leptonf, const double EBeam)
{
  double mp = conf::kProtonMass;
  double nu = utils::GetEnergyTransfer(leptonf, EBeam);
  TVector3 q3 = utils::GetRecoq3(leptonf, EBeam);
  double W2 = std::pow(mp + nu, 2) - q3.Mag2();
  if (W2 < 0)
    return 0;
  return TMath::Sqrt(W2);
}

TVector3 utils::GetPT(const TVector3 p)
{
  TVector3 beam_dir(0, 0, 1);
  double vect_parallel = p.Dot(beam_dir);

  // Calculate transverse vector:
  TVector3 vect_T = p - (vect_parallel * beam_dir);
  return vect_T;
}

double utils::DeltaAlphaT(const TVector3 p1 /*out electron*/, const TVector3 p2 /*proton*/)
{
  TVector3 P1T_dir = utils::GetPT(p1).Unit();
  TVector3 DeltaPT_dir = utils::DeltaPT(p1, p2).Unit();

  return acos(-P1T_dir.Dot(DeltaPT_dir)) * TMath::RadToDeg();
}

double utils::DeltaAlphaT(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TVector3 P1T_dir = utils::GetPT(out_electron.Vect()).Unit();
  TVector3 DeltaPT_dir = utils::DeltaPT(out_electron, hadrons).Unit();

  return acos(-P1T_dir.Dot(DeltaPT_dir)) * TMath::RadToDeg();
}

TVector3 utils::DeltaPT(const TVector3 p1, const TVector3 p2)
{
  TVector3 P1_T = utils::GetPT(p1);
  TVector3 P2_T = utils::GetPT(p2);

  return P1_T + P2_T;
}

TVector3 utils::DeltaPT(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TVector3 P1_T = utils::GetPT(out_electron.Vect());
  TLorentzVector tot_hadron = utils::TotHadron(hadrons);
  TVector3 P2_T = utils::GetPT(tot_hadron.Vect());

  return P1_T + P2_T;
}

double utils::DeltaPTx(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TVector3 deltapt = utils::DeltaPT(out_electron, hadrons);
  double dalphat = utils::DeltaAlphaT(out_electron, hadrons);
  return deltapt.Mag() * TMath::Sin(dalphat);
}

double utils::DeltaPTy(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TVector3 deltapt = utils::DeltaPT(out_electron, hadrons);
  double dalphat = utils::DeltaAlphaT(out_electron, hadrons);
  return deltapt.Mag() * TMath::Cos(dalphat);
}

double utils::DeltaPhiT(const TVector3 p1, const TVector3 p2)
{
  TVector3 P1T_dir = utils::GetPT(p1).Unit();
  TVector3 P2T_dir = utils::GetPT(p2).Unit();
  return acos(-P1T_dir.Dot(P2T_dir)) * TMath::RadToDeg();
}

double utils::DeltaPhiT(const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TVector3 P1T_dir = utils::GetPT(out_electron.Vect()).Unit();
  TLorentzVector tot_hadron = utils::TotHadron(hadrons);
  TVector3 P2T_dir = utils::GetPT(tot_hadron.Vect()).Unit();
  return acos(-P1T_dir.Dot(P2T_dir)) * TMath::RadToDeg();
}

double utils::HadSystemMass(const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TLorentzVector tot_hadron = utils::TotHadron(hadrons);
  return tot_hadron.Mag();
}

TLorentzVector utils::Missing4Momenta(const double EBeam, const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons, const int tgt)
{
  // k^mu - k'^mu = p_p^mu + p_pi^mu + p_miss^{mu} - (Mp-Be, 0)
  TLorentzVector beam(0, 0, EBeam, EBeam);
  TLorentzVector q = beam - out_electron;

  // Initial nucleon at rest
  TLorentzVector in_nucleon(0, 0, 0, utils::GetParticleMass(conf::kPdgProton)-utils::GetBindingEnergy(tgt));
  TLorentzVector tot_hadron;
  for (auto it = hadrons.begin(); it != hadrons.end(); ++it)
  {
    for (unsigned int i = 0; i < (it->second).size(); ++i)
    {
      if( it->first == conf::kPdgProton ) tot_hadron += (it->second)[i] - in_nucleon ;
      else tot_hadron += (it->second)[i];
    }
  }

  return ( q - tot_hadron );
}

double utils::InferedNucleonMom(const double EBeam, const TLorentzVector out_electron, const std::map<int, std::vector<TLorentzVector>> hadrons, const int tgt)
{
  double BindingE = utils::GetBindingEnergy(tgt);
  double Mn = 0.5 * (conf::kProtonMass + conf::kNeutronMass);
  double MA = utils::GetTargetMass(tgt);
  double Ep = out_electron.E();
  TLorentzVector tot_hadron = utils::TotHadron(hadrons);
  double E_had = tot_hadron.E();
  TLorentzVector beam(0, 0, EBeam, EBeam);
  TVector3 beam_dir = beam.Vect().Unit();
  double pl_e = (out_electron.Vect()).Dot(beam_dir);
  double pl_had = (tot_hadron.Vect()).Dot(beam_dir);
  double R = MA + pl_e + pl_had - Ep - E_had;
  double MA_f = MA - Mn + BindingE;
  TVector3 vect_dpt = utils::DeltaPT(out_electron, hadrons);
  double dpl = 0.5 * (R - (pow(MA_f, 2) + vect_dpt.Mag2()) / R);
  double pn = sqrt(vect_dpt.Mag2() + pow(dpl, 2));
  return pn;
}

double utils::Angle(const TVector3 p1, const TVector3 p2)
{
  return p1.Angle(p2);
}

double utils::AngleBeamDir(TVector3 vector)
{
  TVector3 beam_dir(0, 0, 1);
  return vector.Angle(beam_dir);
}

TLorentzVector utils::TotHadron(const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TLorentzVector tot_hadron;
  for (auto it = hadrons.begin(); it != hadrons.end(); ++it)
  {
    for (unsigned int i = 0; i < (it->second).size(); ++i)
    {
      tot_hadron += (it->second)[i];
    }
  }
  return tot_hadron;
}

TLorentzVector utils::VectorInHadFrame(TLorentzVector vector, const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TLorentzVector tot_hadron = utils::TotHadron(hadrons);
  vector.Boost(-tot_hadron.BoostVector());
  return vector;
}

TLorentzVector utils::FindParticle(const unsigned int particle_pdg, const std::map<int, std::vector<TLorentzVector>> hadrons)
{
  TLorentzVector particle(0, 0, 0, 0);
  for (auto it = hadrons.begin(); it != hadrons.end(); ++it)
  {
    // Find particle pdg
    if ((it->second).size() != 1)
      continue;
    if (abs(it->first) == particle_pdg)
    {
      for (unsigned int i = 0; i < (it->second).size(); ++i)
      {
        particle = (it->second)[i];
      }
    }
  }
  return particle;
}

double utils::GetAdlerAngleTheta(const double EBeam, const TLorentzVector leptonf, const std::map<int, std::vector<TLorentzVector>> hadrons, const unsigned int particle_pdg)
{
  // https://arxiv.org/pdf/1511.00501.pdf
  TLorentzVector beam(0, 0, EBeam, EBeam);
  TLorentzVector q = beam - leptonf;
  TLorentzVector particle = utils::FindParticle(particle_pdg, hadrons);

  // In tot hadron reference frame
  TVector3 z_axis = utils::VectorInHadFrame(q, hadrons).Vect().Unit();
  TVector3 particle_dir_hadframe = utils::VectorInHadFrame(particle, hadrons).Vect().Unit();

  return particle_dir_hadframe.Angle(z_axis); // theta
}

double utils::GetAdlerAnglePhi(const double EBeam, const TLorentzVector leptonf, const std::map<int, std::vector<TLorentzVector>> hadrons, const unsigned int particle_pdg)
{
  // https://arxiv.org/pdf/1511.00501.pdf

  TLorentzVector particle = utils::FindParticle(particle_pdg, hadrons);
  TLorentzVector beam(0, 0, EBeam, EBeam);
  TLorentzVector q = beam - leptonf;

  // Boost to tot hadron frame
  TVector3 particle_dir = utils::VectorInHadFrame(particle, hadrons).Vect();
  TVector3 leptonf_hadframe = utils::VectorInHadFrame(leptonf, hadrons).Vect();

  // Axis definition
  TVector3 z_axis = utils::VectorInHadFrame(q, hadrons).Vect().Unit();
  TVector3 y_axis = z_axis.Cross(leptonf_hadframe).Unit();
  TVector3 x_axis = y_axis.Cross(z_axis);

  TVector3 particle_perp(z_axis.Cross(particle_dir.Cross(z_axis).Unit()));

  return particle_perp.Angle(x_axis); // phi
}

double utils::GetRecoEvPionProduction(const TLorentzVector out_electron, const TLorentzVector out_pion)
{
  double elMass2 = pow(conf::kElectronMass, 2);
  double piMMass2 = pow(conf::kPiMMass, 2);
  double nucMass = conf::kNucleonMass;

  double elEnergy = out_electron.E();
  double piMEnergy = out_pion.E();
  double E_rec = 0;

  if (elEnergy > 0 && piMEnergy > 0)
  {
    double elMag = out_electron.Vect().Mag();
    double piMMag = out_pion.Vect().Mag();

    TLorentzVector beam(0, 0, 1, 1);
    double elAngle = TMath::Cos(out_electron.Vect().Angle(beam.Vect()));
    double piMAngle = TMath::Cos(out_pion.Vect().Angle(beam.Vect()));

    E_rec = 2 * elMass2 + piMMass2 - 2 * nucMass * (elEnergy + piMEnergy) + 2 * out_electron.Dot(out_pion);
    E_rec /= (2 * (elEnergy + piMEnergy - piMMag * piMAngle - elMag * elAngle - nucMass));
  }
  if (E_rec < 0)
    return 0;
  return E_rec;
}

double utils::GetRecoWPionProduction(const TLorentzVector out_electron, const TLorentzVector out_pion)
{
  double E_rec = utils::GetRecoEvPionProduction(out_electron, out_pion);

  TLorentzVector ev_4vect_rec;
  ev_4vect_rec.SetPxPyPzE(0, 0, E_rec, E_rec);

  TLorentzVector nucl_4vect;
  nucl_4vect.SetPxPyPzE(0, 0, 0, conf::kNucleonMass);

  TLorentzVector q_4vect = ev_4vect_rec - out_electron;

  double W_rec = (nucl_4vect + q_4vect).Mag();

  return W_rec;
}
