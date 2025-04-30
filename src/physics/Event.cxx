// _______________________________________________
/*
 * Event Interface base class
 *
 */
#include <iostream>
#include "physics/Event.h"
#include "conf/ParticleI.h"
#include "utils/DetectorUtils.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"

using namespace e4nu ;

Event::Event() {
  this->Initialize() ;
}

Event::~Event() {
  this->Clear();
}

void Event::SetOutLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fOutLepton.SetPxPyPzE( px, py, pz, E ) ;
  return ;
}

void Event::SetInLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fInLepton.SetPxPyPzE( px, py, pz, E ) ;
  return ;
}

void Event::SetFinalParticle( const int pdg, const double E, const double px, const double py, const double pz ) {
  TLorentzVector mom;
  mom.SetPxPyPzE( px, py, pz, E ) ;
  if( fFinalParticles.find(pdg) == fFinalParticles.end() ) {
    std::vector<TLorentzVector> vct ;
    vct.push_back(mom);
    fFinalParticles.insert( std::pair<int,std::vector<TLorentzVector>>(pdg, vct) ) ;
  } else {
    fFinalParticles[pdg].push_back( mom ) ;
  }
}

void Event::SetOutUnCorrLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fOutLeptonUnCorr.SetPxPyPzE( px, py, pz, E ) ;
  return ;
}

void Event::SetInUnCorrLeptonKinematics( const double E, const double px, const double py, const double pz ) {
  fInLeptonUnCorr.SetPxPyPzE( px, py, pz, E ) ;
  return ;
}

void Event::SetFinalParticleUnCorr( const int pdg, const double E, const double px, const double py, const double pz ) {
  TLorentzVector mom;
  mom.SetPxPyPzE( px, py, pz, E ) ;
  if( fFinalParticlesUnCorr.find(pdg) == fFinalParticlesUnCorr.end() ) {
    std::vector<TLorentzVector> vct ;
    vct.push_back(mom);
    fFinalParticlesUnCorr.insert( std::pair<int,std::vector<TLorentzVector>>(pdg, vct) ) ;
  } else {
    fFinalParticlesUnCorr[pdg].push_back( mom ) ;
  }
}

void Event::StoreAnalysisRecord( unsigned int analysis_step ) {
  double weight = this->GetTotalWeight() ;
  std::map<int,std::vector<TLorentzVector>> part_map = this->GetFinalParticles4Mom() ;
  std::vector<int> pdg_list ;
  for( auto it = part_map.begin() ; it != part_map.end() ; ++it ) {
    for( unsigned int i = 0 ; i < part_map[it->first].size() ; ++i ) pdg_list.push_back( it->first ) ;
  }
  std::pair<std::vector<int>,double> pair ( pdg_list, weight ) ;
  fAnalysisRecord[analysis_step] = pair ;
  return ;
}

unsigned int Event::GetEventMultiplicity( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) const {
  // return number of charged particles in event
  unsigned int multiplicity = 0 ;
  for( auto it = hadronic_system.begin() ; it != hadronic_system.end() ; ++it ) {
    multiplicity += std::abs( utils::GetParticleCharge( it->first ) );
  }
  return multiplicity ;
}

unsigned int Event::GetNSignalParticles( std::map<int,std::vector<TLorentzVector>> hadronic_system, const std::map<int,unsigned int> topology ) const {
  unsigned int N_signal = 0 ;
  for( auto it = hadronic_system.begin() ; it != hadronic_system.end() ; ++it ) {
    if( it->first == conf::kPdgElectron ) continue ;
    if( topology.find(it->first) != topology.end() ) {
      N_signal += hadronic_system[it->first].size() ;
    }
  }

  return N_signal ;
}

int Event::GetEventTotalVisibleCharge( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) const {
  // return number of charged particles in event
  unsigned int charge = 0 ;
  for( auto it = hadronic_system.begin() ; it != hadronic_system.end() ; ++it ) {
    charge += utils::GetParticleCharge( it->first ) ;
  }
  return charge ;
}

double Event::GetObservable( const std::string observable ) {

  int ID = GetEventID();
  int TargetPdg = GetTargetPdg();
  int InLeptonPdg = GetInLeptPdg();
  int OutLeptonPdg = GetOutLeptPdg();
  double BeamE = GetInLepton4Mom().E();

  TLorentzVector out_mom = GetOutLepton4Mom();
  double Efl = out_mom.E();
  double pfl = out_mom.P();
  double pflx = out_mom.Px();
  double pfly = out_mom.Py();
  double pflz = out_mom.Pz();
  double pfl_theta = out_mom.Theta() * TMath::RadToDeg();
  out_mom.SetPhi(out_mom.Phi() + TMath::Pi());
  unsigned int ElectronSector = utils::GetSector(out_mom.Phi());
  double pfl_phi = out_mom.Phi() * TMath::RadToDeg();

  double RecoQELEnu = utils::GetQELRecoEnu(out_mom, TargetPdg);
  double RecoEnergyTransfer = utils::GetEnergyTransfer(out_mom, BeamE);
  double Recoq3 = utils::GetRecoq3(out_mom, BeamE).Mag();
  double RecoQ2 = utils::GetRecoQ2(out_mom, BeamE);
  double RecoXBJK = utils::GetRecoXBJK(out_mom, BeamE);
  double RecoW = utils::GetRecoW(out_mom, BeamE);

  // flip hadrons phi
  std::map<int, std::vector<TLorentzVector>> hadron_map = GetFinalParticles4Mom();
  for (auto it = hadron_map.begin(); it != hadron_map.end(); ++it)
  {
    for (unsigned int i = 0; i < (it->second).size(); ++i)
    {
      TLorentzVector had_4mom = (it->second)[i];
      had_4mom.SetPhi((it->second)[i].Phi() + TMath::Pi());
      (it->second)[i] = had_4mom;
    }
  }

  TLorentzVector p_max(0, 0, 0, 0);
  double max_mom = 0;
  for (unsigned int i = 0; i < hadron_map[conf::kPdgProton].size(); ++i)
    {
      if (hadron_map[conf::kPdgProton][i].P() > max_mom)
	{
	  max_mom = hadron_map[conf::kPdgProton][i].P();
	  p_max = hadron_map[conf::kPdgProton][i];
	}
    }

  double HadronsAngle = 0;
  std::vector<TLorentzVector> particles;
  for (auto it = hadron_map.begin(); it != hadron_map.end(); ++it)
  {
    if ((it->second).size() != 1)
      continue;
    for (unsigned int i = 0; i < (it->second).size(); ++i)
    {
      particles.push_back((it->second)[i]);
    }
  }
  if (particles.size() == 2)
  {
    HadronsAngle = utils::Angle(particles[0].Vect(), particles[1].Vect()) * TMath::RadToDeg();
  }

  double proton_E = p_max.E();
  double proton_mom = p_max.P();
  double proton_momx = p_max.Px();
  double proton_momy = p_max.Py();
  double proton_momz = p_max.Pz();
  double proton_theta = p_max.Theta() * TMath::RadToDeg();
  double proton_phi = p_max.Phi() * TMath::RadToDeg();
  double ECal = utils::GetECal(out_mom.E(), GetFinalParticles4Mom(), TargetPdg);
  double DiffECal = utils::GetECal(out_mom.E(), GetFinalParticles4Mom(), TargetPdg) - BeamE;
  double AlphaT = utils::DeltaAlphaT(out_mom.Vect(), p_max.Vect());
  double DeltaPT = utils::DeltaPT(out_mom.Vect(), p_max.Vect()).Mag();
  double DeltaPhiT = utils::DeltaPhiT(out_mom.Vect(), p_max.Vect());
  double HadAlphaT = utils::DeltaAlphaT(out_mom, hadron_map);
  double HadDeltaPT = utils::DeltaPT(out_mom, hadron_map).Mag();
  double HadDeltaPTx = utils::DeltaPTx(out_mom, hadron_map);
  double HadDeltaPTy = utils::DeltaPTy(out_mom, hadron_map);
  double HadDeltaPhiT = utils::DeltaPhiT(out_mom, hadron_map);
  double InferedNucleonMom = utils::InferedNucleonMom(BeamE, out_mom, hadron_map, TargetPdg);

  TLorentzVector pip_max(0, 0, 0, 0);
  max_mom = 0;
  for (unsigned int i = 0; i < hadron_map[conf::kPdgPiP].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiP][i].P() > max_mom)
	{
	  max_mom = hadron_map[conf::kPdgPiP][i].P();
	  pip_max = hadron_map[conf::kPdgPiP][i];
	}
    }

  double pip_E = pip_max.E();
  double pip_mom = pip_max.P();
  double pip_momx = pip_max.Px();
  double pip_momy = pip_max.Py();
  double pip_momz = pip_max.Pz();
  double pip_theta = pip_max.Theta() * TMath::RadToDeg();
  double pip_phi = pip_max.Phi() * TMath::RadToDeg();

  TLorentzVector pim_max(0, 0, 0, 0);
  max_mom = 0;
  for (unsigned int i = 0; i < hadron_map[conf::kPdgPiM].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiM][i].P() > max_mom)
	{
	  max_mom = hadron_map[conf::kPdgPiM][i].P();
	  pim_max = hadron_map[conf::kPdgPiM][i];
	}
    }


  // Adler angles
  double AdlerAngleThetaP = utils::GetAdlerAngleTheta(BeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  double AdlerAnglePhiP = utils::GetAdlerAnglePhi(BeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  double AdlerAngleThetaPi = utils::GetAdlerAngleTheta(BeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg();
  double AdlerAnglePhiPi = utils::GetAdlerAnglePhi(BeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg();

  // Angle between q and had system
  double Angleqvshad = utils::Angle(utils::GetRecoq3(out_mom, BeamE), utils::TotHadron(hadron_map).Vect()) * TMath::RadToDeg();

  double HadSystemMass = utils::HadSystemMass(hadron_map);
  double pim_E = pim_max.E();
  double pim_mom = pim_max.P();
  double pim_momx = pim_max.Px();
  double pim_momy = pim_max.Py();
  double pim_momz = pim_max.Pz();
  double pim_theta = pim_max.Theta() * TMath::RadToDeg();
  double pim_phi = pim_max.Phi() * TMath::RadToDeg();

  double MissingEnergy = utils::Missing4Momenta(BeamE, out_mom, hadron_map).E();
  double MissingMomentum = utils::Missing4Momenta(BeamE, out_mom, hadron_map).P();
  double MissingAngle = utils::Missing4Momenta(BeamE, out_mom, hadron_map).Theta() * TMath::RadToDeg();

  TLorentzVector pi_mom(0, 0, 0, 0);
  max_mom = 0;
  for (unsigned int i = 0; i < hadron_map[conf::kPdgPiP].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiP][i].P() > max_mom)
	{
	  max_mom = hadron_map[conf::kPdgPiP][i].P();
	  pi_mom = hadron_map[conf::kPdgPiP][i];
	}
    }
  max_mom = 0;
  for (unsigned int i = 0; i < hadron_map[conf::kPdgPiM].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiM][i].P() > max_mom)
	{
	  max_mom = hadron_map[conf::kPdgPiM][i].P();
	  pi_mom = hadron_map[conf::kPdgPiM][i];
	}
    }

  double RecoEvPion = utils::GetRecoEvPionProduction(out_mom, pi_mom);
  double RecoWPion = utils::GetRecoWPionProduction(out_mom, pi_mom);
  double ElectronPT = utils::GetPT(out_mom.Vect()).Mag();
  double PionPT = utils::GetPT(pi_mom.Vect()).Mag();

  int TrueNProtons = GetTrueNProtons();
  int TrueNNeutrons = GetTrueNNeutrons();
  int TrueNPiP = GetTrueNPiP();
  int TrueNPiM = GetTrueNPiM();
  int TrueNPi0 = GetTrueNPi0();
  int TrueNCh = GetTrueNCh();

  if( observable == "Efl" ) return Efl ;
  else if ( observable == "pfl" ) return pfl ;
  else if ( observable == "pflx") return pflx ;
  else if ( observable == "pfly") return pfly ;
  else if ( observable == "pflz") return pflz ;
  else if ( observable == "pfl_theta") return pfl_theta ;
  else if ( observable == "pfl_phi") return pfl_phi ;
  else if ( observable == "RecoQELEnu") return RecoQELEnu ;
  else if ( observable == "RecoEnergyTransfer") return RecoEnergyTransfer ;
  else if ( observable == "Recoq3") return Recoq3 ;
  else if ( observable == "RecoQ2") return RecoQ2 ;
  else if ( observable == "RecoW") return RecoW;
  else if ( observable == "RecoXBJK") return RecoXBJK;
  else if ( observable == "HadSystemMass") return HadSystemMass;
  else if ( observable == "ElectronSector") return ElectronSector;
  else if ( observable == "MissingEnergy") return MissingEnergy;
  else if ( observable == "MissingMomentum") return MissingMomentum;
  else if ( observable == "MissingAngle") return MissingAngle;
  else if ( observable == "ECal") return ECal;
  else if ( observable == "DiffECal") return DiffECal;
  else if ( observable == "HadronsAngle") return HadronsAngle;
  else if ( observable == "Angleqvshad") return Angleqvshad;
  else if ( observable == "RecoEvPion") return RecoEvPion;
  else if ( observable == "RecoWPion") return RecoWPion;
  else if ( observable == "ElectronPT") return ElectronPT;
  else if ( observable == "PionPT") return PionPT;
  else if ( observable == "proton_E") return proton_E;
  else if ( observable == "proton_mom") return proton_mom;
  else if ( observable == "proton_momx") return proton_momx;
  else if ( observable == "proton_momy") return proton_momy;
  else if ( observable == "proton_momz") return proton_momz;
  else if ( observable == "proton_theta") return proton_theta;
  else if ( observable == "proton_phi") return proton_phi;
  else if ( observable == "AlphaT") return AlphaT;
  else if ( observable == "DeltaPT") return DeltaPT;
  else if ( observable == "DeltaPhiT") return DeltaPhiT;
  else if ( observable == "AdlerAngleThetaP") return AdlerAngleThetaP;
  else if ( observable == "AdlerAnglePhiP") return AdlerAnglePhiP;
  else if ( observable == "pip_E") return pip_E;
  else if ( observable == "pip_mom") return pip_mom;
  else if ( observable == "pip_momx") return pip_momx;
  else if ( observable == "pip_momy") return pip_momy;
  else if ( observable == "pip_momz") return pip_momz;
  else if ( observable == "pip_theta") return pip_theta;
  else if ( observable == "pip_phi") return pip_phi;
  else if ( observable == "pim_E") return pim_E;
  else if ( observable == "pim_mom") return pim_mom;
  else if ( observable == "pim_momx") return pim_momx;
  else if ( observable == "pim_momy") return pim_momy;
  else if ( observable == "pim_momz") return pim_momz;
  else if ( observable == "pim_theta") return pim_theta;
  else if ( observable == "pim_phi") return pim_phi;
  else if ( observable == "AdlerAngleThetaPi" ) return AdlerAngleThetaPi ;
  else if ( observable == "AdlerAnglePhiPi" ) return AdlerAnglePhiPi ;
  else if ( observable == "HadAlphaT") return HadAlphaT;
  else if ( observable == "HadDeltaPT") return HadDeltaPT;
  else if ( observable == "HadDeltaPTx") return HadDeltaPTx;
  else if ( observable == "HadDeltaPTy") return HadDeltaPTy;
  else if ( observable == "HadDeltaPhiT") return HadDeltaPhiT;
  else if ( observable == "InferedNucleonMom") return InferedNucleonMom;
  else if ( observable == "TrueNProtons") return TrueNProtons;
  else if ( observable == "TrueNNeutrons") return TrueNNeutrons;
  else if ( observable == "TrueNPiP") return TrueNPiP;
  else if ( observable == "TrueNPiM") return TrueNPiM;
  else if ( observable == "TrueNPi0") return TrueNPi0;
  else if ( observable == "TrueNCh") return TrueNCh;

  std::cout << observable << " is NOT defined " << std::endl;
  return 0 ;
}

void Event::SetTargetPdg( int target_pdg ) {
  if( target_pdg == 2212 ) target_pdg = 1000010010 ;
  fTargetPdg = target_pdg ;
}

void Event::SetMottXSecWeight(void) {
  // Set Mott XSec
  double reco_Q2 = utils::GetRecoQ2( this->GetOutLepton4Mom(), this->GetInLepton4Mom().E() ) ;
  fMottXSecWght = std::pow( reco_Q2, 2 ) ;
}

TVector3 Event::GetRecoq3() const {
  return utils::GetRecoq3( fOutLepton, fInLepton.E() ) ;
}

// Radiative correction utils
void Event::SetInCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) {
  fInCorrLepton.SetPxPyPzE( px, py, pz, energy ) ;
  return ;
}
void Event::SetOutCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) {
  fOutCorrLepton.SetPxPyPzE( px, py, pz, energy ) ;
  return ;
}

//

void Event::Initialize() {
  fFinalParticles.clear() ;
  fFinalParticlesUnCorr.clear() ;
  fIsMC = false ;
  fEventID = 0 ;
  fWeight = 0 ;
  fEventID = 0 ;
  fTargetPdg = 0 ;
  fInLeptPdg = 11 ;
  fOutLeptPdg = 11 ;
  fNP = 0 ;
  fNN = 0 ;
  fNPiP = 0 ;
  fNPiM = 0 ;
  fNPi0 = 0 ;
  fNKM = 0 ;
  fNKP = 0 ;
  fNK0 = 0 ;
  fNEM = 0 ;
  fNOther = 0 ;

  fIsMC = true ;
  fIsEM = false ;
  fIsCC = false ;
  fIsNC = false ;
  fIsQEL = false ;
  fIsRES = false ;
  fIsMEC = false ;
  fIsDIS = false ;

  fTrueQ2s = 0 ;
  fTrueWs = 0 ;
  fTruexs = 0 ;
  fTrueys = 0 ;
  fTrueQ2 = 0 ;
  fTrueW = 0 ;
  fTruex = 0 ;
  fTruey = 0 ;
  fresid = -1000;

  fInLepton.SetPxPyPzE( 0,0,0,0 ) ;
  fOutLepton.SetPxPyPzE( 0,0,0,0 ) ;
  fInCorrLepton.SetPxPyPzE( 0,0,0,0 ) ;
  fOutCorrLepton.SetPxPyPzE( 0,0,0,0 ) ;
  fInLeptonUnCorr.SetPxPyPzE( 0,0,0,0 ) ;
  fOutLeptonUnCorr.SetPxPyPzE( 0,0,0,0 ) ;

  fAnalysisRecord.clear();
}

void Event::Clear() {

  fFinalParticles.clear() ;
  fFinalParticlesUnCorr.clear() ;
  fAnalysisRecord.clear();
}
