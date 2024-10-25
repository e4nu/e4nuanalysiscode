/*
 * Analysis Interface base class
 *
 */
#include <iostream>
#include "TFile.h"
#include "TDirectoryFile.h"
#include "analysis/CLAS6AnalysisI.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"
#include "utils/ParticleUtils.h"
#include "utils/DetectorUtils.h"
#include "conf/FiducialCutI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/CLAS6ConstantsI.h"
#include "conf/ConstantsI.h"

using namespace e4nu;

CLAS6AnalysisI::CLAS6AnalysisI()
{
  kMult_signal = GetNTopologyParticles();
  this->Initialize();
}

CLAS6AnalysisI::~CLAS6AnalysisI()
{
  delete fData;
}

bool CLAS6AnalysisI::LoadData(void)
{
  if (!IsConfigured())
    return false;

  std::string file = GetInputFile();
  double nevents = GetNEventsToRun();
  double first_event = GetFirstEventToRun();

  if (!kIsDataLoaded)
  {
    fData = new CLAS6EventHolder(file, first_event, nevents);
    kNEvents = fData->GetNEvents();
    kIsDataLoaded = true;
  }
  return kIsDataLoaded;
}

Event *CLAS6AnalysisI::GetEvent(const unsigned int event_id)
{
  return fData->GetEvent(event_id);
}

Event *CLAS6AnalysisI::GetValidEvent(const unsigned int event_id)
{

  Event *event = fData->GetEvent(event_id);
  if (!event)
  {
    delete event;
    return nullptr;
  }

  // Apply Generic analysis cuts
  if (!AnalysisI::Analyse(*event))
  {
    delete event;
    return nullptr;
  }

  // No further code is needed
  // Fiducial cuts are already taken care of
  // No Need to apply them again

  return event;
}

unsigned int CLAS6AnalysisI::GetNEvents(void) const
{
  return (unsigned int)fData->GetNEvents();
}

void CLAS6AnalysisI::Initialize()
{
  fData = nullptr;
}

bool CLAS6AnalysisI::Finalise(std::map<int, std::vector<e4nu::Event>> &event_holder)
{

  if (!AnalysisI::Finalise())
    return false;

  // Store corrected background in event sample
  unsigned int min_mult = GetMinBkgMult();
  for (unsigned int k = 0; k < event_holder[min_mult].size(); ++k)
  {
    StoreTree(event_holder[min_mult][k]);

    double norm_weight = 1;
    if (ApplyCorrWeights())
    {
      norm_weight = event_holder[min_mult][k].GetTotalWeight();
    }

    // Store in histogram(s)
    for (unsigned int j = 0; j < kHistograms.size(); ++j)
    {
      kHistograms[j]->Fill(event_holder[min_mult][k].GetObservable(GetObservablesTag()[j]), norm_weight);
    }
  }

  // Normalize
  unsigned int tgt_pdg = GetConfiguredTarget();
  double EBeam = GetConfiguredEBeam();

  // Get constants
  unsigned int MassNumber = utils::GetMassNumber(tgt_pdg);
  double IntegratedCharge = conf::GetIntegratedCharge(tgt_pdg, EBeam);
  double TargetLength = conf::GetTargetLength(tgt_pdg);
  double TargetDensity = conf::GetTargetDensity(tgt_pdg);

  if (NormalizeHist())
  {
    for (unsigned int j = 0; j < kHistograms.size(); ++j)
    {
      double NBins = kHistograms[j]->GetNbinsX();

      for (int k = 1; k <= NBins; k++)
      {
        double content = kHistograms[j]->GetBinContent(k);
        double error = kHistograms[j]->GetBinError(k);
        double width = kHistograms[j]->GetBinWidth(k);
        double newcontent = content / width;
        double newerror = error / width;
        kHistograms[j]->SetBinContent(k, newcontent);
        kHistograms[j]->SetBinError(k, newerror);
      }
      kHistograms[j]->Scale(kConversionFactorCm2ToMicroBarn * MassNumber / (IntegratedCharge * TargetLength * TargetDensity * kOverallUnitConversionFactor));
    }
  }

  return true;
}

bool CLAS6AnalysisI::StoreTree(Event event)
{
  static bool n = true;
  int ID = event.GetEventID();
  int TargetPdg = event.GetTargetPdg();
  int InLeptonPdg = event.GetInLeptPdg();
  int OutLeptonPdg = event.GetOutLeptPdg();
  double TotWeight = event.GetTotalWeight();
  double EventWght = event.GetEventWeight();
  double BeamE = event.GetInLepton4Mom().E();

  unsigned int RecoNProtons = event.GetRecoNProtons();
  unsigned int RecoNNeutrons = event.GetRecoNNeutrons();
  unsigned int RecoNPiP = event.GetRecoNPiP();
  unsigned int RecoNPiM = event.GetRecoNPiM();
  unsigned int RecoNPi0 = event.GetRecoNPi0();
  unsigned int RecoNKP = event.GetRecoNKP();
  unsigned int RecoNKM = event.GetRecoNKM();
  unsigned int RecoNK0 = event.GetRecoNK0();
  unsigned int RecoNEM = event.GetRecoNEM();

  TLorentzVector out_mom = event.GetOutLepton4Mom();
  double Efl = out_mom.E();
  double pfl = out_mom.P();
  double pflx = out_mom.Px();
  double pfly = out_mom.Py();
  double pflz = out_mom.Pz();
  double pfl_theta = out_mom.Theta() * TMath::RadToDeg();
  double pfl_phi = out_mom.Phi();
  unsigned int ElectronSector = utils::GetSector(pfl_phi);
  pfl_phi *= TMath::RadToDeg();

  double RecoQELEnu = utils::GetQELRecoEnu(out_mom, TargetPdg);
  double RecoEnergyTransfer = utils::GetEnergyTransfer(out_mom, BeamE);
  double Recoq3 = utils::GetRecoq3(out_mom, BeamE).Mag();
  double RecoQ2 = utils::GetRecoQ2(out_mom, BeamE);
  double RecoXBJK = utils::GetRecoXBJK(out_mom, BeamE);
  double RecoW = utils::GetRecoW(out_mom, BeamE);

  std::map<int, unsigned int> topology = GetTopology();
  static bool topology_has_protons = false;
  if (topology.count(conf::kPdgProton) != 0)
  {
    if (topology[conf::kPdgProton] != 0)
      topology_has_protons = true;
  }
  static bool topology_has_pip = false;
  if (topology.count(conf::kPdgPiP) != 0)
  {
    if (topology[conf::kPdgPiP] != 0)
      topology_has_pip = true;
  }

  static bool topology_has_pim = false;
  if (topology.count(conf::kPdgPiM) != 0)
  {
    if (topology[conf::kPdgPiM] != 0)
      topology_has_pim = true;
  }

  unsigned int TopMult = GetNTopologyParticles();
  std::map<int, std::vector<TLorentzVector>> hadron_map = event.GetFinalParticles4Mom();
  TLorentzVector p_max(0, 0, 0, 0);
  if (topology_has_protons)
  {
    double max_mom = 0;
    for (unsigned int i = 0; i < hadron_map[conf::kPdgProton].size(); ++i)
    {
      if (hadron_map[conf::kPdgProton][i].P() > max_mom)
      {
        max_mom = hadron_map[conf::kPdgProton][i].P();
        p_max = hadron_map[conf::kPdgProton][i];
      }
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
  double ECal = utils::GetECal(out_mom.E(), event.GetFinalParticles4Mom(), TargetPdg);
  double DiffECal = utils::GetECal(out_mom.E(), event.GetFinalParticles4Mom(), TargetPdg) - BeamE;
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
  if (topology_has_pip)
  {
    double max_mom = 0;
    for (unsigned int i = 0; i < hadron_map[conf::kPdgPiP].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiP][i].P() > max_mom)
      {
        max_mom = hadron_map[conf::kPdgPiP][i].P();
        pip_max = hadron_map[conf::kPdgPiP][i];
      }
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
  if (topology_has_pim)
  {
    double max_mom = 0;
    for (unsigned int i = 0; i < hadron_map[conf::kPdgPiM].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiM][i].P() > max_mom)
      {
        max_mom = hadron_map[conf::kPdgPiM][i].P();
        pim_max = hadron_map[conf::kPdgPiM][i];
      }
    }
  }

  double pim_E = pim_max.E();
  double pim_mom = pim_max.P();
  double pim_momx = pim_max.Px();
  double pim_momy = pim_max.Py();
  double pim_momz = pim_max.Pz();
  double pim_theta = pim_max.Theta() * TMath::RadToDeg();
  double pim_phi = pim_max.Phi() * TMath::RadToDeg();
  double HadSystemMass = utils::HadSystemMass(hadron_map);

  double MissingEnergy = utils::Missing4Momenta(BeamE, out_mom, hadron_map).E();
  double MissingMomentum = utils::Missing4Momenta(BeamE, out_mom, hadron_map).P();
  double MissingAngle = utils::Missing4Momenta(BeamE, out_mom, hadron_map).Theta() * TMath::RadToDeg();

  unsigned int MassNumber = utils::GetMassNumber(TargetPdg);
  double IntegratedCharge = conf::GetIntegratedCharge(TargetPdg, BeamE);
  double TargetLength = conf::GetTargetLength(TargetPdg);
  double TargetDensity = conf::GetTargetDensity(TargetPdg);
  unsigned int InitialNEvents = kNEvents;
  double ConversionFactor = kConversionFactorCm2ToMicroBarn / kOverallUnitConversionFactor;
  double DataNormalization = kConversionFactorCm2ToMicroBarn * MassNumber / (IntegratedCharge * TargetLength * TargetDensity * kOverallUnitConversionFactor);

  // Adler angles
  double AdlerAngleThetaP = utils::GetAdlerAngleTheta(BeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  double AdlerAnglePhiP = utils::GetAdlerAnglePhi(BeamE, out_mom, hadron_map, conf::kPdgProton) * TMath::RadToDeg();
  double AdlerAngleThetaPi = utils::GetAdlerAngleTheta(BeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg();
  double AdlerAnglePhiPi = utils::GetAdlerAnglePhi(BeamE, out_mom, hadron_map, abs(conf::kPdgPiM)) * TMath::RadToDeg();

  double Angleqvshad = utils::Angle(utils::GetRecoq3(out_mom, BeamE), utils::TotHadron(hadron_map).Vect()) * TMath::RadToDeg();

  TLorentzVector pi_mom(0, 0, 0, 0);
  if (topology_has_pip)
  {
    double max_mom = 0;
    for (unsigned int i = 0; i < hadron_map[conf::kPdgPiP].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiP][i].P() > max_mom)
      {
        max_mom = hadron_map[conf::kPdgPiP][i].P();
        pi_mom = hadron_map[conf::kPdgPiP][i];
      }
    }
  }
  else if (topology_has_pim)
  {
    double max_mom = 0;
    for (unsigned int i = 0; i < hadron_map[conf::kPdgPiM].size(); ++i)
    {
      if (hadron_map[conf::kPdgPiM][i].P() > max_mom)
      {
        max_mom = hadron_map[conf::kPdgPiM][i].P();
        pi_mom = hadron_map[conf::kPdgPiM][i];
      }
    }
  }

  double RecoEvPion = utils::GetRecoEvPionProduction(BeamE, out_mom, pi_mom);
  double RecoWPion = utils::GetRecoWPionProduction(BeamE, out_mom, pi_mom);
  double ElectronPT = utils::GetPT(out_mom.Vect()).Mag();
  double PionPT = utils::GetPT(pi_mom.Vect()).Mag();

  // Store name of files used
  const char *InputROOTFile = kInputFile.c_str();
  const char *OutputROOTFile = kOutputFile.c_str();
  const char *InputXSecFile = kXSecFile.c_str();

  bool IsBkg = event.IsBkg();
  if (n == true)
  {
    kAnalysisTree->Branch("InputROOTFile", &InputROOTFile, "InputROOTFile/C");
    kAnalysisTree->Branch("OutputROOTFile", &OutputROOTFile, "OutputROOTFile/C");
    kAnalysisTree->Branch("InputXSecFile", &InputXSecFile, "InputXSecFile/C");
    kAnalysisTree->Branch("ID", &ID, "ID/I");
    kAnalysisTree->Branch("TargetPdg", &TargetPdg, "TargetPdg/I");
    kAnalysisTree->Branch("InLeptonPdg", &InLeptonPdg, "InLeptonPdg/I");
    kAnalysisTree->Branch("OutLeptonPdg", &OutLeptonPdg, "OutLeptonPdg/I");
    kAnalysisTree->Branch("BeamE", &BeamE, "BeamE/D");
    kAnalysisTree->Branch("TotWeight", &TotWeight, "TotWeight/D");
    kAnalysisTree->Branch("EventWght", &EventWght, "EventWght/D");
    kAnalysisTree->Branch("IsBkg", &IsBkg, "IsBkg/O");
    kAnalysisTree->Branch("RecoNProtons", &RecoNProtons, "RecoNProtons/I");
    kAnalysisTree->Branch("RecoNNeutrons", &RecoNNeutrons, "RecoNNeutrons/I");
    kAnalysisTree->Branch("RecoNPiP", &RecoNPiP, "RecoNPiP/I");
    kAnalysisTree->Branch("RecoNPiM", &RecoNPiM, "RecoNPiM/I");
    kAnalysisTree->Branch("RecoNPi0", &RecoNPi0, "RecoNPi0/I");
    kAnalysisTree->Branch("RecoNKP", &RecoNKP, "RecoNKP/I");
    kAnalysisTree->Branch("RecoNKM", &RecoNKM, "RecoNKM/I");
    kAnalysisTree->Branch("RecoNK0", &RecoNK0, "RecoNK0/I");
    kAnalysisTree->Branch("RecoNEM", &RecoNEM, "RecoNEM/I");
    kAnalysisTree->Branch("TopMult", &TopMult, "TopMult/I");
    kAnalysisTree->Branch("Efl", &Efl, "Efl/D");
    kAnalysisTree->Branch("pfl", &pfl, "pfl/D");
    kAnalysisTree->Branch("pflx", &pflx, "pflx/D");
    kAnalysisTree->Branch("pfly", &pfly, "pfly/D");
    kAnalysisTree->Branch("pflz", &pflz, "pflz/D");
    kAnalysisTree->Branch("pfl_theta", &pfl_theta, "pfl_theta/D");
    kAnalysisTree->Branch("pfl_phi", &pfl_phi, "pfl_phi/D");
    kAnalysisTree->Branch("RecoQELEnu", &RecoQELEnu, "RecoQELEnu/D");
    kAnalysisTree->Branch("RecoEnergyTransfer", &RecoEnergyTransfer, "RecoEnergyTransfer/D");
    kAnalysisTree->Branch("Recoq3", &Recoq3, "Recoq3/D");
    kAnalysisTree->Branch("RecoQ2", &RecoQ2, "RecoQ2/D");
    kAnalysisTree->Branch("RecoW", &RecoW, "RecoW/D");
    kAnalysisTree->Branch("RecoXBJK", &RecoXBJK, "RecoXBJK/D");
    kAnalysisTree->Branch("ElectronSector", &ElectronSector, "ElectronSector/I");
    kAnalysisTree->Branch("MissingEnergy", &MissingEnergy, "MissingEnergy/D");
    kAnalysisTree->Branch("MissingMomentum", &MissingMomentum, "MissingMomentum/D");
    kAnalysisTree->Branch("MissingAngle", &MissingAngle, "MissingAngle/D");
    kAnalysisTree->Branch("ECal", &ECal, "ECal/D");
    kAnalysisTree->Branch("DiffECal", &DiffECal, "DiffECal/D");
    kAnalysisTree->Branch("HadronsAngle", &HadronsAngle, "HadronsAngle/D");
    kAnalysisTree->Branch("Angleqvshad", &Angleqvshad, "Angleqvshad/D");
    kAnalysisTree->Branch("RecoEvPion", &RecoEvPion, "RecoEvPion/D");
    kAnalysisTree->Branch("RecoWPion", &RecoWPion, "RecoWPion/D");
    kAnalysisTree->Branch("ElectronPT", &ElectronPT, "ElectronPT/D");
    kAnalysisTree->Branch("PionPT", &PionPT, "PionPT/D");

    if (topology_has_protons)
    {
      kAnalysisTree->Branch("proton_E", &proton_E, "proton_E/D");
      kAnalysisTree->Branch("proton_mom", &proton_mom, "proton_mom/D");
      kAnalysisTree->Branch("proton_momx", &proton_momx, "proton_momx/D");
      kAnalysisTree->Branch("proton_momy", &proton_momy, "proton_momy/D");
      kAnalysisTree->Branch("proton_momz", &proton_momz, "proton_momz/D");
      kAnalysisTree->Branch("proton_theta", &proton_theta, "proton_theta/D");
      kAnalysisTree->Branch("proton_phi", &proton_phi, "proton_phi/D");
      kAnalysisTree->Branch("AlphaT", &AlphaT, "AlphaT/D");
      kAnalysisTree->Branch("DeltaPT", &DeltaPT, "DeltaPT/D");
      kAnalysisTree->Branch("DeltaPhiT", &DeltaPhiT, "DeltaPhiT/D");
      kAnalysisTree->Branch("AdlerAngleThetaP", &AdlerAngleThetaP, "AdlerAngleThetaP/D");
      kAnalysisTree->Branch("AdlerAnglePhiP", &AdlerAnglePhiP, "AdlerAnglePhiP/D");
    }

    if (topology_has_pip)
    {
      kAnalysisTree->Branch("pip_E", &pip_E, "pip_E/D");
      kAnalysisTree->Branch("pip_mom", &pip_mom, "pip_mom/D");
      kAnalysisTree->Branch("pip_momx", &pip_momx, "pip_momx/D");
      kAnalysisTree->Branch("pip_momy", &pip_momy, "pip_momy/D");
      kAnalysisTree->Branch("pip_momz", &pip_momz, "pip_momz/D");
      kAnalysisTree->Branch("pip_theta", &pip_theta, "pip_theta/D");
      kAnalysisTree->Branch("pip_phi", &pip_phi, "pip_phi/D");
    }

    if (topology_has_pim)
    {
      kAnalysisTree->Branch("pim_E", &pim_E, "pim_E/D");
      kAnalysisTree->Branch("pim_mom", &pim_mom, "pim_mom/D");
      kAnalysisTree->Branch("pim_momx", &pim_momx, "pim_momx/D");
      kAnalysisTree->Branch("pim_momy", &pim_momy, "pim_momy/D");
      kAnalysisTree->Branch("pim_momz", &pim_momz, "pim_momz/D");
      kAnalysisTree->Branch("pim_theta", &pim_theta, "pim_theta/D");
      kAnalysisTree->Branch("pim_phi", &pim_phi, "pim_phi/D");
    }

    if (topology_has_pip || topology_has_pim)
    {
      kAnalysisTree->Branch("AdlerAngleThetaPi", &AdlerAngleThetaPi, "AdlerAngleThetaPi/D");
      kAnalysisTree->Branch("AdlerAnglePhiPi", &AdlerAnglePhiPi, "AdlerAnglePhiPi/D");
    }

    kAnalysisTree->Branch("HadAlphaT", &HadAlphaT, "HadAlphaT/D");
    kAnalysisTree->Branch("HadDeltaPT", &HadDeltaPT, "HadDeltaPT/D");
    kAnalysisTree->Branch("HadDeltaPTx", &HadDeltaPTx, "HadDeltaPTx/D");
    kAnalysisTree->Branch("HadDeltaPTy", &HadDeltaPTy, "HadDeltaPTy/D");
    kAnalysisTree->Branch("HadDeltaPhiT", &HadDeltaPhiT, "HadDeltaPhiT/D");
    kAnalysisTree->Branch("HadSystemMass", &HadSystemMass, "HadSystemMass/D");
    kAnalysisTree->Branch("InferedNucleonMom", &InferedNucleonMom, "InferedNucleonMom/D");

    // Adding Normalization information
    kAnalysisTree->Branch("MassNumber", &MassNumber, "MassNumber/I");
    kAnalysisTree->Branch("IntegratedCharge", &IntegratedCharge, "IntegratedCharge/D");
    kAnalysisTree->Branch("TargetLength", &TargetLength, "TargetLength/D");
    kAnalysisTree->Branch("TargetDensity", &TargetDensity, "TargetDensity/D");
    kAnalysisTree->Branch("InitialNEvents", &InitialNEvents, "InitialNEvents/I");
    kAnalysisTree->Branch("ConversionFactor", &ConversionFactor, "ConversionFactor/D");
    kAnalysisTree->Branch("DataNormalization", &DataNormalization, "DataNormalization/D");

    n = false;
  }

  kAnalysisTree->Fill();
  return true;
}
