/*
 * MC Analysis analysis
 * This class deals with analysis features specific for MC data
 *
 */
#include <iostream>
#include "TFile.h"
#include "TDirectoryFile.h"
#include "analysis/MCCLAS6AnalysisI.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"
#include "utils/DetectorUtils.h"
#include "conf/AnalysisCutsI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/ConstantsI.h"

using namespace e4nu;

MCCLAS6AnalysisI::MCCLAS6AnalysisI()
{
  kMult_signal = GetNTopologyParticles();
  this->Initialize();
}

MCCLAS6AnalysisI::~MCCLAS6AnalysisI()
{
  delete fData;
}

bool MCCLAS6AnalysisI::LoadData(void)
{
  if (!IsConfigured())
    return false;

  std::string file = GetInputFile();
  double nevents = GetNEventsToRun();
  double first_event = GetFirstEventToRun();

  if (!kIsDataLoaded)
  {
    fData = new MCEventHolder(file, first_event, nevents);
    kNEvents = fData->GetNEvents();
    kIsDataLoaded = true;
  }
  return kIsDataLoaded;
}

Event *MCCLAS6AnalysisI::GetEvent(const unsigned int event_id)
{
  return fData->GetEvent(event_id);
}

Event *MCCLAS6AnalysisI::GetValidEvent(const unsigned int event_id)
{

  Event *event;
  if (IsNoFSI())
  {
    // This function will load the full event using pre FSI nucleon kinematics
    event = fData->GetEventNoFSI(event_id);
  }
  else
    event = fData->GetEvent(event_id);

  if (!event)
  {
    delete event;
    return nullptr;
  }

  // Apply Generic analysis cuts (0-3)
  if (!AnalysisI::Analyse(*event))
  {
    delete event;
    return nullptr;
  }

  // Step 3 : smear particles momentum
  if (ApplyReso())
  {
    this->SmearParticles(*event);
  }

  // Step 4 : Remove true Bkg events if requested :
  if (IsTrueSignal())
  {
    // Apply theta cut on hadrons:
    if (!this->ApplyFiducialCut(*event, false))
    {
      delete event;
      return nullptr;
    }
    std::map<int, unsigned int> Topology = GetTopology();
    std::map<int, std::vector<TLorentzVector>> hadrons = event->GetFinalParticles4Mom();
    for (auto it = Topology.begin(); it != Topology.end(); ++it)
    {
      if (it->first == conf::kPdgElectron)
        continue;
      if (hadrons[it->first].size() != it->second)
      {
        delete event;
        return nullptr;
      }
    }
  }

  // Step 5 : Apply fiducial
  // Moved to general Analysis as it is till used for the data when we apply a systematic shift to Phi
  // to compute geometric acceptance systematic.

  // Step 5: Apply Acceptance Correction (Efficiency correction)
  // This takes into account the efficiency detection of each particle in theta and phi
  // We want to apply it at the end to correctly account for the acceptance in background substracted events
  this->ApplyAcceptanceCorrection(*event);

  // Store analysis record after fiducial cut and acceptance correction (2):
  event->StoreAnalysisRecord(kid_fid);
  return event;
}

void MCCLAS6AnalysisI::SmearParticles(Event &event)
{
  double EBeam = GetConfiguredEBeam();
  TLorentzVector out_mom = event.GetOutLepton4Mom();

  utils::ApplyResolution(conf::kPdgElectron, out_mom, EBeam);
  event.SetOutLeptonKinematics(out_mom);

  // Apply for other particles
  std::map<int, std::vector<TLorentzVector>> part_map = event.GetFinalParticles4Mom();
  for (std::map<int, std::vector<TLorentzVector>>::iterator it = part_map.begin(); it != part_map.end(); ++it)
  {
    std::vector<TLorentzVector> vtemp;
    for (unsigned int i = 0; i < (it->second).size(); ++i)
    {
      TLorentzVector temp = (it->second)[i];
      utils::ApplyResolution(it->first, temp, EBeam);
      vtemp.push_back(temp);
    }
    part_map[it->first] = vtemp;
  }
  event.SetFinalParticlesKinematics(part_map);
}

unsigned int MCCLAS6AnalysisI::GetNEvents(void) const
{
  return (unsigned int)fData->GetNEvents();
}

void MCCLAS6AnalysisI::Initialize()
{
  if (IsData())
    return;
  fData = nullptr;

  // Get run configurables
  double EBeam = GetConfiguredEBeam();
  unsigned int Target = GetConfiguredTarget();

  if (NormalizeHist())
  {
    // Get xsec from cross section spline
    std::string target_tag = "";
    if (GetConfiguredTarget() == conf::kPdgC12)
      target_tag = "e-_C12";
    if (GetConfiguredTarget() == conf::kPdgHe3)
      target_tag = "e-_He3";
    if (GetConfiguredTarget() == conf::kPdgHe4)
      target_tag = "e-_He4";
    if (GetConfiguredTarget() == conf::kPdgFe56)
      target_tag = "e-_Fe56";
    if (GetConfiguredTarget() == conf::kPdgO16)
      target_tag = "e-_O16";

    std::cout << " *********************************** " << GetXSecFile()<< std::endl;
    std::unique_ptr<TFile> xsec_file = std::unique_ptr<TFile>(new TFile((GetXSecFile()).c_str(), "READ"));
    if (!xsec_file)
    {
      std::cout << " ERROR: Xsec file does not exist: " << GetXSecFile() << std::endl;
      kIsConfigured = false;
      return;
    }

    TDirectoryFile *xsec_dir = (TDirectoryFile *)xsec_file->Get(target_tag.c_str());
    if (!xsec_dir)
    {
      std::cout << " ERROR: Xsec dir does not exist: " << target_tag << std::endl;
      kIsConfigured = false;
      return;
    }

    TGraph *gxsec = (TGraph *)xsec_dir->Get("tot_em");
    if (!gxsec)
    {
      std::cout << " ERROR: Cannot create graph for " << "tot_em" << std::endl;
      kIsConfigured = false;
      return;
    }

    kXSec = gxsec->Eval(GetConfiguredEBeam());
    std::cout << " Total XSec (" << GetConfiguredEBeam() << ") = " << kXSec << std::endl;
    xsec_file->Close();
  }
}

bool MCCLAS6AnalysisI::Finalise(std::map<int, std::vector<e4nu::Event>> &event_holder)
{
  if (!AnalysisI::Finalise())
    return false;

  auto tags = GetObservablesTag() ;
  unsigned int min_mult = GetMinBkgMult();
  for (unsigned int k = 0; k < event_holder[min_mult].size(); ++k)
  {
    // Store corrected background in event sample
    StoreTree(event_holder[min_mult][k]);
    double norm_weight = event_holder[min_mult][k].GetTotalWeight();

    // Store in histogram(s)
    for (unsigned int obs_id = 0; obs_id < tags.size(); ++obs_id)
    {
      if( !kHistograms[tags[obs_id]][0] ) continue ; 
      kHistograms[tags[obs_id]][0]->Fill(event_holder[min_mult][k].GetObservable(tags[obs_id]), norm_weight);
    }

    PlotBkgInformation(event_holder[min_mult][k]);
  }

  // Normalize
  if (NormalizeHist())
  {
    for (unsigned int obs_id = 0; obs_id < tags.size(); ++obs_id)
    {
      if( !kHistograms[tags[obs_id]][0] ) continue ;
      double NBins = kHistograms[tags[obs_id]][0]->GetNbinsX();

      for (int k = 1; k <= NBins; k++)
      {
        double content = kHistograms[tags[obs_id]][0]->GetBinContent(k);
        double error = kHistograms[tags[obs_id]][0]->GetBinError(k);
        double width = kHistograms[tags[obs_id]][0]->GetBinWidth(k);
        double newcontent = content / width;
        double newerror = error / width;
        kHistograms[tags[obs_id]][0]->SetBinContent(k, newcontent);
        kHistograms[tags[obs_id]][0]->SetBinError(k, newerror);
      }

      kHistograms[tags[obs_id]][0]->Scale(kXSec * kConversionFactorCm2ToMicroBarn * TMath::Power(10., -38.) / GetNEventsToRun());
    }
  }

  return true;
}

bool MCCLAS6AnalysisI::StoreTree(Event event)
{
  static bool n = true;
  int ID = event.GetEventID();
  int TargetPdg = event.GetTargetPdg();
  int InLeptonPdg = event.GetInLeptPdg();
  int OutLeptonPdg = event.GetOutLeptPdg();
  double TotWeight = event.GetTotalWeight();
  double AccWght = event.GetAccWght();
  double EventWght = event.GetEventWeight();
  double BeamE = event.GetInLepton4Mom().E();

  bool CC = event.IsCC();
  bool NC = event.IsNC();
  bool EM = event.IsEM();
  bool QEL = event.IsQEL();
  bool RES = event.IsRES();
  bool MEC = event.IsMEC();
  bool DIS = event.IsDIS();

  unsigned int TrueNProtons = event.GetTrueNProtons();
  unsigned int TrueNNeutrons = event.GetTrueNNeutrons();
  unsigned int TrueNPiP = event.GetTrueNPiP();
  unsigned int TrueNPiM = event.GetTrueNPiM();
  unsigned int TrueNPi0 = event.GetTrueNPi0();
  unsigned int TrueNKP = event.GetTrueNKP();
  unsigned int TrueNKM = event.GetTrueNKM();
  unsigned int TrueNK0 = event.GetTrueNK0();
  unsigned int TrueNEM = event.GetTrueNEM();
  unsigned int TrueNOther = event.GetTrueNOther();

  unsigned int RecoNProtons = event.GetRecoNProtons();
  unsigned int RecoNNeutrons = event.GetRecoNNeutrons();
  unsigned int RecoNPiP = event.GetRecoNPiP();
  unsigned int RecoNPiM = event.GetRecoNPiM();
  unsigned int RecoNPi0 = event.GetRecoNPi0();
  unsigned int RecoNKP = event.GetRecoNKP();
  unsigned int RecoNKM = event.GetRecoNKM();
  unsigned int RecoNK0 = event.GetRecoNK0();
  unsigned int RecoNEM = event.GetRecoNEM();

  double TrueQ2s = event.GetTrueQ2s();
  double TrueWs = event.GetTrueWs();
  double Truexs = event.GetTruexs();
  double Trueys = event.GetTrueys();
  double TrueQ2 = event.GetTrueQ2();
  double TrueW = event.GetTrueW();
  double Truex = event.GetTruex();
  double Truey = event.GetTruey();
  int resid = event.GetRESID();

  TLorentzVector out_mom = event.GetOutLepton4Mom();
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

  double MottXSecScale = event.GetMottXSecWeight();

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

  // flip hadrons phi
  std::map<int, std::vector<TLorentzVector>> hadron_map = event.GetFinalParticles4Mom();
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

  double RecoEvPion = utils::GetRecoEvPionProduction(out_mom, pi_mom);
  double RecoWPion = utils::GetRecoWPionProduction(out_mom, pi_mom);
  double ElectronPT = utils::GetPT(out_mom.Vect()).Mag();
  double PionPT = utils::GetPT(pi_mom.Vect()).Mag();

  bool IsBkg = event.IsBkg();
  unsigned int InitialNEvents = GetNEventsToRun();
  double ConversionFactor = kConversionFactorCm2ToMicroBarn * TMath::Power(10., -38.);
  double TotalXSec = kXSec;
  double MCNormalization = TotalXSec * ConversionFactor / InitialNEvents;

  // Store name of files used
  const char *InputROOTFile = kInputFile.c_str();
  const char *OutputROOTFile = kOutputFile.c_str();
  const char *InputXSecFile = kXSecFile.c_str();

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
    kAnalysisTree->Branch("CC", &CC, "CC/O");
    kAnalysisTree->Branch("NC", &NC, "NC/O");
    kAnalysisTree->Branch("EM", &EM, "EC/O");
    kAnalysisTree->Branch("QEL", &QEL, "QEL/O");
    kAnalysisTree->Branch("RES", &RES, "RES/O");
    kAnalysisTree->Branch("MEC", &MEC, "MEC/O");
    kAnalysisTree->Branch("DIS", &DIS, "DIS/O");
    kAnalysisTree->Branch("TrueQ2s", &TrueQ2s, "TrueQ2s/D");
    kAnalysisTree->Branch("TrueWs", &TrueWs, "TrueWs/D");
    kAnalysisTree->Branch("Truexs", &Truexs, "Truexs/D");
    kAnalysisTree->Branch("Trueys", &Trueys, "Trueys/D");
    kAnalysisTree->Branch("TrueQ2", &TrueQ2, "TrueQ2/D");
    kAnalysisTree->Branch("TrueW", &TrueW, "TrueW/D");
    kAnalysisTree->Branch("Truex", &Truex, "Truex/D");
    kAnalysisTree->Branch("Truey", &Truey, "Truey/D");
    kAnalysisTree->Branch("TotWeight", &TotWeight, "TotWeight/D");
    kAnalysisTree->Branch("EventWght", &EventWght, "EventWght/D");
    kAnalysisTree->Branch("AccWght", &AccWght, "AccWght/D");
    kAnalysisTree->Branch("MottXSecScale", &MottXSecScale, "MottXSecScale/D");
    kAnalysisTree->Branch("IsBkg", &IsBkg, "IsBkg/O");
    kAnalysisTree->Branch("TrueNProtons", &TrueNProtons, "TrueNProtons/I");
    kAnalysisTree->Branch("TrueNNeutrons", &TrueNNeutrons, "TrueNNeutrons/I");
    kAnalysisTree->Branch("TrueNPiP", &TrueNPiP, "TrueNPiP/I");
    kAnalysisTree->Branch("TrueNPiM", &TrueNPiM, "TrueNPiM/I");
    kAnalysisTree->Branch("TrueNPi0", &TrueNPi0, "TrueNPi0/I");
    kAnalysisTree->Branch("TrueNKP", &TrueNKP, "TrueNKP/I");
    kAnalysisTree->Branch("TrueNKM", &TrueNKM, "TrueNKM/I");
    kAnalysisTree->Branch("TrueNK0", &TrueNK0, "TrueNK0/I");
    kAnalysisTree->Branch("TrueNEM", &TrueNEM, "TrueNEM/I");
    kAnalysisTree->Branch("TrueNOther", &TrueNOther, "TrueNOther/I");
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
    kAnalysisTree->Branch("HadSystemMass", &HadSystemMass, "HadSystemMass/D");
    kAnalysisTree->Branch("ElectronSector", &ElectronSector, "ElectronSector/I");
    kAnalysisTree->Branch("MissingEnergy", &MissingEnergy, "MissingEnergy/D");
    kAnalysisTree->Branch("MissingMomentum", &MissingMomentum, "MissingMomentum/D");
    kAnalysisTree->Branch("MissingAngle", &MissingAngle, "MissingAngle/D");
    kAnalysisTree->Branch("ECal", &ECal, "ECal/D");
    kAnalysisTree->Branch("DiffECal", &DiffECal, "DiffECal/D");
    kAnalysisTree->Branch("resid", &resid, "resid/I");
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
    kAnalysisTree->Branch("InferedNucleonMom", &InferedNucleonMom, "InferedNucleonMom/D");

    // Add normaization information
    kAnalysisTree->Branch("InitialNEvents", &InitialNEvents, "InitialNEvents/I");
    kAnalysisTree->Branch("ConversionFactor", &ConversionFactor, "ConversionFactor/D");
    kAnalysisTree->Branch("TotalXSec", &TotalXSec, "TotalXSec/D");
    kAnalysisTree->Branch("MCNormalization", &MCNormalization, "MCNormalization/D");

    n = false;
  }

  kAnalysisTree->Fill();

  return true;
}
