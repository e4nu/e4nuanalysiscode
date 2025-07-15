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

  // (*) Step 4 : Smear before appling mom cuts 
  if( ApplyReso() && !IsData() ) {
    this -> SmearParticles( *event ) ; 
  }
  
  // Step 5: Apply Acceptance Correction (Efficiency correction)
  // This takes into account the efficiency detection of each particle in theta and phi
  // We want to apply it at the end to correctly account for the acceptance in background substracted events
  this->ApplyAcceptanceCorrection(*event);

  // Store analysis record after fiducial cut and acceptance correction (2):
  event->StoreAnalysisRecord(kid_fid);

  return event;
}

unsigned int MCCLAS6AnalysisI::GetNEvents(void) const
{
  return (unsigned int)fData->GetNEvents();
}

void MCCLAS6AnalysisI::Initialize()
{
  if (IsData()) return;
  fData = nullptr;
  
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

  if (NormalizeHist()) { 
    for (unsigned int obs_id = 0; obs_id < tags.size(); ++obs_id) { 
      for( unsigned int i = 0 ; i < kHistograms[tags[obs_id]].size() ; ++i ) {
	if( !kHistograms[tags[obs_id]][i] ) continue ;
      
	kHistograms[tags[obs_id]][i]->Scale(kXSec * kConversionFactorCm2ToMicroBarn * TMath::Power(10., -38.) / GetNEventsToRun());
      }
    }
  }
  return true;
}

bool MCCLAS6AnalysisI::StoreTree(Event event)
{
  static bool isFirstCall = true;

  // Storing general variables:
  double ConversionFactor = kConversionFactorCm2ToMicroBarn * TMath::Power(10., -38.);
  double TotalXSec = kXSec;
  double MCNormalization = TotalXSec * ConversionFactor / GetNEventsToRun();

  // Storing additional GENIE information:
  bool CC = event.IsCC();
  bool NC = event.IsNC();
  bool EM = event.IsEM();
  bool QEL = event.IsQEL();
  bool RES = event.IsRES();
  bool MEC = event.IsMEC();
  bool DIS = event.IsDIS();
  double TrueQ2s = event.GetTrueQ2s();
  double TrueWs = event.GetTrueWs();
  double Truexs = event.GetTruexs();
  double Trueys = event.GetTrueys();
  double TrueQ2 = event.GetTrueQ2();
  double TrueW = event.GetTrueW();
  double Truex = event.GetTruex();
  double Truey = event.GetTruey();
  int resid = event.GetRESID();

  unsigned int TrueNProtons = event.GetTrueNProtons();
  unsigned int TrueNNeutrons = event.GetTrueNNeutrons();
  unsigned int TrueNPiP = event.GetTrueNPiP();
  unsigned int TrueNPiM = event.GetTrueNPiM();
  unsigned int TrueNPi0 = event.GetTrueNPi0();
  unsigned int TrueNCh = event.GetTrueNCh();
  unsigned int TrueNKP = event.GetTrueNKP();
  unsigned int TrueNKM = event.GetTrueNKM();
  unsigned int TrueNK0 = event.GetTrueNK0();
  unsigned int TrueNEM = event.GetTrueNEM();
  unsigned int TrueNOther = event.GetTrueNOther();

  // Analysis specifics:
  double AccWght = event.GetAccWght();
  bool IsBkg = event.IsBkg();

  // // Rotate φ to match CLAS convention (only for MC events)
  TLorentzVector out_mom = event.GetOutLepton4Mom();
  out_mom.SetPhi(out_mom.Phi() + TMath::Pi());
  event.SetOutLeptonKinematics(out_mom);

  //Rotate φ to match CLAS convention (only for MC events)
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
  event.SetFinalParticlesKinematics(hadron_map);

  // Store name of files used
  const char *InputXSecFile = kXSecFile.c_str();

  if (isFirstCall == true)
  {
    kAnalysisTree->Branch("InputXSecFile", &InputXSecFile, "InputXSecFile/C");
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
    kAnalysisTree->Branch("AccWght", &AccWght, "AccWght/D");
    kAnalysisTree->Branch("IsBkg", &IsBkg, "IsBkg/O");
    kAnalysisTree->Branch("TrueNProtons", &TrueNProtons, "TrueNProtons/I");
    kAnalysisTree->Branch("TrueNNeutrons", &TrueNNeutrons, "TrueNNeutrons/I");
    kAnalysisTree->Branch("TrueNPiP", &TrueNPiP, "TrueNPiP/I");
    kAnalysisTree->Branch("TrueNPiM", &TrueNPiM, "TrueNPiM/I");
    kAnalysisTree->Branch("TrueNPi0", &TrueNPi0, "TrueNPi0/I");
    kAnalysisTree->Branch("TrueNCh", &TrueNCh, "TrueNCh/I");
    kAnalysisTree->Branch("TrueNKP", &TrueNKP, "TrueNKP/I");
    kAnalysisTree->Branch("TrueNKM", &TrueNKM, "TrueNKM/I");
    kAnalysisTree->Branch("TrueNK0", &TrueNK0, "TrueNK0/I");
    kAnalysisTree->Branch("TrueNEM", &TrueNEM, "TrueNEM/I");
    kAnalysisTree->Branch("TrueNOther", &TrueNOther, "TrueNOther/I");
    kAnalysisTree->Branch("resid", &resid, "resid/I");

    // Add normaization information
    kAnalysisTree->Branch("TotalXSec", &TotalXSec, "TotalXSec/D");
    kAnalysisTree->Branch("MCNormalization", &MCNormalization, "MCNormalization/D");

    isFirstCall = false;
  }

  // Store other branches which are generic to CLAS data also:
  AnalysisI::StoreTree(event);

  return true;
}
