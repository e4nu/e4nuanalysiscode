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
#include "utils/TargetUtils.h"
#include "utils/DetectorUtils.h"
#include "conf/FiducialCutI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/CLAS6ConstantsI.h"
#include "conf/ConstantsI.h"
#include "conf/CLAS6RunRequirements.h"

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

  // check run is valid. This is used to manually remove runs if needed 
  if ( ! IsRunValid( event->GetEventRunNumber() ) ) { 
    delete event;
    return nullptr;
  }
  
  // Double checking we have the run in the root file (as it is easy to miss one by eye...)
  SetIsRunInFile( event->GetEventRunNumber() ) ; 


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

  // Remove runs which do  not have the right beam energy or target
  RetrieveAllRuns( GetConfiguredEBeam(), utils::GetTargetName(GetConfiguredTarget()) ) ; 
  // Remove additional runs due to issues in the run
  NeglectRuns( conf::GetInvalidRuns( GetConfiguredEBeam(), utils::GetTargetName(GetConfiguredTarget()) ) ) ;  

}

bool CLAS6AnalysisI::Finalise(std::map<int, std::vector<e4nu::Event>> &event_holder)
{

  if (!AnalysisI::Finalise())
    return false;

  auto tags = GetObservablesTag();
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
    for (unsigned int j = 0; j < tags.size(); ++j)
    {
      if( !kHistograms[tags[j]][0] ) continue ;
      kHistograms[tags[j]][0]->Fill(event_holder[min_mult][k].GetObservable(tags[j]), norm_weight);
    }
  }

  // Normalize
  unsigned int tgt_pdg = GetConfiguredTarget();
  // Get constants
  unsigned int MassNumber = utils::GetMassNumber(tgt_pdg);
  // Now we compute directly from the table; It was computed incorrectly in old method !! conf::GetIntegratedCharge(tgt_pdg, EBeam);
  double IntegratedCharge = GetTotalIntegratedCharge(); 
  std::cout << " Total integrated charge " << IntegratedCharge << std::endl;
  double TargetLength = conf::GetTargetLength(tgt_pdg);
  double TargetDensity = conf::GetTargetDensity(tgt_pdg);

  if (NormalizeHist())
  {
    for (unsigned int j = 0; j < tags.size(); ++j)
    {
      if( !kHistograms[tags[j]][0] ) continue ;
      double NBins = kHistograms[tags[j]][0]->GetNbinsX();

      for (int k = 1; k <= NBins; k++)
      {
        double content = kHistograms[tags[j]][0]->GetBinContent(k);
        double error = kHistograms[tags[j]][0]->GetBinError(k);
        double width = kHistograms[tags[j]][0]->GetBinWidth(k);
        double newcontent = content / width;
        double newerror = error / width;
        kHistograms[tags[j]][0]->SetBinContent(k, newcontent);
        kHistograms[tags[j]][0]->SetBinError(k, newerror);
      }
      kHistograms[tags[j]][0]->Scale(kConversionFactorCm2ToMicroBarn * MassNumber / (IntegratedCharge * TargetLength * TargetDensity * kOverallUnitConversionFactor));
    }
  }

  return true;
}

bool CLAS6AnalysisI::StoreTree(Event event)
{
  static bool isFirstCall = true;
  // Storing general variables:
  double BeamE = event.GetInLepton4Mom().E();
  int TargetPdg = event.GetTargetPdg();
  unsigned int MassNumber = utils::GetMassNumber(TargetPdg);
  double IntegratedCharge = GetTotalIntegratedCharge();
  double TargetLength = conf::GetTargetLength(TargetPdg);
  double TargetDensity = conf::GetTargetDensity(TargetPdg);
  double ConversionFactor = kConversionFactorCm2ToMicroBarn / kOverallUnitConversionFactor;
  double DataNormalization = kConversionFactorCm2ToMicroBarn * MassNumber / (IntegratedCharge * TargetLength * TargetDensity * kOverallUnitConversionFactor);
  unsigned int RunNumber = event.GetEventRunNumber();
  
  if( isFirstCall ){
    kAnalysisTree->Branch("MassNumber", &MassNumber, "MassNumber/I");
    kAnalysisTree->Branch("IntegratedCharge", &IntegratedCharge, "IntegratedCharge/D");
    kAnalysisTree->Branch("TargetLength", &TargetLength, "TargetLength/D");
    kAnalysisTree->Branch("TargetDensity", &TargetDensity, "TargetDensity/D");
    kAnalysisTree->Branch("ConversionFactor", &ConversionFactor, "ConversionFactor/D");
    kAnalysisTree->Branch("DataNormalization", &DataNormalization, "DataNormalization/D");
    kAnalysisTree->Branch("RunNumber", &RunNumber, "RunNumber/I");
    isFirstCall = false;
  }
  // Store other branches which are generic to CLAS data also:
  AnalysisI::StoreTree(event);

  return true;
}
