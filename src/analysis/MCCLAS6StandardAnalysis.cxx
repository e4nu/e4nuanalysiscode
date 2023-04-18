/*
 * MC Analysis standard analysis
 * This class is a template for specific e4nu analysis
 * wich require additional constrains from the standard implementation
 * implemented in MCCLAS6AnalysisI
 * 
 */
#include <iostream>
#include "TFile.h"
#include "TDirectoryFile.h"
#include "analysis/MCCLAS6StandardAnalysis.h"
#include "utils/ParticleUtils.h"
#include "utils/KinematicUtils.h"
#include "conf/ParticleI.h"
#include "conf/TargetI.h"
#include "utils/DetectorUtils.h"
#include "conf/FiducialCutI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/ConstantsI.h"

using namespace e4nu ; 

MCCLAS6StandardAnalysis::MCCLAS6StandardAnalysis() {}

MCCLAS6StandardAnalysis::~MCCLAS6StandardAnalysis() {}

EventI * MCCLAS6StandardAnalysis::GetValidEvent( const unsigned int event_id ) {
  MCEvent * event = (MCEvent*) MCCLAS6AnalysisI::GetValidEvent(event_id) ;

  // Add additional constrains here ...
  // Operations ...

  return event; 
}


