// _______________________________________________
/*
 * Main E4Nu analysis code
 * It deals with the distinction between data and MC
 * It is also responsible for the signal/bkg selection
 * And instiantates the Background substraction method
 * 
 */
#include <iostream>
#include "analysis/E4NuAnalysis.h"
#include "conf/ParticleI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/Utils.h"

using namespace e4nu ; 

E4NuAnalysis::E4NuAnalysis() {this->Initialize();}

E4NuAnalysis::E4NuAnalysis( const std::string conf_file ) : AnalysisI(conf_file), MCAnalysisI(), CLAS6AnalysisI() { this->Initialize();}

E4NuAnalysis::E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) : AnalysisI(EBeam, TargetPdg), MCAnalysisI(), CLAS6AnalysisI() { this->Initialize();}

E4NuAnalysis::~E4NuAnalysis() {;}

bool E4NuAnalysis::LoadData(void) {
  if( !IsConfigured() ) {
    std::cout << "ERROR: Configuration failed" <<std::endl;
    return false ;
  }
  if( IsData() ) return CLAS6AnalysisI::LoadData();
  return MCAnalysisI::LoadData() ; 
}

EventI * E4NuAnalysis::GetValidEvent( const unsigned int event_id ) {
  if( IsData() ) return CLAS6AnalysisI::GetValidEvent( event_id ) ; 
  return MCAnalysisI::GetValidEvent( event_id ) ; 
}

unsigned int E4NuAnalysis::GetNEvents( void ) const {
  if( IsData() ) return CLAS6AnalysisI::GetNEvents() ;
  return MCAnalysisI::GetNEvents() ;
}

bool E4NuAnalysis::Analyse(void) {
  unsigned int total_nevents = GetNEvents() ;
  // Loop over events
  for( unsigned int i = 0 ; i < total_nevents ; ++i ) {
    //Print percentage
    if( i==0 || i % 100000 == 0 ) utils::PrintProgressBar(i, total_nevents);
  
    // Get valid event after analysis
    // It returns cooked event, with detector effects
    EventI * event ;
    if( IsData() ) event = (EventI*) CLAS6AnalysisI::GetValidEvent(i) ; 
    else event = (EventI*) MCAnalysisI::GetValidEvent(i) ; 
  
    if( ! event ) {
      continue ;
    }

    // Classify events as signal or Background
    bool is_signal = true ;
    std::map<int,std::vector<TLorentzVector>> part_map = event -> GetFinalParticles4Mom() ;
    std::map<int,unsigned int> Topology = GetTopology();
    //Topology ID
    for( auto it = Topology.begin() ; it != Topology.end() ; ++it ) {
      if( it->first == conf::kPdgElectron ) continue ; 
      else if( part_map[it->first].size() != it->second ) {
	is_signal = false ; 
	break ; 
      } 
    }

    unsigned int signal_mult = GetMinBkgMult() ;  
    if( is_signal ) {
      // Storing in background the signal events
      if( kAnalysedEventHolder.find(signal_mult) == kAnalysedEventHolder.end() ) {
	std::vector<EventI*> temp ( 1, event ) ;
	kAnalysedEventHolder[signal_mult] = temp ; 
      } else { 
	kAnalysedEventHolder[signal_mult].push_back( event ) ; 
      }
    } else { // BACKGROUND 
      event->SetIsBkg(true); 

      // Get Number of signal particles, "multiplicity"
      unsigned int mult_bkg = event->GetNSignalParticles( part_map, Topology ) ; 

      // Check events are above minumum particle multiplicity 
      bool is_signal_bkg = true ; 
      std::map<int,std::vector<TLorentzVector>> hadrons = event->GetFinalParticles4Mom() ;
      for( auto it = Topology.begin(); it!=Topology.end();++it){
	if( it->first == conf::kPdgElectron ) continue ; 
	for( auto part = hadrons.begin() ; part != hadrons.end() ; ++part ) {
	  if( hadrons.find(it->first) != hadrons.end() && hadrons[it->first].size() < it->second ) { is_signal_bkg = false ; break ; } 
	  if( hadrons.find(it->first) == hadrons.end() &&  it->second != 0 ) { is_signal_bkg = false ; break ; }
	}
      }

      if( !is_signal_bkg ) {
	continue ; 
      }      
      // Only store background events with multiplicity > mult_signal
      // Also ignore background events above the maximum multiplicity
      if( mult_bkg > signal_mult && mult_bkg <= GetMaxBkgMult() ) {
	if( kAnalysedEventHolder.find(mult_bkg) == kAnalysedEventHolder.end() ) {
	  std::vector<EventI*> temp ( 1, event ) ;
	  kAnalysedEventHolder[mult_bkg] = temp ; 
	} else { 
	  kAnalysedEventHolder[mult_bkg].push_back( event ) ; 
	}
      }
      continue ; 
    }
  }  
  return true ; 
}

bool E4NuAnalysis::SubtractBackground() {
  
  if( ! BackgroundI::NewBackgroundSubstraction( kAnalysedEventHolder ) ) return false ;  
  //  if( ! BackgroundI::AcceptanceCorrection( kAnalysedEventHolder ) ) return false ; 

  return true ; 
} 

bool E4NuAnalysis::Finalise( ) {
  
  bool is_ok = true ; 
  if( IsData() ) is_ok = CLAS6AnalysisI::Finalise(kAnalysedEventHolder) ; 
  is_ok = MCAnalysisI::Finalise(kAnalysedEventHolder) ; 

  if( is_ok ) { 
    for( unsigned int i = 0 ; i < kHistograms.size() ; ++i ) {
      kHistograms[i]->GetXaxis()->SetTitle(GetObservablesTag()[i].c_str()) ; 
      if( NormalizeHist() ) kHistograms[i]->GetYaxis()->SetTitle(("d#sigma/d"+GetObservablesTag()[i]).c_str()) ; 
      else {
	kHistograms[i]->GetYaxis()->SetTitle("Weighted Events * weight") ;  
      }
      kHistograms[i]->SetStats(false); 
      kHistograms[i]->Write() ; 
    }
  }

  kAnalysisTree->Write() ; 

  kOutFile->Close() ;
  std::string out_file = GetOutputFile()+".txt";

  return is_ok ; 
}

void E4NuAnalysis::Initialize(void) {
  kOutFile = std::unique_ptr<TFile>( new TFile( (GetOutputFile()+".root").c_str(),"RECREATE") );

  for( unsigned int i = 0 ; i < GetObservablesTag().size() ; ++i ) {
    kHistograms.push_back( new TH1D( GetObservablesTag()[i].c_str(),GetObservablesTag()[i].c_str(), GetNBins()[i], GetRange()[i][0], GetRange()[i][1] ) ) ; 
  }  

}

