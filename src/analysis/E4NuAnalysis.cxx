// _______________________________________________
/*
 * Main E4Nu analysis code
 * It deals with the distinction between data and MC
 * It is also responsible for the signal/bkg selection
 * And instiantates the Background substraction method
 * 
 * Analysis Type ID : 
 * 0 => Generic analysis
 * Add new id list here...
 */
#include <iostream>
#include "analysis/E4NuAnalysis.h"
#include "conf/ParticleI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/ParticleUtils.h"
#include "utils/Utils.h"

using namespace e4nu ; 

E4NuAnalysis::E4NuAnalysis() {this->Initialize();}

E4NuAnalysis::E4NuAnalysis( const std::string conf_file ) : AnalysisI(conf_file), MCCLAS6StandardAnalysis(), CLAS6StandardAnalysis() { this->Initialize();}

E4NuAnalysis::E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) : AnalysisI(EBeam, TargetPdg), MCCLAS6StandardAnalysis(), CLAS6StandardAnalysis() { this->Initialize();}

E4NuAnalysis::~E4NuAnalysis() {;}

bool E4NuAnalysis::LoadData(void) {
  if( !IsConfigured() ) {
    std::cout << "ERROR: Configuration failed" <<std::endl;
    return false ;
  }
  // Include new analysis classes with the corresponding analysis ID:
  if( IsCLAS6Analysis() ) { 
    if( IsData() ) {
      if( GetAnalysisTypeID() == 0 ) return CLAS6StandardAnalysis::LoadData();
    }else {
      if( GetAnalysisTypeID() == 0 ) return MCCLAS6StandardAnalysis::LoadData() ; 
      else return false ; 
    } if( IsCLAS12Analysis() ) return false ;  
    return false ; 
  }
  return true ; 
}

EventI * E4NuAnalysis::GetValidEvent( const unsigned int event_id ) {
  // Include new analysis classes with the corresponding analysis ID:
  if( IsCLAS6Analysis() ) { 
    if( IsData() ) {
      if( GetAnalysisTypeID() == 0 ) return CLAS6StandardAnalysis::GetValidEvent( event_id ) ; 
    } else{ 
      if( GetAnalysisTypeID() == 0 ) return MCCLAS6StandardAnalysis::GetValidEvent( event_id ) ; 
      else return nullptr ; 
    }
  } else if ( IsCLAS12Analysis() ) return nullptr ; 
  return nullptr; 
}

unsigned int E4NuAnalysis::GetNEvents( void ) const {
  // Include new analysis classes with the corresponding analysis ID:
  if( IsCLAS6Analysis() ) {
    if( IsData() ) {
      if( GetAnalysisTypeID() == 0 ) return CLAS6StandardAnalysis::GetNEvents() ;
    } else { 
      if( GetAnalysisTypeID() == 0 ) return MCCLAS6StandardAnalysis::GetNEvents() ;
      else return false ; 
    }
  } 
  return 0 ; 
}

bool E4NuAnalysis::Analyse(void) {
  unsigned int total_nevents = GetNEvents() ;
  // Loop over events
  for( unsigned int i = 0 ; i < total_nevents ; ++i ) {
    //Print percentage
    if( i==0 || i % 100000 == 0 ) utils::PrintProgressBar(i, total_nevents);
  
    // Get valid event after analysis
    // It returns cooked event, with detector effects
    EventI * event = nullptr ;
    // Include new analysis classes with the corresponding analysis ID:
    if( IsCLAS6Analysis() ) {
      if( IsData() ) {
	if( GetAnalysisTypeID() == 0 ) event = (EventI*) CLAS6StandardAnalysis::GetValidEvent(i) ; 
      } else {
	if( GetAnalysisTypeID() == 0 ) event = (EventI*) MCCLAS6StandardAnalysis::GetValidEvent(i) ; 
	else event = nullptr ; 
      }
    } 

    if( ! event ) {
      continue ;
    }

    this->ClassifyEvent( event ) ; // Classify events as signal or Background

  }  
  return true ; 
}

void E4NuAnalysis::ClassifyEvent( EventI * event ) { 
  // Classify as signal or background based on topology
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

  // Store in AnalysedEventHolder
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

    if( !is_signal_bkg ) return ; 

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
  }
  return ; 
}

bool E4NuAnalysis::SubtractBackground() {
  
  if( ! BackgroundI::BackgroundSubstraction( kAnalysedEventHolder ) ) return false ;  
  // if( ! BackgroundI::HadronsAcceptanceCorrection( kAnalysedEventHolder ) ) return false ; 
  //  if( ! BackgroundI::ElectronAcceptanceCorrection( kAnalysedEventHolder ) ) return false ; 

  return true ; 
} 

bool E4NuAnalysis::Finalise( ) {
  
  bool is_ok = true ; 
  if( IsCLAS6Analysis() ) {
    if( IsData() ) is_ok = CLAS6AnalysisI::Finalise(kAnalysedEventHolder) ; 
    else {
      if( GetAnalysisTypeID() == 0 ) is_ok = MCCLAS6StandardAnalysis::Finalise(kAnalysedEventHolder) ; 
    }
  }

  unsigned int hist_size = GetObservablesTag().size() ; 
  if( GetDebugBkg() ) hist_size = kHistograms.size() ; 

  if( is_ok ) { 
    for( unsigned int i = 0 ; i < hist_size ; ++i ) {
      if( !kHistograms[i] ) continue ; 
      if( i < GetObservablesTag().size() ){ 
	kHistograms[i]->GetXaxis()->SetTitle(GetObservablesTag()[i].c_str()) ; 
	if( NormalizeHist() ) kHistograms[i]->GetYaxis()->SetTitle(("d#sigma/d"+GetObservablesTag()[i]).c_str()) ; 
	else {
	  kHistograms[i]->GetYaxis()->SetTitle("Weighted Events") ;  
	}
      } else { 
	kHistograms[i]->GetXaxis()->SetTitle("ECal") ;
	kHistograms[i]->GetYaxis()->SetTitle("Weighted Events") ;
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

  unsigned int ECal_id = 0 ;
  for( unsigned int i = 0 ; i < GetObservablesTag().size() ; ++i ) {
    kHistograms.push_back( new TH1D( GetObservablesTag()[i].c_str(),GetObservablesTag()[i].c_str(), GetNBins()[i], GetRange()[i][0], GetRange()[i][1] ) ) ; 
    if( GetObservablesTag()[i] == "ECal" ) ECal_id = i ; 
  }  

  if( GetNBins()[ECal_id] != 0 && GetDebugBkg() ) {
    // These histograms are used to debug the background
    // It compares the true background distribution to the estimated background distribution
    // Different contributions are considered, depending on the multiplicity or the topology
    // You can find some examples here: https://docdb.lns.mit.edu/e4nudb/0000/000060/001/NewBackgroundMethod_test.pdf

    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_OnlySignal").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_signal = kHistograms.size() -1 ;
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_SignalAccCorr").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_acccorr = kHistograms.size() -1 ; 

    // True Background -> Signal 
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotTrueBkg").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_tottruebkg = kHistograms.size() -1 ;

    unsigned int min_mult = GetMinBkgMult() ; 
    unsigned int max_mult = GetMaxBkgMult() ; 
    unsigned int mult = min_mult + 1 ;
    for( unsigned int j = 0 ; j < max_mult - min_mult ; ++j ) { 
      kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotTrueBkg_mult_"+std::to_string(mult)).c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
      mult += 1 ; 
    }
   
    // Add plots for min_mult + 1, in terms of particle content
    // Hardcoded for simplicity ... Adding few cases
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotTrueBkg_mult_2_2p0pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_2p0pitruebkg = kHistograms.size() -1 ;
    
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotTrueBkg_mult_2_1p1pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ;
    kid_1p1pitruebkg = kHistograms.size() -1 ; 
    
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotTrueBkg_mult_3_2p1pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_2p1pitruebkg = kHistograms.size() -1 ;
    
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotTrueBkg_mult_3_1p2pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_1p2pitruebkg = kHistograms.size() -1 ;
  
  
    // Estimated background correction from background events
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotEstBkg").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_totestbkg = kHistograms.size() -1 ; 

    mult = min_mult + 1 ;
    for( unsigned int j = 0 ; j < max_mult - min_mult ; ++j ) { 
      kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotEstBkg_mult_"+std::to_string(mult)).c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
      mult += 1 ; 
    }

    // Add plots for min_mult + 1, in terms of particle content
    // Hardcoded for simplicity ... Adding few cases
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotEstBkg_mult_2_2p0pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_2p0piestbkg = kHistograms.size() -1 ;
    
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotEstBkg_mult_2_1p1pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ;
    kid_1p1piestbkg = kHistograms.size() -1 ; 
    
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotEstBkg_mult_3_2p1pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_2p1piestbkg = kHistograms.size() -1 ;
    
    kHistograms.push_back( new TH1D( (GetObservablesTag()[ECal_id]+"_TotEstBkg_mult_3_1p2pi").c_str(),GetObservablesTag()[ECal_id].c_str(), GetNBins()[ECal_id], GetRange()[ECal_id][0], GetRange()[ECal_id][1] ) ) ; 
    kid_1p2piestbkg = kHistograms.size() -1 ;
  }
}
  
  
