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
#include "plotting/PlottingUtils.h"

using namespace e4nu ;

E4NuAnalysis::E4NuAnalysis() {}

E4NuAnalysis::E4NuAnalysis( const std::string conf_file ) : AnalysisI(conf_file), MCCLAS6StandardAnalysis(), CLAS6StandardAnalysis() {}

E4NuAnalysis::E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) : AnalysisI(EBeam, TargetPdg), MCCLAS6StandardAnalysis(), CLAS6StandardAnalysis() {}

E4NuAnalysis::~E4NuAnalysis() { }

bool E4NuAnalysis::LoadData(void) {
  if( !IsConfigured() ) {
    std::cout << "ERROR: Configuration failed" <<std::endl;
    return false ;
  }
  // Include new analysis classes with the corresponding analysis ID:
  if( IsCLAS6Analysis() ) {
    if( IsData() ) {
      if( GetAnalysisTypeID() == 0 ) return CLAS6StandardAnalysis::LoadData();
    } else {
      if( GetAnalysisTypeID() == 0 ) return MCCLAS6StandardAnalysis::LoadData() ;
      else return false ;
    } if( IsCLAS12Analysis() ) return false ;
    return false ;
  }
  return true ;
}

Event * E4NuAnalysis::GetValidEvent( const unsigned int event_id ) {
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
      else return 0 ;
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
    Event * event = nullptr ;
    // Include new analysis classes with the corresponding analysis ID:
    if( IsCLAS6Analysis() ) {
      if( IsData() ) {
        if( GetAnalysisTypeID() == 0 ) event = (Event*) CLAS6StandardAnalysis::GetValidEvent(i) ;
      } else {
        if( GetAnalysisTypeID() == 0 ) event = (Event*) MCCLAS6StandardAnalysis::GetValidEvent(i) ;
        else event = nullptr ;
      }
    }

    if( ! event ) continue ;

    this->ClassifyEvent( *event ) ; // Classify events as signal or Background

    delete event ;
  }
  return true ;
}

void E4NuAnalysis::ClassifyEvent( Event event ) {
  // Classify as signal or background based on topology
  bool is_signal = true ;
  std::map<int,std::vector<TLorentzVector>> part_map = event.GetFinalParticles4Mom() ;
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
      std::vector<Event> temp ( 1, event ) ;
      kAnalysedEventHolder[signal_mult] = temp ;
    } else {
      kAnalysedEventHolder[signal_mult].push_back( event ) ;
    }
  } else { // BACKGROUND
    // No need to classify it if we don't apply fiducial
    if ( !ApplyFiducial() || ! GetSubtractBkg() ) return ;
    // Tag event as Background
    event.SetIsBkg(true);

    // Get Number of signal particles, "multiplicity"
    unsigned int mult_bkg = event.GetNSignalParticles( part_map, Topology ) ;

    // Check events are above minumum particle multiplicity
    bool is_signal_bkg = true ;
    std::map<int,std::vector<TLorentzVector>> hadrons = event.GetFinalParticles4Mom() ;
    for( auto it = Topology.begin(); it!=Topology.end();++it){
      if( it->first == conf::kPdgElectron ) continue ;
      if( hadrons.find(it->first) != hadrons.end() && hadrons[it->first].size() < it->second ) { is_signal_bkg = false ; break ; }
      if( hadrons.find(it->first) == hadrons.end() &&  it->second != 0 ) { is_signal_bkg = false ; break ; }
    }

    if( !is_signal_bkg ) return ;

    // Only store background events with multiplicity > mult_signal
    // Also ignore background events above the maximum multiplicity
    if( mult_bkg > signal_mult && mult_bkg <= GetMaxBkgMult() ) {
      if( kAnalysedEventHolder.find(mult_bkg) == kAnalysedEventHolder.end() ) {
        std::vector<Event> temp ( 1, event ) ;
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

  auto tags = GetObservablesTag() ;
  for( unsigned int obs_id = 0 ; obs_id < tags.size() ; ++obs_id ) {

    unsigned int hist_size = kHistograms[tags[obs_id]].size() ;

    if( is_ok ) {
      for( unsigned int i = 0 ; i < hist_size ; ++i ) {
        if( !kHistograms[tags[obs_id]][i] ) continue ;
    
	double NBins = kHistograms[tags[obs_id]][i]->GetNbinsX();
    
	for (int k = 1; k <= NBins; k++) {
	  double content = kHistograms[tags[obs_id]][i]->GetBinContent(k);
	  double error = kHistograms[tags[obs_id]][i]->GetBinError(k);
	  double width = kHistograms[tags[obs_id]][i]->GetBinWidth(k);
	  double newcontent = content / width;
	  double newerror = error / width;
	  kHistograms[tags[obs_id]][i]->SetBinContent(k, newcontent);
	  kHistograms[tags[obs_id]][i]->SetBinError(k, newerror);
	}

	kHistograms[tags[obs_id]][i]->GetXaxis()->SetTitle(tags[obs_id].c_str()) ;
        if( NormalizeHist() ) {
	  kHistograms[tags[obs_id]][i]->GetYaxis()->SetTitle(("d#sigma/d"+tags[obs_id]).c_str()) ;
        } else {
          kHistograms[tags[obs_id]][i]->GetYaxis()->SetTitle("Weighted Events") ;
        }
        kHistograms[tags[obs_id]][i]->SetStats(false);
        kHistograms[tags[obs_id]][i]->Write() ;
      }
    }
  }
  kAnalysisTree->Write() ;

  kOutFile->Close() ;
  std::string out_file = GetOutputFile()+".txt";

  return is_ok ;
}

void E4NuAnalysis::Initialize(void) {
  kOutFile = std::unique_ptr<TFile> ( new TFile( (GetOutputFile()+".root").c_str(),"RECREATE") ) ;
  auto tags = GetObservablesTag() ;
  for( unsigned int obs_id = 0 ; obs_id < tags.size() ; ++obs_id ) {
    std::vector<double> binning = plotting::GetBinning( tags[obs_id], GetConfiguredEBeam(), GetAnalysisKey() );
    if( binning.size() == 0 ) {
      std::cout << " Issue with the binning and validation plots" << std::endl;
      continue;
    }

    if (kHistograms.find(tags[obs_id]) == kHistograms.end()) {
      kHistograms[tags[obs_id]] = { new TH1D( tags[obs_id].c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) };
    } else {
      kHistograms[tags[obs_id]].push_back( new TH1D( tags[obs_id].c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
    }

    if( GetDebugBkg() ) {

      // These histograms are used to debug the background
      // It compares the true background distribution to the estimated background distribution
      // Different contributions are considered, depending on the multiplicity or the topology
      // You can find some examples here: https://docdb.lns.mit.edu/e4nudb/0000/000060/001/NewBackgroundMethod_test.pdf

      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_OnlySignal").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_signal = kHistograms[tags[obs_id]].size() -1 ;
      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_SignalAccCorr").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_acccorr = kHistograms[tags[obs_id]].size() -1 ;

      // True Background -> Signal
      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotTrueBkg").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_tottruebkg = kHistograms[tags[obs_id]].size() -1 ;

      unsigned int min_mult = GetMinBkgMult() ;
      unsigned int max_mult = GetMaxBkgMult() ;
      unsigned int mult = min_mult + 1 ;
      for( unsigned int j = 0 ; j < max_mult - min_mult ; ++j ) {
        kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotTrueBkg_mult_"+std::to_string(mult)).c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
        mult += 1 ;
      }

      // Add plots for min_mult + 1, in terms of particle content
      // Hardcoded for simplicity ... Adding few cases
      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotTrueBkg_mult_2_2p0pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_2p0pitruebkg = kHistograms[tags[obs_id]].size() -1 ;

      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotTrueBkg_mult_2_1p1pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_1p1pitruebkg = kHistograms[tags[obs_id]].size() -1 ;

      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotTrueBkg_mult_3_2p1pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_2p1pitruebkg = kHistograms[tags[obs_id]].size() -1 ;

      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotTrueBkg_mult_3_1p2pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_1p2pitruebkg = kHistograms[tags[obs_id]].size() -1 ;


      // Estimated background correction from background events
      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotEstBkg").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_totestbkg = kHistograms[tags[obs_id]].size() -1 ;

      mult = min_mult + 1 ;
      for( unsigned int j = 0 ; j < max_mult - min_mult ; ++j ) {
        kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotEstBkg_mult_"+std::to_string(mult)).c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
        mult += 1 ;
      }

      // Add plots for min_mult + 1, in terms of particle content
      // Hardcoded for simplicity ... Adding few cases
      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotEstBkg_mult_2_2p0pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_2p0piestbkg = kHistograms[tags[obs_id]].size() -1 ;

      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotEstBkg_mult_2_1p1pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_1p1piestbkg = kHistograms[tags[obs_id]].size() -1 ;

      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotEstBkg_mult_3_2p1pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_2p1piestbkg = kHistograms[tags[obs_id]].size() -1 ;

      kHistograms[tags[obs_id]].push_back( new TH1D( (tags[obs_id]+"_TotEstBkg_mult_3_1p2pi").c_str(),tags[obs_id].c_str(), binning.size()-1, &binning[0] ) ) ;
      kid_1p2piestbkg = kHistograms[tags[obs_id]].size() -1 ;
    }

    for( unsigned int i = 0 ; i < kHistograms[tags[obs_id]].size() ; ++i ) {
      kHistograms[tags[obs_id]][i]->Sumw2();
    }
  }

}
