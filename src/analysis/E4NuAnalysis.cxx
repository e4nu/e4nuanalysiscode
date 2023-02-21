// _______________________________________________
/*
 * Analysis Interface base class
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

E4NuAnalysis::E4NuAnalysis( const std::string conf_file ) : AnalysisI(conf_file), MCAnalysisI(), CLASAnalysisI() { this->Initialize();}

E4NuAnalysis::E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) : AnalysisI(EBeam, TargetPdg), MCAnalysisI(), CLASAnalysisI() { this->Initialize();}

E4NuAnalysis::~E4NuAnalysis() {
  delete fElectronFit ; 
}

bool E4NuAnalysis::LoadData(void) {
  if( !IsConfigured() ) {
    std::cout << "ERROR: Configuration failed" <<std::endl;
    return false ;
  }
//  if( IsData() ) return CLASAnalysisI::LoadData(file);
  return MCAnalysisI::LoadData() ; 
}

EventI * E4NuAnalysis::GetValidEvent( const unsigned int event_id ) {
  //if( IsData() ) CLASAnalysisI::GetEvent( event_id ) ; 
  return MCAnalysisI::GetValidEvent( event_id ) ; 
}

unsigned int E4NuAnalysis::GetNEvents( void ) const {
  return MCAnalysisI::GetNEvents() ;
}

bool E4NuAnalysis::Analyse(void) {
  unsigned int total_nevents = GetNEvents() ;
  // Loop over events
  for( unsigned int i = 0 ; i < total_nevents ; ++i ) {
    //Print percentage
    if( i==0 || i % 100000 == 0 ) utils::PrintProgressBar(i, total_nevents);
  
    std::unique_ptr<EventI> event = std::unique_ptr<EventI>( (EventI*) MCAnalysisI::GetValidEvent(i) ); 
    if( ! event ) {
      continue ;
    }
  }
  
  return true ; 
}

bool E4NuAnalysis::SubtractBackground() {
  return MCAnalysisI::SubtractBackground() ; 
  // Add for class
} 

bool E4NuAnalysis::Finalise( ) {
  
  bool is_ok = MCAnalysisI::Finalise() ; 

  // Store histograms 
  for( unsigned int i = 0 ; i < kHistograms.size() ; ++i ) {
    kHistograms[i]->GetXaxis()->SetTitle(GetObservablesTag()[i].c_str()) ; 
    if( NormalizeHist() ) kHistograms[i]->GetYaxis()->SetTitle(("d#sigma/d"+GetObservablesTag()[i]).c_str()) ; 
    else {
      if( ApplyCorrWeights() ) kHistograms[i]->GetYaxis()->SetTitle("Weighted Events * weight") ;  
      else kHistograms[i]->GetYaxis()->SetTitle("UnWeighted #Events") ;  
    }
    kHistograms[i]->SetStats(false); 
    kHistograms[i]->Write() ; 
  }
  kAnalysisTree->Write() ; 
  kOutFile->Close() ;
  std::string out_file = GetOutputFile()+".txt";

  return is_ok ; 
}

void E4NuAnalysis::Initialize(void) {
  kOutFile = std::unique_ptr<TFile>( new TFile( (GetOutputFile()+".root").c_str(),"RECREATE") );
  double Ebeam = GetConfiguredEBeam() ; 

  for( unsigned int i = 0 ; i < GetObservablesTag().size() ; ++i ) {
    kHistograms.push_back( new TH1D( GetObservablesTag()[i].c_str(),GetObservablesTag()[i].c_str(), GetNBins()[i], GetRange()[i][0], GetRange()[i][1] ) ) ; 
  }  

  fElectronFit = new TF1( "myElectronFit", "[0]+[1]/x",0.,0.5);
  if( Ebeam == 1.161 ) { fElectronFit -> SetParameters(17,7) ; }
  if( Ebeam == 2.261 ) { fElectronFit -> SetParameters(16,10.5) ; }
  if( Ebeam == 4.461 ) { fElectronFit -> SetParameters(13.5,15) ; }

}

