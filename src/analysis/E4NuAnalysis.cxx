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

E4NuAnalysis::E4NuAnalysis( const std::string conf_file ) : ConfigureI(conf_file), MCAnalysisI(), CLASAnalysisI() { this->Initialize();}

E4NuAnalysis::E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) : ConfigureI(EBeam, TargetPdg), MCAnalysisI(), CLASAnalysisI() { this->Initialize();}

E4NuAnalysis::~E4NuAnalysis() {

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

    TLorentzVector in_mom = event -> GetInLepton4Mom() ;
    TLorentzVector out_mom = event -> GetOutLepton4Mom() ;
    double EBeam = in_mom.E() ; 

    if( out_mom.Theta() * 180 / TMath::Pi() < GetElectronMinTheta( out_mom ) ) continue ;
    ++fNEventsAfterEThetaCut ;     

    if( ApplyOutElectronCut() ){
      if( out_mom.P() < conf::GetMinMomentumCut( conf::kPdgElectron, EBeam ) ) continue ; 
    }
    ++fNEventsAfterEMomCut ;     

    if( ApplyThetaSlice() ) {
      if( out_mom.Theta() * 180./TMath::Pi() < conf::kMinEThetaSlice ) continue ; 
      if( out_mom.Theta() * 180./TMath::Pi() > conf::kMaxEThetaSlice ) continue ; 
      ++fNEventsAfterThetaCut ; 
    }

    if( ApplyPhiOpeningAngle() ) {
      double phi = out_mom.Phi() ; 
      if ( ! IsData() ) phi += TMath::Pi() ; 
      if ( ! conf::ValidPhiOpeningAngle( phi ) ) continue ;  
    }
    ++fNEventsAfterPhiOpeningAngleCut ; 

    if( ApplyGoodSectorPhiSlice() ) {
      double phi = out_mom.Phi() ; 
      if ( ! IsData() ) phi += TMath::Pi() ; 
      if ( ! conf::GoodSectorPhiSlice( phi ) ) continue ; 
    }
    ++fNEventsAfterPhiCut ; 

    double reco_Q2 = utils::GetRecoQ2( out_mom, EBeam ) ; 
    double W_var = utils::GetRecoW( out_mom, EBeam ) ;

    if( ApplyQ2Cut() ) {
      double MaxQ2 = 0 ; 
      if( conf::GetQ2Cut( MaxQ2, EBeam ) ) {
	if( reco_Q2 < MaxQ2 ) continue ; 
      }
      ++fNEventsAfterQ2Cut ; 
    }

    if( ApplyWCut() ) {
      double MinW = 0 ; 
      if( conf::GetWCut( MinW, EBeam ) ) {
	if( W_var > MinW ) continue ; 
      }
      ++fNEventsAfterWCut ; 
    }
    
    for( unsigned int j = 0 ; j < kHistograms.size() ; ++j ) {
      kHistograms[j]-> Fill( event-> GetObservable( GetObservablesTag()[j] ) , event->GetTotalWeight() ) ;
    }
  }
  
  return true ; 
}

bool E4NuAnalysis::SubstractBackground(void) {
  unsigned int max_mult = GetMaxBkgMult(); 
  unsigned int min_mult = GetMinBkgMult(); // Signal multiplicity
  std::map<int,unsigned int> Topology = GetTopology();
  
  unsigned int m = max_mult ;
  while ( m >= min_mult ) {
    if( fBkg.find(m) != fBkg.end() ) {
      std::cout<< " Number of events with multiplicity " << m << " = " << fBkg[m].size() <<std::endl; 
      // Calculate weight for events with multiplicity m->m-1
      // Add events in fBkg(m-1)
      for( auto i = 0 ; i < fBkg[m].size() ; ++i ) {
	std::cout << " m = " << m << " id = " << fBkg[m][i].GetEventID() << std::endl;
      }
    }
    --m ; 
  }
  
} 


bool E4NuAnalysis::Finalise( ) {
  
  bool is_ok = MCAnalysisI::Finalise() ; 

  // Store histograms 
  for( unsigned int i = 0 ; i < kHistograms.size() ; ++i ) {
    kHistograms[i]->GetXaxis()->SetTitle(GetObservablesTag()[i].c_str()) ; 
    if( NormalizeHist() ) kHistograms[i]->GetYaxis()->SetTitle(("d#sigma/d"+GetObservablesTag()[i]).c_str()) ; 
    else kHistograms[i]->GetYaxis()->SetTitle("NEvents * weight") ;  
    kHistograms[i]->SetStats(false); 
    kHistograms[i]->Write() ; 
  }
  kAnalysisTree->Write() ; 
  kOutFile->Close() ;
  std::string out_file = GetOutputFile()+".txt";

  std::cout << " Events after electron momentum cut = " << fNEventsAfterEMomCut << std::endl;
  std::cout << " Events after electron theta cut = " << fNEventsAfterEThetaCut << std::endl;
  if( ApplyQ2Cut() ) std::cout << " Events after Q2 cut = " << fNEventsAfterQ2Cut << std::endl;
  if( ApplyWCut() ) std::cout << " Events after W cut = "<< fNEventsAfterWCut <<std::endl;
  if( ApplyThetaSlice() ) std::cout << " Events after theta cut = " << fNEventsAfterThetaCut << std::endl;
  if( ApplyPhiOpeningAngle() ) std::cout << " Events after phi opening angle cut = " << fNEventsAfterPhiOpeningAngleCut << std::endl;
  if( ApplyGoodSectorPhiSlice() ) std::cout << " Events after phi cut = " << fNEventsAfterPhiCut << std::endl;

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

double E4NuAnalysis::GetElectronMinTheta( TLorentzVector emom ) {
  return fElectronFit ->Eval(emom.P()) ; 
}
