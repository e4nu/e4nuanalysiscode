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
  delete fRotation ;
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
    }
    --m; 
  }

  std::map<int,std::vector<TLorentzVector>> particles ;
  std::map<int,std::vector<TLorentzVector>> particles_uncorr ; 
  TLorentzVector V4_el ;
  unsigned int bkg_mult = min_mult + 1 ;
 
  // remove multiplicity 2 contribution to signal...
  if ( fBkg.find(bkg_mult) != fBkg.end() ) {
     if ( fBkg[bkg_mult].size() != 0 ) {
      for ( unsigned int i = 0 ; i < fBkg[bkg_mult].size() ; ++i ) {
	particles = fBkg[bkg_mult][i].GetFinalParticles4Mom();
	particles_uncorr = fBkg[bkg_mult][i].GetFinalParticlesUnCorr4Mom(); // This map needs to change with the cuts as well...
	
	// 2p0pi -> 1p0pi ( multiplicity 2 -> multiplicity 1 ) 
	if( particles[conf::kPdgProton].size() == 2 
	    && particles[conf::kPdgPiP].size() == 0 
	    && particles[conf::kPdgPiM].size() == 0 
	    && particles[conf::kPdgPi0].size() == 0 
	    && particles[conf::kPdgPhoton].size() == 0 ) {

	  V4_el = fBkg[bkg_mult][i].GetOutLepton4Mom();
	  
	  double E_tot_2p[bkg_mult]={0};
	  double p_perp_tot_2p[bkg_mult]={0};
	  double N_prot_both = 0;
	  double P_N_2p[bkg_mult]={0};
	  
	  TVector3 V3_2prot_corr[bkg_mult];
	  for ( unsigned k = 0 ; k < bkg_mult ; ++k ) {
	    V3_2prot_corr[k] = particles[conf::kPdgProton][k].Vect() ; 
	  }
	  
	  TVector3 V3_2prot_uncorr[bkg_mult];
	  for ( unsigned k = 0 ; k < bkg_mult ; ++k ) {
	    V3_2prot_uncorr[k] = particles_uncorr[conf::kPdgProton][k].Vect() ; 
	  }

	  fRotation->SetQVector( fBkg[bkg_mult][i].GetRecoq3() );	
	  fRotation->prot2_rot_func( V3_2prot_corr, V3_2prot_uncorr, V4_el, E_tot_2p, p_perp_tot_2p, P_N_2p , &N_prot_both);
	  
	  for( unsigned int j = 0 ; j < bkg_mult ; ++j ) {
	   fBkg[bkg_mult][i].SetEventWeight( -P_N_2p[j] ) ; 
	   fBkg[bkg_mult][i].SetIsBkg(false); // For now, we change to signal... with negative weight
	   if ( fBkg.find(min_mult) != fBkg.end() ) {
	     fBkg[min_mult].push_back( fBkg[bkg_mult][i] ) ; 
	   } else {
	     std::vector<e4nu::EventI> temp = { fBkg[bkg_mult][i] } ;
	     fBkg[min_mult] = temp ; 
	   }
	  }
	}
      } 
     particles.clear() ;
     particles_uncorr.clear() ;
     }
  }

  // Store corrected background in event sample

  for( unsigned int k = 0 ; k < fBkg[min_mult].size() ; ++k ) {
    // if( IsData() ) 
    // MCAnalysisI::StoreTree( (MCEventI) fBkg[min_mult][k] ) ;  
    for( unsigned int j = 0 ; j < kHistograms.size() ; ++j ) {
      // Store in histogram
      kHistograms[j]-> Fill( fBkg[min_mult][k].GetObservable( GetObservablesTag()[j] ) , fBkg[min_mult][k].GetTotalWeight() ) ;
    }
  }

  return true ; 
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

  fRotation = new Subtraction();
  fRotation->InitSubtraction( Ebeam, GetConfiguredTarget(), GetNRotations(), kFiducialCut);
  fRotation->ResetQVector(); //Resets q vector to (0,0,0)
}

