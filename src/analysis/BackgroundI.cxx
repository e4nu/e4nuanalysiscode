/** 
 * \info These parameters are configurable 
 * the default values are set here
 **/
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "analysis/BackgroundI.h"
#include "conf/FiducialCutI.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/AccpetanceMapsI.h"
#include "conf/AnalysisCutsI.h"
#include "utils/KinematicUtils.h"
#include "utils/Utils.h"

using namespace e4nu; 

BackgroundI::BackgroundI( ) {

  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 

  if( kIsConfigured && fFiducialCut ) {
    fRotation = new Subtraction();
    fRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), fFiducialCut);
    fRotation->ResetQVector(); //Resets q vector to (0,0,0) 
  }
}

BackgroundI::BackgroundI( const std::string input_file ) : ConfigureI( input_file ) {

  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 

  if( kIsConfigured && fFiducialCut ) {
    fRotation = new Subtraction();
    fRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), fFiducialCut);
    fRotation->ResetQVector(); //Resets q vector to (0,0,0)
  }
}

BackgroundI::BackgroundI( const double EBeam, const unsigned int TargetPdg ) : ConfigureI( EBeam, TargetPdg ) {

  if( kIsConfigured ) kIsConfigured = InitializeFiducial() ; 

  if( kIsConfigured && fFiducialCut ) {
    fRotation = new Subtraction();
    fRotation->InitSubtraction( GetConfiguredEBeam(), GetConfiguredTarget(), GetNRotations(), fFiducialCut);
    fRotation->ResetQVector(); //Resets q vector to (0,0,0)
  }
}    

BackgroundI::~BackgroundI() {
  fBkg.clear() ; 
  delete fRotation ;
  delete fFiducialCut ;
}

bool BackgroundI::InitializeFiducial(void) {
  double EBeam = GetConfiguredEBeam() ; 
  unsigned int Target = GetConfiguredTarget() ;

  if( ApplyFiducial() ) {
    // Initialize fiducial for this run
    fFiducialCut = new Fiducial() ; 
    fFiducialCut -> InitPiMinusFit( EBeam ) ; 
    fFiducialCut -> InitEClimits(); 
    fFiducialCut -> up_lim1_ec -> Eval(60) ;
    fFiducialCut -> SetConstants( conf::GetTorusCurrent( EBeam ), Target , EBeam ) ;
    fFiducialCut -> SetFiducialCutParameters( EBeam ) ;
    if( !fFiducialCut ) return false ; 
  } else { return true ; }

  return true ; 
}

bool BackgroundI::SubstractBackground(void) {

  unsigned int max_mult = GetMaxBkgMult(); 

  if( GetSubstractBkg() == true ) {

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
	    if ( fBkg.find(min_mult) != fBkg.end() ) {
	      fBkg[min_mult].push_back( fBkg[bkg_mult][i] ) ; 
	    } else {
	      std::vector<e4nu::MCEvent> temp = { fBkg[bkg_mult][i] } ;
	      fBkg[min_mult] = temp ; 
	    }
	  }
	
	} 
	particles.clear() ;
	particles_uncorr.clear() ;
      }
    }
  }
  return true ; 
} 
