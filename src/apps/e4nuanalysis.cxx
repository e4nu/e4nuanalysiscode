// _____________________________________________________________
/* This app is used to run a generic e4nu analysis            */
// _____________________________________________________________

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include "TH1D.h"
#include "analysis/E4NuAnalysis.h"
#include "plotting/PlottingUtils.h"

using namespace std; 
using namespace e4nu;
using namespace e4nu::plotting;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////
// --conf-file) Configuration file location                    //   
// --root-file) InputFile                                      // 
// --output-file) OutputFile                                   //
// --analysis-type) Type of analysis:                          // 
//    - ComputeTrueAccCorr                                     // 
//    - ComputeTrueRecoAccCorr                                 //
//    - IsData                                                 //   
// --rad-corr bool ; used to change the output name of file    //
//                   when radiative corr are used in MC        //
// --xsec-file) XSecFile (only for MC)                         //
// --bkg-mult) Add multiplicity                                //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "E4Nu analysis ongoing..." << std::endl;

  // This object can be initialized with a configuration file which contains information on the event run, 
  // cuts and analysis requirements, and output file location
  char * env = std::getenv("E4NUANALYSIS") ; 
  std::string path( env ) ; 
  path += "/" ;
  std::string config_file = "ConfFiles/example_configuration.txt" ;

  if( argc > 1 && ExistArg("conf-file",argc,argv)) {
    config_file = GetArg("conf-file",argc,argv) ;
  }

  E4NuAnalysis * analysis = new E4NuAnalysis((path+config_file).c_str()) ;
  if( ! analysis ) return 0 ; 
  // Compute acceptance correction files
  bool compute_trueacc = analysis -> ComputeTrueAccCorrection() ; 
  bool compute_truerecoacc = analysis -> ComputeTrueRecoAccCorrection() ; 
  bool is_data = analysis -> IsData() ;
 
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("root-file",argc,argv)) {
      analysis -> SetInputFile( GetArg("root-file",argc,argv)); 
    }
    if( ExistArg("analysis-type",argc,argv)) {
      if ( GetArg("analysis-type",argc,argv) == "ComputeTrueAccCorr" ) { 
	std::cout << " ComputeTrueAccCorr = True " <<std::endl;
	compute_trueacc = true ; compute_truerecoacc = false ; 
      }
      else if ( GetArg("analysis-type",argc,argv) == "ComputeTrueRecoAccCorr" ) {
	std::cout << " ComputeRecoAccCorr = True " <<std::endl;
	compute_trueacc = false ; compute_truerecoacc = true ; 
      }
      else if ( GetArg("analysis-type",argc,argv) == "IsData" ) { is_data = true ; compute_trueacc = false ; compute_truerecoacc = false ; }
    }
    if( ExistArg("rad-corr",argc,argv) && !is_data ) { 
      analysis -> SetRadCorr( true ) ; 
    }
    if( ExistArg("xsec-file",argc,argv) && !is_data ) { 
      analysis -> SetXSecFile( GetArg("xsec-file",argc,argv) ) ; 
    } 
    if( ExistArg("bkg-mult",argc,argv) && is_data ) {
      analysis -> SetMaxBkgMult( atoi(GetArg("bkg-mult",argc,argv).c_str()) ) ;
    }
    if( ExistArg("output-file",argc,argv)) {
      std::string final_name = GetArg("output-file",argc,argv) ; 
      if( is_data ) {
	unsigned int max_mult = analysis->GetMaxBkgMult() ; 
	final_name += "_"+std::to_string(max_mult)+"MaxBkgMult";
      }
      analysis -> SetOutputFile( final_name );
    }
  }

  if ( is_data ) {   
    compute_trueacc = false ; 
    compute_truerecoacc = false ; 
    analysis -> SetApplyFiducial( true ) ; 
    analysis -> SetApplyAccWeights( true ) ; 
    analysis -> SetApplyReso( true ) ;  
    std::string OutputFile_data = analysis->GetOutputFile() + "_clas6data" ;
    analysis -> SetOutputFile( OutputFile_data ) ; 
  }

  if ( compute_trueacc ) {
    analysis -> SetTrueSignal( true ) ;
    analysis -> SetApplyFiducial( false ) ;
    analysis -> SetApplyAccWeights( false ) ;
    //    analysis -> SetApplyReso( false ) ;
    analysis -> SetUseAllSectors( true ) ;
    analysis -> EnableAllSectors( true ) ;
    std::string OutputFile_true = analysis->GetOutputFile() + "_true" ;
    if( analysis -> IsRadiated() ) OutputFile_true += "_radcorr";
    analysis -> SetOutputFile( OutputFile_true ) ;
    std::cout << " Computing true analysis distributions for acceptance correction..."<<std::endl;
  } else if ( compute_truerecoacc ) { 
    analysis -> SetTrueSignal( true ) ; 
    analysis -> SetApplyFiducial( true ) ; 
    analysis -> SetApplyAccWeights( true ) ; 
    analysis -> SetApplyReso( true ) ; 
    analysis -> SetSubtractBkg( false ) ; 
    std::string OutputFile_reco = analysis->GetOutputFile() + "_truereco" ;
    if( analysis -> IsRadiated() ) OutputFile_reco += "_radiated";
    analysis -> SetOutputFile( OutputFile_reco ) ; 
    std::cout << " Computing true reconstructed analysis distributions for acceptance correction..."<<std::endl;
  }
  analysis -> PrintConfiguration() ;
  analysis -> Initialize() ;

  if( ! analysis -> LoadData() ) return 0 ;  
  // This first steps deals with smearing effects, acceptance weights, fiducial cuts, etc. 
  // It also classifies events as signal or background
  analysis -> Analyse() ; 

  // SubstractBacgkround calculates the background probabilities from the identified background events
  // and stores the background events with the substracted probabilities
  // For the stored histograms, the background is substracted
  analysis -> SubtractBackground() ; 
  
  // Stores all the information in a TTree file
  // If requested, it also stores the requested histograms in an output root file
  analysis -> Finalise();

  delete analysis ;
  return 0 ; 
}
