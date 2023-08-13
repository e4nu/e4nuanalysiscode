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

using namespace std; 
using namespace e4nu;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////
// 1) Configuration file location                              //   
// 2) InputFile                                                // 
// 3) OutputFile                                               //
// 4) Type of analysis:                                        // 
//    - ComputeTrueAccCorr                                     // 
//    - ComputeTrueRecoAccCorr                                 //
//    - IsData                                                 //   
// 5) XSecFile (only for MC)                                   //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "E4Nu analysis ongoing..." << std::endl;

  // This object can be initialized with a configuration file which contains information on the event run, 
  // cuts and analysis requirements, and output file location
  char * env = std::getenv("E4NUANALYSIS") ; 
  std::string path( env ) ; 
  path += "/" ;
  std::string config_file = "ConfFiles/example_configuration.txt" ;
  if ( argv[1] ) config_file = argv[1] ;

  E4NuAnalysis * analysis = new E4NuAnalysis((path+config_file).c_str()) ;
  if( ! analysis ) return 0 ; 

  if( argv[2] ) { 
    analysis -> SetInputFile( argv[2] ) ;
  } else { 
    std::cout << " ERR: Input file not specified. Abort... " << std::endl;
    return 0 ; 
  }

  if( argv[3] ) { 
    analysis ->SetOutputFile( argv[3] ) ; 
  } else { 
    std::cout << " ERR: Output file not specified. Abort... " << std::endl;
    return 0 ; 
  }

  // Compute acceptance correction files
  bool compute_trueacc = analysis -> ComputeTrueAccCorrection() ; 
  bool compute_truerecoacc = analysis -> ComputeTrueRecoAccCorrection() ; 
  bool is_data = analysis -> IsData() ; 

  if( argv[4] ) {
    
    if ( strcmp( argv[4], "ComputeTrueAccCorr" ) == 0 ) { 
      std::cout << " ComputeTrueAccCorr = True " <<std::endl;
      compute_trueacc = true ; compute_truerecoacc = false ; 
    }
    else if ( strcmp( argv[4], "ComputeTrueRecoAccCorr" ) == 0 ) { 
      std::cout << " ComputeRecoAccCorr = True " <<std::endl;
      compute_trueacc = false ; compute_truerecoacc = true ; 
    }
    else if ( strcmp( argv[4], "IsData" ) == 0 ) { is_data = true ; compute_trueacc = false ; compute_truerecoacc = false ; }
    else { 
      std::cout << " The second argument defines the type of run: (1) ComputeTrueAccCorr: Analyse using true information,\n" 
		<< "  (2) ComputeTrueRecoAccCorr: Analyse using True events + detector effects,\n"
		<< "  (3) IsData: Analyse data. \n Set one of the correct options\n"
		<< " You can also simply disable this option. In this case, the default from the configuration file will be used.\n";
      return 0 ; 
    }
  }

  if( argv[5] && !is_data ) { 
    analysis -> SetXSecFile( argv[5] ) ; 
  } else { 
    std::cout << " Input file not specified. Abort... " << std::endl;
    return 0 ; 
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
    analysis -> SetApplyReso( false ) ;
    analysis -> SetUseAllSectors( true ) ;
    analysis -> EnableAllSectors( true ) ;
    std::string OutputFile_true = analysis->GetOutputFile() + "_true" ;
    analysis -> SetOutputFile( OutputFile_true ) ;
    std::cout << " Computing true analysis distributions for acceptance correction..."<<std::endl;
  } else if ( compute_truerecoacc ) { 
    analysis -> SetTrueSignal( true ) ; 
    analysis -> SetApplyFiducial( true ) ; 
    analysis -> SetApplyAccWeights( true ) ; 
    analysis -> SetApplyReso( true ) ; 
    analysis -> SetSubtractBkg( false ) ; 
    std::string OutputFile_reco = analysis->GetOutputFile() + "_truereco" ;
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
