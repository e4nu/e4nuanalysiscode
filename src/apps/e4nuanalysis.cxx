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

  // Compute acceptance correction files
  bool compute_trueacc = analysis -> ComputeTrueAccCorrection() ; 
  bool compute_truerecoacc = analysis -> ComputeTrueRecoAccCorrection() ; 

  if ( analysis -> IsData() ) {   
      compute_trueacc = false ; 
      compute_truerecoacc = false ; 
      analysis -> SetApplyFiducial( true ) ; 
      analysis -> SetApplyAccWeights( true ) ; 
      analysis -> SetApplyReso( true ) ;  
      std::string OutputFile_data = analysis->GetOutputFile() + "_clas6data" ;
      analysis -> SetOutputFile( OutputFile_data ) ; 
  }

  if( compute_trueacc && compute_truerecoacc ) { 
    std::cout << " Due to ROOT related issues, the code can only compute one at the time... Abort." <<std::endl;
    return 0 ; 
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
