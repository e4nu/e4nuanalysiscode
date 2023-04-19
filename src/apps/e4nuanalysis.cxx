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

int main( void ) {
  std::cout << "E4Nu analysis ongoing..." << std::endl;

  // This object can be initialized with a configuration file which contains information on the event run, 
  // cuts and analysis requirements, and output file location
  char * env = std::getenv("E4NUANALYSIS") ; 
  std::string path( env ) ; 
  path += "/ConfFiles/" ;
  E4NuAnalysis * analysis = new E4NuAnalysis((path+"example_configuration.txt").c_str()) ;
  if( ! analysis ) return 0 ; 
  
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
