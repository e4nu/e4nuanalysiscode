// _____________________________________________________________
/* This app is used to run a generic e4nu analysis            */
// _____________________________________________________________

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "TH1D.h"
#include "analysis/E4NuAnalysis.h"

using namespace std; 
using namespace e4nu;

int main( void ) {
  std::cout << "E4Nu analysis ongoing..." << std::endl;
  
  // This object can be initialized with a configuration file which contains information on the event run, 
  // cuts and analysis_1p1pim requirements, and output file location

  E4NuAnalysis * analysis = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/ConfFiles/e4nu_1p1pim_conf.txt") ;
  if( ! analysis ) return 0 ; 
  
  if( ! analysis -> LoadData() ) return 0 ;  

  analysis -> Analyse() ; 

  analysis -> SubtractBackground() ; 

  analysis -> Finalise();

  delete analysis ;

}
