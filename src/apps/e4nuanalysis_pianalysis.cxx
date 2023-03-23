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
  E4NuAnalysis * analysis_1p1pim = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/ConfFiles/e4nu_1p1pim_conf.txt") ;
  if( ! analysis_1p1pim ) return 0 ; 
  
  if( ! analysis_1p1pim -> LoadData() ) return 0 ;  

  analysis_1p1pim -> Analyse() ; 

  analysis_1p1pim -> SubtractBackground() ; 

  analysis_1p1pim -> Finalise();

  delete analysis_1p1pim ;

  E4NuAnalysis * analysis_1p1pip = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/ConfFiles/e4nu_1p1pip_conf.txt") ;
  if( ! analysis_1p1pip ) return 0 ; 
  
  if( ! analysis_1p1pip -> LoadData() ) return 0 ;  

  analysis_1p1pip -> Analyse() ; 

  analysis_1p1pip -> SubtractBackground() ; 

  analysis_1p1pip -> Finalise();

  delete analysis_1p1pip ;

  E4NuAnalysis * analysis_Inclusive1pim = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/ConfFiles/e4nu_1pimInclusive_conf.txt") ;
  if( ! analysis_Inclusive1pim ) return 0 ; 
  
  if( ! analysis_Inclusive1pim -> LoadData() ) return 0 ;  

  analysis_Inclusive1pim -> Analyse() ; 

  analysis_Inclusive1pim -> SubtractBackground() ; 

  analysis_Inclusive1pim -> Finalise();

  delete analysis_Inclusive1pim ;

  E4NuAnalysis * analysis_Inclusive1pip = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/ConfFiles/e4nu_1pipInclusive_conf.txt") ;
  if( ! analysis_Inclusive1pip ) return 0 ; 
  
  if( ! analysis_Inclusive1pip -> LoadData() ) return 0 ;  

  analysis_Inclusive1pip -> Analyse() ; 

  analysis_Inclusive1pip -> SubtractBackground() ; 

  analysis_Inclusive1pip -> Finalise();

  delete analysis_Inclusive1pip ;

  return 0 ; 
}
