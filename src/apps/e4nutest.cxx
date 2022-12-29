// _____________________________________________________________
/* This app is used to test aspects of the e4nu software only */
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
  std::cout << "Test ongoing..." << std::endl;
  E4NuAnalysis * analysis = new E4NuAnalysis("/genie/app/users/jtenavid/e4v/E4NuAnalysis/Source/vfork/ConfFiles/example_configuration.txt") ;
  if( ! analysis ) return 0 ; 
  
  if( ! analysis -> LoadData() ) return 0 ;  

  analysis -> Analyse() ; 
  analysis -> SubstractBackground() ; 
  analysis -> Finalise();

  delete analysis ;
  return 0 ; 
}
