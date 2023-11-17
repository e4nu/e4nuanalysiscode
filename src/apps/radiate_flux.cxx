// _________________________________________________________________________________
/* This app is used to generate the radiated flux to be used for event generation */
// _________________________________________________________________________________
#include <iostream>
#include <vector>
#include "utils/RadiativeCorrUtils.h"
#include "plotting/PlottingUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"
#include "TFile.h"
#include "TH1D.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf; 
using namespace e4nu::plotting;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Generating radiated flux..." << std::endl;
  string output_file = "radiatedflux.root";
  double EBeam = 1 ; 
  int tgt = 1000060120 ;
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("output-file",argc,argv)) {
      output_file = GetArg("output-file",argc,argv); 
    }
    if( ExistArg("ebeam",argc,argv)) {
      EBeam = stod(GetArg("ebeam",argc,argv)); 
    }
    if( ExistArg("target",argc,argv)) {
      tgt = stoi(GetArg("target",argc,argv)); 
    }
  }

  std::unique_ptr<TFile> myFile( TFile::Open(output_file.c_str(), "RECREATE") );
  TH1D * hradflux = new TH1D( "hradflux", "Radiated Flux", 50, 0.75, EBeam+0.02) ;   

  TLorentzVector V4_beam(0,0,EBeam,EBeam);
  unsigned int nentries = 10000; 
  for( unsigned int i = 0 ; i < nentries ; ++i ) { 
    double egamma =  SIMCEnergyLoss( EBeam, V4_beam, 11, tgt ) ;
    hradflux -> Fill( EBeam - egamma ) ; 
  }
  hradflux->Scale(1./hradflux->GetEntries());
  myFile->WriteObject(hradflux,"hradflux");
  return 0 ;
}
