// _________________________________________________________________________________
/* This app is used to generate the radiated flux to be used for event generation */
// _________________________________________________________________________________
#include <iostream>
#include <vector>
#include "utils/RadiativeCorrUtils.h"
#include "plotting/PlottingUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"
#include "conf/RadConstants.h"
#include "TFile.h"
#include "TH1D.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf; 
using namespace e4nu::plotting;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
// --output-file : name and path of file where to store output //
// --ebeam : beam energy of your experiment                    //
// --target : target pdg                                       //
// --thickness : thickness of your experiment target           //
// --Nbins : number of bins for radiated flux histogram        //
// --Emin : minimum energy for your histogram axis             //
// --Emax : maximum energy for your histogram axis             //
//                                                             //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Generating radiated flux..." << std::endl;
  string output_file = "radiatedflux.root";
  double EBeam = 1 ; 
  int tgt = 1000060120 ;
  int nbins = 50 ; 
  double Emin = 0.75 ;
  double Emax = EBeam+0.02 ;
  double thickness = e4nu::conf::GetThickness(tgt); // Defaulted to CLAS6
  if( argc > 1 ) { // configure rest of analysis
    if( ExistArg("output-file",argc,argv)) {
      output_file = GetArg("output-file",argc,argv); 
    }
    if( ExistArg("ebeam",argc,argv)) {
      EBeam = stod(GetArg("ebeam",argc,argv)); 
      Emax = EBeam+0.02 ; 
    }
    if( ExistArg("target",argc,argv)) {
      tgt = stoi(GetArg("target",argc,argv)); 
      thickness = e4nu::conf::GetThickness(tgt); // Defaulted to CLAS6
    }
    if( ExistArg("thickness",argc,argv)) {
      thickness = stoi(GetArg("thickness",argc,argv)); 
    }
    if( ExistArg("Nbins",argc,argv)) {
      nbins = stoi(GetArg("Nbins",argc,argv)); 
    }
    if( ExistArg("Emin",argc,argv)) {
      Emin = stod(GetArg("Emin",argc,argv)); 
    }
    if( ExistArg("Emax",argc,argv)) {
      Emax = stod(GetArg("Emax",argc,argv)); 
    }
  }

  std::unique_ptr<TFile> myFile( TFile::Open(output_file.c_str(), "RECREATE") );
  TH1D * hradflux = new TH1D( "hradflux", "Radiated Flux", nbins, Emin, Emax) ;   

  TLorentzVector V4_beam(0,0,EBeam,EBeam);
  unsigned int nentries = 10000; 
  for( unsigned int i = 0 ; i < nentries ; ++i ) { 
    double egamma = SIMCEnergyLoss( EBeam, V4_beam, 11, tgt, thickness ) ;
    hradflux -> Fill( EBeam - egamma ) ; 
  }
  hradflux->Scale(1./hradflux->GetEntries());
  myFile->WriteObject(hradflux,"hradflux");

  std::cout << " Flux stored in " << output_file << std::endl;
  return 0 ;
}
