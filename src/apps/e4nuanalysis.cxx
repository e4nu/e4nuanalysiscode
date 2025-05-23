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

/////////////////////////////////////////////////////////////////////////
// Options:                                                            //
/////////////////////////////////////////////////////////////////////////
// --conf-file) Configuration file location                            //
// --root-file) Input MC File (GENIE gst format or CLAS6 data)         //
// --output-file) Output file which contains event information         //
// --analysis-type) Type of analysis:                                  //
//    - ComputeTrueAccCorr (output file named as "_true.root")         //
//    - ComputeTrueRecoAccCorr(output file named as "_truereco .root") //
//    - IsData                                                         //
//    - ClosureTest                                                    //
// --rad-corr bool ; used to change the output name of file            //
//                   when radiative corr are used in MC                //
// --xsec-file) XSecFile (only for MC)                                 //
// --bkg-mult) Maximum multiplicity used in bkg subtraction method     //
// --phi-shift) Shift on phy applied to reduce fiducial volume         //
/////////////////////////////////////////////////////////////////////////

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
  bool compute_closuretest = analysis -> GetDebugBkg() ;
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
      else if ( GetArg("analysis-type",argc,argv) == "ClosureTest" ) {
        std::cout << " Implementing closure test = True " <<std::endl;
        compute_trueacc = false ; compute_truerecoacc = false ; compute_closuretest = true ;
      }
      else if ( GetArg("analysis-type",argc,argv) == "IsData" ) { is_data = true ; compute_trueacc = false ; compute_truerecoacc = false ; }
    }
    if( ExistArg("rad-corr",argc,argv) && !is_data ) {
      analysis -> SetRadCorr( true ) ;
    }
    if( ExistArg("xsec-file",argc,argv) && !is_data ) {
      if ( !analysis -> SetXSecFile( GetArg("xsec-file",argc,argv) ) ) return false;
    }
    if( ExistArg("bkg-mult",argc,argv) && is_data ) {
      analysis -> SetMaxBkgMult( atoi(GetArg("bkg-mult",argc,argv).c_str()) ) ;
    }
    if( ExistArg("phi-shift",argc,argv) ) {
      analysis -> SetFidAngleShift( stod(GetArg("phi-shift",argc,argv).c_str()) ) ;
    }

    if( ExistArg("output-file",argc,argv)) {
      std::string final_name = GetArg("output-file",argc,argv) ;
      if( is_data ) {
        unsigned int max_mult = analysis->GetMaxBkgMult() ;
        final_name += "_"+std::to_string(max_mult)+"MaxBkgMult";
      }
      if( ExistArg("phi-shift",argc,argv) ) {
        double shift = stod(GetArg("phi-shift",argc,argv).c_str()) ;
        if ( shift != 0 ) {
          final_name += "_Shift_"+std::to_string(shift)+"deg";
        }
      }

      analysis -> SetOutputFile( final_name );
    }
  }

  if ( is_data ) {
    compute_trueacc = false ;
    compute_truerecoacc = false ;
    // We need to apply the fiducial in the data as well for technical reasons
    // & also when we want to compute the systematics associated to the geometrical shift
    // For these reasons, we leave it in. It should not affect at all the results when we do not shift the fiducial
    // As these are already accounted for when we process the data
    analysis -> SetApplyFiducial( true ) ;
    analysis -> SetApplyAccWeights( false ) ;
    analysis -> SetApplyReso( false ) ;
    analysis -> SetSubtractBkg( true ) ;
    std::string OutputFile_data = analysis->GetOutputFile() + "_clas6data" ;
    analysis -> SetOutputFile( OutputFile_data ) ;
  }

  if ( compute_trueacc ) {
    analysis -> SetTrueSignal( true ) ;
    analysis -> SetApplyFiducial( false ) ;
    analysis -> SetApplyAccWeights( false ) ;
    // Correcting for smearing gives unphysical ECal peak. We opt to not correct for this effect
    // MC Generators will have to smear their results to compare against data
    analysis -> SetApplyReso( true ) ;
    // !!!!!!!
    //analysis -> SetUseAllSectors( true ) ;
    //analysis -> EnableAllSectors( true ) ;
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
  } else if ( compute_closuretest ) {
    analysis -> SetTrueSignal( false ) ;
    analysis -> SetApplyFiducial( true ) ;
    analysis -> SetApplyAccWeights( true ) ;
    analysis -> SetApplyReso( true ) ;
    analysis -> SetSubtractBkg( true ) ;
    analysis -> SetOutputFile( analysis->GetOutputFile() + "_closuretest"  ) ;
    analysis -> SetDebugBkg(true) ;
    std::cout << " Computing Closure test..."<<std::endl;
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
  // The format of the output root file is set in the MCCLAS6AnalysisI or CLAS6AnalysisI StoreTree function
  analysis -> Finalise();

  delete analysis ;
  return 0 ;
}
