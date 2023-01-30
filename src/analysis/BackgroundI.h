/**
 * \info These parameters are configurable 
 * the default values are set here
 **/

#ifndef _BACKGROUND_I_H_
#define _BACKGROUND_I_H_

#include <vector>
#include <map>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "analysis/ConfigureI.h"
#include "physics/EventI.h"
#include "physics/MCEvent.h"
#include "utils/Subtraction.h"

namespace e4nu { 

  class BackgroundI : public ConfigureI {

  public: 
    // Default constructor
    BackgroundI() ;
    BackgroundI( const std::string input_file ) ;
    BackgroundI( const double EBeam, const unsigned int TargetPdg ) ;

    bool SubstractBackground(void) ; 
      
  protected:
    bool InitializeFiducial(void) ;
    virtual ~BackgroundI();

    // Background definition
    //    std::map<int,std::vector<e4nu::EventI>> fBkg;
    std::map<int,std::vector<e4nu::MCEvent>> fBkg;

    Fiducial * fFiducialCut ;
    Subtraction * fRotation;

  };
}

#endif
