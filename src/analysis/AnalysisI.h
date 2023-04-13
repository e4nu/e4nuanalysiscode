/**
 * This class deals with common analysis features which are present in data and MC
 **/

#ifndef _ANALYSIS_I_H_
#define _ANALYSIS_I_H_

#include <vector>
#include <map>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "physics/EventI.h"
#include "utils/Fiducial.h"
#include "analysis/BackgroundI.h"

namespace e4nu { 

  class AnalysisI : public BackgroundI {

  public: 
    // Default constructor
    AnalysisI() ;
    AnalysisI( const std::string input_file ) ;
    AnalysisI( const double EBeam, const unsigned int TargetPdg ) ;

    bool Analyse( EventI * event ) ; 
    bool Finalise(void) const ; 

  protected : 
    void CookEvent( EventI * event ) ;
    void PlotBkgInformation( EventI * event ) ;

    // ID for Background historams
    unsigned int kid_signal, kid_tottruebkg, kid_totestbkg, kid_acccorr;
    unsigned int kid_2p0pitruebkg, kid_1p1pitruebkg, kid_2p1pitruebkg, kid_1p2pitruebkg ;
    unsigned int kid_2p0piestbkg, kid_1p1piestbkg, kid_2p1piestbkg, kid_1p2piestbkg ;

    virtual ~AnalysisI();
      
  };
}

#endif
