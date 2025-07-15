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
#include "utils/Fiducial.h"
#include "analysis/BackgroundI.h"

namespace e4nu {

  class AnalysisI : public BackgroundI {

  public:
    // Default constructor
    AnalysisI() ;
    AnalysisI( const std::string input_file ) ;
    AnalysisI( const double EBeam, const unsigned int TargetPdg ) ;

    bool Analyse( Event & event ) ;
    void Initialize(void) ;
    bool Finalise(void) ;

  protected :
    void CookEvent( Event & event ) ;
    void PlotBkgInformation( const Event event ) ;
    void ApplyMomentumCut( Event & event ) ;
    bool ApplyFiducialCutExtra( Event & event ) ;
    bool ApplyFiducialCut( Event & event, bool apply_fiducial ) ;
    double GetElectronMinTheta( TLorentzVector emom ) ;
    bool StoreTree(Event event);

    // ID for Background historams
    unsigned int kid_signal, kid_tottruebkg, kid_totestbkg, kid_acccorr;
    unsigned int kid_2p0pitruebkg, kid_1p1pitruebkg, kid_2p1pitruebkg, kid_1p2pitruebkg ;
    unsigned int kid_2p0piestbkg, kid_1p1piestbkg, kid_2p1piestbkg, kid_1p2piestbkg ;

    virtual ~AnalysisI();

  private :
    TF1 * kElectronFit = nullptr ;

  };
}

#endif
