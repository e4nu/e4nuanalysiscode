/**
 * This class  is the interface for any analysis
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _ANALYSIS_I_H_
#define _ANALYSIS_I_H_

#include <iostream>
#include "conf/ConfigureI.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class AnalysisI : public ConfigureI {
  public : 

    // Main Analyse function
    bool Analyse(void) const ; 

    // Store events on root file for further analysis
    bool StoreAnalysis( const std::string out_file ) { ; }

    // Create histograms and store in file 
    bool Finalize( const std::string out_file ) { ; }

    virtual ~AnalysisI();

  protected : 
    AnalysisI(); 
    AnalysisI( const std::string conf_file ) ;
    AnalysisI( const double EBeam, const unsigned int TargetPdg ) ; 

    // Load Data from root file:
    virtual bool LoadData( const std::string file ) = 0 ; 
    virtual bool LoadData( const std::string file, const unsigned int nmax ) = 0 ; 
    virtual void GetEvent( const unsigned int event_id ) const = 0 ;

    bool fIsDataLoaded ;

  private :
    void Initialize(void) ;

  };
}

#endif
