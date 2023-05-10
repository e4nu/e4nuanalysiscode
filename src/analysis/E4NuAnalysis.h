/**
 * This class  is the interface for any analysis
 * Summary: This class deals with a generic E4nu analysis. 
 * It will treat the data according to it's configuration as MC or CLAS6 data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _E4NUANALYSIS_H_
#define _E4NUANALYSIS_H_

#include <iostream>
#include "TF1.h"
// Include here new analysis classes
#include "analysis/MCCLAS6StandardAnalysis.h"
#include "analysis/CLAS6StandardAnalysis.h"

using namespace e4nu::conf ; 

namespace e4nu {
  class E4NuAnalysis : 
    // Include here new analysis classes 
    public MCCLAS6StandardAnalysis, 
    public CLAS6StandardAnalysis { 
  public : 
    E4NuAnalysis(); 
    E4NuAnalysis( const std::string conf_file ) ;
    E4NuAnalysis( const double EBeam, const unsigned int TargetPdg ) ;

    bool LoadData(void) ; 
    bool Analyse(void) ; 
    void ClassifyEvent( Event event ) ;
    bool SubtractBackground( void ) ;
    bool Finalise(void);

    virtual ~E4NuAnalysis();

  private : 

    e4nu::Event * GetValidEvent( const unsigned int event_id ) ;
    unsigned int GetNEvents( void ) const ;

    // Event Holder for signal and background
    std::map<int,std::vector<e4nu::Event>> kAnalysedEventHolder;

    void Initialize(void) ; 
    
  };
}

#endif
