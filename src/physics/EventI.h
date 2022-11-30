/**
 * This class contains all the information of an event
 * It is an interface - therefore it is independent of the release type, which can be either MC or data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _EVENT_I_H_
#define _EVENT_I_H_

#include <iostream>
#include "TLorentzVector.h"

namespace e4nu {
  class EventI {
  public : 
    virtual ~EventI();

    bool IsMC(void) { return fIsMC ;}
    const TLorentzVector GetInLepton4Mom(void) const { return fInLepton ; }
    const TLorentzVector GetOutLepton4Mom(void) const { return fOutLepton ; }
    const int GetTargetPdg(void) const { return ftarget_pdg ; }
    const int GetHitNuclPdg(void) const { return fhit_nucl ; }

  protected : 
    EventI(); 

    // Common Functionalities    
    void SetTargetPDG( const int target_pdg ) { ftarget_pdg = target_pdg ; } 
    void SetHitNucleon( const int hit_nucl ) { fhit_nucl = hit_nucl ; }
    
    void SetOutLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetInLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 

    void SetOutLeptonKinematics( TLorentzVector tlvect ) { fOutLepton = tlvect ; }
    void SetInLeptonKinematics( TLorentzVector tlvect ) { fInLepton = tlvect ; }
    
    // Common funtionalities which depend on MC or data

    bool fIsMC ;
  private :
   
    int ftarget_pdg ; 
    int fhit_nucl ; 

    TLorentzVector fInLepton ; 
    TLorentzVector fOutLepton ; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
