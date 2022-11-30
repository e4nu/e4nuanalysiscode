/**
 * This class contains all the information of an event
 * It is an interface - therefore it is independent of the release type, which can be either MC or data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _EVENT_MC_H_
#define _EVENT_MC_H_

#include <iostream>
#include "physics/EventI.h"

namespace e4nu {
  class MCEvent : EventI {
  public : 
    MCEvent(); 
    virtual ~MCEvent();

    bool IsQEL(void) const { return fIsQEL; }
    bool IsRES(void) const { return fIsRES; } 
    bool IsMEC(void) const { return fIsMEC; }

  protected : 
    void SetIsQEL( bool qel ) { fIsQEL = qel ; } 
    void SetIsRES( bool res ) { fIsRES = res ; } 
    void SetIsMEC( bool mec ) { fIsMEC = mec ; }
    void SetISDIS( bool dis ) { fIsDIS = dis ; } 

    // Flip phi with respect to GENIE 
    // GENIE Coordinate system is flipped with respect to class
    void SetOutLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetInLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 

    // Common funtionalities which depend on MC or data = definition 

  private :
    bool fIsQEL ;
    bool fIsRES ; 
    bool fIsMEC ; 
    bool fIsDIS ; 
    
  };
}

#endif
