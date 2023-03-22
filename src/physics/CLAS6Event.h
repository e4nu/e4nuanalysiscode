/**
 * This class contains all the information of an event
 * It is an interface - therefore it is independent of the release type, which can be either CLAS6 or data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _EVENT_CLAS6_H_
#define _EVENT_CLAS6_H_

#include <iostream>
#include "physics/EventI.h"
#include "physics/CLAS6EventHolder.h"

namespace e4nu {
  class CLAS6Event : public EventI {
  public : 
    CLAS6Event(); 
    virtual ~CLAS6Event();

    TLorentzVector GetVertex(void) const { return fVertex ; }

    friend class CLAS6EventHolder ; 

  protected : 
 
    void SetVertex(const double vx, const double vy, const double vz, const double t) { fVertex.SetXYZT(vx, vy, vz, t) ; }

    // Flip phi with respect to GENIE 
    // GENIE Coordinate system is flipped with respect to class
    void SetOutLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetInLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 


  private :

    TLorentzVector fVertex ; 
  };
}

#endif
