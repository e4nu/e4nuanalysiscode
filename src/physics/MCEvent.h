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
  class MCEvent : public EventI {
  public : 
    MCEvent(); 
    virtual ~MCEvent();

    bool IsEM(void) const { return fIsEM; }
    bool IsCC(void) const { return fIsCC; }
    bool IsNC(void) const { return fIsNC; }
    bool IsQEL(void) const { return fIsQEL; }
    bool IsRES(void) const { return fIsRES; } 
    bool IsMEC(void) const { return fIsMEC; }
    bool IsDIS(void) const { return fIsDIS; }
    
    double GetTrueQ2s(void) const { return fTrueQ2s ; }
    double GetTrueWs(void) const { return fTrueWs; }
    double GetTruexs(void) const { return fTruexs ; }
    double GetTrueys(void) const { return fTrueys ; }
    double GetTrueQ2(void) const { return fTrueQ2 ; }
    double GetTrueW(void) const { return fTrueW ; }
    double GetTruex(void) const { return fTruex ; }
    double GetTruey(void) const { return fTruey ; }

    TLorentzVector GetVertex(void) const { return fVertex ; }
    
    bool GetAccWght(void) const { return fAccWght ; }
    void SetAccWght( bool wght ) { fAccWght = wght ; }
    friend class MCEventHolder ; 

  protected : 
    void SetIsEM( const bool em ) { fIsEM = em ; }
    void SetIsCC( const bool cc ) { fIsCC = cc ; }
    void SetIsNC( const bool nc ) { fIsNC = nc ; }
    void SetIsQEL( const bool qel ) { fIsQEL = qel ; } 
    void SetIsRES( const bool res ) { fIsRES = res ; } 
    void SetIsMEC( const bool mec ) { fIsMEC = mec ; }
    void SetIsDIS( const bool dis ) { fIsDIS = dis ; } 

    void SetTrueQ2s( const double Q2s ) { fTrueQ2s = Q2s ; }
    void SetTrueWs( const double Ws ) { fTrueWs = Ws ; }
    void SetTruexs( const double xs ) { fTruexs = xs ; }
    void SetTrueys( const double ys ) { fTrueys = ys ; }
    void SetTrueQ2( const double Q2 ) { fTrueQ2 = Q2 ; }
    void SetTrueW( const double W ) { fTrueW = W ; }
    void SetTruex( const double x ) { fTruex = x ; }
    void SetTruey( const double y ) { fTruey = y ; } 

    void SetVertex(const double vx, const double vy, const double vz, const double t) { fVertex.SetXYZT(vx, vy, vz, t) ; }

    // Flip phi with respect to GENIE 
    // GENIE Coordinate system is flipped with respect to class
    void SetOutLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetInLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 
    void SetFinalParticle( const int pdg, const double E, const double px, const double py, const double pz ); 

  private :
    bool fIsEM ; 
    bool fIsCC ; 
    bool fIsNC ; 
    bool fIsQEL ;
    bool fIsRES ; 
    bool fIsMEC ; 
    bool fIsDIS ;
    
    double fTrueQ2s ; 
    double fTrueWs ; 
    double fTruexs ; 
    double fTrueys ; 
    double fTrueQ2 ; 
    double fTrueW ; 
    double fTruex ; 
    double fTruey ; 

    bool fAccWght = 1 ;
    TLorentzVector fVertex ; 

  };
}

#endif
