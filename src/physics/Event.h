/**
 * This class contains all the information of an event
 * It is an interface - therefore it is independent of the release type, which can be either MC or data
 * \date October 2022                                                                                                                                                                                              
 **/

#ifndef _EVENT_I_H_
#define _EVENT_I_H_

#include <iostream>
#include <map> 
#include "TLorentzVector.h"
#include "conf/ParticleI.h"

namespace e4nu {
  class Event {
  public : 
    Event(); 
    virtual ~Event();

    bool IsMC(void) { return fIsMC ;}
    unsigned int GetEventID(void) const { return fEventID ; } 
    int GetTargetPdg(void) const { return fTargetPdg ; }
    int GetInLeptPdg(void) const { return fInLeptPdg ; }
    int GetOutLeptPdg(void) const { return fOutLeptPdg ; }
    TLorentzVector GetInLepton4Mom(void) const { return fInLepton ; }
    TLorentzVector GetOutLepton4Mom(void) const { return fOutLepton ; }
    std::map<int,std::vector<TLorentzVector>> GetFinalParticles4Mom(void) const { return fFinalParticles ; }
    TLorentzVector GetInLeptonUnCorr4Mom(void) const { return fInLeptonUnCorr ; }
    TLorentzVector GetOutLeptonUnCorr4Mom(void) const { return fOutLeptonUnCorr ; }
    std::map<int,std::vector<TLorentzVector>> GetFinalParticlesUnCorr4Mom(void) const { return fFinalParticlesUnCorr ; }
    TLorentzVector GetVertex(void) const { return fVertex ; } 
    TVector3 GetRecoq3(void) const ; 

    void SetEventID( const unsigned int id ) { fEventID = id ; }
    void SetTargetPdg( int target_pdg ) ;
    void SetInLeptPdg( const int pdg ) { fInLeptPdg = pdg ; }
    void SetOutLeptPdg( const int pdg ) { fOutLeptPdg = pdg ; }
    void SetInLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 
    void SetOutLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetFinalParticle( const int pdg, const double E, const double px, const double py, const double pz ) ; 
    void SetOutUnCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetInUnCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 
    void SetFinalParticleUnCorr( const int pdg, const double E, const double px, const double py, const double pz ) ; 
    void SetOutLeptonKinematics( const TLorentzVector & tlvect ) { fOutLepton = tlvect ; }
    void SetInLeptonKinematics( const TLorentzVector & tlvect ) { fInLepton = tlvect ; }
    void SetFinalParticlesKinematics( const std::map<int,std::vector<TLorentzVector>> part_map ) { fFinalParticles = part_map ; }
    void SetOutLeptonUnCorrKinematics( const TLorentzVector & tlvect ) { fOutLeptonUnCorr = tlvect ; }
    void SetFinalParticlesUnCorrKinematics( const std::map<int,std::vector<TLorentzVector>> part_map ) { fFinalParticlesUnCorr = part_map ; }
    void SetVertex(const double vx, const double vy, const double vz, const double t) { fVertex.SetXYZT(vx, vy, vz, t) ; }

    void SetNProtons( const unsigned int n ) { fNP = n ; }
    void SetNNeutrons( const unsigned int n ) { fNN = n ; }
    void SetNPiP( const unsigned int n ) { fNPiP = n ; }
    void SetNPiM( const unsigned int n ) { fNPiM = n ; }
    void SetNPi0( const unsigned int n ) { fNPi0 = n ; }
    void SetNKP( const unsigned int n ) { fNKP = n ; }
    void SetNKM( const unsigned int n ) { fNKM = n ; }
    void SetNK0( const unsigned int n ) { fNK0 = n ; }
    void SetNEM( const unsigned int n ) { fNEM = n ; }
    void SetNOther( const unsigned int n ) { fNOther = n ; }

    // Number of particle getters:
    unsigned int GetRecoNProtons(void) { return fFinalParticles[conf::kPdgProton].size() ; }
    unsigned int GetRecoNNeutrons(void) { return fFinalParticles[conf::kPdgNeutron].size() ; }
    unsigned int GetRecoNPiP(void) { return fFinalParticles[conf::kPdgPiP].size() ; }
    unsigned int GetRecoNPiM(void) { return fFinalParticles[conf::kPdgPiM].size() ; }
    unsigned int GetRecoNPi0(void) { return fFinalParticles[conf::kPdgPi0].size() ; }
    unsigned int GetRecoNKP(void) { return fFinalParticles[conf::kPdgKP].size() ; }
    unsigned int GetRecoNKM(void) { return fFinalParticles[conf::kPdgKM].size() ; }
    unsigned int GetRecoNK0(void) { return fFinalParticles[conf::kPdgK0].size() ; }
    unsigned int GetRecoNEM(void) { return fFinalParticles[conf::kPdgPhoton].size() ; }
    
    // Analysed event properties
    double GetTotalWeight(void) const { return fWeight * fAccWght * fMottXSecWght ; }
    double GetEventWeight(void) const { return fWeight ; }
    void SetEventWeight( double wght ) { fWeight = wght ; }
    void SetMottXSecWeight(void) ; 
    double GetMottXSecWeight(void) const { return fMottXSecWght ; }
    double GetAccWght(void) const { return fAccWght ; }
    void SetAccWght( const double wght ) { fAccWght = wght ; }
    bool IsBkg(void) const{ return fIsBkg ; }
    void SetIsBkg( const bool bkg ) { fIsBkg = bkg ; }

    // Get analised event information
    double GetObservable( const std::string observable ) ;
    unsigned int GetEventMultiplicity( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) const ;
    unsigned int GetNSignalParticles( std::map<int,std::vector<TLorentzVector>> hadronic_system, const std::map<int,unsigned int> topology ) const ;
    int GetEventTotalVisibleCharge( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) const ;

    // Background debugging methods
    // This method returns the pdg of the visible particles before and after fiducial cuts
    std::map<unsigned int,std::pair<std::vector<int>,double>> GetAnalysisRecord(void) const { return fAnalysisRecord; }
    void StoreAnalysisRecord( unsigned int analysis_step ) ; 

    // GENIE Specific variables
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
    int GetRESID(void) const { return fresid ; }

    unsigned int GetTrueNProtons(void) const { return fNP ; }
    unsigned int GetTrueNNeutrons(void) const { return fNN ; }
    unsigned int GetTrueNPiP(void) const { return fNPiP ; }
    unsigned int GetTrueNPiM(void) const { return fNPiM ; }
    unsigned int GetTrueNPi0(void) const { return fNPi0 ; }
    unsigned int GetTrueNKP(void) const { return fNKP ; }
    unsigned int GetTrueNKM(void) const { return fNKM ; }
    unsigned int GetTrueNK0(void) const { return fNK0 ; }
    unsigned int GetTrueNEM(void) const { return fNEM ; }
    unsigned int GetTrueNOther(void) const { return fNOther ; } 

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
    void SetRESID( const int resid ) { fresid = resid ; } 
    
    // Radiative correction utils
    TLorentzVector GetInCorrLepton4Mom(void) const { return fInCorrLepton ; }
    TLorentzVector GetOutCorrLepton4Mom(void) const { return fOutCorrLepton ; }
    void SetInCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 
    void SetOutCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;

  protected : 
    // Common funtionalities which depend on MC or data 
    bool fIsMC ;
    TLorentzVector fInLepton ; 
    TLorentzVector fInCorrLepton ;  // Only used in radiative corrections app 
    TLorentzVector fOutLepton ; 
    TLorentzVector fOutCorrLepton ; // Only used in radiative corrections app 
    std::map<int,std::vector<TLorentzVector>> fFinalParticles ; 

    // Store uncorrected kinematics
    TLorentzVector fInLeptonUnCorr ; 
    TLorentzVector fOutLeptonUnCorr ; 
    std::map<int,std::vector<TLorentzVector>> fFinalParticlesUnCorr ; 

    unsigned int fNP, fNN, fNPiP, fNPiM, fNPi0, fNKP, fNKM, fNK0, fNEM, fNOther ; 

    double fWeight = 1. ; 
    double fAccWght = 1. ;
    double fMottXSecWght = 1. ;

  private :

    unsigned int fEventID ; 
    int fTargetPdg ; 
    int fInLeptPdg ; 
    int fOutLeptPdg ; 
    TLorentzVector fVertex ; 

    // GENIE Specific Variables
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
    int fresid ; 

    // Background ID
    bool fIsBkg = false ; 
    
    // Analysis Record, used to store analized events as a function of multiplicity
    std::map<unsigned int,std::pair<std::vector<int>,double>> fAnalysisRecord; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
