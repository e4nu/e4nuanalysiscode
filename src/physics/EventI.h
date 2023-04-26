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
  class EventI {
  public : 
    EventI(); 
    virtual ~EventI();

    bool IsMC(void) { return fIsMC ;}
    unsigned int GetEventID(void) const { return fEventID ; } 
    TLorentzVector GetInLepton4Mom(void) const { return fInLepton ; }
    TLorentzVector GetOutLepton4Mom(void) const { return fOutLepton ; }
    std::map<int,std::vector<TLorentzVector>> GetFinalParticles4Mom(void) const { return fFinalParticles ; }
    TLorentzVector GetInLeptonUnCorr4Mom(void) const { return fInLeptonUnCorr ; }
    TLorentzVector GetOutLeptonUnCorr4Mom(void) const { return fOutLeptonUnCorr ; }
    std::map<int,std::vector<TLorentzVector>> GetFinalParticlesUnCorr4Mom(void) const { return fFinalParticlesUnCorr ; }
 
    int GetTargetPdg(void) const { return fTargetPdg ; }
    int GetInLeptPdg(void) const { return fInLeptPdg ; }
    int GetOutLeptPdg(void) const { return fOutLeptPdg ; }

    unsigned int GetRecoNProtons(void) { return fFinalParticles[conf::kPdgProton].size() ; }
    unsigned int GetRecoNNeutrons(void) { return fFinalParticles[conf::kPdgNeutron].size() ; }
    unsigned int GetRecoNPiP(void) { return fFinalParticles[conf::kPdgPiP].size() ; }
    unsigned int GetRecoNPiM(void) { return fFinalParticles[conf::kPdgPiM].size() ; }
    unsigned int GetRecoNPi0(void) { return fFinalParticles[conf::kPdgPi0].size() ; }
    unsigned int GetRecoNKP(void) { return fFinalParticles[conf::kPdgKP].size() ; }
    unsigned int GetRecoNKM(void) { return fFinalParticles[conf::kPdgKM].size() ; }
    unsigned int GetRecoNK0(void) { return fFinalParticles[conf::kPdgK0].size() ; }
    unsigned int GetRecoNEM(void) { return fFinalParticles[conf::kPdgPhoton].size() ; }

    double GetTotalWeight(void) const { return fWeight * fAccWght * fMottXSecWght ; }
    double GetEventWeight(void) const { return fWeight ; }
    void SetEventWeight( double wght ) { fWeight = wght ; }
    void SetMottXSecWeight(void) ; 
    double GetMottXSecWeight(void) const { return fMottXSecWght ; }
    double GetAccWght(void) const { return fAccWght ; }
    bool IsBkg(void) const{ return fIsBkg ; }
    void SetIsBkg( const bool bkg ) { fIsBkg = bkg ; }

    void SetOutLeptonKinematics( const TLorentzVector & tlvect ) { fOutLepton = tlvect ; }
    void SetInLeptonKinematics( const TLorentzVector & tlvect ) { fInLepton = tlvect ; }
    void SetFinalParticlesKinematics( const std::map<int,std::vector<TLorentzVector>> part_map ) { fFinalParticles = part_map ; }

    void SetOutLeptonUnCorrKinematics( const TLorentzVector & tlvect ) { fOutLeptonUnCorr = tlvect ; }
    void SetFinalParticlesUnCorrKinematics( const std::map<int,std::vector<TLorentzVector>> part_map ) { fFinalParticlesUnCorr = part_map ; }
    
    double GetObservable( const std::string observable ) ;
    unsigned int GetEventMultiplicity( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) ;
    unsigned int GetNSignalParticles( std::map<int,std::vector<TLorentzVector>> hadronic_system, const std::map<int,unsigned int> topology ) ;
    int GetEventTotalVisibleCharge( const std::map<int,std::vector<TLorentzVector>> hadronic_system ) ;
    TVector3 GetRecoq3(void) const ; 

    // Background debugging methods
    // This method returns the pdg of the visible particles before and after fiducial cuts
    std::map<unsigned int,std::pair<std::vector<int>,double>> GetAnalysisRecord(void) { return fAnalysisRecord; }
    void StoreAnalysisRecord( unsigned int analysis_step ) ; 

  protected : 

    // Common Functionalities    
    void SetEventID( const unsigned int id ) { fEventID = id ; }
    void SetTargetPdg( const int target_pdg ) { fTargetPdg = target_pdg ; } 
    void SetInLeptPdg( const int pdg ) { fInLeptPdg = pdg ; }
    void SetOutLeptPdg( const int pdg ) { fOutLeptPdg = pdg ; }
    
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

    void SetOutLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetInLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 
    void SetFinalParticle( const int pdg, const double E, const double px, const double py, const double pz ) ; 

    void SetOutUnCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) ;
    void SetInUnCorrLeptonKinematics( const double energy, const double px, const double py, const double pz ) ; 
    void SetFinalParticleUnCorr( const int pdg, const double E, const double px, const double py, const double pz ) ; 
    
    // Common funtionalities which depend on MC or data 
    bool fIsMC ;
    TLorentzVector fInLepton ; 
    TLorentzVector fOutLepton ; 
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

    bool fIsBkg = false ; 
    
    std::map<unsigned int,std::pair<std::vector<int>,double>> fAnalysisRecord; 

    void Initialize(void) ;
    void Clear(void); 

  };
}

#endif
