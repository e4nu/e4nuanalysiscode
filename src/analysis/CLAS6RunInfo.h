/**
 * This class reads the run information
 * \date June 2025                                                                                                                                                                                              
 **/

#ifndef _CLAS6RUNINFO_H_
#define _CLAS6RUNINFO_H_

#include <iostream>
#include <map>
#include <vector>
#include <string>

using namespace std ;

namespace e4nu {
  class CLAS6RunInfo {
  public : 
    CLAS6RunInfo(); 

  protected : 
    virtual ~CLAS6RunInfo();

    void LoadData(void);
    void RetrieveAllRuns( const double beam, const string target ) ; 
    void NeglectRuns( const std::vector<unsigned int> run_ids ) ;
    double GetTotalIntegratedCharge(void) ; 
    bool IsRunValid( const unsigned int run_id ) ;
    void SetIsRunInFile( const unsigned int run_id ) ;
    bool IsEnabledRun( const unsigned int run_id) { return kEnabledRun[run_id] ; } 
    double GetRunIntegratedCharge( const unsigned int run_id ) { return kLiveC[run_id]*kFcupA[run_id]; }
    double GetRunBeam( const unsigned int run_id ) {return kBeam[run_id]; }
    string GetRunTarget( const unsigned int run_id ) {return kTarget[run_id]; }
    string GetRunDate( const unsigned int run_id ) {return kDate[run_id]; }
    string GetRunTime( const unsigned int run_id ) {return kTime[run_id]; }
    unsigned int GetRunNFiles(const unsigned int run_id ) {return kNFiles[run_id]; }
    unsigned int GetRunEvents(const unsigned int run_id ) {return kNRunEvents[run_id]; }
    double GetRunI( const unsigned int run_id ) {return kI[run_id] ; }
    int GetRunTorus( const unsigned int run_id ) {return kTorus[run_id]; } 
    unsigned int GetRunECHI( const unsigned int run_id ) {return kechi[run_id]; }  
    string GetRunComments( const unsigned int run_id ) {return kComments[run_id]; } 

  private :
    unsigned int kRunID = 0;
    std::map<unsigned int,bool> kEnabledRun ; 
    std::map<unsigned int,bool> kIsRunInFile ; 
    std::map<unsigned int,string> kDate ;
    std::map<unsigned int,string> kTime ;
    std::map<unsigned int,unsigned int> kNFiles ;
    std::map<unsigned int,unsigned int> kNRunEvents ; 
    std::map<unsigned int,string> kTarget ;
    std::map<unsigned int,double> kLiveC ;
    std::map<unsigned int,double> kBeam ;
    std::map<unsigned int,double> kI;
    std::map<unsigned int,int> kTorus;
    std::map<unsigned int,int> kMini;
    std::map<unsigned int,unsigned int> kechi ;
    std::map<unsigned int,double> kFcupA ;
    std::map<unsigned int,string> kComments;

  };
}

#endif
