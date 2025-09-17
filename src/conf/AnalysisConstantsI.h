/**                                                                                                                                                                                                                
 * This file contains general constants to be used in any analysis
 **/

#ifndef _ANALYSISCONST_I_H_
#define _ANALYSISCONST_I_H_

namespace e4nu {
  namespace conf {
    // Phase space limits
    static const double kQ2Max = 1.75 ; 
    static const double kQ2Min = 0.4 ; 
    static const double kMinTheta = 10 ;
    static const double kMaxTheta = 60 ;
    static const double kMinEePrime = 0.5 ; 
    static const double kMaxEePrime = 2.5 ;
    static const double kPhiOpeningAngle = 6 ; 
    
    static const double kMinThetaElectron = 15.;
    static const double kMaxThetaElectron = 45.;
    static const double kMaxThetaHadrons = 140.;
    static const double kMinThetaProton = 10.;
    static const double kMinThetaPiPlus = 10.;
    static const double kMinThetaPiMinus = 22.;
    static const double kMinThetaGamma = 8.;

    static const double kMinEThetaSlice = 36 ;
    static const double kMaxEThetaSlice = 39 ; 

    static const double kCenterFirstSector = 30;
    static const double kCenterSecondSector = 90;
    static const double kCenterThirdSector = 150;
    static const double kCenterFourthSector = 210;
    static const double kCenterFifthSector = 270;
    static const double kCenterSixthSector = 330;

    static const double kPhotonRadCut = 40 ; 
    static const double kPhotonEPhiDiffCut = 30 ; 
  }
}

#endif
