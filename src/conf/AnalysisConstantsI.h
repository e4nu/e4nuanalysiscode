/**                                                                                                                                                                                                                  
 * This file contains general constants to be used in any analysis
 **/

#ifndef _ANALYSISCONST_I_H_
#define _ANALYSISCONST_I_H_

namespace e4nu {
  namespace conf {
    // Phase space limits
    const double kQ2Max = 1.75 ; 
    const double kQ2Min = 0.4 ; 
    const double kMinTheta = 10 ;
    const double kMaxTheta = 60 ;
    const double kMinEePrime = 0.5 ; 
    const double kMaxEePrime = 2.5 ;
    const double kPhiOpeningAngle = 6 ; 
    
    const double kMinThetaProton = 12.;
    const double kMinThetaPiPlus = 12.;
    const double kMinThetaPiMinus = 0.;
    const double kMinThetaGamma = 8.;

    const double kMinEThetaSlice = 36 ;
    const double kMaxEThetaSlice = 39 ; 

    const double kCenterFirstSector = 30;
    const double kCenterSecondSector = 90;
    const double kCenterThirdSector = 150;
    const double kCenterFourthSector = 210;
    const double kCenterFifthSector = 270;
    const double kCenterSixthSector = 330;

    const double kPhotonRadCut = 40 ; 
    const double kPhotonEPhiDiffCut = 30 ; 
  }
}

#endif
