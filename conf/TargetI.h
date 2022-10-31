/** 
 * \info This script contains the constant binding energy parameters
 **/


#ifndef _TARGET_I_H_
#define _TARGET_I_H_

namespace e4nu {
  namespace TargetI { 

    const unsigned int kPdgHe3 = 1000020030; 
    const unsigned int kPdgHe4 = 1000020040; 
    const unsigned int kPdgC12 = 1000060120 ; 
    const unsigned int kPdgFe56 = 1000260560 ; 
    const unsigned int kPdgCH2 ; //= 100ZZZAAAI; 
    const unsigned int kPdgD = 1000010020 ; 
    const unsigned int kPdgO16 = 1000080160 ; 
    const unsigned int kPdgFreeP = 1000010010 ;
    const unsigned int kPdgFreeN = 1000000010 ;

    // Binding energy in GeV
    const unsigned double kBEH = 0.008481 ; 
    const unsigned double kBEHe3 = 0.0077 ;
    const unsigned double kBEHe4 = 0.0283 ; 
    const unsigned double kBeD2 = 0.00222 ; 
    const unsigned double kBEC12 = 0.09215 ; 
    const unsigned double kBEFe56 = 0.49226 ; 
    const unsigned double kBEB = 0.0762 ;
    const unsigned double kBEMn = 0.4820764 ; 

    // ECalOffset
    const unsigned double kECalOffsetHe3 = 0.004 ; 
    const unsigned double kECalOffsetHe4 = 0.005 ; 
    const unsigned double kECalOffsetC12 = 0.005 ; 
    const unsigned double kECalOffsetFe56 = 0.011 ; 
  }
}
#endif _TARGET_I_H_
