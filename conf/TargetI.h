/** 
 * \info This script contains the constant binding energy parameters
 **/


#ifndef _TARGET_I_H_
#define _TARGET_I_H_

namespace e4nu {
  namespace conf { 

    const unsigned int kPdgHe3 = 1000020030; 
    const unsigned int kPdgHe4 = 1000020040; 
    const unsigned int kPdgC12 = 1000060120 ; 
    const unsigned int kPdgFe56 = 1000260560 ;
    const unsigned int kPdgD = 1000010020 ; 
    const unsigned int kPdgO16 = 1000080160 ; 
    const unsigned int kPdgFreeP = 1000010010 ;
    const unsigned int kPdgFreeN = 1000000010 ;

    // Binding energy in GeV
    const double kBEH = 0.008481 ; 
    const double kBEHe3 = 0.0077 ;
    const double kBEHe4 = 0.0283 ; 
    const double kBED2 = 0.00222 ; 
    const double kBEC12 = 0.09215 ; 
    const double kBEFe56 = 0.49226 ; 
    const double kBEB = 0.0762 ;
    const double kBEMn = 0.4820764 ; 

    // ECalOffset
    const double kECalOffsetHe3 = 0.004 ; 
    const double kECalOffsetHe4 = 0.005 ; 
    const double kECalOffsetC12 = 0.005 ; 
    const double kECalOffsetFe56 = 0.011 ; 
  }
}
#endif
