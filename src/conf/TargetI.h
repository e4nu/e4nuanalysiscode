/** 
 * \info This script contains the constant binding energy parameters
 **/


#ifndef _TARGET_I_H_
#define _TARGET_I_H_

namespace e4nu {
  namespace conf { 

    static const unsigned int kPdgH   = 1000010010; 
    static const unsigned int kPdgHe3 = 1000020030; 
    static const unsigned int kPdgHe4 = 1000020040; 
    static const unsigned int kPdgC12 = 1000060120 ; 
    static const unsigned int kPdgFe56 = 1000260560 ;
    static const unsigned int kPdgD = 1000010020 ; 
    static const unsigned int kPdgO16 = 1000080160 ; 
    static const unsigned int kPdgFreeP = 1000010010 ;
    static const unsigned int kPdgFreeN = 1000000010 ;

    // Binding energy in GeV
    static const double kBEH = 0.008481 ; 
    static const double kBEHe3 = 0.0077 ;
    static const double kBEHe4 = 0.0283 ; 
    static const double kBED2 = 0.00222 ; 
    static const double kBEC12 = 0.09215 ; 
    static const double kBEFe56 = 0.49226 ; 
    static const double kBEB = 0.0762 ;
    static const double kBEMn = 0.4820764 ; 

    // ECalOffset
    static const double kECalOffsetHe3 = 0.004 ; 
    static const double kECalOffsetHe4 = 0.005 ; 
    static const double kECalOffsetC12 = 0.005 ; 
    static const double kECalOffsetFe56 = 0.011 ; 
  }
}
#endif
