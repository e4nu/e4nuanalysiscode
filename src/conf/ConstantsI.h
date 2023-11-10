/**                                     
 * \info This script contains general information on generic constants
 **/

#ifndef _CONST_I_H_
#define _CONST_I_H_

namespace e4nu {
  namespace conf {
    static const double kConversionFactorCm2ToMicroBarn = TMath::Power(10.,30.); // cm^2 to μbarn
    // 1e -> 1.6x10^-19 C
    // 1C -> 6.25x10^18 e
    // 1mC -> 6.25x10^15 e
    // Thus the numbers above have to be multiplied by this number to make sure that we refer to electrons and not charge
    static const double kConversionFactorChargeToElectrons = 6.25*TMath::Power(10.,15.);
    // Avogadro constant: 6x10^23
    // number of atoms in 12 grams of the isotope 12C
    // 1 gr -> 6x10^23 / 12 = 5x10^22 atoms
    static const double kAvogadroNumber = 6*TMath::Power(10.,23);
    static const double kOverallUnitConversionFactor = kConversionFactorChargeToElectrons * kAvogadroNumber;
    static const double kAem = 1./137.03599976; // EM coupling const, dimensionless 
    static const double kAem2  = TMath::Power(kAem,2);
    static const double kPi = TMath::Pi(); 
    
  }
}
#endif
