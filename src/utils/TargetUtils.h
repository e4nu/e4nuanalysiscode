/**
 * This file contains utils specific for target properties
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _TARGET_UTILS_H_
#define _TARGET_UTILS_H_

namespace e4nu { 
  namespace utils
    {
      double GetECalOffset( const unsigned int target_pdg ) ;
      double GetBindingEnergy( const unsigned int target_pdg ) ; 
      double GetTargetMass( const unsigned int target_pdg ) ; 
      double GetResidualTargetMass( const unsigned int target_pdg ) ; 
      unsigned int GetTargetNProtons( const unsigned int target_pdg ); 
      unsigned int GetTargetNNeutrons( const unsigned int target_pdg );
      unsigned int GetMassNumber( const unsigned int target_pdg );
    }
}

#endif
