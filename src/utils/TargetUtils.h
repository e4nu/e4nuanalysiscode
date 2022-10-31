/**
 * This file contains utils specific for target properties
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _TARGET_UTILS_H_
#define _TARGET_UTILS_H_

namespace e4nu { 
  namespace TargetUtils 
    {
      unsigned double GetECalOffset( const unsigned int target_pdg ) ;
      unsigned double GetBindingEnergy( const unsigned int target_pdg ) ; 
      unsigned double GetTargetMass( const unsigned int target_pdg ) ; 
      unsigned double GetResidualTargetMass( const unsigned int target_pdg ) ; 
      unsigned int GetTargetNProtons( const unsigned int target_pdg ); 
      unsigned int GetTargetNProtons( const unsigned int target_pdg );
    }
}

#endif _TARGET_UTILS_H_
