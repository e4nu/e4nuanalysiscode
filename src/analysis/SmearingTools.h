/**
 * This file contains smearing tools
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date October 2022
 **/

#ifndef _SMEARING_TOOLS_H_
#define _SMEARING_TOOLS_H_

#include "TLorentzVector.h"

namespace e4nu { 
  namespace analysis
    {
      bool ApplySmearing( TLorentzVector & particle, const int pdg, const bool is_MC = true ) ;
     }
}

#endif 
