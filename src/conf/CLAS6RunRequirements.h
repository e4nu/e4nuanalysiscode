/**
 * This file contains cuts on CLAS6 runs 
 * \author Julia Tena Vidal \at Tel Aviv University
 * \date Jun 2025
 **/

#ifndef _RUNREQ_I_H_
#define _RUNREQ_I_H_

#include <vector>
#include <string>
namespace e4nu { 
  namespace conf {

    std::vector<unsigned int> GetInvalidRuns( const double Beam, const std::string Target ) ;
    
  }
}

#endif
