/**
 * This file contains utils specific for detector
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date October 2022                                                                                                                                                                                              
 **/
#include <iostream>
#include "utils/Utils.h"

using namespace e4nu;

void utils::PrintProgressBar( const unsigned int curr_event, const unsigned int total_events ) {

  double progress = (double) curr_event / (double) total_events ; 
  if( progress > 1 ) return ; 

  int barWidth = 70;
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] Event: " << curr_event << " Percentage = " << int(progress * 100.0) << " %\r";
  std::cout.flush();

  std::cout << std::endl;
}
