#ifndef rungenieanalysis
#define rungenieanalysis
#include "genie_analysis.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{


  if( argc < 5 ){
    std::cout<<"Please specify the target (3He, 56Fe, C12, 4He), the beam energy (2261 or 4461), the data type (CLAS=0 or simulation=1) "<<std::endl;
    std::cout<<"and the number of rotations"<<std::endl;
    std::cout<<"================= Usage ==============="<<std::endl;
    std::cout<<"./genie_analysis target beam_energy 0/1 #rot"<<std::endl;
    exit(1);
  }


  std::string target  = argv[1];
  std::string beam_en = argv[2];
  int choice = atoi(argv[3]);
  int rotations = atoi(argv[4]);

  if (choice != 1 && choice != 0) {
    std::cout << "Unknown option for parameter 3. It should be either 0 or 1. The given value is " << choice << std::endl;
    return 0;
  }

  if (rotations <= 0) {
    std::cout << "Not a valid number for the number of rotations (parameter 4). The given value is " << rotations << std::endl;
    return 0;
  }

  genie_analysis  t(target,beam_en, rotations);
  t.Loop(choice);


  return 0;
}
#endif
