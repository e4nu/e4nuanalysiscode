#ifndef rungenieanalysis
#define rungenieanalysis
#include "genie_analysis.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{


  if( argc < 4 ){
    std::cout<<"Please specify the target (3He, 56Fe, C12, 4He), the beam energy (2261 or 4461) and the data type (CLAS=0 or simulation=1)"<<std::endl;
    std::cout<<"================= Usage ==============="<<std::endl;
    std::cout<<"./genie_analysis target beam_energy 0/1"<<std::endl;
    exit(1);
  }


  std::string target  = argv[1];
  std::string beam_en = argv[2];
  int choice = atoi(argv[3]);

  if (choice != 1 && choice != 0) {
    std::cout << "Unknown option for parameter 3. It should be either 0 or 1. The given value is " << choice << std::endl;
    return 0;
  }


  genie_analysis  t(target,beam_en);
  if (choice ==1){
    t.Loop();
  }
  if (choice ==0){
    t.LoopCLAS();
  }

  return 0;
}
#endif
