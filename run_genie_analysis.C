#ifndef rungenieanalysis
#define rungenieanalysis
#include "genie_analysis.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{


  if( argc < 3 ){
    cout<<"Please specify the target (3He, 56Fe, C12, 4He) and the beam energy (2261 or 4461)"<<endl;
    cout<<"================= Usage ==============="<<endl;
    cout<<"./genie_analysis target beam_energy"<<endl;
    exit(1);
  }


  std::string target  = argv[1];
  std::string beam_en = argv[2];


  genie_analysis  t(target,beam_en);
  t.Loop();

  return 0;
}
#endif
