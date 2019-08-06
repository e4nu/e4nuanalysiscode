#include "e2a_ep_neutrino6_united4_radphot.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{


  if( argc < 3 ){
    cout<<"Please specify the target (3He, 56Fe, C12, 4He) and the beam energy (2261 or 4461)"<<endl;
    cout<<"================= Usage ==============="<<endl;
    cout<<"./run_e2a_ep_neutrino6_united4_radphot.cc target beam_energy"<<endl;
    exit(1);
  }


  std::string target  = argv[1];
  std::string beam_en = argv[2];

  e2a_ep_neutrino6_united4_radphot  t(target,beam_en);

  t.Loop();

  return 0;
}
