#ifndef rune2aep
#define rune2aep
#include "e2a_ep_neutrino6_united4_radphot.h"
// apapadop
#include "FilterData.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{

	if( argc < 4 ){
		cout<<"Please specify the target (3He, 56Fe, C12, 4He), the beam energy (2261 or 4461) and what you want to do (Create Plots = 0, Filter = 1)"<<endl;
		cout<<"================= Usage ==============="<<endl;
		cout<<"./run_e2a_ep_neutrino6_united4_radphot.cc target beam_energy 0/1 (Create Plots = 0, Filter = 1)"<<endl;
		exit(1);
	}

	std::string target	= argv[1];
	std::string beam_en = argv[2];
	int choice = atoi(argv[3]);

	if (choice != 1 && choice != 0) {
		std::cout << "Unknown option for parameter 3. It should be either 0 or 1. The given value is " << choice << std::endl;
		return 0;
	}

	e2a_ep_neutrino6_united4_radphot	t(target,beam_en);
	FilterData	filter(target,beam_en);
	if (choice == 0) {
		t.Loop();
	}

	if (choice == 1) {

		cout << "Filtering sample" << endl;
		filter.Loop();
	}

	return 0;
}
#endif
