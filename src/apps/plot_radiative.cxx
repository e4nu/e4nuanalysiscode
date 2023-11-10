// __________________________________________________________________________
/* This app is used to plot quantities related to radiative corrections    */
// __________________________________________________________________________
#include <iostream>
#include <vector>
#include "utils/RadiativeCorrUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"

using namespace std;
using namespace e4nu;
using namespace e4nu::utils;
using namespace e4nu::conf; 

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) {

  std::cout << "Plotting e4nu analysis radiative corrections..." << std::endl;

  double EBeam = 1; // To be configurable 
  double W = kNucleonMass;
  double maxQ2 = 2*kNucleonMass*(EBeam-kElectronMass);
  double minQ2 = 0.01;
  double nbins = 100; 
  double binwidth= (maxQ2-minQ2)/nbins;
  vector<double> Q2bins, QELVertex, QELVaumm ; 
  for( unsigned int i = 0 ; i < nbins ; ++i ) { 
    double Q2 = minQ2 + i*binwidth ;
    Q2bins.push_back(Q2);
    QELVertex.push_back(utils::RadCorrQELVertex(Q2));
    QELVaumm.push_back(utils::RadCorrQELVacumm(Q2));
  }

  TGraph * gRadCorrQELVertex = new TGraph( nbins, &Q2bins[0], &QELVertex[0] );
  TGraph * gRadCorrQELVertex = new TGraph( nbins, &Q2bins[&QEL], x[0] );

  return 0 ;
}
