/**
 * This file contains utils which aim to study the effect of radiative corrections
 * qualitively before included in the GENIE event generator
 * It computes the correction factors to the cross section due to vertex, vacum or radiative effects
 * The methods used depend on the exclusive final state measured
 * References are provided for each case
 * \author Julia Tena Vidal \at Tel Aviv University                                                                                                                                                                 
 * \date Nov 2023                                                                                                                                                                                              
 **/
#include <iostream>
#include <TMath.h>
#include <TF1.h>
#include "utils/RadiativeCorrUtils.h"
#include "utils/TargetUtils.h"
#include "utils/ParticleUtils.h"
#include "conf/ConstantsI.h"
#include "conf/ParticleI.h"

using namespace e4nu;
using namespace conf;

double utils::VanderhagenELoss( const double Q2 , const double Ee ) {
  double e_gamma_min = 1E-25;
  double e_gamma_max = 0.2*Ee ;
  TF1 *f = new TF1("f","([0]/x)*TMath::Power(x/[1],[0])",e_gamma_min,e_gamma_max);
  double a = (kAem/kPi)*(TMath::Log(Q2)/pow(kElectronMass,2) - 1.);
  f->SetParameter(0,a);
  f->SetParameter(1,Ee);
  double energyLoss = f->GetRandom();
  delete f;
  return energyLoss ; 
}

double utils::RadCorrQELVertex( const double Q2 ) {
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.62.025501
  // A67 equation
  double log = TMath::Log(Q2/pow(kElectronMass,2)) ;
  double dvertex = (3./2.)*log-2-0.5*pow(log,2)+pow(kPi,2)/6.;
  dvertex *= (kAem/kPi);
  return dvertex;
}
 
double utils::RadCorrQELVacumm( const double Q2 ) {
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.62.025501
  // A69 equation
  double log = TMath::Log(Q2/pow(kElectronMass,2)) ;
  double dvaccumm = log - 5./3.;
  dvaccumm *= (kAem/kPi)*2./3.;
  return dvaccumm;  
}

double utils::RadCorrQELRealRad( const double Q2, const double E, const double Ep, const double theta ) {
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.62.025501
  // A65 equation
  double e_gamma_min = 1E-25;
  TF1 * fsp = new TF1("fsp","TMath::Log(1-x)* (1./x)");
  double SP = -1*fsp->Integral(e_gamma_min,pow(cos(theta)/2,2.));
  delete fsp ; 

  double log = TMath::Log(Q2/pow(kElectronMass,2)) ;
  double deltaE = VanderhagenELoss( Q2, E );
  double dreal =TMath::Log(pow(deltaE,2)/(E*Ep))*(log - 1); 
  dreal -= 0.5 * pow(TMath::Log(E/Ep),2);
  dreal += 0.5 * pow(log,2);
  dreal -= pow(kPi,2)/3. ;
  dreal += SP;
  dreal *= kAem/kPi;
  return dreal ; 
}

double utils::SIMCBFactor( const double tgt_pdg ) { 
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.054610
  double Z = utils::GetTargetNProtons(tgt_pdg);
  double L1 = TMath::Log(184.15) - (1./3)*TMath::Log(Z);
  double L2 = TMath::Log(1194.) - (2./3)*TMath::Log(Z);
  if( Z ==1 ) { 
    L1 = 5.31;
    L2 = 6.144;
  }
  double b = (1./9)*(12 + (Z+1)/(Z*L1 + L2));
  return b ;
}

double utils::SIMCEnergyLoss(const double EBeam, const TLorentzVector particle, const int p_pdg, const double tgt_pdg, const double thickness, const double max_Ephoton ) {
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.054610
  double b = SIMCBFactor( tgt_pdg );
  double lambda = (kAem/kPi)*( 2*TMath::Log(2*EBeam)/kElectronMass - 1 );//+ TMath::Log(0.5*(1-particle.CosTheta())) ) ;
  if( p_pdg == kPdgProton ) lambda = (kAem/kPi)*( TMath::Log((particle.E()+particle.P())/(particle.E()-particle.P())) - 2 ) ;
  lambda += b*thickness;

  double e_gamma_max = max_Ephoton*EBeam ;
  double e_gamma_min = 1E-25;
  double power_hi = pow(e_gamma_max,lambda);
  double power_lo  = pow(e_gamma_min,lambda);

  TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
  f->SetParameter(0,lambda);
  f->SetParameter(1,power_hi - power_lo);
  double energyLoss = f->GetRandom();
  delete f;
  return energyLoss ; 
}

double utils::SimpleEnergyLoss(const double EBeam, const double tgt_pdg, const double thickness, const double max_Ephoton ) {
  // Reference https://github.com/adishka/Generator/blob/adi_radcorr/src/Physics/Common/RadiativeCorrector.cxx
  double b = SIMCBFactor( tgt_pdg );  
  double lambda = (kAem/kPi)*( 2*TMath::Log(2*EBeam)/kElectronMass - 1 ) + b*thickness;

  double e_gamma_max = max_Ephoton*EBeam ;
  double e_gamma_min = 1E-25;
  double power_hi = pow(e_gamma_max,lambda);
  double power_lo  = pow(e_gamma_min,lambda);

  TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
  f->SetParameter(0,lambda);
  f->SetParameter(1,pow(EBeam,-1.*lambda));
  double energyLoss = f->GetRandom();
  delete f;
  return energyLoss ; 
}


double utils::SIMCRadCorrQELRadInElectron( const double EBeam, const TLorentzVector electron, const int tgt, const double thickness, const double max_Ephoton ) {
  // This takes into account the radiation weight due to external and internal radiation of the incoming electron
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.64.054610
  double EPhoton = EBeam - electron.E();
  double b = SIMCBFactor( tgt );
  double lambda_e = (kAem/kPi)*( 2*TMath::Log(2*EBeam/kElectronMass) - 1 ) ;
  double g = lambda_e + b*thickness;
  double C = g/(TMath::Gamma(1+b*thickness)*pow(EBeam,b*thickness)*pow(electron.P()*EBeam,lambda_e/2)); 
  double e_gamma_max = max_Ephoton*EBeam ;
  double e_gamma_min = 1E-25;
  double power_hi = pow(e_gamma_max,lambda_e);
  double power_lo  = pow(e_gamma_min,lambda_e);
  double W_rad_e = (C/g)*(power_hi-power_lo);
  double Phi_ext_e = 1. - ( (b*thickness) / EBeam / g * EPhoton) ;
  double radcor_weight = W_rad_e*Phi_ext_e;
  return radcor_weight ;
} 


double utils::SIMCRadCorrQELRadOutElectron( const TLorentzVector electron_vertex, TLorentzVector & out_electron, const double EBeam, const double Q2, const int tgt, const double thickness, const double max_Ephoton ) {
  // This takes into account the radiation weight due to external and internal radiation of the incoming electron
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.64.054610
  double delta_hard = -1.*(kAem/kPi)*(-28./9.+(13./6.)*TMath::Log(Q2/pow(kElectronMass,2)));
  // decay electron here 

  TLorentzVector photon = electron_vertex - out_electron ;

  double lambda_e = (kAem/kPi)*( 2*TMath::Log(2*electron_vertex.P()/kElectronMass) - 1 + TMath::Log(0.5*(1-electron_vertex.CosTheta()))) ;
  double b = SIMCBFactor( tgt );
  double g = lambda_e + b*thickness;
  double C = g/(TMath::Gamma(1+b*thickness)*pow(electron_vertex.P(),b*thickness)*pow(electron_vertex.P()*out_electron.P(),lambda_e/2)); 
  double e_gamma_max = max_Ephoton*EBeam ;
  double e_gamma_min = 1E-25;
  double power_hi = pow(e_gamma_max,lambda_e);
  double power_lo  = pow(e_gamma_min,lambda_e);
  double W_rad_el = (C/g)*(power_hi-power_lo);
  double Phi_ext_el = 1. - ( (b*thickness) / electron_vertex.E() / g * photon.E()) ;
  
  return W_rad_el*Phi_ext_el*(1-delta_hard);
} 
