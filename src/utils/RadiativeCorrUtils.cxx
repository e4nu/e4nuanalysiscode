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

double utils::SIMCEnergyLoss(const double EBeam, const TLorentzVector particle, const int p_pdg, const double tgt_pdg, const double thickness, const double max_Ephoton ) {
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.054610
  double Z = utils::GetTargetNProtons(tgt_pdg);
  double L1 = TMath::Log(184.15) - (1./3)*TMath::Log(Z);
  double L2 = TMath::Log(1194.) - (2./3)*TMath::Log(Z);
  if( Z ==1 ) { 
    L1 = 5.31;
    L2 = 6.144;
  }
  double b = (1./9)*(12 + (Z+1)/(Z*L1 + L2));
  
  double lambda = (kAem/kPi)*( 2*TMath::Log(2*TMath::Sqrt(pow(EBeam,2)-pow(kElectronMass,2))/utils::GetParticleMass(p_pdg)) - 1 );//+ TMath::Log(0.5*(1-particle.CosTheta())) ) ;
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
  double Z = utils::GetTargetNProtons(tgt_pdg);
  double L1 = TMath::Log(184.15) - (1./3)*TMath::Log(Z);
  double L2 = TMath::Log(1194.) - (2./3)*TMath::Log(Z);
  if( Z ==1 ) { 
    L1 = 5.31;
    L2 = 6.144;
  }
  double b = (1./9)*(12 + (Z+1)/(Z*L1 + L2));
  
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
  std::cout << energyLoss << std::endl;
  return energyLoss ; 
}

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

double utils::RadCorrQELRealProtonD1( const double Q2, const double E, const double Ep, const double theta ) {
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.62.025501  
  double rho = Q2 + 4 * pow( kNucleonMass,2 ) ; 
  double x   = pow( sqrt(Q2) + rho ,2 ) * 0.25 / pow( kNucleonMass,2 ) ;
  double ELab = (kNucleonMass*E-0.5*Q2)/kNucleonMass; // A43 // Eetilde
  double EpLab = (kNucleonMass*Ep+0.5*Q2)/kNucleonMass; // A43 // Eetilde

  // Epel denotes the elastic scattered electron lab energy to distinguish from Ep
  double deltaEp = VanderhagenELoss( Q2, Ep );  
  double Epel = Ep + deltaEp ;
  double deltaELab = E - Ep - (E-Epel)*Ep/Epel; // A47, in terms of lab quantities 
  double eta = deltaELab / (Epel-Ep); //A48
  double DeltaEs = eta*( Epel - Ep ); 

  double e_gamma_min = 1E-25;
  TF1 * fsp = new TF1("fsp","TMath::Log(1-x)* (1./x)");
  double Sp = -1*fsp->Integral(e_gamma_min,pow(cos(theta)/2,2.));
  delete fsp ;

  double delta1 = TMath::Log(4*pow(DeltaEs,2)/Q2/x)*TMath::Log(eta) + Sp * (1-eta/x) - Sp * ( 1-1/eta/x) ; 
  delta1 *= 2*kAem/kPi;
  return delta1; 
}

double utils::RadCorrQELRealProtonD20( const double Q2, const double E, const double Ep, const double theta, const double deltaE, const TLorentzVector nucleon ) {
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.62.025501  
  double rho = Q2 + 4 * pow( kNucleonMass,2 ) ; 
  double x   = pow( sqrt(Q2) + rho ,2 ) * 0.25 / pow( kNucleonMass,2 ) ;
  double ELab = (kNucleonMass*E-0.5*Q2)/kNucleonMass; // A43 // Eetilde
  double EpLab = (kNucleonMass*Ep+0.5*Q2)/kNucleonMass; // A43 // Eetilde

  // Epel denotes the elastic scattered electron lab energy to distinguish from Ep
  double deltaEp = VanderhagenELoss( Q2, Ep );  
  double Epel = Ep + deltaEp ;
  double deltaELab = E-Ep - (E-Epel)*Ep/Epel; // A47, in terms of lab quantities 
  double eta = deltaELab / (Epel-Ep); //A48
  double DeltaEs = eta*( Epel - Ep ); 

  double e_gamma_min = 1E-25;
  TF1 * fsp = new TF1("fsp","TMath::Log(1-x)* (1./x)");
  double Sp = -1*fsp->Integral(e_gamma_min,pow(cos(theta)/2,2.));
  delete fsp ;

  double delta20 = TMath::Log(4*pow(DeltaEs,2)/pow(kNucleonMass,2))*((nucleon.E()/nucleon.Mag())*TMath::Log(x)-1) + 1 ;
  delta20 += (nucleon.E()/nucleon.Mag()) *(-0.5*pow(TMath::Log(x),2)- TMath::Log(x) * TMath::Log(pow(rho,2)/pow(kNucleonMass,2)) + TMath::Log(x) - Sp *( 1- 1/pow(x,2)) -2*Sp/x +pow(kPi,2)/6. );  
  delta20 *= kAem/kPi;
  return delta20;
}
 
double utils::RadCorrQELRealProtonD21( const double Q2, const double E, const double Ep, const double theta) {
  return 0; 
}
