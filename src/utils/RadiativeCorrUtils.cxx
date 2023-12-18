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

using namespace std;
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
  double lambda = 2*TMath::Log(2*TMath::Sqrt(pow(EBeam,2)-pow(kElectronMass,2))/utils::GetParticleMass(p_pdg)) - 1 ;
  if( particle.Pz() != particle.E() ) lambda += TMath::Log(0.5*(1-particle.CosTheta())) ;
  if( p_pdg == kPdgProton ) lambda = ( TMath::Log((particle.E()+particle.P())/(particle.E()-particle.P())) - 2 ) ;
  lambda *= (kAem/kPi) ;
  lambda += b*thickness;
  if( lambda < 0 ) return 0; 

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

  TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
  f->SetParameter(0,lambda);
  f->SetParameter(1,pow(EBeam,-1.*lambda));
  double energyLoss = f->GetRandom();
  delete f;
  return energyLoss ; 
}

TLorentzVector utils::SIMCRadCorrQELRadOutElectron( const TLorentzVector electron_vertex, TLorentzVector & out_electron, const int tgt, const double thickness, const double max_Ephoton, const string model ) {
  // Compute true detected outgoing electron kinematics with energy loss method
  double egamma = 0 ; 
  if( model == "simc" ) egamma = SIMCEnergyLoss( electron_vertex.E(),electron_vertex, 11, tgt, thickness, max_Ephoton ) ;
  else if ( model == "simple" ) egamma = SimpleEnergyLoss( electron_vertex.E(), tgt, thickness, max_Ephoton ) ; 
  if( egamma < 0 ) egamma = 0 ;
  TLorentzVector OutGamma = GetEmittedHardPhoton( electron_vertex, egamma ) ;
  if( OutGamma.E() < 0 )  OutGamma.SetPxPyPzE(0,0,0,0);
  out_electron = electron_vertex - OutGamma ; 
  return OutGamma;
} 

double utils::SIMCRadCorrWeight( const TLorentzVector corr_electron, const TLorentzVector detected_electron, const double EBeam, const double Q2, const int tgt, const double thickness, const double max_Ephoton, const std::string model ) {

  double weight = 1; 
  if( model == "simple") { 
    // Reference ?
    weight = 1 + (2*kAem /kPi) * ( (13./12.)* (TMath::Log(Q2/pow(kElectronMass,2)) - 1) - (17./36.)
				   - (1./4.) * pow(TMath::Log(corr_electron.E()*detected_electron.E()),2)
				   - (1./2.) * ( (pow(kPi,2)/6) -  TMath::DiLog(TMath::Power(TMath::Cos(0.5*detected_electron.Theta()),2.)) ) );
  } else if ( "simc" ) { 
    
    // This takes into account the radiation weight due to external and internal radiation of the incoming electron
    // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.64.054610
    double EPhoton = EBeam - corr_electron.E();
    double b = SIMCBFactor( tgt );
    double e_mom = TMath::Sqrt(pow(EBeam,2)-pow(kElectronMass,2)) ; 
    double lambda_e = (kAem/kPi)*( 2*TMath::Log(2*e_mom/kElectronMass) - 1 ) ;
    double bt = b*thickness;
    double g = lambda_e + bt;
    double C = 1/(TMath::Gamma(1+bt)*pow(e_mom,bt)*pow(e_mom*corr_electron.P(),lambda_e/2)); 
    double e_gamma_max = max_Ephoton*EBeam ;
    double e_gamma_min = 1E-25;
    double power_hi = pow(e_gamma_max,g);
    double power_lo  = pow(e_gamma_min,g);
    double W_rad_e = C*(power_hi-power_lo);
    double Phi_ext_e = 1. ;
    if( EPhoton != 0 ) Phi_ext_e -= ( (b*thickness) / EBeam / g * EPhoton) ;
  
    std::cout <<  " W_rad_e = " << W_rad_e <<  " Phi_ext_e = " << Phi_ext_e << std::endl;

    // This takes into account the radiation weight due to external and internal radiation of the incoming electron
    // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.64.054610
    //double delta_hard = -1.*(kAem/kPi)*(-28./9.+(13./6.)*TMath::Log(Q2/pow(kElectronMass,2))); // Is this ok?
    double delta_hard = 2.*(kAem/kPi)*(-3./4*TMath::Log(Q2/pow(kElectronMass,2))+1+5./9.-1./3.*TMath::Log(Q2/pow(kElectronMass,2))); // Is this ok?
    std::cout << (1-delta_hard)<<std::endl;
    // Compute weight
    lambda_e = (kAem/kPi)*( 2*TMath::Log(2*corr_electron.P()/kElectronMass) - 1 + TMath::Log(0.5*(1-corr_electron.CosTheta()))) ;
    g = lambda_e + b*thickness;
    C = g/(TMath::Gamma(1+b*thickness)*pow(corr_electron.P(),b*thickness)*pow(corr_electron.P()*detected_electron.P(),lambda_e/2)); 
    double W_rad_el = (C/g)*(power_hi-power_lo);
    double Phi_ext_el = 1. ;
    if( EPhoton != 0 ) Phi_ext_el -= ( (b*thickness) / EBeam / g * EPhoton) ;  

    weight= W_rad_e*Phi_ext_e*W_rad_el*Phi_ext_el*(1-delta_hard);

  }
  return weight;
}

TLorentzVector utils::GetEmittedHardPhoton( TLorentzVector electron, double eloss ) {

  double ptLoss;
  double pzLoss;
  if (electron.Pz()==electron.E()) // for the z direction going probe theta = -nan 
    {
      ptLoss = 0.;
      pzLoss = eloss;
    }
  else
    {
      ptLoss = eloss*sin(electron.Theta()); 
      pzLoss = eloss*cos(electron.Theta());
    }

  TLorentzVector p4RadGamma;
  p4RadGamma.SetPxPyPzE(ptLoss*cos(electron.Phi()),ptLoss*sin(electron.Phi()),pzLoss,eloss);
  return p4RadGamma ;
} 
