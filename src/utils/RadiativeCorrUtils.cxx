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

double utils::SIMCEnergyLoss(const TLorentzVector particle, const int p_pdg, const double tgt_pdg, const double thickness, const double max_Ephoton ) {
  // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.054610
  double b = SIMCBFactor( tgt_pdg );
  double lambda = TMath::Log(4*pow(particle.P(),2)/pow(utils::GetParticleMass(p_pdg),2)) - 1 ;
  if( particle.Pz() != particle.E() ) lambda += TMath::Log(0.5*(1-particle.CosTheta())) ;
  lambda += 2*TMath::Log(4.325/particle.E());
  lambda *= (kAem/kPi) ;
  lambda += b*thickness;
  if( lambda < 0 ) return 0; 

  double e_gamma_max = max_Ephoton*particle.E() ;
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

double utils::SimpleEnergyLoss(const TLorentzVector electron, const double tgt_pdg, const double thickness, const double max_Ephoton ) {
  // Reference https://github.com/adishka/Generator/blob/adi_radcorr/src/Physics/Common/RadiativeCorrector.cxx
  double b = SIMCBFactor( tgt_pdg ); 
  double lambda = (kAem/kPi)*( 2*TMath::Log(2*electron.E()/kElectronMass) - 1 ) + b*thickness;
  //  if( particle.Pz() != particle.E() ) lambda += TMath::Log(0.5*(1-particle.CosTheta())) ;
  double e_gamma_max = max_Ephoton*electron.E() ;
  double e_gamma_min = 1E-25;

  TF1 *f = new TF1("f","[0]*pow(x,[0]-1)/[1]",e_gamma_min,e_gamma_max);
  f->SetParameter(0,lambda);
  f->SetParameter(1,pow(electron.E(),-1.*lambda));
  double energyLoss = f->GetRandom();
  delete f;
  return energyLoss ; 
}

TLorentzVector utils::RadOutElectron( const TLorentzVector electron_vertex, TLorentzVector & out_electron, const int tgt, const double thickness, const double max_Ephoton, const string model ) {
  // Compute true detected outgoing electron kinematics with energy loss method
  double egamma = 0 ; 
  //if( model == "simc" || model == "schwinger" || model == "vanderhaeghen" || model == "motsai" || model == "myversion" ) {
  egamma = SIMCEnergyLoss( electron_vertex, 11, tgt, thickness, max_Ephoton ) ;
  //}
  //  else if ( model == "simple" ) egamma = SimpleEnergyLoss( electron_vertex, tgt, thickness, max_Ephoton ) ; 
  if( egamma < 0 ) egamma = 0 ;
  TLorentzVector OutGamma = GetEmittedHardPhoton( electron_vertex, egamma ) ;
  if( OutGamma.E() < 0 )  OutGamma.SetPxPyPzE(0,0,0,0);
  out_electron = electron_vertex - OutGamma ; 
  return OutGamma;
} 

double utils::SIMCRadCorrWeight( const e4nu::Event & event, const double thickness, const double max_Ephoton, const std::string model ) {

  double weight = 1; 
  TLorentzVector Beam  = event.GetInLepton4Mom() ;
  TLorentzVector Detected = event.GetOutLepton4Mom() ;
  TLorentzVector InRad = event.GetInCorrLepton4Mom() ; 
  TLorentzVector OutRad = event.GetOutCorrLepton4Mom();
  double Q2 = event.GetTrueQ2s();
  unsigned int tgt = event.GetTargetPdg();
  double Emax = 0.2*Beam.E();
  double Emin = 1E-15;
  TF1 * fsp = new TF1("fsp","TMath::Log(1-x)* (1./x)"); // Spence function

  if( model == "simple") { 
    // Reference ?
    weight = 1 + (2*kAem /kPi) * ( (13./12.)* (TMath::Log(Q2/pow(kElectronMass,2)) - 1) - (17./36.)
				   - (1./4.) * pow(TMath::Log(OutRad.E()*Detected.E()),2)
				   - (1./2.) * ( (pow(kPi,2)/6) -  TMath::DiLog(TMath::Power(TMath::Cos(0.5*Detected.Theta()),2.)) ) );
  } else if ( "schwinger" ) { 
    // 10.1103/PhysRev.76.790
    weight = 1 + (2*kAem / kPi) *( ( TMath::Log( Beam.P() / Emax ) - (13./12.))* (TMath::Log(Q2/pow(kElectronMass,2)) - 1) + (17./36)); 
  } else if ( "vanderhaeghen" ) { 
    // 10.1103/physrevc.62.025501
    double e_gamma_min = 1E-25;
    double SP = -1*fsp->Integral(Emin,pow(TMath::Cos(Detected.Theta())/2,2.));
    double delta = TMath::Log(pow(Emax,2)/(Beam.E()*Detected.E()))* (TMath::Log(Q2/pow(kElectronMass,2)) - 1) + 13./9.*TMath::Log(Q2/pow(kElectronMass,2)) ;
    delta -= ( 29./9. + 0.5 * pow(TMath::Log(Beam.E()/Detected.E()),2) + pow(kPi,2)/6.);
    delta += SP; 
    delta *= (kAem/kPi);
    weight = 1+delta;

  } else if ( "motsai" ) {
    // 10.1103/RevModPhys.41.205
    // https://inspirehep.net/files/1fcaa81f63f50d7bf56a22ce2c6b8b58 Equation II.2
    double SP = fsp->Integral(Emin,pow(-TMath::Sin(Detected.Theta())/2,2.));
    double Fth = TMath::Log(pow(TMath::Sin(Detected.Theta()/2),2.))*TMath::Log(pow(TMath::Cos(Detected.Theta()/2),2.))*SP;
    double delta = (TMath::Log(Beam.E()/Emax) - 13./12.) * (TMath::Log(Q2/pow(kElectronMass,2)) - 1) + 17./36. + 0.5*Fth;
    delta *= - 2*(kAem/kPi);
    weight = 1 - delta ; 
  } else if ( "myversion" ) { 
    // 10.1103/PhysRevC.64.054610
    double delta_hard = 2.*(kAem/kPi)*(-3./4*TMath::Log(Q2/pow(kElectronMass,2)+1)+5./9.-1./3*TMath::Log(Q2/pow(kElectronMass,2))); 
    weight = 1 - delta_hard ; 
  } else if ( "simc" ) { 
    // This takes into account the radiation weight due to external and internal radiation of the incoming electron
    // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.64.054610
    double delta_hard = 2.*(kAem/kPi)*(-3./4*TMath::Log(Q2/pow(kElectronMass,2))+1+5./9.-1./3.*TMath::Log(Q2/pow(kElectronMass,2))); 

    double b = SIMCBFactor( tgt );
    double bt = b*thickness;

    // Calculate correction for external radiation
    //double lambda_0 = (kAem/kPi)*(2*TMath::Log(Beam.P()/Detected.P())+TMath::Log(0.5*(1-TMath::Cos(Detected.Theta()))));
    double lambda_0_e = (kAem/kPi)*(2*TMath::Log(Beam.P()/InRad.P()));
    //+TMath::Log(0.5*(1-TMath::Cos(InRad.Theta()))));
    double lambda_0_el = (kAem/kPi)*(2*TMath::Log(OutRad.P()/Detected.P())+TMath::Log(0.5*(1-TMath::Cos(Detected.Theta()))));
    double lambda_e = (kAem/kPi)*(TMath::Log(pow(2*Beam.P()/kElectronMass,2)-1))+lambda_0_e;
    double lambda_el = (kAem/kPi)*(TMath::Log(pow(2*Detected.P()/kElectronMass,2)-1))+lambda_0_el;
    double g_e = lambda_e + bt;
    double g_el = lambda_el + bt;

    double Phi_ext_e = 1. ;
    double Phi_ext_el = 1. ;
    double EPhoton_i = Beam.E() - InRad.E() ; // should be energy of decayed electron
    if( EPhoton_i != 0 ) Phi_ext_e -= bt * EPhoton_i / Beam.P() / g_e ;  
    double EPhoton_f = OutRad.E() - Detected.E();
    if( EPhoton_f != 0 ) Phi_ext_el -= bt * EPhoton_f / OutRad.P() / g_el ;  

    double e_gamma_max_e = max_Ephoton*Beam.E() ;
    double e_gamma_max_el = max_Ephoton*OutRad.E() ;
    double e_gamma_min = 1E-25;
    double power_hi_e = pow(e_gamma_max_e,g_e);
    double power_lo_e  = pow(e_gamma_min,g_e);
    double power_hi_el = pow(e_gamma_max_el,g_el);
    double power_lo_el  = pow(e_gamma_min,g_el);

    //    double C_e = g_e/(TMath::Gamma(1+bt)*pow(Beam.P(),bt)*pow(Beam.P()*Detected.P(),lambda_e/2)); 
    double C_e = g_e/(TMath::Gamma(1+bt)*pow(Beam.P(),bt)*pow(Beam.P()*InRad.P(),lambda_e/2)); 
    double C_el = g_el/(TMath::Gamma(1+bt)*pow(OutRad.P(),bt)*pow(OutRad.P()*Detected.P(),lambda_el/2)); 

    double W_e = 1;//(C_e/g_e)*(power_hi_e-power_lo_e);
    double W_el = 1;//(C_el/g_el)*(power_hi_el-power_lo_el);

    weight= W_e*W_el*Phi_ext_e*Phi_ext_el*(1-delta_hard);
  }
  delete fsp ;  

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
