#ifndef SUBTRACTION_H
#define SUBTRACTION_H

#include <TVectorT.h>
#include <TVector3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include "Fiducial.h"


struct Subtraction {
  //initialize Beam Energy String, Target Name String, Map of binding energies and Fiducial cuts
  //N_tot of rotations, q vector is set to (0,0,0)
  void InitSubtraction(std::string in_beam_en, std::string in_target_name, std::map<std::string,double> in_bind_en, int in_nrot, Fiducial *in_fiducial) {
      std::cout << " Subtraction start of initialize " << std::endl;
      fbeam_en = in_beam_en;
      target_name = in_target_name;
      bind_en = in_bind_en;
      fiducialcut = in_fiducial;
      N_tot = in_nrot;
      std::cout << " InitSubtraction: N_tot " << N_tot << " , target_name " << target_name << std::endl;
      std::cout << " Test InitSubtraction Fiducials " << in_fiducial->up_lim1_ec->Eval(60) << std::endl;
      V3q.SetX(0);
      V3q.SetY(0);
      V3q.SetZ(0);
  }

  void  prot2_rot_func(TVector3 V3prot[2],TVector3 V3prot_uncorr[2], TLorentzVector V4el, double Ecal_2pto1p[2],  double pmiss_perp_2pto1p[2],  double  P2pto1p[2],  double *Nboth);
  void  prot3_rot_func(TVector3 V3prot[3],TVector3 V3prot_uncorr[3], TLorentzVector V4el, double Ecal_3pto2p[][2],double pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3], double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det);

  void  prot1_pi1_rot_func(TVector3  V3prot, TVector3 V3pi, int q_pi, double *N_pi_p, double *N_nopi_p);
  void  prot1_pi2_rot_func(TVector3  V3prot, TVector3 V3pi[2], int q_pi[2], double *P_1p0pi, double P_1p1pi[2]);
  void  prot1_pi3_rot_func(TVector3  V3prot, TVector3 V3pi[3], int q_pi[3], double *P_tot);
  void  prot2_pi1_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_1pi, int q_pi, TLorentzVector V4_el, double Ecal_2p1pi_to2p0pi[2],double p_miss_perp_2p1pi_to2p0pi[2],double P_2p1pito2p0pi[2],double P_2p1pito1p1pi[2],double P_2p1pito1p0pi[2],double *P_tot);
  void  prot2_pi2_rot_func(TVector3 V3_2prot_corr[2],TVector3 V3_2prot_uncorr[2],TVector3 V3_2pi[2], int q_pi[2], TLorentzVector V4_el, double Ecal_2p2pi[2],double p_miss_perp_2p2pi[2],double P_tot_2p[2]);
  void  prot3_pi1_rot_func(TVector3 V3_3prot_corr[3],TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, int q_pi, TLorentzVector V4_el, double Ecal_3p1pi[3],double p_miss_perp_3p1pi[3],double P_tot_3p[3]);

  void  pi1_rot_func(TVector3 V3_pi, int q_pi,  double *P_pi);
  void  pi2_rot_func(TVector3 V3_pi[2], int q_pi[2], double *P_0pi,double P_1pi[2]);
  void  pi3_rot_func(TVector3 V3_pi[3], int q_pi[3], double *P_0pi,double P_1pi[3],double P_320[3],double P_3210[][2]);
  void  pi4_rot_func(TVector3 V3_pi[4], int q_pi[4], double *P_0pi,double *P_410,double *P_420,double *P_4210,double *P_430,double *P_4310,double *P_4320,double *P_43210);

  void  SetQVector(TVector3 qin) {
    V3q.SetX(qin.X());
    V3q.SetY(qin.Y());
    V3q.SetZ(qin.Z());
  }

  void  ResetQVector() {
    V3q.SetX(0);
    V3q.SetY(0);
    V3q.SetZ(0);
  }

  void  PrintQVector() {
    std::cout << "Subtraction Class stored q vector: ( " << V3q.X() << " , " << V3q.Y() << " , " << V3q.Z() << " ) " << std::endl;
  }

  Bool_t EFiducialCut(std::string beam_en, TVector3 momentum) {
    return fiducialcut->EFiducialCut(beam_en, momentum);
  }
  Bool_t PFiducialCut(std::string beam_en, TVector3 momentum) {
    return fiducialcut->PFiducialCut(beam_en, momentum);
  }
  bool Pi_phot_fid_united(std::string beam_en, TVector3 V3_pi_phot, int q_pi_phot) {
    return fiducialcut->Pi_phot_fid_united(beam_en,V3_pi_phot, q_pi_phot);
  }

  Fiducial *fiducialcut;
  TVector3 V3q;
  std::string fbeam_en;
  std::string target_name;
  std::map<std::string,double> bind_en;
  int N_tot;



};

#endif
