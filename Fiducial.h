#ifndef FIDUCIAL_H
#define FIDUCIAL_H

namespace Fiducial {

  using namespace std;

  //4.4 GeV parameters
  //e- parameters
  Float_t fgPar_4Gev_2250_Efid_t0_p[6][2];// 4GeV e- fiducial cut parameters
  Float_t fgPar_4Gev_2250_Efid_t1_p[6][6];
  Float_t fgPar_4Gev_2250_Efid_b_p[6][2][6];
  Float_t fgPar_4Gev_2250_Efid_a_p[6][2][6];
  Float_t fgPar_Efid_Theta_S5_extra[8][4];
  Float_t fgPar_Efid_Theta_S4_extra[2][4];
  Float_t fgPar_Efid_Theta_S3_extra[4][4];
  Float_t fgPar_Efid_Theta_S5[8][8];
  Float_t fgPar_Efid_Theta_S4[2][8];
  Float_t fgPar_Efid_Theta_S3[4][8];
  //proton parameter
  Float_t fgPar_Pfid_ScpdS2[2][6];
  Float_t fgPar_Pfid_ScpdS3[8][6];
  Float_t fgPar_Pfid_ScpdS4[4][6];
  Float_t fgPar_Pfid_ScpdS5[8][6];
  Float_t fgPar_Pfid_ScpdS2_extra[2][4];
  Float_t fgPar_Pfid_ScpdS3_extra[8][4];
  Float_t fgPar_Pfid_ScpdS4_extra[4][4];
  Float_t fgPar_Pfid_ScpdS5_extra[8][4];
  Float_t fgPar_4Gev_2250_Pfidft1l[6][6];
  Float_t fgPar_4Gev_2250_Pfidft1r[6][6];
  Float_t fgPar_4Gev_2250_Pfidft2l[6][6];
  Float_t fgPar_4Gev_2250_Pfidft2r[6][6];
  Float_t fgPar_4Gev_2250_Pfidbt1l[6][6];
  Float_t fgPar_4Gev_2250_Pfidbt1r[6][6];
  Float_t fgPar_4Gev_2250_Pfidbt2l[6][6];
  Float_t fgPar_4Gev_2250_Pfidbt2r[6][6];
  Float_t fgPar_4Gev_2250_Pfidbl[6][6];
  Float_t fgPar_4Gev_2250_Pfidbr[6][6];

  //1.1 GeV parameters
   //electron parameters
  double fgPar_1gev_750_Efid[6][5][6];
  double fgPar_1gev_750_Efid_Theta_S3[4][8];
  double fgPar_1gev_750_Efid_Theta_S4[2][8];
  double fgPar_1gev_750_Efid_Theta_S5[8][8];
  double fgPar_1gev_1500_Efid[6][5][6];
  double fgPar_1gev_1500_Efid_Theta_S3[4][8];
  double fgPar_1gev_1500_Efid_Theta_S4[2][8];
  double fgPar_1gev_1500_Efid_Theta_S5[8][8];

  //proton parameters
  double fgPar_1gev_750_Pfid[6][5][6];
  double fgPar_1gev_750_Pfid_ScpdS2[2][6];
  double fgPar_1gev_750_Pfid_ScpdS3[8][6];
  double fgPar_1gev_750_Pfid_ScpdS4[4][6];
  double fgPar_1gev_750_Pfid_ScpdS5[8][6];
  double fgPar_1gev_1500_Pfid[6][5][6];
  double fgPar_1gev_1500_Pfid_ScpdS2[2][6];
  double fgPar_1gev_1500_Pfid_ScpdS3[8][6];
  double fgPar_1gev_1500_Pfid_ScpdS4[4][6];
  double fgPar_1gev_1500_Pfid_ScpdS5[8][6];

  //Pion Plus parameter?
  double fgPar_1gev_750_Piplfid[6][5][6];
  double fgPar_1gev_750_Piplfid_ScpdS2[2][6];
  double fgPar_1gev_750_Piplfid_ScpdS3[8][6];
  double fgPar_1gev_750_Piplfid_ScpdS4[4][6];
  double fgPar_1gev_750_Piplfid_ScpdS5[8][6];
  double fgPar_1gev_1500_Piplfid[6][5][6];
  double fgPar_1gev_1500_Piplfid_ScpdS2[2][6];
  double fgPar_1gev_1500_Piplfid_ScpdS3[8][6];
  double fgPar_1gev_1500_Piplfid_ScpdS4[4][6];
  double fgPar_1gev_1500_Piplfid_ScpdS5[8][6];

  //Pion Minus parameter?
  double fgPar_1gev_750_Pimfid[6][5][6];
  double fgPar_1gev_750_Pimfid_Theta_S3[4][8];
  double fgPar_1gev_750_Pimfid_Theta_S4[2][8];
  double fgPar_1gev_750_Pimfid_Theta_S5[8][8];
  double fgPar_1gev_750_Pimfid_Theta_S3_extra[4][4];
  double fgPar_1gev_750_Pimfid_Theta_S4_extra[2][4];
  double fgPar_1gev_750_Pimfid_Theta_S5_extra[8][4];
  double fgPar_1gev_1500_Pimfid[6][5][6];
  double fgPar_1gev_1500_Pimfid_Theta_S3[4][8];
  double fgPar_1gev_1500_Pimfid_Theta_S4[2][8];
  double fgPar_1gev_1500_Pimfid_Theta_S5[8][8];
  double fgPar_1gev_1500_Pimfid_Theta_S3_extra[4][4];
  double fgPar_1gev_1500_Pimfid_Theta_S4_extra[2][4];
  double fgPar_1gev_1500_Pimfid_Theta_S5_extra[8][4];


}
#endif
