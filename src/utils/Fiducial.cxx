#ifndef FIDUCIAL_CXX
#define FIDUCIAL_CXX

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TVectorT.h>
#include <TVector3.h>
#include <TF1.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include "utils/Fiducial.h"
#include "conf/AnalysisConstantsI.h"
#include "conf/ParticleI.h"

using namespace e4nu ; 

void Fiducial::InitPiMinusFit(const double beam_en)
{
  if (beam_en == 1.161) { myPiMinusFit = std::unique_ptr<TF1>(new TF1("myPiMinusFit","17.+4./TMath::Power(x,1.)",0,5.)); }
  if (beam_en == 2.261) { myPiMinusFit = std::unique_ptr<TF1>(new TF1("myPiMinusFit","(x<0.35)*(25.+7./TMath::Power(x,1.)) + (x>0.35)*(16.+10/TMath::Power(x,1.))",0,5.)); }
  if (beam_en == 4.461) { myPiMinusFit = std::unique_ptr<TF1>(new TF1("myPiMinusFit","(x<0.35)*(25.+7./TMath::Power(x,1.)) + (x>0.35)*(16.+10/TMath::Power(x,1.))",0,5.)); }
}

void Fiducial::InitEClimits()
{
  up_lim1_ec =std::unique_ptr<TF1>( new TF1("up_lim1_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  up_lim2_ec =std::unique_ptr<TF1>( new TF1("up_lim2_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  up_lim3_ec =std::unique_ptr<TF1>( new TF1("up_lim3_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  up_lim4_ec =std::unique_ptr<TF1>( new TF1("up_lim4_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  up_lim5_ec =std::unique_ptr<TF1>( new TF1("up_lim5_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  up_lim6_ec =std::unique_ptr<TF1>( new TF1("up_lim6_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  low_lim1_ec=std::unique_ptr<TF1>( new TF1("low_lim1_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  low_lim2_ec=std::unique_ptr<TF1>( new TF1("low_lim2_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  low_lim3_ec=std::unique_ptr<TF1>( new TF1("low_lim3_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  low_lim4_ec=std::unique_ptr<TF1>( new TF1("low_lim4_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  low_lim5_ec=std::unique_ptr<TF1>( new TF1("low_lim5_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  low_lim6_ec=std::unique_ptr<TF1>( new TF1("low_lim6_ec","[0]+(x-[1])*(x-[1])*[2]",0,360));
  rightside_lim1_ec=std::unique_ptr<TF1>( new TF1("rightside_lim1_ec","[0]*(x+[1])+[2]",0,360));
  leftside_lim1_ec=std::unique_ptr<TF1>( new TF1("leftside_lim1_ec","[0]*(x+[1])+[2]",0,360));
  rightside_lim2_ec=std::unique_ptr<TF1>( new TF1("rightside_lim2_ec","[0]*(x+[1])+[2]",0,360));
  leftside_lim2_ec=std::unique_ptr<TF1>( new TF1("leftside_lim2_ec","[0]*(x+[1])+[2]",0,360));
  rightside_lim3_ec=std::unique_ptr<TF1>( new TF1("rightside_lim3_ec","[0]*(x+[1])+[2]",0,360));
  leftside_lim3_ec=std::unique_ptr<TF1>( new TF1("leftside_lim3_ec","[0]*(x+[1])+[2]",0,360));
  rightside_lim4_ec=std::unique_ptr<TF1>( new TF1("rightside_lim4_ec","[0]*(x+[1])+[2]",0,360));
  leftside_lim4_ec=std::unique_ptr<TF1>( new TF1("leftside_lim4_ec","[0]*(x+[1])+[2]",0,360));
  rightside_lim5_ec=std::unique_ptr<TF1>( new TF1("rightside_lim5_ec","[0]*(x+[1])+[2]",0,360));
  leftside_lim5_ec=std::unique_ptr<TF1>( new TF1("leftside_lim5_ec","[0]*(x+[1])+[2]",0,360));
  rightside_lim6_ec=std::unique_ptr<TF1>( new TF1("rightside_lim6_ec","[0]*(x+[1])+[2]",0,360));
  leftside_lim6_ec=std::unique_ptr<TF1>( new TF1("leftside_lim6_ec","[0]*(x+[1])+[2]",0,360));


  up_lim1_ec->SetParameters(0.995,30,-0.0001);
  up_lim2_ec->SetParameters(0.995,90,-0.0001);
  up_lim3_ec->SetParameters(0.995,150,-0.0001);
  up_lim4_ec->SetParameters(0.995,210,-0.0001);
  up_lim5_ec->SetParameters(0.995,270,-0.0001);
  up_lim6_ec->SetParameters(0.995,330,-0.0001);
  low_lim1_ec->SetParameters(0.7,30,-0.00005);
  low_lim2_ec->SetParameters(0.7,90,-0.00005);
  low_lim3_ec->SetParameters(0.7,150,-0.00005);
  low_lim4_ec->SetParameters(0.7,210,-0.00005);
  low_lim5_ec->SetParameters(0.7,270,-0.00005);
  low_lim6_ec->SetParameters(0.7,330,-0.00005);
  leftside_lim1_ec->SetParameters(0.11,0,0.03);
  rightside_lim1_ec->SetParameters(-0.11,-60,0.03);
  leftside_lim2_ec->SetParameters(0.11,-60,0.03);
  rightside_lim2_ec->SetParameters(-0.11,-120,0.03);
  leftside_lim3_ec->SetParameters(0.11,-120,0.03);
  rightside_lim3_ec->SetParameters(-0.11,-180,0.03);
  leftside_lim4_ec->SetParameters(0.11,-180,0.03);
  rightside_lim4_ec->SetParameters(-0.11,-240,0.03);
  leftside_lim5_ec->SetParameters(0.11,-240,0.03);
  rightside_lim5_ec->SetParameters(-0.11,-300,0.03);
  leftside_lim6_ec->SetParameters(0.11,-300,0.03);
  rightside_lim6_ec->SetParameters(-0.11,-360,0.03);
}

void Fiducial::SetConstants(int in_TorusCurrent, int target_pdg, double in_en_beam) {
  fTorusCurrent = in_TorusCurrent;
  ftarget_pdg = target_pdg;
  en_beam = in_en_beam;
}

bool Fiducial::SetFiducialCutParameters(double beam_en){

  // reads from a file the parameters of the fiducial cut functions
  // Please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu(UNH)
  std::string local_dir = std::getenv("E4NUANALYSIS");
  std::string fbeam_en = std::to_string((int)(beam_en*1000));

  if( beam_en>4. && beam_en<5.){    //
    // reads FC parameters for 4.4GeV , e- and p fiducial cut parameters at 4GeV
    //
    std::ifstream param_file((local_dir+"/data/FiducialsCorrections/CLAS6/FCP_"+fbeam_en+"_"+std::to_string(fTorusCurrent)+".dat").c_str());
    std::ifstream param_file2((local_dir+"/data/FiducialsCorrections/CLAS6/PFID_"+fbeam_en+"_"+std::to_string(fTorusCurrent)+".dat").c_str());

    //	std::ifstream param_file("./FCP_4461_2250.dat");
    int param_type, sector;
    double data[6];
    while ( (sector!=6 || param_type!=21))
      {
	param_file >> param_type;
	param_file >> sector >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];
	// Test the type of parameter and assign it to the proper data array
	//  std::cout << param_type << " " << sector << std::endl;
	switch (param_type)
	  {
	  case  0:
	    for(int k=0; k<2; k++) fgPar_4Gev_2250_Efid_t0_p[sector-1][k] = data[k];
	    break;
	  case  1:
	    for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_t1_p[sector-1][k] = data[k];
	    break;
	  case 10:
	    for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[sector-1][0][k] = data[k];
	    break;
	  case 11:
	    for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[sector-1][1][k] = data[k];
	    break;
	  case 20:
	    for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[sector-1][0][k] = data[k];
	    break;
	  case 21:
	    for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[sector-1][1][k] = data[k];
	    break;
	  default:
	    printf("Error in Efid parameter file!\nReceived parameter type %d, which is not found.\nAborting!\n\n\n",param_type);
	    break;
	  }
      } // Done reading in Fiducial Region Parameters

    // ---
    for(int i = 0 ; i < 4 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_4Gev_2250_Efid_Theta_S3[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 2 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_4Gev_2250_Efid_Theta_S4[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 8 ; j++){
	param_file >> fgPar_4Gev_2250_Efid_Theta_S5[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 4 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file >> fgPar_4Gev_2250_Efid_Theta_S3_extra[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 2 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file >> fgPar_4Gev_2250_Efid_Theta_S4_extra[i][j];
      }
    }
    // ---
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file >> fgPar_4Gev_2250_Efid_Theta_S5_extra[i][j];
      }
    }
    param_file.close();



    for(int i = 0 ; i < 6 ; i++){
      for(int j = 0 ; j < 6 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfidft1l[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidft1r[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidft2l[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidft2r[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidbt1l[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidbt1r[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidbt2l[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidbt2r[i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidbl  [i][j];
	param_file2 >> fgPar_4Gev_2250_Pfidbr  [i][j];
      }
    }

    for(int i = 0 ; i < 2 ; i++){//reading the proton bad TOF cuts at 4GeV
      for(int j = 0 ; j < 6 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS2[i][j];
      }
    }
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 6 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS3[i][j];
      }
    }
    for(int i = 0 ; i < 4 ; i++){
      for(int j = 0 ; j < 6 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS4[i][j];
      }
    }
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 6 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS5[i][j];
      }
    }
    for(int i = 0 ; i < 2 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][j];
      }
    }
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][j];
      }
    }
    for(int i = 0 ; i < 4 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][j];
      }
    }
    for(int i = 0 ; i < 8 ; i++){
      for(int j = 0 ; j < 4 ; j++){
	param_file2 >> fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][j];
      }
    }
  } else if( beam_en > 1. && beam_en < 2.){
    std::ifstream param_file((local_dir+"/data/FiducialsCorrections/CLAS6/FCP_"+fbeam_en+"_"+std::to_string(fTorusCurrent)+".dat").c_str());
    std::ifstream param_file2((local_dir+"/data/FiducialsCorrections/CLAS6/PFID_"+fbeam_en+"_"+std::to_string(fTorusCurrent)+".dat").c_str());
    std::ifstream param_file3((local_dir+"/data/FiducialsCorrections/CLAS6/PIPFID_"+fbeam_en+"_"+std::to_string(fTorusCurrent)+".dat").c_str());
    std::ifstream param_file4((local_dir+"/data/FiducialsCorrections/CLAS6/PIMFID_"+fbeam_en+"_"+std::to_string(fTorusCurrent)+".dat").c_str());

    //
    // reads FC parameters for 1.1GeV , e- fiducial cut parameters at 1GeV
    //

    if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
      {

	for(Int_t sector=0;sector<6;sector++)
	  {
	    for(Int_t thetapar=0;thetapar<5;thetapar++)
	      {
		for(Int_t mompar=0;mompar<6;mompar++)
		  {
		    param_file >> fgPar_1gev_1500_Efid[sector][thetapar][mompar];
		  }
	      }
	  }
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file >> fgPar_1gev_1500_Efid_Theta_S3[i][j];
	  }
	}
	// ---
	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file >> fgPar_1gev_1500_Efid_Theta_S4[i][j];
	  }
	}
	// ---
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file >> fgPar_1gev_1500_Efid_Theta_S5[i][j];
	  }
	}

      }


    if ( fTorusCurrent < 760 && fTorusCurrent > 740){

      for(Int_t sector=0;sector<6;sector++)
	{
	  for(Int_t thetapar=0;thetapar<5;thetapar++)
	    {
	      for(Int_t mompar=0;mompar<6;mompar++)
		{
		  param_file >> fgPar_1gev_750_Efid[sector][thetapar][mompar];
		}
	    }
	}
      for(int i = 0 ; i < 4 ; i++){
	for(int j = 0 ; j < 8 ; j++){
	  param_file >> fgPar_1gev_750_Efid_Theta_S3[i][j];
	}
      }
      // ---
      for(int i = 0 ; i < 2 ; i++){
	for(int j = 0 ; j < 8 ; j++){
	  param_file >> fgPar_1gev_750_Efid_Theta_S4[i][j];
	}
      }
      // ---
      for(int i = 0 ; i < 8 ; i++){
	for(int j = 0 ; j < 8 ; j++){
	  param_file >> fgPar_1gev_750_Efid_Theta_S5[i][j];
	}
      }
    }

    param_file.close();



    //
    // reads FC parameters for 1.1GeV , p fiducial cut parameters at 1GeV
    //

    if ( fTorusCurrent< 1510 && fTorusCurrent > 1490)
      {
	for(Int_t sector=0;sector<6;sector++)
	  {
	    for(Int_t phipar=0;phipar<5;phipar++)
	      {
		for(Int_t mompar=0;mompar<6;mompar++)
		  {
		    param_file2 >> fgPar_1gev_1500_Pfid[sector][phipar][mompar];
		    //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
		    //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
		  }
	      }
	  }
	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_1500_Pfid_ScpdS2[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_1500_Pfid_ScpdS3[i][j];
	  }
	}
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_1500_Pfid_ScpdS4[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_1500_Pfid_ScpdS5[i][j];
	  }
	}
      }
    if (fTorusCurrent< 760 && fTorusCurrent > 740)
      {
	for(Int_t sector=0;sector<6;sector++)
	  {
	    for(Int_t phipar=0;phipar<5;phipar++)
	      {
		for(Int_t mompar=0;mompar<6;mompar++)
		  {
		    param_file2 >> fgPar_1gev_750_Pfid[sector][phipar][mompar];
		    //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
		    //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
		  }
	      }
	  }
	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_750_Pfid_ScpdS2[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_750_Pfid_ScpdS3[i][j];
	  }
	}
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_750_Pfid_ScpdS4[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file2 >> fgPar_1gev_750_Pfid_ScpdS5[i][j];
	  }
	}
      }


    param_file2.close();



    //reads pimi fiducial cut parameters at 1GeV

    if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
      {
	for(Int_t sector=0;sector<6;sector++)
	  {
	    for(Int_t thetapar=0;thetapar<5;thetapar++)
	      {
		for(Int_t mompar=0;mompar<6;mompar++)
		  {
		    param_file4 >> fgPar_1gev_1500_Pimfid[sector][thetapar][mompar];
		  }
	      }
	  }
	// ---
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S3[i][j];
	  }
	}
	// ---
	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S4[i][j];
	  }
	}
	// ---
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S5[i][j];
	  }
	}
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 4 ; j++){
	    param_file4>> fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][j];
	  }
	}
	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 4 ; j++){
	    param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 4 ; j++){
	    param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][j];
	  }
	}
      }
    if (fTorusCurrent< 760 && fTorusCurrent > 740)
      {
	for(Int_t sector=0;sector<6;sector++)
	  {
	    for(Int_t thetapar=0;thetapar<5;thetapar++)
	      {
		for(Int_t mompar=0;mompar<6;mompar++)
		  {
		    param_file4 >> fgPar_1gev_750_Pimfid[sector][thetapar][mompar];
		  }
	      }
	  }
	// ---
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file4 >> fgPar_1gev_750_Pimfid_Theta_S3[i][j];
	  }
	}
	// ---
	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file4 >> fgPar_1gev_750_Pimfid_Theta_S4[i][j];
	  }
	}
	// ---
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 8 ; j++){
	    param_file4 >> fgPar_1gev_750_Pimfid_Theta_S5[i][j];
	  }
	}
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 4 ; j++){
	    param_file4 >> fgPar_1gev_750_Pimfid_Theta_S3_extra[i][j];
	  }
	}
	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 4 ; j++){
	    param_file4 >> fgPar_1gev_750_Pimfid_Theta_S4_extra[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 4 ; j++){
	    param_file4 >> fgPar_1gev_750_Pimfid_Theta_S5_extra[i][j];
	    //	  std::cout << fgPar_1gev_750_Pimfid_Theta_S5_extra[i][j] << std::endl;
	  }
	}
      }
    param_file4.close();


    //reads fiducial cut parameters for pi+

    if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
      {
	for(Int_t sector=0;sector<6;sector++)
	  {
	    for(Int_t phipar=0;phipar<5;phipar++)
	      {
		for(Int_t mompar=0;mompar<6;mompar++)
		  {
		    param_file3 >> fgPar_1gev_1500_Piplfid[sector][phipar][mompar];
		    //  std::cout << "PFID " << fgPar_1gev_1500_Pfid[sector][phipar][mompar]  << std::endl;
		  }
	      }
	  }


	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS2[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS3[i][j];
	  }
	}
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS4[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS5[i][j];
	  }
	}


      }


    if ( fTorusCurrent < 760 && fTorusCurrent > 740)
      {
	for(Int_t sector=0;sector<6;sector++)
	  {
	    for(Int_t phipar=0;phipar<5;phipar++)
	      {
		for(Int_t mompar=0;mompar<6;mompar++)
		  {
		    param_file3 >> fgPar_1gev_750_Piplfid[sector][phipar][mompar];
		    //  std::cout << "PFID " << fgPar_1gev_750_Pfid[sector][phipar][mompar]  << std::endl;
		  }
	      }
	  }

	for(int i = 0 ; i < 2 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_750_Piplfid_ScpdS2[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_750_Piplfid_ScpdS3[i][j];
	  }
	}
	for(int i = 0 ; i < 4 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_750_Piplfid_ScpdS4[i][j];
	  }
	}
	for(int i = 0 ; i < 8 ; i++){
	  for(int j = 0 ; j < 6 ; j++){
	    param_file3 >> fgPar_1gev_750_Piplfid_ScpdS5[i][j];
	  }
	}
      }
    param_file3.close();
  }
  else printf("There are no fiducial cut parameters to be read at %3.1f GeV!\n", beam_en);
  return true ; 
}

Bool_t Fiducial::GetEPhiLimits(double beam_en, Float_t momentum, Float_t theta, Int_t sector,Float_t *EPhiMin, Float_t *EPhiMax){
  //Begin_Html
  /*</pre>
    Information for electron fiducial cut,
    returns the minimum and maximum phi accepted for a given momentum, theta and sector
    momentum is in GeV/c, theta is in degrees, 0 <= sector <= 5
    EPhiMin and EPhiMax are in degrees
    Function returns False if inputs are out of bounds
    1.1 GeV not implemented yet
    tested against EFiducialCut to make sure the limits are identical
    2.2 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.1 < Ef < 2.261
    2 inconsistent events out of 10^6
    4.4 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.3 < Ef < 4.461
    0 inconsistent events out of 10^6
    Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
    For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
    <pre>
  */
  //End_Html
  std::string fbeam_en = std::to_string((int)(beam_en*1000));
  if (sector < 0 || sector > 5) return kFALSE;    // bad input

  if(beam_en>4. && beam_en<5. && fTorusCurrent>2240 && fTorusCurrent<2260){// 4.4GeV fiducial cuts by protopop@jlab.org
    if ((theta < 15.) || (momentum < 0.9)) return kFALSE;         // out of range
    Float_t t0, t1, b[2], a[2];

    if (momentum > 3.7) momentum = 3.7; // don't extrapolate past the data


    // uncomment this if you want 100MeV energy bins
    //Enrgy = 0.100*int(Enrgy/0.100);


    // calculates parameters of cut functions for this energy
    t0 = fgPar_4Gev_2250_Efid_t0_p[sector][0]/pow(momentum, fgPar_4Gev_2250_Efid_t0_p[sector][1]);
    t1 = 0.; for(int k=0; k<6; k++) t1 += (fgPar_4Gev_2250_Efid_t1_p[sector][k]*pow(momentum, k));
    for(int l=0; l<2; l++){
      b[l] = 0.; for(int k=0; k<6; k++) b[l] += (fgPar_4Gev_2250_Efid_b_p[sector][l][k]*pow(momentum, k));
      a[l] = 0.; for(int k=0; k<6; k++) a[l] += (fgPar_4Gev_2250_Efid_a_p[sector][l][k]*pow(momentum, k));
    }



    // adjust upper limit according to hardware
    if(t1 < 45.) t1 = 45.;
    if(t0 < theta && theta < t1){

      *EPhiMin = 60.*sector - b[0]*(1. - 1/((theta - t0)/(b[0]/a[0]) + 1.));
      *EPhiMax = 60.*sector + b[1]*(1. - 1/((theta - t0)/(b[1]/a[1]) + 1.));
      // if(momentum<1.65 && momentum>1.60)cout<<sector<<"  "<<a[0]<<"    "<<a[1]<<"    "<<a[2]<<endl;
    }
    else {
      *EPhiMin = 60.*sector;
      *EPhiMax = 60.*sector;
    }


  }   // 4.4 GeV e2a
  else {
    return kFALSE;     // wrong beam energy/torus
  }
  return kTRUE;
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

Bool_t Fiducial::FiducialCutExtra( const int pdg, const double beam_en, TVector3 momentum, const bool is_data ){
  if( !is_data ) {
    // Electron fiducial cut, return kTRUE if pass or kFALSE if not
    momentum.SetPhi( momentum.Phi() + TMath::Pi() ) ;
  }

  // Apply first angle cuts in both MC and data
  bool extra = true ;   
  if ( pdg == conf::kPdgElectron ) {
    if( momentum.Theta() * TMath::RadToDeg() < conf::kMinThetaElectron || momentum.Theta() * TMath::RadToDeg() > conf::kMaxThetaElectron ) extra = false ; 
  }
  else if ( pdg == conf::kPdgProton ) extra *= PFiducialCutExtra( beam_en, momentum ) ;
  else if ( pdg == conf::kPdgPiP ) extra *= PiplFiducialCutExtra( beam_en, momentum ) ;
  else if ( pdg == conf::kPdgPiM ) extra *= PimiFiducialCutExtra(beam_en, momentum) ;
  else if ( pdg == conf::kPdgPhoton ) extra *= Phot_fidExtra(momentum) ;    
  
  if ( extra == false ) return false ; 
			     
  return true ;
}

  // --------------------------------------------------------------------------------------------------------------------------------------------------------
  
Bool_t Fiducial::FiducialCut( const int pdg, const double beam_en, TVector3 momentum, const bool is_data, 
			      const bool apply_fiducial ) {
  
  if( !is_data ) {
    // Electron fiducial cut, return kTRUE if pass or kFALSE if not
    momentum.SetPhi( momentum.Phi() + TMath::Pi() ) ;
  }

  if( !FiducialCutExtra( pdg, beam_en, momentum, is_data ) ) return false;
			     
  if( !apply_fiducial ) return true ; 
  
  // Apply cuts to include real detector fiducial
  if ( pdg == conf::kPdgElectron ) return EFiducialCut( beam_en, momentum ) ; 
  else if ( pdg == conf::kPdgProton ) return PFiducialCut( beam_en, momentum ) ; 
  else if ( pdg == conf::kPdgPiP ) return Pi_phot_fid_united( beam_en, momentum, 1 ) ; 
  else if ( pdg == conf::kPdgPiM ) return Pi_phot_fid_united (beam_en, momentum, -1 ) ; 
  else if ( pdg == conf::kPdgPhoton ) return Pi_phot_fid_united( beam_en, momentum, 0 ) ; 
  
  return true ; 
}


// --------------------------------------------------------------------------------------------------------------------------------------------------------

Bool_t Fiducial::EFiducialCut(double beam_en, TVector3 momentum) {

  Bool_t status = kTRUE;
  std::string fbeam_en = std::to_string((int)(beam_en*1000));
  bool SCpdcut = true;

  // 1.1 GeV electron fiducials & 750 torus field
  if ( beam_en > 1. && beam_en < 2. && fTorusCurrent > 740 && fTorusCurrent < 1510) {

    Float_t mom = momentum.Mag();
    Float_t phi = momentum.Phi() * TMath::RadToDeg() ;
    if (phi < -30.) phi += 360.;
    Float_t theta = momentum.Theta() * TMath::RadToDeg() ;
    Int_t sector = (Int_t)((phi+30.)/60.);
    if(sector < 0) sector = 0;
    if(sector > 5) sector = 5; // to match array index

    phi -= sector * 60;
    Double_t elmom = (momentum.Mag())*1000;
    Double_t thetapars[5]={0,0,0,0,0};

    for ( Int_t mompar = 0; mompar < 6; mompar++) {
      for (Int_t thetapar = 0; thetapar < 5; thetapar++) {
	if ( (fTorusCurrent > 1490) && (fTorusCurrent < 1510) ) {
	  // 1500A torus current
	  thetapars[thetapar]+=fgPar_1gev_1500_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
	}

	if ( (fTorusCurrent > 740) && (fTorusCurrent < 760) ) {
	  // 750A torus current
	  thetapars[thetapar]+=fgPar_1gev_750_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
	}
      }
    }

    Int_t uplow;
    Double_t thetacutoff;
    Float_t p_thetae = mom, thetamax_e = 0;

    if (p_thetae > 1.05) { p_thetae = 1.05; }
    else if (p_thetae < 0.4)  { p_thetae = 0.4; }

    for (int i = 4; i >= 0;i--) { thetamax_e = thetamax_e*p_thetae + el_thetamax1[i]; }

    if (phi <= 0) {
      uplow = 1;
      thetacutoff = ((phi*(thetapars[0]-(thetapars[1]/thetapars[2])))+
		     (double(uplow)*thetapars[2]*thetapars[0]))/(phi+(double(uplow)*thetapars[2]));
    } else {
      uplow = -1;
      thetacutoff = ( (phi*(thetapars[0]-(thetapars[3]/thetapars[4]))) + (double(uplow)*thetapars[4]*thetapars[0]))/(phi+(double(uplow)*thetapars[4]) );
    }

    status = (theta>thetacutoff) && (thetacutoff>=thetapars[0]) && (elmom>300) && (elmom<=1100)  && theta<=thetamax_e;

    if (SCpdcut && (fTorusCurrent>1490) && (fTorusCurrent<1510) ) {  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
        
      if (status) {

	int tsector = sector + 1;

	// sector 3 has two bad paddles

	if (tsector == 3) {

	  float badpar3[4]; // 4 parameters to determine the positions of the two theta gaps

	  for (int i = 0; i < 4; i++) {

	    badpar3[i] = 0;

	    // calculate the parameters using pol7

	    for (int d=7; d>=0; d--) { badpar3[i] = badpar3[i]*mom + fgPar_1gev_1500_Efid_Theta_S3[i][d]; }
					
	  }

	  for (int ipar = 0; ipar < 2;ipar++) { status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]); }

	}

	// sector 4 has one bad paddle

	else if (tsector == 4) {

	  float badpar4[2]; // 2 parameters to determine the position of the theta gap
		    
	  for (int i = 0; i < 2; i++) {

	    badpar4[i] = 0;

	    // calculate the parameters using pol7

	    for (int d=7; d>=0; d--) { badpar4[i] = badpar4[i]*mom + fgPar_1gev_1500_Efid_Theta_S4[i][d]; }
		    
	  }

	  status = !(theta>badpar4[0] && theta<badpar4[1]);

	}

	// sector 5 has four bad paddles

	else if (tsector == 5) {

	  Float_t badpar5[8]; // 8 parameters to determine the positions of the four theta gaps

	  for (Int_t i = 0; i < 8; i++) {

	    badpar5[i] = 0;

	    // calculate the parameters using pol7

	    for (Int_t d = 7; d >= 0; d--) { badpar5[i] = badpar5[i]*mom + fgPar_1gev_1500_Efid_Theta_S5[i][d]; }

	  }

	  if (mom < 1.25) { badpar5[0] = 23.4 * 1500 / 2250; }

	  if (mom < 1.27) { badpar5[1] = 24.0 * 1500 / 2250; } // some dummy constants. see fiducial cuts webpage.

	  for (Int_t ip = 0; ip < 4; ip++) { status = status && !(theta>badpar5[2*ip] && theta<badpar5[2*ip+1]); }

	}

      }
    }

    if (SCpdcut && (fTorusCurrent>740) && (fTorusCurrent<760) ) {  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

      if (status) {

	int tsector = sector + 1;
	mom = momentum.Mag();

	//sector 2 has one gap

	if (tsector == 2) {

	  double parsec2_l, parsec2_h;
	  if (mom < 0.4) mom = 0.4;
	  parsec2_l = fid_1gev_750_efid_S2[0][0]+fid_1gev_750_efid_S2[0][1]/mom +fid_1gev_750_efid_S2[0][2]/(mom*mom) +fid_1gev_750_efid_S2[0][3]/(mom*mom*mom);
	  parsec2_h = fid_1gev_750_efid_S2[1][0]+fid_1gev_750_efid_S2[1][1]/mom +fid_1gev_750_efid_S2[1][2]/(mom*mom) +fid_1gev_750_efid_S2[1][3]/(mom*mom*mom);
	  status=status && !(theta>parsec2_l && theta<parsec2_h);

	}

	// sector 3 has four gaps, the last two appear only at low momenta (p<0.3) and affect only pimi
		  
	if (tsector == 3) {

	  double parsec3_l[4],parsec3_h[4];
		    
	  for (int d = 0; d < 4; d++) {

	    mom = momentum.Mag();
	    if ( (d==2 || d==3) && mom>0.3 ) { mom = 0.3; }
	    else if ( d < 2 && mom < 0.4) { mom = 0.4; }
	    parsec3_l[d] = fid_1gev_750_efid_S3[d][0][0]+fid_1gev_750_efid_S3[d][0][1]/mom +fid_1gev_750_efid_S3[d][0][2]/(mom*mom) +fid_1gev_750_efid_S3[d][0][3]/(mom*mom*mom);
	    parsec3_h[d]= fid_1gev_750_efid_S3[d][1][0]+fid_1gev_750_efid_S3[d][1][1]/mom +fid_1gev_750_efid_S3[d][1][2]/(mom*mom) +fid_1gev_750_efid_S3[d][1][3]/(mom*mom*mom);
	    status=status && !(theta>parsec3_l[d] && theta<parsec3_h[d]);
		    
	  }

	}

	//sector 4 has two gaps , second gap appears only at p<0.3 and theta>105 and affects only pimi
			  
	else if (tsector == 4) {

	  double parsec4_l[2],parsec4_h[2];
			    
	  for (int d = 0; d < 2;d++) {

	    mom = momentum.Mag();
	    if( d == 0 && mom < 0.775 ) { mom = 0.775; }
	    else if ( d == 1 && mom > 0.3 ) { mom = 0.3; }
	    parsec4_l[d] = fid_1gev_750_efid_S4[d][0][0]+fid_1gev_750_efid_S4[d][0][1]/mom +fid_1gev_750_efid_S4[d][0][2]/(mom*mom) +fid_1gev_750_efid_S4[d][0][3]/(mom*mom*mom);
	    parsec4_h[d] = fid_1gev_750_efid_S4[d][1][0]+fid_1gev_750_efid_S4[d][1][1]/mom +fid_1gev_750_efid_S4[d][1][2]/(mom*mom) +fid_1gev_750_efid_S4[d][1][3]/(mom*mom*mom);
	    status=status && !(theta>parsec4_l[d] && theta<parsec4_h[d]);
			    
	  }

	}

	//sector 5 has three gaps

	else if (tsector == 5) {

	  double parsec5_l[3],parsec5_h[3];
			    
	  for (int d = 0; d < 3; d++ ) {

	    mom = momentum.Mag();
	    if ( d == 0 && mom > 0.3) { mom = 0.3; } //first one shows up only for pimi at p<0.3 and theta~128
	    else if ( d > 0 && mom < 0.5 ) { mom = 0.5; }
	    parsec5_l[d] = fid_1gev_750_efid_S5[d][0][0]+fid_1gev_750_efid_S5[d][0][1]/mom +fid_1gev_750_efid_S5[d][0][2]/(mom*mom) +fid_1gev_750_efid_S5[d][0][3]/(mom*mom*mom);
	    parsec5_h[d] = fid_1gev_750_efid_S5[d][1][0]+fid_1gev_750_efid_S5[d][1][1]/mom +fid_1gev_750_efid_S5[d][1][2]/(mom*mom) +fid_1gev_750_efid_S5[d][1][3]/(mom*mom*mom);
	    status=status && !(theta>parsec5_l[d] && theta<parsec5_h[d]);
	  }
	}
      }
    }
    return status;
  }

  // -------------------------------------------------------------------------------------------------
  // 2GeV fiducials
  if ( beam_en > 2. &&  beam_en < 3. && fTorusCurrent > 2240 && fTorusCurrent < 2260) {
    Float_t phi = momentum.Phi() * TMath::RadToDeg() ;
    if ( phi < -30. ) { phi += 360.; }
    Int_t sector = (Int_t)( (phi+30.) / 60.);
    if ( sector < 0 ) { sector = 0; }
    if ( sector > 5 ) { sector = 5; }
    phi -= sector*60;
    Float_t theta = momentum.Theta() * 180./TMath::Pi();
    Float_t mom = momentum.Mag();
    Float_t par[6];  // six parameters to determine the outline of Theta vs Phi

    for (Int_t i = 0; i < 6; i++ ) {

      par[i] = 0;

      for (Int_t d = 8; d >= 0; d--) {

	par[i] = par[i]*mom + fgPar_2GeV_2250_Efid[sector][i][d];

      } // calculate the parameters using pol8

    }

    if (phi < 0) {

      Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
      status = (theta > tmptheta && tmptheta>=par[0] && theta<par[1]);

    } else {

      Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
      status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);

    }

    // by now, we have checked if the electron is within the outline of theta vs phi plot

    if (SCpdcut) { // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap

      if (status) {

	Int_t tsector = sector + 1;

	if (tsector == 3) { // sector 3 has two bad paddles

	  Float_t badpar3[4]; // 4 parameters to determine the positions of the two theta gaps

	  for (Int_t i = 0; i < 4; i++) {

	    badpar3[i] = 0;

	    for (Int_t d = 7; d >= 0; d--) {

	      badpar3[i] = badpar3[i]*mom +  fgPar_2GeV_2250_EfidTheta_S3[i][d];

	    } // calculate the parameters using pol7

	  }

	  for (Int_t ipar = 0; ipar < 2; ipar++ ) { status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]); }

	}

	else if (tsector == 4) { // sector 4 has one bad paddle

	  Float_t badpar4[2]; // 2 parameters to determine the position of the theta gap

	  for ( Int_t i = 0; i < 2; i++) {

	    badpar4[i] = 0;

	    for (Int_t d = 7; d >= 0; d--) {

	      badpar4[i] = badpar4[i]*mom +  fgPar_2GeV_2250_EfidTheta_S4[i][d];

	    }  // calculate the parameters using pol7

	  }

	  status = !(theta>badpar4[0] && theta<badpar4[1]);

	}

	else if (tsector == 5) { // sector 5 has four bad paddles

	  Float_t badpar5[8];  // 8 parameters to determine the positions of the four theta gaps

	  for (Int_t i = 0; i < 8; i++) {

	    badpar5[i] = 0;

	    for (Int_t d = 7; d >= 0; d--) {

	      badpar5[i] = badpar5[i]*mom +  fgPar_2GeV_2250_EfidTheta_S5[i][d];

	    } // calculate the parameters using pol7

	  }

	  if (mom<1.25) badpar5[0] = 23.4;
	  if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.

	  for(Int_t ipar = 0; ipar < 4; ipar++) { status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]); }
	}
      }
    }
  }
  // ------------------------------------------------------------------------------------------------
  // 4GeV electron fiducials
  if ( beam_en > 4. && beam_en < 5. && fTorusCurrent > 2240 && fTorusCurrent < 2260) {

    //		//Begin_Html
    //		</pre>
    //		Electron fiducial cut, return kTRUE if the electron is in the fiducial volume
    //		modified 14 May 2001 lbw
    //		Now calls GetEPhiLimits for 2.2 and 4.4 GeV
    //		tested against EFiducialCut for both 2.2 (with and without bad scintillator cuts) and 4.4 GeV
    //		discrepancy less than 2 in 10^6 events
    //		Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
    //		For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
    //		Please refer to <a href="http://www.jlab.org/Hall-B/secure/e2/stevenmc/FiducialCuts/index.html">1.1 GeV fiducial cuts</a> -- Steven McLauchlan (GU).
    //		<pre>
    //		//End_Html

    Float_t phiMin, phiMax;
    Float_t mom = momentum.Mag();
    Float_t phi = momentum.Phi()*180./TMath::Pi();
    if(phi<-30.) phi += 360.;
    Float_t theta = momentum.Theta()*180./TMath::Pi();
    Int_t  sector = (Int_t)((phi+30.)/60.);
    if(sector < 0) sector = 0;
    if(sector > 5) sector = 5; // to match array index
    // all the work is now done in GetEPhiLimits

    status = GetEPhiLimits(beam_en,mom, theta, sector, &phiMin, &phiMax);

    if (status) {
      status = status && (phi > phiMin) && (phi < phiMax);
    }

    if(mom <= 2.0) {

      SCpdcut = true;

      if (SCpdcut) {  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

	if (status) {

	  int tsector = sector + 1;

	  // sector 3 has two bad paddles

	  if (tsector == 3) {

	    float badpar3[4]; // 4 parameters to determine the positions of the two theta gaps

	    for (int i = 0; i < 4; i++) {

	      badpar3[i] = 0;

	      // calculate the parameters using pol7

	      for (int d = 7; d >= 0; d-- ) { badpar3[i] = badpar3[i]*mom + fgPar_4Gev_2250_Efid_Theta_S3[i][d]; }

	    }

	    for(int ipar=0;ipar<2;ipar++) { status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]); }

	  }

	  // sector 4 has one bad paddle

	  else if ( tsector == 4) {

	    float badpar4[2]; // 2 parameters to determine the position of the theta gap

	    for (int i = 0; i < 2; i++) {

	      badpar4[i] = 0;

	      // calculate the parameters using pol7

	      for (int d = 7; d >= 0; d--) { badpar4[i] = badpar4[i]*mom + fgPar_4Gev_2250_Efid_Theta_S4[i][d]; }

	    }

	    status = !(theta>badpar4[0] && theta<badpar4[1]);

	  }

	  // sector 5 has four bad paddles

	  else if (tsector == 5) {

	    Float_t badpar5[8]; // 8 parameters to determine the positions of the four theta gaps

	    for (Int_t i = 0; i < 8; i++) {

	      badpar5[i] = 0;

	      // calculate the parameters using pol7

	      for (Int_t d = 7; d >= 0; d--) { badpar5[i] = badpar5[i]*mom + fgPar_4Gev_2250_Efid_Theta_S5[i][d];}
			
	    }

	    if (mom<1.25) badpar5[0] = 23.4;
	    if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
			
	    for (Int_t ipar = 0; ipar < 4; ipar++) { status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]); }

	  }

	}

      }

      return (status && (phi < phiMax) && (phi>phiMin));

    } else {

      SCpdcut = true;
			
      if (SCpdcut) { // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

	if (status) {

	  int tsector = sector + 1;

	  // sector 3 has two bad paddles

	  if (tsector == 3) {

	    float badpar3[4]; // 4 parameters to determine the positions of the two theta gaps

	    for (int i = 0; i < 4; i++) {

	      badpar3[i] = 0;

	      // calculate the parameters using 1/p

	      badpar3[i] = fgPar_4Gev_2250_Efid_Theta_S3_extra[i][0] + fgPar_4Gev_2250_Efid_Theta_S3_extra[i][1]/mom + fgPar_4Gev_2250_Efid_Theta_S3_extra[i][2]/(mom*mom) + fgPar_4Gev_2250_Efid_Theta_S3_extra[i][3]/(mom*mom*mom);

	    }

	    for(int ipar=0;ipar<2;ipar++) { status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]); }

	  }

	  // sector 4 has one bad paddle

	  else if (tsector == 4) {

	    float badpar4[2]; // 2 parameters to determine the position of the theta gap

	    for (int i = 0; i < 2; i++) {

	      badpar4[i] = 0;

	      // calculate the parameters using 1/p

	      badpar4[i] = fgPar_4Gev_2250_Efid_Theta_S4_extra[i][0] + fgPar_4Gev_2250_Efid_Theta_S4_extra[i][1]/mom + fgPar_4Gev_2250_Efid_Theta_S4_extra[i][2]/(mom*mom) + fgPar_4Gev_2250_Efid_Theta_S4_extra[i][3]/(mom*mom*mom);

	    }

	    status = !(theta>badpar4[0] && theta<badpar4[1]);
				
	  }

	  // sector 5 has four bad paddles

	  else if (tsector == 5) {

	    Float_t badpar5[8]; // 8 parameters to determine the positions of the four theta gaps

	    for ( Int_t i = 0; i < 8; i++) {

	      badpar5[i] = 0;

	      // calculate the parameters using 1/p

	      badpar5[i] = fgPar_4Gev_2250_Efid_Theta_S5_extra[i][0] + fgPar_4Gev_2250_Efid_Theta_S5_extra[i][1]/mom + fgPar_4Gev_2250_Efid_Theta_S5_extra[i][2]/(mom*mom) + fgPar_4Gev_2250_Efid_Theta_S5_extra[i][3]/(mom*mom*mom);

	    }

	    if (mom < 1.25) badpar5[0] = 23.4;
	    if (mom < 1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.

	    for ( Int_t ipar = 0; ipar < 4; ipar++) { status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]); }

	  }

	}

      }
			
      return (status && (phi < phiMax) && (phi>phiMin));

    }

  }

  return status;

}

// ------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------

// apapadop // Nov 11 2020 // Narrow band 30 deg in phi and either accepting ALL theta or theta_pos > 12 deg (piplus & protons) and theta_pi- > 30

double Fiducial::GetPhi(TVector3 momentum) {

  double phi = momentum.Phi() * TMath::RadToDeg() + 30.;
  if (phi < 0) { phi += 360; }
  if (phi > 360) { phi -= 360; }

  return phi;

}

// ------------------------------------------------------------------------------------------------------------------------------------------

double Fiducial::GetTheta(TVector3 momentum) {

  double theta = momentum.Theta() * TMath::RadToDeg() ;
  if (theta < 0) { theta += 180; }
  if (theta > 180) { theta -= 180; }

  return theta;

}

// ------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------------------------

Bool_t Fiducial::PFiducialCut(double beam_en, TVector3 momentum){
  //Positive Hadron Fiducial Cut
  //Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).

  if (beam_en == 0) {

    bool status = true;

    double phi = GetPhi(momentum);
    double theta = GetTheta(momentum);

    if ( !(TMath::Abs(phi - conf::kCenterFirstSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterSecondSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterThirdSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterFourthSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterFifthSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterSixthSector) < conf::kPhiOpeningAngle) ) { status = false; }

    if (theta < conf::kMinThetaProton) { status = false; }

    return status;

  } else {

    Bool_t status = kTRUE;
    std::string fbeam_en = std::to_string((int)(beam_en*1000));

    if ( beam_en>1. && beam_en<2.) {

      Float_t theta = momentum.Theta()*180/M_PI;
      Float_t phi = momentum.Phi()  *180/M_PI;
      if(phi<-30) phi+=360;
      Int_t sector = Int_t ((phi+30)/60);
      if(sector<0) sector=0;
      if(sector>5) sector=5;
      phi -= sector*60;
      Float_t p = momentum.Mag();

      if ( fTorusCurrent < 1510 && fTorusCurrent > 1490) {

	Double_t phipars[5] = {0,0,0,0,0};
	status = true;
	bool SCpdcut = true;
	if (p < .3) { return false; }
	if (p > 1) { p = 1; }

	for (Int_t mompar=0;mompar<6;mompar++) {

	  for(Int_t phipar=0;phipar<5;phipar++) {
          
	    phipars[phipar]+=fgPar_1gev_1500_Pfid[sector][phipar][mompar]*pow(p,mompar);
	    //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << std::endl;
	    //std::cout << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;

	  }

	}

	Double_t phicutoff;

	if (phi<=0) {

	  phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
	  //std::cout << "bottom " << theta << std::endl;
	  status = ((phi>phicutoff) && (theta>phipars[4]));

	} else {

	  phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
	  //std::cout << "top " << phicutoff << std::endl;
	  status = ((phi<phicutoff) && (theta>phipars[4]));

	}

	if (status && SCpdcut) { // cut bad scintillator paddles

	  Int_t tsector = sector + 1;
	  Float_t mom_scpd = p;          // Momentum for bad sc paddles cuts
	  if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
	  if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c

	  if (tsector==2) { // sector 2 has one bad paddle

	    Float_t badpar2[2]; // 2 parameters to determine the position of the theta gap

	    for (Int_t i=0; i<2; i++) {

	      badpar2[i] = 0;
  					
	      for (Int_t d=5; d>=0; d--) {

		badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS2[i][d];

	      } // calculate the parameters using pol5

	    }

	    status = status && !(theta>badpar2[0]&&theta<badpar2[1]);

	  } else if(tsector==3) { // sector 3 has four bad paddles

	    Float_t badpar3[8]; // 8 parameters to determine the positions of the theta gaps

	    for (Int_t i=0; i<8; i++){
	      badpar3[i] = 0;
	      for (Int_t d=5; d>=0; d--){
		badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS3[i][d];
	      }                // calculate the parameters using pol5
	    }
	    for (Int_t ipar=0;ipar<4;ipar++){
	      status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	    }
	  }
	  else if(tsector==4){ // sector 4 has two bad paddles
	    Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	    for (Int_t i=0; i<4; i++){
	      if (i==0 || i==1)
		if (mom_scpd > .65)
		  mom_scpd = .65;
	      badpar4[i] = 0;
	      for (Int_t d=5; d>=0; d--){
		badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS4[i][d];
	      }                // calculate the parameters using pol5
	    }
	    for (Int_t ipar=0;ipar<2;ipar++){
	      status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	    }
	  }
	  else if(tsector==5){ // sector 5 has four bad paddles
	    Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	    for (Int_t i=0; i<8; i++){
	      badpar5[i] = 0;
	      for (Int_t d=5; d>=0; d--){
		badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS5[i][d];
	      }                // calculate the parameters using pol5
	    }
	    for (Int_t ipar=0;ipar<4;ipar++){
	      status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	    }
	  }
	}

	// ---------------------------------------------------------------------------------------------

	// apapadop: Nov 23 2020, adding a lower theta bound of 12 deg for the protons
	if (!PFiducialCutExtra(beam_en, momentum)) { status = false; } 

	// ---------------------------------------------------------------------------------------------

	return status;
      }
      if (fTorusCurrent < 760 && fTorusCurrent > 740){
	Double_t phipars[5]={0,0,0,0,0};
	status = true;
	bool SCpdcut = true;
	if (p < .3)
	  return false;
	if (p > 1)
	  p = 1;
	for(Int_t mompar=0;mompar<6;mompar++) {
	  for(Int_t phipar=0;phipar<5;phipar++) {
	    phipars[phipar]+=fgPar_1gev_750_Pfid[sector][phipar][mompar]*pow(p,mompar);
	    //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
	  }
	}

	Double_t phicutoff;
	Float_t p_theta=p, thetamax_p = 0;
	if (p_theta>0.95)  p_theta = 0.95;
	else if(p_theta<0.1)   p_theta = 0.1;
	for(int i=4;i>=0;i--) thetamax_p = thetamax_p*p_theta+pipl_thetamax1[i];

	if(phi<=0) {
	  phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
	  //std::cout << "bottom " << theta << std::endl;
	  status = ((phi>phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
	}
	else {
	  phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
	  //std::cout << "top " << phicutoff << std::endl;
	  status = ((phi<phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
	}
	if(status && SCpdcut){ // cut bad scintillator paddles
	  Int_t tsector = sector + 1;
	  Float_t mom_scpd = momentum.Mag();     // momentum for bad sc paddles cuts
	  //in Marianas new code it is Float_t mom_scpd = momentum.Mag(), previoulsy it was mom_scpd = p;
	  //thus if condition of p 40 lines above is not used. Is this correct F.H. 10/31/19
	  //NOT USED if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
	  if (mom_scpd<0.15)mom_scpd=0.15; // momentum smaller than 150 MeV/c, use 150 MeV/c

	  //sector 1 has two gaps , 1st gap appears only at smaller momenta and theta>110 and affects only pipl
	  if(tsector == 1){
	    double parsec1_l[2],parsec1_h[2];
	    for(int d=0;d<2;d++){
	      //mom has to be reset for the two "d"-values in the for-loop F.H. 31/10/19
	      mom_scpd =momentum.Mag();
	      if(d==0 && mom_scpd>0.45 )mom_scpd=0.45;
	      else if(d==0 && mom_scpd<0.15 )mom_scpd=0.15;
	      else if(d==1 && mom_scpd>0.55) mom_scpd=0.55;
	      else if(d==1 && mom_scpd<0.3) mom_scpd=0.3;
	      parsec1_l[d]= fid_1gev_750_pfid_S1[d][0][0]+fid_1gev_750_pfid_S1[d][0][1]/mom_scpd +fid_1gev_750_pfid_S1[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
	      parsec1_h[d]= fid_1gev_750_pfid_S1[d][1][0]+fid_1gev_750_pfid_S1[d][1][1]/mom_scpd +fid_1gev_750_pfid_S1[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
	      status=status && !(theta>parsec1_l[d] && theta<parsec1_h[d]);
	    }
	  }

	  //sector 2 has 1 gap1
	  else if(tsector == 2){
            double parsec2_l,parsec2_h;
            parsec2_l= fid_1gev_750_pfid_S2[0][0]+fid_1gev_750_pfid_S2[0][1]/mom_scpd +fid_1gev_750_pfid_S2[0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[0][3]/(mom_scpd*mom_scpd*mom_scpd);
            parsec2_h= fid_1gev_750_pfid_S2[1][0]+fid_1gev_750_pfid_S2[1][1]/mom_scpd +fid_1gev_750_pfid_S2[1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[1][3]/(mom_scpd*mom_scpd*mom_scpd);
            status=status && !(theta>parsec2_l && theta<parsec2_h);
	  }

	  //sector 3 has four gaps
	  else if(tsector == 3){
	    double parsec3_l[4],parsec3_h[4];
	    for(int d=0;d<4;d++){
	      parsec3_l[d]= fid_1gev_750_pfid_S3[d][0][0]+fid_1gev_750_pfid_S3[d][0][1]/mom_scpd +fid_1gev_750_pfid_S3[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
	      parsec3_h[d]= fid_1gev_750_pfid_S3[d][1][0]+fid_1gev_750_pfid_S3[d][1][1]/mom_scpd +fid_1gev_750_pfid_S3[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
	      status=status && !(theta>parsec3_l[d] && theta<parsec3_h[d]);
	    }
	  }

	  //sector 4 has two gaps
	  else  if(tsector == 4){
	    double parsec4_l[2],parsec4_h[2];
	    for(int d=0;d<2;d++){
	      parsec4_l[d]= fid_1gev_750_pfid_S4[d][0][0]+fid_1gev_750_pfid_S4[d][0][1]/mom_scpd +fid_1gev_750_pfid_S4[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
	      parsec4_h[d]= fid_1gev_750_pfid_S4[d][1][0]+fid_1gev_750_pfid_S4[d][1][1]/mom_scpd +fid_1gev_750_pfid_S4[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
	      status=status && !(theta>parsec4_l[d] && theta<parsec4_h[d]);
	    }
	  }

	  //sector 5 has four gaps
	  else if(tsector == 5){//the fourth bad TOF pd. can be seen only below b=0.3 and so there are just three bad TOFs for p
	    double parsec5_l[3],parsec5_h[3];
	    for(int d=0;d<3;d++){
	      //  if(d==0 && d==1 && mom_scpd>0.6)mom_scpd=0.6;
	      mom_scpd = momentum.Mag();
	      if(d==2 && mom_scpd<0.5)mom_scpd=0.5;
	      parsec5_l[d]= fid_1gev_750_pfid_S5[d][0][0]+fid_1gev_750_pfid_S5[d][0][1]/mom_scpd +fid_1gev_750_pfid_S5[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
	      parsec5_h[d]= fid_1gev_750_pfid_S5[d][1][0]+fid_1gev_750_pfid_S5[d][1][1]/mom_scpd +fid_1gev_750_pfid_S5[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
	      status=status && !(theta>parsec5_l[d] && theta<parsec5_h[d]);
	    }
	  }

	  //sector 6 has two gaps
	  else if(tsector == 6){
	    double parsec6_l[2],parsec6_h[2];
	    for(int d=0;d<2;d++){
	      mom_scpd = momentum.Mag();
	      if(mom_scpd>0.6 )mom_scpd=0.6;
	      else if(mom_scpd<0.3)mom_scpd=0.3;
	      parsec6_l[d]= fid_1gev_750_pfid_S6[d][0][0]+fid_1gev_750_pfid_S6[d][0][1]/mom_scpd +fid_1gev_750_pfid_S6[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
	      parsec6_h[d]= fid_1gev_750_pfid_S6[d][1][0]+fid_1gev_750_pfid_S6[d][1][1]/mom_scpd +fid_1gev_750_pfid_S6[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
	      status=status && !(theta>parsec6_l[d] && theta<parsec6_h[d]);
	    }
	  }

	}


	// ---------------------------------------------------------------------------------------------

	// apapadop: Nov 23 2020, adding a lower theta bound of 12 deg for the protons
	if (!PFiducialCutExtra(beam_en, momentum)) { status = false; } 

	// ---------------------------------------------------------------------------------------------

	return status;
      }

    }


    if ( beam_en>2. && beam_en<3. && fTorusCurrent>2240 && fTorusCurrent<2260){
      bool SCpdcut = true;
      Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
      Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum.Theta()*180/TMath::Pi();
      Float_t p = momentum.Mag();
      Float_t mom_for = p;              // momentum for forward constraints
      if (mom_for<0.3) mom_for = 0.3;   // momentum smaller than 300 MeV/c, use 300 MeV/c
      if (mom_for>1.6) mom_for = 1.6;   // momentum greater than 1.6 GeV/c, use 1.6 GeV/c
      Float_t mom_bak = p;              // momentum for backward constraints
      if (mom_bak<0.2) mom_bak = 0.2;   // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_bak>1.0) mom_bak = 1.0;   // momentum greater than 1.0 GeV/c, use 1.0 GeV/c
      Float_t theta0 = 8.5;
      Float_t phi_lower = -24.0;
      Float_t phi_upper = 24.0;
      Float_t par_for[4], par_bak[4];
      for (Int_t i=0; i<4; i++){
        par_for[i] = 0; par_bak[i] = 0;
        for (Int_t d=6; d>=0; d--){
          par_for[i] = par_for[i]*mom_for +  fgPar_2GeV_2250_Pfid_For[sector][i][d];
          par_bak[i] = par_bak[i]*mom_bak +  fgPar_2GeV_2250_Pfid_Bak[sector][i][d];
        }
      }
      if (phi < 0) {
        Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
        status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
      }
      else {
        Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
        status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);
      }                     // now the forward constrains are checked
      if ( status ) {       // now check the backward constrains
        if(theta>par_bak[0]) status = kFALSE;
        else if(theta>par_bak[1]) status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]);
      }

      if(status && SCpdcut){ // cut bad scintillator paddles

        Int_t tsector = sector + 1;
        Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
        if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
        if(tsector==2){      // sector 2 has one bad paddle
	  Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
	  for (Int_t i=0; i<2; i++){
	    badpar2[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar2[i] = badpar2[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS2[i][d];
	    }                // calculate the parameters using pol5
	  }
	  status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
	}
	else if(tsector==3){ // sector 3 has four bad paddles
	  Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar3[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar3[i] = badpar3[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS3[i][d];
	    }                // calculate the parameters using pol5
	  }
	  for (Int_t ipar=0;ipar<4;ipar++){
	    status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
	}
	else if(tsector==4){ // sector 4 has two bad paddles
	  Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	  for (Int_t i=0; i<4; i++){
	    badpar4[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar4[i] = badpar4[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS4[i][d];
	    }                // calculate the parameters using pol5
	  }
	  for (Int_t ipar=0;ipar<2;ipar++){
	    status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	  }
	}
	else if(tsector==5){ // sector 5 has four bad paddles
	  Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar5[i] = badpar5[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS5[i][d];
	    }                // calculate the parameters using pol5
	  }
	  for (Int_t ipar=0;ipar<4;ipar++){
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	  }
	}
      }
    }

    if (beam_en>4. && beam_en<5. && fTorusCurrent>2240 && fTorusCurrent<2260){//4 GeV Fiducial Cut Rustam Niyazov

      Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
      Int_t sector = Int_t ((phi+30)/60); if(sector<0)sector=0; if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum.Theta()*180/TMath::Pi();
      Float_t p = momentum.Mag();

      Float_t parfidl[3];for(Int_t i=0; i<3; i++){parfidl[i]=0;}
      Float_t parfidr[3];for(Int_t i=0; i<3; i++){parfidr[i]=0;}
      Float_t parfidbl[2];for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
      Float_t parfidbr[2];for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
      Float_t cphil=0;Float_t cphir=0;
      Float_t phi45l=0; Float_t phi45r=0;
      Float_t phi60l=0; Float_t phi60r=0;
      Float_t theta_min=11;

      bool Forward=kFALSE; //defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
      Int_t thetab=45; //this variable defines the edge point for Forward<->Backward regions
      Float_t p1=0.575; //last bin momentum for region p<0.6 GeV/c
      Float_t theta_max=140;
      if(p<0.2)p=0.2; //momentum less than 0.2 GeV/c, use 0.2 GeV/c
      if(p>4.4)p=4.4; //momentum greater than 4.4 GeV/c, use 4.4 GeV/c

      //get parametrized values of theta_max for p<0.6 GeV/c region
      if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
      //get parametrized values of theta_max for p>0.6 GeV/c region
      else{theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p1;}

      //Get the momentum dependent parameters for Forward Region (theta <45 deg)
      Forward=kTRUE;
      if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
        //parameters for hyperbolic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;
        }
      }
      else{//forward2 defines  regions of momenta and p>0.6 GeV/c
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;
        }
      }
      phi45l=parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
      phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
      if(theta>thetab){//backward region defined by theta >45 deg.
        if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
        if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

        //Get the momentum dependent parameters for Backward Region

        Forward=kFALSE;
        if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
          //parameters for quadratic function
          for (Int_t i=0; i<3; i++){
            Int_t j=2*i;
            parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
            parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;
          }
          //these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
          for (Int_t i=0; i<2; i++){
            Int_t j=2*i;
            parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
            parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;
          }
        }
        else{//backward2 defines  regions of momenta p>0.6 GeV/c
          //parameters for quadratic function
          for (Int_t i=0; i<3; i++){
            Int_t j=2*i;
            parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
            parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;
          }
          //these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
          for (Int_t i=0; i<2; i++){
            Int_t j=2*i;
            parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
            parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;
          }
        }
      }

      if(Forward){//Forward region
        if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.
	cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
	cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
      }
      else{//Backward region
        phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
        phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

        if(theta<60){
          cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
          cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
        }
        Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0;
        //dr and er are theta_flat and phi_edge parameters for phi>0;
        dl=parfidbl[0];el=parfidbl[1];
        dr=parfidbr[0];er=parfidbr[1];

        if(theta>45&&theta<60){ //BackwardA region
          //try to match parametrized values from Forward region to Backward region parameters
          if(cphil>phi45l)cphil=phi45l;
          if(cphir<phi45r)cphir=phi45r;
        }
        //BackwardB region & phi<0
        else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant
        else if(theta>dl&&theta<=theta_max){
          cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line
        else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
        //BackwardB region & phi>0
        if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant
        else if(theta>dr&&theta<=theta_max){
          cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line
        else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
      }//Backward Region

      if(phi<0) status=(phi>cphil); //check the constrains
      else if(phi>=0) {status=(phi<cphir);
      }

      if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min

      if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

      //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for
      //some range of momentum, where edge does not look good.
      bool s1s4=(theta<11.7&&(sector==0||sector==3));
      bool s5=(theta<12.2&&sector==4);
      bool s6=(theta<11.4&&sector==5);
      if(p>=0.6&&p<1.5&&(s1s4||s5||s6)) status=kFALSE;



      bool SCpdcut = true;

      if(status && SCpdcut){ // cut bad scintillator paddles
        if(p < 1.0){
          Int_t tsector = sector + 1;
          Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
          if (mom_scpd<0.3)mom_scpd=0.3; // momentum smaller than 200 MeV/c, use 200 MeV/c
          if(tsector==2){      // sector 2 has one bad paddle
            Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
            for (Int_t i=0; i<2; i++){
              badpar2[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar2[i] = badpar2[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2[i][d];
              }                // calculate the parameters using pol5
            }
            status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
          }
          else if(tsector==3){ // sector 3 has four bad paddles
            Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
            for (Int_t i=0; i<8; i++){
              badpar3[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar3[i] = badpar3[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3[i][d];
              }                // calculate the parameters using pol5
            }
            for (Int_t ipar=0;ipar<4;ipar++){
              status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
            }
          }
          else if(tsector==4){ // sector 4 has two bad paddles
            Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
            for (Int_t i=0; i<4; i++){
              badpar4[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar4[i] = badpar4[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4[i][d];
              }                // calculate the parameters using pol5
            }
            for (Int_t ipar=0;ipar<2;ipar++){
              status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
            }
          }
          else if(tsector==5){ // sector 5 has four bad paddles
            Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
            for (Int_t i=0; i<8; i++){
              badpar5[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar5[i] = badpar5[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5[i][d];
              }                // calculate the parameters using pol5
            }
            for (Int_t ipar=0;ipar<4;ipar++){
              status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
            }
          }
        }
        else{
          int tsector = sector + 1;
          double mom_scpd =p;
          // sector 2 has one bad paddles
          if (tsector == 2){
            float badpar2[2];            // 4 parameters to determine the positions of the two theta gaps
            for (int i=0; i<2; i++){
              badpar2[i] = 0;
              // calculate the parameters using 1/p
              badpar2[i] = fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(int ipar=0;ipar<1;ipar++)
	      status = status && !(theta>badpar2[2*ipar] && theta<badpar2[2*ipar+1]);
	    // status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
          if (tsector == 3){
            float badpar3[8];            // 4 parameters to determine the positions of the two theta gaps
            for (int i=0; i<8; i++){
              badpar3[i] = 0;
              // calculate the parameters using 1/p
              badpar3[i] = fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(int ipar=0;ipar<4;ipar++)
              status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
          // sector 4 has two bad paddle
          else if (tsector == 4){
            float badpar4[4];     // 2 parameters to determine the position of the theta gap
            for (int i=0; i<4; i++){
              badpar4[i] = 0;
              // calculate the parameters using 1/p
              badpar4[i] = fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(int ipar=0;ipar<2;ipar++)
	      status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
          // sector 5 has four bad paddles
          else if (tsector == 5){
            Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
            for (Int_t i=0; i<8; i++){
              badpar5[i] = 0;
              // calculate the parameters using 1/p
              badpar5[i] = fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(Int_t ipar=0;ipar<4;ipar++)
              status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
          }
        }
      }



    }

    // ---------------------------------------------------------------------------------------------

    // apapadop: Nov 23 2020, adding a lower theta bound of 12 deg for the protons
    if (!PFiducialCutExtra(beam_en, momentum)) { status = false; } 

    // ---------------------------------------------------------------------------------------------

    return status;

  } // end of the if statement for non-"" beam energy

}

// ---------------------------------------------------------------------------------------------

Bool_t Fiducial::PiplFiducialCut(double beam_en, TVector3 momentum, Float_t *philow, Float_t *phiup){
  //Positive Hadron Fiducial Cut
  //Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).

  if (beam_en == 0) {

    bool status = true;

    double phi = GetPhi(momentum);
    double theta = GetTheta(momentum);

    if ( !(TMath::Abs(phi - conf::kCenterFirstSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterSecondSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterThirdSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterFourthSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterFifthSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterSixthSector) < conf::kPhiOpeningAngle) ) { status = false; }

    if (theta < conf::kMinThetaPiPlus) { status = false; }

    return status;

  } else {
    Bool_t status = kTRUE;

    if(beam_en>1. && beam_en<2.){

      Float_t theta = momentum.Theta()*180/M_PI;
      Float_t phi   = momentum.Phi()  *180/M_PI;
      if(phi<-30) phi+=360;
      Int_t sector = Int_t ((phi+30)/60);
      if(sector<0) sector=0;
      if(sector>5) sector=5;
      phi -= sector*60;
      Float_t p = momentum.Mag();


      if (fTorusCurrent < 1510 && fTorusCurrent > 1490){
	Double_t phipars[5]={0,0,0,0,0};
	status = true;
	bool SCpdcut = true;
	if (p < 0.15)p=0.15;
	if (p > 1)
	  p = 1;
	for(Int_t mompar=0;mompar<6;mompar++) {
	  for(Int_t phipar=0;phipar<5;phipar++) {
	    phipars[phipar]+=fgPar_1gev_1500_Pfid[sector][phipar][mompar]*pow(p,mompar);
	    //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
	  }
	}

	Double_t phicutoff;
	if(phi<=0) {
	  phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
	  //std::cout << "bottom " << theta << std::endl;
	  status = ((phi>phicutoff) && (theta>phipars[4]));
	}
	else {
	  phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
	  //std::cout << "top " << phicutoff << std::endl;
	  status = ((phi<phicutoff) && (theta>phipars[4]));
	}
	if(status && SCpdcut){ // cut bad scintillator paddles
	  Int_t tsector = sector + 1;
	  Float_t mom_scpd = p;          // Momentum for bad sc paddles cuts
	  if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
	  if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
	  if(tsector==2){      // sector 2 has one bad paddle
	    Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
	    for (Int_t i=0; i<2; i++){
	      badpar2[i] = 0;
	      for (Int_t d=5; d>=0; d--){
		badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS2[i][d];
	      }                // calculate the parameters using pol5
	    }
	    status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
	  }
	  else if(tsector==3){ // sector 3 has four bad paddles
	    Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	    for (Int_t i=0; i<8; i++){
	      badpar3[i] = 0;
	      for (Int_t d=5; d>=0; d--){
		badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS3[i][d];
	      }                // calculate the parameters using pol5
	    }
	    for (Int_t ipar=0;ipar<4;ipar++){
	      status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	    }
	  }
	  else if(tsector==4){ // sector 4 has two bad paddles
	    Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	    for (Int_t i=0; i<4; i++){
	      if (i==0 || i==1)
		if (mom_scpd > .65)
		  mom_scpd = .65;
	      badpar4[i] = 0;
	      for (Int_t d=5; d>=0; d--){
		badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS4[i][d];
	      }                // calculate the parameters using pol5
	    }
	    for (Int_t ipar=0;ipar<2;ipar++){
	      status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	    }
	  }
	  else if(tsector==5){ // sector 5 has four bad paddles
	    Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	    for (Int_t i=0; i<8; i++){
	      badpar5[i] = 0;
	      for (Int_t d=5; d>=0; d--){
		badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS5[i][d];
	      }                // calculate the parameters using pol5
	    }
	    for (Int_t ipar=0;ipar<4;ipar++){
	      status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	    }
	  }
	}

	// ---------------------------------------------------------------------------------------------

	// apapadop: Nov 23 2020, adding a lower theta bound of 12 deg for the pi plus
	if (!PiplFiducialCutExtra(beam_en, momentum)) { status = false; } 

	// ---------------------------------------------------------------------------------------------

	return status;
      }
      if (fTorusCurrent < 760 && fTorusCurrent > 740){
	Double_t phipars[5]={0,0,0,0,0};
	status = true;
	bool SCpdcut = true;
	if (p < 0.15)p=0.15;
	if (p > 1)
	  p = 1;
	for(Int_t mompar=0;mompar<6;mompar++) {
	  for(Int_t phipar=0;phipar<5;phipar++) {
	    phipars[phipar]+=fgPar_1gev_750_Pfid[sector][phipar][mompar]*pow(p,mompar);
	    //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
	  }
	}

	Double_t phicutoff;
	Float_t p_theta=p, thetamax_p = 0;
	if (p_theta>0.95)  p_theta = 0.95;
	else if(p_theta<0.1)   p_theta = 0.1;
	for(int i=4;i>=0;i--)thetamax_p = thetamax_p*p_theta+pipl_thetamax1[i];

	if(phi<=0) {
	  phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
	  //std::cout << "bottom " << theta << std::endl;
	  status = ((phi>phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
	}
	else {
	  phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
	  //std::cout << "top " << phicutoff << std::endl;
	  status = ((phi<phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
	}
	if(status && SCpdcut){ // cut bad scintillator paddles
	  Int_t tsector = sector + 1;
	  Float_t mom_scpd = momentum.Mag();          // momentum for bad sc paddles cuts
	  //F.H. 10/31/19 Mariana's update mom_scpd =momentum.Mag() compared to mom_scpd= p before; but this skips if conditions above on "p"
	  //NOT USED by Mariana's update  if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
	  if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c

	  //sector 1 has two gaps , 1st gap appears only at smaller momenta and theta>110 and affects only pipl
	  if(tsector == 1){
	    double parsec1_l[2],parsec1_h[2];
	    for(int d=0;d<2;d++){

	      mom_scpd =momentum.Mag();
	      if(d==0 && mom_scpd>0.45 )mom_scpd=0.45;
	      else if(d==0 && mom_scpd<0.15 )mom_scpd=0.15;
	      else if(d==1 && mom_scpd>0.55) mom_scpd=0.55;
	      else if(d==1 && mom_scpd<0.3) mom_scpd=0.3;
	      parsec1_l[d]= fid_1gev_750_pfid_S1[d][0][0]+fid_1gev_750_pfid_S1[d][0][1]/mom_scpd +fid_1gev_750_pfid_S1[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
	      parsec1_h[d]= fid_1gev_750_pfid_S1[d][1][0]+fid_1gev_750_pfid_S1[d][1][1]/mom_scpd +fid_1gev_750_pfid_S1[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
	      status=status && !(theta>parsec1_l[d] && theta<parsec1_h[d]);
	    }
	  }
	  //sector 2 has 1 gap
	  else  if(tsector == 2){
	    double parsec2_l,parsec2_h;
	    parsec2_l= fid_1gev_750_pfid_S2[0][0]+fid_1gev_750_pfid_S2[0][1]/mom_scpd +fid_1gev_750_pfid_S2[0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[0][3]/(mom_scpd*mom_scpd*mom_scpd);
	    parsec2_h= fid_1gev_750_pfid_S2[1][0]+fid_1gev_750_pfid_S2[1][1]/mom_scpd +fid_1gev_750_pfid_S2[1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[1][3]/(mom_scpd*mom_scpd*mom_scpd);
	    status=status && !(theta>parsec2_l && theta<parsec2_h);
	  }
	  //sector 3 has four gaps
	  else if(tsector == 3){
            double parsec3_l[4],parsec3_h[4];
            for(int d=0;d<4;d++){
              parsec3_l[d]= fid_1gev_750_pfid_S3[d][0][0]+fid_1gev_750_pfid_S3[d][0][1]/mom_scpd +fid_1gev_750_pfid_S3[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
              parsec3_h[d]= fid_1gev_750_pfid_S3[d][1][0]+fid_1gev_750_pfid_S3[d][1][1]/mom_scpd +fid_1gev_750_pfid_S3[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
              status=status && !(theta>parsec3_l[d] && theta<parsec3_h[d]);
            }
	  }
	  //sector 4 has two gaps
	  else if(tsector == 4){
            double parsec4_l[2],parsec4_h[2];
            for(int d=0;d<2;d++){
              parsec4_l[d]= fid_1gev_750_pfid_S4[d][0][0]+fid_1gev_750_pfid_S4[d][0][1]/mom_scpd +fid_1gev_750_pfid_S4[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
              parsec4_h[d]= fid_1gev_750_pfid_S4[d][1][0]+fid_1gev_750_pfid_S4[d][1][1]/mom_scpd +fid_1gev_750_pfid_S4[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
              status=status && !(theta>parsec4_l[d] && theta<parsec4_h[d]);
            }
	  }
	  //sector 5 has four gaps
	  else if(tsector == 5){
            double parsec5_l[4],parsec5_h[4];
            for(int d=0;d<4;d++){
              mom_scpd=momentum.Mag();
              //  if(d==0 && d==1 && mom_scpd>0.6)mom_scpd=0.6;
              if(d==2 && mom_scpd<0.5)mom_scpd=0.5;
              if(d==3 && mom_scpd>0.3)mom_scpd=0.3;
              parsec5_l[d]= fid_1gev_750_pfid_S5[d][0][0]+fid_1gev_750_pfid_S5[d][0][1]/mom_scpd +fid_1gev_750_pfid_S5[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
              parsec5_h[d]= fid_1gev_750_pfid_S5[d][1][0]+fid_1gev_750_pfid_S5[d][1][1]/mom_scpd +fid_1gev_750_pfid_S5[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
              status=status && !(theta>parsec5_l[d] && theta<parsec5_h[d]);
            }
	  }
	  //sector 6 has two gaps
	  else if(tsector == 6){
            double parsec6_l[2],parsec6_h[2];
            for(int d=0;d<2;d++){
	      mom_scpd = momentum.Mag();
	      if(mom_scpd>0.6 )mom_scpd=0.6;
	      else if(mom_scpd<0.3)mom_scpd=0.3;
	      parsec6_l[d]= fid_1gev_750_pfid_S6[d][0][0]+fid_1gev_750_pfid_S6[d][0][1]/mom_scpd +fid_1gev_750_pfid_S6[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
	      parsec6_h[d]= fid_1gev_750_pfid_S6[d][1][0]+fid_1gev_750_pfid_S6[d][1][1]/mom_scpd +fid_1gev_750_pfid_S6[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
	      status=status && !(theta>parsec6_l[d] && theta<parsec6_h[d]);
            }
	  }
	}

	// ---------------------------------------------------------------------------------------------

	// apapadop: Nov 23 2020, adding a lower theta bound of 12 deg for the pi plus
	if (!PiplFiducialCutExtra(beam_en, momentum)) { status = false; } 

	// ---------------------------------------------------------------------------------------------

	return status;
      }
    }

    if (beam_en>2. && beam_en<3. && fTorusCurrent>2240 && fTorusCurrent<2260){
      bool SCpdcut = true;
      Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
      Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum.Theta()*180/TMath::Pi();
      Float_t p = momentum.Mag();
      Float_t mom_for = p;              // momentum for forward constraints
      if (mom_for<0.3) mom_for = 0.3;   // momentum smaller than 300 MeV/c, use 300 MeV/c
      if (mom_for>1.6) mom_for = 1.6;   // momentum greater than 1.6 GeV/c, use 1.6 GeV/c
      Float_t mom_bak = p;              // momentum for backward constraints
      if (mom_bak<0.2) mom_bak = 0.2;   // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_bak>1.0) mom_bak = 1.0;   // momentum greater than 1.0 GeV/c, use 1.0 GeV/c
      Float_t theta0 = 8.5;
      Float_t phi_lower = -24.0;
      Float_t phi_upper = 24.0;
      Float_t phimin, phimax;
      Float_t par_for[4], par_bak[4];
      for (Int_t i=0; i<4; i++){
        par_for[i] = 0; par_bak[i] = 0;
        for (Int_t d=6; d>=0; d--){
          par_for[i] = par_for[i]*mom_for +  fgPar_2GeV_2250_Pfid_For[sector][i][d];
          par_bak[i] = par_bak[i]*mom_bak +  fgPar_2GeV_2250_Pfid_Bak[sector][i][d];
        }
      }
      if (phi < 0) {
        Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
	phimin = par_for[1]/((theta-theta0)+par_for[1]/par_for[0])-par_for[0];
	phimax = par_for[0]-par_for[1]/((theta-theta0)+par_for[1]/par_for[0]);
	*philow = phimin;
	*phiup = phimax;
        status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
      }
      else {
        Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
        phimin = par_for[3]/(theta-theta0+par_for[3]/par_for[2])-par_for[2];
        phimax = par_for[2]-par_for[3]/(theta-theta0+par_for[3]/par_for[2]);
        *phiup = phimax;
        *philow = phimin;
        status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);
      }                     // now the forward constrains are checked
      if ( status ) {       // now check the backward constrains
        if(theta>par_bak[0]) status = kFALSE;
        else if(theta>par_bak[1]) status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]);
      }

      if(status && SCpdcut){ // cut bad scintillator paddles

        Int_t tsector = sector + 1;
        Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
        //if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c //pion threshold is 150 mev, gaps have been now been updated to be valid to lower momentum
        if(tsector==2){      // sector 2 has one bad paddle
          Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
          for (Int_t i=0; i<2; i++){
            badpar2[i] = 0;
            for (Int_t d=5; d>=0; d--){
	      badpar2[i] = badpar2[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS2[i][d];
            }                // calculate the parameters using pol5
          }
          status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
        }
        else if(tsector==3){ // sector 3 has four bad paddles
	  Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar3[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar3[i] = badpar3[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS3[i][d];
	    }                // calculate the parameters using pol5
	  }
	  for (Int_t ipar=0;ipar<4;ipar++){
	    status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	  }
        }
        else if(tsector==4){ // sector 4 has two bad paddles
	  Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	  for (Int_t i=0; i<4; i++){
	    badpar4[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar4[i] = badpar4[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS4[i][d];
	    }                // calculate the parameters using pol5
	  }
	  for (Int_t ipar=0;ipar<2;ipar++){
	    status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	  }
        }
        else if(tsector==5){ // sector 5 has four bad paddles
	  Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar5[i] = badpar5[i]*mom_scpd +  fgPar_2GeV_2250_Pfid_ScpdS5[i][d];
	    }                // calculate the parameters using pol5
	  }
	  for (Int_t ipar=0;ipar<4;ipar++){
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	  }
        }
      }

    }

    if (beam_en>4. && beam_en<5. && fTorusCurrent>2240 && fTorusCurrent<2260){//4 GeV Fiducial Cut Rustam Niyazov

      Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
      Int_t sector = Int_t ((phi+30)/60); if(sector<0)sector=0; if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum.Theta()*180/TMath::Pi();
      Float_t p = momentum.Mag();

      Float_t parfidl[3];for(Int_t i=0; i<3; i++){parfidl[i]=0;}
      Float_t parfidr[3];for(Int_t i=0; i<3; i++){parfidr[i]=0;}
      Float_t parfidbl[2];for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
      Float_t parfidbr[2];for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
      Float_t cphil=0;Float_t cphir=0;
      Float_t phi45l=0; Float_t phi45r=0;
      Float_t phi60l=0; Float_t phi60r=0;
      Float_t theta_min=11;

      bool Forward=kFALSE; //defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
      Int_t thetab=45; //this variable defines the edge point for Forward<->Backward regions
      Float_t p1=0.575; //last bin momentum for region p<0.6 GeV/c
      Float_t theta_max=140;
      if(p<0.2)p=0.2; //momentum less than 0.2 GeV/c, use 0.2 GeV/c
      if(p>4.4)p=4.4; //momentum greater than 4.4 GeV/c, use 4.4 GeV/c

      //get parametrized values of theta_max for p<0.6 GeV/c region
      if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
      //get parametrized values of theta_max for p>0.6 GeV/c region
      else{theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p1;}

      //Get the momentum dependent parameters for Forward Region (theta <45 deg)
      Forward=kTRUE;
      if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
        //parameters for hyperbolic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;
        }
      }
      else{//forward2 defines  regions of momenta and p>0.6 GeV/c
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;
        }
      }
      phi45l=parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
      phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
      if(theta>thetab){//backward region defined by theta >45 deg.
        if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
        if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

        //Get the momentum dependent parameters for Backward Region

        Forward=kFALSE;
        if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
          //parameters for quadratic function
          for (Int_t i=0; i<3; i++){
            Int_t j=2*i;
            parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
            parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;
          }
          //these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
          for (Int_t i=0; i<2; i++){
            Int_t j=2*i;
            parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
            parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;
          }
        }
        else{//backward2 defines  regions of momenta p>0.6 GeV/c
          //parameters for quadratic function
          for (Int_t i=0; i<3; i++){
            Int_t j=2*i;
            parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
            parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;
          }
          //these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
          for (Int_t i=0; i<2; i++){
            Int_t j=2*i;
            parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
            parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;
          }
        }
      }

      if(Forward){//Forward region
        if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.
	cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
	cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
      }
      else{//Backward region
        phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
        phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

        if(theta<60){
          cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
          cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
        }
        Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0;
        //dr and er are theta_flat and phi_edge parameters for phi>0;
        dl=parfidbl[0];el=parfidbl[1];
        dr=parfidbr[0];er=parfidbr[1];

        if(theta>45&&theta<60){ //BackwardA region
          //try to match parametrized values from Forward region to Backward region parameters
          if(cphil>phi45l)cphil=phi45l;
          if(cphir<phi45r)cphir=phi45r;
        }
        //BackwardB region & phi<0
        else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant
        else if(theta>dl&&theta<=theta_max){
          cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line
        else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
        //BackwardB region & phi>0
        if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant
        else if(theta>dr&&theta<=theta_max){
          cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line
        else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
      }//Backward Region


      if(phi<0) status=(phi>cphil); //check the constrains
      else if(phi>=0) {status=(phi<cphir);
      }

      if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min

      if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

      //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for
      //some range of momentum, where edge does not look good.
      bool s1s4=(theta<11.7&&(sector==0||sector==3));
      bool s5=(theta<12.2&&sector==4);
      bool s6=(theta<11.4&&sector==5);
      if(p>=0.6&&p<1.5&&(s1s4||s5||s6)) status=kFALSE;

      *philow = cphil;
      *phiup = cphir;



      bool SCpdcut = true;
      if(status && SCpdcut){ // cut bad scintillator paddles
        if(p < 1.0){
          Int_t tsector = sector + 1;
          Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
          if (mom_scpd<0.3)mom_scpd=0.3; // momentum smaller than 200 MeV/c, use 200 MeV/c
          if(tsector==2){      // sector 2 has one bad paddle
            Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
            for (Int_t i=0; i<2; i++){
              badpar2[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar2[i] = badpar2[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2[i][d];
              }                // calculate the parameters using pol5
            }
            status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
          }
          else if(tsector==3){ // sector 3 has four bad paddles
            Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
            for (Int_t i=0; i<8; i++){
              badpar3[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar3[i] = badpar3[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3[i][d];
              }                // calculate the parameters using pol5
            }
            for (Int_t ipar=0;ipar<4;ipar++){
              status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
            }
          }
          else if(tsector==4){ // sector 4 has two bad paddles
            Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
            for (Int_t i=0; i<4; i++){
              badpar4[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar4[i] = badpar4[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4[i][d];
              }                // calculate the parameters using pol5
            }
            for (Int_t ipar=0;ipar<2;ipar++){
              status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
            }
          }
          else if(tsector==5){ // sector 5 has four bad paddles
            Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
            for (Int_t i=0; i<8; i++){
              badpar5[i] = 0;
              for (Int_t d=5; d>=0; d--){
                badpar5[i] = badpar5[i]*mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5[i][d];
              }                // calculate the parameters using pol5
            }
            for (Int_t ipar=0;ipar<4;ipar++){
              status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
            }
          }
        }
        else{
          int tsector = sector + 1;
          double mom_scpd =p;
          // sector 2 has one bad paddles
          if (tsector == 2){
            float badpar2[2];            // 4 parameters to determine the positions of the two theta gaps
            for (int i=0; i<2; i++){
              badpar2[i] = 0;
              // calculate the parameters using 1/p
              badpar2[i] = fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS2_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(int ipar=0;ipar<1;ipar++)
              status = status && !(theta>badpar2[2*ipar] && theta<badpar2[2*ipar+1]);
          }
          if (tsector == 3){
            float badpar3[8];            // 4 parameters to determine the positions of the two theta gaps
            for (int i=0; i<8; i++){
              badpar3[i] = 0;
              // calculate the parameters using 1/p
              badpar3[i] = fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS3_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(int ipar=0;ipar<4;ipar++)
              status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
          // sector 4 has two bad paddle
          else if (tsector == 4){
            float badpar4[4];     // 2 parameters to determine the position of the theta gap
            for (int i=0; i<4; i++){
              badpar4[i] = 0;
              // calculate the parameters using 1/p
              badpar4[i] = fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS4_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(int ipar=0;ipar<2;ipar++)
	      status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
          // sector 5 has four bad paddles
          else if (tsector == 5){
            Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
            for (Int_t i=0; i<8; i++){
              badpar5[i] = 0;
              // calculate the parameters using 1/p
              badpar5[i] = fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][0] + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][1]/mom_scpd + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_4Gev_2250_Pfid_ScpdS5_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
            }
            for(Int_t ipar=0;ipar<4;ipar++)
              status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
          }
        }
      }

    }

    // ---------------------------------------------------------------------------------------------

    // apapadop: Nov 23 2020, adding a lower theta bound of 12 deg for the pi plus
    if (!PiplFiducialCutExtra(beam_en, momentum)) { status = false; } 

    // ---------------------------------------------------------------------------------------------

    return status;

  } // end of the if statement for the non-"" case

}

// ----------------------------------------------------------------------------------------------------------------------------------

Bool_t Fiducial::PFiducialCutExtra(double beam_en, TVector3 momentum) {

  double theta = momentum.Theta() * TMath::RadToDeg() ;

  if (theta < conf::kMinThetaProton) return false;
  if (theta >= conf::kMaxThetaHadrons ) return false; 
  return true ;

}

// ----------------------------------------------------------------------------------------------------------------------------------

Bool_t Fiducial::PiplFiducialCutExtra(double beam_en, TVector3 momentum) {

  double theta = momentum.Theta() * TMath::RadToDeg() ;

  if (theta < conf::kMinThetaPiPlus) return false;
  if (theta >= conf::kMaxThetaHadrons ) return false;
  return true;
}

// ----------------------------------------------------------------------------------------------------------------------------------

Bool_t Fiducial::Phot_fidExtra(TVector3 momentum) {

  double theta = momentum.Theta() * TMath::RadToDeg() ;
  
  if (theta < conf::kMinThetaGamma) return false;
  if (theta >= conf::kMaxThetaHadrons ) return false;
  return true;

}

// ----------------------------------------------------------------------------------------------------------------------------------

Bool_t Fiducial::PimiFiducialCutExtra(double beam_en, TVector3 momentum) {
  double theta = momentum.Theta() * TMath::RadToDeg() ;
  double mom = momentum.Mag();
  double theta_min = myPiMinusFit->Eval(mom);
  if (theta < theta_min) return false;
  if (theta >= conf::kMaxThetaHadrons ) return false;
  return true;
}

// ----------------------------------------------------------------------------------------------------------------------------------

//using two GeV pimi fiducial cuts for both 2 and 4 GeV analysis

//PimiFiducialCut - Fiducial cuts function for e2a analysis
//                  Two step process: 1) Cut on phi fiducial edges
//                                    2) Remove bad TOF paddles
//                  Uses 2 GeV pimi fiducial cuts for both 2 and 4 GeV analysis, from electron cuts plus low momentum extension
//                  October 27th, 2020 - Theta gap parameters at 1 GeV (750A), and 2/4 GeV (2250A) updated to replace reuse of
//                                       electron gap functions, which are too messy to reuse accurately and account for CC cuts inappropriate to pions
//                                       This update also includes Florian's update to loosen the constraint on upper limit of theta, which at high
//                                       momentum is still much too tight
//                                       A future update will implement new phi fiducial cuts at all momenta, but the current bug fixes are considered
//                                       adequate for the zero pion analyses

Bool_t Fiducial::PimiFiducialCut(double beam_en, TVector3 momentum, Float_t *pimi_philow, Float_t *pimi_phiup){
  
  if (beam_en == 0) {

    bool status = true;

    double phi = GetPhi(momentum);
    double theta = GetTheta(momentum);

    if ( !(TMath::Abs(phi - conf::kCenterFirstSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterSecondSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterThirdSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterFourthSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterFifthSector) < conf::kPhiOpeningAngle || TMath::Abs(phi - conf::kCenterSixthSector) < conf::kPhiOpeningAngle) ) { status = false; }

    if (theta < conf::kMinThetaPiMinus) { status = false; }

    return status;

  } else {

    // Electron fiducial cut, return kTRUE if pass or kFALSE if not
    //--------------------------------------------------------------
    //Preamble, set up common variables

    Bool_t status = kTRUE;  //set initial return status true, logic below will decide if that changes

    //Variables used by the cuts
    TVector3 mom = momentum;
    double phi = mom.Phi();
    if (phi < -M_PI/6.) phi+= 2.*M_PI;
    int sector = (phi+M_PI/6.)/(M_PI/3.);
    sector = sector%6;

    int tsector = sector + 1;   //temporary variable for counting sector number, given arrays start at 0

    double phi_deg = phi * 180./M_PI;
    phi_deg -= sector*60;  //fiducial cuts defined using individual sector distributions, redrawn to be centred on 0 degrees in phi. This angle allows those cut definitions to be applied

    double theta = mom.Theta();
    double theta_deg = theta * 180./M_PI;

    double mom_pi = mom.Mag();

    Float_t p_theta = mom_pi, thetamax = 0;
    //--------------------------------------------------------------


    //--------------------------------------------------------------
    //1.1GeV - uses electron fiducial cuts
    if(beam_en>1. &&  beam_en<2.) {

      Double_t phipars[5]={0,0,0,0,0};
      Double_t phicutoff;

      status = true;
      //piminus cuts at 1 GeV only valid in momentum range 150 MeV - 1.1 GeV
      //There's no data outside this momentum range at 1 GeV, it's an excessive condition, but does no harm
      if (mom_pi < .15){
	mom_pi = .15;
      }
      if (mom_pi > 1.1){
	mom_pi = 1.1;
      }

      //Torus setting of 1500 A
      //Not used in Mariana's analysis, so are these still valid?
      if( fTorusCurrent>1490 && fTorusCurrent<1510){

	for(Int_t mompar=0;mompar<6;mompar++) {
	  for(Int_t phipar=0;phipar<5;phipar++) {
	    phipars[phipar]+=fgPar_1gev_1500_Pimfid[sector][phipar][mompar]*pow(mom_pi,mompar);
	    //std::cout << mom_e << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
	  }
	}

	//cut on phi fiducial edges
	if(phi_deg<=0) {
	  phicutoff = phipars[1]*(1.-(1./((theta_deg-phipars[4])/phipars[3]+1.)));
	  status = ((phi_deg>phicutoff) && (theta_deg>phipars[4]));
	}
	else {
	  phicutoff = phipars[0]*(1.-(1./((theta_deg-phipars[4])/phipars[2]+1.)));
	  status = ((phi_deg<phicutoff) && (theta_deg>phipars[4]));
	}

	//////////////////////////////////Remove bad TOF paddles ///////////////////////////////////
	//This takes the form of a momentum-dependent cut on theta, defined from stored parameters//
	//There are two sets of gap arameters - above and below 300 MeV momentum. Clearly an      //
	//attempt has been made to extend electron cuts at this field setting too                 //
	////////////////////////////////////////////////////////////////////////////////////////////
	if (mom_pi >= .3) //the electron cuts which are reused for pi minus are only defined from 300 MeV
	  {
	    bool SCpdcut = true;
	    if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
	      if (status){
		// sector 3 has two bad paddles
		if (tsector == 3){
		  float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
		  for (int i=0; i<4; i++){
		    badpar3[i] = 0;
		    // calculate the parameters using pol7
		    for (int d=7; d>=0; d--){
		      badpar3[i] = badpar3[i]*mom_pi + fgPar_1gev_1500_Pimfid_Theta_S3[i][d];
		    }
		  }
		  for(int ipar=0;ipar<2;ipar++)
		    status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
		}
		// sector 4 has one bad paddle
		else if (tsector == 4){
		  float badpar4[2];     // 2 parameters to determine the position of the theta gap
		  for (int i=0; i<2; i++){
		    badpar4[i] = 0;
		    // calculate the parameters using pol7
		    for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_pi + fgPar_1gev_1500_Pimfid_Theta_S4[i][d];}
		  }
		  status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
		}
		// sector 5 has four bad paddles
		else if (tsector == 5){
		  Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
		  for (Int_t i=0; i<8; i++){
		    badpar5[i] = 0;
		    // calculate the parameters using pol7
		    for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_pi + fgPar_1gev_1500_Pimfid_Theta_S5[i][d];}
		  }
		  if (mom_pi<1.25) badpar5[0] = 23.4*1500/2250;
		  if (mom_pi<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.
		  for(Int_t ipar=0;ipar<4;ipar++)
		    status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
		}
	      }
	    }

	    // ---------------------------------------------------------------------------------------------

	    // apapadop: Nov 23 2020, adding a lower theta bound of A + B/P for the pi minus
	    if (!PimiFiducialCutExtra(beam_en, momentum)) { status = false; } 

	    // ---------------------------------------------------------------------------------------------

	    return status;
	  }
	else{     //momentum below 300 MeV
	  bool SCpdcut = true;
	  if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
	    if (status){
	      // sector 3 has two bad paddles
	      if (tsector == 3){
		float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
		for (int i=0; i<4; i++){
		  badpar3[i] = 0;
		  // calculate the parameters using 1/p
		  badpar3[i] = fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][1]/mom_pi + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][2]/(mom_pi*mom_pi) + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][3]/(mom_pi*mom_pi*mom_pi);
		}
		for(int ipar=0;ipar<2;ipar++)
		  status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
	      }
	      // sector 4 has one bad paddle
	      else if (tsector == 4){
		float badpar4[2];     // 2 parameters to determine the position of the theta gap
		for (int i=0; i<2; i++){
		  badpar4[i] = 0;
		  // calculate the parameters using 1/p
		  badpar4[i] = fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][1]/mom_pi + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][2]/(mom_pi*mom_pi) + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][3]/(mom_pi*mom_pi*mom_pi);
		}
		status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
	      }
	      // sector 5 has four bad paddles
	      else if (tsector == 5){
		Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
		for (Int_t i=0; i<8; i++){
		  badpar5[i] = 0;
		  // calculate the parameters using 1/p
		  badpar5[i] = fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][1]/mom_pi + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][2]/(mom_pi*mom_pi) + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][3]/(mom_pi*mom_pi*mom_pi);
		}
		if (mom_pi<1.25) badpar5[0] = 23.4*1500/2250;
		if (mom_pi<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.
		for(Int_t ipar=0;ipar<4;ipar++)
		  status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
	      }
	    }
	  }

	  // ---------------------------------------------------------------------------------------------

	  // apapadop: Nov 23 2020, adding a lower theta bound of A + B/P for the pi minus
	  if (!PimiFiducialCutExtra(beam_en, momentum)) { status = false; } 

	  // ---------------------------------------------------------------------------------------------

	  return (status);
	}
      }
      //------------------------------
      //Torus current setting of 750 A
      //------------------------------
      if ( fTorusCurrent>740 && fTorusCurrent<760){

	for(Int_t mompar=0;mompar<6;mompar++) {
	  for(Int_t phipar=0;phipar<5;phipar++) {
	    phipars[phipar]+=fgPar_1gev_750_Pimfid[sector][phipar][mompar]*pow(mom_pi,mompar);
	    //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
	  }
	}

	//Temporary momentum parameter to define a maximum limit on theta in the range 0.1-0.7 GeV pion momentum
	//Very few events 700 MeV at 1.1 GeV beam energy, and the momentum threshold for pion detection is 150 MeV
	//All this does is remove some insignificant scruff
	if (p_theta>0.7)  p_theta = 0.7;
	else if(p_theta<0.1)   p_theta = 0.1;

	//Calculate theta limit from stored parameters
	for(int i=4;i>=0;i--){
	  thetamax = thetamax*p_theta+pimi_thetamax1[i];
	}

	//cut on phi fiducial edges
	if(phi_deg<=0) {
	  phicutoff = phipars[1]*(1.-(1./((theta_deg-phipars[4])/phipars[3]+1.)));
	  status = ((phi_deg>phicutoff) && (theta_deg>phipars[4]) && theta_deg<=thetamax);
	}
	else {
	  phicutoff = phipars[0]*(1.-(1./((theta_deg-phipars[4])/phipars[2]+1.)));
	  status = ((phi_deg<phicutoff) && (theta_deg>phipars[4]) && theta_deg<=thetamax);
	}

	///////////////////////////////////Remove bad TOF paddles //////////////////////////////////
	//This takes the form of a momentum-dependent cut on theta, defined from stored parameters//
	////////////////////////////////////////////////////////////////////////////////////////////
	bool SCpdcut = true;
	if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
	  if (status){
	    mom_pi = mom.Mag();
	    //Reset ignores mom_e > 1.1 cut from above. Is this okay? F.H. 10/31/19  - Yes, as it applies the theta cut functions at all momenta SF 10/11/20

	    //could probably replace parsec[1,6]_{l/h} with
	    //double thetaCutLo, thetaCutHi;
	    //and eventually propagate to all theta gap determinations

	    //Using a set of fid_1gev_750_pimifid_S[1,6] parameters, hard coded in Fiducial.h, to avoid confusion with electron, althoug some parameters are the same
	    //sector 1 has one gap
	    if(tsector == 1){
	      double parsec1_l,parsec1_h;
	      //Evaluate theta gap boundaries from saved parameters
	      parsec1_l= fid_1gev_750_pimifid_S1[0][0]+fid_1gev_750_pimifid_S1[0][1]/mom_pi +fid_1gev_750_pimifid_S1[0][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S1[0][3]/(mom_pi*mom_pi*mom_pi);
	      parsec1_h= fid_1gev_750_pimifid_S1[1][0]+fid_1gev_750_pimifid_S1[1][1]/mom_pi +fid_1gev_750_pimifid_S1[1][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S1[1][3]/(mom_pi*mom_pi*mom_pi);

	      status = status && !(theta_deg>parsec1_l && theta_deg<parsec1_h);
	    }
	    // //sector 2 has one gap, however, it doesn't seem to be present in pi minus, it may be a CC gap and only applicable to electrons
	    //           if(tsector == 2){
	    //             double parsec2_l,parsec2_h;
	    //             if(mom_pi<0.4) mom_pi=0.4;
	    // 	    //Evaluate theta gap boundaries from saved parameters
	    //             parsec2_l= fid_1gev_750_pimifid_S2[0][0]+fid_1gev_750_pimifid_S2[0][1]/mom_pi +fid_1gev_750_pimifid_S2[0][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S2[0][3]/(mom_pi*mom_pi*mom_pi);
	    //             parsec2_h= fid_1gev_750_pimifid_S2[1][0]+fid_1gev_750_pimifid_S2[1][1]/mom_pi +fid_1gev_750_pimifid_S2[1][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S2[1][3]/(mom_pi*mom_pi*mom_pi);

	    //             status = status && !(theta_deg>parsec2_l && theta_deg<parsec2_h);
	    //           }
	    //sector 3 has three gaps, removing CC gap from electron parameters
	    else if(tsector == 3){
	      double parsec3_l[3],parsec3_h[3]; //do these really need to be arrays?
	      for(int d=0;d<3;d++){
		//Evaluate theta gap boundaries from saved parameters
		parsec3_l[d]= fid_1gev_750_pimifid_S3[d][0][0]+fid_1gev_750_pimifid_S3[d][0][1]/mom_pi +fid_1gev_750_pimifid_S3[d][0][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S3[d][0][3]/(mom_pi*mom_pi*mom_pi);
		parsec3_h[d]= fid_1gev_750_pimifid_S3[d][1][0]+fid_1gev_750_pimifid_S3[d][1][1]/mom_pi +fid_1gev_750_pimifid_S3[d][1][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S3[d][1][3]/(mom_pi*mom_pi*mom_pi);

		status = status && !(theta_deg>parsec3_l[d] && theta_deg<parsec3_h[d]);
	      }
	    }
	    //sector 4 has two gaps, however, the first gap is probably a CC gap, and not applicable to the pi minus, an additional TOF gap replaces this
	    else if(tsector == 4){
	      double parsec4_l[2],parsec4_h[2];
	      for(int d=0;d<2;d++){
		//Evaluate theta gap boundaries from saved parameters
		parsec4_l[d]= fid_1gev_750_pimifid_S4[d][0][0]+fid_1gev_750_pimifid_S4[d][0][1]/mom_pi +fid_1gev_750_pimifid_S4[d][0][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S4[d][0][3]/(mom_pi*mom_pi*mom_pi);
		parsec4_h[d]= fid_1gev_750_pimifid_S4[d][1][0]+fid_1gev_750_pimifid_S4[d][1][1]/mom_pi +fid_1gev_750_pimifid_S4[d][1][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S4[d][1][3]/(mom_pi*mom_pi*mom_pi);

		status = status && !(theta_deg>parsec4_l[d] && theta_deg<parsec4_h[d]);
	      }
	    }
	    //sector 5 has two gaps, omitting the CC gaps from electrons and adding a new pi minus gap
	    else if(tsector == 5){
	      double parsec5_l[2],parsec5_h[2];
	      for(int d=0;d<2;d++){
		//Evaluate theta gap boundaries from saved parameters
		parsec5_l[d]= fid_1gev_750_pimifid_S5[d][0][0]+fid_1gev_750_pimifid_S5[d][0][1]/mom_pi +fid_1gev_750_pimifid_S5[d][0][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S5[d][0][3]/(mom_pi*mom_pi*mom_pi);
		parsec5_h[d]= fid_1gev_750_pimifid_S5[d][1][0]+fid_1gev_750_pimifid_S5[d][1][1]/mom_pi +fid_1gev_750_pimifid_S5[d][1][2]/(mom_pi*mom_pi) +fid_1gev_750_pimifid_S5[d][1][3]/(mom_pi*mom_pi*mom_pi);

		status = status && !(theta_deg>parsec5_l[d] && theta_deg<parsec5_h[d]);
	      }
	    }
	  }
	}

	// ---------------------------------------------------------------------------------------------

	// apapadop: Nov 23 2020, adding a lower theta bound of A + B/P for the pi minus
	if (!PimiFiducialCutExtra(beam_en, momentum)) { status = false; } 

	// ---------------------------------------------------------------------------------------------

	return (status);
      }
    }


    //--------------------------------------------------------------
    //2.2 and 4.4 GeV - uses 2.2 GeV electron cuts, plus a weird fudge at pion momenta below 350 MeV //This MUST be replaced by bespoke cuts
    //--------------------------------------------------------------
    if ( beam_en>2. &&  beam_en<5. && fTorusCurrent>2240 && fTorusCurrent<2260){

      Float_t phimin, phimax;

      //        Float_t par[6];               // six parameters to determine the outline of Theta vs Phi  //why six parameters? it's five everywhere else...
      //p_theta == mom_pi
      //theta maximum limit functions only defined in the momentum range 0.1 - 2.075 //is this an appropriate maximum at 4 GeV?
      if (p_theta>2.075){
	p_theta = 2.075;
      }
      else if(p_theta<0.1){
	p_theta = 0.1;
      }

      //calculate theta limit from loaded parameters (is this thetamax[0] + mom*thetamax[1] + mom^2*thetamax[2] + mom^3*thetamax[3] + mom^4*thetamax[4]?)
      for(int i=4;i>=0;i--){
	thetamax = thetamax*p_theta + pimi_thetamax2and4[i]; //upper theta limit for pi- at different momentum
      }
      //std::cout << "pion minus: mom " << p_theta << " , phi " << phi << " , phi_deg " << phi_deg << " , sector " << sector << " , thetamax " << thetamax << std::endl;
      //---These are the reused electron parameters
      if(mom_pi > 0.35){   //theta vs phi outline for high p region obtained by Bin

        if(mom_pi > 2.){
          mom_pi = 2.; //to extrapolate the cut to higher momenta for pimi //(to badly extrapolate) S.F. October 2020
        }

        Float_t par[6];               // six parameters to determine the outline of Theta vs Phi  //why six parameters? it's five everywhere else...

        for (Int_t i=0; i<6; i++){
          par[i] = 0;
          for (Int_t d=8; d>=0; d--){
            par[i] = par[i]*mom_pi +  fgPar_2GeV_2250_Efid[sector][i][d];
          }                          // calculate the parameters using pol8
        }
        double phi_min_limit = par[3]/((50.0-par[0])+par[3]/par[2])-par[2]; //calculation of phi limit from theta at 50 degree
        double phi_max_limit = par[2]-par[3]/((50.0-par[0])+par[3]/par[2]); //calculation of phi limit from theta at 50 degree

	//  std::cout << "par[0] " << par[0] << " , par[1] " << par[1] <<  std::endl;
        if (phi_deg < 0) {
          Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi_deg);
          phimin =  par[3]/((theta_deg-par[0])+par[3]/par[2])-par[2];
          phimax =  par[2]-par[3]/((theta_deg-par[0])+par[3]/par[2]);
          *pimi_philow = phimin;
          *pimi_phiup = phimax;
          //Modification by F.H. 10/20/20
          //theta range was limited previously by a cut  on theta_deg<par[1]
          //cut is removed and a cut on phi to phi_limit is added similar to piplus at 2.2 where limits are hardcoded
          //also a cut for thetamax is added
          status = (theta_deg>tmptheta && tmptheta>=par[0] && phi_deg>=phi_min_limit && theta_deg<=thetamax);
	  //  std::cout << "status " << status << " , theta_deg " << theta_deg << ", tmptheta " << tmptheta << " , phimin " << phimin << " , phimax " << phimax << " , phi_deg " << phi_deg << " , sector " << sector << " , thetamax " << thetamax << std::endl;
        }
        else {
          Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi_deg);
          phimin =  par[5]/((theta_deg-par[0])+par[5]/par[4])-par[4];
          phimax =  par[4]-par[5]/((theta_deg-par[0])+par[5]/par[4]);
          *pimi_philow = phimin;
          *pimi_phiup = phimax;
          //Modification by F.H. 10/20/20
          //theta range was limited previously by a cut  on theta_deg<par[1]
          //cut is removed and a cut on phi to phi_limit is added similar to piplus at 2.2 where limits are hardcoded
          //also a cut for thetamax is added
          status = (theta_deg>tmptheta && tmptheta>=par[0] && phi_deg<=phi_max_limit && theta_deg<=thetamax);
        }
      }//end of high momentum cut

      //---The pi- parameters obtained only at low momentum
      else if(mom_pi<=0.35){     //theta vs phi outline for low p obtained by Mariana

	//valid only in the range 0.125 - 0.325 GeV
        if(mom_pi>0.325)mom_pi=0.325; //This causes the notch at the lower theta limit between 325 and 350 MeV
        else if (mom_pi<0.125)mom_pi=0.125;//This shouldn't do anything at all, as pions are selected from 150 MeV

        Float_t params[6];               // six parameters to determine the outline of Theta vs Phi, like the higher momentum region above.  Can change to par[6] and share declaration with above
        for (Int_t i=0; i<6; i++){
          params[i] = 0;
          for (Int_t d=4; d>=0; d--){
            params[i] = params[i]*mom_pi + fid_2gev_2250_pimifid_outline[sector][i][d];
          }                          // calculate the parameters using pol4
        }
        if (phi_deg < 0) {
          phimin =  params[3]/((theta_deg-params[0])+params[3]/params[1])-params[1];  //why is phimin calculated differently in the low momentum?
	  //where's phimax?
          status = (phi_deg>phimin && theta_deg>params[0] && theta_deg<=thetamax);
        }
        else {
	  //where's phimin?
          phimax = params[2]-params[4]/((theta_deg-params[0])+params[4]/params[2]);   //why is phimax calculated differently in the low momentum?
          status = (phi_deg<phimax && theta_deg>params[0]  && theta_deg<=thetamax);
	}
      }//end of low momentum cut

      else{ std::cout << "pi- fiducial outline not determined" << std::endl;}

      ////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////Remove bad TOF paddles ///////////////////////////////////
      //This takes the form of a momentum-dependent cut on theta, defined from stored parameters//
      ////////////////////////////////////////////////////////////////////////////////////////////

      bool SCpdcut=true;

      // by now, we have checked if the electron is within the outline of theta vs phi plot
      if (SCpdcut){  // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
	if (status){ //if event not already rejected

	  //---Further gaps obtained with pi- ------
	  //Gaps obtained with pi-. Includes redefined versions of some electron gaps from fgPar_2GeV_2250_EfidThets_S{1-6}. S.F Oct 26, 2020
	  mom_pi = momentum.Mag();  //reset momentum to event value, just in case the electron paddle nonsense causes problems
	  //sector 1 has two gaps
	  if(tsector == 1){
	    double parsec1_l[2],parsec1_h[2];   //theta gap paramaters
	    for(int d=0;d<2;d++){
	      //-----Define theta limits from stored parameters
	      parsec1_l[d]= fid_2gev_2250_pimifid_S1[d][0][0]+fid_2gev_2250_pimifid_S1[d][0][1]/mom_pi +fid_2gev_2250_pimifid_S1[d][0][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S1[d][0][3]/(mom_pi*mom_pi*mom_pi);
	      parsec1_h[d]= fid_2gev_2250_pimifid_S1[d][1][0]+fid_2gev_2250_pimifid_S1[d][1][1]/mom_pi +fid_2gev_2250_pimifid_S1[d][1][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S1[d][1][3]/(mom_pi*mom_pi*mom_pi);
	      status=status && !(theta_deg>parsec1_l[d] && theta_deg<parsec1_h[d]);
	    }
	  }
	  //sector 3 has 3 gaps
	  else if(tsector == 3){
	    double parsec3_l[3],parsec3_h[3];
	    for(int d=0;d<3;d++){
	      //-----Define theta limits from stored parameters
	      parsec3_l[d]= fid_2gev_2250_pimifid_S3[d][0][0]+fid_2gev_2250_pimifid_S3[d][0][1]/mom_pi +fid_2gev_2250_pimifid_S3[d][0][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S3[d][0][3]/(mom_pi*mom_pi*mom_pi);
	      parsec3_h[d]= fid_2gev_2250_pimifid_S3[d][1][0]+fid_2gev_2250_pimifid_S3[d][1][1]/mom_pi +fid_2gev_2250_pimifid_S3[d][1][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S3[d][1][3]/(mom_pi*mom_pi*mom_pi);

	      status=status && !(theta_deg>parsec3_l[d] && theta_deg<parsec3_h[d]);
	    }
	  }
	  //sector 4 has two gaps, one defined on pi minus, one reimplemented from electron
	  else if(tsector == 4){
	    double parsec4_l[2],parsec4_h[2];
	    for(int d=0;d<2;d++){
	      if(d==0){
		if(mom_pi<0.15)mom_pi=0.15;     //not really needed, below pion detection threshold
		else if(mom_pi>0.5)mom_pi=0.5;  //not really needed, very little data here
	      }
	      //-----Define theta limits from stored parameters
	      parsec4_l[d]= fid_2gev_2250_pimifid_S4[d][0][0]+fid_2gev_2250_pimifid_S4[d][0][1]/mom_pi +fid_2gev_2250_pimifid_S4[d][0][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S4[d][0][3]/(mom_pi*mom_pi*mom_pi);
	      parsec4_h[d]= fid_2gev_2250_pimifid_S4[d][1][0]+fid_2gev_2250_pimifid_S4[d][1][1]/mom_pi +fid_2gev_2250_pimifid_S4[d][1][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S4[d][1][3]/(mom_pi*mom_pi*mom_pi);
	      status=status && !(theta_deg>parsec4_l[d] && theta_deg<parsec4_h[d]);
	    }
	  }
	  //sector 5 has two gaps, plus two moved from electron gaps
	  else if(tsector == 5){
	    double parsec5_l[3],parsec5_h[3];
	    for(int d=0;d<4;d++){
	      //-----Define theta limits from stored parameters
	      parsec5_l[d]= fid_2gev_2250_pimifid_S5[d][0][0]+fid_2gev_2250_pimifid_S5[d][0][1]/mom_pi +fid_2gev_2250_pimifid_S5[d][0][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S5[d][0][3]/(mom_pi*mom_pi*mom_pi);
	      parsec5_h[d]= fid_2gev_2250_pimifid_S5[d][1][0]+fid_2gev_2250_pimifid_S5[d][1][1]/mom_pi +fid_2gev_2250_pimifid_S5[d][1][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S5[d][1][3]/(mom_pi*mom_pi*mom_pi);
	      status=status && !(theta_deg>parsec5_l[d] && theta_deg<parsec5_h[d]);
	    }
	  }
	  //sector 6 has one gap, a gap that could be replaced by a box cut between 58 and 65 degrees to some momentum limit (probably a DC region failure)
	  if(tsector == 6){
	    double parsec6_l,parsec6_h;
	    if(mom_pi<0.175)mom_pi=0.175;
	    else if(mom_pi>0.275)mom_pi=0.275;
	    parsec6_l= fid_2gev_2250_pimifid_S6[0][0]+fid_2gev_2250_pimifid_S6[0][1]/mom_pi +fid_2gev_2250_pimifid_S6[0][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S6[0][3]/(mom_pi*mom_pi*mom_pi);
	    parsec6_h= fid_2gev_2250_pimifid_S6[1][0]+fid_2gev_2250_pimifid_S6[1][1]/mom_pi +fid_2gev_2250_pimifid_S6[1][2]/(mom_pi*mom_pi) +fid_2gev_2250_pimifid_S6[1][3]/(mom_pi*mom_pi*mom_pi);
	    status=status && !(theta_deg>parsec6_l && theta_deg<parsec6_h);
	  }
	  ////end of gaps with pi-
	}
      }
    } //end of 2 and 4GEV pi minus fiducial

    // ---------------------------------------------------------------------------------------------

    // apapadop: Nov 23 2020, adding a lower theta bound of A + B/P for the pi minus
    if (!PimiFiducialCutExtra(beam_en, momentum)) { status = false; } 

    // ---------------------------------------------------------------------------------------------

    return status;

  } // end of the if statemnet for the non-"" case

}

// --------------------------------------------------------------------------------------------------------------------------------------------------


bool Fiducial::Phot_fid(TVector3 V3_phot){

  bool status = true;
  //  double costheta=V3_phot.Pz();
  double costheta=V3_phot.CosTheta();
  //double theta_deg=TMath::ACos(V3_phot.Pz())*TMath::RadToDeg(); //extraneous variable
  double phi_deg=TMath::ATan2(V3_phot.Py(),V3_phot.Px())*TMath::RadToDeg()+30;
  if(phi_deg<0)phi_deg=phi_deg+360;
  bool hot_spot=(phi_deg>185 && phi_deg<191 && costheta<0.71 && costheta>0.67) || (phi_deg>221 && phi_deg<236 && costheta<0.73 && costheta>0.67); // used to kill the two hot spots in sector 4


  if((costheta>low_lim1_ec->Eval(phi_deg) && costheta<up_lim1_ec->Eval(phi_deg) && costheta<rightside_lim1_ec->Eval(phi_deg) && costheta<leftside_lim1_ec->Eval(phi_deg))  ||
     (costheta>low_lim2_ec->Eval(phi_deg) && costheta<up_lim2_ec->Eval(phi_deg) && costheta<rightside_lim2_ec->Eval(phi_deg) && costheta<leftside_lim2_ec->Eval(phi_deg))       ||
     (costheta>low_lim3_ec->Eval(phi_deg) && costheta<up_lim3_ec->Eval(phi_deg) && costheta<rightside_lim3_ec->Eval(phi_deg) && costheta<leftside_lim3_ec->Eval(phi_deg))       ||
     (costheta>low_lim4_ec->Eval(phi_deg) && costheta<up_lim4_ec->Eval(phi_deg) && costheta<rightside_lim4_ec->Eval(phi_deg) && costheta<leftside_lim4_ec->Eval(phi_deg) && !hot_spot)       ||
     (costheta>low_lim5_ec->Eval(phi_deg) && costheta<up_lim5_ec->Eval(phi_deg) && costheta<rightside_lim5_ec->Eval(phi_deg) && costheta<leftside_lim5_ec->Eval(phi_deg))       ||
     (costheta>low_lim6_ec->Eval(phi_deg) && costheta<up_lim6_ec->Eval(phi_deg) && costheta<rightside_lim6_ec->Eval(phi_deg) && costheta<leftside_lim6_ec->Eval(phi_deg))
     ){ //EC only
    status=true;}
  else {
    status=false;}

  return status;
}

// ----------------------------------------------------------------------------------------------------

bool Fiducial::Pi_phot_fid_united(double beam_en, TVector3 V3_pi_phot, int q_pi_phot){

  bool status = false;
  Float_t pi_cphil=0,pi_cphir=0,pi_phimin=0,pi_phimax=0;

  if(q_pi_phot==0) status=Phot_fid(V3_pi_phot);
  if(q_pi_phot>0) status=PiplFiducialCut(beam_en, V3_pi_phot, &pi_cphil, &pi_cphir);
  if(q_pi_phot<0) status=PimiFiducialCut(beam_en, V3_pi_phot, &pi_phimin, &pi_phimax);
  return status;

}

// ----------------------------------------------------------------------------------------------------

bool Fiducial::Pi_phot_fid_unitedExtra(double beam_en, TVector3 V3_pi_phot, int q_pi_phot){

  bool status = true;
  if(q_pi_phot==0) status=Phot_fidExtra(V3_pi_phot);
  if(q_pi_phot>0) status=PiplFiducialCutExtra(beam_en, V3_pi_phot);
  if(q_pi_phot<0) status=PimiFiducialCutExtra(beam_en, V3_pi_phot);
  return status;
}

// ----------------------------------------------------------------------------------------------------
#endif
