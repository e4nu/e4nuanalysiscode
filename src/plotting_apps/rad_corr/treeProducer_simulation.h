#ifndef treeProducer_simulation_h
#define treeProducer_simulation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>

#include <sstream>
#include <iostream>

using namespace std;

TString FSI = "_FSI";
TString E = "4_325";

TString Target = "1H";        
TString FSIModel = "hA2018_LFG"; 
TString WhichSample = "EMradcorrLarryWeights";
TString graphstring = "";

TString FullPath = "./mySamples";

char * env = std::getenv("E4NUANALYSIS") ; 
std::string path = env ;
std::string FullInputName = path+"/data/1H_4_325_Data_Plots_FSI_em.root"; 
//TString file_name = "/genie/app/users/jtenavid/Software/e4v/E4NuAnalysis/Source/e4nuanalysiscode/src/plotting_apps/rad_corr/RadFlux_G18_10a_H_4320MeV.gst.root";
TString file_name = "/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/RadCorr_simple/master-routine_validation_01-eScattering/radcorr/G18_10a_Radcorr_simple_H_4325MeV.gst.root";
//"/genie/app/users/jtenavid/Software/e4v/E4NuAnalysis/Source/e4nuanalysiscode/e_on_1000010010_4325MeV_H_G18_10a_radiated.gst.root";
//"/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/MonoFlux/master-routine_validation_01-eScattering/G18_10a_H_4325MeV_MonoFlux.gst.root";
//
//"/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/G18_10a/MonoFlux/master-routine_validation_01-eScattering/e_on_1000010010_4325MeV_H_G18_10a_radiated.gst.root";
//"/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/MonoFlux/master-routine_validation_01-eScattering/G18_10a_H_4325MeV_MonoFlux.gst.root"; 
//"/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/RadCorr/master-routine_validation_01-eScattering/G18_10a_H_4325MeV_RadFlux.gst.root";
  //
//"/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/AdiMethod_simple/master-routine_validation_01-eScattering/G18_10a_H_4325MeV_AdiSimple.gst.root";
  //"/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/AdiMethod/master-routine_validation_01-eScattering/G18_10a_Adi_simc_H_4325MeV.gst.root";
  //

  //"/pnfs/genie/scratch/users/jtenavid/GENIE_e4nu_Generations/2024Generation/G18_10a/RadCorr/master-routine_validation_01-eScattering/G18_10a_H_4325MeV_RadFlux.gst.root";
  //
  //
//TString file_name = "/genie/app/users/jtenavid/Software/genie-v3/GeneratorRad/e_on_1000010010_4325MeV_simple.gst.root";


bool SCpdcut = true;

double H3_bind_en = 0.008481, He4_bind_en = 0.0283, C12_bind_en = 0.09215, B_bind_en = 0.0762, He3_bind_en = 0.0077, D2_bind_en = 0.00222, Fe_bind_en = 0.49226, Mn_bind_en = 0.4820764;

float fgPar_4Gev_2250_Phi[6][3], fgPar_4Gev_2250_Theta[6][4], fgPar_4Gev_2250_Efid_t0_p[6][2];
float fgPar_4Gev_2250_Efid_t1_p[6][6], fgPar_4Gev_2250_Efid_b_p[6][2][6], fgPar_4Gev_2250_Efid_a_p[6][2][6];
double NRotations = 100;

const Float_t par_EcUVW[6][3] = {{60, 360, 400}, {55, 360, 400}, {50, 363, 400}, {52, 365, 396}, {60, 360, 398}, {50, 362, 398}};

const Float_t fgPar_2Gev_2250_Phi[6][3] = {{0.9903, -3.926E-4, 1.318E-5}, {0.9803, -3.177E-4, 1.706E-5}, {0.9922, 1.836E-4, 1.474E-5}, {0.9898, 5.259E-5, 1.45E-5}, {0.9906, -1.604E-4, 1.66E-5}, {0.9925, 1.902E-4, 1.985E-5}};
const Float_t fgPar_2Gev_2250_Theta[3][4] = {{8.56526E-4, 7.89140E+1, -8.41321E-1, 1.00082}, {6.10625E-1, 8.30600E-1, -4.40544E-1, 0}, {-5.02481, 1.29011E+1, -6.90397, 0}}; //corrections for electron momentum obtained with e1c 2.5Gev 2250A data set (Run 16719 and 16720)
const Float_t fgPar_2GeV_2250_Efid[6][6][9] = {{{62.2935, -92.5133, 87.0360, -38.4696, 6.3177, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {-140.845, 1381.30, -4499.99, 7557.27, -7140.27, 3828.75, -1086.21, 126.468, 0}, {497.951, -1846.42, 2759.58, -1634.71, 345.006, 0, 0, 0, 0}, {9.40986, 180.752, -646.771, 1055.14, -909.094, 424.435, -99.8368, 9.02086, 0}, {288.485, -1016.03, 1463.72, -859.231, 185.976, 0, 0, 0, 0}}, {{61.1474, -88.768, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {21.3087, 138.975, -672.710, 1324.20, -1326.12, 714.866, -197.531, 21.9144, 0}, {375.091, -1411.50, 2082.58, -1192.17, 239.685, 0, 0, 0, 0}, {-121.816, 1182.59, -3800.98, 6319.82, -5937.33, 3179.37, -903.954, 105.764, 0}, {-4781.96, 43165.9, -159567, 318502, -376469, 271207, -116893, 27698.9, -2775.61}}, {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {73.7620, -34.6321, -41.8796, 117.543, -81.2043, 17.1718, 0, 0, 0}, {157.046, -765.472, 1735.21, -2053.86, 1371.34, -515.214, 101.081, -8.07402, 0}, {-608.740, 4827.18, -13239.6, 17742.4, -12420.0, 4369.11, -607.877, 0, 0}, {-274.278, 2380.63, -7560.19, 12582.3, -11924.5, 6464.66, -1863.44, 221.134, 0}, {-1240.72, 8096.04, -19407.0, 23942.9, -16052.3, 5559.32, -776.123, 0, 0}}, {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {-71.2528, 879.668, -3027.37, 5226.61, -4999.19, 2689.35, -761.206, 88.1242, 0}, {-1269.89, 9486.25, -26103.8, 35581.2, -25373.0, 9062.87, -1277.60, 0, 0}, {-186.640, 1811.85, -6032.01, 10283.3, -9808.11, 5285.35, -1501.87, 174.799, 0}, {-530.826, 4643.56, -13864.2, 20580.2, -15898.0, 6106.69, -916.365, 0, 0}}, {{61.6665, -90.4268, 84.5606, -37.2240, 6.09207, 0, 0, 0, 0}, {78.5134, -58.5975, 3.30928, 77.4749, -64.3984, 14.4860, 0, 0, 0}, {-1.53910, 216.936, -701.057, 1167.26, -1111.92, 615.364, -183.854, 22.8595, 0}, {-19.7415, 454.317, -1250.51, 1512.52, -762.408, 137.695, 0, 0, 0}, {-55.9612, 657.449, -2049.73, 3295.30, -2995.85, 1553.68, -427.764, 48.4324, 0}, {-522.682, 3356.77, -7535.50, 8756.49, -5518.61, 1795.60, -235.144, 0, 0}}, {{61.1474, -88.7680, 82.6446, -36.2780, 5.92310, 0, 0, 0, 0}, {73.7620, -34.6321, -41.8796, 117.543, -81.2043, 17.1718, 0, 0, 0}, {-82.0368, 883.261, -2828.84, 4621.53, -4223.56, 2185.52, -598.218, 67.2908, 0}, {608.323, -2743.56, 4942.01, -4045.58, 1558.07, -226.240, 0, 0, 0}, {4.07203, 138.882, -321.983, 282.702, -12.9566, -129.159, 74.5884, -12.9994, 0}, {-866.737, 5984.13, -15129.6, 19134.6, -12757.7, 4276.79, -566.056, 0, 0}}};

const Float_t fgPar_2GeV_2250_EfidTheta_S3[4][8] = {{74.4893, -158.720, 251.241, -200.000, 52.2984, 25.4188, -18.8692, 3.27217}, {90.8413, -226.800, 358.487, -259.260, 30.9359, 68.7248, -38.9760, 6.47933}, {117.102, -429.455, 1208.29, -1922.72, 1791.40, -965.135, 277.459, -32.8536}, {55.0676, 91.0959, -444.252, 791.284, -717.492, 350.325, -87.3235, 8.68087}};

const Float_t fgPar_2GeV_2250_EfidTheta_S4[2][8] = {{77.7940, -192.492, 361.852, -394.127, 246.499, -84.6133, 13.9182, -0.713846}, {41.2902, 110.603, -586.690, 1130.70, -1137.27, 633.345, -185.038, 22.1482}};

const Float_t fgPar_2GeV_2250_EfidTheta_S5[8][8] = {{-12998.3, 57694.0, -109085, 114102, -71303.6, 26616.9, -5494.45, 483.756}, {17842.9, -74659.9, 133869, -133170, 79380.0, -28352.3, 5617.88, -476.314}, {65.5364, -99.8689, 88.5645, 24.8299, -121.327, 102.818, -37.8275, 5.28492}, {66.4049, -76.4096, -20.8674, 230.072, -318.905, 206.721, -66.3286, 8.48753}, {100.262, -358.882, 957.267, -1495.42, 1396.73, -765.881, 226.791, -27.9341}, {50.4447, 48.3032, -315.976, 580.141, -525.583, 252.075, -59.9294, 5.34805}, {78.5845, -155.728, 320.528, -420.296, 341.899, -164.626, 42.5274, -4.50224}, {95.9430, -221.787, 391.495, -350.033, 131.391, 13.2965, -24.0460, 4.92253}};

const Float_t fgPar_2GeV_2250_Pfid_For[6][4][7] = 
  {{{60.2165, -189.720, 446.990, -523.122, 320.721, -97.8518, 11.5258}, 
    {-1457.16, 13814.2, -43182.7, 66646.0, -54355.1, 22423.5, -3683.76}, 
    {17.1086, 54.2974, -103.464, 111.325, -70.7673, 27.2551, -5.02858},
    {-2547.86, 22143.1, -66326.6, 101105.0, -82187.8, 33959.7, -5607.59}}, 
   {{65.7242, -246.922, 759.745, -1198.32, 1007.05, -428.060, 72.2644}, 
    {3384.16, -19353.1, 54083.5, -79843.4, 63870.2, -26079.2, 4250.29}, 
    {85.2489, -441.821, 1327.52, -1978.53, 1567.84, -633.530, 102.928}, 
    {411.998, -533.572, 599.925, 2099.52, -5061.48, 3701.58, -891.843}},
   {{110.022, -558.044, 1512.96, -2098.53, 1579.55, -613.478, 96.3279}, 
    {3937.29, -23745.1, 59651.0, -76988.6, 54276.0, -19900.2, 2974.95}, 
    {35.8488, -46.9595, 107.492, -93.9141, 10.5845, 26.1910, -9.89460}, 
    {-326.838, 4634.99, -11155.2, 11811.4, -5405.80, 554.030, 175.526}}, 
   {{38.9338, -62.8663, 118.218, -56.6953, -40.5083, 46.1782, -11.5822}, 
    {1864.83, -11735.6, 34175.4, -48928.5, 37315.8, -14496.1, 2254.05}, 
    {23.6892, 9.69854, 94.4521, -270.119, 288.132, -140.031, 25.9272}, 
    {-261.086, 4863.13, -11760.4, 13791.1, -8983.19, 3136.52, -457.183}}, 
   {{-11.0252, 348.901, -1172.63, 1980.73, -1759.08, 786.043, -139.299}, 
    {-2231.41, 23477.1, -78229.3, 129238.0, -111761.0, 48561.4, -8370.65}, 
    {104.415, -548.464, 1506.70, -2064.10, 1507.55, -561.677, 83.9247}, 
    {1402.87, -9008.78, 25660.0, -37543.3, 29860.8, -12238.4, 2019.03}}, 
   {{20.4577, 66.1373, -205.218, 372.864, -366.625, 177.596, -33.1168}, 
    {2059.77, -14468.3, 46492.9, -72168.2, 58275.9, -23615.8, 3800.60}, 
    {-18.9897, 392.519, -1234.31, 1950.24, -1623.01, 681.260, -113.806},
    {-3478.50, 32840.9, -104381.0, 167656.0, -143070.0, 61909.3, -10690.1}}};

const Float_t fgPar_2GeV_2250_Pfid_Bak[6][4][7] =
  {{{110.007, 121.302, 97.8380, -1679.71, 4022.73, -3973.09, 1422.42}, 
    {69.7305, 359.843, -876.383, 649.612, 600.059, -1155.43, 472.866}, 
    {13.9334, -236.587, 810.783, -1614.65, 1851.97, -1125.48, 280.069}, 
    {10.1644, 51.7943, -527.843, 2071.12, -3480.34, 2663.52, -768.498}}, 
   {{161.555, -263.801, 770.924, -902.814, 503.641, -319.619, 171.147},
    {154.660, -619.711, 3444.65, -8994.29, 12253.9, -8439.82, 2321.14}, 
    {117.461, -1429.96, 6117.79, -13492.3, 16142.2, -9965.40, 2490.47}, 
    {7.77411, -17.3501, 279.462, -876.326, 1398.82, -1137.49, 365.383}}, 
   {{-31.1460, 1942.49, -9193.97, 21731.0, -26961.3, 16701.7, -4067.85}, 
    {154.660, -654.420, 3774.08, -9920.36, 13333.7, -8953.68, 2386.32}, 
    {63.2709, -867.859, 4000.97, -9557.57, 12215.1, -7926.91, 2052.90}, 
    {-28.1127, 484.636, -2665.71, 7484.94, -10740.7, 7561.79, -2076.70}}, 
   {{172.853, -656.312, 3768.76, -10243.0, 14600.3, -10616.3, 3095.27}, 
    {270.076, -1938.46, 9276.01, -21861.1, 27363.7, -17479.9, 4490.05}, 
    {32.2327, -432.593, 1666.57, -3491.43, 4031.58, -2406.30, 579.944},
    {-44.9153, 638.112, -2971.77, 7223.13, -9328.99, 6080.46, -1576.13}}, 
   {{45.7403, 875.133, -3646.85, 7848.52, -8905.36, 4914.78, -1010.91}, 
    {138.000, -449.485, 2806.13, -7725.44, 10777.3, -7482.95, 2056.80}, 
    {72.7551, -944.002, 4200.92, -9776.76, 12316.6, -7955.78, 2066.50},
    {-9.59531, 180.519, -795.797, 2124.85, -2978.29, 2040.14, -541.811}}, 
   {{77.5100, 494.571, -1625.99, 2397.48, -1177.99, -574.604, 530.446}, 
    {117.869, -56.8761, 330.252, -715.276, 807.257, -497.124, 133.989}, 
    {7.66164, -208.001, 996.883, -2772.33, 4100.81, -3008.90, 864.126}, 
    {-25.3497, 346.501, -1458.46, 3513.62, -4625.70, 3088.01, -818.696}}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS2[2][6] = 
  {{-28.1486, 425.124, -935.693, 1065.39, -608.526, 137.658}, 
   {-15.2084, 345.466, -697.657, 751.738, -419.288, 95.2206}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS3[8][6] = 
  {{17.1490, 294.605, -640.590, 707.758, -386.730, 83.2529},
   {35.9318, 204.580, -404.489, 413.240, -209.580, 41.7819}, 
   {47.6825, 274.777, -754.725, 1117.80, -846.816, 255.607}, 
   {44.7484, 344.543, -872.200, 1113.89, -694.736, 168.061}, 
   {-205.978, 828.617, -1199.65, 875.482, -317.846, 45.6938},
   {-240.595, 961.068, -1370.34, 977.625, -345.743, 48.3834},
   {-136.104, 479.276, -593.135, 374.730, -118.350, 14.7923},
   {-196.773, 700.974, -894.540, 577.460, -185.690, 23.6201}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS4[4][6] = 
  {{81.8115, 139.810, -445.130, 804.212, -821.194, 364.924},
   {79.5053, 317.287, -1582.80, 3987.05, -4880.55, 2305.63}, 
   {-137.480, 633.288, -954.383, 721.057, -269.140, 39.4822}, 
   {-145.605, 697.662, -1088.74, 853.855, -330.883, 50.3421}};

const Float_t fgPar_2GeV_2250_Pfid_ScpdS5[8][6] =
  {{-29.9426, 370.963, -714.697, 707.343, -348.995, 67.7647}, 
   {-27.4173, 372.536, -693.341, 652.792, -302.559, 54.7761}, 
   {-47.1617, 132.967, -104.776, 41.7673, -7.68238, 0.404311}, 
   {-54.5895, 149.685, -111.590, 41.2556, -6.93943, 0.301087}, 
   {-79.1386, 275.678, -341.972, 218.907, -69.5520, 8.66381}, 
   {-97.5794, 352.616, -468.487, 322.829, -111.159, 15.0975}, 
   {22.5823, -182.064, 365.317, -294.653, 108.779, -15.2712}, 
   {-7.59521, 2.91795, 31.6773, -28.3085, 10.5943, -1.57966}};


const Float_t fgPar_4Gev_2250_Pfidft1l[6][6] = { {26.2564,0.441269,-29.7632,94.5137,7.71903,2.10915}, {29.7455,-0.826489,4.09596,91.8187,8.38108,1.5016},  					
						 {29.5399,-0.878321,43.1909,64.9772,11.1844,0.825411},{28.5857,0.4061,98.6296,95.5022,13.7297,0.415071},
						 {31.9803,0.341766,257.124,103.504,14.2357,0.43387},{29.2846,-0.257616,51.1709,84.3207,10.2963,1.69991}};
const Float_t fgPar_4Gev_2250_Pfidft1r[6][6] = { {34.7359,-1.45301,660.653,-79.1375,11.3239,1.05352}, {30.6992,0.71858,442.087,4.20897,3.62722,3.35155},
						 {19.1518,3.71404,-197.134,177.828,9.63173,1.35402}, {23.9897,1.52101,23.9288,71.4476,8.89464,1.69512},
						 {22.6619,2.4697,-54.5174,112.22,11.2561,0.687839}, {20.9859,3.86504,-56.5229,230.635,13.6587,0.270987}}; 
const Float_t fgPar_4Gev_2250_Pfidft2l[6][6] = { {24.683,0.470268,124.501,-9.04329,8.60129,1.66063}, {26.2736,-0.591497,182.954,-51.059,7.65701,2.29757},
						 {24.8681,1.15526,111.322,22.2304,9.46319,1.6834}, {29.3639,1.307,282.797,89.5863,11.7162,0.376266},
						 {36.8099,-0.785452,655.368,46.4935,12.0443,0.500522}, {25.8401,0.899645,141.723,27.6687,9.62103,1.7379}}; 
const Float_t fgPar_4Gev_2250_Pfidft2r[6][6] = { {32.9905,-0.580968,464.263,30.5379,11.7414,0.320415}, {26.8867,0.748481,150.349,51.4182,8.70942,1.51013},
						 {26.0729,0.357197,136.456,24.1839,6.70568,0.820883}, {25.8339,1.018,149.648,38.7987,6.56928,0.527773},
						 {27.997,0.0685368,268.87,-45.3343,5.26386,3.08026}, {30.3568,1.60206,359.39,197.047,11.1523,0.451219}}; 
const Float_t fgPar_4Gev_2250_Pfidbt1l[6][6] = { {-24.4118,4.20154,-0.0480933,-0.0800641,0.000311929,0.000511191}, {-34.5523,8.81812,0.221281,-0.203846,-0.00115322,0.00119883},
						 {-29.4962,6.57417,0.0830637,-0.142094,-0.000271087,0.000801481}, {-29.5177,6.23458,0.183415,-0.160458,-0.00121912,0.0010282},
						 {-19.8091,4.37431,-0.046672,-0.124147,-7.21454e-05,0.000931229}, {-38.1865,10.6462,0.363126,-0.267793,-0.00212252,0.00162732}}; 
const Float_t fgPar_4Gev_2250_Pfidbt1r[6][6] = { {-15.6987,3.34818,-0.155291,-0.102923,0.000736214,0.000775517}, {-15.9442,1.75807,-0.196246,-0.0524198,0.00118102,0.000398854},
						 {-14.4453,1.65733,-0.269699,-0.0423913,0.00187485,0.000274252}, {-18.5972,1.41622,-0.144491,-0.0369631,0.000874762,0.000326006},
						 {-17.1008,0.577868,-0.173353,-0.021315,0.00108238,0.000189545}, {2.21904,-3.38706,-0.636698,0.0953525,0.0038789,-0.000559086}}; 
const Float_t fgPar_4Gev_2250_Pfidbt2l[6][6] = { {-13.7253,-1.53789,-0.296133,0.0648705,0.00269427,-0.000928492}, {-12.356,-2.62192,-0.366191,0.115155,0.0033624,-0.00137599},
						 {-2.52638,-9.6591,-0.743505,0.380195,0.0067055,-0.00369404}, {-34.5804,15.3815,0.417723,-0.489802,-0.00337546,0.00370894},
						 {1.87747,-7.70598,-0.919924,0.376373,0.00776553,-0.00354661}, {-12.3968,-2.37408,-0.367352,0.114661,0.00352523,-0.00148841}}; 
const Float_t fgPar_4Gev_2250_Pfidbt2r[6][6] = { {-29.5895,10.9088,0.248994,-0.326966,-0.00154954,0.00202508}, {-7.20087,-6.19132,-0.568426,0.257971,0.00476513,-0.00236084},
						 {-10.0076,-3.66545,-0.468027,0.163446,0.00421363,-0.00175242}, {-9.03582,-5.14009,-0.515592,0.221044,0.00482855,-0.00237549},
						 {-8.55955,-5.27785,-0.504058,0.201472,0.00404296,-0.00175892}, {-21.122,5.19264,-0.0761427,-0.0826774,0.0018747,-0.000390706}}; 
const Float_t fgPar_4Gev_2250_Pfidbl[6][6] = { {131.839,-6.64199,-22.8623,4.91185,126.5,20}, {132.055,-5.2283,2.20945,-1.57951,128.429,11.4286},
					       {137.945,-7.90553,-12.8716,3.94534,119.857,22.8571}, {124.743,-3.54503,-22.8263,5.62231,130.429,11.4286},
					       {136.455,-7.59559,-18.6847,4.52149,123.5,20}, {126.556,-4.02284,-22.2328,5.23298,124.857,22.8571}};
const Float_t fgPar_4Gev_2250_Pfidbr[6][6] = { {97.3917,2.99764,26.7715,-5.95695,126.5,20}, {132.154,-6.60261,0.000146616,1.53542,128.429,11.4286},
					       {113.746,-1.24667,32.0728,-9.35241,119.857,22.8571}, {118.596,-2.44983,22.2973,-5.40976,130.429,11.4286},
					       {125.129,-3.96273,21.6178,-5.86908,123.5,20}, {111.201,-0.178015,25.1267,-6.55928,124.857,22.8571}}; 


class treeProducer_simulation {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Int_t           iev;
  Int_t           neu;
  Int_t           fspl;
  Int_t           tgt;
  Int_t           Z;
  Int_t           A;
  Int_t           hitnuc;
  Int_t           hitqrk;
  Int_t           resid;
  Bool_t          sea;
  Bool_t          qel;
  Bool_t          mec;
  Bool_t          res;
  Bool_t          dis;
  Bool_t          coh;
  Bool_t          dfr;
  Bool_t          imd;
  Bool_t          imdanh;
  Bool_t          singlek;
  Bool_t          nuel;
  Bool_t          em;
  Bool_t          cc;
  Bool_t          nc;
  Bool_t          charm;
  Int_t           neut_code;
  Int_t           nuance_code;
  Double_t        wght;
  Double_t        xs;
  Double_t        ys;
  Double_t        ts;
  Double_t        Q2s;
  Double_t        Ws;
  Double_t        x;
  Double_t        y;
  Double_t        t;
  Double_t        Q2;
  Double_t        W;
  Double_t        EvRF;
  Double_t        Ev;
  Double_t        pxv;
  Double_t        pyv;
  Double_t        pzv;
  Double_t        En;
  Double_t        pxn;
  Double_t        pyn;
  Double_t        pzn;
  Double_t        El;
  Double_t        pxl;
  Double_t        pyl;
  Double_t        pzl;
  Double_t        pl;
  Double_t        cthl;
  Int_t           nfp;
  Int_t           nfn;
  Int_t           nfpip;
  Int_t           nfpim;
  Int_t           nfpi0;
  Int_t           nfkp;
  Int_t           nfkm;
  Int_t           nfk0;
  Int_t           nfem;
  Int_t           nfother;
  Int_t           nip;
  Int_t           nin;
  Int_t           nipip;
  Int_t           nipim;
  Int_t           nipi0;
  Int_t           nikp;
  Int_t           nikm;
  Int_t           nik0;
  Int_t           niem;
  Int_t           niother;
  Int_t           ni;
  Int_t           pdgi[2];   //[ni]
  Int_t           resc[1];   //[ni]
  Double_t        Ei[2];   //[ni]
  Double_t        pxi[2];   //[ni]
  Double_t        pyi[2];   //[ni]
  Double_t        pzi[2];   //[ni]
  Int_t           nf;
  Int_t           pdgf[120];   //[nf]
  Double_t        Ef[120];   //[nf]
  Double_t        pxf[120];   //[nf]
  Double_t        pyf[120];   //[nf]
  Double_t        pzf[120];   //[nf]
  Double_t        pf[120];   //[nf]
  Double_t        cthf[120];   //[nf]
  Double_t        vtxx;
  Double_t        vtxy;
  Double_t        vtxz;
  Double_t        vtxt;
  Double_t        sumKEf;
  Double_t        calresp0;

  // List of branches
  TBranch        *b_iev;   //!
  TBranch        *b_neu;   //!
  TBranch        *b_fspl;   //!
  TBranch        *b_tgt;   //!
  TBranch        *b_Z;   //!
  TBranch        *b_A;   //!
  TBranch        *b_hitnuc;   //!
  TBranch        *b_hitqrk;   //!
  TBranch        *b_resid;   //!
  TBranch        *b_sea;   //!
  TBranch        *b_qel;   //!
  TBranch        *b_mec;   //!
  TBranch        *b_res;   //!
  TBranch        *b_dis;   //!
  TBranch        *b_coh;   //!
  TBranch        *b_dfr;   //!
  TBranch        *b_imd;   //!
  TBranch        *b_imdanh;   //!
  TBranch        *b_singlek;   //!
  TBranch        *b_nuel;   //!
  TBranch        *b_em;   //!
  TBranch        *b_cc;   //!
  TBranch        *b_nc;   //!
  TBranch        *b_charm;   //!
  TBranch        *b_neut_code;   //!
  TBranch        *b_nuance_code;   //!
  TBranch        *b_wght;   //!
  TBranch        *b_xs;   //!
  TBranch        *b_ys;   //!
  TBranch        *b_ts;   //!
  TBranch        *b_Q2s;   //!
  TBranch        *b_Ws;   //!
  TBranch        *b_x;   //!
  TBranch        *b_y;   //!
  TBranch        *b_t;   //!
  TBranch        *b_Q2;   //!
  TBranch        *b_W;   //!
  TBranch        *b_EvRF;   //!
  TBranch        *b_Ev;   //!
  TBranch        *b_pxv;   //!
  TBranch        *b_pyv;   //!
  TBranch        *b_pzv;   //!
  TBranch        *b_En;   //!
  TBranch        *b_pxn;   //!
  TBranch        *b_pyn;   //!
  TBranch        *b_pzn;   //!
  TBranch        *b_El;   //!
  TBranch        *b_pxl;   //!
  TBranch        *b_pyl;   //!
  TBranch        *b_pzl;   //!
  TBranch        *b_pl;   //!
  TBranch        *b_cthl;   //!
  TBranch        *b_nfp;   //!
  TBranch        *b_nfn;   //!
  TBranch        *b_nfpip;   //!
  TBranch        *b_nfpim;   //!
  TBranch        *b_nfpi0;   //!
  TBranch        *b_nfkp;   //!
  TBranch        *b_nfkm;   //!
  TBranch        *b_nfk0;   //!
  TBranch        *b_nfem;   //!
  TBranch        *b_nfother;   //!
  TBranch        *b_nip;   //!
  TBranch        *b_nin;   //!
  TBranch        *b_nipip;   //!
  TBranch        *b_nipim;   //!
  TBranch        *b_nipi0;   //!
  TBranch        *b_nikp;   //!
  TBranch        *b_nikm;   //!
  TBranch        *b_nik0;   //!
  TBranch        *b_niem;   //!
  TBranch        *b_niother;   //!
  TBranch        *b_ni;   //!
  TBranch        *b_pdgi;   //!
  TBranch        *b_resc;   //!
  TBranch        *b_Ei;   //!
  TBranch        *b_pxi;   //!
  TBranch        *b_pyi;   //!
  TBranch        *b_pzi;   //!
  TBranch        *b_nf;   //!
  TBranch        *b_pdgf;   //!
  TBranch        *b_Ef;   //!
  TBranch        *b_pxf;   //!
  TBranch        *b_pyf;   //!
  TBranch        *b_pzf;   //!
  TBranch        *b_pf;   //!
  TBranch        *b_cthf;   //!
  TBranch        *b_vtxx;   //!
  TBranch        *b_vtxy;   //!
  TBranch        *b_vtxz;   //!
  TBranch        *b_vtxt;   //!
  TBranch        *b_sumKEf;   //!
  TBranch        *b_calresp0;   //!

  treeProducer_simulation(TTree *tree=0);
  virtual ~treeProducer_simulation();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  virtual void SetFiducialCutParameters();
  virtual bool ElectronFiducialCut(TVector3* momentum);
  virtual bool GetEPhiLimits(Float_t momentum, Float_t theta, Int_t sector,Float_t *EPhiMin, Float_t *EPhiMax);
  virtual bool ProtonFiducialCut(TVector3* momentum, Float_t *ptr_cphil, Float_t *ptr_cphir);
  virtual bool PimiFiducialCut(TVector3* momentum, Float_t *pimi_philow, Float_t *pimi_phiup);
  virtual void prot_rot_func_2p(TVector3 V3q, TVector3  V3prot[2],TLorentzVector V4el,double Ecal_2pto1p[2],double  pmiss_perp_2pto1p[2],double  P2pto1p[2], double *Nboth);
  virtual void  prot_rot_func_3p(TVector3 V3q, TVector3  V3prot[3],TLorentzVector V4el,double Ecal_3pto2p[][2],double  pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3],double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det);

};

#endif

#ifdef treeProducer_simulation_cxx
treeProducer_simulation::treeProducer_simulation(TTree *tree) : fChain(0) 
{

  TFile * mcfile = new TFile( file_name,"ROOT");
  tree = (TTree*)mcfile->Get("gst");

  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {

    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(FullInputName.c_str());
    if (!f || !f->IsOpen()) f = new TFile(FullInputName.c_str());
    std::cout<<FullInputName<<std::endl;
    f->GetObject("gst",tree);

  }
  Init(tree);
}

treeProducer_simulation::~treeProducer_simulation()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t treeProducer_simulation::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t treeProducer_simulation::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void treeProducer_simulation::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("iev", &iev, &b_iev);
  fChain->SetBranchAddress("neu", &neu, &b_neu);
  fChain->SetBranchAddress("fspl", &fspl, &b_fspl);
  fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
  fChain->SetBranchAddress("Z", &Z, &b_Z);
  fChain->SetBranchAddress("A", &A, &b_A);
  fChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
  fChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
  fChain->SetBranchAddress("resid", &resid, &b_resid);
  fChain->SetBranchAddress("sea", &sea, &b_sea);
  fChain->SetBranchAddress("qel", &qel, &b_qel);
  fChain->SetBranchAddress("mec", &mec, &b_mec);//go down from here
  fChain->SetBranchAddress("res", &res, &b_res);
  fChain->SetBranchAddress("dis", &dis, &b_dis);
  fChain->SetBranchAddress("coh", &coh, &b_coh);
  fChain->SetBranchAddress("dfr", &dfr, &b_dfr);
  fChain->SetBranchAddress("imd", &imd, &b_imd);
  fChain->SetBranchAddress("imdanh", &imdanh, &b_imdanh);
  fChain->SetBranchAddress("singlek", &singlek, &b_singlek);
  fChain->SetBranchAddress("nuel", &nuel, &b_nuel);
  fChain->SetBranchAddress("em", &em, &b_em);
  fChain->SetBranchAddress("cc", &cc, &b_cc);
  fChain->SetBranchAddress("nc", &nc, &b_nc);
  fChain->SetBranchAddress("charm", &charm, &b_charm);
  fChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
  fChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
  fChain->SetBranchAddress("wght", &wght, &b_wght);
  fChain->SetBranchAddress("xs", &xs, &b_xs);
  fChain->SetBranchAddress("ys", &ys, &b_ys);
  fChain->SetBranchAddress("ts", &ts, &b_ts);
  fChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
  fChain->SetBranchAddress("Ws", &Ws, &b_Ws);
  fChain->SetBranchAddress("x", &x, &b_x);
  fChain->SetBranchAddress("y", &y, &b_y);
  fChain->SetBranchAddress("t", &t, &b_t);
  fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
  fChain->SetBranchAddress("W", &W, &b_W);
  fChain->SetBranchAddress("EvRF", &EvRF, &b_EvRF);
  fChain->SetBranchAddress("Ev", &Ev, &b_Ev);
  fChain->SetBranchAddress("pxv", &pxv, &b_pxv);
  fChain->SetBranchAddress("pyv", &pyv, &b_pyv);
  fChain->SetBranchAddress("pzv", &pzv, &b_pzv);
  fChain->SetBranchAddress("En", &En, &b_En);
  fChain->SetBranchAddress("pxn", &pxn, &b_pxn);
  fChain->SetBranchAddress("pyn", &pyn, &b_pyn);
  fChain->SetBranchAddress("pzn", &pzn, &b_pzn);
  fChain->SetBranchAddress("El", &El, &b_El);
  fChain->SetBranchAddress("pxl", &pxl, &b_pxl);
  fChain->SetBranchAddress("pyl", &pyl, &b_pyl);
  fChain->SetBranchAddress("pzl", &pzl, &b_pzl);
  fChain->SetBranchAddress("pl", &pl, &b_pl);
  fChain->SetBranchAddress("cthl", &cthl, &b_cthl);
  fChain->SetBranchAddress("nfp", &nfp, &b_nfp);
  fChain->SetBranchAddress("nfn", &nfn, &b_nfn);
  fChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
  fChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
  fChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
  fChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
  fChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
  fChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
  fChain->SetBranchAddress("nfem", &nfem, &b_nfem);
  fChain->SetBranchAddress("nfother", &nfother, &b_nfother);
  fChain->SetBranchAddress("nip", &nip, &b_nip);
  fChain->SetBranchAddress("nin", &nin, &b_nin);
  fChain->SetBranchAddress("nipip", &nipip, &b_nipip);
  fChain->SetBranchAddress("nipim", &nipim, &b_nipim);
  fChain->SetBranchAddress("nipi0", &nipi0, &b_nipi0);
  fChain->SetBranchAddress("nikp", &nikp, &b_nikp);
  fChain->SetBranchAddress("nikm", &nikm, &b_nikm);
  fChain->SetBranchAddress("nik0", &nik0, &b_nik0);
  fChain->SetBranchAddress("niem", &niem, &b_niem);
  fChain->SetBranchAddress("niother", &niother, &b_niother);
  fChain->SetBranchAddress("ni", &ni, &b_ni);
  fChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
  fChain->SetBranchAddress("resc", resc, &b_resc);
  fChain->SetBranchAddress("Ei", Ei, &b_Ei);
  fChain->SetBranchAddress("pxi", pxi, &b_pxi);
  fChain->SetBranchAddress("pyi", pyi, &b_pyi);
  fChain->SetBranchAddress("pzi", pzi, &b_pzi);
  fChain->SetBranchAddress("nf", &nf, &b_nf);
  fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
  fChain->SetBranchAddress("Ef", Ef, &b_Ef);
  fChain->SetBranchAddress("pxf", pxf, &b_pxf);
  fChain->SetBranchAddress("pyf", pyf, &b_pyf);
  fChain->SetBranchAddress("pzf", pzf, &b_pzf);
  fChain->SetBranchAddress("pf", pf, &b_pf);
  fChain->SetBranchAddress("cthf", cthf, &b_cthf);
  fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
  fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
  fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
  fChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
  fChain->SetBranchAddress("sumKEf", &sumKEf, &b_sumKEf);
  fChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);
  Notify();
}

Bool_t treeProducer_simulation::Notify()
{
  return kTRUE;
}

void treeProducer_simulation::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t treeProducer_simulation::Cut(Long64_t entry)
{
  return 1;
}
//_____________________________________________________________________________________________________________________________________________________________________________________________________

void treeProducer_simulation::SetFiducialCutParameters() {

  double Ebeam = -99.;

  if (E == "2_261") { Ebeam = 2.261; }
  if (E == "4_461") { Ebeam = 4.461; }

  if (Ebeam > 4. && Ebeam < 5.) {

    FILE *fiducial_par; 
    Int_t ptype;
    Int_t ci;   
    Float_t par[6];
    Char_t Filename[100], ParLocation[80]; 
    sprintf(ParLocation,".");

    //printf("Reading fiducial cut parameters for 4.4GeV/2250A ...\n");
    sprintf(Filename,"%s/FCP_4GeV.par", ParLocation);
    fiducial_par = fopen(Filename,"r");

    while(!feof(fiducial_par)) {

      if(fscanf(fiducial_par, "%i   %i   %f   %f   %f   %f  %f   %f", &ptype, &ci, &par[0], &par[1], &par[2], &par[3], &par[4], &par[5])) {};

      switch (ptype) {

      case  0: for(int k=0; k<2; k++) fgPar_4Gev_2250_Efid_t0_p[ci-1][k] = par[k];
	break; 
      case  1: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_t1_p[ci-1][k] = par[k];
	break;
      case 10: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[ci-1][0][k] = par[k];
	break;
      case 11: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[ci-1][1][k] = par[k];
	break;
      case 20: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[ci-1][0][k] = par[k];
	break;
      case 21: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[ci-1][1][k] = par[k];	
	break;
      default: printf("error!\n"); 

      }

    }

    fclose(fiducial_par);

  } else printf("Don't know how to set fiducial cut parameters for %3.1f GeV!\n", Ebeam);

}
// __________________________________________________________________________________________________________________________________________________________________________________________________

bool treeProducer_simulation::ElectronFiducialCut(TVector3* momentum) {

  double Ebeam = -99., fTorusCurrent = -99.;

  // Torus Current & Beam Energy

  if (E == "2_261") { 

    Ebeam = 2.261, fTorusCurrent = 2250.; 

    // Electron fiducial cut, return kTRUE if pass or kFALSE if not

    Bool_t status = kTRUE;
    if (Ebeam > 2. && Ebeam < 3. && fTorusCurrent > 2240. && fTorusCurrent < 2260.) {

      Float_t phi = momentum->Phi()*180./TMath::Pi() + 180 - 30.; 
      if(phi<-30.) phi+=360.;
      Int_t sector = (Int_t)((phi+30.)/60.); 
      if(sector<0) sector=0; 
      if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum->Theta()*180./TMath::Pi();
      Float_t mom = momentum->Mag();
      Float_t par[6]; // six parameters to determine the outline of Theta vs Phi

      for (Int_t i=0; i<6; i++) {

	par[i] = 0;

	for (Int_t d=8; d>=0; d--) {

	  par[i] = par[i]*mom + fgPar_2GeV_2250_Efid[sector][i][d];

	}// calculate the parameters using pol8

      }

      if (phi < 0) {

	Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
	status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);

      } else {

	Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
	status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);

      }

      // by now, we have checked if the electron is within the outline of theta vs phi plot

      if (SCpdcut) { // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

	if (status) {

	  Int_t tsector = sector + 1;

	  if (tsector == 3){ // sector 3 has two bad paddles

	    Float_t badpar3[4];// 4 parameters to determine the positions of the two theta gaps

	    for (Int_t i=0; i<4; i++){

	      badpar3[i] = 0;

	      for (Int_t d=7; d>=0; d--) {

		badpar3[i] = badpar3[i]*mom + fgPar_2GeV_2250_EfidTheta_S3[i][d];

	      } // calculate the parameters using pol7
	    }

	    for(Int_t ipar=0;ipar<2;ipar++)
	      status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);

	  } else if (tsector == 4) { // sector 4 has one bad paddle

	    Float_t badpar4[2]; // 2 parameters to determine the position of the theta gap

	    for (Int_t i=0; i<2; i++) {
	
	      badpar4[i] = 0;
		
	      for (Int_t d=7; d>=0; d--) {

		badpar4[i] = badpar4[i]*mom + fgPar_2GeV_2250_EfidTheta_S4[i][d];

	      } // calculate the parameters using pol7

	    }

	    status = !(theta>badpar4[0] && theta<badpar4[1]);

	  } else if (tsector == 5) { // sector 5 has four bad paddles

	    Float_t badpar5[8]; // 8 parameters to determine the positions of the four theta gaps

	    for (Int_t i=0; i<8; i++) {

	      badpar5[i] = 0;

	      for (Int_t d=7; d>=0; d--) {

		badpar5[i] = badpar5[i]*mom + fgPar_2GeV_2250_EfidTheta_S5[i][d];

	      } // calculate the parameters using pol7

	    }

	    if (mom<1.25) badpar5[0] = 23.4;
	    if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	    for(Int_t ipar=0;ipar<4;ipar++)
	      status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);

	  }

	}

      }

    }

    return status;

  }

  if (E == "4_461") { 

    Ebeam = 4.461, fTorusCurrent = 2250.; 

    Bool_t status = kTRUE;  
    Float_t phiMin, phiMax;
    Float_t mom = momentum->Mag();
    Float_t phi = momentum->Phi()*180./TMath::Pi() + 180 - 30;
    if(phi<-30.) phi += 360.;
    Float_t theta = momentum->Theta()*180./TMath::Pi();
    Int_t  sector = (Int_t)((phi+30.)/60.); 	
    if(sector < 0) sector = 0; 
    if(sector > 5) sector = 5; 

    // all the work is now done in GetEPhiLimits
	
    status = GetEPhiLimits(mom, theta, sector, &phiMin, &phiMax);

    return status;

  }

  return kTRUE;

}
// __________________________________________________________________________________________________________________________________________________________________________________________________

bool treeProducer_simulation::GetEPhiLimits(Float_t momentum, Float_t theta, Int_t sector,Float_t *EPhiMin, Float_t *EPhiMax) { 

  double Ebeam = -99., fTorusCurrent = -99.;

  // Torus Current & Beam Energy

  if (E == "2_261") { 

    Ebeam = 2.261, fTorusCurrent = 2250.; 

    // Information for electron fiducial cut, 
    // returns the minimum and maximum phi accepted for a given momentum, theta and sector
    // momentum is in GeV/c, theta is in degrees, 0 <= sector <= 5

    if (sector < 0 || sector > 5) return kFALSE; // bad input
	  
    if(Ebeam>4. && Ebeam<5. && fTorusCurrent>2240. && fTorusCurrent<2260.) {

      if ((theta < 15.) || (momentum < 0.9)) return kFALSE; // out of range
      Float_t t0, t1, b[2], a[2];
      if (momentum > 3.7) momentum = 3.7; // don't extrapolate past the data
	   
      // calculates parameters of cut functions for this energy
      t0 = fgPar_4Gev_2250_Efid_t0_p[sector][0]/pow(momentum, fgPar_4Gev_2250_Efid_t0_p[sector][1]);
      t1 = 0.; for(int k=0; k<6; k++) t1 += (fgPar_4Gev_2250_Efid_t1_p[sector][k]*pow(momentum, k)); 
	
      for (int l=0; l<2; l++) {

	b[l] = 0.; for(int k=0; k<6; k++) b[l] += (fgPar_4Gev_2250_Efid_b_p[sector][l][k]*pow(momentum, k)); 
	a[l] = 0.; for(int k=0; k<6; k++) a[l] += (fgPar_4Gev_2250_Efid_a_p[sector][l][k]*pow(momentum, k));

      }
	   
      // adjust upper limit according to hardware
      if(t1 < 45.) t1 = 45.; 
      if(t0 < theta && theta < t1) { 

	*EPhiMin = 60.*sector - b[0]*(1. - 1/((theta - t0)/(b[0]/a[0]) + 1.));
	*EPhiMax = 60.*sector + b[1]*(1. - 1/((theta - t0)/(b[1]/a[1]) + 1.));

      } else {

	*EPhiMin = 60.*sector;
	*EPhiMax = 60.*sector;

      }

    } else {

      return kFALSE; // wrong beam energy/torus
    }

    return kTRUE;

  }

  if (E == "4_461") { 

    Ebeam = 4.461, fTorusCurrent = 2250.; 

    // Information for electron fiducial cut, 
    // returns the minimum and maximum phi accepted for a given momentum, theta and sector
    // momentum is in GeV/c, theta is in degrees, 0 <= sector <= 5

    if (sector < 0 || sector > 5) return kFALSE; // bad input
	  
    if(Ebeam>4. && Ebeam<5. && fTorusCurrent>2240. && fTorusCurrent<2260.) {

      if ((theta < 15.) || (momentum < 0.9)) return kFALSE; // out of range
      Float_t t0, t1, b[2], a[2];
      if (momentum > 3.7) momentum = 3.7; // don't extrapolate past the data
	   
      // calculates parameters of cut functions for this energy
      t0 = fgPar_4Gev_2250_Efid_t0_p[sector][0]/pow(momentum, fgPar_4Gev_2250_Efid_t0_p[sector][1]);
      t1 = 0.; for(int k=0; k<6; k++) t1 += (fgPar_4Gev_2250_Efid_t1_p[sector][k]*pow(momentum, k)); 
	
      for (int l=0; l<2; l++) {

	b[l] = 0.; for(int k=0; k<6; k++) b[l] += (fgPar_4Gev_2250_Efid_b_p[sector][l][k]*pow(momentum, k)); 
	a[l] = 0.; for(int k=0; k<6; k++) a[l] += (fgPar_4Gev_2250_Efid_a_p[sector][l][k]*pow(momentum, k));

      }
	   
      // adjust upper limit according to hardware
      if(t1 < 45.) t1 = 45.; 
      if(t0 < theta && theta < t1) { 

	*EPhiMin = 60.*sector - b[0]*(1. - 1/((theta - t0)/(b[0]/a[0]) + 1.));
	*EPhiMax = 60.*sector + b[1]*(1. - 1/((theta - t0)/(b[1]/a[1]) + 1.));

      } else {

	*EPhiMin = 60.*sector;
	*EPhiMax = 60.*sector;

      }

    } else {

      return kFALSE; // wrong beam energy/torus
    }

    return kTRUE;

  }

  return kTRUE;

}

// __________________________________________________________________________________________________________________________________________________________________________________________________

bool treeProducer_simulation::ProtonFiducialCut(TVector3* momentum, Float_t *ptr_cphil, Float_t *ptr_cphir) { // Proton fiducial cuts

  double Ebeam = -99., fTorusCurrent = -99.;

  // Torus Current & Beam Energy

  if (E == "2_261") { 

    Ebeam = 2.261, fTorusCurrent = 2250.; 

    Bool_t status = kTRUE;

    if (Ebeam > 2. && Ebeam < 3. && fTorusCurrent > 2240. && fTorusCurrent < 2260.) {

      Float_t phi = momentum->Phi()*180/TMath::Pi() + 180 - 30; 
      if(phi<-30) phi+=360;
      Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum->Theta()*180/TMath::Pi();
      Float_t p = momentum->Mag();
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

      for (Int_t i=0; i<4; i++) {

	par_for[i] = 0; par_bak[i] = 0;
				
	for (Int_t d=6; d>=0; d--) {

	  par_for[i] = par_for[i]*mom_for + fgPar_2GeV_2250_Pfid_For[sector][i][d];
	  par_bak[i] = par_bak[i]*mom_bak + fgPar_2GeV_2250_Pfid_Bak[sector][i][d];

	}

      }

      if (phi < 0) {

	Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
	status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
			
      } else {

	Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
	status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);

      }                     // now the forward constrains are checked

      if ( status ) {       // now check the backward constrains

	if(theta>par_bak[0]) status = kFALSE;
	else if(theta>par_bak[1]) 
	  { status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && 
	      (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]); }

      }
		
      if (status && SCpdcut) { // cut bad scintillator paddles
		 
	Int_t tsector = sector + 1;
	Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
	if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
	if(tsector==2){      // sector 2 has one bad paddle
	  Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
	  for (Int_t i=0; i<2; i++){
	    badpar2[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar2[i] = badpar2[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS2[i][d];
	    }                // calculate the parameters using pol5
	  }
	  status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
	}
	else if(tsector==3){ // sector 3 has four bad paddles
	  Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar3[i] = 0;
	    for (Int_t d=5; d>=0; d--){
	      badpar3[i] = badpar3[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS3[i][d];
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
	      badpar4[i] = badpar4[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS4[i][d];
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
	      badpar5[i] = badpar5[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS5[i][d];
	    }                // calculate the parameters using pol5
	  }
	  for (Int_t ipar=0;ipar<4;ipar++){
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	  }
	}
      }
    }
	
    return status;

  }


  if (E == "4_461") { 

    Ebeam = 4.461, fTorusCurrent = 2250.; 

    bool status = kTRUE;

    //float cphil = 0, cphir = 0, pimi_phimin = 0, pimi_phimax = 0;

    if (Ebeam > 4. && Ebeam < 5. && fTorusCurrent > 2240. && fTorusCurrent < 2260.) { // 4 GeV Fiducial Cut

      Float_t phi = momentum->Phi()*180/TMath::Pi() + 180 - 30;
      if (phi < - 30) phi += 360;
      Int_t sector = Int_t ((phi+30)/60); if (sector<0) sector=0; if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum->Theta()*180/TMath::Pi();
      Float_t p = momentum->Mag();

      Float_t parfidl[3]; for(Int_t i = 0; i < 3; i++){parfidl[i] = 0;} 
      Float_t parfidr[3]; for(Int_t i = 0; i < 3; i++){parfidr[i] = 0;}
      Float_t parfidbl[2]; for(Int_t i = 0; i < 2; i++){parfidbl[i] = 0;}
      Float_t parfidbr[2]; for(Int_t i = 0; i < 2; i++){parfidbr[i] = 0;}
      Float_t cphil = 0; Float_t cphir = 0;
      Float_t phi45l = 0; Float_t phi45r = 0;
      Float_t phi60l = 0; Float_t phi60r = 0;
      Float_t theta_min = 11;

      bool Forward = kFALSE; // Defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
      Int_t thetab = 45; // This variable defines the edge point for Forward<->Backward regions
      Float_t p1 = 0.575; // Last bin momentum for region p<0.6 GeV/c
      Float_t theta_max = 140;
      if (p < 0.2) p = 0.2; // Momentum less than 0.2 GeV/c, use 0.2 GeV/c
      if (p > 4.4) p = 4.4; // Momentum greater than 4.4 GeV/c, use 4.4 GeV/c

      //get parametrized values of theta_max for p<0.6 GeV/c region
      if (p < 0.6) { theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p; }
      //get parametrized values of theta_max for p>0.6 GeV/c region
      else {theta_max = fgPar_4Gev_2250_Pfidbl[sector][4] + fgPar_4Gev_2250_Pfidbl[sector][5]*p1; } 
			 
      //Get the momentum dependent parameters for Forward Region (theta < 45 deg)   
      Forward = kTRUE;
      if (p < 0.6) { //forward1 defines  regions of momenta p<0.6 GeV/c
	//parameters for hyperbolic function

	for (Int_t i = 0; i < 3; i++) {

	  Int_t j = 2*i;
	  parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
	  parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;

	}

      } else { //forward2 defines  regions of momenta and p>0.6 GeV/c

	for (Int_t i=0; i<3; i++) {

	  Int_t j=2*i;
	  parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
	  parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;

	}
      }

      phi45l = parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); // Parametrized value of phi at theta = 45 deg.
      phi45r = -parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));

      if(theta > thetab) { // backward region defined by theta > 45 deg. 

	if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
	if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

	//Get the momentum dependent parameters for Backward Region
		  
	Forward=kFALSE;
	if(p<0.6) { //backward1 defines  regions of momenta p<0.6 GeV/c

	  //parameters for quadratic function

	  for (Int_t i=0; i<3; i++) {

	    Int_t j=2*i;
	    parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
	    parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;

	  }

	  //these parameters determine theta_flat and phi_edge at p<0.6 GeV/c

	  for (Int_t i=0; i<2; i++) {

	    Int_t j=2*i;
	    parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
	    parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;

	  }

	} else { //backward2 defines  regions of momenta p>0.6 GeV/c

	  //parameters for quadratic function

	  for (Int_t i=0; i<3; i++) {

	    Int_t j=2*i;
	    parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
	    parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;

	  }

	  //these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum

	  for (Int_t i=0; i<2; i++) {

	    Int_t j=2*i;
	    parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
	    parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;

	  }

	}

      }
	  
      if(Forward) { //Forward region

	if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.   
	cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
	cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));

      } else { //Backward region

	phi60l = parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.; // Parametrized value of phi at theta = 60 deg.
	phi60r = -(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

	if (theta < 60) {

	  cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
	  cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);

	}  

	Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0; 
	//dr and er are theta_flat and phi_edge parameters for phi>0;  
	dl = parfidbl[0]; el = parfidbl[1];
	dr = parfidbr[0]; er = parfidbr[1];

	if (theta > 45 && theta < 60) { // BackwardA region

					// Try to match parametrized values from Forward region to Backward region parameters
	  if (cphil > phi45l) cphil = phi45l; 
	  if (cphir < phi45r) cphir = phi45r;

	}

	// BackwardB region & phi<0
	else if (theta >= 60 && theta <= dl) { cphil = phi60l; } // phi = constant 
	else if (theta > dl && theta <= theta_max) { cphil = (140-theta)*(phi60l-el) / (140-dl) + el;} // phi = stright line 
	else if (theta > theta_max) { cphil = 0; } // cut out if theta > theta_max

	//BackwardB region & phi>0
	if(theta >= 60 && theta <= dr) { cphir = phi60r; } // phi=constant 

	else if (theta > dr && theta <= theta_max) { cphir=(140-theta)*(phi60r-er)/(140-dr) + er;} // phi=stright line 
	else if (theta > theta_max) { cphir = 0; } // cut out if theta > theta_max

      } // Backward Region

      if (phi < 0) status = (phi > cphil); // Check the constrains 
      else if (phi >= 0) { status = (phi < cphir); }

      if (theta < theta_min) status = kFALSE; // Cutting out events below theta_min

      if(Forward && p<0.6 && theta<20.6-11.4*p) status = kFALSE; // Function defines cut of the edge at low theta for p<0.6 GeV/c

      //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for 
      //some range of momentum, where edge does not look good.
      bool s1s4 = (theta < 11.7 && (sector == 0 || sector == 3));
      bool s5 = (theta < 12.2 && sector == 4);
      bool s6 = (theta < 11.4 && sector == 5);
      if (p >= 0.6 && p < 1.5 && (s1s4 || s5 || s6)) status = kFALSE; 

      *ptr_cphil = cphil;
      *ptr_cphir = cphir;

    }

    return status;

  }

  return kTRUE;

}
// __________________________________________________________________________________________________________________________________________________________________________________________________

Bool_t treeProducer_simulation::PimiFiducialCut(TVector3* momentum, Float_t *pimi_philow, Float_t *pimi_phiup) {

  double Ebeam = -99., fTorusCurrent = -99.;

  // Torus Current & Beam Energy

  if (E == "2_261") { 

    Ebeam = 2.261, fTorusCurrent = 2250.; 

    // Electron fiducial cut, return kTRUE if pass or kFALSE if not
    Bool_t status = kTRUE;

    if ( Ebeam > 2. &&Ebeam<3. && fTorusCurrent>2240. && fTorusCurrent<2260.) {

      Float_t phi = momentum->Phi()*180./TMath::Pi() - 30.; 
      if (phi < -30.) phi += 360.;
      Int_t sector = (Int_t)((phi+30.)/60.); 
      if(sector<0)sector=0; 
      if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum->Theta()*180./TMath::Pi();
      if (theta>45) theta = 45;   //to extrapolate the cut to higher theta for pimi
      Float_t mom = momentum->Mag(), phimin, phimax;
      Float_t par[6]; // six parameters to determine the outline of Theta vs Phi

      for (Int_t i=0; i<6; i++) {

	par[i] = 0;

	for (Int_t d=8; d>=0; d--) {

	  par[i] = par[i]*mom + fgPar_2GeV_2250_Efid[sector][i][d];

	} // Calculate the parameters using pol8

      }

      if (phi < 0) {

	Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
	phimin =  par[3]/((theta-par[0])+par[3]/par[2])-par[2];
	phimax =  par[2]-par[3]/((theta-par[0])+par[3]/par[2]);
	*pimi_philow = phimin;
	*pimi_phiup = phimax;
	status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);
				
      } else {

	Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
	phimin =  par[5]/((theta-par[0])+par[5]/par[4])-par[4];
	phimax =  par[4]-par[5]/((theta-par[0])+par[5]/par[4]);
	*pimi_philow = phimin;
	*pimi_phiup = phimax;
	status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);

      }

      // by now, we have checked if the electron is within the outline of theta vs phi plot

      if (SCpdcut) {  // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

	if (status) {

	  Int_t tsector = sector + 1;

	  if (tsector == 3) { // sector 3 has two bad paddles

	    Float_t badpar3[4];            // 4 parameters to determine the positions of the two theta gaps

	    for (Int_t i=0; i<4; i++) {

	      badpar3[i] = 0;

	      for (Int_t d=7; d>=0; d--) {

		badpar3[i] = badpar3[i]*mom + fgPar_2GeV_2250_EfidTheta_S3[i][d];

	      } // Calculate the parameters using pol7

	    }

	    for (Int_t ipar=0;ipar<2;ipar++) { status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]); }

	  } else if (tsector == 4) { // sector 4 has one bad paddle

	    Float_t badpar4[2]; // 2 parameters to determine the position of the theta gap

	    for (Int_t i=0; i<2; i++) {

	      badpar4[i] = 0;

	      for (Int_t d=7; d>=0; d--) {

		badpar4[i] = badpar4[i]*mom + fgPar_2GeV_2250_EfidTheta_S4[i][d];

	      } // calculate the parameters using pol7

	    }

	    status = !(theta>badpar4[0] && theta<badpar4[1]);

	  } else if (tsector == 5) { // sector 5 has four bad paddles

	    Float_t badpar5[8]; // 8 parameters to determine the positions of the four theta gaps

	    for (Int_t i=0; i<8; i++) {

	      badpar5[i] = 0;

	      for (Int_t d=7; d>=0; d--) {

		badpar5[i] = badpar5[i]*mom + fgPar_2GeV_2250_EfidTheta_S5[i][d];

	      } // calculate the parameters using pol7

	    }

	    if (mom<1.25) badpar5[0] = 23.4;
	    if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	    for(Int_t ipar=0;ipar<4;ipar++) { status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]); }

	  }

	}

      }

    }

    return status;

  }

  if (E == "4_461") { 

    Ebeam = 4.461, fTorusCurrent = 2250.; 

    // Pion fiducial cut, return kTRUE if pass or kFALSE if not
    bool status = kTRUE;

    if (fTorusCurrent > 2240. && fTorusCurrent < 2260.) {

      Float_t phi=momentum->Phi()*180./TMath::Pi() + 180 - 30; 
      if(phi<-30.) phi+=360.;
      Int_t sector = (Int_t)((phi+30.)/60.); 
      if(sector<0)sector=0; 
      if(sector>5) sector=5;
      phi -= sector*60;
      Float_t theta = momentum->Theta()*180./TMath::Pi();
      if (theta>45)theta=45;   //to extrapolate the cut to higher theta for pimi
      Float_t mom = momentum->Mag(), phimin, phimax;
      Float_t par[6];   // six parameters to determine the outline of Theta vs Phi

      for (Int_t i=0; i<6; i++) {

	par[i] = 0;

	for (Int_t d=8; d>=0; d--) {

	  par[i] = par[i]*mom + fgPar_2GeV_2250_Efid[sector][i][d];
			  
	}  // calculate the parameters using pol8

      }

      if (phi < 0) {

	Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
	phimin =  par[3]/((theta-par[0])+par[3]/par[2])-par[2];
	phimax =  par[2]-par[3]/((theta-par[0])+par[3]/par[2]);
	*pimi_philow = phimin;
	*pimi_phiup = phimax;
	status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);

      } else {

	Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
	phimin =  par[5]/((theta-par[0])+par[5]/par[4])-par[4];
	phimax =  par[4]-par[5]/((theta-par[0])+par[5]/par[4]);
	*pimi_philow = phimin;
	*pimi_phiup = phimax;
	status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);

      }

      // by now, we have checked if the electron is within the outline of theta vs phi plot

      if (SCpdcut) {  // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

	if (status) {

	  Int_t tsector = sector + 1;
	  if (tsector == 3) {   // sector 3 has two bad paddles

	    Float_t badpar3[4];// 4 parameters to determine the positions of the two theta gaps

	    for (Int_t i=0; i<4; i++) {

	      badpar3[i] = 0;

	      for (Int_t d=7; d>=0; d--) {

		badpar3[i] = badpar3[i]*mom + fgPar_2GeV_2250_EfidTheta_S3[i][d];

	      } // calculate the parameters using pol7
		
	    }

	    for (Int_t ipar=0;ipar<2;ipar++)

	      status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);

	  } else if (tsector == 4) { // sector 4 has one bad paddle

	    Float_t badpar4[2];   // 2 parameters to determine the position of the theta gap
			
	    for (Int_t i=0; i<2; i++) {

	      badpar4[i] = 0;
			
	      for (Int_t d=7; d>=0; d--) {

		badpar4[i] = badpar4[i]*mom + fgPar_2GeV_2250_EfidTheta_S4[i][d];

	      }   // calculate the parameters using pol7

	    }

	    status = !(theta>badpar4[0] && theta<badpar4[1]);

	  } else if (tsector == 5) { // sector 5 has four bad paddles

	    Float_t badpar5[8];   // 8 parameters to determine the positions of the four theta gaps
	
	    for (Int_t i=0; i<8; i++) {

	      badpar5[i] = 0;

	      for (Int_t d=7; d>=0; d--) {

		badpar5[i] = badpar5[i]*mom + fgPar_2GeV_2250_EfidTheta_S5[i][d];

	      }   // calculate the parameters using pol7

	    }

	    if (mom<1.25) badpar5[0] = 23.4;
	    if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	    for(Int_t ipar=0;ipar<4;ipar++)
			
	      status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
		
	  }
			
	}

      }

    }

    return status;

  } 

  return kTRUE;

}

// __________________________________________________________________________________________________________________________________________________________________________________________________

void treeProducer_simulation::prot_rot_func_2p(TVector3 V3q, TVector3  V3prot[2],TLorentzVector V4el,double Ecal_2pto1p[2],double  pmiss_perp_2pto1p[2],double  P2pto1p[2], double *Nboth) {

  double BE = -99.;

  if (Target == "3He") { BE = He3_bind_en-D2_bind_en; }
  if (Target == "4He") { BE = He4_bind_en-H3_bind_en; }
  if (Target == "12C") { BE = C12_bind_en - B_bind_en; }
  if (Target == "56Fe") { BE = Fe_bind_en-Mn_bind_en; }

  const int N2 = 2;                                 
  double rot_angle, N_p2to1[N2]={0},Ntot=100,m_prot=0.9382720813;
  TVector3 V3_2prot[N2], V3_prot_el_2pto1p[N2];
  double N_2=0;

  float cphil = 0, cphir = 0, pimi_phimin = 0, pimi_phimax = 0;  
     
  for (int g1=0; g1<Ntot; g1++) {

    rot_angle=gRandom->Uniform(0,2*TMath::Pi());

    V3_2prot[0]=V3prot[0];
    V3_2prot[1]=V3prot[1];
    V3_2prot[0].Rotate(rot_angle,V3q);
    V3_2prot[1].Rotate(rot_angle,V3q); 

    if(ProtonFiducialCut(&V3_2prot[0], &cphil, &cphir)  && !ProtonFiducialCut(&V3_2prot[1], &cphil, &cphir)) N_p2to1[0]=N_p2to1[0]+1;
    if(!ProtonFiducialCut(&V3_2prot[0], &cphil, &cphir) && ProtonFiducialCut(&V3_2prot[1], &cphil, &cphir))  N_p2to1[1]=N_p2to1[1]+1;
    if(ProtonFiducialCut(&V3_2prot[0], &cphil, &cphir)  && ProtonFiducialCut(&V3_2prot[1], &cphil, &cphir))  N_2=N_2+1;

  }

  //-----------------------------------------  2p to 1p  -----------------------------------------------------------------------

  if (N_2 != 0) {
    
    V3_prot_el_2pto1p[0]=V4el.Vect()+ V3prot[0];
    Ecal_2pto1p[0]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[0].Mag()*V3prot[0].Mag())-m_prot+BE;
    pmiss_perp_2pto1p[0]=TMath::Sqrt(V3_prot_el_2pto1p[0].Px()*V3_prot_el_2pto1p[0].Px()+V3_prot_el_2pto1p[0].Py()*V3_prot_el_2pto1p[0].Py());
    P2pto1p[0]=N_p2to1[0]/N_2;
	    
    V3_prot_el_2pto1p[1]=V4el.Vect()+ V3prot[1];
    Ecal_2pto1p[1]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[1].Mag()*V3prot[1].Mag())-m_prot+BE;
    pmiss_perp_2pto1p[1]=TMath::Sqrt(V3_prot_el_2pto1p[1].Px()*V3_prot_el_2pto1p[1].Px()+V3_prot_el_2pto1p[1].Py()*V3_prot_el_2pto1p[1].Py());
    P2pto1p[1]=N_p2to1[1]/N_2;

  }
  
  *Nboth=N_2;

}

// __________________________________________________________________________________________________________________________________________________________________________________________________

void  treeProducer_simulation::prot_rot_func_3p(TVector3 V3q, TVector3  V3prot[3],TLorentzVector V4el,double Ecal_3pto2p[][2],double  pmiss_perp_3pto2p[][2],double  P3pto2p[][2],double N_p1[3],double Ecal_3pto1p[3],double  pmiss_perp_3pto1p[3], double *N_p3det) {

  double BE = -99.;

  if (Target == "3He") { BE = He3_bind_en-D2_bind_en; }
  if (Target == "4He") { BE = He4_bind_en-H3_bind_en; }
  if (Target == "12C") { BE = C12_bind_en - B_bind_en; }
  if (Target == "56Fe") { BE = Fe_bind_en-Mn_bind_en; }	

  double N_tot=100,m_prot=0.9382720813;
  const int N_3p=3, N_2p=2;
  double N_p2[N_3p]={0},N_p1det[3][N_3p]={0};
  TVector3 V3_3p_rot[N_3p],V3_2p_rot[N_3p],V3_prot_el_3pto2p[N_2p],V3_prot_el_3pto1p[N_3p];
  double rot_angle;
  double N_pthree=0;
  bool prot_stat[N_3p];
  int count =0;

  float cphil = 0, cphir = 0, pimi_phimin = 0, pimi_phimax = 0;  

  for (int g=0; g<N_tot; g++) {

    rot_angle=gRandom->Uniform(0,2*TMath::Pi());

    for(int i=0;i<N_3p;i++) {

      V3_3p_rot[i]= V3prot[i];
      V3_3p_rot[i].Rotate(rot_angle,V3q);

    }   

    for(int ind_p=0;ind_p<N_3p;ind_p++) prot_stat[ind_p]=ProtonFiducialCut(&V3_3p_rot[ind_p], &cphil, &cphir);
  
    if( prot_stat[0]  && !prot_stat[1]  && !prot_stat[2])  N_p1[0]=N_p1[0]+1;
    if(!prot_stat[0] &&   prot_stat[1]  && !prot_stat[2])  N_p1[1]=N_p1[1]+1;
    if(!prot_stat[0] &&  !prot_stat[1]  &&  prot_stat[2])  N_p1[2]=N_p1[2]+1;
    if( prot_stat[0]  &&  prot_stat[1]  &&  prot_stat[2])  N_pthree=N_pthree+1;

    if(prot_stat[0] && prot_stat[1] && !prot_stat[2])   N_p2[0]=N_p2[0]+1;
    if(prot_stat[0] && !prot_stat[1] && prot_stat[2])   N_p2[1]=N_p2[1]+1;
    if(!prot_stat[0] && prot_stat[1] && prot_stat[2])   N_p2[2]=N_p2[2]+1;

  } //Loop of 3p rotations ends
   
  //-----------------------------------------  3p to 1p  -----------------------------------------------------------------------

  for (int j = 0; j < N_3p; j++) { // Looping over 1p combinations out of 3protons
  
    V3_prot_el_3pto1p[j]=V4el.Vect()+ V3prot[j];
    Ecal_3pto1p[j]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[j].Mag()*V3prot[j].Mag())-m_prot+BE;
    pmiss_perp_3pto1p[j]=TMath::Sqrt(V3_prot_el_3pto1p[j].Px()*V3_prot_el_3pto1p[j].Px()+V3_prot_el_3pto1p[j].Py()*V3_prot_el_3pto1p[j].Py());

  }

  //-----------------------------------------  3p to 2p->1p  -----------------------------------------------------------------------

  for (int ind1=0;ind1<N_3p;ind1++){ // Looping over 2p combinations out of 3p

    for (int ind2=0;ind2<N_3p;ind2++) {

      if (ind1!=ind2 && ind1<ind2) {

	for (int g1 = 0; g1 < N_tot; g1++) {

	  rot_angle=gRandom->Uniform(0,2*TMath::Pi());
     
	  V3_2p_rot[ind1]=V3prot[ind1];
	  V3_2p_rot[ind2]=V3prot[ind2];
	  V3_2p_rot[ind1].Rotate(rot_angle,V3q);
	  V3_2p_rot[ind2].Rotate(rot_angle,V3q);       

	  if(ProtonFiducialCut(&V3_2p_rot[ind1], &cphil, &cphir)  && !ProtonFiducialCut(&V3_2p_rot[ind2], &cphil, &cphir)) N_p1det[count][0]=N_p1det[count][0]+1;
	  if(!ProtonFiducialCut(&V3_2p_rot[ind1], &cphil, &cphir) && ProtonFiducialCut(&V3_2p_rot[ind2], &cphil, &cphir))  N_p1det[count][1]=N_p1det[count][1]+1;
	  if(ProtonFiducialCut(&V3_2p_rot[ind1], &cphil, &cphir)  && ProtonFiducialCut(&V3_2p_rot[ind2], &cphil, &cphir))  N_p1det[count][2]=N_p1det[count][2]+1;
	     
	}

	if (N_p1det[count][2] != 0 && N_pthree != 0 ) { 


	  V3_prot_el_3pto2p[0]=V4el.Vect()+ V3prot[ind1];
	  Ecal_3pto2p[count][0]=V4el.E()+ TMath::Sqrt(m_prot*m_prot+V3prot[ind1].Mag()*V3prot[ind1].Mag())-m_prot+BE;
	  pmiss_perp_3pto2p[count][0]=TMath::Sqrt(V3_prot_el_3pto2p[0].Px()*V3_prot_el_3pto2p[0].Px()+V3_prot_el_3pto2p[0].Py()*V3_prot_el_3pto2p[0].Py());
	  P3pto2p[count][0]=N_p1det[count][0]/N_p1det[count][2]*(N_p2[count]/N_pthree);

	  V3_prot_el_3pto2p[1]=V4el.Vect()+ V3prot[ind2];
	  Ecal_3pto2p[count][1]=V4el.E()+TMath::Sqrt(m_prot*m_prot+V3prot[ind2].Mag()*V3prot[ind2].Mag())-m_prot+BE;
	  pmiss_perp_3pto2p[count][1]=TMath::Sqrt(V3_prot_el_3pto2p[1].Px()*V3_prot_el_3pto2p[1].Px()+V3_prot_el_3pto2p[1].Py()*V3_prot_el_3pto2p[1].Py());
	  P3pto2p[count][1]=N_p1det[count][1]/N_p1det[count][2]*(N_p2[count]/N_pthree);

	}
	
	count=count +1;

      }

    }

  }

  *N_p3det=N_pthree;

}

// __________________________________________________________________________________________________________________________________________________________________________________________________













// __________________________________________________________________________________________________________________________________________________________________________________________________














// __________________________________________________________________________________________________________________________________________________________________________________________________

#endif // #ifdef treeProducer_simulation_cxx
