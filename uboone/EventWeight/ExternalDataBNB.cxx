// Script to fill external data fit and covariance values for the BNB beam
//
//Values come from the different fit methods and meson type:
//      - Sanford-Wang: Pi+, Pi-, K0
//      - Extended Sanford-Wang: Pi+
//      - Feymann Scaling: K+
//      - Splines: Pi+
//External data used comes from:
//      - HARP
//      - E910 @6.4 GeV & @12.3 GeV
/////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <vector> 
using namespace std;

#include <cmath>
//#include <cstdlib>
//#include <iomanip>
#include <fstream>
#include <iostream>
#include "TMatrixD.h"
 
//#ifdef __MAKECINT__
//#pragma link C++ class std::vector < std::vector<float> >+;   
//#endif 

//__________________________________________________________________________
int ExternalDataBNB(){

  // Create a new ROOT binary machine.
  // This file may contain any kind of ROOT objects
  // This file is now becoming the current directory.

  std::cout<<" antes de abrir hfile"<<std::endl;

  TFile hfile("BNBExternalData_uBooNE.root","RECREATE","ROOT file with the result from the fit and their covariance, coming from the external data samples and the different fit methods");
  //  TTree *tree = new TTree("T","Tree containing the fit results and their covariances");
  
  /*
   * Make the TFile with the following structure:
   *
   * SW/  directory with all information from Sanford-Wang
   *    PiPlus/
   *       Fit
   *       Cov
   *    PiMinus/
   *       Fit
   *       Cov   
   *    K0/
   *       Fit
   *       Cov
   *       
   * FS/  directory with all information from Feymann scaling
   *    K+/
   *       Fit
   *       Cov
   */

  TDirectory* pSWDir = hfile.mkdir("SW", "Sanford-Wang parametrizations");
  
  TDirectory* pSWPiPlus = pSWDir->mkdir("PiPlus", "Sanford-Wang #pi^{+} parametrization");
  
  pSWPiPlus->cd();
  // prepare the data

  TArrayD SWPiPlusFitVal(9);

  SWPiPlusFitVal[0] = 220.7;
  SWPiPlusFitVal[1] = 1.080;
  SWPiPlusFitVal[2] = 1.000;
  SWPiPlusFitVal[3] = 1.978;
  SWPiPlusFitVal[4] = 1.32;
  SWPiPlusFitVal[5] = 5.572;
  SWPiPlusFitVal[6] = 0.0868;
  SWPiPlusFitVal[7] = 9.686;
  SWPiPlusFitVal[8] = 1.000;

  pSWPiPlus->WriteObject(&SWPiPlusFitVal, "SWPiPlusFitVal");
 
  std::cout<<" despues de llenar primer array"<<std::endl;
  
  TMatrixD SWPiPlusFitCov(9,9);

  // covariance values [row][column]
  
  SWPiPlusFitCov[0][0] = 1707.2;
  SWPiPlusFitCov[0][1] = 1.146;
  SWPiPlusFitCov[0][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[0][3] = -17.646;
  SWPiPlusFitCov[0][4] = -15.968;
  SWPiPlusFitCov[0][5] = -8.81;
  SWPiPlusFitCov[0][6] = -0.7347;
  SWPiPlusFitCov[0][7] = -60.816;
  SWPiPlusFitCov[0][8] = 0.00001;///no entries? cannot be <=0/// not in miniboone code!

  std::cout<<" despues de llenar primera linea primera matriz"<<std::endl;

  SWPiPlusFitCov[1][0] = 1.146;
  SWPiPlusFitCov[1][1] = 0.0396;
  SWPiPlusFitCov[1][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[1][3] = -0.1072;
  SWPiPlusFitCov[1][4] = -0.0993;
  SWPiPlusFitCov[1][5] = 0.0325;
  SWPiPlusFitCov[1][6] = 0.0007;
  SWPiPlusFitCov[1][7] = -0.0777;
  SWPiPlusFitCov[1][8] = 0.00001;///no entries? cannot be <=0///not in miniboone

  SWPiPlusFitCov[2][0] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[2][1] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[2][2] = 0.000001;///no entries? cannot be <=0
  SWPiPlusFitCov[2][3] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[2][4] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[2][5] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[2][6] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[2][7] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[2][8] = 0.;///no entries? cannot be <=0/// not in miniboone

  SWPiPlusFitCov[3][0] = -17.646;
  SWPiPlusFitCov[3][1] = -0.1072;
  SWPiPlusFitCov[3][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[3][3] = 0.5945;
  SWPiPlusFitCov[3][4] = 0.5049;
  SWPiPlusFitCov[3][5] = 0.0655;
  SWPiPlusFitCov[3][6] = 0.0025;
  SWPiPlusFitCov[3][7] = 0.198;
  SWPiPlusFitCov[3][8] = 0.;///no entries? cannot be <=0/// not in miniboone

  SWPiPlusFitCov[4][0] = -15.968;
  SWPiPlusFitCov[4][1] = -0.0993;
  SWPiPlusFitCov[4][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[4][3] = 0.5049;
  SWPiPlusFitCov[4][4] = 0.4411;
  SWPiPlusFitCov[4][5] = 0.0568;
  SWPiPlusFitCov[4][6] = 0.0025;
  SWPiPlusFitCov[4][7] = 0.2271;
  SWPiPlusFitCov[4][8] = 0.;///no entries? cannot be <=0

  SWPiPlusFitCov[5][0] = -8.81;
  SWPiPlusFitCov[5][1] = 0.0325;
  SWPiPlusFitCov[5][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[5][3] = 0.0655;
  SWPiPlusFitCov[5][4] = 0.0568;
  SWPiPlusFitCov[5][5] = 0.2066;
  SWPiPlusFitCov[5][6] = 0.0047;
  SWPiPlusFitCov[5][7] = 0.1031;
  SWPiPlusFitCov[5][8] = 0.;///no entries? cannot be <=0

  SWPiPlusFitCov[6][0] = -0.7347;
  SWPiPlusFitCov[6][1] = 0.0007;
  SWPiPlusFitCov[6][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[6][3] = 0.0025;
  SWPiPlusFitCov[6][4] = 0.0025;
  SWPiPlusFitCov[6][5] = 0.0047;
  SWPiPlusFitCov[6][6] = 0.0005;
  SWPiPlusFitCov[6][7] = 0.0641;
  SWPiPlusFitCov[6][8] = 0.;///no entries? cannot be <=0

  SWPiPlusFitCov[7][0] = -60.816;
  SWPiPlusFitCov[7][1] = -0.0777;
  SWPiPlusFitCov[7][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[7][3] = 0.198;
  SWPiPlusFitCov[7][4] = 0.2271;
  SWPiPlusFitCov[7][5] = 0.1031;
  SWPiPlusFitCov[7][6] = 0.0641;
  SWPiPlusFitCov[7][7] = 16.0189;
  SWPiPlusFitCov[7][8] = 0.;///no entries? cannot be <=0

  SWPiPlusFitCov[8][0] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][1] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][2] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][3] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][4] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][5] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][6] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][7] = 0.;///no entries? cannot be <=0
  SWPiPlusFitCov[8][8] = 0.;///no entries? cannot be <=0

  std::cout<<" despues de llenar primera matriz"<<std::endl;
  
  pSWPiPlus->WriteObject(&SWPiPlusFitCov, "SWPiPlusFitCov");


   TH1F *h1 = new TH1F("h1", "Pi+", 10,0.,221 );

  for (int i = 0; i<9; i++){
    std::cout<<" fit pi plus   "<<  SWPiPlusFitVal[i]<<std::endl;
    h1->Fill(SWPiPlusFitVal[i]);
      }
  h1->Write();


  TDirectory* pSWPiMinus = pSWDir->mkdir("PiMinus", "Sanford-Wang #pi^{-} parametrization");
  
  pSWPiMinus->cd();
  // prepare the data

  TArrayD SWPiMinusFitVal(9);

  SWPiMinusFitVal[0] = 213.7;
  SWPiMinusFitVal[1] = 0.9379;
  SWPiMinusFitVal[2] = 5.454;
  SWPiMinusFitVal[3] = 1.210;
  SWPiMinusFitVal[4] = 1.284;
  SWPiMinusFitVal[5] = 4.781;
  SWPiMinusFitVal[6] = 0.07338;
  SWPiMinusFitVal[7] = 8.329;
  SWPiMinusFitVal[8] = 1.000;

  pSWPiMinus->WriteObject(&SWPiMinusFitVal, "SWPiMinusFitVal");
 
  TMatrixD SWPiMinusFitCov(9,9);

  SWPiMinusFitCov[0][0] = 3688.9;
  SWPiMinusFitCov[0][1] = 7.61;
  SWPiMinusFitCov[0][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[0][3] = -15.666;
  SWPiMinusFitCov[0][4] = -17.48;
  SWPiMinusFitCov[0][5] = -11.329;
  SWPiMinusFitCov[0][6] = -0.9925;
  SWPiMinusFitCov[0][7] = -91.4;
  SWPiMinusFitCov[0][8] = 0.00001;///no entries? =1? cannot be <=0
  
  SWPiMinusFitCov[1][0] = 7.61;
  SWPiMinusFitCov[1][1] = 0.0388;
  SWPiMinusFitCov[1][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[1][3] = -0.0437;
  SWPiMinusFitCov[1][4] = -0.0509;
  SWPiMinusFitCov[1][5] = 0.0102;
  SWPiMinusFitCov[1][6] = -0.0009;
  SWPiMinusFitCov[1][7] = -0.1957;
  SWPiMinusFitCov[1][8] = 0.00001;///no entries? =1? cannot be <=0

  SWPiMinusFitCov[2][0] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][1] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][3] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][4] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][5] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][6] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][7] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[2][8] = 0.00001;///no entries? =1? cannot be <=0

  SWPiMinusFitCov[3][0] = -15.666;
  SWPiMinusFitCov[3][1] = -0.0437;
  SWPiMinusFitCov[3][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[3][3] = 0.0841;
  SWPiMinusFitCov[3][4] = 0.0895;
  SWPiMinusFitCov[3][5] = 0.0301;
  SWPiMinusFitCov[3][6] = 0.0029;
  SWPiMinusFitCov[3][7] = 0.2588;
  SWPiMinusFitCov[3][8] = 0.00001;///no entries? =1? cannot be <=0

  SWPiMinusFitCov[4][0] = -17.48;
  SWPiMinusFitCov[4][1] = -0.0509;
  SWPiMinusFitCov[4][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[4][3] = 0.0895;
  SWPiMinusFitCov[4][4] = 0.0986;
  SWPiMinusFitCov[4][5] = 0.0375;
  SWPiMinusFitCov[4][6] = 0.0033;
  SWPiMinusFitCov[4][7] = 0.3141;
  SWPiMinusFitCov[4][8] = 0.00001;///no entries? =1? cannot be <=0
  
  SWPiMinusFitCov[5][0] = -11.329;
  SWPiMinusFitCov[5][1] = 0.0102;
  SWPiMinusFitCov[5][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[5][3] = 0.0301;
  SWPiMinusFitCov[5][4] = 0.0375;
  SWPiMinusFitCov[5][5] = 0.1595;
  SWPiMinusFitCov[5][6] = 0.0051;
  SWPiMinusFitCov[5][7] = 0.1933;
  SWPiMinusFitCov[5][8] = 0.00001;///no entries? =1? cannot be <=0

  SWPiMinusFitCov[6][0] = -0.9925;
  SWPiMinusFitCov[6][1] = -0.0009;
  SWPiMinusFitCov[6][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[6][3] = 0.0029;
  SWPiMinusFitCov[6][4] = 0.0033;
  SWPiMinusFitCov[6][5] = 0.0051;
  SWPiMinusFitCov[6][6] = 0.0005;
  SWPiMinusFitCov[6][7] = 0.064;
  SWPiMinusFitCov[6][8] = 0.00001;///no entries? =1? cannot be <=0
  
  SWPiMinusFitCov[7][0] = -91.4;
  SWPiMinusFitCov[7][1] = -0.1957;
  SWPiMinusFitCov[7][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[7][3] = 0.2588;
  SWPiMinusFitCov[7][4] = 0.3141;
  SWPiMinusFitCov[7][5] = 0.1933;
  SWPiMinusFitCov[7][6] = 0.064;
  SWPiMinusFitCov[7][7] = 17.242;
  SWPiMinusFitCov[7][8] = 0.00001;///no entries? =1? cannot be <=0

  SWPiMinusFitCov[8][0] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][1] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][2] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][3] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][4] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][5] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][6] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][7] = 0.00001;///no entries? =1? cannot be <=0
  SWPiMinusFitCov[8][8] = 0.00001;///no entries? =1? cannot be <=0
  
  pSWPiMinus->WriteObject(&SWPiMinusFitCov, "SWPiMinusFitCov");

  TDirectory* pSWK0s = pSWDir->mkdir("K0s", "Sanford-Wang #K^{0} parametrization");
  
  pSWK0s->cd();
  // prepare the data

  TArrayD SWK0sFitVal(9);

  SWK0sFitVal[0] = 15.130;
  SWK0sFitVal[1] = 1.975;
  SWK0sFitVal[2] = 4.084;
  SWK0sFitVal[3] = 0.928;
  SWK0sFitVal[4] = 0.731;
  SWK0sFitVal[5] = 4.362;
  SWK0sFitVal[6] = 0.048;
  SWK0sFitVal[7] = 13.3;
  SWK0sFitVal[8] = 1.278;

  pSWK0s->WriteObject(&SWK0sFitVal, "SWK0sFitVal");

  TMatrixD SWK0sFitCov(9,9);

  SWK0sFitCov[0][0] = 32.3;
  SWK0sFitCov[0][1] = -0.09687;
  SWK0sFitCov[0][2] = 0.8215;
  SWK0sFitCov[0][3] = -0.1018;
  SWK0sFitCov[0][4] = -0.2124;
  SWK0sFitCov[0][5] = -0.8902;
  SWK0sFitCov[0][6] = -0.1333;
  SWK0sFitCov[0][7] = 16.55;
  SWK0sFitCov[0][8] = -1.789;

  SWK0sFitCov[1][0] = -0.09687;
  SWK0sFitCov[1][1] = 0.09574;
  SWK0sFitCov[1][2] = 0.03248;
  SWK0sFitCov[1][3] = 0.00131;
  SWK0sFitCov[1][4] = -0.01303;
  SWK0sFitCov[1][5] = 0.08836;
  SWK0sFitCov[1][6] = -0.00031;
  SWK0sFitCov[1][7] = -1.536;
  SWK0sFitCov[1][8] = -0.2156;

  SWK0sFitCov[2][0] = 0.8215;
  SWK0sFitCov[2][1] = 0.03248;
  SWK0sFitCov[2][2] = 0.5283;
  SWK0sFitCov[2][3] = -0.01922;
  SWK0sFitCov[2][4] = 0.02267;
  SWK0sFitCov[2][5] = -0.0033;
  SWK0sFitCov[2][6] = -0.00236;
  SWK0sFitCov[2][7] = 0.0391;
  SWK0sFitCov[2][8] = -0.08017;

  SWK0sFitCov[3][0] = -0.1018;
  SWK0sFitCov[3][1] = 0.00131;
  SWK0sFitCov[3][2] = -0.01922;
  SWK0sFitCov[3][3] = 0.08442;
  SWK0sFitCov[3][4] = 0.00405;
  SWK0sFitCov[3][5] = 0.00071;
  SWK0sFitCov[3][6] = -0.00037;
  SWK0sFitCov[3][7] = -0.01443;
  SWK0sFitCov[3][8] = -0.07301;

  SWK0sFitCov[4][0] = -0.2124;
  SWK0sFitCov[4][1] = -0.01303;
  SWK0sFitCov[4][2] = 0.02267;
  SWK0sFitCov[4][3] = 0.00405;
  SWK0sFitCov[4][4] = 0.00982;
  SWK0sFitCov[4][5] = 0.00287;
  SWK0sFitCov[4][6] = 0.00028;
  SWK0sFitCov[4][7] = -0.05777;
  SWK0sFitCov[4][8] = 0.02966;

  SWK0sFitCov[5][0] = -0.8902;
  SWK0sFitCov[5][1] = 0.08836;
  SWK0sFitCov[5][2] = -0.0033;
  SWK0sFitCov[5][3] = 0.00071;
  SWK0sFitCov[5][4] = 0.00287;
  SWK0sFitCov[5][5] = 0.3599;
  SWK0sFitCov[5][6] = 0.00385;
  SWK0sFitCov[5][7] = -4.751;
  SWK0sFitCov[5][8] = -0.1577;

  SWK0sFitCov[6][0] = -0.1333;
  SWK0sFitCov[6][1] = -0.00031;
  SWK0sFitCov[6][2] = -0.00236;
  SWK0sFitCov[6][3] = -0.00037;
  SWK0sFitCov[6][4] = 0.00028;
  SWK0sFitCov[6][5] = 0.00385;
  SWK0sFitCov[6][6] = 0.00105;
  SWK0sFitCov[6][7] = 0.05806;
  SWK0sFitCov[6][8] = 0.00686;

  SWK0sFitCov[7][0] = 16.55;
  SWK0sFitCov[7][1] = -1.536;
  SWK0sFitCov[7][2] = 0.0391;
  SWK0sFitCov[7][3] = -0.01443;
  SWK0sFitCov[7][4] = -0.05777;
  SWK0sFitCov[7][5] = -4.751;
  SWK0sFitCov[7][6] = 0.05806;
  SWK0sFitCov[7][7] = 130.2;
  SWK0sFitCov[7][8] = 1.222;

  SWK0sFitCov[8][0] = -1.798;
  SWK0sFitCov[8][1] = -0.2156;
  SWK0sFitCov[8][2] = -0.08017;
  SWK0sFitCov[8][3] = -0.07301;
  SWK0sFitCov[8][4] = 0.02966;
  SWK0sFitCov[8][5] = -0.1577;
  SWK0sFitCov[8][6] = 0.00686;
  SWK0sFitCov[8][7] = 1.222;
  SWK0sFitCov[8][8] = 2.948;
 
  //  SWK0sFitCov->Write();

  pSWK0s->WriteObject(&SWK0sFitCov, "SWK0sFitCov");

  TDirectory* pFSDir = hfile.mkdir("FS", "Feynman scaling parametrizations");
  
  TDirectory* pFSKPlus = pFSDir->mkdir("KPlus", "Feynmann K^{+} parametrization");
  
  pFSKPlus->cd();
  // prepare the data

  TArrayD FSKPlusFitVal(7);

  FSKPlusFitVal[0] = 11.70;
  FSKPlusFitVal[1] = 0.88;
  FSKPlusFitVal[2] = 4.77;
  FSKPlusFitVal[3] = 1.51;
  FSKPlusFitVal[4] = 2.21;
  FSKPlusFitVal[5] = 2.17;
  FSKPlusFitVal[6] = 1.51;

  //  FSKPlusFitVal->Write();

  pFSKPlus->WriteObject(&FSKPlusFitVal, "FSKPlusFitVal");

  TMatrixD FSKPlusFitCov(7,7);
   
  FSKPlusFitCov[0][0] = 1.094;
  FSKPlusFitCov[0][1] = 0.0502;
  FSKPlusFitCov[0][2] = 0.00299;
  FSKPlusFitCov[0][3] = -0.0332;
  FSKPlusFitCov[0][4] = -0.0375;
  FSKPlusFitCov[0][5] = 0.125;
  FSKPlusFitCov[0][6] = 0.0743;

  FSKPlusFitCov[1][0] = 0.0502;
  FSKPlusFitCov[1][1] = 0.0161;
  FSKPlusFitCov[1][2] = 0.00139;
  FSKPlusFitCov[1][3] = -0.00144;
  FSKPlusFitCov[1][4] = -0.0126;
  FSKPlusFitCov[1][5] = 0.0322;
  FSKPlusFitCov[1][6] = 0.022;

  FSKPlusFitCov[2][0] = 0.00299;
  FSKPlusFitCov[2][1] = 0.00139;
  FSKPlusFitCov[2][2] = 0.00747;
  FSKPlusFitCov[2][3] = 0.00206;
  FSKPlusFitCov[2][4] = 0.00193;
  FSKPlusFitCov[2][5] = 0.0135;
  FSKPlusFitCov[2][6] = -0.00334;

  FSKPlusFitCov[3][0] = -0.0332;
  FSKPlusFitCov[3][1] = -0.00144;
  FSKPlusFitCov[3][2] = 0.00206;
  FSKPlusFitCov[3][3] = 0.00346;
  FSKPlusFitCov[3][4] = 0.00203;
  FSKPlusFitCov[3][5] = -0.00411;
  FSKPlusFitCov[3][6] = -0.00628;

  FSKPlusFitCov[4][0] = -0.0375;
  FSKPlusFitCov[4][1] = -0.0126;
  FSKPlusFitCov[4][2] = 0.00193;
  FSKPlusFitCov[4][3] = 0.00203;
  FSKPlusFitCov[4][4] = 0.0146;
  FSKPlusFitCov[4][5] = -0.0154;
  FSKPlusFitCov[4][6] = -0.0244;

  FSKPlusFitCov[5][0] = 0.125;
  FSKPlusFitCov[5][1] = 0.0322;
  FSKPlusFitCov[5][2] = 0.0135;
  FSKPlusFitCov[5][3] = -0.00411;
  FSKPlusFitCov[5][4] = -0.0154;
  FSKPlusFitCov[5][5] = 0.182;
  FSKPlusFitCov[5][6] = 0.126;

  FSKPlusFitCov[6][0] = 0.0743;
  FSKPlusFitCov[6][1] = 0.022;
  FSKPlusFitCov[6][2] = -0.00334;
  FSKPlusFitCov[6][3] = -0.00628;
  FSKPlusFitCov[6][4] = -0.0244;
  FSKPlusFitCov[6][5] = 0.126;
  FSKPlusFitCov[6][6] = 0.159;

  //  FSKPlusFitCov->Write();

  pFSKPlus->WriteObject(&FSKPlusFitCov, "FSKPlusFitCov");

  // covariance values [row][column]

//      tree->Fill();

  std::cout<<" Fit results Pi+"<<std::endl;


  // std::cout<<" c1  : "<<SWPiPlusFitVal[0]<<" c2 : "<<  SWPiPlusFitVal[1]<<" c3 : "<<  SWPiPlusFitVal[2]<<" c4 : "<<  SWPiPlusFitVal[3]<<" c5 : "<<  SWPiPlusFitVal[4]<<" c6 : "<<  SWPiPlusFitVal[5]<<" c7 : "<<  SWPiPlusFitVal[6]<<" c8 : "<<  SWPiPlusFitVal[7]<<" c9 : "<<  SWPiPlusFitVal[8]<<std::endl;

  //std::cout<<" Covariance Pi+"<<std::endl;

  //std::cout<<" SWPiPlusFitCov[0][0]  : "<<SWPiPlusFitCov[0][0]<<" SWPiPlusFitCov[0][1] : "<<  SWPiPlusFitCov[0][1]<<" SWPiPlusFitCov[0][2] : "<<  SWPiPlusFitCov[0][2]<<" SWPiPlusFitCov[0][3] : "<<  SWPiPlusFitCov[0][3]<<" SWPiPlusFitCov[0][4] : "<<  SWPiPlusFitCov[0][4]<<" SWPiPlusFitCov[0][5] : "<<  SWPiPlusFitCov[0][5]<<" SWPiPlusFitCov[0][6] : "<<  SWPiPlusFitCov[0][6]<<" SWPiPlusFitCov[0][7] : "<<  SWPiPlusFitCov[0][7]<<" SWPiPlusFitCov[0][8] : "<<  SWPiPlusFitCov[0][8]<<std::endl;

  //std::cout<<" SWPiPlusFitCov[1][0]  : "<<SWPiPlusFitCov[1][0]<<" SWPiPlusFitCov[1][1] : "<<  SWPiPlusFitCov[1][1]<<" SWPiPlusFitCov[1][2] : "<<  SWPiPlusFitCov[1][2]<<" SWPiPlusFitCov[1][3] : "<<  SWPiPlusFitCov[1][3]<<" SWPiPlusFitCov[1][4] : "<<  SWPiPlusFitCov[1][4]<<" SWPiPlusFitCov[1][5] : "<<  SWPiPlusFitCov[1][5]<<" SWPiPlusFitCov[1][6] : "<<  SWPiPlusFitCov[1][6]<<" SWPiPlusFitCov[1][7] : "<<  SWPiPlusFitCov[1][7]<<" SWPiPlusFitCov[0][8] : "<<  SWPiPlusFitCov[1][8]<<std::endl;

  ///////
  
  //std::cout<<" Fit results Pi-"<<std::endl;
  
  //std::cout<<" c1  : "<<SWPiMinusFitVal[0]<<" c2 : "<<  SWPiMinusFitVal[1]<<" c3 : "<<  SWPiMinusFitVal[2]<<" c4 : "<<  SWPiMinusFitVal[3]<<" c5 : "<<  SWPiMinusFitVal[4]<<" c6 : "<<  SWPiMinusFitVal[5]<<" c7 : "<<  SWPiMinusFitVal[6]<<" c8 : "<<  SWPiMinusFitVal[7]<<" c9 : "<<  SWPiMinusFitVal[8]<<std::endl;

  //std::cout<<" Covariance Pi-"<<std::endl;

  //std::cout<<" SWPiMinusFitCov[0][0]  : "<<SWPiMinusFitCov[0][0]<<" SWPiMinusFitCov[0][1] : "<<  SWPiMinusFitCov[0][1]<<" SWPiMinusFitCov[0][2] : "<<  SWPiMinusFitCov[0][2]<<" SWPiMinusFitCov[0][3] : "<<  SWPiMinusFitCov[0][3]<<" SWPiMinusFitCov[0][4] : "<<  SWPiMinusFitCov[0][4]<<" SWPiMinusFitCov[0][5] : "<<  SWPiMinusFitCov[0][5]<<" SWPiMinusFitCov[0][6] : "<<  SWPiMinusFitCov[0][6]<<" SWPiMinusFitCov[0][7] : "<<  SWPiMinusFitCov[0][7]<<" SWPiMinusFitCov[0][8] : "<<  SWPiMinusFitCov[0][8]<<std::endl;

  //std::cout<<" SWPiMinusFitCov[1][0]  : "<<SWPiMinusFitCov[1][0]<<" SWPiMinusFitCov[1][1] : "<<  SWPiMinusFitCov[1][1]<<" SWPiMinusFitCov[1][2] : "<<  SWPiMinusFitCov[1][2]<<" SWPiMinusFitCov[1][3] : "<<  SWPiMinusFitCov[1][3]<<" SWPiMinusFitCov[1][4] : "<<  SWPiMinusFitCov[1][4]<<" SWPiMinusFitCov[1][5] : "<<  SWPiMinusFitCov[1][5]<<" SWPiMinusFitCov[1][6] : "<<  SWPiMinusFitCov[1][6]<<" SWPiMinusFitCov[1][7] : "<<  SWPiMinusFitCov[1][7]<<" SWPiMinusFitCov[0][8] : "<<  SWPiMinusFitCov[1][8]<<std::endl;

  //std::cout<<" Fit results K+"<<std::endl;

  //std::cout<<" c1  : "<<FSKPlusFitVal[0]<<" c2 : "<<  FSKPlusFitVal[1]<<" c3 : "<<  FSKPlusFitVal[2]<<" c4 : "<<  FSKPlusFitVal[3]<<" c5 : "<<  FSKPlusFitVal[4]<<" c6 : "<<  FSKPlusFitVal[5]<<" c7 : "<<  FSKPlusFitVal[6]<<std::endl;

  //std::cout<<" Covariance K+"<<std::endl;

  //std::cout<<" FSKPlusFitCov[0][0]  : "<<FSKPlusFitCov[0][0]<<" FSKPlusFitCov[0][1] : "<<  FSKPlusFitCov[0][1]<<" FSKPlusFitCov[0][2] : "<<  FSKPlusFitCov[0][2]<<" FSKPlusFitCov[0][3] : "<<  FSKPlusFitCov[0][3]<<" FSKPlusFitCov[0][4] : "<<  FSKPlusFitCov[0][4]<<" FSKPlusFitCov[0][5] : "<<  FSKPlusFitCov[0][5]<<" FSKPlusFitCov[0][6] : "<<  FSKPlusFitCov[0][6]<<std::endl;

  //std::cout<<" FSKPlusFitCov[1][0]  : "<<FSKPlusFitCov[1][0]<<" FSKPlusFitCov[1][1] : "<<  FSKPlusFitCov[1][1]<<" FSKPlusFitCov[1][2] : "<<  FSKPlusFitCov[1][2]<<" FSKPlusFitCov[1][3] : "<<  FSKPlusFitCov[1][3]<<" FSKPlusFitCov[1][4] : "<<  FSKPlusFitCov[1][4]<<" FSKPlusFitCov[1][5] : "<<  FSKPlusFitCov[1][5]<<" FSKPlusFitCov[1][6] : "<<  FSKPlusFitCov[1][6]<<std::endl;

  //std::cout<<" Fit results K0_s"<<std::endl;

  //std::cout<<" c1  : "<<SWK0sFitVal[0]<<" c2 : "<<  SWK0sFitVal[1]<<" c3 : "<<  SWK0sFitVal[2]<<" c4 : "<<  SWK0sFitVal[3]<<" c5 : "<<  SWK0sFitVal[4]<<" c6 : "<<  SWK0sFitVal[5]<<" c7 : "<<  SWK0sFitVal[6]<<" c8 : "<<  SWK0sFitVal[7]<<" c9 : "<<  SWK0sFitVal[8]<<std::endl;

  //std::cout<<" Covariance K+"<<std::endl;

  //std::cout<<" SWK0sFitCov[0][0]  : "<<SWK0sFitCov[0][0]<<" SWK0sFitCov[0][1] : "<<  SWK0sFitCov[0][1]<<" SWK0sFitCov[0][2] : "<<  SWK0sFitCov[0][2]<<" SWK0sFitCov[0][3] : "<<  SWK0sFitCov[0][3]<<" SWK0sFitCov[0][4] : "<<  SWK0sFitCov[0][4]<<" SWK0sFitCov[0][5] : "<<  SWK0sFitCov[0][5]<<" SWK0sFitCov[0][6] : "<<  SWK0sFitCov[0][6]<<" SWK0sFitCov[0][7] : "<<  SWK0sFitCov[0][7]<<" SWK0sFitCov[0][8] : "<<  SWK0sFitCov[0][8]<<std::endl;
  
  //  hfile.Write();



  std::cout<<" SWPiPlusFitCov[0][0]  : "<<SWPiPlusFitCov[0][0]<<" SWPiPlusFitCov[0][1] : "<<  SWPiPlusFitCov[0][1]<<" SWPiPlusFitCov[1][0] : "<<  SWPiPlusFitCov[1][0]<<" SWPiPlusFitCov[1][1] : "<<  SWPiPlusFitCov[1][1]<<std::endl;




  //  delete hfile;
  return 0;
}
