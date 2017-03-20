////////////////////////////////////////////////////////////////////////
// \file SpaceChargeMicroBooNE.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for MicroBooNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// LArSoft includes
#include "uboone/SpaceCharge/SpaceChargeMicroBooNE.h"

// Framework includes
#include "cetlib/exception.h"

double spacecharge::SpaceChargeMicroBooNE::fPos[3];
double spacecharge::SpaceChargeMicroBooNE::parA[6][7];
double spacecharge::SpaceChargeMicroBooNE::parB[6];

//-----------------------------------------------
spacecharge::SpaceChargeMicroBooNE::SpaceChargeMicroBooNE(
  fhicl::ParameterSet const& pset
)
{
  Configure(pset);
}

void spacecharge::SpaceChargeMicroBooNE::ParseRepresentationType(std::string s)
{
  if(s=="Parametric")
    fRepresentationType = SpaceChargeRepresentation_t::kParametric;
  else
    fRepresentationType = SpaceChargeRepresentation_t::kUnknown;   
}

//------------------------------------------------
bool spacecharge::SpaceChargeMicroBooNE::Configure(fhicl::ParameterSet const& pset)
{  
  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  fEnableCorrSCE = pset.get<bool>("EnableCorrSCE");

  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
  {
    ParseRepresentationType(pset.get<std::string>("RepresentationType"));

    fInputFilename = pset.get<std::string>("InputFilename");

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename,fname);

    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()) throw cet::exception("SpaceChargeMicroBooNE") << "Could not find the space charge effect file '" << fname << "'!\n";

    if(fRepresentationType==SpaceChargeRepresentation_t::kParametric)
    {      
      for(int i = 0; i < 5; i++)
      {
        g1_x[i] = (TGraph*)infile->Get(Form("deltaX/g1_%d",i));
        g2_x[i] = (TGraph*)infile->Get(Form("deltaX/g2_%d",i));
        g3_x[i] = (TGraph*)infile->Get(Form("deltaX/g3_%d",i));   
        g4_x[i] = (TGraph*)infile->Get(Form("deltaX/g4_%d",i));
        g5_x[i] = (TGraph*)infile->Get(Form("deltaX/g5_%d",i));

        g1_y[i] = (TGraph*)infile->Get(Form("deltaY/g1_%d",i));
        g2_y[i] = (TGraph*)infile->Get(Form("deltaY/g2_%d",i));
        g3_y[i] = (TGraph*)infile->Get(Form("deltaY/g3_%d",i));   
        g4_y[i] = (TGraph*)infile->Get(Form("deltaY/g4_%d",i));
        g5_y[i] = (TGraph*)infile->Get(Form("deltaY/g5_%d",i));
        g6_y[i] = (TGraph*)infile->Get(Form("deltaY/g6_%d",i));

        g1_z[i] = (TGraph*)infile->Get(Form("deltaZ/g1_%d",i));
        g2_z[i] = (TGraph*)infile->Get(Form("deltaZ/g2_%d",i));
        g3_z[i] = (TGraph*)infile->Get(Form("deltaZ/g3_%d",i));   
        g4_z[i] = (TGraph*)infile->Get(Form("deltaZ/g4_%d",i));

	g1_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g1_%d",i));
	g2_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g2_%d",i));
	g3_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g3_%d",i));
	g4_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g4_%d",i));
	g5_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g5_%d",i));

	g1_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g1_%d",i));
	g2_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g2_%d",i));
	g3_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g3_%d",i));
	g4_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g4_%d",i));
	g5_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g5_%d",i));
	g6_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g6_%d",i));

	g1_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g1_%d",i));
	g2_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g2_%d",i));
	g3_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g3_%d",i));
	g4_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g4_%d",i));
      }

      g1_x[5] = (TGraph*)infile->Get("deltaX/g1_5");
      g2_x[5] = (TGraph*)infile->Get("deltaX/g2_5");
      g3_x[5] = (TGraph*)infile->Get("deltaX/g3_5");   
      g4_x[5] = (TGraph*)infile->Get("deltaX/g4_5");
      g5_x[5] = (TGraph*)infile->Get("deltaX/g5_5");

      g1_y[5] = (TGraph*)infile->Get("deltaY/g1_5");
      g2_y[5] = (TGraph*)infile->Get("deltaY/g2_5");
      g3_y[5] = (TGraph*)infile->Get("deltaY/g3_5");   
      g4_y[5] = (TGraph*)infile->Get("deltaY/g4_5");
      g5_y[5] = (TGraph*)infile->Get("deltaY/g5_5");
      g6_y[5] = (TGraph*)infile->Get("deltaY/g6_5");
      
      g1_x[6] = (TGraph*)infile->Get("deltaX/g1_6");
      g2_x[6] = (TGraph*)infile->Get("deltaX/g2_6");
      g3_x[6] = (TGraph*)infile->Get("deltaX/g3_6");
      g4_x[6] = (TGraph*)infile->Get("deltaX/g4_6");
      g5_x[6] = (TGraph*)infile->Get("deltaX/g5_6");

      g1_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g1_5");
      g2_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g2_5");
      g3_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g3_5");
      g4_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g4_5");
      g5_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g5_5");

      g1_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g1_5");
      g2_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g2_5");
      g3_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g3_5");
      g4_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g4_5");
      g5_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g5_5");
      g6_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g6_5");

      g1_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g1_6");
      g2_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g2_6");
      g3_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g3_6");
      g4_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g4_6");
      g5_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g5_6");
    }

    infile->Close();
  }

  if(fEnableCorrSCE == true)
  {
    // Grab other parameters from pset  
  }

  return true;
}

//------------------------------------------------
bool spacecharge::SpaceChargeMicroBooNE::Update(uint64_t ts) 
{
  if (ts == 0) return false;

  return true;
}

//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// spatial distortions
bool spacecharge::SpaceChargeMicroBooNE::EnableSimSpatialSCE() const
{
  return fEnableSimSpatialSCE;
}

//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// E field distortions
bool spacecharge::SpaceChargeMicroBooNE::EnableSimEfieldSCE() const
{
  return fEnableSimEfieldSCE;
}

//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeMicroBooNE::EnableCorrSCE() const
{
  return fEnableCorrSCE;
}

void spacecharge::SpaceChargeMicroBooNE::SetPosition(double xVal, double yVal, double zVal) const
{
  fPos[0] = TransformX(xVal); fPos[1] = TransformY(yVal); fPos[2] = TransformZ(zVal);
}

//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in ionization electron drift
void spacecharge::SpaceChargeMicroBooNE::GetPosOffsets(double xVal, double yVal, double zVal,
						     std::vector<double> & posOffsets) const
{
  posOffsets.resize(3,0.0);
  if(!IsInsideBoundaries(xVal,yVal,zVal))
    return;

  SetPosition(xVal,yVal,zVal);
  
  if(fRepresentationType==SpaceChargeRepresentation_t::kParametric){
    posOffsets[0] = GetXPosOffsetParametric();
    posOffsets[1] = GetYPosOffsetParametric();
    posOffsets[2] = GetZPosOffsetParametric();
  }

}

//----------------------------------------------------------------------------
/// Provides X position offset using a parametric representation
double spacecharge::SpaceChargeMicroBooNE::GetXPosOffsetParametric() const
{      

  for(int j = 0; j < 7; ++j)
    {
      parA[0][j] = g1_x[j]->Eval(fPos[2]);
      parA[1][j] = g2_x[j]->Eval(fPos[2]);
      parA[2][j] = g3_x[j]->Eval(fPos[2]);
      parA[3][j] = g4_x[j]->Eval(fPos[2]);
      parA[4][j] = g5_x[j]->Eval(fPos[2]);
    }
  
  f1_x->SetParameters(parA[0]);
  f2_x->SetParameters(parA[1]);
  f3_x->SetParameters(parA[2]);
  f4_x->SetParameters(parA[3]);
  f5_x->SetParameters(parA[4]);
  
  parB[0] = f1_x->Eval(fPos[1]);
  parB[1] = f2_x->Eval(fPos[1]);
  parB[2] = f3_x->Eval(fPos[1]);
  parB[3] = f4_x->Eval(fPos[1]);
  parB[4] = f5_x->Eval(fPos[1]);
  
  fFinal_x->SetParameters(parB);
  return 100.0*fFinal_x->Eval(fPos[0]);
  
}


//----------------------------------------------------------------------------
/// Provides y position offset using a parametric representation
double spacecharge::SpaceChargeMicroBooNE::GetYPosOffsetParametric() const
{      
  
  for(int j = 0; j < 6; j++)
    {
      parA[0][j] = g1_y[j]->Eval(fPos[2]);
      parA[1][j] = g2_y[j]->Eval(fPos[2]);
      parA[2][j] = g3_y[j]->Eval(fPos[2]);
      parA[3][j] = g4_y[j]->Eval(fPos[2]);
      parA[4][j] = g5_y[j]->Eval(fPos[2]);
      parA[5][j] = g6_y[j]->Eval(fPos[2]);
    }
  
  f1_y->SetParameters(parA[0]);
  f2_y->SetParameters(parA[1]);
  f3_y->SetParameters(parA[2]);
  f4_y->SetParameters(parA[3]);
  f5_y->SetParameters(parA[4]);
  f6_y->SetParameters(parA[5]);

  parB[0] = f1_y->Eval(fPos[0]);
  parB[1] = f2_y->Eval(fPos[0]);
  parB[2] = f3_y->Eval(fPos[0]);
  parB[3] = f4_y->Eval(fPos[0]);
  parB[4] = f5_y->Eval(fPos[0]);
  parB[5] = f6_y->Eval(fPos[0]);
  
  fFinal_y->SetParameters(parB);
  return 100.0*fFinal_y->Eval(fPos[1]);
}

//----------------------------------------------------------------------------
/// Provides z position offset using a parametric representation
double spacecharge::SpaceChargeMicroBooNE::GetZPosOffsetParametric() const
{      
  static double parA[6][7];
  static double parB[6];
  
  for(int j = 0; j < 5; j++)
    {
      parA[0][j] = g1_z[j]->Eval(fPos[2]);
      parA[1][j] = g2_z[j]->Eval(fPos[2]);
      parA[2][j] = g3_z[j]->Eval(fPos[2]);
      parA[3][j] = g4_z[j]->Eval(fPos[2]);
    }
  
  f1_z->SetParameters(parA[0]);
  f2_z->SetParameters(parA[1]);
  f3_z->SetParameters(parA[2]);
  f4_z->SetParameters(parA[3]);
  
  parB[0] = f1_z->Eval(fPos[1]);
  parB[1] = f2_z->Eval(fPos[1]);
  parB[2] = f3_z->Eval(fPos[1]);
  parB[3] = f4_z->Eval(fPos[1]);
  
  fFinal_z->SetParameters(parB);
  return 100.0*fFinal_z->Eval(fPos[0]);
}

//----------------------------------------------------------------------------
/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.)
void spacecharge::SpaceChargeMicroBooNE::GetEfieldOffsets(double xVal, double yVal, double zVal,
							std::vector<double> & efieldOffsets) const
{
  efieldOffsets.resize(3,0.0);

  SetPosition(xVal,yVal,zVal);
  
  if(fRepresentationType==SpaceChargeRepresentation_t::kParametric){
    efieldOffsets[0] = -1.*GetXEfieldOffsetParametric();
    efieldOffsets[1] = -1.*GetYEfieldOffsetParametric();
    efieldOffsets[2] = -1.*GetZEfieldOffsetParametric();
  }
  
}

//----------------------------------------------------------------------------
/// Provides X position offset using a parametric representation
double spacecharge::SpaceChargeMicroBooNE::GetXEfieldOffsetParametric() const
{      

  for(int j = 0; j < 7; ++j)
    {
      parA[0][j] = g1_Ex[j]->Eval(fPos[2]);
      parA[1][j] = g2_Ex[j]->Eval(fPos[2]);
      parA[2][j] = g3_Ex[j]->Eval(fPos[2]);
      parA[3][j] = g4_Ex[j]->Eval(fPos[2]);
      parA[4][j] = g5_Ex[j]->Eval(fPos[2]);
    }
  
  f1_Ex->SetParameters(parA[0]);
  f2_Ex->SetParameters(parA[1]);
  f3_Ex->SetParameters(parA[2]);
  f4_Ex->SetParameters(parA[3]);
  f5_Ex->SetParameters(parA[4]);
  
  parB[0] = f1_Ex->Eval(fPos[1]);
  parB[1] = f2_Ex->Eval(fPos[1]);
  parB[2] = f3_Ex->Eval(fPos[1]);
  parB[3] = f4_Ex->Eval(fPos[1]);
  parB[4] = f5_Ex->Eval(fPos[1]);
  
  fFinal_Ex->SetParameters(parB);
  return fFinal_Ex->Eval(fPos[0]);
  
}


//----------------------------------------------------------------------------
/// Provides y position offset using a parametric representation
double spacecharge::SpaceChargeMicroBooNE::GetYEfieldOffsetParametric() const
{

  for(int j = 0; j < 6; j++)
    {
      parA[0][j] = g1_Ey[j]->Eval(fPos[2]);
      parA[1][j] = g2_Ey[j]->Eval(fPos[2]);
      parA[2][j] = g3_Ey[j]->Eval(fPos[2]);
      parA[3][j] = g4_Ey[j]->Eval(fPos[2]);
      parA[4][j] = g5_Ey[j]->Eval(fPos[2]);
      parA[5][j] = g6_Ey[j]->Eval(fPos[2]);
    }
  
  f1_Ey->SetParameters(parA[0]);
  f2_Ey->SetParameters(parA[1]);
  f3_Ey->SetParameters(parA[2]);
  f4_Ey->SetParameters(parA[3]);
  f5_Ey->SetParameters(parA[4]);
  f6_Ey->SetParameters(parA[5]);

  parB[0] = f1_Ey->Eval(fPos[0]);
  parB[1] = f2_Ey->Eval(fPos[0]);
  parB[2] = f3_Ey->Eval(fPos[0]);
  parB[3] = f4_Ey->Eval(fPos[0]);
  parB[4] = f5_Ey->Eval(fPos[0]);
  parB[5] = f6_Ey->Eval(fPos[0]);
  
  fFinal_Ey->SetParameters(parB);
  return fFinal_Ey->Eval(fPos[1]);
}

//----------------------------------------------------------------------------
/// Provides z position offset using a parametric representation
double spacecharge::SpaceChargeMicroBooNE::GetZEfieldOffsetParametric() const
{
  
  for(int j = 0; j < 5; j++)
    {
      parA[0][j] = g1_Ez[j]->Eval(fPos[2]);
      parA[1][j] = g2_Ez[j]->Eval(fPos[2]);
      parA[2][j] = g3_Ez[j]->Eval(fPos[2]);
      parA[3][j] = g4_Ez[j]->Eval(fPos[2]);
    }
  
  f1_Ez->SetParameters(parA[0]);
  f2_Ez->SetParameters(parA[1]);
  f3_Ez->SetParameters(parA[2]);
  f4_Ez->SetParameters(parA[3]);
  
  parB[0] = f1_Ez->Eval(fPos[1]);
  parB[1] = f2_Ez->Eval(fPos[1]);
  parB[2] = f3_Ez->Eval(fPos[1]);
  parB[3] = f4_Ez->Eval(fPos[1]);
  
  fFinal_Ez->SetParameters(parB);
  return fFinal_Ez->Eval(fPos[0]);
}

//----------------------------------------------------------------------------
/// Transform X to SCE X coordinate:  [2.56,0.0] --> [0.0,2.50]
double spacecharge::SpaceChargeMicroBooNE::TransformX(double xVal) const
{
  return 2.50 - (2.50/2.56)*(xVal/100.0) - 1.25;
}

//----------------------------------------------------------------------------
/// Transform Y to SCE Y coordinate:  [-1.165,1.165] --> [0.0,2.50]
double spacecharge::SpaceChargeMicroBooNE::TransformY(double yVal) const
{
  return (2.50/2.33)*((yVal/100.0)+1.165) - 1.25;
}

//----------------------------------------------------------------------------
/// Transform Z to SCE Z coordinate:  [0.0,10.37] --> [0.0,10.0]
double spacecharge::SpaceChargeMicroBooNE::TransformZ(double zVal) const
{
  return (10.0/10.37)*(zVal/100.0);
}

//----------------------------------------------------------------------------
/// Check to see if point is inside boundaries of map (allow to go slightly out of range)
bool spacecharge::SpaceChargeMicroBooNE::IsInsideBoundaries(double xVal, double yVal, double zVal) const
{
  if((xVal < 0.0) || (xVal > 260.0) || (yVal < -120.0) || (yVal > 120.0) || (zVal < 0.0) || (zVal > 1050.0))
    return false;

  return true;
}
