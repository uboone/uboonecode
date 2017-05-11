////////////////////////////////////////////////////////////////////////
// \file SpaceChargeMicroBooNE.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for MicroBooNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <string>
#include <cassert>

// LArSoft includes
#include "uboone/SpaceCharge/SpaceChargeMicroBooNE.h"

// Framework includes
#include "cetlib/exception.h"

// ROOT includes
#include "TFile.h"


namespace {
  
  /// TGraph wrapper with RAII.
  class SortedTGraphFromFile {
    TGraph* graph = nullptr;
    
      public:
    SortedTGraphFromFile(TFile& file, char const* graphName)
      : graph(static_cast<TGraph*>(file.Get(graphName)))
      {
        if (!graph) {
          throw cet::exception("SpaceChargeMicroBooNE")
            << "Graph '" << graphName << "' not found in '" << file.GetPath()
            << "'\n";
        }
        graph->Sort();
      }
    
    ~SortedTGraphFromFile() { delete graph; }
    
    TGraph const& operator* () const { return *graph; }
    
  }; // class SortedTGraphFromFile
  
} // local namespace


//-----------------------------------------------
spacecharge::SpaceChargeMicroBooNE::SpaceChargeMicroBooNE(
  fhicl::ParameterSet const& pset
)
{
  Configure(pset);
}

//------------------------------------------------
bool spacecharge::SpaceChargeMicroBooNE::Configure(fhicl::ParameterSet const& pset)
{  
  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  fEnableCorrSCE = pset.get<bool>("EnableCorrSCE");

  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
  {
    auto const reprTypeString = pset.get<std::string>("RepresentationType");
    fRepresentationType = ParseRepresentationType(reprTypeString);
    if (fRepresentationType == SpaceChargeRepresentation_t::kUnknown) {
      throw cet::exception("SpaceChargeMicroBooNE")
        << "Unknown space charge representation type: '" << reprTypeString
        << "'\n";
    }
    fInputFilename = pset.get<std::string>("InputFilename");

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename,fname);

    TFile infile(fname.c_str(), "READ");
    if(!infile.IsOpen()) throw cet::exception("SpaceChargeMicroBooNE") << "Could not find the space charge effect file '" << fname << "'!\n";

    switch (fRepresentationType) {
      case SpaceChargeRepresentation_t::kParametric:
        for(int i = 0; i < 5; i++)
        {
          g1_x[i] = makeInterpolator(infile, Form("deltaX/g1_%d",i));
          g2_x[i] = makeInterpolator(infile, Form("deltaX/g2_%d",i));
          g3_x[i] = makeInterpolator(infile, Form("deltaX/g3_%d",i));   
          g4_x[i] = makeInterpolator(infile, Form("deltaX/g4_%d",i));
          g5_x[i] = makeInterpolator(infile, Form("deltaX/g5_%d",i));

          g1_y[i] = makeInterpolator(infile, Form("deltaY/g1_%d",i));
          g2_y[i] = makeInterpolator(infile, Form("deltaY/g2_%d",i));
          g3_y[i] = makeInterpolator(infile, Form("deltaY/g3_%d",i));   
          g4_y[i] = makeInterpolator(infile, Form("deltaY/g4_%d",i));
          g5_y[i] = makeInterpolator(infile, Form("deltaY/g5_%d",i));
          g6_y[i] = makeInterpolator(infile, Form("deltaY/g6_%d",i));

          g1_z[i] = makeInterpolator(infile, Form("deltaZ/g1_%d",i));
          g2_z[i] = makeInterpolator(infile, Form("deltaZ/g2_%d",i));
          g3_z[i] = makeInterpolator(infile, Form("deltaZ/g3_%d",i));   
          g4_z[i] = makeInterpolator(infile, Form("deltaZ/g4_%d",i));

          g1_Ex[i] = makeInterpolator(infile, Form("deltaExOverE/g1_%d",i));
          g2_Ex[i] = makeInterpolator(infile, Form("deltaExOverE/g2_%d",i));
          g3_Ex[i] = makeInterpolator(infile, Form("deltaExOverE/g3_%d",i));
          g4_Ex[i] = makeInterpolator(infile, Form("deltaExOverE/g4_%d",i));
          g5_Ex[i] = makeInterpolator(infile, Form("deltaExOverE/g5_%d",i));
          
          g1_Ey[i] = makeInterpolator(infile, Form("deltaEyOverE/g1_%d",i));
          g2_Ey[i] = makeInterpolator(infile, Form("deltaEyOverE/g2_%d",i));
          g3_Ey[i] = makeInterpolator(infile, Form("deltaEyOverE/g3_%d",i));
          g4_Ey[i] = makeInterpolator(infile, Form("deltaEyOverE/g4_%d",i));
          g5_Ey[i] = makeInterpolator(infile, Form("deltaEyOverE/g5_%d",i));
          g6_Ey[i] = makeInterpolator(infile, Form("deltaEyOverE/g6_%d",i));
          
          g1_Ez[i] = makeInterpolator(infile, Form("deltaEzOverE/g1_%d",i));
          g2_Ez[i] = makeInterpolator(infile, Form("deltaEzOverE/g2_%d",i));
          g3_Ez[i] = makeInterpolator(infile, Form("deltaEzOverE/g3_%d",i));
          g4_Ez[i] = makeInterpolator(infile, Form("deltaEzOverE/g4_%d",i));
        }

        g1_x[5] = makeInterpolator(infile, "deltaX/g1_5");
        g2_x[5] = makeInterpolator(infile, "deltaX/g2_5");
        g3_x[5] = makeInterpolator(infile, "deltaX/g3_5");   
        g4_x[5] = makeInterpolator(infile, "deltaX/g4_5");
        g5_x[5] = makeInterpolator(infile, "deltaX/g5_5");

        g1_y[5] = makeInterpolator(infile, "deltaY/g1_5");
        g2_y[5] = makeInterpolator(infile, "deltaY/g2_5");
        g3_y[5] = makeInterpolator(infile, "deltaY/g3_5");   
        g4_y[5] = makeInterpolator(infile, "deltaY/g4_5");
        g5_y[5] = makeInterpolator(infile, "deltaY/g5_5");
        g6_y[5] = makeInterpolator(infile, "deltaY/g6_5");
        
        g1_x[6] = makeInterpolator(infile, "deltaX/g1_6");
        g2_x[6] = makeInterpolator(infile, "deltaX/g2_6");
        g3_x[6] = makeInterpolator(infile, "deltaX/g3_6");
        g4_x[6] = makeInterpolator(infile, "deltaX/g4_6");
        g5_x[6] = makeInterpolator(infile, "deltaX/g5_6");

        g1_Ex[5] = makeInterpolator(infile, "deltaExOverE/g1_5");
        g2_Ex[5] = makeInterpolator(infile, "deltaExOverE/g2_5");
        g3_Ex[5] = makeInterpolator(infile, "deltaExOverE/g3_5");
        g4_Ex[5] = makeInterpolator(infile, "deltaExOverE/g4_5");
        g5_Ex[5] = makeInterpolator(infile, "deltaExOverE/g5_5");
        
        g1_Ey[5] = makeInterpolator(infile, "deltaEyOverE/g1_5");
        g2_Ey[5] = makeInterpolator(infile, "deltaEyOverE/g2_5");
        g3_Ey[5] = makeInterpolator(infile, "deltaEyOverE/g3_5");
        g4_Ey[5] = makeInterpolator(infile, "deltaEyOverE/g4_5");
        g5_Ey[5] = makeInterpolator(infile, "deltaEyOverE/g5_5");
        g6_Ey[5] = makeInterpolator(infile, "deltaEyOverE/g6_5");
        
        g1_Ex[6] = makeInterpolator(infile, "deltaExOverE/g1_6");
        g2_Ex[6] = makeInterpolator(infile, "deltaExOverE/g2_6");
        g3_Ex[6] = makeInterpolator(infile, "deltaExOverE/g3_6");
        g4_Ex[6] = makeInterpolator(infile, "deltaExOverE/g4_6");
        g5_Ex[6] = makeInterpolator(infile, "deltaExOverE/g5_6");
        break; // kParametric
      case SpaceChargeRepresentation_t::kUnknown:
        assert(false); // logic error
        break;
    } // switch

    infile.Close();
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

//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in ionization electron drift
geo::Vector_t spacecharge::SpaceChargeMicroBooNE::GetPosOffsets(geo::Point_t const& point) const
{
  if (!EnableSimSpatialSCE()) return {}; // no correction, zero displacement
  if(!IsInsideBoundaries(point)) return {}; // zero-initialised
  switch (fRepresentationType) {
    
    case SpaceChargeRepresentation_t::kParametric:
      return GetPosOffsetsParametric(point);
    
    case SpaceChargeRepresentation_t::kUnknown:
      assert(false); // logic error: can't be unknown
      return {}; // zero-initialised
    
  } // switch

  assert(false); // logic error: can't be here
  return {};
}

//----------------------------------------------------------------------------
/// Provides position offsets using a parametric representation
geo::Vector_t spacecharge::SpaceChargeMicroBooNE::GetPosOffsetsParametric
  (geo::Point_t const& point) const
{
  geo::Point_t transformedPoint = Transform(point);
  return {
    GetOnePosOffsetParametricX(transformedPoint),
    GetOnePosOffsetParametricY(transformedPoint),
    GetOnePosOffsetParametricZ(transformedPoint)
    };
}

//----------------------------------------------------------------------------
double spacecharge::SpaceChargeMicroBooNE::GetOnePosOffsetParametricX
  (geo::Point_t const& point) const
{
  
  double parA[5][7];
  double const zValNew = point.Z();
  for(size_t j = 0; j < 7; j++) 
  {
    parA[0][j] = g1_x[j].Eval(zValNew);
    parA[1][j] = g2_x[j].Eval(zValNew);
    parA[2][j] = g3_x[j].Eval(zValNew);
    parA[3][j] = g4_x[j].Eval(zValNew);
    parA[4][j] = g5_x[j].Eval(zValNew);
  }

  f1_x_poly_t f1_x;
  f2_x_poly_t f2_x;
  f3_x_poly_t f3_x;
  f4_x_poly_t f4_x;
  f5_x_poly_t f5_x;
  
  f1_x.SetParameters(parA[0]);
  f2_x.SetParameters(parA[1]);
  f3_x.SetParameters(parA[2]);
  f4_x.SetParameters(parA[3]);
  f5_x.SetParameters(parA[4]);

  double const aValNew = point.Y();
  double const parB[] = {
    f1_x.Eval(aValNew),
    f2_x.Eval(aValNew),
    f3_x.Eval(aValNew),
    f4_x.Eval(aValNew),
    f5_x.Eval(aValNew)
  };
  
  double const bValNew = point.X();
  return 100.0*fFinal_x_poly_t::Eval(bValNew, parB);
}

//----------------------------------------------------------------------------
/// Provides one position offset using a parametric representation, for Y axis
double spacecharge::SpaceChargeMicroBooNE::GetOnePosOffsetParametricY
  (geo::Point_t const& point) const
{
  
  double parA[6][6];
  double const zValNew = point.Z();
  for(size_t j = 0; j < 6; j++)
  {
    parA[0][j] = g1_y[j].Eval(zValNew);
    parA[1][j] = g2_y[j].Eval(zValNew);
    parA[2][j] = g3_y[j].Eval(zValNew);
    parA[3][j] = g4_y[j].Eval(zValNew);
    parA[4][j] = g5_y[j].Eval(zValNew);
    parA[5][j] = g6_y[j].Eval(zValNew);
  }
  
  f1_y_poly_t f1_y;
  f2_y_poly_t f2_y;
  f3_y_poly_t f3_y;
  f4_y_poly_t f4_y;
  f5_y_poly_t f5_y;
  f6_y_poly_t f6_y;
  
  f1_y.SetParameters(parA[0]);
  f2_y.SetParameters(parA[1]);
  f3_y.SetParameters(parA[2]);
  f4_y.SetParameters(parA[3]);
  f5_y.SetParameters(parA[4]);
  f6_y.SetParameters(parA[5]);
  
  double const aValNew = point.X();
  
  double const parB[] = {
    f1_y.Eval(aValNew),
    f2_y.Eval(aValNew),
    f3_y.Eval(aValNew),
    f4_y.Eval(aValNew),
    f5_y.Eval(aValNew),
    f6_y.Eval(aValNew)
  };
  
  double const bValNew = point.Y();
  return 100.0*fFinal_y_poly_t::Eval(bValNew, parB);
} // spacecharge::SpaceChargeMicroBooNE::GetOnePosOffsetParametricY()


//----------------------------------------------------------------------------
/// Provides one position offset using a parametric representation, for a given
/// axis
double spacecharge::SpaceChargeMicroBooNE::GetOnePosOffsetParametricZ
  (geo::Point_t const& point) const
{
  
  double parA[4][5];
  double const zValNew = point.Z();
  for(size_t j = 0; j < 5; j++)
  {
    parA[0][j] = g1_z[j].Eval(zValNew);
    parA[1][j] = g2_z[j].Eval(zValNew);
    parA[2][j] = g3_z[j].Eval(zValNew);
    parA[3][j] = g4_z[j].Eval(zValNew);
  }

  f1_z_poly_t f1_z;
  f2_z_poly_t f2_z;
  f3_z_poly_t f3_z;
  f4_z_poly_t f4_z;
  
  f1_z.SetParameters(parA[0]);
  f2_z.SetParameters(parA[1]);
  f3_z.SetParameters(parA[2]);
  f4_z.SetParameters(parA[3]);
  
  double const aValNew = point.Y();
  double const parB[] = {
    f1_z.Eval(aValNew),
    f2_z.Eval(aValNew),
    f3_z.Eval(aValNew),
    f4_z.Eval(aValNew)
  };
  
  double const bValNew = point.X();
  return 100.0*fFinal_z_poly_t::Eval(bValNew, parB);
  
} // spacecharge::SpaceChargeMicroBooNE::GetOnePosOffsetParametricZ()


//----------------------------------------------------------------------------
/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.)
geo::Vector_t spacecharge::SpaceChargeMicroBooNE::GetEfieldOffsets(geo::Point_t const& point) const
{
  if (!EnableSimEfieldSCE()) return {}; // no correction, zero distortion
  if(!IsInsideBoundaries(point)) return {}; // zero-initialised
  
  switch (fRepresentationType) {
    
    case SpaceChargeRepresentation_t::kParametric:
      return -GetEfieldOffsetsParametric(point);
    
    case SpaceChargeRepresentation_t::kUnknown:
      assert(false); // logic error: can't be unknown
      return {};
    
  } // switch
  
  assert(false); // logic error: can't get here
  return {};
}

//----------------------------------------------------------------------------
/// Provides E field offsets using a parametric representation
geo::Vector_t spacecharge::SpaceChargeMicroBooNE::GetEfieldOffsetsParametric
  (geo::Point_t const& point) const
{
  geo::Point_t const transformedPoint = Transform(point);
  return {
    GetOneEfieldOffsetParametricX(transformedPoint),
    GetOneEfieldOffsetParametricY(transformedPoint),
    GetOneEfieldOffsetParametricZ(transformedPoint)
    };
}

//----------------------------------------------------------------------------
/// Provides one E field offset using a parametric representation, for x axis
double spacecharge::SpaceChargeMicroBooNE::GetOneEfieldOffsetParametricX
  (geo::Point_t const& point) const
{
  
  double parA[5][7];
  double const zValNew = point.Z();
  for(size_t j = 0; j < 7; j++)
  {
    parA[0][j] = g1_Ex[j].Eval(zValNew);
    parA[1][j] = g2_Ex[j].Eval(zValNew);
    parA[2][j] = g3_Ex[j].Eval(zValNew);
    parA[3][j] = g4_Ex[j].Eval(zValNew);
    parA[4][j] = g5_Ex[j].Eval(zValNew);
  }

  f1_Ex_poly_t f1_Ex;
  f2_Ex_poly_t f2_Ex;
  f3_Ex_poly_t f3_Ex;
  f4_Ex_poly_t f4_Ex;
  f5_Ex_poly_t f5_Ex;
  
  f1_Ex.SetParameters(parA[0]);
  f2_Ex.SetParameters(parA[1]);
  f3_Ex.SetParameters(parA[2]);
  f4_Ex.SetParameters(parA[3]);
  f5_Ex.SetParameters(parA[4]);
  
  double const aValNew = point.Y();
  double const parB[] = {
    f1_Ex.Eval(aValNew),
    f2_Ex.Eval(aValNew),
    f3_Ex.Eval(aValNew),
    f4_Ex.Eval(aValNew),
    f5_Ex.Eval(aValNew)
  };
  
  double const bValNew = point.X();
  return fFinal_Ex_poly_t::Eval(bValNew, parB);
  
} // spacecharge::SpaceChargeMicroBooNE::GetOneEfieldOffsetParametricX()

//----------------------------------------------------------------------------
/// Provides one E field offset using a parametric representation, for y axis
double spacecharge::SpaceChargeMicroBooNE::GetOneEfieldOffsetParametricY
  (geo::Point_t const& point) const
{
  
  double parA[6][6];
  double const zValNew = point.Z();
  for(size_t j = 0; j < 6; j++)
  {
    parA[0][j] = g1_Ey[j].Eval(zValNew);
    parA[1][j] = g2_Ey[j].Eval(zValNew);
    parA[2][j] = g3_Ey[j].Eval(zValNew);
    parA[3][j] = g4_Ey[j].Eval(zValNew);
    parA[4][j] = g5_Ey[j].Eval(zValNew);
    parA[5][j] = g6_Ey[j].Eval(zValNew);
  }

  f1_Ey_poly_t f1_Ey;
  f2_Ey_poly_t f2_Ey;
  f3_Ey_poly_t f3_Ey;
  f4_Ey_poly_t f4_Ey;
  f5_Ey_poly_t f5_Ey;
  f6_Ey_poly_t f6_Ey;
  
  f1_Ey.SetParameters(parA[0]);
  f2_Ey.SetParameters(parA[1]);
  f3_Ey.SetParameters(parA[2]);
  f4_Ey.SetParameters(parA[3]);
  f5_Ey.SetParameters(parA[4]);
  f6_Ey.SetParameters(parA[5]);

  double const aValNew = point.X();
  double const parB[] = {
    f1_Ey.Eval(aValNew),
    f2_Ey.Eval(aValNew),
    f3_Ey.Eval(aValNew),
    f4_Ey.Eval(aValNew),
    f5_Ey.Eval(aValNew),
    f6_Ey.Eval(aValNew)
  };
  
  double const bValNew = point.Y();
  return fFinal_Ey_poly_t::Eval(bValNew, parB);
  
} // spacecharge::SpaceChargeMicroBooNE::GetOneEfieldOffsetParametricY()


//----------------------------------------------------------------------------
/// Provides one E field offset using a parametric representation, for z axis
double spacecharge::SpaceChargeMicroBooNE::GetOneEfieldOffsetParametricZ
  (geo::Point_t const& point) const
{
  
  double parA[4][5];
  double const zValNew = point.Z();
  for(size_t j = 0; j < 5; j++)
  {
    parA[0][j] = g1_Ez[j].Eval(zValNew);
    parA[1][j] = g2_Ez[j].Eval(zValNew);
    parA[2][j] = g3_Ez[j].Eval(zValNew);
    parA[3][j] = g4_Ez[j].Eval(zValNew);
  }

  f1_Ez_poly_t f1_Ez;
  f2_Ez_poly_t f2_Ez;
  f3_Ez_poly_t f3_Ez;
  f4_Ez_poly_t f4_Ez;
  
  f1_Ez.SetParameters(parA[0]);
  f2_Ez.SetParameters(parA[1]);
  f3_Ez.SetParameters(parA[2]);
  f4_Ez.SetParameters(parA[3]);

  double const aValNew = point.Y();
  double const parB[] = {
    f1_Ez.Eval(aValNew),
    f2_Ez.Eval(aValNew),
    f3_Ez.Eval(aValNew),
    f4_Ez.Eval(aValNew)
  };
  
  double const bValNew = point.X();
  return fFinal_Ez_poly_t::Eval(bValNew, parB);
} // spacecharge::SpaceChargeMicroBooNE::GetOneEfieldOffsetParametricZ()


//----------------------------------------------------------------------------
/// Transform X to SCE X coordinate:  [2.56,0.0] --> [0.0,2.50]
double spacecharge::SpaceChargeMicroBooNE::TransformX(double xVal) const
{
  double xValNew;
  xValNew = 2.50 - (2.50/2.56)*(xVal/100.0);
  xValNew -= 1.25;

  return xValNew;
}

//----------------------------------------------------------------------------
/// Transform Y to SCE Y coordinate:  [-1.165,1.165] --> [0.0,2.50]
double spacecharge::SpaceChargeMicroBooNE::TransformY(double yVal) const
{
  double yValNew;
  yValNew = (2.50/2.33)*((yVal/100.0)+1.165);
  yValNew -= 1.25;

  return yValNew;
}

//----------------------------------------------------------------------------
/// Transform Z to SCE Z coordinate:  [0.0,10.37] --> [0.0,10.0]
double spacecharge::SpaceChargeMicroBooNE::TransformZ(double zVal) const
{
  double zValNew;
  zValNew = (10.0/10.37)*(zVal/100.0);

  return zValNew;
}

//----------------------------------------------------------------------------
geo::Point_t spacecharge::SpaceChargeMicroBooNE::Transform
  (geo::Point_t const& point) const
{
  return
    { TransformX(point.X()), TransformY(point.Y()), TransformZ(point.Z()) };
}

//----------------------------------------------------------------------------
/// Check to see if point is inside boundaries of map (allow to go slightly out of range)
bool spacecharge::SpaceChargeMicroBooNE::IsInsideBoundaries(geo::Point_t const& point) const
{
  return !(
       (point.X() <    0.0) || (point.X() >  260.0)
    || (point.Y() < -120.0) || (point.Y() >  120.0)
    || (point.Z() <    0.0) || (point.Z() > 1050.0)
    );
}


//----------------------------------------------------------------------------

spacecharge::SpaceChargeMicroBooNE::SpaceChargeRepresentation_t
spacecharge::SpaceChargeMicroBooNE::ParseRepresentationType
  (std::string repr_str)
{
  
  if (repr_str == "Parametric") return SpaceChargeRepresentation_t::kParametric;
  else                          return SpaceChargeRepresentation_t::kUnknown;
  
} // spacecharge::SpaceChargeMicroBooNE::ParseRepresentationType()


//----------------------------------------------------------------------------
gsl::Interpolator spacecharge::SpaceChargeMicroBooNE::makeInterpolator
  (TFile& file, char const* graphName)
{
  SortedTGraphFromFile graph(file, graphName);
  return gsl::Interpolator(*graph);
} // spacecharge::SpaceChargeMicroBooNE::makeInterpolator()


//----------------------------------------------------------------------------
