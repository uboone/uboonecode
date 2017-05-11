////////////////////////////////////////////////////////////////////////
// \file SpaceChargeMicroBooNE.h
//
// \brief header of class for storing/accessing space charge distortions for MicroBooNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGE_SPACECHARGEMICROBOONE_H
#define SPACECHARGE_SPACECHARGEMICROBOONE_H

// LArSoft libraries
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"

// GNU Scientific library
#include "gsl/gsl_interp.h" // gsl_interp, gsl_interp_linear, ...
#include "gsl/gsl_poly.h" // gsl_poly_eval()

// ROOT includes
#include "TGraph.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <array>
#include <algorithm> // std::copy_n()


// forward declarations
class TFile;


namespace gsl {
  
  /// Wrapper to GSL interpolator, with RAII and data extraction from TGraph.
  class Interpolator {
    
    std::vector<double> xa, ya; ///< Input data.
    gsl_interp* fInterp = nullptr; ///< Pointer to the GSL interpolator.
    
    void swap(Interpolator& with)
      {
        std::swap(xa, with.xa);
        std::swap(ya, with.ya);
        std::swap(fInterp, with.fInterp);
      }
    
      public:
    
    /// Default constructor: object will be useless until reassigned.
    Interpolator() = default;
    
    /// Initializes from the data of a TGraph. The graph is assumed sorted.
    Interpolator(TGraph const& graph)
      : xa(graph.GetX(), graph.GetX() + graph.GetN())
      , ya(graph.GetY(), graph.GetY() + graph.GetN())
      , fInterp(gsl_interp_alloc(gsl_interp_linear, xa.size()))
      {
        // assert(graph.TestBit(TGraph::kIsSortedX));
        gsl_interp_init(fInterp, xa.data(), ya.data(), xa.size());
      }
    
    ~Interpolator()
      { if (fInterp) { gsl_interp_free(fInterp); fInterp = nullptr; } }
    
    Interpolator(Interpolator const&) = delete;
    Interpolator(Interpolator&& from) { swap(from); }
    Interpolator& operator= (Interpolator const&) = delete;
    Interpolator& operator= (Interpolator&& from)
      { swap(from); return *this; }
    
    /// Interpolates to `x`, capping it into x range limits.
    double Eval(double x) const
      {
        if (x <= xa.front()) return ya.front();
        if (x >= xa.back()) return ya.back();
        // use no accelerator
        return gsl_interp_eval(fInterp, xa.data(), ya.data(), x, nullptr);
      }
    
  }; // class Interpolator
  
  
  /// GSL polynomial evaluation wrapper.
  template <unsigned int N>
  struct PolynomialBase {
    static constexpr unsigned int Degree = N; ///< Degree of the polynomial.
    
    /// Number of parameters.
    static constexpr unsigned int NParams = Degree + 1;
    
    using Params_t = std::array<double, NParams>; ///< Set of parameters.
    
    /// Evaluates the polynomial with specified parameters at point `x`.
    static double Eval(double x, double const* params);
    
    /// Evaluates the polynomial with specified parameters at point `x`.
    static double Eval(double x, Params_t const& params)
      { return Eval(x, params.data()); }
    
  }; // class PolynomialBase<>
  
  
  /// Implementation of polynomial evaluations mimicking some of TF1 I/F.
  template <unsigned int N>
  class Polynomial: public PolynomialBase<N> {
    using PolyImpl_t = PolynomialBase<N>;
      public:
    static constexpr unsigned int Degree = PolyImpl_t::Degree;
    static constexpr unsigned int NParams = PolyImpl_t::NParams;
    
    using Params_t = typename PolyImpl_t::Params_t;
    
    /// Default constructor: coefficients are *not initialised*.
    Polynomial() = default;
    
    /// Copies the specified parameters into this object.
    void SetParameters(double const* params)
      { std::copy_n(params, NParams, fParams.begin()); }
    
    
    using PolynomialBase<N>::Eval;
    
    /// Evaluates the polynomial at the specified point `x`.
    double Eval(double x) const
      { return PolyImpl_t::Eval(x, fParams.data()); }
    
      protected:
    Params_t fParams; ///< Parameters of the polynomial.
    
  }; // class Polynomial<>
  
} // namespace gsl


namespace spacecharge {
  
  
  
  class SpaceChargeMicroBooNE : public SpaceCharge {
 
    public:

      typedef enum {
        kParametric,
        kUnknown
      } SpaceChargeRepresentation_t;
    
      explicit SpaceChargeMicroBooNE(fhicl::ParameterSet const& pset);
      SpaceChargeMicroBooNE(SpaceChargeMicroBooNE const&) = delete;
      virtual ~SpaceChargeMicroBooNE() = default;
      
      bool Configure(fhicl::ParameterSet const& pset);
      bool Update(uint64_t ts=0);
      
      bool EnableSimSpatialSCE() const override;
      bool EnableSimEfieldSCE() const override;
      bool EnableCorrSCE() const override;
      geo::Vector_t GetPosOffsets(geo::Point_t const& point) const override;
      geo::Vector_t GetEfieldOffsets(geo::Point_t const& point) const override;
      
    protected:

      static SpaceChargeRepresentation_t ParseRepresentationType
        (std::string repr_str);
      
      geo::Vector_t GetPosOffsetsParametric(geo::Point_t const& point) const;
      double GetOnePosOffsetParametricX(geo::Point_t const& point) const;
      double GetOnePosOffsetParametricY(geo::Point_t const& point) const;
      double GetOnePosOffsetParametricZ(geo::Point_t const& point) const;

      geo::Vector_t GetEfieldOffsetsParametric(geo::Point_t const& point) const;
      double GetOneEfieldOffsetParametricX(geo::Point_t const& point) const;
      double GetOneEfieldOffsetParametricY(geo::Point_t const& point) const;
      double GetOneEfieldOffsetParametricZ(geo::Point_t const& point) const;
      double TransformX(double xVal) const;
      double TransformY(double yVal) const;
      double TransformZ(double zVal) const;
      geo::Point_t Transform(geo::Point_t const& point) const;

      bool IsInsideBoundaries(geo::Point_t const& point) const;

      bool fEnableSimSpatialSCE;
      bool fEnableSimEfieldSCE;
      bool fEnableCorrSCE;
      
      SpaceChargeRepresentation_t fRepresentationType = kUnknown;
      std::string fInputFilename;
      
      std::array<gsl::Interpolator, 7U> g1_x;
      std::array<gsl::Interpolator, 7U> g2_x;
      std::array<gsl::Interpolator, 7U> g3_x;
      std::array<gsl::Interpolator, 7U> g4_x;
      std::array<gsl::Interpolator, 7U> g5_x;
      
      std::array<gsl::Interpolator, 7U> g1_y;
      std::array<gsl::Interpolator, 7U> g2_y;
      std::array<gsl::Interpolator, 7U> g3_y;
      std::array<gsl::Interpolator, 7U> g4_y;
      std::array<gsl::Interpolator, 7U> g5_y;
      std::array<gsl::Interpolator, 7U> g6_y;
      
      std::array<gsl::Interpolator, 7U> g1_z;
      std::array<gsl::Interpolator, 7U> g2_z;
      std::array<gsl::Interpolator, 7U> g3_z;
      std::array<gsl::Interpolator, 7U> g4_z;
      
      typedef gsl::Polynomial<6> f1_x_poly_t;
      typedef gsl::Polynomial<6> f2_x_poly_t;
      typedef gsl::Polynomial<6> f3_x_poly_t;
      typedef gsl::Polynomial<6> f4_x_poly_t;
      typedef gsl::Polynomial<6> f5_x_poly_t;
      typedef gsl::PolynomialBase<4> fFinal_x_poly_t;
      
      typedef gsl::Polynomial<5> f1_y_poly_t;
      typedef gsl::Polynomial<5> f2_y_poly_t;
      typedef gsl::Polynomial<5> f3_y_poly_t;
      typedef gsl::Polynomial<5> f4_y_poly_t;
      typedef gsl::Polynomial<5> f5_y_poly_t;
      typedef gsl::Polynomial<5> f6_y_poly_t;
      typedef gsl::PolynomialBase<5> fFinal_y_poly_t;
      
      typedef gsl::Polynomial<4> f1_z_poly_t;
      typedef gsl::Polynomial<4> f2_z_poly_t;
      typedef gsl::Polynomial<4> f3_z_poly_t;
      typedef gsl::Polynomial<4> f4_z_poly_t;
      typedef gsl::PolynomialBase<3> fFinal_z_poly_t;
    
      std::array<gsl::Interpolator, 7U> g1_Ex;
      std::array<gsl::Interpolator, 7U> g2_Ex;
      std::array<gsl::Interpolator, 7U> g3_Ex;
      std::array<gsl::Interpolator, 7U> g4_Ex;
      std::array<gsl::Interpolator, 7U> g5_Ex;
      
      std::array<gsl::Interpolator, 7U> g1_Ey;
      std::array<gsl::Interpolator, 7U> g2_Ey;
      std::array<gsl::Interpolator, 7U> g3_Ey;
      std::array<gsl::Interpolator, 7U> g4_Ey;
      std::array<gsl::Interpolator, 7U> g5_Ey;
      std::array<gsl::Interpolator, 7U> g6_Ey;
      
      std::array<gsl::Interpolator, 7U> g1_Ez;
      std::array<gsl::Interpolator, 7U> g2_Ez;
      std::array<gsl::Interpolator, 7U> g3_Ez;
      std::array<gsl::Interpolator, 7U> g4_Ez;
      
      typedef gsl::Polynomial<6> f1_Ex_poly_t;
      typedef gsl::Polynomial<6> f2_Ex_poly_t;
      typedef gsl::Polynomial<6> f3_Ex_poly_t;
      typedef gsl::Polynomial<6> f4_Ex_poly_t;
      typedef gsl::Polynomial<6> f5_Ex_poly_t;
      typedef gsl::PolynomialBase<4> fFinal_Ex_poly_t;
      
      typedef gsl::Polynomial<5> f1_Ey_poly_t;
      typedef gsl::Polynomial<5> f2_Ey_poly_t;
      typedef gsl::Polynomial<5> f3_Ey_poly_t;
      typedef gsl::Polynomial<5> f4_Ey_poly_t;
      typedef gsl::Polynomial<5> f5_Ey_poly_t;
      typedef gsl::Polynomial<5> f6_Ey_poly_t;
      typedef gsl::PolynomialBase<5> fFinal_Ey_poly_t;
      
      typedef gsl::Polynomial<4> f1_Ez_poly_t;
      typedef gsl::Polynomial<4> f2_Ez_poly_t;
      typedef gsl::Polynomial<4> f3_Ez_poly_t;
      typedef gsl::Polynomial<4> f4_Ez_poly_t;
      typedef gsl::PolynomialBase<3> fFinal_Ez_poly_t;
    
    
      /// Returns a new GSL interpolator with data from TGraph in specified file.
      static gsl::Interpolator makeInterpolator
        (TFile& file, char const* graphName);
    
  }; // class SpaceChargeMicroBooNE
} //namespace spacecharge



//------------------------------------------------------------------------------
//--- Template implementation
//---
template <unsigned int N>
double gsl::PolynomialBase<N>::Eval
  (double x, double const* params)
{
  // GSL implementation
  return gsl_poly_eval(params, NParams, x);
#if 0
  // hand-craft implementation
  
  // gave up using GSL because I cant figure out how to convince CET build to link to it
  double const* param = params + NParams;
  double p = *--param;
  while (param != params) {
    p *= x;
    p += *--param;
  }
  return p;
#endif // 0
} // gsl::PolynomialBase<N>::Eval()


//------------------------------------------------------------------------------

#endif // SPACECHARGE_SPACECHARGEMICROBOONE_H
