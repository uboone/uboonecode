#ifndef UBOONETPCENERGYCALIBPROVIDER_CXX
#define UBOONETPCENERGYCALIBPROVIDER_CXX

#include "UbooneTPCEnergyCalibProvider.h"

// art/LArSoft libraries
#include "cetlib/exception.h"

#include <fstream>

namespace lariov {

  //constructor
  UbooneTPCEnergyCalibProvider::UbooneTPCEnergyCalibProvider(fhicl::ParameterSet const& p) :
    fXNormProvider(p.get<fhicl::ParameterSet>("XNormCorrectionProvider")),
    fdEdxProvider(p.get<fhicl::ParameterSet>("dEdxCorrectionProvider")) {   

    fXShapeProvider.resize(3, NULL);    
    fXShapeProvider[0].reset(new UboonedqdxCorrectionProvider(p.get<fhicl::ParameterSet>("XShapeCorrectionProvider_Plane0")));
    fXShapeProvider[1].reset(new UboonedqdxCorrectionProvider(p.get<fhicl::ParameterSet>("XShapeCorrectionProvider_Plane1")));
    fXShapeProvider[2].reset(new UboonedqdxCorrectionProvider(p.get<fhicl::ParameterSet>("XShapeCorrectionProvider_Plane2")));
    
    fYZProvider.resize(3, NULL);
    fYZProvider[0].reset(new UboonedqdxCorrectionProvider(p.get<fhicl::ParameterSet>("YZCorrectionProvider_Plane0")));
    fYZProvider[1].reset(new UboonedqdxCorrectionProvider(p.get<fhicl::ParameterSet>("YZCorrectionProvider_Plane1")));
    fYZProvider[2].reset(new UboonedqdxCorrectionProvider(p.get<fhicl::ParameterSet>("YZCorrectionProvider_Plane2")));
  }
  
  void UbooneTPCEnergyCalibProvider::Reconfigure(fhicl::ParameterSet const& p) {
    
    fXShapeProvider[0]->Reconfigure(p.get<fhicl::ParameterSet>("XShapeCorrectionProvider_Plane0"));
    fXShapeProvider[1]->Reconfigure(p.get<fhicl::ParameterSet>("XShapeCorrectionProvider_Plane1"));
    fXShapeProvider[2]->Reconfigure(p.get<fhicl::ParameterSet>("XShapeCorrectionProvider_Plane2"));
    
    fYZProvider[0]->Reconfigure(p.get<fhicl::ParameterSet>("YZCorrectionProvider_Plane0"));
    fYZProvider[1]->Reconfigure(p.get<fhicl::ParameterSet>("YZCorrectionProvider_Plane1"));
    fYZProvider[2]->Reconfigure(p.get<fhicl::ParameterSet>("YZCorrectionProvider_Plane2"));
    
    fXNormProvider.Reconfigure(p.get<fhicl::ParameterSet>("XNormCorrectionProvider"));
    fdEdxProvider.Reconfigure(p.get<fhicl::ParameterSet>("dEdxCorrectionProvider"));
  }
  
  bool UbooneTPCEnergyCalibProvider::Update(DBTimeStamp_t ts) {
    
    //get around compiler optimization shenanigans by putting each update call in a separate line before returning
    bool return_val1 = fXShapeProvider[0]->Update(ts);
    bool return_val2 = fXShapeProvider[1]->Update(ts);
    bool return_val3 = fXShapeProvider[2]->Update(ts);
    bool return_val4 = fYZProvider[0]->Update(ts);
    bool return_val5 = fYZProvider[1]->Update(ts);
    bool return_val6 = fYZProvider[2]->Update(ts);
    bool return_val7 = fXNormProvider.Update(ts);
    bool return_val8 = fdEdxProvider.Update(ts);
    
    return return_val1 || return_val2 || return_val3 || return_val4 || return_val5 || return_val6 || return_val7 || return_val8;
    
  }
  
  float UbooneTPCEnergyCalibProvider::YZdqdxCorrection(int plane, float y, float z) const {
    std::vector<float> vec(2,y);
    vec[1] = z;
    return fYZProvider.at(plane)->Correction(vec);
  }
  
  float UbooneTPCEnergyCalibProvider::YZdqdxCorrectionErr(int plane, float y, float z) const {
    std::vector<float> vec(2,y);
    vec[1] = z;
    return fYZProvider.at(plane)->CorrectionErr(vec);
  }
  
  float UbooneTPCEnergyCalibProvider::DriftdqdxCorrection(int plane, float drift_coord) const {
    return this->XdqdxCorrection(plane, drift_coord);
  }
  
  float UbooneTPCEnergyCalibProvider::DriftdqdxCorrectionErr(int plane, float drift_coord) const {
    return this->XdqdxCorrectionErr(plane, drift_coord);
  }

  float UbooneTPCEnergyCalibProvider::XShapedqdxCorrection(int plane, float x) const {
    std::vector<float> vec(1, x);
    return fXShapeProvider.at(plane)->Correction(vec);
  }
  
  float UbooneTPCEnergyCalibProvider::XShapedqdxCorrectionErr(int plane, float x) const {
    std::vector<float> vec(1, x);
    return fXShapeProvider.at(plane)->CorrectionErr(vec);
  }
  
  float UbooneTPCEnergyCalibProvider::XNormdqdxCorrection(int plane) const {
    return fXNormProvider.Correction(plane);
  }
  
  float UbooneTPCEnergyCalibProvider::XNormdqdxCorrectionErr(int plane) const {
    return fXNormProvider.CorrectionErr(plane);
  }
  
  float UbooneTPCEnergyCalibProvider::XdqdxCorrection(int plane, float x) const {
    return (this->XShapedqdxCorrection(plane,x))*(this->XNormdqdxCorrection(plane));
  }
  
  float UbooneTPCEnergyCalibProvider::XdqdxCorrectionErr(int plane, float x) const {
    float norm_frac_err  = ( this->XNormdqdxCorrectionErr(plane) )     / (this->XNormdqdxCorrection(plane)     );
    float shape_frac_err = ( this->XShapedqdxCorrectionErr(plane, x) ) / (this->XShapedqdxCorrection(plane, x) );
    return (this->XdqdxCorrection(plane,x))*sqrt(norm_frac_err*norm_frac_err + shape_frac_err*shape_frac_err);
  }
  
  float UbooneTPCEnergyCalibProvider::YdqdxCorrection(int plane, float y) const {
    std::cout<<"WARNING: MicroBooNE dqdx calibration does not factorize into separate y and z dependencies"<<std::endl;
    return 1.0;
  }
  
  float UbooneTPCEnergyCalibProvider::YdqdxCorrectionErr(int plane, float y) const {
    std::cout<<"WARNING: MicroBooNE dqdx calibration does not factorize into separate y and z dependencies"<<std::endl;
    return 0.0;
  }
  
  float UbooneTPCEnergyCalibProvider::ZdqdxCorrection(int plane, float z) const {
    std::cout<<"WARNING: MicroBooNE dqdx calibration does not factorize into separate y and z dependencies"<<std::endl;
    return 1.0;
  }
  
  float UbooneTPCEnergyCalibProvider::ZdqdxCorrectionErr(int plane, float y) const {
    std::cout<<"WARNING: MicroBooNE dqdx calibration does not factorize into separate y and z dependencies"<<std::endl;
    return 0.0;
  }
    
  float UbooneTPCEnergyCalibProvider::ThetadqdxCorrection(int plane, float theta) const {
    std::cout<<"WARNING: theta dqdx Correction Not implemented!"<<std::endl;
    return 1.0;
  }
  
  float UbooneTPCEnergyCalibProvider::ThetadqdxCorrectionErr(int plane, float theta) const {
    std::cout<<"WARNING: theta dqdx Correction Not implemented!"<<std::endl;
    return 0.0;
  }
  
  float UbooneTPCEnergyCalibProvider::PhidqdxCorrection(int plane, float phi) const {
    std::cout<<"WARNING: MicroBooNE dqdx calibration does not have a phi dependence"<<std::endl;
    return 1.0;
  }
  
  float UbooneTPCEnergyCalibProvider::PhidqdxCorrectionErr(int plane, float phi) const {
    std::cout<<"WARNING: MicroBooNE dqdx calibration does not have a phi dependence"<<std::endl;
    return 0.0;
  }
  
  float UbooneTPCEnergyCalibProvider::TotaldqdxCorrection(int plane, float x, float y, float z, float theta /*=0.0*/, float phi /*=0.0*/) const {
    return ( this->XdqdxCorrection(plane, x) )*( this->YZdqdxCorrection(plane, y, z) );
  }
  
  float UbooneTPCEnergyCalibProvider::TotaldqdxCorrectionErr(int plane, float x, float y, float z, float theta /*=0.0*/, float phi /*=0.0*/) const {
    float x_frac_err =     ( this->XdqdxCorrectionErr(plane, x) )         / ( this->XdqdxCorrection(plane, x) );
    float yz_frac_err =    ( this->YZdqdxCorrectionErr(plane, y, z) )     / ( this->YZdqdxCorrection(plane, y, z) );
    
    
    return (this->TotaldqdxCorrection(plane, x, y, z, theta, phi))*sqrt(x_frac_err*x_frac_err + yz_frac_err*yz_frac_err);
  }
  
  float UbooneTPCEnergyCalibProvider::dEdxCorrection(int plane) const {
    return fdEdxProvider.Correction(plane);
  }
  
  float UbooneTPCEnergyCalibProvider::dEdxCorrectionErr(int plane) const {
    return fdEdxProvider.CorrectionErr(plane);
  }

} //end namespace lariov

#endif
