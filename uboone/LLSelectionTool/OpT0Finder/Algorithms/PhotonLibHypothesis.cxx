#ifndef PHOTONLIBHYPOTHESIS_CXX
#define PHOTONLIBHYPOTHESIS_CXX

#include "PhotonLibHypothesis.h"
//#include "OpT0Finder/PhotonLibrary/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

namespace flashana {

  PhotonLibHypothesis::PhotonLibHypothesis(const std::string name)
    : BaseFlashHypothesis(name)
  {}

  //void PhotonLibHypothesis::Configure(const ::fcllite::PSet &pset)
  void PhotonLibHypothesis::Configure(const ::fhicl::ParameterSet &pset)
  {}
  
  void PhotonLibHypothesis::FillEstimate(const QCluster_t& trk,
					 Flash_t &flash) const
  {
    art::ServiceHandle<phot::PhotonVisibilityService> vis;
    double xyz[3]={0.};
    
    size_t n_pmt = BaseAlgorithm::NOpDets();//n_pmt returns 0 now, needs to be fixed
    
    for ( auto& v : flash.pe_v ) v = 0;
    
    for ( size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {

      for ( size_t ipt = 0; ipt < trk.size(); ++ipt) {
	
        auto const& pt = trk[ipt];
	
        double q = pt.q;
	
        //q *= ::phot::PhotonVisibilityService::GetME().GetVisibility( pt.x, pt.y, pt.z, ipmt)*0.0093;
	xyz[0] = pt.x;
	xyz[1] = pt.y;
	xyz[2] = pt.z;
	q *= vis->GetVisibility(xyz,ipmt) * 0.0093;
        flash.pe_v[ipmt] += q;
	//std::cout << "PMT : " << ipmt << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q << std::endl;
	
      }
    }

    return;
  }
}
#endif
