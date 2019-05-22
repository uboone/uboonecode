/**
 * \file NuMuCCNpEventSelection.h
 *
 * \ingroup UBXSec
 * 
 * \brief Class def header for a class NuMuCCNpEventSelection
 *
 * @author Simone Marcocci
 */

/** \addtogroup UBXSec

    @{*/
#ifndef NUMUCCNPEVENTSELECTION_H
#define NUMUCCNPEVENTSELECTION_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "uboone/UBXSec/DataTypes/UBXSecEvent.h"

namespace ubana{
  
  /**
   \class NuMuCCNpEventSelection
   User defined class NuMuCCNpEventSelection ... these comments are used to generate
   doxygen documentation!
 */

  class NuMuCCNpEventSelection {
    
  public:
    
    /// Default constructor
    NuMuCCNpEventSelection();

    /// Default destructor
    ~NuMuCCNpEventSelection(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Print the current configuration
    void PrintConfig();

    /// Set the UBXSecEvent
    void SetEvent(UBXSecEvent*);

    /// Returns true if this event is selected
    bool IsSelected(int & slice_index, std::map<std::string,bool> & failure_map);

    /// Initialises the failure map
    void InitializeFailureMap( std::map<std::string,bool>&  );

    /// This function returns if a 3D point is within the fiducial volume
    bool inCV(float x, float y, float z);
  protected:

    UBXSecEvent* _ubxsec_event;
    bool _event_is_set = false;
    bool _configured = false;

    bool _verbose;

    float _FVx = 256.35;
    float _FVy = 233;
    float _FVz = 1036.8;


    double _FVborderx;
    double _FVbordery;
    double _FVborderz;
    int _ntracks;
    bool _showerastrack;
    unsigned _collection_hits;
    double _proton_p_threshold;
    double _chi2_cut;

  };
}

#endif
/** @} */ // end of doxygen group 

