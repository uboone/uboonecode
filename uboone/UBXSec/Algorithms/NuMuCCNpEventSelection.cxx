#ifndef NUMUCCNPEVENTSELECTION_CXX
#define NUMUCCNPEVENTSELECTION_CXX

#include "NuMuCCNpEventSelection.h"
#include <iostream>

namespace ubana {

  NuMuCCNpEventSelection::NuMuCCNpEventSelection() {
    _configured = false;
  }

  void NuMuCCNpEventSelection::Configure(fhicl::ParameterSet const& pset)
  {
    
    _FVborderx                         = pset.get< double > ( "FVborderx", 10. );
    _FVbordery                         = pset.get< double > ( "FVborderx", 20. );
    _FVborderz                         = pset.get< double > ( "FVborderx", 10. );
    _ntracks                           = pset.get< int    > ( "NTracks" );
    _showerastrack                     = pset.get< bool   > ( "ShowerAsTrack", true );
    _collection_hits                   = pset.get< unsigned    > ( "MinCollectionHits", 5 );
    _proton_p_threshold                = pset.get< double > ( "ProtonPThreshold", 0.3 );
    _chi2_cut                          = pset.get< double > ( "Chi2Cut", 88 );

    _verbose                           = pset.get< bool > ( "Verbose" );

    _configured = true;
  }

  void NuMuCCNpEventSelection::PrintConfig() {

    std::cout << "--- NuMuCCNpEventSelection configuration:" << std::endl;
    std::cout << "---   _FVborderx                         = " << _FVborderx << std::endl;
    std::cout << "---   _FVbordery                         = " << _FVbordery << std::endl;
    std::cout << "---   _FVborderz                         = " << _FVborderz << std::endl;
    std::cout << "---   _tracks                            = " << _ntracks << std::endl;
    std::cout << "---   _showerastrack                     = " << _showerastrack << std::endl;
    std::cout << "---   _collection_hits                   = " << _collection_hits << std::endl;
    std::cout << "---   _proton_p_threshold                = " << _proton_p_threshold << std::endl;
    std::cout << "---   _chi2_cut                          = " << _chi2_cut << std::endl;
    
    std::cout << "---   _verbose                           = " << _verbose << std::endl;
  }

  void NuMuCCNpEventSelection::SetEvent(UBXSecEvent* ubxsec_event) {

    _ubxsec_event = ubxsec_event;
    _event_is_set = true;

  }

  void NuMuCCNpEventSelection::InitializeFailureMap( std::map<std::string,bool> & failure_map ) {

	    failure_map["n_tracks"] = false;
	    failure_map["uncontained_protons"] = false;
	    failure_map["no_protons"] = false;
	    failure_map["low_collection_hits"] = false;
	    failure_map["low_p_proton"] = false;
	    failure_map["fail_chi2"] = false;

  }

  bool NuMuCCNpEventSelection::inCV(float x, float y, float z) { //this is NOT the CC inclusive FV
    if(x < (_FVx - _FVborderx) && (x > _FVborderx) && (y < (_FVy/2. - _FVbordery)) && (y > (-_FVy/2. + _FVbordery)) && (z < (_FVz - _FVborderz)) && (z > _FVborderz)) return true;
    else return false;
  }

  bool NuMuCCNpEventSelection::IsSelected(int & slice_index, std::map<std::string,bool> & failure_map) {
   
    if (_verbose) std::cout << "[NuMuCCNpEventSelection] Starts. " << std::endl;
    
    //set failure map
    InitializeFailureMap( failure_map );

    if (!_configured || !_event_is_set) {
      std::cout << "Call to NuMuCCNpEventSelection::IsSelected() without having done configuration or having set the UBXSecEvent" << std::endl;
      throw std::exception();
    }

    slice_index = _ubxsec_event->selected_slice;
    std::string reason = "no_failure";

    // ************ 
    // Number of tracks
    // ************

    bool return_value = true;
    int number = 0;
    if(_showerastrack) number =  _ubxsec_event->num_pfp;
    else number = _ubxsec_event->num_pfp_tracks;

    if (number < _ntracks) return_value = false;

    if ( !return_value ) {
	    reason = "n_tracks";
      
      if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection FAILED. Number of tracks is " << number << "<" << _ntracks << "." << std::endl;
      	return return_value;
    
    } else {
	    failure_map["n_tracks"] = true;
    }
    
    
    // ************ 
    // All protons contained in FV
    // ************

    int uncontained_proton=0;
    
    // loop over all the proton candidates and select the events with all the proton candidates 
    for( size_t n_trk_pfp=0; n_trk_pfp<_ubxsec_event->pfp_reco_startx.size(); n_trk_pfp++){
       if ( _ubxsec_event->pfp_reco_ismuoncandidate[n_trk_pfp] ) continue;
       
       if(_showerastrack){
         if( _ubxsec_event->pfp_reco_numtracks[n_trk_pfp] !=1 ) continue;
       } else{
         if( !_ubxsec_event->pfp_reco_istrack[n_trk_pfp] ) continue;  
       }
       if(!inCV(_ubxsec_event->pfp_reco_startx[n_trk_pfp], _ubxsec_event->pfp_reco_starty[n_trk_pfp], _ubxsec_event->pfp_reco_startz[n_trk_pfp]) ||
           !inCV(_ubxsec_event->pfp_reco_endx[n_trk_pfp],   _ubxsec_event->pfp_reco_endy[n_trk_pfp],   _ubxsec_event->pfp_reco_endz[n_trk_pfp]))   {
           uncontained_proton +=1;
       }
    }
   
    if( uncontained_proton>0 ) {
	    reason = "uncontained_protons";
      
      if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection FAILED. Number of uncontained protons is " << uncontained_proton << "." << std::endl;
      	return false;
    
    } else {
	    failure_map["uncontained_protons"] = true;
    }
   
    // ************ 
    // Find the leading proton candidate and perform the minCol cut
    // ************
    int pind=-999;
    float temp_length=-999.0;
    
    for(size_t np=0; np<_ubxsec_event->pfp_reco_length.size(); np++){
        if( _ubxsec_event->pfp_reco_ismuoncandidate[np] )  continue;
        if( !_ubxsec_event->pfp_reco_istrack[np] && !_ubxsec_event->pfp_reco_isshower[np] ) continue;      
        
        if(_showerastrack){
           if(_ubxsec_event->pfp_reco_numtracks[np] !=1) continue;
        } else {
           if( !_ubxsec_event->pfp_reco_istrack[np] ) continue; 
        }
        if(_ubxsec_event->pfp_reco_length[np] > temp_length) {
           temp_length=_ubxsec_event->pfp_reco_length[np];
           pind=np;
        }
    }
     
    if( pind < 0 ) {
	    reason = "no_protons";
      
      if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection FAILED. There are no candidate protons in the event." << std::endl;
      	return false;
    
    } else {
	    failure_map["no_protons"] = true;
    }
   
    if( _ubxsec_event->pfp_reco_dEdx[pind].size()< _collection_hits ) { //_collection_hits = 5
	    reason = "low_collection_hits";
      
      if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection FAILED. Number of collection hits is " << _ubxsec_event->pfp_reco_dEdx[pind].size()  << " < " << _collection_hits << std::endl;
      	return false;
    
    } else {
	    failure_map["low_collection_hits"] = true;
    }
   
    
    // ************ 
    // Threshold on leading proton momentum
    // ************

    if( _ubxsec_event->pfp_reco_Mom_proton[pind]<_proton_p_threshold ) { //_proton_p_threshold = 300 GeV/c
	    reason = "low_p_proton";
      
      if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection FAILED. The momentum of the leading proton is " << _ubxsec_event->pfp_reco_Mom_proton[pind]  << " < " << _proton_p_threshold << std::endl;
      	return false;
    
    } else {
	    failure_map["low_p_proton"] = true;
    }
   
    
    // ************ 
    // Perform the chi2 cut for selecting protons
    // ************

     
    Int_t npcand_fail_chi2=0;
    for(size_t ntrk=0; ntrk<_ubxsec_event->pfp_reco_chi2_proton.size(); ntrk++){
        if( _ubxsec_event->pfp_reco_ismuoncandidate[ntrk] ) continue;
        if( !_ubxsec_event->pfp_reco_istrack[ntrk] && !_ubxsec_event->pfp_reco_isshower[ntrk]) continue;
        
	if(_showerastrack){
          if(_ubxsec_event->pfp_reco_numtracks[ntrk] !=1) continue;
        } else {
          if( !_ubxsec_event->pfp_reco_istrack[ntrk] ) continue;
        }
        
	if(_ubxsec_event->pfp_reco_dEdx[ntrk].size()>=_collection_hits && _ubxsec_event->pfp_reco_chi2_proton[ntrk]>_chi2_cut) {
		npcand_fail_chi2++;
	 	if (_verbose) std::cout << "[NuMuCCNpEventSelection] Chi2 is " << _ubxsec_event->pfp_reco_chi2_proton[ntrk] << std::endl;
	} //chi2_cut=88
    }
    
    if( npcand_fail_chi2>0 ) { 
	    reason = "fail_chi2";
      
      if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection FAILED. Number of candidate protons not passing chi2: " << npcand_fail_chi2 << std::endl;
      	return false;
    
    } else {
	    failure_map["fail_chi2"] = true;
    }

    for (auto iter : failure_map) {
      if (!iter.second) {
        if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection FAILED. Did not pass " << iter.first << " cut." << std::endl;
        return false;
      }
    }

    if (_verbose) std::cout << "[NuMuCCNpEventSelection] Selection PASSED. All cuts are satisfied." << std::endl;

    return true;

  }
}


#endif
