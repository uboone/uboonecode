////////////////////////////////////////////////////////////////////////
//
//  MCShowerRecoPart source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCSHOWERRECOPART_CC
#define MCSHOWERRECOPART_CC

#include "MCShowerRecoPart.h"

namespace larreco {

  const unsigned int MCShowerRecoPart::kINVALID_UINT = std::numeric_limits<unsigned int>::max();
  const int MCShowerRecoPart::kINVALID_INT = std::numeric_limits<int>::max();

  //##################################################################
  MCShowerRecoPart::MCShowerRecoPart(fhicl::ParameterSet const& pset)
  //##################################################################
  {
    _debug_mode = pset.get<bool>("DebugMode");
    _comb_time_cut  = pset.get<double>("TimeCut");
    _comb_dist2_cut = pset.get<double>("Dist2Cut");
    _merge = pset.get<bool>("EnableMerge");
  }

  //----------------------------------------------------------------------------------------------- 
  void MCShowerRecoPart::ConstructShower(const art::Handle<std::vector<simb::MCParticle> > mcpArray)
  //-----------------------------------------------------------------------------------------------
  {
    if(!mcpArray.isValid()) 

      throw cet::exception(__FUNCTION__) << "Invalid particle array provided. Nothing done!";
    
    AddParticles(mcpArray);

    ConstructGranularShower();

    if(_merge)

      CombineGranularShower();
    
  }

  //----------------------------------------------------------------
  void MCShowerRecoPart::ClearAndReservePartArray(size_t n)
  //----------------------------------------------------------------
  {
    _track_index.clear();
    _track_id.clear();
    _mother.clear();
    _pdgcode.clear();
    _daughters.clear();
    _shower_id.clear();
    _start_vtx.clear();
    _start_mom.clear();

    _track_id.reserve(n);
    _mother.reserve(n);
    _pdgcode.reserve(n);
    _start_vtx.reserve(n);
    _start_mom.reserve(n);
    _daughters.reserve(n);
    _shower_id.reserve(n);
  }

  //----------------------------------------------------------------
  void MCShowerRecoPart::AddParticle(unsigned int track_id,
				    unsigned int mother_track_id,
				    int pdgcode,
				    const TLorentzVector &start_vtx,
				    const TLorentzVector &start_mom,
				    const std::set<unsigned int> &daughters)
  //----------------------------------------------------------------
  {
    unsigned int index = _track_index.size();

    _track_index.insert(std::pair<unsigned int, unsigned int>(track_id,index));
    _track_id.push_back(track_id);
    _mother.push_back(mother_track_id);
    _pdgcode.push_back(pdgcode);
    _start_vtx.push_back(std::vector<double>(4,0.));
    _start_mom.push_back(std::vector<double>(4,0.));
    _daughters.push_back(std::set<unsigned int>(daughters));

    _start_vtx[index][0]=start_vtx.X();
    _start_vtx[index][1]=start_vtx.Y();
    _start_vtx[index][2]=start_vtx.Z();
    _start_vtx[index][3]=start_vtx.T();

    _start_mom[index][0]=start_mom.X();
    _start_mom[index][1]=start_mom.Y();
    _start_mom[index][2]=start_mom.Z();
    _start_mom[index][3]=start_mom.T();
    
  }
  
  //--------------------------------------------------------------------------------------------
  void MCShowerRecoPart::AddParticles(const art::Handle<std::vector<simb::MCParticle> > mcpArray)
  //--------------------------------------------------------------------------------------------
  {

    //art::Handle<std::vector<simb::MCParticle> > mcpArray;
    //evt.getByLabel(fG4ModName,mcpArray);

    // Read in all particles' information
    ClearAndReservePartArray(mcpArray->size());
    for(size_t i=0; i < mcpArray->size(); ++i) {

      const art::Ptr<simb::MCParticle> mcp_ptr(mcpArray,i);


      std::set<unsigned int> daughters;
      for(size_t i=0; i<(size_t)(mcp_ptr->NumberDaughters()); ++i)
	daughters.insert(mcp_ptr->Daughter(i));

      this->AddParticle((unsigned int)(mcp_ptr->TrackId()),
			(unsigned int)(mcp_ptr->Mother()),
			mcp_ptr->PdgCode(),
			mcp_ptr->Position(),
			mcp_ptr->Momentum(),
			daughters);
    }
  }

  void MCShowerRecoPart::GetTrackStartInfo(const unsigned int &index,
					  double &start_x,
					  double &start_y,
					  double &start_z,
					  double &start_time)
  {
    if(index > _track_id.size())

      throw cet::exception(__FUNCTION__) << Form("Particle index %d not found!",index);

    start_x = _start_vtx.at(index).at(0);
    start_y = _start_vtx.at(index).at(1);
    start_z = _start_vtx.at(index).at(2);
    start_time = _start_vtx.at(index).at(3);

  }

  //----------------------------------------------------------------
  void MCShowerRecoPart::ConstructGranularShower()
  //----------------------------------------------------------------
  {

    if(!_mother.size()) return;

    _shower_id.clear();
    _shower_id.resize(_start_vtx.size(),-1);
    _shower_index.clear();
    _shower_daughters.clear();

    // Make shower-parent candidate look up table
    std::set<unsigned int> shower_parent_id;
    for(size_t i=0; i<_track_id.size(); ++i) {

      if( _pdgcode.at(i) == 22 ||
	  _pdgcode.at(i) == 11 ||
	  _pdgcode.at(i) == -11 ) {

	shower_parent_id.insert(_track_id.at(i));
      }

    }

    std::vector<std::map<Double_t,std::set<unsigned int> > > ordered_shower_daughters;
    for(size_t i=0; i<_mother.size(); ++i) {
      
      unsigned int parent_id  = 0;
      unsigned int grandma_id = _mother.at(i);
      double       time       = _start_vtx.at(i).at(3);

      while(1) {

	if( shower_parent_id.find(grandma_id) == shower_parent_id.end() )

	  break;
	  
	else {

	  parent_id = grandma_id;

	  auto grandma_iter = _track_index.find(grandma_id);

	  if(grandma_iter == _track_index.end()) break;

	  grandma_id = _mother.at((*grandma_iter).second);

	}

      }

      if(!parent_id) continue;
      unsigned int parent_index = (*(_track_index.find(parent_id))).second;

      auto shower_index_iter = _shower_index.find(parent_index);
      size_t shower_index = 0;
      if(shower_index_iter == _shower_index.end()) {

	shower_index = _shower_index.size();
	_shower_index.insert(std::pair<unsigned int, unsigned int>(parent_index,shower_index));
	ordered_shower_daughters.push_back(std::map<double,std::set<unsigned int> >());
      }
      else shower_index = (*shower_index_iter).second;

      // Add to a shower daughter... if there's already a particle @ same time T, shift by 1e-12
      if(ordered_shower_daughters.at(shower_index).find(time) == ordered_shower_daughters.at(shower_index).end() )

	ordered_shower_daughters[shower_index].insert(std::pair<double,std::set<unsigned int> >(time,std::set<unsigned int>()));

      ordered_shower_daughters[shower_index][time].insert(i);
      _shower_id[i]=shower_index;

      // In case this particle has daughters, add them
      for(auto const daughter_track : _daughters.at(i)) {

	auto const daughter_index_iter = _track_index.find(daughter_track);
	if( daughter_index_iter == _track_index.end() ) continue;
	
	ordered_shower_daughters[shower_index][time].insert((*daughter_index_iter).second);
	_shower_id[(*daughter_index_iter).second]=shower_index;
      }
    }

    // Store ordered daughters' list
    _shower_daughters.clear();
    _shower_daughters.reserve(ordered_shower_daughters.size());
    for(auto time_daughters : ordered_shower_daughters) {
      _shower_daughters.push_back(std::vector<unsigned int>());
      auto shower_daughters = _shower_daughters.rbegin();
      shower_daughters->reserve(time_daughters.size());
      std::set<unsigned int> unique_daughter_set;
      for(auto const time_daughter_pair : time_daughters){
	for(auto const daughter_index : time_daughter_pair.second){
	  if(unique_daughter_set.find(daughter_index)==unique_daughter_set.end()) {
	    unique_daughter_set.insert(daughter_index);
	    shower_daughters->push_back(daughter_index);
	  }
	}
      }
    }

    if(_debug_mode) {

      std::cout << std::endl << Form("Found %zu granular showers...",_shower_index.size())<<std::endl;
      
      for(auto mother_iter = _shower_index.begin();
	  mother_iter != _shower_index.end();
	  ++mother_iter) {
	
	unsigned int mother_index = (*mother_iter).first;
	unsigned int shower_index = (*mother_iter).second;
	unsigned int mother_track = _track_id.at(mother_index);
	unsigned int ndaughters = _shower_daughters.at(shower_index).size();
	double x,y,z,t;
	GetTrackStartInfo(mother_index,x,y,z,t);

	std::cout<<Form("Shower track ID = %d, PDG=%d,  @(%g, %g, %g, %g) ...  %d daughters ",
			mother_track,
			_pdgcode.at(mother_index),
			x,y,z,t,
			ndaughters)<<std::endl;
	
	unsigned int first_daughter_index = (*(_shower_daughters.at(shower_index).begin()));
	unsigned int daughter_pdgcode = _pdgcode.at(first_daughter_index);
	
	GetTrackStartInfo(first_daughter_index,x,y,z,t);
	std::cout << Form("  Daughter %d starting @ (%g, %g, %g, %g) =>",daughter_pdgcode,x,y,z,t) << std::endl;
	
      }
      std::cout<<std::endl<<"End ConstructGranularShower ..."<<std::endl<<std::endl;
    }
  }

  void MCShowerRecoPart::CombineGranularShower()
  {
    // Let's get start & end track ID of each shower
    std::vector<unsigned int> track_start_v;
    track_start_v.reserve(_shower_daughters.size());
    std::vector<unsigned int> track_end_v;
    track_end_v.reserve(_shower_daughters.size());

    for(auto shower_iter = _shower_index.begin();
	shower_iter != _shower_index.end();
	++shower_iter) {
      
      unsigned int shower_index  = (*shower_iter).second;
      unsigned int primary_index = (*shower_iter).first;

      track_start_v.push_back(primary_index);
      track_end_v.push_back  ((*(_shower_daughters.at(shower_index).rbegin())));

      if(_debug_mode) std::cout << Form("Shower Track ID %d => %d", 
					_track_id.at(primary_index),
					track_end_v.at(track_end_v.size()-1))
				<<std::endl;
    }
    
    std::vector<std::set<unsigned int> > combination_v;
    combination_v.reserve(track_start_v.size());
    
    for(size_t i=0; i<track_start_v.size(); ++i) {

      std::set<unsigned int> combination;
      double tstart_i, xstart_i, ystart_i, zstart_i;
      GetTrackStartInfo(track_start_v.at(i),
			xstart_i,ystart_i,zstart_i,tstart_i);

      for(size_t j=0; j<track_start_v.size(); ++j) {
	
	// Check if shower "i" belongs to "j".
	// Skip if "i"=="j"
	if(i==j) continue;

	// Skip if the combination (i,j) is already inspected and found to be combined
	if(combination_v.size()>j && combination_v.at(j).find(i) != combination_v.at(j).end())
	  continue;

	// Skip if "j" start time is not inside "i" shower time, which is 1st daughter to last daughter time
	double tstart_j, xstart_j, ystart_j, zstart_j;
	double tend_j, xend_j, yend_j, zend_j;
	GetTrackStartInfo( (*(_shower_daughters.at(j).begin())),
			   xstart_j,ystart_j,zstart_j,tstart_j);
	GetTrackStartInfo( (*(_shower_daughters.at(j).rbegin())),
			   xend_j,yend_j,zend_j,tend_j);
	
	if( tstart_j > (tstart_i + _comb_time_cut) || (tend_j + _comb_time_cut) < tstart_i )
	  continue;

	if(_debug_mode) std::cout<<Form("Shower %zu starts within %zu! Inspecting for merging...",i,j)<<std::endl;

	bool combine = false;
	double dT = 0;
	double dL2 = 0;
	if( tstart_i < tstart_j ){
	// Case 1: shower "i" is just before shower "j" ... use 1st daughter
	  dT = tstart_i - tstart_j;
	  double x,y,z,t;
	  auto const daughter_iter = (_shower_daughters.at(j).begin());
	  GetTrackStartInfo((*daughter_iter), x, y, z, t);
	  dL2 = (pow(xstart_i - x,2) + pow(ystart_i - y,2) + pow(zstart_i - z,2));
	  combine = dL2 < _comb_dist2_cut;
	}else if( tend_j < tstart_i){
	// Case 2: shower "i" is just after shower "j"
	  dT = tstart_i - tend_j;
	  double x,y,z,t;
	  auto const daughter_iter = (_shower_daughters.at(j).rbegin());
	  GetTrackStartInfo((*daughter_iter), x, y, z, t);
	  dL2 = (pow(xstart_i - x,2) + pow(ystart_i - y,2) + pow(zstart_i - z,2));
	  combine = dL2 < _comb_dist2_cut;
	}else{
	// Case 3: shower "i" is within shower "j"

	// Loop over "j"'s daughters and see if i's mother is generated nearby
	  double x,y,z,t;
	  for(auto daughter_iter = _shower_daughters.at(j).begin();
	      daughter_iter != _shower_daughters.at(j).end();
	      ++daughter_iter) {

	    double       step_time  = _start_vtx.at((*daughter_iter)).at(3);
	    unsigned int step_index = (*daughter_iter);
	    
	    if( (step_time < tstart_i) && (tstart_i - step_time) > _comb_time_cut ) continue;
	    if( (step_time > tstart_i) && (step_time - tstart_i) > _comb_time_cut ) break;

	    if(_debug_mode) std::cout<<Form(" Inspecting @ %d...",_track_id.at(step_index))<<std::endl;
	    
	    GetTrackStartInfo(step_index,x,y,z,t);
	    double dist2 = (pow(xstart_i - x,2) + pow(ystart_i - y,2) + pow(zstart_i - z,2));
	    combine = dist2 < _comb_dist2_cut;
	    if(combine) { dL2 = dist2; dT = tstart_i - step_time; break; }
	  }

	}

	if(combine) { 
	 
	  if(_debug_mode) std::cout<<Form("    Found merging point! dT = %g, dX^2 = %g", dT,dL2)<<std::endl;

	  combination.insert(j);

	  if(j < i) combination_v.at(j).insert(i);
	}
      }
      combination.insert(i);
      combination_v.push_back(combination);
    }

    // Now let's find all possible pairs.
    // Loop over each vector content of combinations
    for(size_t i=0; i<combination_v.size(); ++i) {

      // Keep looking for more showers to be merged until
      // we find all. 
      while(combination_v.at(i).size() < combination_v.size()) {
	size_t my_size = combination_v.at(i).size();

	// Loop over other sets of shower index
	for(size_t j=0; j<combination_v.size(); ++j) {
	  // Don't bother the same shower index.
	  if(i==j) continue;

	  // Loop over MY shower indexes and see if any of them
	  // is in THIS shower indexes. If found so, break and
	  // merge.
	  bool do_merge = false;
	  for(auto merged_index : combination_v.at(i))

	    if(combination_v.at(j).find(merged_index) != combination_v.at(j).end())
	      {do_merge = true; break;}

	  if(do_merge)

	    for(auto index : combination_v.at(j))
	      
	      combination_v.at(i).insert(index);

	}
	
	if(my_size == combination_v.at(i).size())
	  break;
      }
    }

    if(_debug_mode){
      std::cout<<std::endl<<"Found combinations..."<<std::endl;
      for(size_t i=0; i<combination_v.size(); ++i) {
	std::cout<<Form("    Index %zu : ",i);
	for(auto comb : combination_v.at(i))

	  std::cout<<comb <<" ";
	std::cout<<std::endl;
      }
    }

    // Now make a combined shower index
    std::vector<std::set<unsigned int> > sorted_combination_v;

    for(size_t i=0; i<combination_v.size(); ++i) {

      int sorted_index = -1;

      // First, attempt to find "i" and associated indexes in sorted_combination_v contents.
      for(size_t j=0; j<sorted_combination_v.size(); ++j) {

	if(sorted_combination_v.at(j).find(i) != sorted_combination_v.at(j).end()) {
	  sorted_index = j; 
	  break;
	}
	for(auto const index : combination_v.at(i)) {

	  if(sorted_combination_v.at(j).find(index) != sorted_combination_v.at(j).end()) {
	    sorted_index = j;
	    break;
	  }

	}
	if(sorted_index>0) break;
      }
      // If not found, create a new set and insert
      if(sorted_index<0) {
	std::set<unsigned int> sorted_combination(combination_v.at(i));
	sorted_combination.insert(i);
	sorted_combination_v.push_back(sorted_combination);
      }
      // Else combine sets
      else{

	sorted_combination_v[sorted_index].insert(i);
	for(auto const index : combination_v.at(i))

	  sorted_combination_v[sorted_index].insert(index);
      }
    }

    if(_debug_mode){
      std::cout<<std::endl<<"Found combinations..."<<std::endl;
      for(size_t i=0; i<sorted_combination_v.size(); ++i) {
	std::cout<<Form("    Index %zu : ",i);
	for(auto comb : sorted_combination_v.at(i))

	  std::cout<<comb <<" ";
	std::cout<<std::endl;
      }
    }

    std::map<unsigned int,unsigned int> combined_shower_index;
    std::vector<std::map<double,std::vector<unsigned int> > > combined_shower_daughters;
    for(size_t i=0; i<sorted_combination_v.size(); ++i) {

      int super_mother_index=-1;
      double super_mother_time=-1;
      std::map<double,std::vector<unsigned int> > daughters;
      for(auto const index : sorted_combination_v.at(i)) {

	// Combine daughters
	for(auto const daughter_index : _shower_daughters.at(index)) {
	  
	  double daughter_time = _start_vtx.at(daughter_index).at(3);
	  if(daughters.find(daughter_time)==daughters.end())
	    daughters.insert(std::pair<double,std::vector<unsigned int> >(daughter_time,
									  std::vector<unsigned int>(1,
												    daughter_index)
									  )
			     );
	  else
	    daughters[daughter_time].push_back(daughter_index);
	  
	}

	// Find super mother
	int mother_index = track_start_v.at(index);
	double mother_time = _start_vtx.at(mother_index).at(3);
	if(super_mother_index<0 || mother_time < super_mother_time) {
	  super_mother_index=mother_index;
	  super_mother_time=mother_time;
	}else if(mother_time == super_mother_time && mother_index < super_mother_index) {
	  super_mother_index = mother_index;
	  super_mother_time=mother_time;
	}
      }
      // Add mothers that is not super mother, to a daughter
      for(auto const index : sorted_combination_v.at(i)) {
	int mother_index = track_start_v.at(index);
	double mother_time = _start_vtx.at(mother_index).at(3);
	if(mother_index!=super_mother_index) {
	  if(daughters.find(mother_time)==daughters.end())
	    daughters.insert(std::pair<double,std::vector<unsigned int> >(mother_time,
									  std::vector<unsigned int>(1,
												    mother_index)
									  )
			     );
	  else
	    daughters[mother_time].push_back(mother_index);
	}
	else
	  combined_shower_index.insert(std::pair<unsigned int,unsigned int>(super_mother_index,combined_shower_daughters.size()));
      }
      combined_shower_daughters.push_back(daughters);
    }

    // Re-set showers
    _shower_index = combined_shower_index; 

    // Update shower id
    _shower_daughters.clear();
    _shower_daughters.reserve(combined_shower_daughters.size());
    for(size_t i=0; i<_shower_id.size(); ++i)
      _shower_id[i]=-1;

    for(auto const daughters_map : combined_shower_daughters) {
      
      unsigned int shower_index = _shower_daughters.size();
      std::vector<unsigned int> daughters_v;
      daughters_v.reserve(daughters_map.size());
      for(auto daughters_map_iter = daughters_map.begin();
	  daughters_map_iter != daughters_map.end();
	  ++daughters_map_iter) {

	for(auto const index : (*daughters_map_iter).second) {
	  _shower_id[index] = shower_index;
	  daughters_v.push_back(index);
	}
      }
      _shower_daughters.push_back(daughters_v);
    }

    for(auto mother_iter = _shower_index.begin();
	mother_iter != _shower_index.end();
	++mother_iter) {

      _shower_id[(*mother_iter).first] = (*mother_iter).second;

      if(_debug_mode) 
	std::cout
	  <<Form("  Combined shower %d ... mother track=%d, pdg=%d with %zu daughters...",
		 (*mother_iter).second,
		 _track_id.at((*mother_iter).first),
		 _pdgcode.at((*mother_iter).first),
		 _shower_daughters.at((*mother_iter).second).size())
	  << std::endl;
    }    

    if(_debug_mode) {
      std::vector<unsigned int> part_count_v(_shower_index.size(),0);
      unsigned int undefined=0;
      std::cout<<std::endl;
      for(size_t i=0; i<_shower_id.size(); ++i) {
	if(_shower_id.at(i)<0) {
	  undefined++;
	  std::cout << Form("  Track %d (PDG=%d) @ (%g, %g, %g, %g) with %zu daughters does not belong to a shower!",
			    _track_id[i],
			    _pdgcode.at(i),
			    _start_vtx.at(i).at(0), _start_vtx.at(i).at(1), _start_vtx.at(i).at(2), _start_vtx.at(i).at(3),
			    _daughters.at(i).size()) 
		    << std::endl;
	}else if((size_t)(_shower_id.at(i))>=part_count_v.size())
	  
	  throw cet::exception(__FUNCTION__) << Form("Track %d PDG %d has ill-defined shower index %d!",
						     _track_id.at(i),_pdgcode.at(i),_shower_id.at(i));
	else
	  part_count_v[_shower_id.at(i)]++;
      }
      std::cout << Form("  %d tracks do not belong to shower...",undefined) << std::endl;
      
      for(size_t i=0; i<part_count_v.size(); ++i)
	
	std::cout << Form("  %d tracks belong to shower %zu...", part_count_v.at(i),i)<<std::endl;
      
    }
    
  } // namespace opdet

}
#endif
