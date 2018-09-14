#ifndef UBOONEDQDXCORRECTIONPROVIDER_CXX
#define UBOONEDQDXCORRECTIONPROVIDER_CXX

#include "UboonedqdxCorrectionProvider.h"

// art/LArSoft libraries
#include "cetlib_except/exception.h"


#include <fstream>

namespace lariov {

  //constructor      
  UboonedqdxCorrectionProvider::UboonedqdxCorrectionProvider(fhicl::ParameterSet const& p) :
    DatabaseRetrievalAlg(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg")) {	
    
    this->Reconfigure(p);
  }
  
  void UboonedqdxCorrectionProvider::Reconfigure(fhicl::ParameterSet const& p) {
    
    this->DatabaseRetrievalAlg::Reconfigure(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg"));
    fData.Clear();
    IOVTimeStamp tmp = IOVTimeStamp::MaxTimeStamp();
    tmp.SetStamp(tmp.Stamp()-1, tmp.SubStamp());
    fData.SetIoV(tmp, IOVTimeStamp::MaxTimeStamp());
    
    fCoord_min.clear();
    fCoord_max.clear();
    fNbins.clear();
    fBinEdgeNames.clear();
    fIsFixedBinSize = true;
    fCoordToBinMap.clear();

    bool UseDB      = p.get<bool>("UseDB", false);
    bool UseFile    = p.get<bool>("UseFile", false);
    std::string fileName = p.get<std::string>("FileName", "");

    //priority:  (1) use db, (2) use table, (3) use defaults
    //If none are specified, use defaults
    if ( UseDB )      fDataSource = DataSource::Database;
    else if (UseFile) fDataSource = DataSource::File;
    else              fDataSource = DataSource::Default;

    if (fDataSource == DataSource::Default) {
      float default_correction     = p.get<float>("DefaultCorrection");
      float default_correction_err = p.get<float>("DefaultCorrectionErr");

      DqdxCorrection defaultCalib(0);
      defaultCalib.SetCorrection(default_correction);
      defaultCalib.SetCorrectionErr(default_correction_err);
      fData.AddOrReplaceRow(defaultCalib);
      
    }
    else if (fDataSource == DataSource::File) {
      cet::search_path sp("FW_SEARCH_PATH");
      std::string abs_fp = sp.find_file(fileName);
      std::cout << "Using dqdx corrections from local file: "<<abs_fp<<"\n";
      std::ifstream file(abs_fp);
      if (!file) {
        throw cet::exception("UboonedqdxCorrectionProvider")
          << "File "<<abs_fp<<" is not found.";
      }
      
      std::string line;
      int line_number = -1;
      unsigned int NDimensions = 0;
      int prev_bin = -1;
      std::vector<std::vector<float> > bin_list;
      std::vector<std::vector<float> > prev_bin_list;
      std::vector<float> prev_low_edge;
      std::vector<float> prev_high_edge;
      while (std::getline(file, line)) {
        line_number++;
	
	
	//-------------------
	//initialization
	//-------------------
	if (line_number==0) {
	  NDimensions = std::count(line.begin(), line.end(), ',') - 2;
	  bin_list.resize(NDimensions);
	  prev_bin_list.resize(NDimensions);
	  prev_low_edge.resize(NDimensions);
	  prev_high_edge.resize(NDimensions);
	}
      
      
        //-------------------
        //get bin number, make sure it is reasonable
	//-------------------
        size_t current_comma = line.find(',');
        DBChannelID_t bin = (DBChannelID_t)std::stoi(line.substr(0, current_comma)); 
	if ( (int)bin != prev_bin+1 || (int)bin != line_number) {
	  throw cet::exception("UboonedqdxCorrectionProvider")
	    << "Bins in file "<<abs_fp<<" are not in order of increasing bin number or do not start at zero";
	}
	prev_bin = bin;
	
	
	//-------------------
	//get bin high and low edges
	//-------------------
	std::vector<float> low_edge;
	std::vector<float> high_edge;
	for (unsigned int n_dim = 0; n_dim != NDimensions; ++n_dim) {
      
          low_edge.push_back( std::stof(line.substr(current_comma+1, line.find(',',current_comma+1))) );        
          current_comma = line.find(',',current_comma+1);
          
	  high_edge.push_back( std::stof(line.substr(current_comma+1, line.find(',',current_comma+1))) );
	  current_comma = line.find(',',current_comma+1);
	  
	  //sanity check
	  if (high_edge[n_dim] <= low_edge[n_dim]) {
	    throw cet::exception("UboonedqdxCorrectionProvider")
	      << "Some bins in file "<<abs_fp<<" have negative size";
	  }
	}
	
	
	//-------------------
	//make sure bins are ordered correctly in file by checking that this current bin low_edge is greater than the previous bin lowedge
	//-------------------
	/*if (line_number != 0 && !UboonedqdxCorrectionProvider::CompareFunction(prev_low_edge, low_edge) ) {
	  throw cet::exception("UboonedqdxCorrectionProvider")
	    << "Bins in file "<<abs_fp<<" are not ordered correctly!";
        }*/
	
	
	//-------------------
	//maintain list of bin edges in each dimension
	//-------------------
	for (unsigned int n_dim = 0; n_dim != NDimensions; ++n_dim) {
	  //normally we add to the list of bin edges if the bin edge increases
	  if ( line_number==0 || low_edge[n_dim] > bin_list[n_dim].back() ) {
	    bin_list[n_dim].push_back(low_edge[n_dim]);
	  }
	  //unless we encounter a decrease in this dimension, which signifies that we looped back to the beginning
	  else if (low_edge[n_dim] < bin_list[n_dim].back()) {
	  
	    //finish current bin list by adding previous high edge
	    bin_list[n_dim].push_back(prev_high_edge[n_dim]);
	    
	    //if we just finished the first time through the bins in this dimension, then initialize the prev_bin_list 
	    if (prev_bin_list[n_dim].empty()) {
	      prev_bin_list[n_dim] = bin_list[n_dim];
	      
	    }
	    //else check for consistency with the previous bin list
	    else if (bin_list[n_dim] != prev_bin_list[n_dim]) {
	      fIsFixedBinSize = false;
	      throw cet::exception("UboonedqdxCorrectionProvider")
		<< "Bins in file "<<abs_fp<<" are not consistent!";
            }
	    
	    //Now that we finished checking the bin list, clear it and start filling again
	    bin_list[n_dim].clear();
	    bin_list[n_dim].push_back(low_edge[n_dim]);
	  }
	  	  
	} //end loop for maintaining the list of bin edges in each dimension
	
	
	//-------------------	
	//update prev_high_edge and prev_low_edge
	//-------------------
	prev_low_edge = low_edge;
	prev_high_edge = high_edge; 
	
	
	//-------------------
	//get the coorection and its error, make sure the number of columns hasn't changed
	//-------------------
        float correction = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));	
	current_comma = line.find(',',current_comma+1);
        float correction_err = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	if ( current_comma == std::string::npos || current_comma != line.find_last_of(',') ) {
	  throw cet::exception("UboonedqdxCorrectionProvider")
          << "Number of columns in calibration file do not match n dimensions";
        }


        //-------------------
        //Add calibration data to fData
	//-------------------
	DqdxCorrection dp(0);
        dp.SetChannel(bin);
	dp.SetCorrection(correction);
        dp.SetCorrectionErr(correction_err);       
        fData.AddOrReplaceRow(dp);
	

        //-------------------
	//fill map of bin number to coordinate low edge
	//-------------------
	fCoordToBinMap.push_back(low_edge);
      }//end loop over file lines
      
      //Set fCoord_min, fCoord_max, fIsFixedSize, and fNbins.
      for (unsigned int n_dim=0; n_dim!=NDimensions; ++n_dim) {
        //First add the last high edge to each of the bin lists 
        bin_list[n_dim].push_back(prev_high_edge[n_dim]);
	
	//Do the last consistency check
	if (n_dim>0 && bin_list[n_dim] != prev_bin_list[n_dim]) {
	  fIsFixedBinSize = false;
	  throw cet::exception("UboonedqdxCorrectionProvider")
	    << "Bins in file "<<abs_fp<<" are not consistent!";
        }
	
	//Do the fixed size check
	if (fIsFixedBinSize) {
	  float first_bin_size = bin_list[n_dim][1] - bin_list[n_dim][0];
	  for (std::vector<float>::const_iterator itV = bin_list[n_dim].begin()+2; itV != bin_list[n_dim].end(); ++itV) {
	    if ( fabs((*itV - (*itV-1)) - first_bin_size)/first_bin_size > 1.0e-6) {
	      fIsFixedBinSize = false;
	      break;
	    }
	  }
	}
      
        //Set these variables using bin_list
        fCoord_min.push_back(bin_list[n_dim].front());
	fCoord_max.push_back(bin_list[n_dim].back());	
	fNbins[n_dim] = bin_list[n_dim].size() - 1;
      }

    } //end if datasource is a file
    else {    
      fBinEdgeNames = p.get<std::vector<std::string> >("BinEdgeNames");
    }
  }



  bool UboonedqdxCorrectionProvider::Update(DBTimeStamp_t ts) {
    
    if (fDataSource != DataSource::Database) return false;
      
    if (!this->UpdateFolder(ts)) return false;

    //DBFolder was updated, so now update the Snapshot
    fData.Clear();
    fData.SetIoV(this->Begin(), this->End());

    std::vector<DBChannelID_t> bins;
    fFolder->GetChannelList(bins);
    
    //prep
    fCoord_min.resize(fBinEdgeNames.size());
    fCoord_max.resize(fBinEdgeNames.size());
    fNbins.resize(fBinEdgeNames.size(),0);
    fIsFixedBinSize = true;
    fCoordToBinMap.clear();
    
    
    std::vector<float> prev_low_edge(fBinEdgeNames.size());
    std::vector<float> prev_high_edge(fBinEdgeNames.size());
    std::vector<std::vector<float> > bin_list(fBinEdgeNames.size());
    std::vector<std::vector<float> > prev_bin_list(fBinEdgeNames.size());  
    for (auto it_bin = bins.begin(); it_bin != bins.end(); ++it_bin) {

      //-------------------
      //make sure bin numbers are increasing by 1
      //-------------------
      if ( (it_bin==bins.begin() && *it_bin != 0) ||
           (it_bin!=bins.begin() && *it_bin != (*(it_bin-1)) + 1) )  {
	throw cet::exception("UboonedqdxCorrectionProvider")
	  << "Bins in TPC dq/dx calib database are not in order of increasing bin number or do not start at zero";
      }
      
      
      //-------------------
      //Get bin edges
      //-------------------
      std::vector<float> low_edge(fBinEdgeNames.size());
      std::vector<float> high_edge(fBinEdgeNames.size());
      double correction, correction_err;
      bool dummy_bin=true;
      for (unsigned int n_dim = 0; n_dim != fBinEdgeNames.size(); ++n_dim) {
      
        std::string low_edge_str = fBinEdgeNames[n_dim]+"_low_edge";
	std::string high_edge_str = fBinEdgeNames[n_dim]+"_high_edge";
	double le, he;
        fFolder->GetNamedChannelData(*it_bin, low_edge_str,   le);
        fFolder->GetNamedChannelData(*it_bin, high_edge_str,  he); 
	low_edge[n_dim] = (float)le;
	high_edge[n_dim] = (float)he;
	
	dummy_bin = dummy_bin && le > 1.0e5;
	
	//sanity check
	if (high_edge[n_dim] <= low_edge[n_dim]) {
	  throw cet::exception("UboonedqdxCorrectionProvider")
	    << "Some bins in TPC dq/dx calib database have negative size";
	}
      }
      
      //-------------------
      //break out of loop over bins if we've reached the last bin
      //i.e. the bin edges have dummy values
      //-------------------
      if (dummy_bin) break;
      
      
      //-------------------
      //make sure bins are ordered correctly in the database by checking that this current bin low_edge is greater than the previous bin lowedge
      //-------------------
      /*if (it_bin != bins.begin() && !UboonedqdxCorrectionProvider::CompareFunction(prev_low_edge, low_edge) ) {
	throw cet::exception("UboonedqdxCorrectionProvider")
	  << "Bins in TPC dq/dx calib database are not ordered correctly! "<<prev_low_edge.front()<<","<<prev_low_edge.back()
	  <<" vs. "<<low_edge.front()<<","<<low_edge.back();
      }*/
	
	
      //-------------------
      //maintain list of bin edges in each dimension
      //-------------------
      for (unsigned int n_dim = 0; n_dim != fBinEdgeNames.size(); ++n_dim) {
	//normally we add to the list of bin edges if the bin edge increases
	if ( it_bin==bins.begin() || low_edge[n_dim] > bin_list[n_dim].back() ) {
	  bin_list[n_dim].push_back(low_edge[n_dim]);
	}
	//unless we encounter a decrease in this dimension, which signifies that we looped back to the beginning
	else if (low_edge[n_dim] < bin_list[n_dim].back()) {

	  //finish current bin list by adding previous high edge
	  bin_list[n_dim].push_back(prev_high_edge[n_dim]);

	  //if we just finished the first time through the bins in this dimension, then initialize the prev_bin_list 
	  if (prev_bin_list[n_dim].empty()) {
	    prev_bin_list[n_dim] = bin_list[n_dim];
	  }
	  //else check for consistency with the previous bin list
	  else if (bin_list[n_dim] != prev_bin_list[n_dim]) {
	    fIsFixedBinSize = false;
	    throw cet::exception("UboonedqdxCorrectionProvider")
	      << "Bins in TPC dq/dx calib database are not consistent sizes!";
          }

	  //Now that we finished checking the bin list, clear it and start filling again
	  bin_list[n_dim].clear();
	  bin_list[n_dim].push_back(low_edge[n_dim]);
	}

      } //end loop for maintaining the list of bin edges in each dimension
	
      
      //-------------------	
      //update prev_high_edge and prev_low_edge
      //-------------------
      prev_low_edge = low_edge;
      prev_high_edge = high_edge; 

      
      //-------------------
      //get the coorection and its error, add to local cache
      //-------------------
      fFolder->GetNamedChannelData(*it_bin, "correction",     correction);
      fFolder->GetNamedChannelData(*it_bin, "correction_err", correction_err); 
      
      DqdxCorrection pg(*it_bin);
      
      pg.SetCorrection( (float)correction );
      pg.SetCorrectionErr( (float)correction_err );

      fData.AddOrReplaceRow(pg);

      //fill map of bin number to coordinate low edge
      fCoordToBinMap.push_back(low_edge);
    }//end loop over bins

    //Set fCoord_min, fCoord_max, fIsFixedSize, and fNbins.
    for (unsigned int n_dim=0; n_dim!=fBinEdgeNames.size(); ++n_dim) {
      //First add the last high edge to each of the bin lists 
      bin_list[n_dim].push_back(prev_high_edge[n_dim]);

      //Do the last consistency check
      if (n_dim>0 && bin_list[n_dim] != prev_bin_list[n_dim]) {
	fIsFixedBinSize = false;
	throw cet::exception("UboonedqdxCorrectionProvider")
	  << "Bins in TPC dq/dx calib database are not consistent sizes!";
      }

      //Do the fixed size check
      if (fIsFixedBinSize) {
	float first_bin_size = bin_list[n_dim][1] - bin_list[n_dim][0];
	for (std::vector<float>::const_iterator itV = bin_list[n_dim].begin()+2; itV != bin_list[n_dim].end(); ++itV) {
	  if ( fabs((*itV - *(itV-1)) - first_bin_size)/first_bin_size > 1.0e-6) {
	    fIsFixedBinSize = false;
	    std::cout<<"Bin size varied: "<<n_dim<<" "<<*itV<<" "<<*(itV-1)<<" "<<first_bin_size<<std::endl;
	    break;
	  }
	}
      }

      //Set these variables using bin_list
      fCoord_min[n_dim] = bin_list[n_dim].front();
      fCoord_max[n_dim] = bin_list[n_dim].back();
      //std::cout<<"STATEMENT: Set min and max bin coords for this dimension to: "<<fCoord_min[n_dim]<<" and "<<fCoord_max[n_dim]<<std::endl;	
      fNbins[n_dim] = bin_list[n_dim].size() - 1;
    }
    
    return true;
  }
  
  
  bool UboonedqdxCorrectionProvider::CompareFunction(const std::vector<float>& coord, const std::vector<float>& bin_lower_bounds) {
    if (bin_lower_bounds.empty()) return false;
    if (coord.empty()) return true;
    
    std::vector<float>::const_iterator it_coord = coord.cbegin();
    std::vector<float>::const_iterator it_bounds = bin_lower_bounds.cbegin();
    
    std::vector<float>::const_iterator end_coord = coord.cend();
    std::vector<float>::const_iterator end_bounds = bin_lower_bounds.cend();
    
    while (it_coord!=end_coord && it_bounds!=end_bounds) {
      if (*it_coord < *it_bounds) return false;
      ++it_coord; ++it_bounds;
    }
    return false;
  }
    
  
  int UboonedqdxCorrectionProvider::GetBinNumber(const std::vector<float>& input_coord) const {
    
    std::vector<float> coord = input_coord;
    if (coord.size() != fCoord_min.size()) {
      throw cet::exception("UboonedqdxCorrectionProvider")
        << "number of input coordinate dimensions ("<<coord.size()<<") doesn't match initialization ("<<fCoord_min.size()<<")";
    }
    
    for (unsigned int n_dim=0; n_dim != coord.size(); ++n_dim) {
      if (coord[n_dim] < fCoord_min[n_dim]) {
        coord[n_dim] = fCoord_min[n_dim];
      }
      else if (coord[n_dim] >= fCoord_max[n_dim]) {
        coord[n_dim] = fCoord_max[n_dim]*(1.0 - 1.0e-5);
      }
    }
    
    //std::cout<<"STATEMENT: bin ranges: "<<coord[0]<<" - "<<fCoord_min[0]<<","<<fCoord_max[0]<<std::endl;
    
    if (fIsFixedBinSize) {
      int ret_val = 0;
      int dim_product = 1;
      for (int n_dim=coord.size()-1; n_dim >= 0; --n_dim) {
        ret_val += ((int)(fNbins[n_dim]*(coord[n_dim]-fCoord_min[n_dim])/(fCoord_max[n_dim]-fCoord_min[n_dim])))*dim_product;
	dim_product *= fNbins[n_dim];
      }
      /*std::cout<<"Chosen bin "<<ret_val<<std::endl;
      for (int n_dim=coord.size()-1; n_dim >= 0; --n_dim) {
        std::cout<<fNbins[n_dim]<<" "<<fCoord_min[n_dim]<<" "<<fCoord_max[n_dim]<<std::endl;
      }*/
      return ret_val;
    }
    
    //not fixed bins
    throw cet::exception("UboonedqdxCorrectionProvider")
      <<"GetBinNumber: Not fixed bins: doesn't work yet"<<std::endl;
      
    /*if (coord.size()==2) {
      int tmp_bin=0;
      for (auto iter = fCoordToBinMap.begin(); iter!= fCoordToBinMap.end(); ++iter, ++tmp_bin) {
        std::cout<<"Bin: "<<tmp_bin<<" coord "<<iter->at(0)<<","<<iter->at(1)<<std::endl;
      }
    }
    
    std::vector<std::vector<float> >::const_iterator it = std::upper_bound(fCoordToBinMap.begin(), fCoordToBinMap.end(), coord, UboonedqdxCorrectionProvider::CompareFunction);
    if (coord.size()==2) std::cout<<"Chosen bin coord: "<<it->at(0)<<","<<it->at(1)<<std::endl;
    std::cout<<"Chosen bin number: "<<std::distance(fCoordToBinMap.cbegin(),it)-1<<std::endl;
    std::cout<<"Compares: "<<this->CompareFunction(coord, *it)<<std::endl;
    if ( it==fCoordToBinMap.end() ) {
      throw cet::exception("UboonedqdxCorrectionProvider")
        << "coordinate is out of search bounds!";
    }
    return std::distance(fCoordToBinMap.cbegin(),it)-1;*/
  }
  
  
  const DqdxCorrection& UboonedqdxCorrectionProvider::DqdxCorrectionObject(const std::vector<float>& coord) const { 
    return fData.GetRow( this->GetBinNumber(coord) );
  }
 
      
  float UboonedqdxCorrectionProvider::Correction(const std::vector<float>& coord) const {
    return this->DqdxCorrectionObject(coord).Correction();
  }
  
  
  float UboonedqdxCorrectionProvider::CorrectionErr(const std::vector<float>& coord) const {
    return this->DqdxCorrectionObject(coord).CorrectionErr();
  }
  
}//end namespace lariov
	
#endif
        
