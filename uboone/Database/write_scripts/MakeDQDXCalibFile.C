void MakeDQDXCalibFile(std::string file_name, std::string out_dir, std::string out_tag) {

  
  TFile* f = new TFile(file_name.c_str());
  std::vector<TH1D*> x_corrections;
  std::vector<TH2D*> yz_corrections;

  std::string base_name_x = "correction_x_plane";
  std::string base_name_yz = "correction_yz_plane";
  for (unsigned int plane=0; plane!=3; ++plane) {
    std::string full_name_x = base_name_x + std::to_string(plane);
    x_corrections.push_back((TH1D*)f->Get(full_name_x.c_str()));
    
    std::string full_name_yz = base_name_yz + std::to_string(plane);
    yz_corrections.push_back((TH2D*)f->Get(full_name_yz.c_str()));
  }

  TCanvas* c = new TCanvas("c","c");
  yz_corrections[2]->Draw("COLZ");
  
  for (unsigned int plane=0; plane!=3; ++plane) {
    std::string out_name_x = out_dir+"/dqdx_xcorrections_plane"+std::to_string(plane)+"_run_NA_"+out_tag+".dat";
    
    std::ofstream out_file_x;
    out_file_x.open(out_name_x);
    
    out_file_x << "# time 1" << std::endl;
    out_file_x << "# channel x_bin_low_edge x_bin_high_edge correction correction_err"<< std::endl;
    
    unsigned int global_bin_counter_x=0;
    for (Int_t binX=1; binX <= x_corrections[plane]->GetNbinsX(); ++binX) {
      out_file_x <<      binX-1 
               <<" "<< x_corrections[plane]->GetBinLowEdge(binX) 
	       <<" "<< x_corrections[plane]->GetBinLowEdge(binX+1)
	       <<" "<< x_corrections[plane]->GetBinContent(binX)
	       <<" "<< x_corrections[plane]->GetBinContent(binX)/20.0
	       <<std::endl;
      global_bin_counter_x = binX;
    }
    
    //add dummy bins
    while (global_bin_counter_x <= 130) {
      out_file_x <<   global_bin_counter_x
                 <<" "<< 1.0e6
		 <<" "<< 1.1e6
		 <<" "<< 1.0
		 <<" "<< 1.0
		 <<std::endl;
      global_bin_counter_x++;
    }
    out_file_x.close();
    
    std::string out_name_yz = out_dir+"/dqdx_yzcorrections_plane"+std::to_string(plane)+"_run_NA_"+out_tag+".dat";
    
    std::ofstream out_file_yz;
    out_file_yz.open(out_name_yz);
    
    out_file_yz << "# time 1" << std::endl;
    out_file_yz << "# channel y_bin_low_edge y_bin_high_edge z_bin_low_edge z_bin_high_edge correction correction_err"<< std::endl;
    
    unsigned int global_bin_counter_yz=0;
    for (Int_t binY=1; binY <= yz_corrections[plane]->GetNbinsY(); ++binY) {
      for (Int_t binZ=1; binZ <= yz_corrections[plane]->GetNbinsX(); ++binZ) {
        out_file_yz << global_bin_counter_yz
	            << " " << yz_corrections[plane]->GetYaxis()->GetBinLowEdge(binY) 
		    << " " << yz_corrections[plane]->GetYaxis()->GetBinLowEdge(binY+1) 
                    << " " << yz_corrections[plane]->GetXaxis()->GetBinLowEdge(binZ) 
		    << " " << yz_corrections[plane]->GetXaxis()->GetBinLowEdge(binZ+1) 
		    << " " << yz_corrections[plane]->GetBinContent(binZ, binY)
	            << " " << yz_corrections[plane]->GetBinContent(binZ, binY)/20.0
		    <<std::endl;
        global_bin_counter_yz++;
      }
    }
    
    //add dummy bins
     while (global_bin_counter_yz <= 10000) {
      out_file_yz <<   global_bin_counter_yz
                 <<" "<< 1.0e6
		 <<" "<< 1.1e6
		 <<" "<< 1.0e6
		 <<" "<< 1.1e6
		 <<" "<< 1.0
		 <<" "<< 1.0
		 <<std::endl;
      global_bin_counter_yz++;
    }
    out_file_yz.close();
  }//end loop over planes    
}
