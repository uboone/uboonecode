

#include <iostream>
#include <vector>
#include <fstream>
#include "../../GetPermutations.h"


int main(int const argc, char const * argv[]) {

  int const jobs = 10;

  double fstart_prox = 4;
  double fstart_prox_min = 2;
  double fstart_prox_max = 10;
  double fstart_prox_inc = 4;
  double fshower_prox = -1;
  double fshower_prox_min = 5;
  double fshower_prox_max = 40;
  double fshower_prox_inc = 7;
  double fmax_bp_dist = -1;
  double fmax_bp_dist_min = 5;
  double fmax_bp_dist_max = 50;
  double fmax_bp_dist_inc = 9;
  double fcpoa_vert_prox = 13;
  double fcpoa_vert_prox_min = 1; 
  double fcpoa_vert_prox_max = 13;
  double fcpoa_vert_prox_inc = 12;
  double fcpoa_trackend_prox = 10;
  double fcpoa_trackend_prox_min = 2;
  double fcpoa_trackend_prox_max = 10;
  double fcpoa_trackend_prox_inc = 4;  

  std::vector<std::string> const set_str = {
    "start_prox",
    "shower_prox",
    "max_bp_dist",
    "cpoa_vert_prox",
    "cpoa_trackend_prox"
  };

  std::vector<double> const parameter_set = {fstart_prox, fshower_prox, fmax_bp_dist, fcpoa_vert_prox, fcpoa_trackend_prox};
  std::vector<double> const parameter_min = {fstart_prox_min, fshower_prox_min, fmax_bp_dist_min, fcpoa_vert_prox_min, fcpoa_trackend_prox_min};
  std::vector<double> const parameter_max = {fstart_prox_max, fshower_prox_max, fmax_bp_dist_max, fcpoa_vert_prox_max, fcpoa_trackend_prox_max};
  std::vector<double> const parameter_inc = {fstart_prox_inc, fshower_prox_inc, fmax_bp_dist_inc, fcpoa_vert_prox_inc, fcpoa_trackend_prox_inc};
  std::vector<double> const parameter_inc_size = {(fstart_prox_max - fstart_prox_min) / fstart_prox_inc,
						  (fshower_prox_max - fshower_prox_min) / fshower_prox_inc,
						  (fmax_bp_dist_max - fmax_bp_dist_min) / fmax_bp_dist_inc,
						  (fcpoa_vert_prox_max - fcpoa_vert_prox_min) / fcpoa_vert_prox_inc,
						  (fcpoa_trackend_prox_max - fcpoa_trackend_prox_min) / fcpoa_trackend_prox_inc};

  std::vector<std::vector<double>> permutation_v;

  GetPermutations(parameter_set,
		  parameter_min,
		  parameter_max,
		  parameter_inc_size,
		  permutation_v);
  CheckPermutations(parameter_set, parameter_inc, permutation_v);

  size_t const job_inc = permutation_v.size() / jobs;
  size_t job_mod = permutation_v.size() % jobs;
  int index = -1;

  std::vector<std::pair<size_t, size_t>> index_v;
  for(size_t i = 0; i < jobs; ++i) {
    ++index;
    size_t min = index;
    index += job_inc - 1;
    if(job_mod > 0) {
      ++index;
      --job_mod;
    }
    index_v.push_back({min, index});
  }

  for(size_t file_i = 0; file_i < jobs; ++file_i) {
  
    std::ofstream ofile("output_fcls/file"+std::to_string(file_i)+".fcl");

    ofile << "#include \"services_microboone.fcl\"\n#include \"messageservice.fcl\"\n#include \"time_memory_tracker_microboone.fcl\"\n\nprocess_name : RunVertexBuilder\n\nservices:\n{\n\n\tTFileService: {fileName: \"RunVertexBuilder.root\"}\n\tTimeTracker:  @local::microboone_time_tracker\n\tMemoryTracker: @local::microboone_memory_tracker\n\t@table::microboone_simulation_services\n\n\n}#end services\n\nsource:\n{\n\n\tmodule_type: RootInput\n\tmaxEvents: -1\n\tfirstRun: 1\n\tfirstEvent: 1\n\n}\n\nphysics:\n{\n\tanalyzers:\n\t{\n\t\tRunVB:\n\t\t{\n\t\t\tmodule_type: \"RunVertexBuilder\"\n\t\t\ttrack_producer: \"pandoraNu\"\n\t\t\tshower_producer: \"pandoraNu\"\n\t\t\thit_producer: \"pandoraCosmicHitRemoval\"\n\t\t\trmcmassociation_producer: \"crHitRemovalTruthMatch\"\n\n";

    //ofile << "\t\t\ttrack_only: \"true\"\n\n";

    for(size_t i = 0; i < parameter_set.size(); ++i) {
      double const d = parameter_set.at(i);
      std::string const par = "\t\t\t" + set_str.at(i);
      if(d != -1) ofile << par << ": \"" << d << "\"\n\n";
      else {
	ofile << par << "_min: \"" << parameter_min.at(i) << "\"\n"
	      << par << "_max: \"" << parameter_max.at(i) << "\"\n"
	      << par << "_inc: \"" << parameter_inc.at(i) << "\"\n\n";
      }
    }

    ofile << "\t\t\tmin_index: \"" << index_v.at(file_i).first << "\"\n"
	  << "\t\t\tmax_index: \"" << index_v.at(file_i).second << "\"\n\n"
	  << "\n\t\t}\n\t}\n\tanalysis: [RunVB]\n\tend_paths: [analysis]\n}\n\nservices.DetectorPropertiesService.NumberTimeSamples:        6400\nservices.DetectorPropertiesService.ReadOutWindowSize:        6400\nservices.DetectorClocksService.InheritClockConfig:           false\nservices.DetectorClocksService.TriggerOffsetTPC:             -0.400e3\nservices.DetectorClocksService.TrigModuleName:               \"daq\"\nservices.DetectorClocksService.InheritClockConfig: false\n";

    ofile.close();

  }

  return 0;

}
