


void GetPermutations(std::vector<double> const & parameter_set,
		     std::vector<double> const & parameter_min,
		     std::vector<double> const & parameter_max,
		     std::vector<double> const & parameter_inc_size,
		     std::vector<std::vector<double>> & permutation_v) {

  bool continue_bool;
  std::vector<double> parameters = parameter_min;
  for(size_t i = 0; i < parameter_set.size(); ++i) {
    if(parameter_set.at(i) != -1) parameters.at(i) = parameter_set.at(i);     
  }
  permutation_v.push_back(parameters);
  do {
    continue_bool = false;
    for(int i = parameters.size() - 1; i >= 0; --i) {
      if(parameter_set.at(i) != -1) continue;
      if(parameters.at(i) < parameter_max.at(i)) {
	parameters.at(i) += parameter_inc_size.at(i);
	for(size_t j = parameters.size() - 1; int(j) > i; --j) {
	  if(parameter_set.at(j) == -1) parameters.at(j) = parameter_min.at(j);
	  else parameters.at(j) = parameter_set.at(j);
	}
	continue_bool = true;
	break;
      }
    } 
    if(continue_bool) permutation_v.push_back(parameters);
  } while(continue_bool);

}



void CheckPermutations(std::vector<double> const & parameter_set,
		       std::vector<double> const & parameter_inc,
		       std::vector<std::vector<double>> const & permutation_v) {

  int permutations = 1;
  for(size_t i = 0; i < parameter_set.size(); ++i) {
    if(parameter_set.at(i) != -1) continue;
    permutations *= (parameter_inc.at(i) + 1);
  }

  for(std::vector<double> const & v : permutation_v) {
    for(double const d : v) std::cout << d << " ";
    std::cout << "\n";
  }

  std::cout << "size: " << permutation_v.size() << " permutations: " << permutations << "\n";

}

