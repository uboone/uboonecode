#ifndef _WEIGHTCALC_H_
#define _WEIGHTCALC_H_
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>
#include <map>

//weight calc base
namespace evwgh {
  typedef std::map<std::string, std::vector<double> > WeightMap_t;

  class WeightCalc
  {    
  public:
    virtual void                Configure(fhicl::ParameterSet const& pset) = 0;
    virtual std::vector<std::vector<double> > GetWeight(art::Event & e) = 0;  
    void                        SetName(std::string name) {fName=name;};
    std::string                 GetName() {return fName;}; 
    
  private:
    std::string fName;
  };
}

#endif // _WEIGHTCALC_H_
