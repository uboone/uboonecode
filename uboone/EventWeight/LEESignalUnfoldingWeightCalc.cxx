#include "WeightCalcCreator.h"
#include "WeightCalc.h"

#include <iostream>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "nutools/RandomUtils/NuRandomService.h"


#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TFile.h"
#include "TH1D.h"

namespace evwgh {
	class LEESignalUnfoldingWeightCalc : public WeightCalc
	{
		public:
			LEESignalUnfoldingWeightCalc();
			void Configure(fhicl::ParameterSet const& pset);
			std::vector<std::vector<double> > GetWeight(art::Event & e);

		private:    
			std::vector<double> fWeightArray;
			std::string fGenieModuleLabel;
			int fNmultisims;
			int fNuType;

			double fWeights[200];

			DECLARE_WEIGHTCALC(LEESignalUnfoldingWeightCalc)
	};
	LEESignalUnfoldingWeightCalc::LEESignalUnfoldingWeightCalc()
	{
	}

	void LEESignalUnfoldingWeightCalc::Configure(fhicl::ParameterSet const& p)
	{    
		fGenieModuleLabel= p.get< std::string > ("genie_module_label");
		fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
		fNmultisims = pset.get<int>("number_of_multisims");
		fNuType     = pset.get<int>("nutype");		
		std::string unfoldedWeights = pset.get< std::string >("unfolded_weights_file");

		TFile fUnfold(unfoldedWeights.c_str(),"read");
		if(!fUnfold.IsOpen()){
			throw cet::exception(__FUNCTION__) << GetName()<<":: Warning: File did not open correctly."<< std::endl;
		}

		for (int ibin=0;ibin<200;ibin++) {
			fWeights[ibin]=(dynamic_cast<TH1D*> (fUnfold.Get("unfolded_ratio")))->GetBinContent(ibin+1);
		}
		fUnfold.Close();

	}

	std::vector<std::vector<double> > LEESignalUnfoldingWeightCalc::GetWeight(art::Event & e)
	{
		//calculate weight(s) here 
		std::vector<std::vector<double>> weight;

		// * MC flux information
		art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
		std::vector<art::Ptr<simb::MCFlux> > fluxlist;
		if (e.getByLabel(fGenieModuleLabel,mcfluxListHandle))
			art::fill_ptr_vector(fluxlist, mcfluxListHandle);
		else{ return weight;}

		// * MC truth information
		art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
		std::vector<art::Ptr<simb::MCTruth> > mclist;
		if (e.getByLabel(fGenieModuleLabel,mctruthListHandle))
			art::fill_ptr_vector(mclist, mctruthListHandle);
		else{return weight;}


		weight.resize(mclist.size());
		for (unsigned int inu=0;inu<mclist.size();inu++) {
			weight[inu].resize(fNmultisims);

			int bin=-9999;

			if( fluxlist[inu]->fntype == fNuType){

				for (int i=0;i<fNmultisims;i++) {
					double enu=mclist[inu]->GetNeutrino().Nu().E();
					bin=enu/0.05;     
					double lee_signal_scaling = 0.0;
					if(bin > 0 && bin < 200) lee_signal_scaling = fWeights[bin]; 
					weight[inu][i] = lee_signal_scaling;
				}

			}else{
				for (int i=0;i<fNmultisims;i++) {
					weight[inu][i]=0.0;
				}
			}


		}
		

		return weight;
	}
	REGISTER_WEIGHTCALC(LEESignalUnfoldingWeightCalc)
}
