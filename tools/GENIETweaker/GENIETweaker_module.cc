// Fake filter module (all events pass) that will manually tweak the GENIE
// configuration. Intended for use in validation of larsim's EventWeight module
// (tweak & generate, compare results with reweighting).
//
// Steven Gardiner <gardiner@fnal.gov>

// Standard library includes
#include <iostream>
#include <limits>
#include <map>

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// GENIE includes (v3 only for now)
//#ifdef GENIE_PRE_R3
// TODO: add GENIE v2 includes
//#else
// Use these includes for GENIE v3
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/EventGeneratorList.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"
#include "RwFramework/GSyst.h"
#include "RwFramework/GSystUncertainty.h"
//#endif

// Helper functions local to this source file
namespace {

  constexpr double BOGUS = std::numeric_limits<double>::lowest();

  genie::Algorithm* find_algorithm(const RgAlg& rg_alg) {
    const genie::Algorithm* alg = genie::AlgFactory::Instance()
      ->GetAlgorithm( rg_alg.name, rg_alg.config );
    return const_cast< genie::Algorithm* >( alg );
  }

  genie::Algorithm* find_algorithm(const std::string& name_slash_config) {
    genie::AlgConfigPool* conf_pool = genie::AlgConfigPool::Instance();
    genie::Registry* gpl = conf_pool->GlobalParameterList();
    RgAlg temp_rg_alg = gpl->GetAlg( name_slash_config );
    return find_algorithm( temp_rg_alg );
  }

  //genie::Algorithm* find_event_generator_module(
  //  const genie::GEVGDriver& evg_driver,
  //  const std::string& ev_gen_config, const std::string& module_name)
  //{
  //  for (const genie::EventGeneratorI* gen : *evg_driver.EventGenerators()) {
  //    // All of the generators in the list will have
  //    // gen->Id().Name() == "genie::EventGenerator"
  //    if ( gen->Id().Config() == ev_gen_config ) {
  //      const genie::Registry& temp_reg = gen->GetConfig();
  //      temp_reg.Print( std::cout );
  //      RgInt num_modules = temp_reg.GetInt( "NModules" );
  //      for (RgInt m = 0; m < num_modules; ++m) {
  //        RgKey module_key = "Module-" + std::to_string(m);
  //        RgAlg alg = temp_reg.GetAlg( module_key );
  //        if (module_name == alg.name) {
  //          return const_cast< genie::Algorithm* >(
  //            genie::AlgFactory::Instance()->GetAlgorithm(alg.name, alg.config)
  //          );
  //        }
  //      }
  //    }
  //  }
  //  return nullptr;
  //}

}

namespace sim {
  class GENIETweaker;
}

class sim::GENIETweaker : public art::EDFilter {
public:
  explicit GENIETweaker(const fhicl::ParameterSet& p);
  virtual ~GENIETweaker();

  bool filter(art::Event& evt) override;
  // Do the tweak in begin run to ensure that it'll always come after
  // GENIEHelper initializes the tune (which is done in GENIEGen::beginJob()).
  // This ensures that GENIEHelper won't overwrite our tweaks.
  bool beginRun(art::Run& run) override;
  //void beginJob() override;

protected:

  genie::Registry* get_registry_for_tweak_dial(genie::rew::GSyst_t dial_label);

  void tweak_parameter(genie::Registry& reg, const std::string& parameter_name,
    genie::rew::GSyst_t twk_dial_label, double twk_dial_value, double& untweaked_value,
    double& tweaked_value);

  void CheckTune(const std::string& tune_name);

  bool fCheckTweaksOnFilter;
  std::vector<genie::rew::GSyst_t> fTweakDialLabels;
  std::vector<double> fTweakDialValues;
  std::vector<double> fUntweakedParameterValues;
  std::vector<double> fTweakedParameterValues;

  static const std::map<genie::rew::GSyst_t, std::string> fDialLabelToParamName;
};


sim::GENIETweaker::GENIETweaker(const fhicl::ParameterSet& pset)
  : art::EDFilter(pset)
{
  //#ifndef GENIE_PRE_R3
  // Check that the GENIE tune is configured (v3+ only)
  this->CheckTune( pset.get<std::string>("tune_name", "${GENIE_XSEC_TUNE}") );
  //#endif

  fCheckTweaksOnFilter = pset.get<bool>("check_on_filter", false);

  // Tweak dials are named using the same codes as grwght1p
  // (or, equivalently, genie::rew::GSyst::FromString())
  auto tweak_dial_names = pset.get< std::vector<std::string> >( "tweak_dial_names" );
  for (const auto& dial_name : tweak_dial_names) {
    // Get the GSyst_t enum value corresponding to the current name
    genie::rew::GSyst_t dial_label = genie::rew::GSyst::FromString( dial_name );
    fTweakDialLabels.push_back( dial_label );
  }

  // Tweak dial values are interpreted as in GENIE's native Reweight package
  // (0 -> untweaked, +1 -> raised by nominal 1-sigma error,
  // -1 -> lowered by nominal 1-sigma error)
  // The 1-sigma errors are obtained from the (sadly hard-coded) values owned
  // by the genie::rew::GSystUncertainty singleton class (just as in grwght1p)
  fTweakDialValues = pset.get< std::vector<double> >( "tweak_dial_values" );

}

sim::GENIETweaker::~GENIETweaker() {
}

bool sim::GENIETweaker::filter(art::Event& /*evt*/) {

  // Don't actually do anything to the event, we just need the tweaking
  // to be done in beginRun() to get GENIE adjusted properly.

  if ( fCheckTweaksOnFilter ) {

    LOG_INFO("sim::GENIETweaker") << "Checking the " << fTweakDialLabels.size()
      << " parameters that were tweaked.";

    for (size_t d = 0; d < fTweakDialLabels.size(); ++d) {
      genie::rew::GSyst_t dial_label = fTweakDialLabels.at( d );
      genie::Registry* reg = get_registry_for_tweak_dial( dial_label );
      std::string param_name = fDialLabelToParamName.at( dial_label );
      double param_value = reg->GetDouble( param_name );

      // If the tweaks made by this module have been reverted or otherwise altered,
      // then reject the current event by returning false, and print an error message.
      double untweaked_value = fUntweakedParameterValues.at( d );
      if ( fTweakDialValues.at( d ) != 0. && param_value == untweaked_value ) {
        LOG_ERROR("sim::GENIETweaker") << "The parameter " << param_name << " controlled by"
        << " the reweight tweak dial " << genie::rew::GSyst::AsString( dial_label )
        << " was reset to its untweaked value " << untweaked_value << '\n';
        return false;
      }

      double tweaked_value = fTweakedParameterValues.at( d );
      if ( param_value != tweaked_value ) {
        LOG_ERROR("sim::GENIETweaker") << "The parameter " << param_name << " controlled by"
        << " the reweight tweak dial " << genie::rew::GSyst::AsString( dial_label )
        << " differs from its tweaked value " << tweaked_value << '\n';
        return false;
      }

      LOG_INFO("sim::GENIETweaker") << "The parameter " << param_name << " controlled"
        << " by the reweight tweak dial " << genie::rew::GSyst::AsString( dial_label )
        << " currently has the value " << param_value << " (changed from its original"
        << " value of " << untweaked_value << ')';
    }
  }

  return true;
}

// The manual tweaking is done here
bool sim::GENIETweaker::beginRun(art::Run& /*run*/) {

  static bool already_tweaked = false;
  if ( already_tweaked ) {
    LOG_INFO("sim::GENIETweaker") << "Values have already been"
      << "tweaked. GENIETweaker is being inappropriately run multiple times.\n";
    return true;
  }

  fTweakedParameterValues.resize( fTweakDialLabels.size() );
  fUntweakedParameterValues.resize( fTweakDialLabels.size() );
  for (size_t d = 0; d < fTweakDialLabels.size(); ++d) {

    genie::rew::GSyst_t dial_label = fTweakDialLabels.at( d );
    double dial_value = fTweakDialValues.at( d );

    genie::Registry* reg = get_registry_for_tweak_dial( dial_label );
    std::string param_name = fDialLabelToParamName.at( dial_label );
    double untweaked_value = BOGUS;
    double tweaked_value = BOGUS;
    tweak_parameter(*reg, param_name, dial_label, dial_value, untweaked_value, tweaked_value);

    fUntweakedParameterValues[ d ] = untweaked_value;
    fTweakedParameterValues[ d ] = tweaked_value;
  }

  already_tweaked = true;

  return true;
}

// Copied from larsim's EventWeight module
// TODO: reduce code duplication here
void sim::GENIETweaker::CheckTune(const std::string& tune_name) {
//// The tune configuration only needs to be checked for GENIE v3+
//#ifndef GENIE_PRE_R3

  std::string fhicl_tune_name = tune_name;

  // The default tune name is ${GENIE_XSEC_TUNE}, which
  // should be converted into the value of the corresponding
  // enviornment variable, as is done below.
  if ( fhicl_tune_name.front() == '$' ) {
    // need to remove ${}'s
    std::string tuneEnvVar = fhicl_tune_name;
    std::string rmchars("$(){} ");
    // std::remove_if removes characters in [first,last) that are found
    //   within the rmchars string. It returns returns a past-the-end
    //   iterator for the new end of the range [funky!]
    // std::string::erase actually trims the string
    tuneEnvVar.erase( std::remove_if(tuneEnvVar.begin(), tuneEnvVar.end(),
      [&rmchars](const char& c) -> bool { return rmchars.find(c) != std::string::npos; }),
      tuneEnvVar.end() );

    const char* tune = std::getenv( tuneEnvVar.c_str() );
    if ( tune ) {
      mf::LogInfo("EventWeight") << "fhicl_tune_name started as '"
        << fhicl_tune_name << "' " << " (env: " << tuneEnvVar << "), "
        << " converted to " << tune;
      fhicl_tune_name = std::string(tune);
    } else {
      mf::LogError("EventWeight") << "fhicl_tune_name started as '"
        << fhicl_tune_name << "', " << " (env: " << tuneEnvVar << "), "
        << " but resolved to a empty string";
      throw cet::exception("UnresolvedTuneName")
        << "can't resolve TuneName: " << fhicl_tune_name;
    }
  }

  // If the XSecSplineList returns a non-empty string as the current tune name,
  // then genie::RunOpt::BuildTune() has already been called.
  std::string current_tune = genie::XSecSplineList::Instance()->CurrentTune();
  if ( current_tune.empty() ) {
    // We need to build the GENIE tune config
    mf::LogInfo("EventWeight") << "Configuring GENIE tune \""
      << fhicl_tune_name << '\"';

    // Constructor automatically calls grunopt->Init();
    genie::RunOpt* grunopt = genie::RunOpt::Instance();
    grunopt->SetTuneName( fhicl_tune_name );
    grunopt->BuildTune();
  }
  else {
    // It has already been built, so just check consistency
    if ( fhicl_tune_name != current_tune) {
      throw cet::exception("TuneNameMismatch") << "Requested GENIE tune \""
        << fhicl_tune_name << "\" does not match previously built tune \""
        << current_tune << '\"';
    }
  }

//#endif
}

void sim::GENIETweaker::tweak_parameter(genie::Registry& reg,
  const std::string& parameter_name, genie::rew::GSyst_t twk_dial_label,
  double twk_dial_value, double& untweaked_value, double& tweaked_value)
{
  untweaked_value = reg.GetDouble( parameter_name );

  int tweak_sign = 1;
  if ( twk_dial_value < 0. ) tweak_sign = -1;

  // Get nominal 1-sigma error for this tweak dial from
  // the genie::rew::GSystUncertainty singleton class
  auto* gsu = genie::rew::GSystUncertainty::Instance();
  double sigma = gsu->OneSigmaErr(twk_dial_label, tweak_sign);
  tweaked_value = untweaked_value + twk_dial_value*sigma;

  reg.Set( parameter_name, tweaked_value );

  LOG_INFO("sim::GENIETweaker") << "Changed " << parameter_name << " from "
    << untweaked_value << " to " << tweaked_value << " (" << (twk_dial_value < 0 ? "" : "+")
    << twk_dial_value << " sigma)";
}


genie::Registry* sim::GENIETweaker::get_registry_for_tweak_dial(
  genie::rew::GSyst_t dial_label)
{
  // Get a non-const pointer to the registry that owns each parameter
  // controlled by a tweak dial
  genie::Registry* reg = nullptr;
  switch (dial_label) {

    // NC quasielastic cross section parameters
    case genie::rew::kXSecTwkDial_MaNCEL:
    case genie::rew::kXSecTwkDial_EtaNCEL:
    {
      auto* nc_xsec = find_algorithm("XSecModel@genie::EventGenerator/QEL-NC");
      reg = nc_xsec->GetOwnedConfig();
    }
    break;

    // CC quasielastic cross section parameters
    case genie::rew::kXSecTwkDial_MaCCQE:
    {
      auto* cc_xsec = find_algorithm("XSecModel@genie::EventGenerator/QEL-CC");
      RgAlg ff_alg_id = cc_xsec->GetConfig().GetAlg( "FormFactorsAlg" );
      auto* ff_alg = find_algorithm( ff_alg_id );
      RgAlg axial_ff_model = ff_alg->GetConfig().GetAlg( "AxialFormFactorModel" );
      auto* axial_ff_model_alg = find_algorithm( axial_ff_model );

      reg = axial_ff_model_alg->GetOwnedConfig();
    }
    break;

    // CC resonance production cross section parameters
    case genie::rew::kXSecTwkDial_MaCCRES:
    case genie::rew::kXSecTwkDial_MvCCRES:
    {
      auto* cc_xsec = find_algorithm("XSecModel@genie::EventGenerator/RES-CC");
      reg = cc_xsec->GetOwnedConfig();
    }
    break;

    // NC resonance production cross section parameters
    case genie::rew::kXSecTwkDial_MaNCRES:
    case genie::rew::kXSecTwkDial_MvNCRES:
    {
      auto* nc_xsec = find_algorithm("XSecModel@genie::EventGenerator/RES-NC");
      reg = nc_xsec->GetOwnedConfig();
    }
    break;

    // Coherent pion production cross section parameters
    case genie::rew::kXSecTwkDial_MaCOHpi:
    case genie::rew::kXSecTwkDial_R0COHpi:
    {
      // TODO: tweak NC as well? reweight calculator appears to assume that CC and NC
      // models use the same parameters
      auto* coh_cc_xsec = find_algorithm("XSecModel@genie::EventGenerator/COH-CC");
      reg = coh_cc_xsec->GetOwnedConfig();
    }
    break;

    // Deep inelastic scattering cross section parameters
    case genie::rew::kXSecTwkDial_AhtBY:
    case genie::rew::kXSecTwkDial_BhtBY:
    case genie::rew::kXSecTwkDial_CV1uBY:
    case genie::rew::kXSecTwkDial_CV2uBY:
    {
      // TODO: tweak NC as well? reweight calculator appears to assume that CC and NC
      // models use the same parameters
      auto* cc_xsec = find_algorithm("XSecModel@genie::EventGenerator/DIS-CC");
      RgAlg sf_alg_id = cc_xsec->GetConfig().GetAlg( "SFAlg" );
      auto* sf_alg = find_algorithm( sf_alg_id );
      reg = sf_alg->GetOwnedConfig();
    }
    break;

    // Non-resonant background parameters
    case genie::rew::kXSecTwkDial_RvpCC1pi:
    case genie::rew::kXSecTwkDial_RvpCC2pi:
    case genie::rew::kXSecTwkDial_RvpNC1pi:
    case genie::rew::kXSecTwkDial_RvpNC2pi:
    case genie::rew::kXSecTwkDial_RvnCC1pi:
    case genie::rew::kXSecTwkDial_RvnCC2pi:
    case genie::rew::kXSecTwkDial_RvnNC1pi:
    case genie::rew::kXSecTwkDial_RvnNC2pi:
    case genie::rew::kXSecTwkDial_RvbarpCC1pi:
    case genie::rew::kXSecTwkDial_RvbarpCC2pi:
    case genie::rew::kXSecTwkDial_RvbarpNC1pi:
    case genie::rew::kXSecTwkDial_RvbarpNC2pi:
    case genie::rew::kXSecTwkDial_RvbarnCC1pi:
    case genie::rew::kXSecTwkDial_RvbarnCC2pi:
    case genie::rew::kXSecTwkDial_RvbarnNC1pi:
    case genie::rew::kXSecTwkDial_RvbarnNC2pi:
    {
      reg = const_cast< genie::Registry* >( genie::AlgConfigPool::Instance()
        ->CommonParameterList("NonResBackground") );
    }
    break;

    // Intranuke hA parameters
    // TODO: add elastic fate fractions for GENIE v2
    case genie::rew::kINukeTwkDial_MFP_N:
    case genie::rew::kINukeTwkDial_FrCEx_N:
    case genie::rew::kINukeTwkDial_FrInel_N:
    case genie::rew::kINukeTwkDial_FrAbs_N:
    case genie::rew::kINukeTwkDial_FrPiProd_N:
    case genie::rew::kINukeTwkDial_MFP_pi:
    case genie::rew::kINukeTwkDial_FrCEx_pi:
    case genie::rew::kINukeTwkDial_FrInel_pi:
    case genie::rew::kINukeTwkDial_FrAbs_pi:
    case genie::rew::kINukeTwkDial_FrPiProd_pi:
    {
      auto* fsi_alg = find_algorithm("HadronTransp-Model");
      reg = fsi_alg->GetOwnedConfig();

      // The scale factors for these are *optional* parameters, so they
      // may need to be inserted into the registry before proceeding
      std::string temp_param_name = fDialLabelToParamName.at( dial_label );
      if ( !reg->Exists( temp_param_name ) ) {
        // The default setting for these should be to scale by 1.
        reg->Set( temp_param_name, 1. );
      }
    }
    break;

    default:
    break;
  }

  return reg;
}

// Lookup table to convert from a genie::rew::GSyst_t tweak dial label
// into the parameter name used to store the value in a genie::Registry.
const std::map<genie::rew::GSyst_t, std::string>
sim::GENIETweaker::fDialLabelToParamName =
{
  { genie::rew::kXSecTwkDial_MaNCEL,   "QEL-Ma"       },
  { genie::rew::kXSecTwkDial_EtaNCEL,  "EL-Axial-Eta" },
  { genie::rew::kXSecTwkDial_MaCCQE,   "QEL-Ma"       },
  { genie::rew::kXSecTwkDial_MaCCRES,  "RES-Ma"       },
  { genie::rew::kXSecTwkDial_MvCCRES,  "RES-Mv"       },
  { genie::rew::kXSecTwkDial_MaNCRES,  "RES-Ma"       },
  { genie::rew::kXSecTwkDial_MvNCRES,  "RES-Mv"       },
  { genie::rew::kXSecTwkDial_MaCOHpi,  "COH-Ma"       },
  { genie::rew::kXSecTwkDial_R0COHpi,  "COH-Ro"       },
  { genie::rew::kXSecTwkDial_RvpCC1pi, "DIS-HMultWgt-vp-CC-m2" },
  { genie::rew::kXSecTwkDial_RvpCC2pi, "DIS-HMultWgt-vp-CC-m3" },
  { genie::rew::kXSecTwkDial_RvpNC1pi, "DIS-HMultWgt-vp-NC-m2" },
  { genie::rew::kXSecTwkDial_RvpNC2pi, "DIS-HMultWgt-vp-NC-m3" },
  { genie::rew::kXSecTwkDial_RvnCC1pi, "DIS-HMultWgt-vn-CC-m2" },
  { genie::rew::kXSecTwkDial_RvnCC2pi, "DIS-HMultWgt-vn-CC-m3" },
  { genie::rew::kXSecTwkDial_RvnNC1pi, "DIS-HMultWgt-vn-NC-m2" },
  { genie::rew::kXSecTwkDial_RvnNC2pi, "DIS-HMultWgt-vn-NC-m3" },
  { genie::rew::kXSecTwkDial_RvbarpCC1pi, "DIS-HMultWgt-vbp-CC-m2" },
  { genie::rew::kXSecTwkDial_RvbarpCC2pi, "DIS-HMultWgt-vbp-CC-m3" },
  { genie::rew::kXSecTwkDial_RvbarpNC1pi, "DIS-HMultWgt-vbp-NC-m2" },
  { genie::rew::kXSecTwkDial_RvbarpNC2pi, "DIS-HMultWgt-vbp-NC-m3" },
  { genie::rew::kXSecTwkDial_RvbarnCC1pi, "DIS-HMultWgt-vbn-CC-m2" },
  { genie::rew::kXSecTwkDial_RvbarnCC2pi, "DIS-HMultWgt-vbn-CC-m3" },
  { genie::rew::kXSecTwkDial_RvbarnNC1pi, "DIS-HMultWgt-vbn-NC-m2" },
  { genie::rew::kXSecTwkDial_RvbarnNC2pi, "DIS-HMultWgt-vbn-NC-m3" },
  { genie::rew::kINukeTwkDial_MFP_N,       "FSI-Nucleon-MFPScale"        },
  { genie::rew::kINukeTwkDial_FrCEx_N,     "FSI-Nucleon-FracCExScale"    },
  { genie::rew::kINukeTwkDial_FrInel_N,    "FSI-Nucleon-FracInelScale"   },
  { genie::rew::kINukeTwkDial_FrAbs_N,     "FSI-Nucleon-FracAbsScale"    },
  { genie::rew::kINukeTwkDial_FrPiProd_N,  "FSI-Nucleon-FracPiProdScale" },
  { genie::rew::kINukeTwkDial_MFP_pi,      "FSI-Pion-MFPScale"           },
  { genie::rew::kINukeTwkDial_FrCEx_pi,    "FSI-Pion-FracCExScale"       },
  { genie::rew::kINukeTwkDial_FrInel_pi,   "FSI-Pion-FracInelScale"      },
  { genie::rew::kINukeTwkDial_FrAbs_pi,    "FSI-Pion-FracAbsScale"       },
  { genie::rew::kINukeTwkDial_FrPiProd_pi, "FSI-Pion-FracPiProdScale"    },
  { genie::rew::kXSecTwkDial_AhtBY,  "BY-A"    },
  { genie::rew::kXSecTwkDial_BhtBY,  "BY-B"    },
  { genie::rew::kXSecTwkDial_CV1uBY, "BY-Cv1U" },
  { genie::rew::kXSecTwkDial_CV2uBY, "BY-Cv2U" }
};

DEFINE_ART_MODULE(sim::GENIETweaker)
