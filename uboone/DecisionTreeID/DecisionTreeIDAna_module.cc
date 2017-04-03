////////////////////////////////////////////////////////////////////////
// Class:       DecisionTreeIDAna
// Plugin Type: analyzer (art v2_05_00)
// File:        DecisionTreeIDAna_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "TString.h"
#include "TTree.h"

#include "uboone/DecisionTreeID/Algorithms/TrackFeatures.h"

#include <memory>
#include <iostream>

class DecisionTreeIDAna;

class DecisionTreeIDAna : public art::EDAnalyzer {
public:
  explicit DecisionTreeIDAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DecisionTreeIDAna(DecisionTreeIDAna const &) = delete;
  DecisionTreeIDAna(DecisionTreeIDAna &&) = delete;
  DecisionTreeIDAna & operator = (DecisionTreeIDAna const &) = delete;
  DecisionTreeIDAna & operator = (DecisionTreeIDAna &&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const & e) override;

private:

  // Pointer to the feature algorithm
  ::dtfeatures::TrackFeatures fTrackFeatures;

  std::string _trackmodulelabel;

  bool        _getlabel;
  bool        _getprediction;
  std::string _dtassoclabel;

  std::string _fcsvname;
  bool        _writecsv;

  TTree* _tree1;
  int _run, _subrun, _event, _trackid;
  float _nhits,_length,_starty,_startz,_endy,_endz;
  float _theta,_phi,_distlenratio;
  float _startdqdx,_enddqdx,_dqdxdiff,_dqdxratio;
  float _totaldqdx,_averagedqdx;
  float _cosmicscore,_coscontscore;
  float _pidpida,_pidchi,_cfdistance;
  float _label;
  float _predict_p,_predict_mu,_predict_pi,_predict_em,_predict_cos;

  std::ofstream fout;

  const anab::CosmicTagID_t TAGID_P  = anab::CosmicTagID_t::kGeometry_YY;
  const anab::CosmicTagID_t TAGID_MU = anab::CosmicTagID_t::kGeometry_YZ;
  const anab::CosmicTagID_t TAGID_PI = anab::CosmicTagID_t::kGeometry_ZZ;
  const anab::CosmicTagID_t TAGID_EM = anab::CosmicTagID_t::kGeometry_XX;
  const anab::CosmicTagID_t TAGID_CS = anab::CosmicTagID_t::kGeometry_XY;
    
};


DecisionTreeIDAna::DecisionTreeIDAna(fhicl::ParameterSet const & p)
  : EDAnalyzer(p) 
{
  _trackmodulelabel = p.get<std::string>("TrackModuleLabel");
  _getlabel         = p.get<bool>("GetClassLabel", false);
  _getprediction    = p.get<bool>("GetDTPrediction", false);
  _dtassoclabel     = p.get<std::string>("DTAssocLabel", "decisiontreeid");
  _fcsvname         = p.get<std::string>("CSVFileOut");
  _writecsv         = p.get<bool>("WriteCSV", false);
   
  fTrackFeatures.Configure(p.get<fhicl::ParameterSet>("FeaturesConfig"));

  art::ServiceHandle<art::TFileService> tfs;
  _tree1 = tfs->make<TTree>("tree","");
  _tree1->Branch("run",    &_run,    "run/I");
  _tree1->Branch("subrun", &_subrun, "subrun/I");
  _tree1->Branch("event",  &_event,  "event/I");
  _tree1->Branch("trackid",&_trackid,"trackid/I");
  _tree1->Branch("nhits",&_nhits,"nhits/F");
  _tree1->Branch("length",&_length,"length/F");
  _tree1->Branch("starty",&_starty,"starty/F");
  _tree1->Branch("startz",&_startz,"startz/F");
  _tree1->Branch("endy",&_endy,"endy/F");
  _tree1->Branch("endz",&_endz,"endz/F");
  _tree1->Branch("theta",&_theta,"theta/F");
  _tree1->Branch("phi",&_phi,"phi/F");
  _tree1->Branch("distlenratio",&_distlenratio,"distlenratio/F");
  _tree1->Branch("startdqdx",&_startdqdx,"startdqdx/F");
  _tree1->Branch("enddqdx",&_enddqdx,"enddqdx/F");
  _tree1->Branch("dqdxdiff",&_dqdxdiff,"dqdxdiff/F");
  _tree1->Branch("dqdxratio",&_dqdxratio,"dqdxratio/F");
  _tree1->Branch("totaldqdx",&_totaldqdx,"totaldqdx/F");
  _tree1->Branch("averagedqdx",&_averagedqdx,"averagedqdx/F");
  _tree1->Branch("cosmicscore",&_cosmicscore,"cosmicscore/F");
  _tree1->Branch("coscontscore",&_coscontscore,"coscontscore/F");
  _tree1->Branch("pidpida",&_pidpida,"pidpida/F");
  _tree1->Branch("pidchi",&_pidchi,"pidchi/F");
  _tree1->Branch("cfdistance",&_cfdistance,"cfdistance/F");
  _tree1->Branch("label",&_label,"label/F");
  _tree1->Branch("predict_p",&_predict_p,"predict_p/F");
  _tree1->Branch("predict_mu",&_predict_mu,"predict_mu/F");
  _tree1->Branch("predict_pi",&_predict_pi,"predict_pi/F");
  _tree1->Branch("predict_em",&_predict_em,"predict_em/F");
  _tree1->Branch("predict_cos",&_predict_cos,"predict_cos/F");

}

void DecisionTreeIDAna::analyze(art::Event const & e)
{

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // open csv for writing 
  if(_writecsv) fout.open(_fcsvname.c_str(),std::ofstream::out | std::ofstream::app);

  // recover handle for tracks that we want to analyze
  art::Handle< std::vector<recob::Track> > trackVecHandle;
  e.getByLabel(_trackmodulelabel, trackVecHandle);

  if(trackVecHandle.isValid())
  {
    art::FindManyP<anab::CosmicTag> dtAssns(trackVecHandle, e, _dtassoclabel); 

    std::vector< std::vector<float> > evtdata = fTrackFeatures.CreateFeatures(e, trackVecHandle);
    std::vector<float> evtlabels;
    if(_getlabel) evtlabels = fTrackFeatures.ClassLabel(e, trackVecHandle);

    // Loop over input tracks
    for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
    {
      art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);

      std::vector<float> trkdata = evtdata.at(trkIdx);

      if(_getlabel) _label = evtlabels.at(trkIdx);
      else _label = -9999.;

      _predict_p   = -9999.;
      _predict_mu  = -9999.;
      _predict_pi  = -9999.;
      _predict_em  = -9999.;
      _predict_cos = -9999.;

      if(_getprediction && dtAssns.isValid())
      {

        std::vector< art::Ptr<anab::CosmicTag> > dtVec = dtAssns.at(trk.key());
        for(auto const& dttag : dtVec)
        {
          if(dttag->CosmicType() == TAGID_P)  _predict_p   = dttag->CosmicScore();
          else if(dttag->CosmicType() == TAGID_MU) _predict_mu  = dttag->CosmicScore();
          else if(dttag->CosmicType() == TAGID_PI) _predict_pi  = dttag->CosmicScore();
          else if(dttag->CosmicType() == TAGID_EM) _predict_em  = dttag->CosmicScore();
          else if(dttag->CosmicType() == TAGID_CS) _predict_cos = dttag->CosmicScore();
        }
      }

      if(_writecsv)
      {
        for(auto it = trkdata.begin(); it != trkdata.end(); ++it)
        {
          if(std::next(it) != trkdata.end()) fout << *it << ",";
          else fout << *it;
        }
        if(_getlabel) fout << "," << _label;
        if(_getprediction) fout << "," << _predict_p  << ","
                                       << _predict_mu << ","
                                       << _predict_pi << ","
                                       << _predict_em << ","
                                       << _predict_cos;

        fout << std::endl;
      }

      _trackid      = trk->ID();
      _nhits        = trkdata[0];
      _length       = trkdata[1];
      _starty       = trkdata[2];
      _startz       = trkdata[3];
      _endy         = trkdata[4];
      _endz         = trkdata[5];
      _theta        = trkdata[6];
      _phi          = trkdata[7];
      _distlenratio = trkdata[8];
      _startdqdx    = trkdata[9];
      _enddqdx      = trkdata[10];
      _dqdxdiff     = trkdata[11];
      _dqdxratio    = trkdata[12];
      _totaldqdx    = trkdata[13];
      _averagedqdx  = trkdata[14];
      _cosmicscore  = trkdata[15];
      _coscontscore = trkdata[16];
      _pidpida      = trkdata[17];
      _pidchi       = trkdata[18];
      _cfdistance   = trkdata[19];

      _tree1->Fill();

    }
  }
  if(_writecsv) fout.close();
}

void DecisionTreeIDAna::beginJob()
{
  if(_writecsv)
  {
    // open csv and add header
    fout.open(_fcsvname.c_str(),std::ofstream::out | std::ofstream::app);

    fout << "nhits,length,starty,startz,endy,endz,theta,phi,"
            "distlenratio,startdqdx,enddqdx,dqdxdiff,dqdxratio,"
            "totaldqdx,averagedqdx,cosmicscore,coscontscore,"
            "pidpida,pidchi,cfdistance";

    if(_getlabel) fout << ",label";
    if(_getprediction)
    {
      fout << ",predict_p,predict_mu,predict_pi,predict_em,predict_cos";
    }
    fout << std::endl;
    fout.close();
  }
}
DEFINE_ART_MODULE(DecisionTreeIDAna)
