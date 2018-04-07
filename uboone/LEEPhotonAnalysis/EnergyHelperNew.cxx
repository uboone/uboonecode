#ifndef ENERGYHELPERNEW_CXX
#define ENERGYHELPERNEW_CXX

#include "EnergyHelperNew.h"

namespace lee
{

EnergyHelperNew::EnergyHelperNew() {
  detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
}


  /*
void EnergyHelperNew::trackResiduals(const art::Event &e,
                                  std::string _pfp_producer,
                                  art::Ptr<recob::Track> candidate_track,
                                  double &mean,
                                  double &std)
{
  lar_pandora::TrackVector allPfParticleTracks;
  lar_pandora::TracksToHits trackToHitsMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, trackToHitsMap);
  
  std::vector<TVector3> hit_v;   // a vec of hits from coll plane
  std::vector<TVector3> track_v; // a vec of hits from coll plane

  // Collect hits
  auto iter = trackToHitsMap.find(candidate_track);
  if (iter != trackToHitsMap.end())
  {
    std::vector<art::Ptr<recob::Hit>> hits = iter->second;
    for (auto hit : hits)
    {
      if (hit->View() == 2)
      {
        TVector3 h(hit->WireID().Wire, hit->PeakTime(), 0);
        //std::cout << "emplacing hit with wire " << h.X() << ", and time " << h.Y() << std::endl;
        hit_v.emplace_back(h);
      }
    }
  }

  // Collect track points
  for (size_t i = 0; i < candidate_track->NumberTrajectoryPoints(); i++)
  {
    try
    {
      if (candidate_track->HasValidPoint(i))
      {
        TVector3 trk_pt = candidate_track->LocationAtPoint(i);
        double wire = geo->NearestWire(trk_pt, 2);
        double time = detprop->ConvertXToTicks(trk_pt.X(), geo::PlaneID(0, 0, 2));
        TVector3 p(wire, time, 0.);
        //std::cout << "emplacing track point on wire " << p.X() << ", and time " << p.Y() << std::endl;
        track_v.emplace_back(p);
      }
    }
    catch (...)
    {
      continue;
    }
  }

  ubana::TrackQuality _track_quality;
  _track_quality.SetTrackPoints(track_v);
  _track_quality.SetHitCollection(hit_v);
  std::pair<double, double> residual_mean_std = _track_quality.GetResiduals();
  
  if (!isnan(residual_mean_std.first) && !isnan(residual_mean_std.second)) {
    mean = residual_mean_std.first;
    std = residual_mean_std.second;
  }

}
  */


void EnergyHelperNew::energyFromHits(recob::PFParticle const &pfparticle,
                                    std::vector<int>    &nHits,
                                    std::vector<double> &pfenergy,
                                    const art::Event &evt,
                                    std::string _pfp_producer)
{
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> cluster_pfparticle_ass(pfparticle_handle, evt,
                                                        _pfp_producer);
  std::vector<art::Ptr<recob::Cluster>> clusters =
      cluster_pfparticle_ass.at(pfparticle.Self());

      std::vector<double> _gain;
      if (evt.isRealData()){
        _gain = _data_gain;
      }  
      else{
        _gain = _mc_gain;
      }

  nHits.resize(3);
  pfenergy.resize(3);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    art::FindManyP<recob::Hit> hit_cluster_ass(cluster_handle, evt,
                                               _pfp_producer);
    std::vector<art::Ptr<recob::Hit>> hits =
        hit_cluster_ass.at(clusters[icl].key());

    for (size_t ihit = 0; ihit < hits.size(); ++ihit)
    {
      auto plane_nr = hits[ihit]->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;

      pfenergy[plane_nr] += hits[ihit]->Integral() * _gain[plane_nr] * work_function / recombination_factor /1000; // convert MeV to GeV
      nHits[plane_nr] ++;
    }
  }
}


  std::vector<double> EnergyHelperNew::trackEnergy_dedx(const art::Ptr<recob::Track> &track,
							const art::Event &evt,
							std::string _pfp_producer)
{

  auto const &track_handle = evt.getValidHandle<std::vector<recob::Track>>(_pfp_producer);
  art::FindManyP<anab::Calorimetry> calo_track_ass(track_handle, evt,"pandoraNucalo");

  const std::vector<art::Ptr<anab::Calorimetry>> calos = calo_track_ass.at(track->ID());

  std::vector<double> result(3, 0);

  for (size_t ical = 0; ical < calos.size(); ++ical)
  {
    if (!calos[ical])
      continue;

    if (!calos[ical]->PlaneID().isValid)
      continue;

    int planenum = calos[ical]->PlaneID().Plane;

    if (planenum < 0 || planenum > 2)
      continue;
    if (result.at(planenum) != 0)
      continue;
    
    result.at(planenum) = calos[ical]->KineticEnergy() / 1000;  // convert to GeV

  }

  return result;

}

void EnergyHelperNew::nHits(size_t pfp_id,
                         const art::Event &evt,
                         std::vector<int> &nHits,
                         std::string _pfp_producer)
{

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);
  nHits.resize(3);

  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, _pfp_producer);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    nHits[clusters[icl]->Plane().Plane] = hits_per_clusters.at(clusters[icl].key()).size();
  }
}

void EnergyHelperNew::PCA(size_t pfp_id,
                       const art::Event &evt,
                       std::vector<std::vector<double>> &pca_planes,
                       std::string _pfp_producer)
{

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);

  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, _pfp_producer);

  double drift = detprop->DriftVelocity() * 1e-3;
  double fromTickToNs = 4.8 / detprop->ReadOutWindowSize() * 1e6;
  double wireSpacing = 0.3;

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    TPrincipal fPrincipal(2, "D");

    std::vector<art::Ptr<recob::Hit>> hits = hits_per_clusters.at(clusters[icl].key());

    for (auto &hit : hits)
    {
      double data[2];
      double w = hit->WireID().Wire * wireSpacing;
      double t = fromTickToNs * drift * hit->PeakTime();
      data[0] = w;
      data[1] = t;
      fPrincipal.AddRow(data);
    }

    fPrincipal.MakePrincipals();
    pca_planes[clusters[icl]->Plane().Plane][0] = (*fPrincipal.GetEigenValues())[0];
    pca_planes[clusters[icl]->Plane().Plane][1] = (*fPrincipal.GetEigenValues())[1];
  }
}

void EnergyHelperNew::dQdx(size_t pfp_id,
                        const art::Event &evt,
                        std::vector<double> &dqdx,
                        std::vector<double> &dqdx_hits,
                        double m_dQdxRectangleLength,
                        double m_dQdxRectangleWidth,
                        std::string _pfp_producer)
{

  std::vector<double> _gain;
  if (evt.isRealData())
    _gain = _data_gain;
  else
    _gain = _mc_gain;

  detinfo::DetectorProperties const *detprop =
      lar::providerFrom<detinfo::DetectorPropertiesService>();
  // art::ServiceHandle<geo::Geometry> geom;

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  TVector3 pfp_dir;

  //For a shower
  if (pfparticle_handle->at(pfp_id).PdgCode() == 11)
  {

    art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const &shower_obj = shower_per_pfpart.at(pfp_id);

    try {
      pfp_dir.SetX(shower_obj->Direction().X());
      pfp_dir.SetY(shower_obj->Direction().Y());
      pfp_dir.SetZ(shower_obj->Direction().Z());
    } catch (...) {
      return;
    }
  }
  // For a track
  else
  {

    art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, _pfp_producer);
    auto const &track_obj = track_per_pfpart.at(pfp_id);

    pfp_dir.SetX(track_obj->StartDirection().X());
    pfp_dir.SetY(track_obj->StartDirection().Y());
    pfp_dir.SetZ(track_obj->StartDirection().Z());
  }


  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt,
                                                     _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt,
                                               _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters =
      clusters_per_pfpart.at(pfp_id);

  double drift = detprop->DriftVelocity() * 1e-3;
  //std::cout << "[dQdx] Clusters size " << clusters.size() << std::endl;

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {

    std::vector<art::Ptr<recob::Hit>> hits =
        hits_per_clusters.at(clusters[icl].key());

    /*
    std::cout << "[dQdx] "
              << "Cluster " << icl << std::endl;

    std::cout << "[dQdx] "
              << "Wire coordinate " << clusters[icl]->StartWire() << std::endl;
    std::cout << "[dQdx] "
              << "Tick coordinate " << clusters[icl]->StartTick() << std::endl;
    */
    // TODO Use variable from detector properties!
    // To get the time in ns -> 4.8 ms / 9600 ticks * 1e6 = 500
    // 0.3 wire spacing

    double fromTickToNs = 4.8 / detprop->ReadOutWindowSize() * 1e6;
    double wireSpacing = 0.3;

    std::vector<double> cluster_start = {
        clusters[icl]->StartWire() * wireSpacing,
        drift * clusters[icl]->StartTick() * fromTickToNs};
    std::vector<double> cluster_end = {clusters[icl]->EndWire() * wireSpacing,
                                       drift * clusters[icl]->EndTick() *
                                           fromTickToNs};

    double cluster_length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) +
                                 pow(cluster_end[1] - cluster_start[1], 2));
    if (cluster_length <= 0)
      continue;

    std::vector<double> cluster_axis = {cos(clusters[icl]->StartAngle()),
                                        sin(clusters[icl]->StartAngle())};

    // Build rectangle 4 x 1 cm around the cluster axis
    std::vector<std::vector<double>> points;
    geoHelper.buildRectangle(m_dQdxRectangleLength, m_dQdxRectangleWidth,
                             cluster_start, cluster_axis, points);

    /*
    std::cout << "[dQdx] Point 1 " << points[0][0] << " " << points[0][1] << std::endl;
    std::cout << "[dQdx] Point 2 " << points[1][0] << " " << points[1][1] << std::endl;
    std::cout << "[dQdx] Point 3 " << points[2][0] << " " << points[2][1] << std::endl;
    std::cout << "[dQdx] Point 4 " << points[3][0] << " " << points[3][1] << std::endl;
    */

    std::vector<double> dqdxs;

    bool first = true;

    for (auto &hit : hits)
    {

      // std::cout << "[PandoraLEE] Hit wire ID " << hit->WireID().Wire <<
      // std::endl;
      // std::cout << "Hit peak time " << hit->PeakTime() << std::endl;

      std::vector<double> hit_pos = {hit->WireID().Wire * wireSpacing,
                                     fromTickToNs * drift * hit->PeakTime()};

      double pitch =
          geoHelper.getPitch(pfp_dir, clusters[icl]->Plane().Plane);

      bool is_within = geoHelper.isInside(hit_pos, points);

      // Hit within the rectangle. The function considers points on the border
      // as outside, so we manually add the first point

      if (is_within || first)
      {
        double q = hit->Integral() * _gain[clusters[icl]->Plane().Plane];
        dqdxs.push_back(q / pitch);
        if (clusters[icl]->Plane().Plane == 2)
        {
          dqdx_hits.push_back(q / pitch);
        }
      }
      first = false;

      // std::cout << "[dQdx] Hit pos " << is_within << " " << hit_pos[0] << " " << hit_pos[1] << std::endl;
    }

    // Get the median
    size_t n = dqdxs.size() / 2;

    std::nth_element(dqdxs.begin(), dqdxs.begin() + n, dqdxs.end());
    if (n > 0)
    {
      /*
      std::cout << "[dQdx] Plane dQdx " << clusters[icl]->Plane().Plane << " "
                << icl << " "
                << dqdxs[n] << " " << dqdxs[n] * (23. / 1e6) / 0.62
                << std::endl;
      */
      dqdx[clusters[icl]->Plane().Plane] = dqdxs[n];
    }
  }
}

void EnergyHelperNew::dEdxFromdQdx(std::vector<double> &dedx,
                                std::vector<double> &dqdx)
{
  for (size_t i = 0; i < dqdx.size(); i++)
  {
    if (dqdx[i] > 0)
      dedx[i] = dqdx[i] * (work_function) / recombination_factor;
  }
}

} // namespace lee

#endif
