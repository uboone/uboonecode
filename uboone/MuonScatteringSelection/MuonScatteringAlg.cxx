#include "MuonScatteringAlg.h"

#include <functional>
#include <unordered_map>
#include "TMath.h"

andy::MuonScatteringAlg::MuonScatteringAlg() : FVx(256.35), FVy(233.), FVz(1036.8), borderx(10), bordery(10), borderz(40) {
  std::cout << "made a new MuonScatteringAlg object" << std::endl;
}

bool andy::MuonScatteringAlg::inFV(double xyz[3]){
  double x = xyz[0];
  double y = xyz[1];
  double z = xyz[2];
  //double FVx = 256.35;
  //double FVy = 233.;
  //double FVz = 1036.8;
  //double borderx = 10.; 
  //double bordery = 10;
  //double borderz = 40.;
  if (x > (FVx - borderx)){ return false;}
  if (x < (borderx)){ return false;}
  if (y > (FVy/2. - bordery)){ return false;}
  if (y < (-1*FVx/2. + bordery)){ return false;}
  if (z > (FVz - borderz)){ return false;}
  if (z< (borderz)){return false;}

  return true;
}

bool andy::MuonScatteringAlg::fullyContained(recob::Track t){
  double xyz_start[3] = {t.Start().X(), t.Start().Y(), t.Start().Z()};
  double xyz_end[3] = {t.End().X(), t.End().Y(), t.End().Z()};
  if (not inFV(xyz_start)){return false;}
  if (not inFV(xyz_end)){return false;}

  return true;
}

double andy::MuonScatteringAlg::closestApproach(recob::Track t, double xyz[3]){
  double min_dist = 99999;
  unsigned int nTrajPoints = t.NumberTrajectoryPoints();
  for (unsigned int point(0); point < nTrajPoints; ++point){
    double track_x = t.TrajectoryPoint(point).position.X();
    double track_y = t.TrajectoryPoint(point).position.Y();
    double track_z = t.TrajectoryPoint(point).position.Z();
    double dist_x = (xyz[0]-track_x);
    double dist_y = (xyz[1]-track_y);
    double dist_z = (xyz[2]-track_z);
    double dist = TMath::Sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
    if (dist < min_dist){
      min_dist = dist;
    }
  }
  return min_dist;
}

double andy::MuonScatteringAlg::absDistance(double * a, double * b){
  double x = a[0] - b[0];
  double y = a[1] - b[1];
  double z = a[2] - b[2];
  double distsq = x*x + y*y + z*z;
  double dist = TMath::Sqrt(distsq);
  return dist;
}

TVector3 andy::MuonScatteringAlg::TPCentryPoint(const simb::MCParticle & part){
  auto traj = part.Trajectory();
  for (unsigned int i(0); i < traj.size() ; i++){
    TVector3 pos = traj.Position(i).Vect();
    double posarray[3] = {pos.X(), pos.Y(), pos.Z()}; 
    if (inFV(posarray)){
      return pos;
    }
  }
  TVector3 y(-999,-999,-999);
  return y;
}

//simb::MCParticle GetMCParticleFromRecoTrack(art::Handle<recob::Track> &t, art::Event const & e){
//  art::InputTag tag("pandoraNu");
//  art::FindManyP<recob::Hit>      fmht(t, e, tag);
//}
