#pragma once

// system include files
#include <cmath>
#include <memory>
#include <vector>

// user include files
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackerRecHit2D/interface/VectorHit2D.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "RecoTracker/PixelSeeding/interface/OrderedHitTriplet.h"
#include "RecoTracker/TkSeedingLayers/interface/OrderedSeedingHits.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoTracker/PixelSeeding/interface/OrderedHitTriplet.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "DataFormats/TrackerRecHit2D/interface/VectorHit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/ClusterParameterEstimator.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkClonerImpl.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
//
// class declaration
//


/// helper struct to represent the seeding grid.
/// For cosmics, bin in cartesian coordinates, sort along y. 
class SeedingGrid{
  public: 
    SeedingGrid (int nBinsX, double xmin, double xmax, int nBinsY, double ymin, double ymax,int nBinsZ, double zmin, double zmax ):
      nBinsX_(nBinsX),
      xmin_(xmin), 
      xmax_(xmax), 
      nBinsY_(nBinsY),
      ymin_(ymin), 
      ymax_(ymax), 
      nBinsZ_(nBinsZ),
      zmin_(zmin), 
      zmax_(zmax){
        recHits_.resize(nBinsX*nBinsY*nBinsZ);
      }

    void addHit(const BaseTrackerRecHit* h){
      recHits_.at(getBin(binX(h->globalPosition().x()), binY(h->globalPosition().y()), binZ(h->globalPosition().z()))).push_back(h);
    }

    int binX(double x)const {
      return std::clamp<int>(std::floor(( x- xmin_) / (xmax_ - xmin_) * nBinsX_),0,nBinsX_-1);
    }
    int binY(double y)const {
      return std::clamp<int>(std::floor(( y- ymin_) / (ymax_ - ymin_) * nBinsY_),0,nBinsY_-1);
    }
    int binZ(double z)const {
      return std::clamp<int>(std::floor(( z- zmin_) / (zmax_ - zmin_) * nBinsZ_),0,nBinsZ_-1);
    }

    const std::vector<const BaseTrackerRecHit*> & getHits(int binX, int binY, int binZ) const{
      if (getBin(binX,binY,binZ) >= (int)recHits_.size()) std::cout << " request bin "<<binX<<" "<<binY<<" "<<binZ<<" mapping to "<<getBin(binX,binY,binZ)<<std::endl; 
      return recHits_.at(getBin(binX,binY,binZ));

    }

    void sort(){
      for (auto & vec : recHits_){
        std::sort(vec.begin(),vec.end(),[](const BaseTrackerRecHit* h1, const BaseTrackerRecHit* h2){return h1->globalPosition().y() < h2->globalPosition().y(); });
      }
    }

    int getBin(int bx, int by, int bz) const{
      // std::cout << " bins "<<bx<<" "<<by<<" "<<bz<<"  map to "<<bz + nBinsZ_ * (by + nBinsY_ * bx);
      return bz + nBinsZ_ * (by + nBinsY_ * bx);
    }

    int nBinsX() const {return nBinsX_;}
    int nBinsY() const {return nBinsY_;}
    int nBinsZ() const {return nBinsZ_;}

    double xmin() const {return  xmin_;}
    double xmax() const {return  xmax_;}
    double ymin() const {return  ymin_;}
    double ymax() const {return  ymax_;}
    double zmin() const {return  zmin_;}
    double zmax() const {return  zmax_;}

  private: 
    int nBinsX_;
    double xmin_; 
    double xmax_; 
    int nBinsY_;
    double ymin_; 
    double ymax_; 
    int nBinsZ_;
    double zmin_; 
    double zmax_; 

    std::vector<std::vector<const BaseTrackerRecHit*>> recHits_ = {};
};


class CosmicGridTripletSeeder : public edm::stream::EDProducer<> {
public:
    explicit CosmicGridTripletSeeder(const edm::ParameterSet& par);
    ~CosmicGridTripletSeeder() override = default;


    using protoSeed = OrderedHitTriplet; 

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void produce(edm::Event &e, const edm::EventSetup &c) override;

private:
    // event state, encapsuled to ensure thread safety 
    struct TripletSeederEventState{
        SeedingGrid grid; 
        const MagneticField *magfield = nullptr;
        const TrackerGeometry *tracker = nullptr;
        TkClonerImpl cloner;  // FIXME
        std::unique_ptr<KFUpdator> theUpdator = nullptr;
        std::unique_ptr<PropagatorWithMaterial> thePropagatorAl = nullptr;
        std::unique_ptr<PropagatorWithMaterial> thePropagatorOp = nullptr;
    };

     /// load conditions required for the event
     TripletSeederEventState initEventState(const edm::EventSetup &c);
     /// populate the seeding grid 
     bool populateGrid(const edm::Event& e, TripletSeederEventState & state); 
     /// using the seeding grid, form triplets
     bool formTriplets(const TripletSeederEventState & state, std::vector<protoSeed> & found);

     /// get the triplets for one particular cell combination
     void formTriplets(const std::vector<const BaseTrackerRecHit*> & topCands, 
                       const std::vector<const BaseTrackerRecHit*> & centerCands, 
                       const std::vector<const BaseTrackerRecHit*> & bottomCands, 
                       const TripletSeederEventState & state,
                       std::vector<protoSeed> & found, 
                       std::unordered_multiset<const BaseTrackerRecHit*> & trackUsage);
     /// fit the triplets with kalman into trajectory seeds
     bool fitTriplets(const TripletSeederEventState & state, const std::vector<protoSeed>& seeds, TrajectorySeedCollection & output);
     /// fit of a single triplet into a trajectory seed 
     bool fitTriplet(const TripletSeederEventState & state, const protoSeed& seed, TrajectorySeedCollection & output);

    /// helper borrowed from SimpleCosmicBONSeeder - possibly refactor into helper class for deployment
    std::pair<GlobalVector, int> pqFromHelixFit(const GlobalPoint &inner,
                                            const GlobalPoint &middle,
                                            const GlobalPoint &outer,
                                            const MagneticField* magfield); 

     /// data dependencies 
     
    edm::EDGetTokenT<VectorHitCollection> vectorHitsToken_;  // vector hit collection
    edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> otRecHitsToken_;  // strip hits in the outer tracker
    edm::EDGetTokenT<SiPixelRecHitCollection> pixelRecHitsToken_;   // pixel hits

    /// condition dependencies 

    // B-field 
    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
    // geometry
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerToken_;

    /// tools 
    const edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> ttrhBuilderToken_;


};

