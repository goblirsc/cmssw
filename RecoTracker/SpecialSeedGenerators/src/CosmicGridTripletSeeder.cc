
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TEllipse.h>
#include <cmath>
#include <memory>
#include <set>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/VectorHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "RecoTracker/SpecialSeedGenerators/interface/CosmicGridTripletSeeder.h"

//
 

CosmicGridTripletSeeder::CosmicGridTripletSeeder(const edm::ParameterSet& iConfig)
    : vectorHitsToken_(consumes<VectorHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vectorHits"))),
      otRecHitsToken_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("OTRecHits"))),
      pixelRecHitsToken_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("PixelRecHits"))),
      magfieldToken_(esConsumes(iConfig.getParameter<edm::ESInputTag>("MagneticFieldRecord"))),
      trackerToken_(esConsumes()),
      ttrhBuilderToken_(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("TTRHBuilder")))) {
  produces<TrajectorySeedCollection>();
}

void CosmicGridTripletSeeder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("vectorHits", edm::InputTag("siPhase2VectorHits:accepted"));
  desc.addUntracked<edm::InputTag>("OTRecHits", edm::InputTag("Phase2TrackerRecHits"));
  desc.addUntracked<edm::InputTag>("PixelRecHits", edm::InputTag("siPixelRecHits"));
  desc.add<std::string>("TTRHBuilder", "WithTrackAngle");
  desc.add<edm::ESInputTag>("MagneticFieldRecord", edm::ESInputTag("", ""));
  descriptions.addWithDefaultLabel(desc);
}

std::pair<GlobalVector, int> CosmicGridTripletSeeder::pqFromHelixFit(const GlobalPoint& inner,
                                                                     const GlobalPoint& middle,
                                                                     const GlobalPoint& outer,
                                                                     const MagneticField* magfield) {
  std::cout << "DEBUG PZ =====" << std::endl;
  FastHelix helix(inner, middle, outer, magfield->nominalValue(), magfield);
  GlobalVector gv = helix.stateAtVertex().momentum();  // status on inner hit
  std::cout << "FastHelix P = " << gv << "\n";
  std::cout << "FastHelix Q = " << helix.stateAtVertex().charge() << "\n";

  // My attempt (with different approx from FastHelix)
  // 1) fit the circle
  FastCircle theCircle(inner, middle, outer);
  double rho = theCircle.rho();

  // 2) Get the PT
  GlobalVector tesla = magfield->inTesla(middle);
  double pt = 0.01 * rho * (0.3 * tesla.z());

  // 3) Get the PX,PY at OUTER hit (VERTEX)
  double dx1 = outer.x() - theCircle.x0();
  double dy1 = outer.y() - theCircle.y0();
  double py = pt * dx1 / rho, px = -pt * dy1 / rho;
  if (px * (middle.x() - outer.x()) + py * (middle.y() - outer.y()) < 0.) {
    px *= -1.;
    py *= -1.;
  }

  // 4) Get the PZ through pz = pT*(dz/d(R*phi)))
  double dz = inner.z() - outer.z();
  double sinphi = (dx1 * (inner.y() - theCircle.y0()) - dy1 * (inner.x() - theCircle.x0())) / (rho * rho);
  double dphi = std::abs(std::asin(sinphi));
  double pz = pt * dz / (dphi * rho);

  int myq = ((theCircle.x0() * py - theCircle.y0() * px) / tesla.z()) > 0. ? +1 : -1;

  std::pair<GlobalVector, int> mypq(GlobalVector(px, py, pz), myq);

  std::cout << "Gio: pt = " << pt << std::endl;
  std::cout << "Gio: dz = " << dz << ", sinphi = " << sinphi << ", dphi = " << dphi
            << ", dz/drphi = " << (dz / dphi / rho) << std::endl;
  std::cout << "Gio's fit P = " << mypq.first << "\n";
  std::cout << "Gio's fit Q = " << myq << "\n";

  return mypq;
}

CosmicGridTripletSeeder::TripletSeederEventState CosmicGridTripletSeeder::initEventState(const edm::EventSetup& c) {
  auto magfield = &c.getData(magfieldToken_);
  auto tracker = &c.getData(trackerToken_);
  auto cloner = dynamic_cast<TkTransientTrackingRecHitBuilder const&>(c.getData(ttrhBuilderToken_)).cloner();
  return TripletSeederEventState{SeedingGrid{8, -120., 120., 6, -120., 120., 8, -280., 280.},
                                 magfield,
                                 tracker,
                                 cloner,
                                 std::make_unique<KFUpdator>(),
                                 std::make_unique<PropagatorWithMaterial>(alongMomentum, 0.1057, magfield),
                                 std::make_unique<PropagatorWithMaterial>(oppositeToMomentum, 0.1057, magfield)};
}

void CosmicGridTripletSeeder::produce(edm::Event& e, const edm::EventSetup& c) {
  // setup the event state object
  TripletSeederEventState state = initEventState(c);

  // populate the seeding grid
  populateGrid(e,state);

  // sort the seeding grid
  state.grid.sort();

  // now form triplets
  std::vector<protoSeed> triplets;
  formTriplets(state, triplets);

  // book the trajectory seed collection
  auto output = std::make_unique<TrajectorySeedCollection>();

  // and fit the triplets
  fitTriplets(state, triplets, *output);

  // put our output on the store
  e.put(std::move(output));
}

bool CosmicGridTripletSeeder::populateGrid(const edm::Event& iEvent, CosmicGridTripletSeeder::TripletSeederEventState & state){
  const auto & otHitCollection = iEvent.get(otRecHitsToken_); 
  const auto & pixelHitCollection = iEvent.get(pixelRecHitsToken_); 
  const auto & vectorHits = iEvent.get(vectorHitsToken_); 

  std::set<const TrackingRecHit*> recHitsSeen{}; 
  std::vector<const VectorHit*> vhSeen{}; 
  /// step 1: Add the vector hits, and remember all raw hits associated to them 
  for (auto  ds : vectorHits){
      for (const VectorHit & vh : ds){
        state.grid.addHit(&vh); 
        vhSeen.push_back(&vh); 
      }
  }
  std::cout << " CGS: Done adding vector hits, now have "<<vhSeen.size()<<" unique hits"<< std::endl; 
  /// step 2: Add remaining OT hits, excluding those already on vector hits 
  for (auto  ds : otHitCollection){
      for (const Phase2TrackerRecHit1D & otHit : ds){
        bool unique = true; 
        for (const auto* vh : vhSeen){
            if (vh->sharesInput(&otHit,TrackingRecHit::some)){
                unique = false; 
                break; 
            }
        }
        if (!unique){
            std::cout << " skip a hit overlapping with VH" << std::endl; 
            continue; 
        }
        if (recHitsSeen.contains(&otHit)){
            std::cout << "this is already on a VH!"<< std::endl; 
            continue;
        }
        state.grid.addHit(&otHit); 
        recHitsSeen.insert(&otHit); 
      }
  }
  std::cout << " CGS: Done adding strip hits, now have "<<recHitsSeen.size() + vhSeen.size()<<" unique hits"<< std::endl; 
  /// step 3: Add pixel hits if desired  
  for (auto  ds : pixelHitCollection){
      for (const SiPixelRecHit & pix : ds){
        state.grid.addHit(&pix); 
      }
  }
  std::cout << " CGS: Done adding pixel hits, now have "<<recHitsSeen.size()<<" unique hits"<< std::endl; 

  // now sort all bins of the grid by ascending global y.
  state.grid.sort(); 

  return true;
}
/// using the seeding grid, form triplets
bool CosmicGridTripletSeeder::formTriplets(const CosmicGridTripletSeeder::TripletSeederEventState & state, std::vector<CosmicGridTripletSeeder::protoSeed> & found){

    return true; 
}
     /// get the triplets for one particular cell combination
void CosmicGridTripletSeeder::formTriplets(const std::vector<const BaseTrackerRecHit*> & topCands, 
                const std::vector<const BaseTrackerRecHit*> & centerCands, 
                const std::vector<const BaseTrackerRecHit*> & bottomCands, 
                const TripletSeederEventState & state,
                std::vector<protoSeed> & found, 
                std::unordered_multiset<const BaseTrackerRecHit*> & trackUsage){

  for (const BaseTrackerRecHit* top : topCands){
    if (trackUsage.count(top) > 2) continue; 
    for (const BaseTrackerRecHit* center: centerCands){
      if (trackUsage.count(center) > 2) continue; 
      for (const BaseTrackerRecHit* bottom: bottomCands){
        if (trackUsage.count(bottom) > 2) continue; 
        if (center == top || top == bottom || center == bottom) continue; 
        if (top->detUnit() == bottom->detUnit() || top->detUnit() == center->detUnit() || bottom->detUnit() == center->detUnit()) continue; 
        if (center->globalPosition().y() > top->globalPosition().y()) continue;
        if (bottom->globalPosition().y() > center->globalPosition().y()) continue;
        protoSeed ps (top, center, bottom);
        trackUsage.insert(top); 
        trackUsage.insert(center); 
        trackUsage.insert(bottom); 
        found.push_back(ps); 

      }
    }
  }
}
/// fit the triplets with kalman into trajectory seeds
bool CosmicGridTripletSeeder::fitTriplets(const CosmicGridTripletSeeder::TripletSeederEventState & state, const std::vector<CosmicGridTripletSeeder::protoSeed>& seeds, TrajectorySeedCollection & output){
    return true; 
}
/// fit of a single triplet into a trajectory seed 
bool CosmicGridTripletSeeder::fitTriplet(const CosmicGridTripletSeeder::TripletSeederEventState & state, const CosmicGridTripletSeeder::protoSeed& seed, TrajectorySeedCollection & output){
    return true; 
}
