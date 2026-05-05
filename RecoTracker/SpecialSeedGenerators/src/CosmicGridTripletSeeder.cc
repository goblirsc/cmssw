
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
    std::unordered_multiset<const BaseTrackerRecHit*> trackUsage; 
    // loop downwards over the starting y bin, top to bottom 
    for (int yBin = state.grid.nBinsY()-1; yBin >0 ; --yBin){
      // loop rectangularily over the x-z grid 
      for (int xBin = 0; xBin < state.grid.nBinsX(); ++xBin){
        for (int zBin = 0; zBin < state.grid.nBinsZ(); ++zBin){
          const auto & topCandidates = state.grid.getHits(xBin,yBin,zBin); 
          if (topCandidates.empty()) continue;         
          formTriplets(xBin,yBin,zBin,state,found,trackUsage); 
        }
      }
    }
    std::cout << " found "<<found.size()<<" new triplets "<< std::endl;
    return true; 
}

void CosmicGridTripletSeeder::formTriplets(int xBin,
                                           int yBin,
                                           int zBin,
                                           const TripletSeederEventState& state,
                                           std::vector<protoSeed>& found,
                                           std::unordered_multiset<const BaseTrackerRecHit*>& trackUsage) {
  const SeedingGrid& grid = state.grid;

  const auto& topHits = grid.getHits(xBin, yBin, zBin);
  // top to bottom navigation: Check a 3x3 grid for the next hit in the chain
  if (yBin == 0)
    return;

  for (int dxCenter = -1; dxCenter < 2; ++dxCenter) {
    if (xBin + dxCenter < 0 || xBin + dxCenter >= grid.nBinsX())
      continue;
    for (int dzCenter = -1; dzCenter < 2; ++dzCenter) {
      if (zBin + dzCenter < 0 || zBin + dzCenter >= grid.nBinsZ())
        continue;
      const auto& middleHits = grid.getHits(xBin + dxCenter, yBin - 1, zBin + dzCenter);
      formTriplets(topHits, topHits, middleHits, state, found, trackUsage);
      formTriplets(topHits, middleHits, middleHits, state, found, trackUsage);

      if (yBin == 1)
        continue;
      for (int dxBottom = -1; dxBottom < 2; ++dxBottom) {
        if (xBin + dxCenter + dxBottom < 0 || xBin + dxCenter + dxBottom >= grid.nBinsX())
          continue;
        if ((dxCenter < 0 && dxBottom > 0) || (dxCenter > 0 && dxBottom < 0))
          continue;
        for (int dzBottom = -1; dzBottom < 2; ++dzBottom) {
          if (zBin + dzCenter + dzBottom < 0 || zBin + dzCenter + dzBottom >= grid.nBinsZ())
            continue;
          if ((dzCenter < 0 && dzBottom > 0) || (dzCenter > 0 && dzBottom < 0))
            continue;
          const auto& bottomHits = grid.getHits(xBin + dxCenter + dxBottom, yBin - 2, zBin + dzCenter + dzBottom);
          formTriplets(topHits, middleHits, bottomHits, state, found, trackUsage);
        }
      }
    }
  }
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
        if ((top->globalPosition().z() - center->globalPosition().z()) * (center->globalPosition().z() - bottom->globalPosition().z())  < 0) continue;
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
bool CosmicGridTripletSeeder::fitTriplets(const CosmicGridTripletSeeder::TripletSeederEventState & state, const std::vector<CosmicGridTripletSeeder::protoSeed>& triplets, TrajectorySeedCollection & output){
    for (const protoSeed & ps : triplets ){
      fitTriplet(state, ps, output); 
    }
    std::cout << " built "<<output.size()<<" seeds from "<<triplets.size()<<" triplets "<< std::endl; 
    return true; 
}
/// fit of a single triplet into a trajectory seed 
bool CosmicGridTripletSeeder::fitTriplet(const CosmicGridTripletSeeder::TripletSeederEventState & state, const CosmicGridTripletSeeder::protoSeed& triplet, TrajectorySeedCollection & output){
  typedef TrajectoryStateOnSurface TSOS;


    OrderedHitTriplet trip = triplet;  


    GlobalPoint inner =
        state.tracker->idToDet((*(trip.inner())).geographicalId())->surface().toGlobal((*(trip.inner())).localPosition());

    GlobalPoint middle =
        state.tracker->idToDet((*(trip.middle())).geographicalId())->surface().toGlobal((*(trip.middle())).localPosition());

    GlobalPoint outer =
        state.tracker->idToDet((*(trip.outer())).geographicalId())->surface().toGlobal((*(trip.outer())).localPosition());

    // std::cout << "Processing triplet " << ": " << inner << " + " << middle << " + " << outer << std::endl;

    if ((outer.y() - inner.y()) * outer.y() < 0) {
      std::swap(inner, outer);
      trip = OrderedHitTriplet(trip.outer(), trip.middle(), trip.inner());

        // std::cout << "The seed was going away from CMS! swapped in <-> out" << std::endl;
        // std::cout << "Processing swapped triplet  : " << inner << " + " << middle << " + " << outer
                  // << std::endl;
    }

    // First use FastHelix out of the box
    std::pair<GlobalVector, int> pq = pqFromHelixFit(inner, middle, outer, state.magfield);
    GlobalVector gv = pq.first;
    float ch = pq.second;
    float Mom = sqrt(gv.x() * gv.x() + gv.y() * gv.y() + gv.z() * gv.z());

    if (Mom > 10000 || edm::isNotFinite(Mom)) {
        // std::cout << "Processing triplet " << ": fail for momentum." << std::endl;
      return false;
    }

    if (gv.perp() < 2.5) {
        // std::cout << "Processing triplet " << ": fail for pt = " << gv.perp() << " < ptMin = 2.5"
                  // << std::endl;
      return false;
    }

    const Propagator *propagator = nullptr;
    if ((outer.y() - inner.y()) > 0) {
        // std::cout << "Processing triplet " << ":  downgoing." << std::endl;
      propagator = state.thePropagatorAl.get();
    } else {
      gv = -1 * gv;
      ch = -1. * ch;
      propagator = state.thePropagatorOp.get();
        // std::cout << "Processing triplet " << ":  upgoing." << std::endl;
    }

    if ((gv.z() * (outer.z() - inner.z()) > 0) && (fabs(outer.z() - inner.z()) > 5) && (fabs(gv.z()) > .01)) {
      // std::cout << "ORRORE: outer.z()-inner.z() = " << (outer.z() - inner.z()) << ", gv.z() = " << gv.z()
      //           << std::endl;
    }

    GlobalTrajectoryParameters Gtp(outer, gv, int(ch), state.magfield);
    FreeTrajectoryState CosmicSeed(Gtp, CurvilinearTrajectoryError(AlgebraicSymMatrix55(AlgebraicMatrixID())));
    CosmicSeed.rescaleError(100);
    // std::cout << "Processing triplet " << ". start from " << std::endl;
    // std::cout << "    X  = " << outer << ", P = " << gv << std::endl;
    // std::cout << "    Cartesian error (X,P) = \n" << CosmicSeed.cartesianError().matrix() << std::endl;

    edm::OwnVector<TrackingRecHit> hits;
    OrderedHitTriplet seedHits(trip.outer(), trip.middle(), trip.inner());
    TSOS propagated, updated;
    bool fail = false;
    for (size_t ih = 0; ih < 3; ++ih) {
      // if ((ih == 2) && seedOnMiddle_) {
      //   if (seedVerbosity_ > 2)
      //     std::cout << "Stopping at middle hit, as requested." << std::endl;
      //   break;
      // }
      if (ih == 0) {
        propagated = propagator->propagate(CosmicSeed, state.tracker->idToDet((*seedHits[ih]).geographicalId())->surface());
      } else {
        propagated = propagator->propagate(updated, state.tracker->idToDet((*seedHits[ih]).geographicalId())->surface());
      }
      if (!propagated.isValid()) {
        // std::cout << "Processing triplet "  << ", hit " << ih << ": failed propagation." << std::endl;
        fail = true;
        break;
      } else {
          // std::cout << "Processing triplet "  << ", hit " << ih << ": propagated state = " << propagated;
      }
      SeedingHitSet::ConstRecHitPointer tthp = seedHits[ih];
      auto newtth = static_cast<SeedingHitSet::RecHitPointer>(state.cloner(*tthp, propagated));
      updated = state.theUpdator->update(propagated, *newtth);
      hits.push_back(newtth);
      if (!updated.isValid()) {
          // std::cout << "Processing triplet "  << ", hit " << ih << ": failed update." << std::endl;
        fail = true;
        break;
      } else {
          // std::cout << "Processing triplet "  << ", hit " << ih << ": updated state = " << updated;
      }
    }
    if (!fail && updated.isValid() && (updated.globalMomentum().perp() < 2.5)) {
        // std::cout << "Processing triplet "  << ": failed for final pt " << updated.globalMomentum().perp() << " < 2.5"
        //           << std::endl;
      fail = true;
    }
    if (!fail && updated.isValid() && (updated.globalMomentum().mag() < 2.5)) {
        // std::cout << "Processing triplet "  << ": failed for final p " << updated.globalMomentum().perp() << " < 2.5"
        //            << std::endl;
      fail = true;
    }
    if (fail) return false; 
    // if (!fail) {
        // if (seedVerbosity_ > 2) {
        //   std::cout << "Processing triplet "  << ", rescale error by " << rescaleError_
        //             << ": state BEFORE rescaling " << updated;
        //   std::cout << "    Cartesian error (X,P) before rescaling= \n"
        //             << updated.cartesianError().matrix() << std::endl;
        // }
        // updated.rescaleError(rescaleError_);
      // }
      // if (seedVerbosity_ > 0) {
      // std::cout << "Processed  triplet "  << ": success (saved as #" << out.size() << ") : " << inner << " + "
                // << middle << " + " << outer << std::endl;
      // std::cout << "    pt = " << updated.globalMomentum().perp() << "    eta = " << updated.globalMomentum().eta()
                // << "    phi = " << updated.globalMomentum().phi() << "    ch = " << updated.charge() << std::endl;
      // if (seedVerbosity_ > 1) {
        // std::cout << "    State:" << updated;
      // } else {
      //   std::cout << "    X  = " << updated.globalPosition() << ", P = " << updated.globalMomentum() << std::endl;
      // }
      // std::cout << "    Cartesian error (X,P) = \n" << updated.cartesianError().matrix() << std::endl;
    // }

    PTrajectoryStateOnDet const &PTraj = trajectoryStateTransform::persistentState(
        // updated, (*(seedOnMiddle_ ? trip.middle() : trip.inner())).geographicalId().rawId());
        updated, (*trip.inner()).geographicalId().rawId());
    output.push_back(TrajectorySeed(PTraj, hits, ((outer.y() - inner.y() > 0) ? alongMomentum : oppositeToMomentum)));
    if (output.size() > size_t(50)) {
      output.clear();
      edm::LogError("TooManySeeds") << "Found too many seeds, bailing out.\n";
      return false;
    }
    return true;

}
