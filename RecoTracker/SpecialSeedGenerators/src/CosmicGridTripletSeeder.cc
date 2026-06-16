
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TEllipse.h>
#include <cmath>
#include <memory>
#include <set>
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/SiStripDetId/interface/SiStripEnums.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/VectorHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "Geometry/CommonTopologies/interface/GeomDetEnumerators.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "RecoTracker/SpecialSeedGenerators/interface/CosmicGridTripletSeeder.h"

//
 

CosmicGridTripletSeeder::CosmicGridTripletSeeder(const edm::ParameterSet& iConfig)
    : vectorHitsToken_(mayConsume<VectorHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vectorHits"))),
      otRecHitsToken_(mayConsume<Phase2TrackerRecHit1DCollectionNew>(iConfig.getUntrackedParameter<edm::InputTag>("OTRecHits"))),
      matchedStripHitsToken_(mayConsume<SiStripMatchedRecHit2DCollection>(iConfig.getUntrackedParameter<edm::InputTag>("matchedStripHits"))),
      rPhiHitsToken_(mayConsume<SiStripRecHit2DCollection>(iConfig.getUntrackedParameter<edm::InputTag>("rPhiHits"))),
      // stereoHitsToken_(mayConsume<SiStripRecHit2DCollection>(iConfig.getUntrackedParameter<edm::InputTag>("stereoHits"))),
      pixelRecHitsToken_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("PixelRecHits"))),
      magfieldToken_(esConsumes(iConfig.getParameter<edm::ESInputTag>("MagneticFieldRecord"))),
      trackerToken_(esConsumes()),
      ttrhBuilderToken_(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("TTRHBuilder")))),
      m_nGridX(iConfig.getParameter<int>("nGridX")),
      m_nGridY(iConfig.getParameter<int>("nGridY")),
      m_nGridZ(iConfig.getParameter<int>("nGridZ")) {
  produces<TrajectorySeedCollection>();
  produces<edm::OwnVector<TrackingRecHit>>();

}

void CosmicGridTripletSeeder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("vectorHits", edm::InputTag("siPhase2VectorHits:accepted"));
  desc.addUntracked<edm::InputTag>("OTRecHits", edm::InputTag("Phase2TrackerRecHits"));
  desc.addUntracked<edm::InputTag>("matchedStripHits", edm::InputTag("siStripMatchedRecHits","matchedRecHit"));
  desc.addUntracked<edm::InputTag>("rPhiHits", edm::InputTag("siStripMatchedRecHits","rphiRecHit"));
  // desc.addUntracked<edm::InputTag>("stereoHits", edm::InputTag("siStripMatchedRecHits","matchedRecHit"));
  desc.addUntracked<edm::InputTag>("PixelRecHits", edm::InputTag("siPixelRecHits"));
  desc.add<std::string>("TTRHBuilder", "WithTrackAngle");
  desc.add<edm::ESInputTag>("MagneticFieldRecord", edm::ESInputTag("", ""));
  desc.add<int>("nGridX",1);
  desc.add<int>("nGridY",1);
  desc.add<int>("nGridZ",1);
  descriptions.addWithDefaultLabel(desc);
}

std::pair<GlobalVector, int> CosmicGridTripletSeeder::pqFromHelixFit(const GlobalPoint& inner,
                                                                     const GlobalPoint& middle,
                                                                     const GlobalPoint& outer,
                                                                     const MagneticField* magfield) {
  // std::cout << "DEBUG PZ =====" << std::endl;
  FastHelix helix(inner, middle, outer, magfield->nominalValue(), magfield);
  // GlobalVector gv = helix.stateAtVertex().momentum();  // status on inner hit
  // std::cout << "FastHelix P = " << gv << "\n";
  // std::cout << "FastHelix Q = " << helix.stateAtVertex().charge() << "\n";

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

  // std::cout << "Gio: pt = " << pt << std::endl;
  // std::cout << "Gio: dz = " << dz << ", sinphi = " << sinphi << ", dphi = " << dphi
  //           << ", dz/drphi = " << (dz / dphi / rho) << std::endl;
  // std::cout << "Gio's fit P = " << mypq.first << "\n";
  // std::cout << "Gio's fit Q = " << myq << "\n";

  return mypq;
}

CosmicGridTripletSeeder::TripletSeederEventState CosmicGridTripletSeeder::initEventState(const edm::EventSetup& c) {
  auto magfield = &c.getData(magfieldToken_);
  auto tracker = &c.getData(trackerToken_);
  auto cloner = dynamic_cast<TkTransientTrackingRecHitBuilder const&>(c.getData(ttrhBuilderToken_)).cloner();

  return TripletSeederEventState{SeedingGrid{m_nGridX, -120., 120., m_nGridY, -120., 120., m_nGridZ, -280., 280.},
                                 {},
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
  auto outtriplets = std::make_unique<edm::OwnVector<TrackingRecHit>>();

  // and fit the triplets
  fitTriplets(state, triplets, *output);

  for (auto & trip : triplets){
    outtriplets->push_back(trip.inner()->cloneHit()); 
    outtriplets->push_back(trip.middle()->cloneHit()); 
    outtriplets->push_back(trip.outer()->cloneHit()); 
  }

  // put our output on the store
  e.put(std::move(output));
  e.put(std::move(outtriplets));
}

bool CosmicGridTripletSeeder::populateGrid(const edm::Event& iEvent, CosmicGridTripletSeeder::TripletSeederEventState & state){
  edm::Handle<VectorHitCollection> vectorHits; 
  edm::Handle<Phase2TrackerRecHit1DCollectionNew> otHitCollection; 
  edm::Handle<SiStripMatchedRecHit2DCollection> matchedStripHits; 
  edm::Handle<SiStripRecHit2DCollection> rPhiStripHits; 

  bool hasOT = iEvent.getByToken(otRecHitsToken_, otHitCollection );
  bool hasVec = iEvent.getByToken(vectorHitsToken_, vectorHits );
  bool hasMatchedStrips = iEvent.getByToken(matchedStripHitsToken_, matchedStripHits );
  bool hasRphiStrips = iEvent.getByToken(rPhiHitsToken_, rPhiStripHits );

  const SiPixelRecHitCollection & pixelHitCollection = iEvent.get(pixelRecHitsToken_);

  std::set<const TrackingRecHit*> recHitsSeen{}; 
  std::vector<const VectorHit*> vhSeen{}; 
  /// step 1: Add the vector hits, and remember all raw hits associated to them 
  if (hasVec){
    for (auto  ds : *vectorHits){
        for (const VectorHit & vh : ds){
          // endcap vector hits have invalid directional information, 
          // hence we skip these. 
          if (state.tracker->idToDet(vh.geographicalId())->subDetector() == GeomDetEnumerators::SubDetector::P2OTEC){
            continue; 
          }
          state.grid.addHit(&vh); 
          vhSeen.push_back(&vh); 
        }
    }
  }
  if (hasOT){
    /// step 2: Add remaining OT hits, excluding those already on vector hits 
    for (auto  ds : *otHitCollection){
        for (const auto & otHit : ds){
          bool unique = true; 
          for (const auto* vh : vhSeen){
              if (vh->sharesInput(&otHit,TrackingRecHit::some)){
                  unique = false; 
                  state.vhConstituents[vh].push_back(&otHit);
                  break; 
              }
          }
          if (!unique){
              continue; 
          }
          state.grid.addHit(&otHit); 
          recHitsSeen.insert(&otHit); 
        }
    }
  }
  /// phase-1 matched strip hits
  if (hasMatchedStrips){
    for (auto  ds : *matchedStripHits){
      for (const auto & stripHit : ds){
        bool unique = true; 
        for (const auto* vh : vhSeen){
            if (vh->sharesInput(&stripHit,TrackingRecHit::some)){
                unique = false; 
                state.vhConstituents[vh].push_back(&stripHit);
                break; 
            }
        }
        if (!unique){
            continue; 
        }
        state.grid.addHit(&stripHit); 
        recHitsSeen.insert(&stripHit); 
      }
    }
  } 

  // phase-1 rphi strip hits 
  if (hasRphiStrips){
    for (auto  ds : *rPhiStripHits){
      for (const auto & stripHit : ds){
        bool unique = true; 
        for (const auto* vh : vhSeen){
            if (vh->sharesInput(&stripHit,TrackingRecHit::some)){
                unique = false; 
                state.vhConstituents[vh].push_back(&stripHit);
                break; 
            }
        }
        if (!unique){
            continue; 
        }
        state.grid.addHit(&stripHit); 
        recHitsSeen.insert(&stripHit); 
      }
    }
  }

  /// step 3: Add pixel hits if desired  
  for (auto  ds : pixelHitCollection){
      for (const auto & pix : ds){
        state.grid.addHit(&pix); 
      }
  }

  // now sort all bins of the grid by ascending global y.
  state.grid.sort(); 

  return true;
}

/// using the seeding grid, form triplets
bool CosmicGridTripletSeeder::formTriplets(const CosmicGridTripletSeeder::TripletSeederEventState & state, std::vector<CosmicGridTripletSeeder::protoSeed> & found){
    std::unordered_multiset<const BaseTrackerRecHit*> trackUsage; 
    // loop downwards over the starting y bin, top to bottom 
    for (int yBin = state.grid.nBinsY()-1; yBin >=0 ; --yBin){
      // loop rectangularily over the x-z grid 
      for (int xBin = 0; xBin < state.grid.nBinsX(); ++xBin){
        for (int zBin = 0; zBin < state.grid.nBinsZ(); ++zBin){
          const auto & topCandidates = state.grid.getHits(xBin,yBin,zBin); 
          if (topCandidates.empty()) continue;         
          formTriplets(xBin,yBin,zBin,state,found,trackUsage); 
        }
      }
    }
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
  // if (yBin  0)
    // return;

  std::vector<const BaseTrackerRecHit *> hitCands; 
  hitCands.insert(hitCands.end(),topHits.begin(), topHits.end() ); 

  for (int dxCenter = -1; dxCenter < 2; ++dxCenter) {
    if (xBin + dxCenter < 0 || xBin + dxCenter >= grid.nBinsX())
      continue;
    for (int dzCenter = -1; dzCenter < 2; ++dzCenter) {
      if (zBin + dzCenter < 0 || zBin + dzCenter >= grid.nBinsZ())
        continue;
      // formTriplets(topHits, topHits, topHits, state, found, trackUsage);
      if (yBin > 0){
        const auto& middleHits = grid.getHits(xBin + dxCenter, yBin - 1, zBin + dzCenter);
        hitCands.insert(hitCands.end(),middleHits.begin(), middleHits.end() ); 
        // // formTriplets(topHits, topHits, middleHits, state, found, trackUsage);
        // // formTriplets(topHits, middleHits, middleHits, state, found, trackUsage);
        // for (int dxBottom = -1; dxBottom < 2; ++dxBottom) {
        //   if (xBin + dxCenter + dxBottom < 0 || xBin + dxCenter + dxBottom >= grid.nBinsX())
        //     continue;
        //   if ((dxCenter < 0 && dxBottom > 0) || (dxCenter > 0 && dxBottom < 0))
        //     continue;
        //   for (int dzBottom = -1; dzBottom < 2; ++dzBottom) {
        //     if (zBin + dzCenter + dzBottom < 0 || zBin + dzCenter + dzBottom >= grid.nBinsZ())
        //       continue;
        //     if ((dzCenter < 0 && dzBottom > 0) || (dzCenter > 0 && dzBottom < 0))
        //       continue;
        //     if (yBin > 1){
        //       const auto& bottomHits = grid.getHits(xBin + dxCenter + dxBottom, yBin - 2, zBin + dzCenter + dzBottom);
        //       formTriplets(topHits, middleHits, bottomHits, state, found, trackUsage);
        //     }
        //   }
        // }
      }
    }
  }
  std::sort(hitCands.begin(),hitCands.end(),[](const BaseTrackerRecHit* h1, const BaseTrackerRecHit* h2){return h1->globalPosition().y() < h2->globalPosition().y(); });
  formTriplets(hitCands, hitCands, hitCands, state, found, trackUsage);
}

     /// get the triplets for one particular cell combination
void CosmicGridTripletSeeder::formTriplets(const std::vector<const BaseTrackerRecHit*> & topCands, 
                const std::vector<const BaseTrackerRecHit*> & centerCands, 
                const std::vector<const BaseTrackerRecHit*> & bottomCands, 
                const TripletSeederEventState & state,
                std::vector<protoSeed> & found, 
                std::unordered_multiset<const BaseTrackerRecHit*> & trackUsage){

  for (const BaseTrackerRecHit* top : topCands){
    if (trackUsage.count(top) > 12) continue; 
    const auto & gTop = top->globalPosition(); 
    for (const BaseTrackerRecHit* center: centerCands){
      const auto & gCenter = center->globalPosition(); 
      if (trackUsage.count(center) > 12) continue; 
      for (const BaseTrackerRecHit* bottom: bottomCands){
        const auto & gBottom = bottom->globalPosition(); 
        if (trackUsage.count(bottom) > 12) continue; 
        if (center == top || top == bottom || center == bottom) continue; 
        if (top->sameDetModule(*bottom) || top->sameDetModule(*center) || center->sameDetModule(*bottom)) continue; 
        if (gCenter.y() > gTop.y()) continue;
        if (gBottom.y() > gCenter.y()) continue;
        // veto "zigzag" in z 
        if ((gTop.z() - gCenter.z()) * (gCenter.z() - gBottom.z())  < 0 && std::abs(gTop.z() - gBottom.z() > 1.0)) continue;
        // veto "zigzag" in x 
        if ((gTop.x() - gCenter.x()) * (gCenter.x() - gBottom.x())  < 0 && std::abs(gTop.x() - gBottom.x() > 1.0)) continue;
        const VectorHit* cenVH = dynamic_cast<const VectorHit*>(center); 
        if (cenVH){
          double dydx = 0.5 * ((gTop.y() - gCenter.y())/(gTop.x() - gCenter.x()) + (gCenter.y() - gBottom.y())/(gCenter.x() - gBottom.x())); 
          double dydx_vh = cenVH->globalDirectionVH().y() / cenVH->globalDirectionVH().x();
          std::cout << " dydx from trip "<<dydx<<" and from vh "<<dydx_vh<<std::endl; 
        }
        // try this:
        // double dyz_top =  (gTop.z() - gCenter.z()) / (gTop.y() - gCenter.y());
        // double dyz_bot =  (gCenter.z() - gBottom.z()) / (gCenter.y() - gBottom.y()); 
        // double dyz_trip =  (gTop.z() - gBottom.z()) / (gTop.y() - gBottom.y()); 
        // double dxy_top =  (gTop.x() - gCenter.x()) / (gTop.y() - gCenter.y());
        // double dxy_bot =  (gCenter.x() - gBottom.x()) / (gCenter.y() - gBottom.y()); 
        // double dxy_trip =  (gTop.x() - gBottom.x()) / (gTop.y() - gBottom.y()); 
        // std::cout << " dyz_top "<<dyz_top<<" "<<" dyz_bottom " <<dyz_bot<< " dyz_trip "<<dyz_trip<< " dTrip " <<5.0 / std::abs((gTop.y() - gBottom.y()))<< std::endl; 
        // std::cout << " dxy_top "<<dxy_top<<" "<<" dxy_bottom " <<dxy_bot<< " dxy_trip "<<dxy_trip<< std::endl; 
        // if (std::abs(dxy_top - dxy_bot) > 0.5 * std::abs(dxy_trip)) continue;



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
    // std::cout << " built "<<output.size()<<" seeds from "<<triplets.size()<<" triplets "<< std::endl; 
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

    if ((outer.y() - inner.y()) * outer.y() < 0) {
      std::swap(inner, outer);
      trip = OrderedHitTriplet(trip.outer(), trip.middle(), trip.inner());
    }

    // First use FastHelix out of the box
    std::pair<GlobalVector, int> pq = pqFromHelixFit(inner, middle, outer, state.magfield);
    GlobalVector gv = pq.first;
    float ch = pq.second;
    float Mom = sqrt(gv.x() * gv.x() + gv.y() * gv.y() + gv.z() * gv.z());

    if (Mom > 1000000 || edm::isNotFinite(Mom)) {
        // std::cout << "Processing triplet " << ": fail for momentum." << std::endl;
      return false;
    }

    if (gv.perp() < 0.5) {
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
    }

    if ((gv.z() * (outer.z() - inner.z()) > 0) && (fabs(outer.z() - inner.z()) > 5) && (fabs(gv.z()) > .01)) {
    }


    edm::OwnVector<TrackingRecHit> hits;
    std::vector<const BaseTrackerRecHit*> seedHits;
    for (const BaseTrackerRecHit* hit : {trip.outer(), trip.middle(), trip.inner()}){
      const VectorHit* vh = dynamic_cast<const VectorHit*>(hit);
      if (vh){
        auto found = state.vhConstituents.find(vh);
        if (found != state.vhConstituents.end()){
          for (auto & component : found->second){
            seedHits.push_back(component);
          }
        }
      }
      else{
        seedHits.push_back(hit);
        // std::cout <<"         "<< seedHits.back()->globalPosition().y()<< std::endl; 
      }
    }
    std::sort(seedHits.begin(),seedHits.end(),[&](const BaseTrackerRecHit* h1, const BaseTrackerRecHit* h2){
      if (outer.y() > 0 ){
        return h1->globalPosition().y() > h2->globalPosition().y();
      }
      else{
        return h1->globalPosition().y() < h2->globalPosition().y();
      }
    }); 
    outer = seedHits.front()->globalPosition(); 
    GlobalTrajectoryParameters Gtp(outer, gv, int(ch), state.magfield);
    FreeTrajectoryState CosmicSeed(Gtp, CurvilinearTrajectoryError(AlgebraicSymMatrix55(AlgebraicMatrixID())));
    CosmicSeed.rescaleError(100);
    // std::cout << "Processing triplet " << ". start from " << std::endl;
    // std::cout << "    X  = " << outer << ", P = " << gv << std::endl;
    // std::cout << "    Cartesian error (X,P) = \n" << CosmicSeed.cartesianError().matrix() << std::endl;
    // std::cout << " start prop from "<<outer << std::endl ;
      
    TSOS propagated, updated;
    bool fail = false;
    for (size_t ih = 0; ih < seedHits.size(); ++ih) {
      // if ((ih == 2) && seedOnMiddle_) {
      //   if (seedVerbosity_ > 2)
      //     std::cout << "Sinnerping at middle hit, as requested." << std::endl;
      //   break;
      // }
      // std::cout <<" about to prop from "<<CosmicSeed.position()<<" to "<<state.tracker->idToDet((*seedHits[ih]).geographicalId())->surface().position()<<" with momentum "<<CosmicSeed.momentum()<<std::endl;
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
          // std::cout << "Processing triplet "  << ", hit " << ih << ": ok"<<std::endl;
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
          // std::cout << "Processing triplet "  << ", hit " << ih << ": ok"<<std::endl; 
      }
    }
    if (!fail && updated.isValid() && (updated.globalMomentum().perp() < 2.5)) {
        // std::cout << "Processing triplet "  << ": failed for final pt " << updated.globalMomentum().perp() << " < 2.5"
                  // << std::endl;
      fail = true;
    }
    if (!fail && updated.isValid() && (updated.globalMomentum().mag() < 2.5)) {
        // std::cout << "Processing triplet "  << ": failed for final p " << updated.globalMomentum().perp() << " < 2.5"
                  //  << std::endl;
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
        // updated.rescaleError(100);
      // }
    //   if (true) {
    //   std::cout << "Processed  triplet "  << ": success (saved as #" << output.size() << ") : " << inner << " + "
    //             << middle << " + " << outer << std::endl;
    //   std::cout << "    pt = " << updated.globalMomentum().perp() << "    eta = " << updated.globalMomentum().eta()
    //             << "    phi = " << updated.globalMomentum().phi() << "    ch = " << updated.charge() << std::endl;
    //   if (true) {
    //     std::cout << "    State:" << updated;
    //   } else {
    //     std::cout << "    X  = " << updated.globalPosition() << ", P = " << updated.globalMomentum() << std::endl;
    //   }
    //   std::cout << "    Cartesian error (X,P) = \n" << updated.cartesianError().matrix() << std::endl;
    // }   

    PTrajectoryStateOnDet const &PTraj = trajectoryStateTransform::persistentState(
        updated, hits.back().geographicalId().rawId());
    // output.push_back(TrajectorySeed(PTraj, hits, ((outer.y() - inner.y() > 0) ? alongMomentum : oppositeToMomentum)));
    output.push_back(TrajectorySeed(PTraj, hits, ((outer.y() - inner.y() > 0) ? alongMomentum : oppositeToMomentum)));

    // std::cout <<"  Wrote seed at "<<state.tracker->idToDet(output.back().startingState().detId())->toGlobal(output.back().startingState().parameters().position())<<" in direction "<<(output.back().direction() == 0 ? "opposite" : "along")<<" with momentum direction "<<state.tracker->idToDet(output.back().startingState().detId())->toGlobal(output.back().startingState().parameters().momentum())<<std::endl;
    if (output.size() > size_t(500)) {
      output.clear();
      edm::LogError("TooManySeeds") << "Found too many seeds, bailing out.\n";
      return false;
    }
    return true;

}
