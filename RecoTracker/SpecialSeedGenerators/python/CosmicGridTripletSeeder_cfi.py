import FWCore.ParameterSet.Config as cms

from RecoTracker.SpecialSeedGenerators.cosmicGridTripletSeeder_cfi import cosmicGridTripletSeeder 
cosmicGridTripletSeeds = cosmicGridTripletSeeder.clone(
   vectorHits    = cms.untracked.InputTag("siPhase2VectorHits:accepted"),
   OTRecHits = cms.untracked.InputTag("siPhase2RecHits"),
   PixelRecHits = cms.untracked.InputTag("siPixelRecHits"),
   TTRHBuilder = cms.string('WithTrackAngle'),
   MagneticFieldRecord = cms.ESInputTag('', ''),
)

