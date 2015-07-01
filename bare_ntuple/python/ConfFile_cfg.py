import FWCore.ParameterSet.Config as cms

process = cms.Process("BARE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "PHYS14_25_V2"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:myfile.root' if it is local
'file:/afs/cern.ch/user/s/shoh/eos/cms/store/relval/CMSSW_7_4_1/RelValADDMonoJet_d3MD3_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/80CF5456-B9EC-E411-93DA-002618FDA248.root'
    )
)

process.bare = cms.EDAnalyzer('bare_ntuple',
                              offlineprimaryvertices = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                              pfj = cms.untracked.InputTag('slimmedJets'), #ak4PFJetsCHS
                              muon = cms.untracked.InputTag('slimmedMuons'), #pt>5 GeV pass PFMuons ID
                              met = cms.untracked.InputTag('slimmedMETs'), #type 1 PFMET
                              electron = cms.untracked.InputTag('slimmedElectrons'), #all gedGsfElectrons
                              tau = cms.untracked.InputTag('slimmedTaus'), # pt>18 GeV passing basic decayModeFinding id 
                              photon = cms.untracked.InputTag('slimmedPhotons') # gedPhoton pt > 14 GeV 
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService", fileName = cms.string('bare_ntuple.root') )

process.p = cms.Path(process.bare)

#process.p=cms.Path(process.bare*process.dump)

#process.p=cms.Path(process.dump)
