import FWCore.ParameterSet.Config as cms

process = cms.Process("BARE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:myfile.root' if it is local
        #'/store/mc/Phys14DR/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02FBEDF5-8373-E411-AEC8-E0CB4E553665.root'
#'file:/afs/cern.ch/user/s/shoh/eos/cms/store/relval/CMSSW_7_4_1/RelValADDMonoJet_d3MD3_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/80CF5456-B9EC-E411-93DA-002618FDA248.root'
'/store/mc/Phys14DR/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02FBEDF5-8373-E411-AEC8-E0CB4E553665.root'
    )
)

process.bare = cms.EDAnalyzer('bare_ntuple',
                              pfj = cms.untracked.InputTag('ak4PFJetsCHS')
                              #PFJetCorToken_ = cms.untracked.InputTag('ak4PFchsL1FastL2L3')
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('bare_ntuple.root') )
path = cms.Path(process.bare)
