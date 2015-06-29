import FWCore.ParameterSet.Config as cms


process = cms.Process("DEMO")
process.load("FWCore.MessageService.MessageLogger_cfi")


era ='PHYS14_V4_MC'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "PHYS14_25_V2::All"
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                                   connect = cms.string('sqlite_file:/afs/cern.ch/work/s/shoh/analysis/monojet/13TeV_phys14/CMSSW_7_4_0/src/bare_ntuple/bare_ntuple/python/PHYS14_V4_MC.db'),
                                                           toGet =  cms.VPSet(
               cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                                label= cms.untracked.string("AK4PF")),
               cms.PSet(record = cms.string("JetCorrectionsRecord"),
                                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                                label= cms.untracked.string("AK4PFchs")),
                )
                                                           )
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

process.load("RecoJets.Configuration.RecoJets_EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         #'root://cmsxrootd.fnal.gov//store/mc/TP2023SHCALDR/QCD_Pt-15to3000_Tune4C_14TeV_pythia8/GEN-SIM-RECO/SHCALJan23_NoPU_PH2_1K_FB_V6-v1/00000/001FF988-5CA3-E411-BE7A-00261894388A.root'          
#'/store/mc/Phys14DR/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02FBEDF5-8373-E411-AEC8-E0CB4E553665.root'
    )
)



from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.Configuration.RecoGenJets_cff import *
from bare_ntuple.bare_ntuple.JetCorrection_cff import *
process.load('bare_ntuple.bare_ntuple.JetCorrection_cff')

ak4GenJets = ak5GenJets.clone ( rParam=0.4 )
ak4GenJetsNoMuNoNu = ak4GenJets.clone ( src = 'genParticlesForJetsNoMuNoNu' )
#This is the collection I actually use
ak4GenJetsNoNu = ak4GenJets.clone ( src = 'genParticlesForJetsNoNu' )

setattr(process,'genParticlesForJetsNoNu',genParticlesForJetsNoNu)
sequence = cms.Sequence(genParticlesForJetsNoNu)
setattr(process,'ak4GenJetsNoNu',ak4GenJetsNoNu) #the second argument is the input tag!
setattr(process,'ak4PFchsJetsL1FastL2L3',ak4PFchsJetsL1FastL2L3)
sequence = cms.Sequence(sequence * ak4GenJetsNoNu * ak4PFchsJetsL1FastL2L3)

process.demo = cms.EDAnalyzer('demo',srcRho = cms.InputTag('kt6PFJets', 'rho'))

sequence = cms.Sequence(sequence * process.demo)
setattr(process,'demoAK4Sequence',sequence)
path = cms.Path( sequence )
setattr(process,'demoAK4Path',path)

