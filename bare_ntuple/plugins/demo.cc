// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include <TMath.h>

#include <memory>
#include <map>
#include <math.h>
#include <vector>

//
// class declaration
//

class demo : public edm::EDAnalyzer {
public:
  explicit demo(const edm::ParameterSet&);
  ~demo();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  const JetCorrector* jetCorrector_;

  edm::InputTag srcRho_;
  //edm::InputTag jetCorrector_;

  unsigned int nGenJet;
  unsigned int nPFJet;
  unsigned int nPFCHSJet;
  unsigned int nPFCHSL1L2L3Jet;

  TFile* rootFile_;

};


//
// constructors and destructor
//
demo::demo(const edm::ParameterSet& iConfig):
srcRho_(iConfig.getParameter<edm::InputTag>("srcRho"))
{
  //now do what ever initialization is needed
  rootFile_ = new TFile("ntuple.root","RECREATE");
  rootFile_->cd();

  jetCorrector_=0;
}


demo::~demo()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
demo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;  
  
  jetCorrector_ = JetCorrector::getJetCorrector("ak4PFchsL1FastL2L3", iSetup);

  /////////// gen Jets //////////////
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(InputTag("ak4GenJetsNoNu"), genJets);
  //iConfig.getParameter<edm::InputTag>("GenJetAlgorithm" )
  nGenJet=0;
  for (GenJetCollection::const_iterator genjets = genJets->begin(); genjets != genJets->end(); genjets++) {
    nGenJet++;  
  }
  
  /////////// ak4pf jets //////////////
  Handle<reco::PFJetCollection> pfjetcoll;
  iEvent.getByLabel(InputTag("ak4PFJets"), pfjetcoll);
  const reco::PFJetCollection *pfjets = pfjetcoll.product();
  reco::PFJetCollection::const_iterator pfjetclus = pfjets->begin();

  nPFJet=0;
  for(pfjetclus = pfjets->begin(); pfjetclus!= pfjets->end(); ++pfjetclus){
    nPFJet++;  
  }
  
  /////////// ak4pfchs jets //////////////
  Handle<reco::PFJetCollection> pfjetchscoll;
  iEvent.getByLabel(InputTag("ak4PFJetsCHS"), pfjetchscoll);
  const reco::PFJetCollection *pfchsjets = pfjetchscoll.product();
  reco::PFJetCollection::const_iterator pfjetchsclus = pfchsjets->begin();

  nPFCHSJet=0;
  for(pfjetchsclus = pfchsjets->begin(); pfjetchsclus!= pfchsjets->end(); ++pfjetchsclus){
    nPFCHSJet++;
  }  
  
  /////////// ak4pfchs L1L2L3 corrected jets //////////////
  Handle<reco::PFJetCollection> pfjetL1L2L3chscoll;
  iEvent.getByLabel(InputTag("ak4PFchsJetsL1FastL2L3"), pfjetL1L2L3chscoll);
  const reco::PFJetCollection *pfchsL1L2L3jets = pfjetL1L2L3chscoll.product();
  reco::PFJetCollection::const_iterator pfjetL1L2L3chsclus = pfchsL1L2L3jets->begin();

  nPFCHSL1L2L3Jet=0;
  for(pfjetL1L2L3chsclus = pfchsL1L2L3jets->begin(); pfjetL1L2L3chsclus!= pfchsL1L2L3jets->end(); ++pfjetL1L2L3chsclus){
    nPFCHSL1L2L3Jet++;
  } 
}


// ------------ method called once each job just before starting event loop  ------------
void 
demo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
demo::endJob() 
{
  // go to *OUR* root file and store histograms
  rootFile_->cd();
  rootFile_->Write();
  rootFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
demo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(demo);
