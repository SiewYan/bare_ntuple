// -*- C++ -*-
//
// Package:    Bare_Ntuple/bare_ntuple
// Class:      bare_ntuple
// 
/**\class bare_ntuple bare_ntuple.cc Bare_Ntuple/bare_ntuple/plugins/bare_ntuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Siew Yan Hoh
//         Created:  Wed, 24 Jun 2015 09:16:10 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "TH1D.h"
//#include "TTree.h"
#include "TNtuple.h"



#define MAXVTX 200

// class declaration

class bare_ntuple : public edm::EDAnalyzer {
   public:
      explicit bare_ntuple(const edm::ParameterSet&);
      ~bare_ntuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::InputTag pfj_;//PFJetCorToken_;

  TTree *mtree;
  
  unsigned int mrun;
  unsigned int mlumi;
  unsigned int mevent;
  unsigned int mnvtx;
  unsigned int mnjets;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
bare_ntuple::bare_ntuple(const edm::ParameterSet& iConfig):
  pfj_(iConfig.getUntrackedParameter<edm::InputTag>("pfj"))
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;

  mtree = fs->make<TTree>( "ntuple" , " monojet_input " );


}


bare_ntuple::~bare_ntuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
bare_ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //using namespace reco;

   mrun   = iEvent.id().run();
   mlumi  = iEvent.luminosityBlock();
   mevent = iEvent.id().event();

   ////////////////Vertex////////////////
   reco::Vertex primaryvtx;
   Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel("offlineSlimmedPrimaryVertices",recVtxs);
   int  pvind=0;
   for(unsigned int ind=0;ind<recVtxs->size();ind++){
     if(!((*recVtxs)[ind].isFake())){
       //mPVx[pvind]      = (*recVtxs)[ind].x();
       //mPVy[pvind]      = (*recVtxs)[ind].y();
       //mPVz[pvind]      = (*recVtxs)[ind].z();
       //mPVchi2[pvind]   = (*recVtxs)[ind].chi2();
       //mPVndof[pvind]   = (*recVtxs)[ind].ndof();
       //mPVntracks[pvind]= (*recVtxs)[ind].tracksSize();
       if(pvind == 0) primaryvtx = (*recVtxs)[ind];
       pvind++;
       if(pvind>=MAXVTX) break;
     }
   }
   mnvtx = pvind;

   /////////////////RecoJet/////////////////////
   unsigned int count=0;
   Handle<pat::JetCollection> pfjetH;
   //ak4PFJetsCHS
   iEvent.getByLabel(pfj_, pfjetH);
   for ( pat::JetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ) {
     //jetpt->Fill(jet->pt());
     count++;
   }
   mnjets = count;
 
   mtree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
bare_ntuple::beginJob()
{

  mtree->Branch( "run" , &mrun ,"run/I" );
  mtree->Branch( "lumi" , &mlumi , "lumi/I" );
  mtree->Branch( "event" , &mevent , "event/I" );
  mtree->Branch( "nvtx" , &mnvtx , "nvtx/I" );
  mtree->Branch( "njets" , &mnjets , "njets/I" );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
bare_ntuple::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
bare_ntuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
bare_ntuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
bare_ntuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
bare_ntuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bare_ntuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bare_ntuple);
