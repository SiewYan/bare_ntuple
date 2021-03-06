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

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TH1D.h"

//
// class declaration
//

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
  edm::InputTag pfj_;

  TH1D *jetpt;
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

  jetpt      = fs->make<TH1D>( "jetpt" , " AK5PFJETS_CHS " , 100 , 0. , 1000. );

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
   using namespace reco;

   Handle<PFJetCollection> pfjetH;
   iEvent.getByLabel("ak4PFJetsCHS", pfjetH);
   for ( PFJetCollection::const_iterator jet = pfjetH->begin(); jet != pfjetH->end(); ++jet ) {
     jetpt->Fill(jet->pt());
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
bare_ntuple::beginJob()
{
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
