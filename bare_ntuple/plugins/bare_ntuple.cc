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
#include <iostream>
#include <memory>
#include <algorithm>
#include <iterator>
#include <vector>
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
#define MAXMUON 30

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
  typedef std::pair<float,double> OrderPair;
  struct IndexByPt {

    const bool operator() (const OrderPair& j1, const OrderPair& j2) const {return j1.second > j2.second;}

  };  

  float deltaPhi(float v1, float v2);
  
  edm::InputTag offlineprimaryvertices_;
  edm::InputTag pfj_;
  edm::InputTag muon_;
  edm::InputTag met_;
  edm::InputTag electron_;
  edm::InputTag tau_;
  edm::InputTag photon_;

  TTree *mtree;

  IndexByPt indexComparator;

  int count_event;
  int cut_1;
  int cut_2;
  int cut_3;
  int cut_4;
  int cut_5;
  int cut_6;
  int cut_7;
  int cut_8;
  int cut_9;


  
  int mrun;
  int mlumi;
  int mevent;
  int mnvtx;
  int mnjets;

  int mnmuonsveto;

  //double  mMuonPt[MAXMUON];
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
  offlineprimaryvertices_(iConfig.getUntrackedParameter<edm::InputTag>("offlineprimaryvertices")),
  pfj_(iConfig.getUntrackedParameter<edm::InputTag>("pfj")),
  //muon_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muon"))),
  muon_(iConfig.getUntrackedParameter<edm::InputTag>("muon")),
  met_(iConfig.getUntrackedParameter<edm::InputTag>("met")),
  electron_(iConfig.getUntrackedParameter<edm::InputTag>("electron")),
  tau_(iConfig.getUntrackedParameter<edm::InputTag>("tau")),
  photon_(iConfig.getUntrackedParameter<edm::InputTag>("photon"))
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;

  count_event=0;

  mtree = fs->make<TTree>( "ntuple" , " monojet_input " );

  cut_1=0;
  cut_2=0;
  cut_3=0;
  cut_4=0;
  cut_5=0;
  cut_6=0;
  cut_7=0;
  cut_8=0;
  cut_9=0;

}

bare_ntuple::~bare_ntuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

float bare_ntuple::deltaPhi(float v1, float v2) 
{
  float diff = fabs(v2 - v1);
  float corr = 2*acos(-1.) - diff;
  if (diff < acos(-1.))
    return diff;
  else
    return corr;
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

   count_event++;

   mrun   = iEvent.id().run();
   mlumi  = iEvent.luminosityBlock();
   mevent = iEvent.id().event();

   ////////////////Vertex////////////////
   reco::Vertex primaryvtx;
   Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel(offlineprimaryvertices_,recVtxs);
   int  pvind=0;
   for(unsigned int ind=0;ind<recVtxs->size();ind++){
     if(!((*recVtxs)[ind].isFake())){
       if(pvind == 0) primaryvtx = (*recVtxs)[ind];
       pvind++;
       if(pvind>=MAXVTX) break;
     }
   }
   mnvtx = pvind;

   bool cut1=false;
   if (pvind > 0){cut_1++; cut1=true;}

   /////////////////RecoJet/////////////////////
   unsigned int count=0;
   Handle<pat::JetCollection> pfjetH;
   iEvent.getByLabel(pfj_, pfjetH);

   std::vector<OrderPair> jet_ordered;

   for ( unsigned int ind=0; ind<(*pfjetH).size() && ind<MAXMUON; ind++ )
     {
       const pat::Jet& jet = (*pfjetH)[ind];

       //Variable for defining a Jet LooseID 
       float NHF = jet.neutralHadronEnergyFraction();

       float NEMF = jet.neutralEmEnergyFraction();

       float CHF = jet.chargedHadronEnergyFraction();

       float MUF = jet.muonEnergyFraction();

       float CEMF = jet.chargedEmEnergyFraction();

       float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();

       float CHM = jet.chargedMultiplicity();

       //float pileupJetId = -999;
       //if ( jet.hasUserFloat("pileupJetId:fullDiscriminant") ) pileupJetId = jet.userFloat("pileupJetId:fullDiscriminant");
       
       //Prescribed LooseID
       if ( (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((abs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jet.eta())>2.4) )
	 {
	   if ( jet.hasUserFloat("pileupJetId:fullDiscriminant") && jet.pt() > 30 && abs(jet.eta()) < 2.5 && (CHF > 0.2 && NHF < 0.7 && NEMF < 0.7) )
	     {jet_ordered.push_back(std::make_pair(jet.phi(),jet.pt())); count++;}
	   else {jet_ordered.push_back(std::make_pair(5.0,-1.0));}
	 }
       else {jet_ordered.push_back(std::make_pair(5.0,-1.0));}
  
     }
   
   std::sort(jet_ordered.begin(), jet_ordered.end(), indexComparator);

   bool cut2 =false;
   bool cut3 =false;

   if ( jet_ordered[0].second > 110 && jet_ordered[0].second > 0 /*&& (CHF > 0.2 && NHF < 0.7 && NEMF < 0.7)*/ ){if (cut1){cut_2++; cut2=true;}}

   if ( (abs(jet_ordered[0].first) < 4 && abs(jet_ordered[1].first) < 4 ) && (deltaPhi(jet_ordered[0].first,jet_ordered[1].first) < 2.5) )
     { if (cut1 && cut2) { cut_3++; cut3=true; } }
   
   if ( jet_ordered[0].second > 110 )
     {
       if ( (abs(jet_ordered[0].first) < 4 && abs(jet_ordered[1].first) < 4 ) && (deltaPhi(jet_ordered[0].first,jet_ordered[1].first) < 2.5) )
	 {
	   if ( count > 0 && ( count == 1 || count ==2 ) ) mnjets = count;
	 }
     }
   /////////////////MET////////////////////////////
   Handle<pat::METCollection> metcorr;
   iEvent.getByLabel(met_,metcorr);

   bool cut4 =false;

   for ( pat::METCollection::const_iterator mety = metcorr->begin(); mety != metcorr->end() ; ++mety )
     {
       if( mety->pt() > 200){ if(cut1 && cut2 && cut3) {cut_4++; cut4=true;} };
     }

   bool cut5 =false;
   if( count < 3){if ( cut1 && cut2 && cut3 && cut4) {cut_5++; cut5=true;}}


   /////////////////PFMuons///////////////////////
   unsigned int muon_count=0;
   Handle<pat::MuonCollection> pfmuon;
   iEvent.getByLabel(muon_, pfmuon);
   //iEvent.getByToken(muon_, pfmuon);

   for ( unsigned int ind=0; ind<(*pfmuon).size() && ind<MAXMUON; ind++ )
     {
       const pat::Muon& muon = (*pfmuon)[ind];
       if ( !(muon.pt() > 10 && abs(muon.eta()) < 2.4 && muon.isLooseMuon()) )
	 muon_count++;
     }
   mnmuonsveto = muon_count;

   bool cut6=false;
   if (muon_count == 0){if ( cut1 && cut2 && cut3 && cut4 && cut5 ) {cut_6++; cut6=true;}}

   ///////////////////PFElectron////////////////////
   unsigned int e_count=0;
   Handle<pat::ElectronCollection> elec;
   iEvent.getByLabel(electron_, elec);                                                                                                               

   for ( unsigned int ind=0; ind<(*elec).size() && ind<MAXMUON; ind++ )
     {
       const pat::Electron& electron = (*elec)[ind];
       if ( !(electron.pt() > 10 && abs(electron.eta()) < 2.5) )
         e_count++;
     }

   bool cut7=false;
   if (e_count == 0){if ( cut1 && cut2 && cut3 && cut4 && cut5 && cut6 ) {cut_7++; cut7=true;}}


   ///////////////////PFTau///////////////
   unsigned int tau_count=0;
   Handle<pat::TauCollection> pftau;
   iEvent.getByLabel(tau_, pftau);
   
   for ( unsigned int ind=0; ind<(*pftau).size() && ind<MAXMUON; ind++ )
     {
       const pat::Tau& tau = (*pftau)[ind];
       if ( !(tau.pt() > 18 && abs(tau.eta()) < 2.3) )
         tau_count++;
     }

   bool cut8=false;
   if (tau_count == 0){if ( cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 ) {cut_8++; cut8=true;}}

   //////////////////////photon//////////////////
   unsigned int photon_count=0;
   Handle<pat::PhotonCollection> pho;
   iEvent.getByLabel(photon_, pho);

   for ( unsigned int ind=0; ind<(*pho).size() && ind<MAXMUON; ind++ )
     {
       const pat::Photon& photon = (*pho)[ind];
       if ( !(photon.pt() > 15 && abs(photon.eta()) < 2.5) )
         photon_count++;
     }

   if (photon_count == 0){if ( cut1 && cut2 && cut3 && cut4 && cut5 && cut6 && cut7 && cut8 ) cut_9++;}



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
  mtree->Branch( "nmuonsveto" , &mnmuonsveto , "nmuonsveto/I" );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
bare_ntuple::endJob() 
{

  std::cout<<"Number of event = "<<count_event<<std::endl;
  std::cout<<"Number of good vertex = "<<cut_1<<std::endl;
  std::cout<<"Number of pt > 110 GeV = "<<cut_2<<std::endl;
  std::cout<<"Number of deltaphi = "<<cut_3<<std::endl;
  std::cout<<"Number of MET > 200 GeV ="<<cut_4<<std::endl;
  std::cout<<"Number of N < 3 ="<<cut_5<<std::endl;
  std::cout<<"Number of Lepton veto ="<<cut_6+cut_7<<std::endl;
  std::cout<<"Number of Tau veto ="<<cut_8<<std::endl;
  std::cout<<"Number of Photon veto ="<<cut_9<<std::endl;

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
