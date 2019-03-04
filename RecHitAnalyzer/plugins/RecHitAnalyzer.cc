// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
// 
//
// Original Author:  Michael Andrews
//         Created:  Sat, 14 Jan 2017 17:45:54 GMT
//
//
#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

//
// constructors and destructor
//
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
{

  //EBRecHitCollectionT_ = iConfig.getParameter<edm::InputTag>("EBRecHitCollection");
  EBRecHitCollectionT_ = iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection");
  //EBDigiCollectionT_ = iConfig.getParameter<edm::InputTag>("selectedEBDigiCollection");
  //EBDigiCollectionT_ = iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  EERecHitCollectionT_ = iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection");
  //EERecHitCollectionT_ = iConfig.getParameter<edm::InputTag>("EERecHitCollection"));
  HBHERecHitCollectionT_ = iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection");
  //TRKRecHitCollectionT_  = iConfig.getParameter<edm::InputTag>("trackRecHitCollection");

  genParticleCollectionT_ = iConfig.getParameter<edm::InputTag>("genParticleCollection");
  photonCollectionT_ = iConfig.getParameter<edm::InputTag>("gedPhotonCollection");
  //photonCollectionT_ = iConfig.getParameter<edm::InputTag>("photonCollection");
  //jetCollectionT_ = iConfig.getParameter<edm::InputTag>("ak4PFJetCollection");
  jetCollectionT_ = iConfig.getParameter<edm::InputTag>("PFJetCollection");
  genJetCollectionT_ = iConfig.getParameter<edm::InputTag>("genJetCollection");
  trackCollectionT_ = iConfig.getParameter<edm::InputTag>("trackCollection");
  pfCandCollectionT_ = iConfig.getParameter<edm::InputTag>("pfCollection");


  siPixelRecHitCollectionT_   = iConfig.getParameter<edm::InputTag>("siPixelRecHitCollection");
  siStripRecHitCollectionT_ = iConfig.getParameter<std::vector<edm::InputTag> >("siStripRecHitCollection");

  mode_      = iConfig.getParameter<std::string>("mode");
  minJetPt_  = iConfig.getParameter<double>("minJetPt");
  maxJetEta_ = iConfig.getParameter<double>("maxJetEta");
  std::cout << " >> Mode set to " << mode_ << std::endl;
  if ( mode_ == "JetLevel" ) {
    doJets_ = true;
    nJets_ = iConfig.getParameter<int>("nJets");
    std::cout << "\t>> nJets set to " << nJets_ << std::endl;
  } else if ( mode_ == "EventLevel" ) {
    doJets_ = false;
  } else {
    std::cout << " >> Assuming EventLevel Config. " << std::endl;
    doJets_ = false;
  }

  // Initialize file writer
  // NOTE: initializing dynamic-memory histograms outside of TFileService
  // will cause memory leaks
  //usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_sel = fs->make<TH1F>("h_sel", "isSelected;isSelected;Events", 2, 0., 2.);

  //////////// TTree //////////

  // These will be use to create the actual images
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  //branchesEvtSel       ( RHTree, fs );
  //branchesEvtSel_jet   ( RHTree, fs );
  if ( doJets_ ) {
    branchesEvtSel_jet( RHTree, fs );
  } else {
    branchesEvtSel( RHTree, fs );
  }
  branchesEB           ( RHTree, fs );
  branchesEE           ( RHTree, fs );
  branchesHBHE         ( RHTree, fs );
  //branchesECALatHCAL   ( RHTree, fs );
  branchesECALstitched ( RHTree, fs );
  branchesHCALatEBEE   ( RHTree, fs );
  branchesTracksAtEBEE(RHTree, fs);
  branchesTracksAtECALstitched( RHTree, fs);
  branchesPFCandsAtECALstitched( RHTree, fs);
  branchesTRKlayersAtEBEE(RHTree, fs);
  branchesTRKlayersAtECALstitched(RHTree, fs);
  //branchesTRKvolumeAtEBEE(RHTree, fs);
  //branchesTRKvolumeAtECAL(RHTree, fs);

  // For FC inputs
  //RHTree->Branch("FC_inputs",      &vFC_inputs_);

} // constructor
//
RecHitAnalyzer::~RecHitAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}
//
// member functions
//
// ------------ method called for each event  ------------
void
RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  nTotal++;

  // Check run-dependent count
  unsigned int run = iEvent.id().run();
  //if ( run == 194533 && !(runCount[0] < runTotal[0]) ) {
  //  return;
  //} else if ( run == 200519 && !(runCount[1] < runTotal[1]) ) {
  //  return;
  //} else if ( run == 206859 && !(runCount[2] < runTotal[2]) ) {
  //  return;
  //}

  // ----- Apply event selection cuts ----- //

  //edm::Handle<reco::PhotonCollection> photons;
  //iEvent.getByLabel(photonCollectionT_, photons);
  //edm::Handle<reco::GenParticleCollection> genParticles;
  //iEvent.getByLabel(genParticleCollectionT_, genParticles);

  //std::cout << "GenCol.size: " << genParticles->size() << std::endl;
  bool passedSelection = false;
  //passedSelection = runEvtSel( iEvent, iSetup );
  //passedSelection = runEvtSel_jet( iEvent, iSetup );

  if ( doJets_ ) {
    passedSelection = runEvtSel_jet( iEvent, iSetup );
  } else {
    passedSelection = runEvtSel( iEvent, iSetup );
  }

  if ( !passedSelection ) {
    h_sel->Fill( 0. );
    return;
  }

  fillEB( iEvent, iSetup );
  fillEE( iEvent, iSetup );
  fillHBHE( iEvent, iSetup );
  //fillECALatHCAL( iEvent, iSetup );
  fillECALstitched( iEvent, iSetup );
  fillHCALatEBEE( iEvent, iSetup );
  fillTracksAtEBEE( iEvent, iSetup );
  for (int i=0;i<3;i++)
  {
    fillTracksAtECALstitched( iEvent, iSetup, i );
  }
  fillPFCandsAtECALstitched( iEvent, iSetup );
  fillTRKlayersAtEBEE( iEvent, iSetup );
  fillTRKlayersAtECALstitched( iEvent, iSetup );
  //fillTRKvolumeAtEBEE( iEvent, iSetup );
  //fillTRKvolumeAtECAL( iEvent, iSetup );

  ////////////// 4-Momenta //////////
  //fillFC( iEvent, iSetup );

  // Fill RHTree
  RHTree->Fill();
  h_sel->Fill( 1. );
  if ( run == 194533 ) {
    runCount[0]++;
  } else if ( run == 200519 ) {
    runCount[1]++;
  } else if ( run == 206859 ) {
    runCount[2]++;
  } else {
    runCount[3]++;
  }

} // analyze()


// ------------ method called once each job just before starting event loop  ------------
void 
RecHitAnalyzer::beginJob()
{
  for ( int i=0; i < 4; i++ ) {
    runCount[i] = 0;
  }
  nTotal = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecHitAnalyzer::endJob() 
{
  for ( int i=0; i < 4; i++ ) {
    std::cout << " >> i:" << i << " runCount:" << runCount[i] << std::endl;
  }
  std::cout << " >> nPassed[0:2]: " << runCount[0]+runCount[1]+runCount[2] << std::endl;
  std::cout << " >> nTotal: " << nTotal << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

/*
//____ Fill FC diphoton variables _____//
void RecHitAnalyzer::fillFC ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  vFC_inputs_.clear();

  int ptOrder[2] = {0, 1};
  if ( vPho_[1].Pt() > vPho_[0].Pt() ) {
      ptOrder[0] = 1;
      ptOrder[1] = 0;
  }
  for ( int i = 0; i < 2; i++ ) {
    vFC_inputs_.push_back( vPho_[ptOrder[i]].Pt()/m0_ );
    vFC_inputs_.push_back( vPho_[ptOrder[i]].Eta() );
  }
  vFC_inputs_.push_back( TMath::Cos(vPho_[0].Phi()-vPho_[1].Phi()) );

} // fillFC() 
*/

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
