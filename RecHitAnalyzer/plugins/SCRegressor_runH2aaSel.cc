#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

struct gen_obj {
  unsigned int idx;
  double pt;
};

std::vector<gen_obj> vAs;
std::vector<unsigned int> vGenAIdxs;

// Initialize branches _____________________________________________________//
void SCRegressor::branchesH2aaSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("mHgen",     &mHgen_);
  //tree->Branch("FC_inputs", &vFC_inputs_);
  //tree->Branch("hltAccept", &hltAccept_);

  tree->Branch("A_mass",    &vA_mass_);
  tree->Branch("A_DR",      &vA_DR_);
  tree->Branch("A_E",       &vA_E_);
  tree->Branch("A_pT",      &vA_pT_);
  tree->Branch("A_eta",     &vA_eta_);
  tree->Branch("A_phi",     &vA_phi_);
  tree->Branch("A_recoIdx", &vA_recoIdx_);

  hdPhidEta = fs->make<TH2F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma)",
              6, 0., 6.*0.0174, 6, 0., 6.*0.0174);
}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runH2aaSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  vAs.clear();
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vH;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );
    if ( std::abs(iGen->pdgId()) != 35 ) continue;
    if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 25 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    if ( iGen->daughter(0)->pdgId() != 22 || iGen->daughter(1)->pdgId() != 22 ) continue;

    gen_obj Gen_obj = { iG, std::abs(iGen->pt()) };
    vAs.push_back( Gen_obj );
    vH += iGen->p4();

  } // gen particles
  if ( vAs.size() != 2 ) return false;
  mHgen_ = vH.mass();

  // Sort As by pT, for abitrary N
  std::sort( vAs.begin(), vAs.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  vGenAIdxs.clear();
  //vPreselPhoIdxs_.clear();
  //bool keepEvt = true;
  //float ptomcut[2] = { 125./3., 125./4. };
  for ( unsigned int iG = 0; iG < vAs.size(); iG++ ) {
    reco::GenParticleRef iGen( genParticles, vAs[iG].idx );
    if ( debug ) std::cout << " >> pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy() << std::endl;
    //if ( std::abs(iGen->eta()) > 1.21 ) keepEvt = false;
    //if ( std::abs(iGen->pt()) < ptomcut[iG] ) keepEvt = false;
    //if ( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) > 0.15 ) keepEvt = false;
    //vPreselPhoIdxs_.push_back( vAs[iG].idx );
    vGenAIdxs.push_back( vAs[iG].idx );
  }
  //if (keepEvt == false) return false;

  /*
  edm::Handle<edm::TriggerResults> trgs;
  iEvent.getByToken( trgResultsT_, trgs );

  const edm::TriggerNames &triggerNames = iEvent.triggerNames( *trgs );
  //std::cout << ">> N triggers:" << trgs->size() << std::endl;
  for ( unsigned int iT = 0; iT < trgs->size(); iT++ ) {
    //std::cout << " name["<<iT<<"]:"<<triggerNames.triggerName(iT)<< std::endl;
  }

  int hltAccept = -1;
  std::string trgName = "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_*_Mass55_v*";
  std::vector< std::vector<std::string>::const_iterator > trgMatches = edm::regexMatch( triggerNames.triggerNames(), trgName );

  if ( !trgMatches.empty() ) {
    hltAccept = 0;
    for ( auto const& iT : trgMatches ) {
      if ( trgs->accept(triggerNames.triggerIndex(*iT)) ) hltAccept = 1;
    }
  }
  hltAccept_ = hltAccept;
  */

  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillH2aaSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  vA_E_.clear();
  vA_pT_.clear();
  vA_eta_.clear();
  vA_phi_.clear();
  vA_mass_.clear();
  vA_DR_.clear();
  vA_recoIdx_.clear();
  float dPhi, dEta, dR, recoDR;
  int recoDR_idx;
  for ( unsigned int iG : vGenAIdxs ) {

    reco::GenParticleRef iGen( genParticles, iG );

    vA_E_.push_back( std::abs(iGen->energy()) );
    vA_pT_.push_back( std::abs(iGen->pt()) );
    vA_eta_.push_back( iGen->eta() );
    vA_phi_.push_back( iGen->phi() );
    vA_mass_.push_back( iGen->mass() );
    vA_DR_.push_back( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) );

    dPhi = reco::deltaPhi( iGen->daughter(0)->phi(), iGen->daughter(1)->phi() );
    dEta = std::abs( iGen->daughter(0)->eta() - iGen->daughter(1)->eta() );
    hdPhidEta->Fill( dPhi, dEta );

    // Get index to dR-matched preselected photon
    recoDR = 2*0.04;
    recoDR_idx = -1;
    for ( unsigned int iP : vRegressPhoIdxs_ ) {
      PhotonRef iPho( photons, iP );
      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());
      if ( dR < recoDR ) {
        recoDR = dR;
        recoDR_idx = iP;
      }
    } // reco pho
    vA_recoIdx_.push_back( recoDR_idx );

  } // gen A

}
