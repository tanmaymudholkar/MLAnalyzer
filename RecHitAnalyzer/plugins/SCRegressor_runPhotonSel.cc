#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

//TODO: Copy pi0 struct implementation in runPiSel()

// Initialize branches _____________________________________________________//
void SCRegressor::branchesPhotonSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  hSC_pT = fs->make<TH1F>("SC_pT", "Pt", 27, 15., 150.);

  tree->Branch("SC_mass",   &vSC_mass_);
  tree->Branch("SC_DR",     &vSC_DR_);
  tree->Branch("SC_pT",     &vSC_pT_);

}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  // Get reco (pT>5GeV) photons associated to gen photons
  float dR;
  float minDR = 100.;
  int minDR_idx = -1;
  mGenPi0_RecoPho.clear();
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );

    if ( iGen->pdgId() != 22 ) continue; //photon
    //if ( iGen->pdgId() != 11 ) continue; //electron
    if ( iGen->status() != 1 ) continue;

    mGenPi0_RecoPho.insert( std::pair<unsigned int, std::vector<unsigned int>>(iG, std::vector<unsigned int>()) );

    // Loop over reco photon collection
    minDR = 100.;
    minDR_idx = -1;
    for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {

      PhotonRef iPho( photons, iP );
      if ( iPho->pt() < 5. ) continue;
      dR = reco::deltaR( iPho->eta(),iPho->phi(), iGen->eta(),iGen->phi() );
      if ( dR > minDR ) continue;

      minDR = dR;
      minDR_idx = iP;
      if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << std::endl;

    } // reco photons
    if ( minDR > 0.04 ) continue;
    mGenPi0_RecoPho[iG].push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << std::endl;

  } // gen photons

  ////////// Apply selection //////////

  float ptCut = 15., etaCut = 1.44;

  if ( debug ) std::cout << " >> PhoCol.size: " << photons->size() << std::endl;
  // Ensure only 1 reco photon associated to each gen pho
  std::vector<unsigned int> vRecoPhoIdxs_;
  for ( auto const& mG : mGenPi0_RecoPho ) {
    if ( mG.second.empty() ) continue;
    if ( mG.second.size() != 1 ) continue;
    vRecoPhoIdxs_.push_back( mG.second[0] );
  }
  if ( debug ) std::cout << " >> RecoPhos.size: " << vRecoPhoIdxs_.size() << std::endl;
  if ( vRecoPhoIdxs_.empty() ) return false;

  // Ensure each pre-selected photon is isolated
  bool isIso;
  vPreselPhoIdxs_.clear();
  for ( unsigned int i = 0; i < vRecoPhoIdxs_.size(); i++ ) {

    PhotonRef iPho( photons, vRecoPhoIdxs_[i] );

    if ( std::abs(iPho->pt()) <= ptCut ) continue;
    if ( std::abs(iPho->eta()) >= etaCut ) continue;
    if ( debug ) std::cout << " >> pT: " << iPho->pt() << " eta: " << iPho->eta() << std::endl;

    /*
    if ( iPho->full5x5_r9() <= 0.5 ) continue;
    if ( iPho->hadTowOverEm() >= 0.08 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;
    //if ( iPho->passElectronVeto() == true ) continue;
    //if ( iPho->userFloat("phoChargedIsolation")/std::abs(iPho->pt()) > 0.3 ) continue;

    if ( iPho->full5x5_r9() <= 0.85 ) {
      if ( iPho->full5x5_sigmaIetaIeta() >= 0.015 ) continue;
      if ( iPho->userFloat("phoPhotonIsolation") >= 4.0 ) continue;
      if ( iPho->trkSumPtHollowConeDR03() >= 6. ) continue;
      //if ( iPho->trackIso() >= 6. ) continue;
    }
    */

    ///*
    // Ensure pre-sel photons are isolated
    isIso = true;
    for ( unsigned int j = 0; j < photons->size(); j++ ) {

      if ( j == vRecoPhoIdxs_[i] ) continue;
      PhotonRef jPho( photons, j );
      //if ( jPho->pt() < 5. ) continue;
      if ( iPho->pt() < 10. ) continue; // pat/miniaod threshold is pt > 10 GeV
      dR = reco::deltaR( iPho->eta(),iPho->phi(), jPho->eta(),jPho->phi() );
      if ( debug ) std::cout << "   >> reco dR:" << dR << std::endl;
      if ( dR > 16*.0174 ) continue;
      isIso = false;
      break;
    } // reco photon j
    if ( !isIso ) continue;
    //*/

    vPreselPhoIdxs_.push_back( vRecoPhoIdxs_[i] );

  } // reco photon i
  if ( debug ) std::cout << " >> PreselPhos.size: " << vPreselPhoIdxs_.size() << std::endl;
  if ( vPreselPhoIdxs_.empty() ) return false;

  if ( debug ) std::cout << " >> Passed selection. " << std::endl;
  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  // Eliminate gen pi0s which are unmatched
  for ( auto const& mG : mGenPi0_RecoPho ) {
    if ( mG.second.empty() ) {
      mGenPi0_RecoPho.erase( mG.first );
      continue;
    }
    if ( std::find(vRegressPhoIdxs_.begin(), vRegressPhoIdxs_.end(), mG.second[0]) == vRegressPhoIdxs_.end() )
      mGenPi0_RecoPho.erase( mG.first );
  }
  if ( debug ) std::cout << " >> mGenPi0.size: " << mGenPi0_RecoPho.size() << std::endl;

  if ( debug ) {
    int iter_idx = 0;
    for ( auto const& mG : mGenPi0_RecoPho ) {
      std::cout << " >> mGenPi0_RecoPho[" << mG.first << "]: " << mG.second[0] << " vs. " << vRegressPhoIdxs_[iter_idx] << std::endl;
      iter_idx++;
    }
  }

  ////////// Store kinematics //////////

  vSC_DR_.clear();
  vSC_pT_.clear();
  vSC_mass_.clear();
  for ( auto const& mG : mGenPi0_RecoPho ) {
    reco::GenParticleRef iGen( genParticles, mG.first );
    vSC_DR_.push_back( 0. );
    vSC_pT_.push_back( iGen->pt() );
    vSC_mass_.push_back( 0. );

    hSC_pT->Fill( iGen->pt() );
  }

}
