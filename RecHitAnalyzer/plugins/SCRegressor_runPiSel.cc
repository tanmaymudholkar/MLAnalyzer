#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesPiSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  //hSC_pT = fs->make<TH1F>("SC_pT", "Pt", 65, 30., 160.);
  //hSC_mass = fs->make<TH1F>("SC_mass", "m_{SC};m_{SC}",50, 0., 0.5);
  //hdR = fs->make<TH1F>("dR_seed_subJet", "#DeltaR(seed,subJet);#DeltaR",50, 0., 50.*0.0174);
  //hdEta = fs->make<TH1F>("dEta_seed_subJet", "#Delta#eta(seed,subJet);#Delta#eta",50, 0., 50.*0.0174);
  //hdPhi = fs->make<TH1F>("dPhi_seed_subJet", "#Delta#phi(seed,subJet);#Delta#phi",50, 0., 50.*0.0174);
  hdPhidEta = fs->make<TH3F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma);m",
      6, 0., 6.*0.0174, 6, 0., 6.*0.0174, 16., 0.,1.6);
  hnPho = fs->make<TH2F>("nPho", "N(m_{#pi},p_{T,#pi})_{reco};m_{#pi^{0}};p_{T,#pi^0}",
      16, 0., 1.6, 17, 15., 100.);

  tree->Branch("SC_mass",   &vSC_mass_);
  tree->Branch("SC_DR",     &vSC_DR_);
  tree->Branch("SC_pT",     &vSC_pT_);

}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runPiSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  // ____________ Gen-level studies ________________ //
  float dR;
  mGenPi0_RecoPho.clear();
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );
    // ID cuts
    if ( debug ) std::cout << " >> pdgId:"<< iGen->pdgId() << " status:" << iGen->status() << " nDaughters:" << iGen->numberOfDaughters() << std::endl;
    if ( std::abs(iGen->pdgId()) != 111 ) continue;
    if ( debug ) std::cout << " >> pdgId:111 nDaughters:" << iGen->numberOfDaughters() << std::endl;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    ///*
    dR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
    if ( dR > 10*.0174 ) continue;
    //*/

    mGenPi0_RecoPho.insert( std::pair<unsigned int, std::vector<unsigned int>>(iG, std::vector<unsigned int>()) );

  } // genParticle loop: count good photons
  if ( debug ) std::cout << " >> mGenPi0_RecoPho.size: " << mGenPi0_RecoPho.size() << std::endl;
  if ( mGenPi0_RecoPho.empty() ) return false;

  ////////// Apply selection //////////

  //float ptCut = 45., etaCut = 1.44;
  float ptCut = 15., etaCut = 1.44;

  // Get reco (pT>5GeV) photons associated to pi0
  float minDR = 100.;
  int minDR_idx = -1;
  // Loop over pi0s
  for ( auto& mG : mGenPi0_RecoPho ) {

    reco::GenParticleRef iGen( genParticles, mG.first );

    if ( debug ) std::cout << " >> pi0[" << mG.first << "]"<< std::endl;
    // Loop over photons from pi0
    for ( unsigned int iD = 0; iD < iGen->numberOfDaughters(); iD++ ) {

      const reco::Candidate* iGenPho = iGen->daughter(iD);

      // Loop over reco photon collection
      minDR = 100.;
      minDR_idx = -1;
      if ( debug ) std::cout << "  >> iD[" << iD << "]" << std::endl;
      for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {

        PhotonRef iPho( photons, iP );
        if ( iPho->pt() < 5. ) continue;
        dR = reco::deltaR( iPho->eta(),iPho->phi(), iGenPho->eta(),iGenPho->phi() );
        if ( dR > minDR ) continue;

        minDR = dR;
        minDR_idx = iP;
        if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << std::endl;

      } // reco photons
      if ( minDR > 0.04 ) continue;
      mG.second.push_back( minDR_idx );
      if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << std::endl;

    } // gen photons
    dR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
    if ( debug ) std::cout << " >> gen dR:" << dR << std::endl;

  } // gen pi0s

  // Ensure only 1 reco photon associated to each gen pi0
  std::vector<unsigned int> vRecoPhoIdxs_;
  for ( auto const& mG : mGenPi0_RecoPho ) {
    if ( debug ) std::cout << " >> pi0[" << mG.first << "] size:" << mG.second.size() << std::endl; 
    // Possibilities:
    // 1) pho_gen1: unmatched, pho_gen2: unmatched => both reco pho failed pT cut or dR matching, reject
    // 2) pho_gen1: reco-matched1, pho_gen1: unmatched => one reco pho failed pT cut or dR matching, can accept
    // 3) pho_gen1: reco-matched1, pho_gen2: reco-matched2
    //    3a) reco-matched1 == reco-matched2 => merged, accept
    //    3b) reco-matched1 != reco-matched2 => resolved, reject
    if ( mG.second.empty() ) continue;
    if ( mG.second.size() == 2 && mG.second[0] != mG.second[1] ) continue; // 2 resolved reco photons 
    vRecoPhoIdxs_.push_back( mG.second[0] );
    //if ( mG.second.size() == 2 && mG.second[0] != mG.second[1] ) vRecoPhoIdxs_.push_back( mG.second[1] );
  } 
  if ( debug ) std::cout << " >> RecoPhos.size: " << vRecoPhoIdxs_.size() << std::endl;
  if ( vRecoPhoIdxs_.empty() ) return false;

  // Ensure each of pi0-matched reco photons passes pre-selection
  // NOTE: at this point there is 1-1 correspondence between pi0 and reco pho,
  // so suffices to check there are as many pre-selected phos as pi0s
  bool isIso;
  vPreselPhoIdxs_.clear();
  if ( debug ) std::cout << " >> PhoCol.size: " << photons->size() << std::endl;
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
      if ( jPho->pt() < 5. ) continue;
      dR = reco::deltaR( iPho->eta(),iPho->phi(), jPho->eta(),jPho->phi() );
      if ( debug ) std::cout << "   >> reco dR:" << dR << std::endl;
      if ( dR > 16*.0174 ) continue;
      isIso = false;
      break;

    } // photon j
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
void SCRegressor::fillPiSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
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
  float dEta, dPhi, dR, mPi0, ptPi0;
  for ( auto const& mG : mGenPi0_RecoPho ) {
    reco::GenParticleRef iGen( genParticles, mG.first );
    mPi0 = iGen->mass();
    ptPi0 = iGen->pt();
    dR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
    dEta = std::abs( iGen->daughter(0)->eta() - iGen->daughter(1)->eta() );
    dPhi = reco::deltaPhi( iGen->daughter(0)->phi(), iGen->daughter(1)->phi() );
    if ( debug ) std::cout << " >> m0:" << mPi0 << " dR:" << dR << " dPhi:" << dPhi << std::endl;

    vSC_DR_.push_back( dR );
    vSC_pT_.push_back( iGen->pt() );
    vSC_mass_.push_back( mPi0 );

    //hPt->Fill( ptPi0 );
    hdPhidEta->Fill( dPhi, dEta, mPi0 );
    //hnPho->Fill( mPi0, iGen->pt() );
    hnPho->Fill( mPi0, ptPi0 );
    //hSC_mass->Fill( mPi0 );
  }

}
