#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesZJetsEleSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("m0",        &m0_);
  tree->Branch("FC_inputs", &vFC_inputs_);
  tree->Branch("hltAccept", &hltAccept_);
  tree->Branch("nRecoPho",  &nRecoPho_);
  tree->Branch("minDR",     &vMinDR_);

  hNpassed_kin      = fs->make<TH1F>("hNpassed_kin", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_presel   = fs->make<TH1F>("hNpassed_presel", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_mGG      = fs->make<TH1F>("hNpassed_mGG", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_nRecoPho = fs->make<TH1F>("hNpassed_nRecoPho", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_hlt      = fs->make<TH1F>("hNpassed_hlt", "isPassed;isPassed;N", 2, 0., 2);
}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runZJetsEleSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<ElectronCollection> electrons;
  iEvent.getByToken(electronCollectionT_, electrons);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  ////////// Apply selection //////////

  //std::cout << " Ele collection size:" << electrons->size() << std::endl;

  std::vector<unsigned int> vRecoEleIdxs;
  for ( unsigned int iE = 0; iE < electrons->size(); iE++ ) {
    ElectronRef iEle( electrons, iE );
    if ( std::abs(iEle->pt()) < 10. ) continue;
    if ( std::abs(iEle->eta()) > 2.4 ) continue;
    if ( std::abs(iEle->eta()) > 1.442 && std::abs(iEle->eta()) < 1.566 ) continue;
    // isLooseEle
    if ( iEle->electronID("cutBasedElectronID-Fall17-94X-V2-loose") == false ) continue;
    vRecoEleIdxs.push_back( iE );
    //vDiEle += iEle->p4();
  }
  if ( vRecoEleIdxs.size() < 2 ) return false;
  //std::cout << " vRecoEleIdxs size:" << vRecoEleIdxs.size() << std::endl;

  int nTightEle = 0;
  bool isInZWindow = false;
  for ( unsigned int iE : vRecoEleIdxs ) {

    // Check ele passes tight ID
    ElectronRef iEle( electrons, iE );
    if ( std::abs(iEle->pt()) < 20. ) continue;
    // isTightEle
    if ( iEle->electronID("cutBasedElectronID-Fall17-94X-V2-medium") == false ) continue;
    nTightEle++;

    // Check is Z-mass compatible with other electrons
    isInZWindow = false;
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vDiEle;
    for ( unsigned int jE : vRecoEleIdxs ) {
      if ( iE == jE ) continue;
      ElectronRef jEle( electrons, jE );
      vDiEle = iEle->p4() + jEle->p4();
      if ( vDiEle.mass() < 60. || vDiEle.mass() > 120 ) continue;
      isInZWindow = true;
      m0_ = vDiEle.mass();
      if ( isInZWindow ) break;
    }// j electrons

    // Stop at first compatible match
    if ( isInZWindow ) break;

  } // i electrons
  if ( nTightEle < 1 ) return false;
  if ( !isInZWindow ) return false;
  //std::cout << " vDiEle.mass():" << m0_<< std::endl;

  //std::cout << " Pho collection size:" << photons->size() << std::endl;
  // Ensure jet passes photon presel
  std::vector<unsigned int> vRecoPhoIdxs;
  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {

    PhotonRef iPho( photons, iP );

    if ( std::abs(iPho->pt()) <= 25. ) continue;
    if ( std::abs(iPho->eta()) >= 1.442 ) continue;

    ///*
    if ( iPho->full5x5_r9() <= 0.5 ) continue;
    if ( iPho->hadTowOverEm() >= 0.08 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;
    //*/
    ///*
    if ( iPho->full5x5_r9() <= 0.85 ) {
      if ( iPho->full5x5_sigmaIetaIeta() >= 0.015 ) continue;
      if ( iPho->userFloat("phoPhotonIsolation") >= 4.0 ) continue;
      //if ( iPho->photonIso() >= 4.0 ) continue;
      if ( iPho->trkSumPtHollowConeDR03() >= 6. ) continue;
      //if ( iPho->trackIso() >= 6. ) continue;
    }
    //*/
    vRecoPhoIdxs.push_back( iP );

  } // photons
  if ( vRecoPhoIdxs.size() < 1 ) return false;
  //std::cout << " vRecoPhoIdxs.size():" << vRecoPhoIdxs.size() << std::endl;

  // Ensure presel photon doesnt overlap with eles
  float dR;
  bool isIsolated = true;
  for ( unsigned int iP : vRecoPhoIdxs ) {
    PhotonRef iPho( photons, iP );
    isIsolated = true;
    for ( unsigned int iE : vRecoEleIdxs ) {
      ElectronRef iEle( electrons, iE );
      dR = reco::deltaR( iPho->eta(),iPho->phi(), iEle->eta(),iEle->phi() );
      //std::cout << " dR to ele["<<iE<<"]:" << dR << std::endl;
      if ( dR < 0.3 ) {
        isIsolated = false;
        break;
      } // dR
    } // vRecoEleIdxs
    if ( !isIsolated ) continue;
    vPreselPhoIdxs_.push_back( iP );
  } // vRecoPhoIdxs
  //std::cout << " n passing ZJetsEle sel:" << vPreselPhoIdxs_.size() << std::endl;

  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillZJetsEleSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  /*
  // Fill kinematic variables
  float dphi[2] = { 0., 0. };
  vFC_inputs_.clear();
  for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
    PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
    vFC_inputs_.push_back( iPho->pt()/m0_ );
    vFC_inputs_.push_back( iPho->eta() );
    dphi[iP] = iPho->phi();
  }
  vFC_inputs_.push_back( TMath::Cos(reco::deltaPhi(dphi[0], dphi[1])) );

  // Get dR of closest reco photon to presel photon
  float minDR, dR;
  vMinDR_.clear();
  for ( unsigned int jP = 0; jP < vRegressPhoIdxs_.size(); jP++ ) {
    PhotonRef jPho( photons, vRegressPhoIdxs_[jP] );
    minDR = 100.;
    for ( unsigned int kP = 0; kP < photons->size(); kP++ ) {
      if ( std::find(vRegressPhoIdxs_.begin(), vRegressPhoIdxs_.end(), kP) != vRegressPhoIdxs_.end() ) continue;
      PhotonRef kPho( photons, kP );
      if ( std::abs(kPho->pt()) < 10. ) continue;
      dR = reco::deltaR( jPho->eta(),jPho->phi(), kPho->eta(),kPho->phi() );
      if ( dR < minDR ) minDR = dR;
    } //k
    vMinDR_.push_back( minDR );
  } //j
  */

}
