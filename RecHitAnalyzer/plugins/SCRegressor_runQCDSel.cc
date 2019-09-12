#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesQCDSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  //tree->Branch("mHgen",     &mHgen_);
  //tree->Branch("FC_inputs", &vFC_inputs_);
  //tree->Branch("hltAccept", &hltAccept_);

  tree->Branch("A_mass",    &vA_mass_);
  tree->Branch("A_DR",      &vA_DR_);
  tree->Branch("A_E",       &vA_E_);
  tree->Branch("A_pT",      &vA_pT_);
  tree->Branch("A_eta",     &vA_eta_);
  tree->Branch("A_phi",     &vA_phi_);

  tree->Branch("A_pdgId",   &vA_pdgId_);
  tree->Branch("A_mothPdgId",   &vA_mothPdgId_);
  tree->Branch("A_jetM",    &vA_jetM_);

  tree->Branch("OutPart_pdgId",   &vOutPart_pdgId_);

  hJetNeuM  = fs->make<TH1F>("hJetNeuM", "jet m_{neu};jet m_{neu};N_{jet}", 48, 0., 1.2);
  hdPhidEta = fs->make<TH2F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma)",
              6, 0., 6.*0.0174, 6, 0., 6.*0.0174);
}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runQCDSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillQCDSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  vA_E_.assign(2, -99.);
  vA_pT_.assign(2, -99.);
  vA_eta_.assign(2, -99.);
  vA_phi_.assign(2, -99.);
  vA_mass_.assign(2, -99.);
  vA_DR_.assign(2, -99.);
  vA_pdgId_.assign(2, -99);
  vA_mothPdgId_.assign(2, -99);
  vA_jetM_.assign(2, -99);

  vOutPart_pdgId_.assign(2, -99);

  float minDR, dR;
  int minDR_idx;
  for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {

    PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
    /*

    //start comment
    // Get "a": mother of gen-matched photon
    minDR = 2*0.04;
    minDR_idx = -1;
    for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

      reco::GenParticleRef iGen( genParticles, iG );

      if ( abs(iGen->eta()) > 1.442 ) continue;
      //if ( abs(iGen->pt())  < 20. ) continue;
      if ( abs(iGen->pdgId()) != 22 ) continue;
      if ( abs(iGen->status()) != 1 ) continue;

      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());
      if ( dR > minDR ) continue;
      minDR = dR;
      minDR_idx = iG;
    } // gen particles

    if ( minDR_idx == -1 ) continue;

    // Fill a kinematics
    reco::GenParticleRef iGen( genParticles, minDR_idx );
    const reco::Candidate* Acand = iGen->mother();

    vA_E_[iP] = std::abs(Acand->energy());
    vA_pT_[iP] = std::abs(Acand->pt());
    vA_eta_[iP] = Acand->eta();
    vA_phi_[iP] = Acand->phi();
    vA_mass_[iP] = Acand->mass();
    vA_DR_[iP] = minDR;
    vA_pdgId_[iP] = Acand->pdgId();
    if ( !(Acand->mother()) ) continue;
    vA_mothPdgId_[iP] = Acand->mother()->pdgId();

    // Get neutral jet mass: all aunts of gen-matched photon that decay to photons within dR < 0.3
    const reco::Candidate* AcandMoth = Acand->mother();
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > jet_p4;
    for ( unsigned int i = 0; i < AcandMoth->numberOfDaughters(); i++ ) {
      const reco::Candidate* AcandSis = AcandMoth->daughter(i);
      if ( AcandSis->numberOfDaughters() == 0 ) continue;
      if ( AcandSis->daughter(0)->pdgId() != 22 ) continue;
      if ( reco::deltaR(AcandSis->eta(),AcandSis->phi(), iPho->eta(),iPho->phi()) >= 0.3 ) continue;
      jet_p4 += AcandSis->p4();
    }// Acand sisters
    hJetNeuM->Fill( jet_p4.mass() );
    vA_jetM_[iP] = jet_p4.mass();
    //end comment
    */

    // Get pdgId of outgoing parton from hardscatter gen-matched to this presel photon
    float minDR_parton = 0.3;
    float minDR_parton_idx = -1;
    for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

      reco::GenParticleRef iGen( genParticles, iG );

      if ( abs(iGen->status()) != 23 ) continue;

      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());
      if ( dR > minDR_parton ) continue;
      minDR_parton = dR;
      minDR_parton_idx = iG;

    }// gen particles
    if ( minDR_parton_idx == -1 ) continue;
    reco::GenParticleRef iGenPart( genParticles, minDR_parton_idx );
    vOutPart_pdgId_[iP] = iGenPart->pdgId();

  } // presel photons

}
