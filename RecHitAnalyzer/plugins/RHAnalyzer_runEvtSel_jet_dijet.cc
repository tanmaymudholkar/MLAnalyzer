#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

using std::vector;

TH1D *h_dijet_jet_pT;
TH1D *h_dijet_jet_E;
TH1D *h_dijet_jet_eta;
TH1D *h_dijet_jet_m0;
TH1D *h_dijet_jet_nJet;
vector<float> vDijet_jet_pT_;
vector<float> vDijet_jet_m0_;
vector<float> vDijet_jet_eta_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet( TTree* tree, edm::Service<TFileService> &fs ) {

  h_dijet_jet_pT    = fs->make<TH1D>("h_jet_pT"  , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_dijet_jet_E     = fs->make<TH1D>("h_jet_E"   , "E;E;Particles"        , 100,  0., 800.);
  h_dijet_jet_eta   = fs->make<TH1D>("h_jet_eta" , "#eta;#eta;Particles"  , 100, -5., 5.);
  h_dijet_jet_nJet  = fs->make<TH1D>("h_jet_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_dijet_jet_m0    = fs->make<TH1D>("h_jet_m0"  , "m0;m0;Events"         , 100,  0., 100.);

  tree->Branch("jetPt",  &vDijet_jet_pT_);
  tree->Branch("jetM",   &vDijet_jet_m0_);
  tree->Branch("jetEta", &vDijet_jet_eta_);

} // branchesEvtSel_jet_dijet()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(jetCollectionT_, jets);
  iEvent.getByLabel(genParticleCollectionT_, genParticles);

  vJetIdxs.clear();
  vDijet_jet_pT_.clear();
  vDijet_jet_m0_.clear();
  vDijet_jet_eta_.clear();

  std::vector<TLorentzVector> had_tops,bdau,wdau;
  for (const auto & p : *genParticles.product())
  {
    int id = p.pdgId();
    if(abs(id) != 6 || p.numberOfDaughters()!=2) continue;
    int iw=-1;
    int ib=-1;
    if (abs(p.daughter(0)->pdgId())==24 && abs(p.daughter(1)->pdgId())==5)
    {
      iw=0;ib=1;
    }
    else
    {
      if(abs(p.daughter(1)->pdgId())==24 && abs(p.daughter(0)->pdgId())==5)
      {
        iw=1;ib=0;
      }
      else continue;
    }
    const reco::Candidate *d = p.daughter(iw);
    const reco::Candidate *b = p.daughter(ib);
    while(d->numberOfDaughters() == 1) d = d->daughter(0);
    if(!(abs(d->daughter(0)->pdgId()) < 10 && abs(d->daughter(1)->pdgId()) < 10)) continue;
    TLorentzVector the_top,the_w,the_b;
    the_top.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
    the_w.SetPtEtaPhiE(d->pt(),d->eta(),d->phi(),d->energy());
    the_b.SetPtEtaPhiE(b->pt(),b->eta(),b->phi(),b->energy());
    had_tops.push_back(the_top);
    wdau.push_back(the_w);
    bdau.push_back(the_b);
  }

  int nJet = 0;
  // Loop over jets
  for ( unsigned ihad=0;ihad<had_tops.size();ihad++)
  {
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ )
    {
      reco::PFJetRef iJet( jets, iJ );
      TLorentzVector vjet;
      vjet.SetPtEtaPhiE(iJet->pt(),iJet->eta(),iJet->phi(),iJet->energy());

      if ( std::abs(iJet->pt()) < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_) continue;
      if (had_tops[ihad].DeltaR(vjet)>0.8) continue;
      if (wdau[ihad].DeltaR(vjet)>0.8) continue;
      if (bdau[ihad].DeltaR(vjet)>0.8) continue;

      if ( debug ) std::cout << " >> jet[" << iJ << "]Pt:" << iJet->pt() << " jetE:" << iJet->energy() << " jetM:" << iJet->mass() << std::endl;

      vJetIdxs.push_back(iJ);

      nJet++;
      //break; This should allow two hardonic tops
    } // jets
    if ( (nJets_ > 0) && (nJet >= nJets_) ) break;
  } // hadronic tops


  if ( debug ) {
    for(int thisJetIdx : vJetIdxs)
      std::cout << " >> vJetIdxs:" << thisJetIdx << std::endl;
  }

  if ( (nJets_ > 0) && (nJet != nJets_) ){
    if ( debug ) std::cout << " Fail jet multiplicity:  " << nJet << " < " << nJets_ << std::endl;
    return false;
  }

  if ( vJetIdxs.size() == 0){
    if ( debug ) std::cout << " No passing jets...  " << std::endl;
    return false;
  }

  if ( debug ) std::cout << " >> has_jet_dijet: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet() 

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);

  h_dijet_jet_nJet->Fill( vJetIdxs.size() );
  // Fill branches and histograms 
  for(int thisJetIdx : vJetIdxs){
    reco::PFJetRef thisJet( jets, thisJetIdx );
    if ( debug ) std::cout << " >> Jet[" << thisJetIdx << "] Pt:" << thisJet->pt() << std::endl;
    h_dijet_jet_pT->Fill( std::abs(thisJet->pt()) );
    h_dijet_jet_E->Fill( thisJet->energy() );
    h_dijet_jet_m0->Fill( thisJet->mass() );
    h_dijet_jet_eta->Fill( thisJet->eta() );
    vDijet_jet_pT_.push_back( std::abs(thisJet->pt()) );
    vDijet_jet_m0_.push_back( thisJet->mass() );
    vDijet_jet_eta_.push_back( thisJet->eta() );
  }

} // fillEvtSel_jet_dijet()
