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



std::vector<std::vector<int> > seljet_genpart_collid;
std::vector<std::vector<int> > seljet_genpart_pdgid;
std::vector<std::vector<int> > seljet_genpart_charge;

std::vector<std::vector<float> > seljet_genpart_px;
std::vector<std::vector<float> > seljet_genpart_py;
std::vector<std::vector<float> > seljet_genpart_pz;
std::vector<std::vector<float> > seljet_genpart_energy;

std::vector<std::vector<int> > seljet_genpart_status;

std::vector<std::vector<int> > seljet_genpart_motherpdgid;
std::vector<std::vector<int> > seljet_genpart_dau1pdgid;
std::vector<std::vector<int> > seljet_genpart_dau2pdgid;

std::vector<float> seljet_px;
std::vector<float> seljet_py;
std::vector<float> seljet_pz;
std::vector<float> seljet_energy;

std::vector<std::vector<float> > seljet_pfcand_px;
std::vector<std::vector<float> > seljet_pfcand_py;
std::vector<std::vector<float> > seljet_pfcand_pz;
std::vector<std::vector<float> > seljet_pfcand_energy;
std::vector<std::vector<int> > seljet_pfcand_type;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet( TTree* tree, edm::Service<TFileService> &fs ) {

  h_dijet_jet_pT    = fs->make<TH1D>("h_jet_pT"  , "p_{T};p_{T};Particles", 300,  0., 1500.);
  h_dijet_jet_E     = fs->make<TH1D>("h_jet_E"   , "E;E;Particles"        , 100,  0., 800.);
  h_dijet_jet_eta   = fs->make<TH1D>("h_jet_eta" , "#eta;#eta;Particles"  , 100, -5., 5.);
  h_dijet_jet_nJet  = fs->make<TH1D>("h_jet_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_dijet_jet_m0    = fs->make<TH1D>("h_jet_m0"  , "m0;m0;Events"         , 100,  0., 100.);

  tree->Branch("jetPt",  &vDijet_jet_pT_);
  tree->Branch("jetM",   &vDijet_jet_m0_);
  tree->Branch("jetEta", &vDijet_jet_eta_);



  tree->Branch("seljet_genpart_collid", &seljet_genpart_collid);
  tree->Branch("seljet_genpart_pdgid", &seljet_genpart_pdgid);
  tree->Branch("seljet_genpart_charge", &seljet_genpart_charge);

  tree->Branch("seljet_genpart_px", &seljet_genpart_px);
  tree->Branch("seljet_genpart_py", &seljet_genpart_py);
  tree->Branch("seljet_genpart_pz", &seljet_genpart_pz);
  tree->Branch("seljet_genpart_energy", &seljet_genpart_energy);

  tree->Branch("seljet_genpart_status", &seljet_genpart_status);

  tree->Branch("seljet_genpart_motherpdgid", &seljet_genpart_motherpdgid);
  tree->Branch("seljet_genpart_dau1pdgid", &seljet_genpart_dau1pdgid);
  tree->Branch("seljet_genpart_dau2pdgid", &seljet_genpart_dau2pdgid);


  tree->Branch("seljet_px", &seljet_px);
  tree->Branch("seljet_py", &seljet_py);
  tree->Branch("seljet_pz", &seljet_pz);
  tree->Branch("seljet_energy", &seljet_energy);

  tree->Branch("seljet_pfcand_px", &seljet_pfcand_px);
  tree->Branch("seljet_pfcand_py", &seljet_pfcand_py);
  tree->Branch("seljet_pfcand_pz", &seljet_pfcand_pz);
  tree->Branch("seljet_pfcand_energy", &seljet_pfcand_energy);
  tree->Branch("seljet_pfcand_type", &seljet_pfcand_type);

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

  int nJet = 0;

  std::vector<TLorentzVector> had_tops,bdau,wdau;
  if (isTTbar_) { //is a ttbar sample
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
      break;// This should allow two hardonic tops
    } // jets
    if ( (nJets_ > 0) && (nJet >= nJets_) ) break;
  } // hadronic tops
  } // isTTbar
  else { //is QCD
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ )
    {
      reco::PFJetRef iJet( jets, iJ );
      if ( std::abs(iJet->pt()) < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      
      if ( debug ) std::cout << " >> jet[" << iJ << "]Pt:" << iJet->pt() << " jetE:" << iJet->energy() << " jetM:" << iJet->mass() << std::endl;

      vJetIdxs.push_back(iJ);
      nJet++;
      if ( (nJets_ > 0) && (nJet >= nJets_) ) break;
    }
  } // is QCD

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

  seljet_genpart_collid.clear();
  seljet_genpart_pdgid.clear();
  seljet_genpart_charge.clear();

  seljet_genpart_px.clear();
  seljet_genpart_py.clear();
  seljet_genpart_pz.clear();
  seljet_genpart_energy.clear();

  seljet_genpart_status.clear();

  seljet_genpart_motherpdgid.clear();
  seljet_genpart_dau1pdgid.clear();
  seljet_genpart_dau2pdgid.clear();

  seljet_px.clear();
  seljet_py.clear();
  seljet_pz.clear();
  seljet_energy.clear();

  seljet_pfcand_px.clear();
  seljet_pfcand_py.clear();
  seljet_pfcand_pz.clear();
  seljet_pfcand_energy.clear();
  seljet_pfcand_type.clear();


  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByLabel( edm::InputTag("genParticles") , genparticles);

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


    seljet_px.push_back(thisJet->px());
    seljet_py.push_back(thisJet->py());
    seljet_pz.push_back(thisJet->pz());
    seljet_energy.push_back(thisJet->energy());

    std::vector<float> pfcand_px;
    std::vector<float> pfcand_py;
    std::vector<float> pfcand_pz;
    std::vector<float> pfcand_energy;
    std::vector<int> pfcand_type;

    for(auto & pfcand : thisJet->getPFConstituents())
    {

      pfcand_px.push_back(pfcand->px());
      pfcand_py.push_back(pfcand->py());
      pfcand_pz.push_back(pfcand->pz());
      pfcand_energy.push_back(pfcand->energy());
      pfcand_type.push_back((int)pfcand->particleId());

    }

    seljet_pfcand_px.push_back(pfcand_px);
    seljet_pfcand_py.push_back(pfcand_py);
    seljet_pfcand_pz.push_back(pfcand_pz);
    seljet_pfcand_energy.push_back(pfcand_energy);
    seljet_pfcand_type.push_back(pfcand_type);

    TLorentzVector TLVJet(thisJet->px(),thisJet->py(),thisJet->pz(),thisJet->energy());
    double cosTheta = TLVJet.CosTheta();
      if (cosTheta*cosTheta >=0)
        TLVJet.SetPx(0.0001);

    std::vector<int> genpart_collid;
    std::vector<int> genpart_pdgid;
    std::vector<int> genpart_charge;

    std::vector<float> genpart_px;
    std::vector<float> genpart_py;
    std::vector<float> genpart_pz;
    std::vector<float> genpart_energy;

    std::vector<int> genpart_status;

    std::vector<int> genpart_motherpdgid;
    std::vector<int> genpart_dau1pdgid;
    std::vector<int> genpart_dau2pdgid;

    std::vector<reco::GenParticle>::const_iterator genpartIterator      = (genparticles.product())->begin();
    std::vector<reco::GenParticle>::const_iterator genpartIteratorEnd   = (genparticles.product())->end();
    for ( ; genpartIterator != genpartIteratorEnd; genpartIterator++)
    {

      TLorentzVector TLVgenpart(genpartIterator->px(),genpartIterator->py(),genpartIterator->pz(),genpartIterator->energy());
      cosTheta = TLVgenpart.CosTheta();
      if (cosTheta*cosTheta >=0)
        TLVgenpart.SetPx(0.0001);
      
      
      if (TLVJet.DeltaR(TLVgenpart)<0.8)
      {
        genpart_collid.push_back(genpartIterator->collisionId());
        genpart_pdgid.push_back(genpartIterator->pdgId());
        genpart_charge.push_back(genpartIterator->charge());

        genpart_px.push_back(genpartIterator->px());
        genpart_py.push_back(genpartIterator->py());
        genpart_pz.push_back(genpartIterator->pz());
        genpart_energy.push_back(genpartIterator->energy());

        genpart_status.push_back(genpartIterator->status());

        if (genpartIterator->numberOfMothers()>0)
        {
          genpart_motherpdgid.push_back(genpartIterator->mother(0)->pdgId());
        }
        else
        {
          genpart_motherpdgid.push_back(-9999);
        }

        switch (genpartIterator->numberOfDaughters())
        {
          case 0:
            genpart_dau1pdgid.push_back(-9999);
            genpart_dau2pdgid.push_back(-9999);
          break;

          case 1:
            genpart_dau1pdgid.push_back(genpartIterator->daughter(0)->pdgId());
            genpart_dau2pdgid.push_back(-9999);
          break;

          default:
            genpart_dau1pdgid.push_back(genpartIterator->daughter(0)->pdgId());
            genpart_dau2pdgid.push_back(genpartIterator->daughter(1)->pdgId());
          break;
        }//switch
      }//DeltaR condition
    }//genpart loop
    seljet_genpart_collid.push_back(genpart_collid);
    seljet_genpart_pdgid.push_back(genpart_pdgid);
    seljet_genpart_charge.push_back(genpart_charge);

    seljet_genpart_px.push_back(genpart_px);
    seljet_genpart_py.push_back(genpart_py);
    seljet_genpart_pz.push_back(genpart_pz);
    seljet_genpart_energy.push_back(genpart_energy);

    seljet_genpart_status.push_back(genpart_status);

    seljet_genpart_motherpdgid.push_back(genpart_motherpdgid);
    seljet_genpart_dau1pdgid.push_back(genpart_dau1pdgid);
    seljet_genpart_dau2pdgid.push_back(genpart_dau2pdgid);

  }//jet loop

} // fillEvtSel_jet_dijet()
