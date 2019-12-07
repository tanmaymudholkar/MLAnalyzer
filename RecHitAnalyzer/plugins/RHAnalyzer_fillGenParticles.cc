#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill Gen Particles ////////////////////////////////

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

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesGenParticles ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for gen particles

  tree->Branch("genpart_collid", &genpart_collid);
  tree->Branch("genpart_pdgid", &genpart_pdgid);
  tree->Branch("genpart_charge", &genpart_charge);

  tree->Branch("genpart_px", &genpart_px);
  tree->Branch("genpart_py", &genpart_py);
  tree->Branch("genpart_pz", &genpart_pz);
  tree->Branch("genpart_energy", &genpart_energy);

  tree->Branch("genpart_status", &genpart_status);

  tree->Branch("genpart_motherpdgid", &genpart_motherpdgid);
  tree->Branch("genpart_dau1pdgid", &genpart_dau1pdgid);
  tree->Branch("genpart_dau2pdgid", &genpart_dau2pdgid);

} // branches()

// Fill genparticles _________________________________________________________________//
void RecHitAnalyzer::fillGenParticles ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  genpart_collid.clear();
  genpart_pdgid.clear();
  genpart_charge.clear();

  genpart_px.clear();
  genpart_py.clear();
  genpart_pz.clear();
  genpart_energy.clear();

  genpart_status.clear();

  genpart_motherpdgid.clear();
  genpart_dau1pdgid.clear();
  genpart_dau2pdgid.clear();

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByLabel( edm::InputTag("genParticles") , genparticles);
  std::vector<reco::GenParticle>::const_iterator genpartIterator      = (genparticles.product())->begin();
  std::vector<reco::GenParticle>::const_iterator genpartIteratorEnd   = (genparticles.product())->end();

  for ( ; genpartIterator != genpartIteratorEnd; genpartIterator++)
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
      
    }

  }


} // fillEB()

