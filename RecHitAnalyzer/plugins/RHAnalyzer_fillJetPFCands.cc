#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill Gen Particles ////////////////////////////////

std::vector<float> pfjet_px;
std::vector<float> pfjet_py;
std::vector<float> pfjet_pz;
std::vector<float> pfjet_energy;

std::vector<std::vector<float> > pfjet_pfcand_px;
std::vector<std::vector<float> > pfjet_pfcand_py;
std::vector<std::vector<float> > pfjet_pfcand_pz;
std::vector<std::vector<float> > pfjet_pfcand_energy;
std::vector<std::vector<int> > pfjet_pfcand_type;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesJetPFCands ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for gen particles


  tree->Branch("pfjet_px", &pfjet_px);
  tree->Branch("pfjet_py", &pfjet_py);
  tree->Branch("pfjet_pz", &pfjet_pz);
  tree->Branch("pfjet_energy", &pfjet_energy);

  tree->Branch("pfjet_pfcand_px", &pfjet_pfcand_px);
  tree->Branch("pfjet_pfcand_py", &pfjet_pfcand_py);
  tree->Branch("pfjet_pfcand_pz", &pfjet_pfcand_pz);
  tree->Branch("pfjet_pfcand_energy", &pfjet_pfcand_energy);
  tree->Branch("pfjet_pfcand_type", &pfjet_pfcand_type);

} // branches()

// Fill genparticles _________________________________________________________________//
void RecHitAnalyzer::fillJetPFCands ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  pfjet_px.clear();
  pfjet_py.clear();
  pfjet_pz.clear();
  pfjet_energy.clear();

  pfjet_pfcand_px.clear();
  pfjet_pfcand_py.clear();
  pfjet_pfcand_pz.clear();
  pfjet_pfcand_energy.clear();
  pfjet_pfcand_type.clear();

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);

  std::vector<reco::PFJet>::const_iterator pfjetIterator      = (jets.product())->begin();
  std::vector<reco::PFJet>::const_iterator pfjetIteratorEnd   = (jets.product())->end();

  for ( ; pfjetIterator != pfjetIteratorEnd; pfjetIterator++)
  {

    pfjet_px.push_back(pfjetIterator->px());
    pfjet_py.push_back(pfjetIterator->py());
    pfjet_pz.push_back(pfjetIterator->pz());
    pfjet_energy.push_back(pfjetIterator->energy());

    std::vector<float> pfcand_px;
    std::vector<float> pfcand_py;
    std::vector<float> pfcand_pz;
    std::vector<float> pfcand_energy;
    std::vector<int> pfcand_type;

    for(auto & pfcand : pfjetIterator->getPFConstituents())
    {

      pfcand_px.push_back(pfcand->px());
      pfcand_py.push_back(pfcand->py());
      pfcand_pz.push_back(pfcand->pz());
      pfcand_energy.push_back(pfcand->energy());
      pfcand_type.push_back((int)pfcand->particleId());

    }

    pfjet_pfcand_px.push_back(pfcand_px);
    pfjet_pfcand_py.push_back(pfcand_py);
    pfjet_pfcand_pz.push_back(pfcand_pz);
    pfjet_pfcand_energy.push_back(pfcand_energy);
    pfjet_pfcand_type.push_back(pfcand_type);


  }


} // fillEB()

