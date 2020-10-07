#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

struct gen_obj {
  unsigned int idx;
  double pt;
};

std::vector<gen_obj> vAs;
std::vector<gen_obj> vHs;
std::vector<unsigned int> vGenAIdxs;
std::vector<unsigned int> vGenHIdxs;
// Initialize branches _____________________________________________________//
void SCRegressor::branchesH2aaSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("mHgen",     &mHgen_);
  //tree->Branch("FC_inputs", &vFC_inputs_);
  //tree->Branch("hltAccept", &hltAccept_);
  tree->Branch("Gen_E",    &vGen_E_);
  tree->Branch("A_mass",    &vA_mass_);
  tree->Branch("A_DR",      &vA_DR_);
  tree->Branch("AA_DR",     &vAA_DR_);
  tree->Branch("A_E",       &vA_E_);
  tree->Branch("A_pT",      &vA_pT_);
  tree->Branch("A_eta",     &vA_eta_);
  tree->Branch("A_phi",     &vA_phi_);
  tree->Branch("A_recoIdx", &vA_recoIdx_);
  tree->Branch("dPhi", &vdPhi_);
  tree->Branch("dEta", &vdEta_);
  hdPhidEta = fs->make<TH2F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma)",
              14, 0., 14*5*0.0174, 14, 0.,14*5.*0.0174);
}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runH2aaSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  vAs.clear();
  vHs.clear();
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vH;
  // std::cout<<std::endl<<std::endl;
  // std::cout<<"No.of Gen particles = "<<genParticles->size()<<std::endl;
   for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
      
      reco::GenParticleRef iGen( genParticles, iG );
     
     // std::cout << " >> pdgId:"<< iGen->pdgId() << " status:" << iGen->status() << " pt:" <<std::abs(iGen->pt())<<" nDaughters:" << iGen->numberOfDaughters() << " eta:" << iGen->eta() <<" phi:" << iGen->phi() <<std::endl;
   
    if ( std::abs(iGen->pdgId()) != 25 && std::abs(iGen->pdgId()) != 35 ) continue;
    //if ( std::abs(iGen->pdgId()) != 25) continue;
    if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 25 ) continue;
    if ( iGen->numberOfDaughters() != 2) continue;
    //if ( iGen->daughter(0)->pdgId() != 22 || iGen->daughter(1)->pdgId() != 22 ) continue;
    // std::cout<<std::endl<<std::endl;
      //std::cout << " >> pdgId:"<< iGen->pdgId() << " status:" << iGen->status() << " pt:" <<std::abs(iGen->pt())<<" nDaughters:" << iGen->numberOfDaughters() << " eta:" << iGen->eta() <<" phi:" << iGen->phi() <<std::endl;
    for ( unsigned int iD = 0; iD < iGen->numberOfDaughters(); iD++ ) {
      const reco::Candidate* iGenPho = iGen->daughter(iD);
        
	if (debug) std::cout << "  >> iD[" << iD <<"]: pdgId:" << iGenPho->pdgId() << " pt:" << iGenPho->pt() << " eta:" << iGenPho->eta() <<" phi:" << iGen->phi() <<std::endl;
    } 
/*    if ( std::abs(iGen->pdgId()) != 35 ) continue;
    if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 25 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    if ( iGen->daughter(0)->pdgId() != 22 || iGen->daughter(1)->pdgId() != 22 ) continue;
*/
    gen_obj Gen_obj = { iG, std::abs(iGen->pt()) };
    if ( std::abs(iGen->pdgId()) == 25){
    	vAs.push_back( Gen_obj );
    	vH += iGen->p4();
     }
    else if ( std::abs(iGen->pdgId()) == 35){
	vHs.push_back( Gen_obj );
			}

  } // gen particles
  //if ( vAs.size() != 3 ) return false;
  mHgen_ = vH.mass();
   
  // Sort As by pT, for abitrary N
  std::sort( vAs.begin(), vAs.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  vGenAIdxs.clear();
  vGenHIdxs.clear();
  //vPreselPhoIdxs_.clear();
  //bool keepEvt = true;
  //float ptomcut[2] = { 125./3., 125./4. };
  //std::cout<<"\n After sorting"<<std::endl;
  for ( unsigned int iG = 0; iG < vAs.size(); iG++ ) {
    reco::GenParticleRef iGen( genParticles, vAs[iG].idx );
    if (debug) std::cout << " >> pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy() << std::endl;
   if (debug)  std::cout<<"vAs[iG].idx = "<<vAs[iG].idx<<std::endl;
    //if ( std::abs(iGen->eta()) > 1.21 ) keepEvt = false;
    //if ( std::abs(iGen->pt()) < ptomcut[iG] ) keepEvt = false;
    //if ( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) > 0.15 ) keepEvt = false;
    //vPreselPhoIdxs_.push_back( vAs[iG].idx );
    vGenAIdxs.push_back( vAs[iG].idx );
  }
   for ( unsigned int iG = 0; iG < vHs.size(); iG++ ) {
   vGenHIdxs.push_back( vHs[iG].idx );
   //std::cout<<"vHs[iG].idx = "<<vHs[iG].idx<<std::endl;
   //std::cout<<vGenHIdxs<<std::endl;   
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
 // std::cout<<"Entered FillH2aaSel"<<std::endl;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
   
  vGen_E_.clear();
  vA_E_.clear();
  vA_pT_.clear();
  vA_eta_.clear();
  vA_phi_.clear();
  vA_mass_.clear();
  vA_DR_.clear();
  vAA_DR_.clear();
  vA_recoIdx_.clear();
  vdPhi_.clear();
  vdEta_.clear();
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
   if (debug)  std::cout<<" dR = "<<reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi())<<std::endl;
    dPhi = std::abs(reco::deltaPhi( iGen->daughter(0)->phi(), iGen->daughter(1)->phi() ));
    dEta = std::abs( iGen->daughter(0)->eta() - iGen->daughter(1)->eta() );
    hdPhidEta->Fill( dPhi, dEta );
    vdPhi_.push_back(dPhi);
    vdEta_.push_back(dEta);
    vGen_E_.push_back(std::abs(iGen->daughter(0)->energy()));
    vGen_E_.push_back(std::abs(iGen->daughter(1)->energy()));
    // Get index to dR-matched preselected photon
    recoDR = 2*0.04;
    recoDR_idx = -1;
    // Want vA_recoIdx_ to store vector index in vRegressPhoIdxs_
    // i.e., vRegressPhoIdxs_[0]:leading reco pho, vRegressPhoIdxs_[1]:sub-leading reco pho
    // not position in original photon collection
    for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
      PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());
      if ( dR < recoDR ) {
        recoDR = dR;
        recoDR_idx = iP;
      }
    } // reco pho
    vA_recoIdx_.push_back( recoDR_idx );

  } // gen A
   // unsigned int iH=vHs[0].idx;
   for ( unsigned int iH : vGenHIdxs ) {
      //  std::cout<<"Enter higgs"<<std::endl;
	//std::cout<<"i = "<<iH<<std::endl; 
   	reco::GenParticleRef iGen( genParticles, iH );
   	vAA_DR_.push_back( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) );
    }
} 
