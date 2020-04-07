#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Fill EB rec hits ////////////////////////////////
// Store event rechits in a vector of length equal
// to number of crystals in EB (ieta:170 x iphi:360)

//TProfile2D *hEB_energy;
//TProfile2D *hEB_time;
//std::vector<float> vEB_energy_;
//std::vector<float> vEB_time_;

// Initialize branches _____________________________________________________//
void SCRegressor::branchesEB ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  //tree->Branch("EB_energy", &vEB_energy_);
  tree->Branch("EB_energyT", &vEB_energyT_);
  tree->Branch("EB_energyZ", &vEB_energyZ_);
  tree->Branch("EB_energyErr", &vEB_energyErr_);
  //tree->Branch("EB_time",   &vEB_time_);

  // Histograms for monitoring
  hEB_energy = fs->make<TProfile2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hEB_energyErr = fs->make<TProfile2D>("EB_energyErr", "Err(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hEB_time = fs->make<TProfile2D>("EB_time", "t(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );

  hSC_energyErr = fs->make<TProfile>("SC_energyErr", "Err vs E;energy;Err",
      100, 0., 100.);

} // branchesEB()

// Fill EB rechits _________________________________________________________________//
void SCRegressor::fillEB ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  float energy_;

  //vEB_energy_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vEB_energyT_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vEB_energyZ_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vEB_energyErr_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  //vEB_time_.assign( EBDetId::kSizeForDenseIndexing, 0. );

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByToken( EBRecHitCollectionT_, EBRecHitsH_);

  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();

  //std::cout << "EB" << std::endl;
  // Fill EB rechits
  for ( EcalRecHitCollection::const_iterator iRHit = EBRecHitsH_->begin();
        iRHit != EBRecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ <= zs ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi() - 1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    // Fill histograms for monitoring
    hEB_energy->Fill( iphi_,ieta_,energy_ );
    hEB_energyErr->Fill( iphi_,ieta_,iRHit->energyError() );
    hEB_time->Fill( iphi_,ieta_,iRHit->time() );
    // Get Hashed Index: provides convenient
    // index mapping from [ieta][iphi] -> [idx]
    idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_

    // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
    //auto cell = caloGeom->getGeometry(ebId);
    auto pos = caloGeom->getPosition(ebId);

    // Fill vectors for images
    //vEB_energy_[idx_] = energy_;
    vEB_energyT_[idx_] = energy_/TMath::CosH(pos.eta());
    vEB_energyZ_[idx_] = energy_*std::abs(TMath::TanH(pos.eta()));
    vEB_energyErr_[idx_] = iRHit->energyError();
    //vEB_time_[idx_] = iRHit->time();

    //std::cout << "idx,ieta,iphi,E:" <<idx_<<","<< ieta_ << "," << iphi_ << "," << iRHit->energy() << std::endl;

  } // EB rechits

} // fillEB()
