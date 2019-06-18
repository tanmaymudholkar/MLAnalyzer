#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesSCaod ( TTree* tree, edm::Service<TFileService> &fs )
{
  hSCaod_energy = fs->make<TProfile2D>("SCaod_energy", "E(i#phi,i#eta);iphi;ieta",
      crop_size, 0, crop_size,
      crop_size, 0, crop_size );
  hSCaod_time = fs->make<TProfile2D>("SCaod_time", "t(i#phi,i#eta);iphi;ieta",
      crop_size, 0, crop_size,
      crop_size, 0, crop_size );

  RHTree->Branch("SCaod_energy",  &vSCaod_energy_);
  RHTree->Branch("SCaod_energyT", &vSCaod_energyT_);
  RHTree->Branch("SCaod_energyZ", &vSCaod_energyZ_);
  RHTree->Branch("SCaod_time",    &vSCaod_time_);

} // branchesSCaod()

// Fill SCaod rechits _________________________________________________________________//
void SCRegressor::fillSCaod ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(AODEBRecHitCollectionT_, EBRecHitsH);

  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();

  vSCaod_energy_.clear();
  vSCaod_energyT_.clear();
  vSCaod_energyZ_.clear();
  vSCaod_time_.clear();
  std::vector<float> SCaod_energy;
  std::vector<float> SCaod_energyT;
  std::vector<float> SCaod_energyZ;
  std::vector<float> SCaod_time;

  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  int iphi_shift, ieta_shift;
  int iphi_crop, ieta_crop;
  for ( unsigned int iP(0); iP < nPho; iP++ ) {

    SCaod_energy.assign(crop_size*crop_size,0.);
    SCaod_energyT.assign(crop_size*crop_size,0.);
    SCaod_energyZ.assign(crop_size*crop_size,0.);
    SCaod_time.assign(crop_size*crop_size,0.);

    iphi_shift = vIphi_Emax_[iP] - 15;
    ieta_shift = vIeta_Emax_[iP] - 15;
    if ( debug ) std::cout << " >> Doing pho img: iphi_Emax,ieta_Emax: " << vIphi_Emax_[iP] << ", " << vIeta_Emax_[iP] << std::endl;

    for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
        iRHit != EBRecHitsH->end();
        ++iRHit) {

      if ( iRHit->energy() < zs ) continue;

      // Convert detector coordinates to ordinals
      EBDetId ebId( iRHit->id() );
      iphi_ = ebId.iphi()-1; // [0,...,359]
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,0,...,84]
      ieta_ += EBDetId::MAX_IETA; // [0,...,169]

      // Convert to [0,...,31][0,...,31]
      ieta_crop = ieta_ - ieta_shift;
      iphi_crop = iphi_ - iphi_shift;
      if ( iphi_crop >= EBDetId::MAX_IPHI ) iphi_crop = iphi_crop - EBDetId::MAX_IPHI; // get wrap-around hits
      if ( iphi_crop < 0 ) iphi_crop = iphi_crop + EBDetId::MAX_IPHI; // get wrap-around hits

      if ( ieta_crop < 0 || ieta_crop > crop_size-1 ) continue;
      if ( iphi_crop < 0 || iphi_crop > crop_size-1 ) continue;

      // Convert to [0,...,32*32-1]
      idx_ = ieta_crop*crop_size + iphi_crop;

      // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
      //auto cell = caloGeom->getGeometry(ebId);
      auto pos = caloGeom->getPosition(ebId);

      // Fill branch arrays
      SCaod_energy[idx_] = iRHit->energy();
      //SCaod_energyT[idx_] = iRHit->energy()/TMath::CosH(vPho_eta_);
      SCaod_energyT[idx_] = iRHit->energy()/TMath::CosH(pos.eta());
      SCaod_energyZ[idx_] = iRHit->energy()*std::abs(TMath::TanH(pos.eta()));
      SCaod_time[idx_] = iRHit->time();
      /*
      vEB_SCenergy_[ebId.hashedIndex()] = iRHit->energy();
      */
      //std::cout << " >> " << iP << ": iphi_,ieta_,E: " << iphi_crop << ", " << ieta_crop << ", " << iRHit->energy() << std::endl; 
      //std::cout << "idx,ieta,iphi,E:" <<idx_<<","<< ieta_crop << "," << iphi_crop << "," << iRHit->energy() << std::endl; 

      // Fill histograms to monitor cumulative distributions
      hSCaod_energy->Fill( iphi_crop,ieta_crop,iRHit->energy() );
      hSCaod_time->Fill( iphi_crop,ieta_crop,iRHit->time() );

    } // EB rechits
    vSCaod_energy_.push_back( SCaod_energy );
    vSCaod_energyT_.push_back( SCaod_energyT );
    vSCaod_energyZ_.push_back( SCaod_energyZ );
    vSCaod_time_.push_back( SCaod_time );

  } // photons
  /*
  for(auto const& e:vSCaod_energy_) {
    std::cout << "array" << std::endl;
    for ( int i =0;i < 32*32;i++) {
      //if (e[i] > 0.) std::cout << "idx,ieta,iphi,E:" <<i<<","<< (int)(i/32) << ","<< i % 32  <<","<< e[i] << std::endl;
    };
  };
  */

} // fillSCaod()
