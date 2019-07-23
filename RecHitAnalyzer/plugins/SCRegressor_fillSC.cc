#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesSC ( TTree* tree, edm::Service<TFileService> &fs )
{
  hSC_energy = fs->make<TProfile2D>("SC_energy", "E(i#phi,i#eta);iphi;ieta",
      crop_size, 0, crop_size,
      crop_size, 0, crop_size );
  hSC_time = fs->make<TProfile2D>("SC_time", "t(i#phi,i#eta);iphi;ieta",
      crop_size, 0, crop_size,
      crop_size, 0, crop_size );

  RHTree->Branch("SC_energy",  &vSC_energy_);
  RHTree->Branch("SC_energyT", &vSC_energyT_);
  RHTree->Branch("SC_energyZ", &vSC_energyZ_);
  RHTree->Branch("SC_time",    &vSC_time_);

} // branchesSC()

// Fill SC rechits _________________________________________________________________//
void SCRegressor::fillSC ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);

  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();

  vSC_energy_.clear();
  vSC_energyT_.clear();
  vSC_energyZ_.clear();
  vSC_time_.clear();
  std::vector<float> SC_energy;
  std::vector<float> SC_energyT;
  std::vector<float> SC_energyZ;
  std::vector<float> SC_time;

  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  int iphi_shift, ieta_shift;
  int iphi_crop, ieta_crop;
  for ( unsigned int iP(0); iP < nPho; iP++ ) {

    SC_energy.assign(crop_size*crop_size,0.);
    SC_energyT.assign(crop_size*crop_size,0.);
    SC_energyZ.assign(crop_size*crop_size,0.);
    SC_time.assign(crop_size*crop_size,0.);

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
      SC_energy[idx_] = iRHit->energy();
      //SC_energyT[idx_] = iRHit->energy()/TMath::CosH(vPho_eta_);
      SC_energyT[idx_] = iRHit->energy()/TMath::CosH(pos.eta());
      SC_energyZ[idx_] = iRHit->energy()*std::abs(TMath::TanH(pos.eta()));
      SC_time[idx_] = iRHit->time();
      /*
      vEB_SCenergy_[ebId.hashedIndex()] = iRHit->energy();
      */
      //std::cout << " >> " << iP << ": iphi_,ieta_,E: " << iphi_crop << ", " << ieta_crop << ", " << iRHit->energy() << std::endl; 
      //std::cout << "idx,ieta,iphi,E:" <<idx_<<","<< ieta_crop << "," << iphi_crop << "," << iRHit->energy() << std::endl; 

      // Fill histograms to monitor cumulative distributions
      hSC_energy->Fill( iphi_crop,ieta_crop,iRHit->energy() );
      hSC_time->Fill( iphi_crop,ieta_crop,iRHit->time() );

    } // EB rechits
    vSC_energy_.push_back( SC_energy );
    vSC_energyT_.push_back( SC_energyT );
    vSC_energyZ_.push_back( SC_energyZ );
    vSC_time_.push_back( SC_time );

  } // photons
  /*
  for(auto const& e:vSC_energy_) {
    std::cout << "array" << std::endl;
    for ( int i =0;i < 32*32;i++) {
      //if (e[i] > 0.) std::cout << "idx,ieta,iphi,E:" <<i<<","<< (int)(i/32) << ","<< i % 32  <<","<< e[i] << std::endl;
    };
  };
  */

} // fillSC()
