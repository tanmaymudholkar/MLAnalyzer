#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Fill ES rec hits /////////////////////////////////////////
// For each endcap, store event rechits in a vector of length
// equal to number of crystals per endcap (ix:100 x iy:100)

static const int nSTRIP = ESDetId::ISTRIP_MAX;// 32
static const int nPLANE = ESDetId::PLANE_MAX;// 2
static const int nZ = ESDetId::IZ_NUM;// 2
static const int nXY = ESDetId::IX_MAX; // 40
static const int nXY_STRIP = nXY*nSTRIP; // 40*32 = 1280
static const int search_window = 5;

TProfile2D *hES_energy[nPLANE][nZ];
TProfile2D *hES_energy_sensor[nPLANE][nZ];
TH2F *hEvt_ES_energy_sensor[nPLANE][nZ];
std::vector<float> vES_energy_[nPLANE][nZ];
std::vector<float> seedIx_, seedIy_;

// Initialize branches _____________________________________________________//
void SCRegressor::branchesES ( TTree* tree, edm::Service<TFileService> &fs ) {

  char hid[10], hname[50], htitle[50];
  for ( int iz(0); iz < nZ; iz++ ) {
    const char *zside_str = (iz > 0) ? "p" : "m";
    for ( int ip(0); ip < nPLANE; ip++ ) {

      const char *plane_str = (ip > 0) ? "Y" : "X"; // X:F:0, Y:R:1
      sprintf(hid, "ES%s%s",zside_str, plane_str);

      // Branches for images
      sprintf(hname, "%s_energy",hid);
      tree->Branch(hname, &vES_energy_[ip][iz]);

      // Histograms for monitoring
      sprintf(hname, "%s_energy_sensor",hid);
      sprintf(htitle,"E(sensorX,sensorY);sensorX;sensorY");
      hES_energy_sensor[ip][iz] = fs->make<TProfile2D>(hname, htitle,
          nXY, 0, nXY,
          nXY, 0, nXY  );
      sprintf(hname, "evt_%s_energy_sensor",hid);
      hEvt_ES_energy_sensor[ip][iz] = new TH2F(hname, htitle,
          nXY, 0, nXY,
          nXY, 0, nXY  );
      sprintf(hname, "%s_energy",hid);
      sprintf(htitle,"E(stripX,stripY);stripX;stripY");
      hES_energy[ip][iz] = fs->make<TProfile2D>(hname, htitle,
          nXY_STRIP, 0, nXY_STRIP,
          nXY_STRIP, 0, nXY_STRIP  );
    } // ip
  } // iz
  tree->Branch("seed_ix", &seedIx_);
  tree->Branch("seed_iy", &seedIy_);

} // branchesES()

// Fill ES rechits _________________________________________________________________//
void SCRegressor::fillES ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int ix_, iy_, iz_, idx_; // rows:iy, columns:ix
  int istrip_, ip_, ibroad_;
  int stripX_, stripY_;
  float energy_;
  GlobalPoint pos;

  ///*
  for ( int iz(0); iz < nZ; iz++ ) {
    for ( int ip(0); ip < nPLANE; ip++ ) {
      vES_energy_[ip][iz].assign( nXY_STRIP*nXY, 0. );
      hEvt_ES_energy_sensor[ip][iz]->Reset();
    }
  }
  //*/

  edm::Handle<ESRecHitCollection> ESRecHitsH_;
  iEvent.getByToken( ESRecHitCollectionT_, ESRecHitsH_ );

  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom;
  caloGeom = caloGeomH_.product();
  //const CaloCellGeometry::RepCorners& repCorners;

  /*
  // Geometry test
  std::cout << "kSize:" << ESDetId::kSizeForDenseIndexing << std::endl;
  for ( int is(0); is < ESDetId::kSizeForDenseIndexing; is++ ) {
    ESDetId esId = ESDetId::detIdFromDenseIndex(is);
    ix_ = esId.six() - 1;
    iy_ = esId.siy() - 1;
    iz_ = (esId.zside() > 0) ? 1 : 0;
    istrip_ = esId.strip() - 1;
    ip_ = esId.plane() - 1;
    //hES_energy_sensor[ip_][iz_]->Fill( ix_, iy_, istrip_ );
    // ES_F: segmented along x
    if ( ip_ == 0 ) {
      stripX_ = ix_*nSTRIP + istrip_;
      stripY_ = iy_*nSTRIP;
      for ( int is(0); is < nSTRIP; is++ ) {
        hES_energy[ip_][iz_]->Fill( stripX_, stripY_+is, stripX_ );
      }
    }
    // ES_R: segmented along y
    else if ( ip_ == 1 ) {
      stripX_ = ix_*nSTRIP;
      stripY_ = iy_*nSTRIP + istrip_;
      for ( int is(0); is < nSTRIP; is++ ) {
        hES_energy[ip_][iz_]->Fill( stripX_+is, stripY_, stripY_ );
      }
    }
  }
  */

  // Fill ES rechits
  for ( ESRecHitCollection::const_iterator iRHit = ESRecHitsH_->begin();
        iRHit != ESRecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    //if ( energy_ <= zs ) continue;
    if ( energy_ <= 0. ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    ESDetId esId( iRHit->id() );
    ix_ = esId.six() - 1;
    iy_ = esId.siy() - 1;
    iz_ = (esId.zside() > 0) ? 1 : 0;
    istrip_ = esId.strip() - 1;
    ip_ = esId.plane() - 1;
    pos  = caloGeom->getPosition( esId );

    /*
    // Get REP corners of cell Id
    const auto repCorners = caloGeom->getGeometry(esId)->getCornersREP();
    // Get min,max phi,eta at plane closest to IP
    // See illustration at bottom
    isBoundary = false;
    minEta_ = repCorners[2].eta();
    maxEta_ = repCorners[0].eta();
    minPhi_ = repCorners[2].phi();
    maxPhi_ = repCorners[0].phi();
    */

    //std::cout << "ix,iy,iz,istrip,iplane:" << ix_ << "," << iy_ << "," << iz_ << "," << istrip_ << "," << ip_
    //          << " | x,y,z:" << pos.x() << "," << pos.y() << "," << pos.z() << std::endl;

    // Fill histograms for monitoring
    hES_energy_sensor[ip_][iz_]->Fill( ix_, iy_, energy_ );
    hEvt_ES_energy_sensor[ip_][iz_]->Fill( ix_, iy_, energy_ );
    ///*
    // ES_F: segmented along x
    if ( ip_ == 0 ) {
      stripX_ = ix_*nSTRIP + istrip_;
      //stripY_ = iy_*nSTRIP;
      stripY_ = iy_;
      // Broadcast along unsegmented direction
      for ( int is(0); is < nSTRIP; is++ ) {
        ibroad_ = (iy_*nSTRIP) + is;
        hES_energy[ip_][iz_]->Fill( stripX_, ibroad_, energy_ );
      }
      idx_ = stripY_*nXY_STRIP + stripX_;
    }
    // ES_R: segmented along y
    else if ( ip_ == 1 ) {
      stripX_ = ix_;
      //stripX_ = ix_*nSTRIP;
      stripY_ = iy_*nSTRIP + istrip_;
      // Broadcast along unsegmented direction
      for ( int is(0); is < nSTRIP; is++ ) {
        //hES_energy[ip_][iz_]->Fill( stripX_+is, stripY_, energy_/float(nSTRIP) );
        ibroad_ = (ix_*nSTRIP) + is;
        hES_energy[ip_][iz_]->Fill( ibroad_, stripY_, energy_ );
      }
      idx_ = stripY_*nXY + stripX_;
    }

    // Fill vectors for images
    vES_energy_[ip_][iz_][idx_] = energy_;
    //*/

  } // ES rechits

  ///*
  const double zESF = 303.846;
  const double zESR = 308.306;
  double phi, theta, radius, z, rcyl;
  GlobalPoint point;
  int siX_, siY_, seedIx, seedIy;
  float binE, seedE;
  const CaloSubdetectorGeometry* esGeom = caloGeom->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  std::cout << "nPho:" << vRegressPhoIdxs_.size() << std::endl;
  for ( int iP : vRegressPhoIdxs_ ) {

    PhotonRef iPho( photons, iP );

    phi = iPho->phi();
    theta = 2.*std::atan( exp(-iPho->eta()) );
    radius = zESF / std::abs(std::cos(theta));

    z = radius * std::cos(theta);
    rcyl = radius * std::sin(theta);
    GlobalPoint point(rcyl * std::cos(phi), rcyl * std::sin(phi), z);
    const DetId detId = esGeom->getClosestCell(point);

    ESDetId esId = ESDetId( detId );
    ix_ = esId.six() - 1;
    iy_ = esId.siy() - 1;
    iz_ = (esId.zside() > 0) ? 1 : 0;
    istrip_ = esId.strip() - 1;
    ip_ = esId.plane() - 1;
    pos  = caloGeom->getPosition( esId );

    //std::cout << "ix,iy,iz,istrip,iplane:" << ix_ << "," << iy_ << "," << iz_ << "," << istrip_ << "," << ip_
    //          << " | x,y,z:" << pos.x() << "," << pos.y() << "," << pos.z() << std::endl;

    seedE = 0.;
    seedIx = -1;
    seedIy = -1;
    // Look for the most energetic HBHE tower deposit within a search window
    for ( int siX = 0; siX < search_window; siX++ ) {

      siX_ = ix_ - (search_window/2)+siX;
      if ( siX_ > nXY-1 ) continue;
      if ( siX_ < 0 ) continue;

      for ( int siY = 0; siY < search_window; siY++ ) {

        siY_ = iy_ - (search_window/2)+siY;
        if ( siY_ > nXY-1 ) continue;
        if ( siY_ < 0 ) continue;

        for ( int siP = 0; siP < nPLANE; siP++ ) {

          binE = hES_energy_sensor[siP][iz_]->GetBinContent( siX_+1, siY_+1 );
          //std::cout << "ix,iy,iz,iplane:" << siX_ << "," << siY_ << "," << iz_ << "," << siP << " | energy:" << binE << std::endl;

          if ( binE <= seedE ) continue;
          seedE = binE;
          seedIx = siX_;
          seedIy = siY_;

        } // siP
      } // siY
    } // siX
    seedIx_.push_back( seedIx );
    seedIy_.push_back( seedIy );

  } // photons
  //*/

} // fillES()
