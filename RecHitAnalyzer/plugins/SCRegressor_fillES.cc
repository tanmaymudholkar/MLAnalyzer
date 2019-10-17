#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Fill ES rec hits /////////////////////////////////////////
// For each endcap, store event rechits in a vector of length
// equal to number of crystals per endcap (ix:100 x iy:100)

static const int nSTRIP = ESDetId::ISTRIP_MAX;// 32
static const int nPLANE = ESDetId::PLANE_MAX;// 2
static const int nZ = ESDetId::IZ_NUM;// 2
static const int nXY = ESDetId::IX_MAX; // 40
static const int nXY_STRIP = nXY*nSTRIP; // 40*32 = 1280

TProfile2D *hES_energy[nPLANE][nZ];
TProfile2D *hES_energy_sensor[nPLANE][nZ];
std::vector<float> vESx_energy_[nZ];
std::vector<float> vESy_energy_[nZ];

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
      // ES_F: segmented along x
      if ( ip == 0 ) {
        tree->Branch(hname, &vESx_energy_[iz]);
      }
      // ES_R: segmented along y
      else if ( ip == 1 ) {
        tree->Branch(hname, &vESy_energy_[iz]);
      }

      // Histograms for monitoring
      sprintf(hname, "%s_energy_sensor",hid);
      sprintf(htitle,"E(sensorX,sensorY);sensorX;sensorY");
      hES_energy_sensor[ip][iz] = fs->make<TProfile2D>(hname, htitle,
          nXY, 0, nXY,
          nXY, 0, nXY  );
      sprintf(hname, "%s_energy",hid);
      sprintf(htitle,"E(stripX,stripY);stripX;stripY");
      hES_energy[ip][iz] = fs->make<TProfile2D>(hname, htitle,
          nXY_STRIP, 0, nXY_STRIP,
          nXY_STRIP, 0, nXY_STRIP  );
    }
  } // iz

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
    vESx_energy_[iz].assign( nXY_STRIP*nXY, 0. );
    vESy_energy_[iz].assign( nXY*nXY_STRIP, 0. );
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

    pos  = caloGeom->getPosition( esId );

    //std::cout << "ix,iy,iz,istrip,iplane:" << ix_ << "," << iy_ << "," << iz_ << "," << istrip_ << "," << ip_
    //          << " | x,y:" << pos.x() << "," << pos.y() << std::endl;

    // Fill histograms for monitoring
    hES_energy_sensor[ip_][iz_]->Fill( ix_, iy_, energy_ );
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
    }
    //*/
    ///*

    // Fill vectors for images
    // Create hashed Index: maps from [iy][ix] -> [idx_]
    idx_ = stripY_*nXY_STRIP + stripX_;
    // ES_F: segmented along x
    if ( ip_ == 0 ) {
      vESx_energy_[iz_][idx_] = energy_;
    }
    // ES_R: segmented along y
    else if ( ip_ == 1 ) {
      vESy_energy_[iz_][idx_] = energy_;
    }
    //*/

  } // ES rechits

} // fillES()
