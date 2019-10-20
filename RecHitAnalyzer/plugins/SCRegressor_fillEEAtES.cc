#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Fill EE rec hits /////////////////////////////////////////
// For each endcap, store event rechits in a vector of length
// equal to number of crystals per endcap (ix:100 x iy:100)

TH2F *hEEatES_energy[nZ];
std::vector<float> vEEatES_energy[nZ];
const int EEAtES_search_window = 3*nSTRIP;

// Initialize branches _____________________________________________________________//
void SCRegressor::branchesEEatES ( TTree* tree, edm::Service<TFileService> &fs ) {

  char hname[50], htitle[50];
  for ( int iz(0); iz < nEE; iz++ ) {
    // Branches for images
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "EE_energy_ES%s",zside);
    tree->Branch(hname,        &vEEatES_energy[iz]);

    // Histograms for monitoring
    hEEatES_energy[iz] = fs->make<TH2F>(hname, htitle,
        nXY_STRIP, 0, nXY_STRIP,
        nXY_STRIP, 0, nXY_STRIP  );
  } // iz

} // branchesEEatES()

ESDetId ESId_from_EtaPhi( float& eta, float& phi, const CaloGeometry* caloGeom ) {

  //const double zESF = 303.846;
  const double zESR = 308.306;
  const CaloSubdetectorGeometry* esGeom = caloGeom->getSubdetectorGeometry( DetId::Ecal, EcalPreshower );

  double theta = 2.*std::atan( exp(-eta) );
  double radius = zESR / std::abs(std::cos(theta));

  double z = radius * std::cos(theta);
  double rcyl = radius * std::sin(theta);
  GlobalPoint point(rcyl * std::cos(phi), rcyl * std::sin(phi), z);
  ESDetId esId( esGeom->getClosestCell(point) );

  return esId;

} // ESId_from_EtaPhi

/*
float minXYcorner( EEDetId &eeId, SCRegressor_fillEEAtES.cc ) {
  const auto xyzCorners = caloGeom->getGeometry(eeId)->getCorners();
  // Get min,max phi,eta at plane closest to IP
  // See illustration at bottom
  //isBoundary = false;
  xCorners = { xyzCorners[0].x(), xyzCorners[1].x(), xyzCorners[2].x(), xyzCorners[3].x() };
  yCorners = { xyzCorners[0].y(), xyzCorners[1].y(), xyzCorners[2].y(), xyzCorners[3].y() };
  //minY_ = xyzCorners[2].y();
  //maxY_ = xyzCorners[0].y();
  minX_ = *std::min_element(xCorners.begin(), xCorners.end());
  maxX_ = *std::max_element(xCorners.begin(), xCorners.end());
  minY_ = *std::min_element(yCorners.begin(), yCorners.end());
  maxY_ = *std::max_element(yCorners.begin(), yCorners.end());
}
*/
// Fill HCAL rechits at EB/EE ______________________________________________________________//
void SCRegressor::fillEEatES ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int ix_, iy_, iz_, idx_; // rows:iy, columns:ix
  int nEExtals_filled;
  float energy_;
  float eta, phi, minY_, minX_, maxY_, maxX_;
  float st_x, st_y;
  float st_minx, st_miny;
  float st_maxx, st_maxy;
  int stX_, stY_;
  int stripX_seed, stripY_seed, st_ix, st_iy, st_xstrip, st_ystrip;
  int nESstrips_in_EExtal;
  GlobalPoint pos;
  bool isBoundary;
  std::vector<int> vESXstrips_in_EExtal;
  std::vector<int> vESYstrips_in_EExtal;

  for ( int iz(0); iz < nEE; iz++ ) {
    vEEatES_energy[iz].assign( nXY_STRIP*nXY_STRIP, 0. );
  }

  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByToken( EERecHitCollectionT_, EERecHitsH_ );
  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom;
  caloGeom = caloGeomH_.product();
  //const CaloCellGeometry::RepCorners& repCorners;

  int ixs[2] = {25, 80};
  int iys[2] = {25, 80};

  std::vector<float> esCornersX{ 0., 0., 0., 0.};
  std::vector<float> esCornersY{ 0., 0., 0., 0.};
  std::vector<float> xCorners{ 0., 0., 0., 0.};
  std::vector<float> yCorners{ 0., 0., 0., 0.};
  for ( EcalRecHitCollection::const_iterator iRHit = EERecHitsH_->begin();
        iRHit != EERecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ <= 0. ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    EEDetId eeId( iRHit->id() );
    //EEDetId eeId(ixs[0], iys[0], 0);
    //ix_ = eeId.ix() - 1;
    //iy_ = eeId.iy() - 1;
    //iz_ = (eeId.zside() > 0) ? 1 : 0;
    /*
    for ( int i=0; i < 2; i++ ) {
      for ( int j=0; j < 2; j++ ) {
        for ( int k=0; k < 2; k++ ) {
          std::cout << "(ix,iy):"<< ixs[i] << ","<<iys[j] <<std::endl;
          EEDetId eeId(ixs[i], iys[j], k);
    */
    /*
    for ( int iC = 0; iC < EEDetId::kSizeForDenseIndexing; iC++ ) {
      // Store EE xtals with centers within corners of HCAL tower
      EEDetId eeId( EEDetId::unhashIndex(iC) );
    */

    // Get REP corners of cell Id
    const auto xyzCorners = caloGeom->getGeometry(eeId)->getCorners();
    // Get min,max phi,eta at plane closest to IP
    // See illustration at bottom
    //isBoundary = false;
    xCorners = { xyzCorners[0].x(), xyzCorners[1].x(), xyzCorners[2].x(), xyzCorners[3].x() };
    yCorners = { xyzCorners[0].y(), xyzCorners[1].y(), xyzCorners[2].y(), xyzCorners[3].y() };
    //minY_ = xyzCorners[2].y();
    //maxY_ = xyzCorners[0].y();
    minX_ = *std::min_element(xCorners.begin(), xCorners.end());
    maxX_ = *std::max_element(xCorners.begin(), xCorners.end());
    minY_ = *std::min_element(yCorners.begin(), yCorners.end());
    maxY_ = *std::max_element(yCorners.begin(), yCorners.end());
    //maxX_ = xyzCorners[0].x();
    //if ( minX_ > maxX_ ) isBoundary = true;
    //std::cout << "(minY:maxY) "<<minY_ <<":"<<maxY_ << ", minX:maxX " << minX_<<":"<<maxX_ <<std::endl;
    /*
    std::cout << "(zside,x0,x1,x2,x3) "<< eeId.zside() <<" : "
     << xyzCorners[0].x() << " : "
     << xyzCorners[1].x() << " : "
     << xyzCorners[2].x() << " : "
     << xyzCorners[3].x() << " : "
     << std::endl;
    std::cout << "(zside,y0,y1,y2,y3) "<< eeId.zside() <<" : "
     << xyzCorners[0].y() << " : "
     << xyzCorners[1].y() << " : "
     << xyzCorners[2].y() << " : "
     << xyzCorners[3].y() << " : "
     << std::endl;
        }
      }
    }
    */
    // Get closest ES id
    pos  = caloGeom->getPosition( eeId );
    eta = pos.eta();
    phi = pos.phi();
    ESDetId nearestESId = ESId_from_EtaPhi( eta, phi, caloGeom );
    ix_ = nearestESId.six() - 1;
    iy_ = nearestESId.siy() - 1;
    iz_ = (nearestESId.zside() > 0) ? 1 : 0;
    stripX_seed = nSTRIP*ix_ + (nSTRIP/2);
    stripY_seed = nSTRIP*iy_ + (nSTRIP/2);
    //std::cout << ESDetId::validDetId(nearestESId.strip(), nearestESId.six(), nearestESId.siy(), nearestESId.plane(), nearestESId.zside()) << std::endl;
    // Search window around closest ES id:
    nESstrips_in_EExtal = 0;
    vESXstrips_in_EExtal.clear();
    vESYstrips_in_EExtal.clear();
    for ( int stX = 0; stX < EEAtES_search_window; stX++ ) {

      stX_ = stripX_seed - (EEAtES_search_window/2)+stX;
      if ( stX_ > nXY_STRIP-1 ) continue;
      if ( stX_ < 0 ) continue;

      for ( int stY = 0; stY < EEAtES_search_window; stY++ ) {

        stY_ = stripY_seed - (EEAtES_search_window/2)+stY;
        if ( stY_ > nXY_STRIP-1 ) continue;
        if ( stY_ < 0 ) continue;

        //std::cout << "stX_,stY_,E:" << stX_ << "," << stY_ << "," << energy_ << std::endl;
        // Define fictional 1.9mm x 1.9mm strip `st` with
        // x-coords of ES_F strip and y-coords of ES_R strip
        st_ix = stX_ / nSTRIP;
        st_xstrip = stX_ % nSTRIP;
        st_iy = stY_ / nSTRIP;
        st_ystrip = stY_ % nSTRIP;
        //std::cout << "ix,xstrip, iy,ystrip:" << st_ix << "," << st_xstrip << " , "<< st_iy << "," << st_ystrip << std::endl;
        //std::cout << ESDetId::validDetId(st_xstrip+1, st_ix+1, st_iy+1, 2, iz_ > 0 ? 1 : -1 ) << std::endl;

        if ( !ESDetId::validDetId(st_xstrip+1, st_ix+1, st_iy+1, 1, iz_ > 0 ? 1 : -1) ) continue;
        ESDetId stId_x  = ESDetId(st_xstrip+1, st_ix+1, st_iy+1, 1, iz_ > 0 ? 1 : -1);
        st_x = caloGeom->getPosition( stId_x ).x(); // ES_F
        // If 1.9mm x 1.9mm centers lie outside of EE edges then do not fill ES image
        if ( st_x < minX_ || st_x > maxX_ ) continue;
        /*
        const auto esCorners_x = caloGeom->getGeometry(stId_x)->getCorners();
        esCornersX = { esCorners_x[0].x(), esCorners_x[1].x(), esCorners_x[2].x(), esCorners_x[3].x() };
        st_minx = *std::min_element(esCornersX.begin(), esCornersX.end());
        st_maxx = *std::max_element(esCornersX.begin(), esCornersX.end());
        if ( st_maxx < minX_ || st_minx > maxX_ ) continue;
        */

        if ( !ESDetId::validDetId(st_ystrip+1, st_ix+1, st_iy+1, 2, iz_ > 0 ? 1 : -1) ) continue;
        ESDetId stId_y  = ESDetId(st_ystrip+1, st_ix+1, st_iy+1, 2, iz_ > 0 ? 1 : -1);
        st_y = caloGeom->getPosition( stId_y ).y(); // ES_R
        if ( st_y < minY_ || st_y > maxY_ ) continue;
        /*
        const auto esCorners_y = caloGeom->getGeometry(stId_y)->getCorners();
        esCornersY = { esCorners_y[0].y(), esCorners_y[1].y(), esCorners_y[2].y(), esCorners_y[3].y() };
        st_miny = *std::min_element(esCornersY.begin(), esCornersY.end());
        st_maxy = *std::max_element(esCornersY.begin(), esCornersY.end());
        if ( st_maxy < minY_ || st_miny > maxY_ ) continue;
        */

        //std::cout << "ix,xstrip, iy,ystrip:" << stId_x.six() << "," << stId_x.strip() << " , "<< stId_y.siy() << "," << stId_y.strip() << std::endl;

        //std::cout << "st_x,st_y:" << st_x << "," << st_y << std::endl;

        // Fill histogram
        //hEEatES_energy[iz_]->SetBinContent( stX_, stY_, hEEatES_energy[iz_]->GetBinContent( stX_, stY_) + energy_ );
        //hEEatES_energy[iz_]->SetBinContent( stX_, stY_, hEEatES_energy[iz_]->GetBinContent( stX_, stY_) + iC );

        vESXstrips_in_EExtal.push_back( stX_ );
        vESYstrips_in_EExtal.push_back( stY_ );

      } // siY
    } // siX
    // loop over zip(stX_,stY_)
    nESstrips_in_EExtal = vESXstrips_in_EExtal.size();
    if ( nESstrips_in_EExtal == 0 ) continue;
    //std::cout << "nESstrips_in_EExtal:" << nESstrips_in_EExtal << std::endl;
    for ( int i = 0; i < nESstrips_in_EExtal; i++ ) {
      //idx_ = vESYstrips_in_EExtal[i]*nXY_STRIP + vESXstrips_in_EExtal[i];
      stX_ = vESXstrips_in_EExtal[i];
      stY_ = vESYstrips_in_EExtal[i];
      idx_ = stY_*nXY_STRIP + stX_;
      vEEatES_energy[iz_][idx_] += energy_/float(nESstrips_in_EExtal);
      hEEatES_energy[iz_]->SetBinContent( stX_, stY_, hEEatES_energy[iz_]->GetBinContent(stX_, stY_) + energy_/float(nESstrips_in_EExtal) );
    }
    // fill image[stX,stY_] vector with energy_/stX_.size();

    //}
    //break;
    /*
    // Loop over all EE crystals
    vEExtals_in_HBHEtower.clear();
    for ( int iC = 0; iC < EEDetId::kSizeForDenseIndexing; iC++ ) {
      // Store EE xtals with centers within corners of HCAL tower
      EEDetId eeId( EEDetId::unhashIndex(iC) );
      pos = caloGeom->getPosition( eeId );
      eta = pos.eta();
      phi = pos.phi();
      if ( eta < minY_ || eta > maxY_ ) continue;
      if ( !isBoundary ) {
        if ( phi < minX_ || phi > maxX_ ) continue;
      } else {
        if ( phi < minX_ && phi > maxX_ ) continue;
      }
      vEExtals_in_HBHEtower.push_back( iC );
      //std::cout << "xtal eta,phi: " << eta << " " << phi << std::endl;
      //if ( iphi > 30 || iphi < 40 ) std::cout << "xtal eta,phi: " << eta << " " << phi << std::endl;
    } // EE
    nEExtals_filled = vEExtals_in_HBHEtower.size();
    // Loop over selected EE xtals
    for ( int iC = 0; iC < nEExtals_filled; iC++ ) {
      // Split HCAL tower energy evenly among xtals
      EEDetId eeId( EEDetId::unhashIndex(vEExtals_in_HBHEtower[iC]) );
      ix_ = eeId.ix() - 1;
      iy_ = eeId.iy() - 1;
      iz_ = (eeId.zside() > 0) ? 1 : 0;
      // Create hashed Index: maps from [iy][ix] -> [idx_]
      idx_ = iy_*EE_MAX_IX + ix_;
      // Fill vector for images
      vHBHE_energy_EE_[iz_][idx_] += ( energy_/float(nEExtals_filled) );
      // Fill histogram for monitoring
      hHBHE_energy_EE_[iz_]->Fill( ix_, iy_, energy_ );
    } // EE, selected
    */

  } // EE rechits

  /*
  for ( int ieta = 19; ieta < 24; ieta++) {
    for ( int iphi = 1; iphi < 7; iphi++) {
      //HcalDetId hId(HcalEndcap, ieta, iphi, 1); //ieta,iphi,d
      HcalDetId hId( iRHit->id() );
      //if ( hId.ieta() < 20 || hId.ieta() > 22 ) continue;

      pos = caloGeom->getGeometry(hId)->getPosition();
      const CaloCellGeometry::RepCorners& repCorners = caloGeom->getGeometry(hId)->getCornersREP();
      for ( unsigned c = 0; c < repCorners.size(); c++ ) {
        std::cout << c << ":rho,eta,phi: " << repCorners[c].rho() << " " << repCorners[c].eta() << " " << repCorners[c].phi() << std::endl;
      }
      //pos = caloGeom->getPosition( hId );
      eta = pos.eta();
      phi = pos.phi();
      //std::cout << "ieta"<<ieta<<hId.ieta()<<",iphi"<<iphi<<hId.iphi()<<": " << eta << " , " << phi << std::endl;
      std::cout << "ieta"<<hId.ieta()<<",iphi"<<hId.iphi()<<": " << eta << " , " << phi << std::endl;
    }
    std::cout << "" <<std::endl;
  }

  // Fill HBHE rechits
  for ( HBHERecHitCollection::const_iterator iRHit = HBHERecHitsH_->begin();
      iRHit != HBHERecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ <= zs ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );

    if ( id.subdetId() == EcalBarrel ) continue;

    // Get detector id and convert to histogram-friendly coordinates
    EEDetId eeId( id );
    ix_ = eeId.ix() - 1;
    iy_ = eeId.iy() - 1;
    iz_ = (eeId.zside() > 0) ? 1 : 0;
    // Fill histograms for monitoring
    hHBHE_energy_EE_[iz_]->Fill( ix_,iy_,energy_ );
    // Create hashed Index: maps from [iy][ix] -> [idx_]
    idx_ = iy_*EE_MAX_IX + ix_;
    // Fill vectors for images
    vHBHE_energy_EE_[iz_][idx_] += energy_;

  //if ( hId.eta() > 20 ) fillEEhit( HcalDetId hId(ieta_,iphi_+1,0) )

  } // HBHE rechits
  */

} // fillEEatES()

// HCAL REPcorners index illustration:
/*
         6_____7   eta
         /|   /|    |
        / |_ /_|   \|/
       / 5  / / 4
      /    / /
     /    / /
    /    / /
  2/___3/ /
   |   | /
   |_ _|/   phi__\
  1     0        /

   / IP
 |/_

*/
