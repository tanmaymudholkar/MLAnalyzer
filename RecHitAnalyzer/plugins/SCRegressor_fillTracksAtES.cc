#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Fill EE rec hits /////////////////////////////////////////
// For each endcap, store event rechits in a vector of length
// equal to number of crystals per endcap (ix:100 x iy:100)

TH2F *hTracksAtES_energy[nZ];
std::vector<float> vTracksAtES_energy[nZ];
const int TrksAtES_search_window = 3*nSTRIP;

// Initialize branches _____________________________________________________________//
void SCRegressor::branchesTracksAtES ( TTree* tree, edm::Service<TFileService> &fs ) {

  char hname[50], htitle[50];
  for ( int iz(0); iz < nEE; iz++ ) {
    // Branches for images
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "TracksPt_ES%s",zside);
    tree->Branch(hname,        &vTracksAtES_energy[iz]);

    // Histograms for monitoring
    sprintf(htitle,"pT(stripX,stripY);stripX;stripY");
    hTracksAtES_energy[iz] = fs->make<TH2F>(hname, htitle,
        nXY_STRIP, 0, nXY_STRIP,
        nXY_STRIP, 0, nXY_STRIP  );
  } // iz

} // branchesTracksAtES()

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

// Fill HCAL rechits at EB/EE ______________________________________________________________//
void SCRegressor::fillTracksAtES ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int ix_, iy_, iz_, idx_; // rows:iy, columns:ix
  float energy_;
  float tkEta, tkPhi, tkPt, minY_, minX_, maxY_, maxX_;
  float st_x, st_y, tk_x, tk_y;
  float st_minx, st_miny;
  float st_maxx, st_maxy;
  int stX_, stY_;
  int stripX_seed, stripY_seed, st_ix, st_iy, st_xstrip, st_ystrip;
  GlobalPoint pos;
  float min_dX, min_dY, min_stX_, min_stY_;

  for ( int iz(0); iz < nEE; iz++ ) {
    vTracksAtES_energy[iz].assign( nXY_STRIP*nXY_STRIP, 0. );
  }

  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom;
  caloGeom = caloGeomH_.product();

  edm::Handle<reco::TrackCollection> tracksH_;
  //edm::Handle<pat::IsolatedTrackCollection> tracksH_;
  iEvent.getByToken( trackCollectionT_, tracksH_ );

  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
  //for ( pat::IsolatedTrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) {
    if ( !(iTk->quality(tkQt_)) ) continue;

    tkEta = iTk->eta();
    tkPhi = iTk->phi();
    tkPt  = iTk->pt();
    tk_x = iTk->vx();
    tk_y = iTk->vy();
    if ( std::abs(tkEta) < 1.7 || std::abs(tkEta) > 2.5 ) continue;

    // Get closest ES id
    ESDetId nearestESId = ESId_from_EtaPhi( tkEta, tkPhi, caloGeom );
    ix_ = nearestESId.six() - 1;
    iy_ = nearestESId.siy() - 1;
    iz_ = (nearestESId.zside() > 0) ? 1 : 0;
    stripX_seed = nSTRIP*ix_ + (nSTRIP/2);
    stripY_seed = nSTRIP*iy_ + (nSTRIP/2);
    //std::cout << ESDetId::validDetId(nearestESId.strip(), nearestESId.six(), nearestESId.siy(), nearestESId.plane(), nearestESId.zside()) << std::endl;
    // Search window around closest ES id:
    min_dX = 100.;
    min_dY = 100.;
    min_stX_ = -1;
    min_stY_ = -1;
    for ( int stX = 0; stX < TrksAtES_search_window; stX++ ) {

      stX_ = stripX_seed - (TrksAtES_search_window/2)+stX;
      if ( stX_ > nXY_STRIP-1 ) continue;
      if ( stX_ < 0 ) continue;

      for ( int stY = 0; stY < TrksAtES_search_window; stY++ ) {

        stY_ = stripY_seed - (TrksAtES_search_window/2)+stY;
        if ( stY_ > nXY_STRIP-1 ) continue;
        if ( stY_ < 0 ) continue;

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

        if ( !ESDetId::validDetId(st_ystrip+1, st_ix+1, st_iy+1, 2, iz_ > 0 ? 1 : -1) ) continue;
        ESDetId stId_y  = ESDetId(st_ystrip+1, st_ix+1, st_iy+1, 2, iz_ > 0 ? 1 : -1);
        st_y = caloGeom->getPosition( stId_y ).y(); // ES_R

        if ( !(std::abs(st_x - tk_x) < min_dX) ) continue;
        if ( !(std::abs(st_y - tk_y) < min_dY) ) continue;
        //std::cout << "ix,xstrip, iy,ystrip:" << stId_x.six() << "," << stId_x.strip() << " , "<< stId_y.siy() << "," << stId_y.strip() << std::endl;

        min_dX = std::abs(st_x - tk_x);
        min_dY = std::abs(st_y - tk_y);
        min_stX_ = stX_;
        min_stY_ = stY_;
        //std::cout << "st_x,st_y:" << st_x << "," << st_y << std::endl;

      } // siY
    } // siX
    // loop over zip(stX_,stY_)
    if ( min_stX_ == -1 || min_stY_ == -1 ) continue;
    idx_ = min_stY_*nXY_STRIP + min_stX_;
    vTracksAtES_energy[iz_][idx_] += tkPt;
    hTracksAtES_energy[iz_]->SetBinContent( min_stX_, min_stY_, hTracksAtES_energy[iz_]->GetBinContent(min_stX_, min_stY_) + tkPt );

  } // tracks

} // fillTracksAtES()
