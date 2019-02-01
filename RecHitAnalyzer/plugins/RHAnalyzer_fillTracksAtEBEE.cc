#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"
/*
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
*/
// Fill EB rec hits ////////////////////////////////
// Store event rechits in a vector of length equal
// to number of crystals in EB (ieta:170 x iphi:360)

TH2F *hTracks_EE[nEE];
TH2F *hTracks_EB;
TH2F *hTracksPt_EE[nEE];
TH2F *hTracksPt_EB;
TH2F *hTracksD0_EE[nEE];
TH2F *hTracksD0_EB;
std::vector<float> vTracksPt_EE_[nEE];
std::vector<float> vTracksD0_EE_[nEE];
std::vector<float> vTracks_EE_[nEE];
std::vector<float> vTracksPt_EB_;
std::vector<float> vTracksD0_EB_;
std::vector<float> vTracks_EB_;

// Initialize branches ____________________________________________________________//
void RecHitAnalyzer::branchesTracksAtEBEE ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("Tracks_EB",   &vTracks_EB_);
  tree->Branch("TracksPt_EB", &vTracksPt_EB_);
  tree->Branch("TracksD0_EB", &vTracksD0_EB_);

  // Histograms for monitoring
  hTracks_EB = fs->make<TH2F>("Tracks_EB", "N(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hTracksPt_EB = fs->make<TH2F>("TracksPt_EB", "pT(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hTracksD0_EB = fs->make<TH2F>("TracksD0_EB", "d0(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );

  char hname[50], htitle[50];
  for ( int iz(0); iz < nEE; iz++ ) {
    // Branches for images
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "Tracks_EE%s",zside);
    tree->Branch(hname,        &vTracks_EE_[iz]);
    sprintf(hname, "TracksPt_EE%s",zside);
    tree->Branch(hname,        &vTracksPt_EE_[iz]);
    sprintf(hname, "TracksD0_EE%s",zside);
    tree->Branch(hname,       &vTracksD0_EE_[iz]);

    // Histograms for monitoring
    sprintf(hname, "Tracks_EE%s",zside);
    sprintf(htitle,"N(ix,iy);ix;iy");
    hTracks_EE[iz] = fs->make<TH2F>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
    sprintf(hname, "TracksPt_EE%s",zside);
    sprintf(htitle,"pT(ix,iy);ix;iy");
    hTracksPt_EE[iz] = fs->make<TH2F>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
    sprintf(hname, "TracksD0_EE%s",zside);
    sprintf(htitle,"D0(ix,ix);ix;iy");
    hTracksD0_EE[iz] = fs->make<TH2F>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
  } // iz

} // branchesEB()

// Fill TRK rechits at EB/EE ______________________________________________________________//
void RecHitAnalyzer::fillTracksAtEBEE ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int ix_, iy_, iz_;
  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  float eta, phi;
  GlobalPoint pos;

  vTracks_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vTracksPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vTracksD0_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  for ( int iz(0); iz < nEE; iz++ ) {
    vTracks_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vTracksPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vTracksD0_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
  }

  edm::Handle<reco::TrackCollection> tracksH_;
  iEvent.getByLabel( trackCollectionT_, tracksH_ );
  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) {
    if ( !(iTk->quality(tkQt_)) ) continue;
    eta = iTk->eta();
    phi = iTk->phi();
    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) {
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      // Fill histograms for monitoring
      hTracks_EB->Fill( iphi_, ieta_ );
      hTracksPt_EB->Fill( iphi_, ieta_, iTk->pt() );
      hTracksD0_EB->Fill( iphi_, ieta_, iTk->d0() );
      idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
      // Fill vectors for images
      vTracks_EB_[idx_] += 1.;
      vTracksPt_EB_[idx_] += iTk->pt();
      vTracksD0_EB_[idx_] += iTk->d0();
    } else if ( id.subdetId() == EcalEndcap ) {
      EEDetId eeId( id );
      ix_ = eeId.ix() - 1;
      iy_ = eeId.iy() - 1;
      iz_ = (eeId.zside() > 0) ? 1 : 0;
      // Fill histograms for monitoring
      hTracks_EE[iz_]->Fill( ix_, iy_ );
      hTracksPt_EE[iz_]->Fill( ix_, iy_, iTk->pt() );
      hTracksD0_EE[iz_]->Fill( ix_, iy_, iTk->d0() );
      // Create hashed Index: maps from [iy][ix] -> [idx_]
      idx_ = iy_*EE_MAX_IX + ix_;
      // Fill vectors for images
      vTracks_EE_[iz_][idx_] += 1.;
      vTracksPt_EE_[iz_][idx_] += iTk->pt();
      vTracksD0_EE_[iz_][idx_] += iTk->d0();
    } 
  } // tracks

  //Get average D0 for each position
  for (unsigned int idx__=0;idx__<vTracks_EB_.size();idx__++) {
    if (vTracks_EB_[idx__] != 0) { vTracksD0_EB_[idx__] = vTracksD0_EB_[idx__] / vTracks_EB_[idx__]; }
  }
  for (int iz__=0;iz__<nEE;iz__++) {
    for (unsigned int idx__=0;idx__<vTracks_EE_[iz__].size();idx__++) {
      if (vTracks_EE_[iz__][idx__] != 0) { vTracksD0_EE_[iz__][idx__] = vTracksD0_EE_[iz__][idx__] / vTracks_EE_[iz__][idx__]; }
    }
  }


} // fillEB()
