#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

// Fill stitched EEm_EB_EEp image /////////////////////............/
// Store all ECAL event rechits into a stitched EEm_EB_EEp image 
// segmented in iphi,ieta spannning the full -3 < eta < 3. 
// Use EB-like granularity giving an extended range ieta=[-140,140].
//
// For endcaps, project EE hits into a helper histogram binned by 
// phi,eta before filling the full extended ECAL(iphi,eta) image.
// For barrel, fill EB hits directly since geometries are 1:1. 
//
// 'ieta_global' keeps track of the global ieta index count used
// for filling the extended image vector vECAL_tracksPt.
// 'ieta_signed' keeps track of the position along [-140,140] used
// for filling the monitoring histogram hECAL_tracksPt.

TH2F *hEvt_EE_tracks[nEE];
TH2F *hEvt_EE_tracksPt[nEE];
TH2F *hEvt_EE_tracksD0[nEE];
TH2F *hEvt_EE_tracksDz[nEE];
TH2F *hEvt_EE_tracksPt_max[nEE];
TH2F *hEvt_EE_tracksD0_max[nEE];
TH2F *hEvt_EE_tracksDz_max[nEE];
TProfile2D *hECAL_tracks;
TProfile2D *hECAL_tracksPt;
TProfile2D *hECAL_tracksD0;
TProfile2D *hECAL_tracksDz;
std::vector<float> vECAL_tracksPt_;
std::vector<float> vECAL_tracksD0_;
std::vector<float> vECAL_tracksDz_;
std::vector<float> vECAL_tracks_;

std::vector<float> vECAL_tracksPt_max_;
std::vector<float> vECAL_tracksD0_max_;
std::vector<float> vECAL_tracksDz_max_;

// TODO Take the D0 and Dz values from the highest pT track for endcap regions

// Initialize branches _______________________________________________________________//
void RecHitAnalyzer::branchesTracksAtECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("ECAL_tracks",      &vECAL_tracks_);
  tree->Branch("ECAL_tracksPt",    &vECAL_tracksPt_);
  tree->Branch("ECAL_tracksD0",    &vECAL_tracksD0_);
  tree->Branch("ECAL_tracksDz",    &vECAL_tracksDz_);
  tree->Branch("ECAL_tracksPt_maxPt",    &vECAL_tracksPt_max_);
  tree->Branch("ECAL_tracksD0_maxPt",    &vECAL_tracksD0_max_);
  tree->Branch("ECAL_tracksDz_maxPt",    &vECAL_tracksDz_max_);
  // Intermediate helper histogram (single event only)
  // Track Intermediate Hists
  hEvt_EE_tracks[0] = new TH2F("evt_EEm_tracks", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracks[1] = new TH2F("evt_EEp_tracks", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  // Track Parameter Intermediate Hists
  hEvt_EE_tracksPt[0] = new TH2F("evt_EEm_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksPt[1] = new TH2F("evt_EEp_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksD0[0] = new TH2F("evt_EEm_tracksD0", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksD0[1] = new TH2F("evt_EEp_tracksD0", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksDz[0] = new TH2F("evt_EEm_tracksDz", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksDz[1] = new TH2F("evt_EEp_tracksDz", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  // Track Max Parameter Intermediate Hists
  hEvt_EE_tracksPt_max[0] = new TH2F("evt_EEm_tracksPt_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksPt_max[1] = new TH2F("evt_EEp_tracksPt_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksD0_max[0] = new TH2F("evt_EEm_tracksD0_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksD0_max[1] = new TH2F("evt_EEp_tracksD0_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksDz_max[0] = new TH2F("evt_EEm_tracksDz_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksDz_max[1] = new TH2F("evt_EEp_tracksDz_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  // Histograms for monitoring
  hECAL_tracks = fs->make<TProfile2D>("ECAL_tracks", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
  hECAL_tracksPt = fs->make<TProfile2D>("ECAL_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
  hECAL_tracksD0 = fs->make<TProfile2D>("ECAL_tracksD0", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
  hECAL_tracksDz = fs->make<TProfile2D>("ECAL_tracksDz", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

} // branchesTracksAtECALstitched()

// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
void fillTracksAtECAL_with_EEproj ( TH2F *hEvt_EE_tracks_, TH2F *hEvt_EE_tracksPt_, TH2F *hEvt_EE_tracksD0_, TH2F *hEvt_EE_tracksDz_,  TH2F *hEvt_EE_tracksD0_max_, TH2F *hEvt_EE_tracksDz_max_, int ieta_global_offset, int ieta_signed_offset ) {

  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float track_;
  float trackPt_;
  float trackD0_; 
  float trackDz_;
  float trackD0_max_;
  float trackDz_max_;
  for (int ieta = 1; ieta < hEvt_EE_tracksPt_->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_tracksPt_->GetNbinsX()+1; iphi++) {

      track_  = hEvt_EE_tracks_->GetBinContent( iphi, ieta );
      trackPt_ = hEvt_EE_tracksPt_->GetBinContent( iphi, ieta );
      trackD0_ = hEvt_EE_tracksD0_->GetBinContent( iphi, ieta );
      trackDz_ = hEvt_EE_tracksDz_->GetBinContent( iphi, ieta );
      trackD0_max_ = hEvt_EE_tracksD0_max_->GetBinContent( iphi, ieta );
      trackDz_max_ = hEvt_EE_tracksDz_max_->GetBinContent( iphi, ieta );
      if ( trackPt_ <= zs ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
      // Fill vector for image
      vECAL_tracks_[idx_] = track_;
      vECAL_tracksPt_[idx_] = trackPt_;
      vECAL_tracksD0_[idx_] = trackD0_;
      vECAL_tracksDz_[idx_] = trackDz_;
      vECAL_tracksD0_max_[idx_] = trackD0_max_;
      vECAL_tracksDz_max_[idx_] = trackDz_max_;
      // Fill histogram for monitoring
      hECAL_tracks->Fill( iphi_, ieta_signed_, track_ );
      hECAL_tracksPt->Fill( iphi_, ieta_signed_, trackPt_ );
      hECAL_tracksD0->Fill( iphi_, ieta_signed_, trackD0_ );
      hECAL_tracksDz->Fill( iphi_, ieta_signed_, trackDz_ );

    } // iphi_
  } // ieta_

} // fillTracksAtECAL_with_EEproj

// Fill stitched EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillTracksAtECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, trackPt_, trackD0_, trackDz_;
  GlobalPoint pos;

  vECAL_tracks_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksD0_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksDz_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksPt_max_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksD0_max_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksDz_max_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  for ( int iz(0); iz < nEE; ++iz ) {
    hEvt_EE_tracks[iz]->Reset();
    hEvt_EE_tracksPt[iz]->Reset();
    hEvt_EE_tracksD0[iz]->Reset();
    hEvt_EE_tracksDz[iz]->Reset();
    hEvt_EE_tracksPt_max[iz]->Reset();
    hEvt_EE_tracksD0_max[iz]->Reset();
    hEvt_EE_tracksDz_max[iz]->Reset();
  }

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByLabel( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByLabel( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  edm::Handle<reco::TrackCollection> tracksH_;
  iEvent.getByLabel( trackCollectionT_, tracksH_ );

  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  int bin;
  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) {
    if ( !(iTk->quality(tkQt_)) ) continue;
    eta = iTk->eta();
    phi = iTk->phi();
    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) continue;
    if ( id.subdetId() == EcalEndcap ) {
      iz_ = (eta > 0.) ? 1 : 0;
      // Fill intermediate helper histogram by eta,phi
      hEvt_EE_tracks[iz_]->Fill( phi, eta );
      hEvt_EE_tracksPt[iz_]->Fill( phi, eta, iTk->pt() );
      hEvt_EE_tracksD0[iz_]->Fill( phi, eta, iTk->d0() );
      hEvt_EE_tracksDz[iz_]->Fill( phi, eta, iTk->dz() );
      bin = hEvt_EE_tracks[iz_]->FindBin( phi, eta );
      if ( iTk->pt() > hEvt_EE_tracksPt_max[iz_]->GetBinContent( bin ) ) {
        hEvt_EE_tracksPt_max[iz_]->SetBinContent( bin, iTk->pt() );
        hEvt_EE_tracksD0_max[iz_]->SetBinContent( bin, iTk->d0() );
        hEvt_EE_tracksDz_max[iz_]->SetBinContent( bin, iTk->dz() );
      }
    }
  } // tracks
 
  // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  ieta_global_offset = 0;
  ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  fillTracksAtECAL_with_EEproj( hEvt_EE_tracks[0], hEvt_EE_tracksPt[0], hEvt_EE_tracksD0[0], hEvt_EE_tracksDz[0], hEvt_EE_tracksD0_max[0], hEvt_EE_tracksDz_max[0], ieta_global_offset, ieta_signed_offset );

  // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  ieta_global_offset = 55;

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) { 
    if ( !(iTk->quality(tkQt_)) ) continue;
    eta = iTk->eta();
    phi = iTk->phi();
    trackPt_ = iTk->pt();
    trackD0_ = iTk->d0();
    trackDz_ = iTk->dz();
    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalEndcap ) continue;
    if ( id.subdetId() == EcalBarrel ) { 
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      if ( trackPt_ <= zs ) continue;
      // Fill vector for image
      ieta_signed = ieta_;
      ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
      idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
      vECAL_tracks_[idx_] += 1.;
      vECAL_tracksPt_[idx_] += trackPt_;
      vECAL_tracksD0_[idx_] += trackD0_;
      vECAL_tracksDz_[idx_] += trackDz_;
      if ( trackPt_ > vECAL_tracksPt_max_[idx_] ) {
        vECAL_tracksPt_max_[idx_] = trackPt_;
        vECAL_tracksD0_max_[idx_] = trackD0_;
        vECAL_tracksDz_max_[idx_] = trackDz_;
      }
      // Fill histogram for monitoring
      hECAL_tracks->Fill( iphi_, ieta_signed );
      hECAL_tracksPt->Fill( iphi_, ieta_signed, trackPt_ );
      hECAL_tracksD0->Fill( iphi_, ieta_signed, trackD0_ );
      hECAL_tracksDz->Fill( iphi_, ieta_signed, trackDz_ );
    }

  } // EB
  // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
  ieta_signed_offset = EB_IETA_MAX;
  fillTracksAtECAL_with_EEproj( hEvt_EE_tracks[1], hEvt_EE_tracksPt[1], hEvt_EE_tracksD0[1], hEvt_EE_tracksDz[1], hEvt_EE_tracksD0_max[1], hEvt_EE_tracksDz_max[1], ieta_global_offset, ieta_signed_offset );

  // Get average D0 and Dz for each position
  for (unsigned int idx_=0;idx_<vECAL_tracks_.size();idx_++) {
    if (vECAL_tracks_[idx_] != 0) {
      vECAL_tracksD0_[idx_] = vECAL_tracksD0_[idx_] / vECAL_tracks_[idx_];
      vECAL_tracksDz_[idx_] = vECAL_tracksDz_[idx_] / vECAL_tracks_[idx_];
    }
  }

} // fillTracksAtECALstitched()
