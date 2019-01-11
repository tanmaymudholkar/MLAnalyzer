#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

// Fill Tracks into stitched EEm_EB_EEp image //////////////////////
// Store all Track positions into a stitched EEm_EB_EEp image 

TH2F *hEvt_EE_EndtracksPt[nEE];
TH2F *hEvt_EE_muonsPt[nEE];
TProfile2D *hECAL_EndtracksPt;
TProfile2D *hECAL_muonsPt;
std::vector<float> vECAL_EndtracksPt_;
std::vector<float> vECAL_muonsPt_;

// Initialize branches _______________________________________________________________//
void RecHitAnalyzer::branchesPFCandsAtECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("ECAL_EndtracksPt",    &vECAL_EndtracksPt_);
  tree->Branch("ECAL_muonsPt",    &vECAL_muonsPt_);

  // Intermediate helper histogram (single event only)
  hEvt_EE_EndtracksPt[0] = new TH2F("evt_EEm_EndtracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_EndtracksPt[1] = new TH2F("evt_EEp_EndtracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  hEvt_EE_muonsPt[0] = new TH2F("evt_EEm_muonsPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_muonsPt[1] = new TH2F("evt_EEp_muonsPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );


  // Histograms for monitoring
  hECAL_EndtracksPt = fs->make<TProfile2D>("ECAL_EndtracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

  hECAL_muonsPt = fs->make<TProfile2D>("ECAL_muonsPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

} // branchesTracksAtECALstitched()

// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
void fillPFCandsAtECAL_with_EEproj ( TH2F *hEvt_EE_EndtracksPt_, TH2F *hEvt_EE_muonsPt_, int ieta_global_offset, int ieta_signed_offset ) {

  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float EndtrackPt_;
  float muonPt_;

  for (int ieta = 1; ieta < hEvt_EE_EndtracksPt_->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_EndtracksPt_->GetNbinsX()+1; iphi++) {

      EndtrackPt_ = hEvt_EE_EndtracksPt_->GetBinContent( iphi, ieta );
      muonPt_ = hEvt_EE_muonsPt_->GetBinContent( iphi, ieta );
      if ( (EndtrackPt_ <= zs)  && (muonPt_ <= zs)) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
      // Fill vector for image
      vECAL_EndtracksPt_[idx_] = EndtrackPt_;
      vECAL_muonsPt_[idx_] = muonPt_;
      // Fill histogram for monitoring
      hECAL_EndtracksPt->Fill( iphi_, ieta_signed_, EndtrackPt_ );
      hECAL_muonsPt->Fill( iphi_, ieta_signed_, muonPt_ );

    } // iphi_
  } // ieta_

} // fillPFCandsAtECAL_with_EEproj

// Fill stitched EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillPFCandsAtECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, trackPt_;
  GlobalPoint pos;

  vECAL_EndtracksPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_muonsPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_EndtracksPt[iz]->Reset();
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_muonsPt[iz]->Reset();

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByLabel( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByLabel( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  edm::Handle<reco::PFCandidateCollection> pfCandsH_;
  iEvent.getByLabel(pfCandCollectionT_, pfCandsH_);

  for ( reco::PFCandidateCollection::const_iterator iPFC = pfCandsH_->begin();
        iPFC != pfCandsH_->end(); ++iPFC ) {
    const reco::Track* thisTrk = &(*iPFC->trackRef());
    if(!thisTrk) continue;

    const math::XYZPointF& ecalPos = iPFC->positionAtECALEntrance();
    eta = ecalPos.eta();
    phi = ecalPos.phi();

    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) continue;
    if ( id.subdetId() == EcalEndcap ) {
      iz_ = (eta > 0.) ? 1 : 0;
      // Fill intermediate helper histogram by eta,phi
      hEvt_EE_EndtracksPt[iz_]->Fill( phi, eta, thisTrk->pt() );
      if(iPFC->particleId() == 3)
	hEvt_EE_muonsPt[iz_]->Fill( phi, eta, thisTrk->pt() );
    }
  } // pfCands



  // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  ieta_global_offset = 0;
  ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  fillPFCandsAtECAL_with_EEproj( hEvt_EE_EndtracksPt[0], hEvt_EE_muonsPt[0], ieta_global_offset, ieta_signed_offset );

  // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  ieta_global_offset = 55;


  for ( reco::PFCandidateCollection::const_iterator iPFC = pfCandsH_->begin();
        iPFC != pfCandsH_->end(); ++iPFC ) {
    const reco::Track* thisTrk = &(*iPFC->trackRef());
    if(!thisTrk) continue;

    const math::XYZPointF& ecalPos = iPFC->positionAtECALEntrance();
    eta = ecalPos.eta();
    phi = ecalPos.phi();

    trackPt_ = thisTrk->pt();

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
      vECAL_EndtracksPt_[idx_] += trackPt_;

      // Fill histogram for monitoring
      hECAL_EndtracksPt->Fill( iphi_, ieta_signed, trackPt_ );
      if(iPFC->particleId() == 3){
	vECAL_muonsPt_[idx_] += trackPt_;
	hECAL_muonsPt->Fill( iphi_, ieta_signed, trackPt_ );
      }
    }

  } // EB PFCands


  // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
  ieta_signed_offset = EB_IETA_MAX;
  fillPFCandsAtECAL_with_EEproj( hEvt_EE_EndtracksPt[1], hEvt_EE_muonsPt[1], ieta_global_offset, ieta_signed_offset );

} // fillPFCandsAtECALstitched()
