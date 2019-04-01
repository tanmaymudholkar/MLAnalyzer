#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
//#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

// Fill adjustable EEm_EB_EEp image /////////////////////............/
// Store all ECAL event rechits into a adjustable EEm_EB_EEp image 
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

TH2F *hEvt_Adj_tracks[Nadjproj];
TH2F *hEvt_Adj_tracksPt[Nadjproj];
TH2F *hEvt_Adj_tracksQPt[Nadjproj];
TH2F *hEvt_Adj_tracksD0[Nadjproj];
TH2F *hEvt_Adj_tracksDz[Nadjproj];
TH2F *hEvt_Adj_tracksPt_max[Nadjproj];
TH2F *hEvt_Adj_tracksD0_max[Nadjproj];
TH2F *hEvt_Adj_tracksDz_max[Nadjproj];

TH2F *hECALadj_tracks[Nadjproj];
TH2F *hECALadj_tracksPt[Nadjproj];
TH2F *hECALadj_tracksQPt[Nadjproj];
TH2F *hECALadj_tracksD0[Nadjproj];
TH2F *hECALadj_tracksDz[Nadjproj];

std::vector<float> vECALadj_tracksPt_[Nadjproj];
std::vector<float> vECALadj_tracksQPt_[Nadjproj];
std::vector<float> vECALadj_tracksD0_[Nadjproj];
std::vector<float> vECALadj_tracksDz_[Nadjproj];
std::vector<float> vECALadj_tracks_[Nadjproj];
std::vector<float> vECALadj_tracksPt_max_[Nadjproj];
std::vector<float> vECALadj_tracksD0_max_[Nadjproj];
std::vector<float> vECALadj_tracksDz_max_[Nadjproj];


// TODO Take the D0 and Dz values from the highest pT track for endcap regions

// Initialize branches _______________________________________________________________//
void RecHitAnalyzer::branchesTracksAtECALadjustable ( TTree* tree, edm::Service<TFileService> &fs ) {
  
  for (unsigned int proj=0; proj<Nadjproj; proj++)
  {
    // Branches for images
    tree->Branch((std::string("ECALadj_tracks")+adj_projections[proj]).c_str(),      &(vECALadj_tracks_[proj]));
    tree->Branch((std::string("ECALadj_tracksPt")+adj_projections[proj]).c_str(),    &(vECALadj_tracksPt_[proj]));
    tree->Branch((std::string("ECALadj_tracksQPt")+adj_projections[proj]).c_str(),    &(vECALadj_tracksPt_[proj]));
    tree->Branch((std::string("ECALadj_tracksD0")+adj_projections[proj]).c_str(),    &(vECALadj_tracksD0_[proj]));
    tree->Branch((std::string("ECALadj_tracksDz")+adj_projections[proj]).c_str(),    &(vECALadj_tracksDz_[proj]));
    tree->Branch((std::string("ECALadj_tracksPt_maxPt")+adj_projections[proj]).c_str(),    &(vECALadj_tracksPt_max_[proj]));
    tree->Branch((std::string("ECALadj_tracksD0_maxPt")+adj_projections[proj]).c_str(),    &(vECALadj_tracksD0_max_[proj]));
    tree->Branch((std::string("ECALadj_tracksDz_maxPt")+adj_projections[proj]).c_str(),    &(vECALadj_tracksDz_max_[proj]));

    std::vector<double> adjEtaBins;
    std::vector<double> adjPhiBins;

    int eta_nbins_HBHE = 2*(hcaldqm::constants::IETA_MAX_HE-1);
    int totalMultiEta[proj] = granularityMultiEta[proj] * granularityMultiECAL[proj];

    for (int i=0; i<eta_nbins_HBHE; i++)
    {
      double step=(eta_bins_HBHE[i+1]-eta_bins_HBHE[i])/totalMultiEta[proj];
      for (int j=0; j<totalMultiEta; j++)
      {
        adjEtaBins.push_back(eta_bins_HBHE[i]+step*j);
      }
    }
    adjEtaBins.push_back(eta_bins_HBHE[eta_nbins_HBHE]);

    totalEtaBins = totalMultiEta*(eta_nbins_HBHE);
    totalPhiBins = granularityMultiPhi * granularityMultiECAL*HBHE_IPHI_NUM;


    // Histograms for monitoring
    hECALadj_tracks[proj] = fs->make<TH2F>((std::string("ECALadj_tracks")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksPt[proj] = fs->make<TH2F>((std::string("ECALadj_tracksPt")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksQPt[proj] = fs->make<TH2F>((std::string("ECALadj_tracksQPt")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksD0[proj] = fs->make<TH2F>((std::string("ECALadj_tracksD0")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksDz[proj] = fs->make<TH2F>((std::string("ECALadj_tracksDz")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );


    hEvt_Adj_tracks[proj] = new TH2F((std::string("evt_Adj_tracks")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hEvt_Adj_tracksPt[proj] = new TH2F((std::string("evt_Adj_tracksPt")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hEvt_Adj_tracksQPt[proj] = new TH2F((std::string("evt_Adj_tracksQPt")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hEvt_Adj_tracksD0[proj] = new TH2F((std::string("evt_Adj_tracksD0")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hEvt_Adj_tracksDz[proj] = new TH2F((std::string("evt_Adj_tracksDz")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hEvt_Adj_tracksPt_max[proj] = new TH2F((std::string("evt_Adj_tracksPt_max")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hEvt_Adj_tracksD0_max[proj] = new TH2F((std::string("evt_Adj_tracksD0_max")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hEvt_Adj_tracksDz_max[proj] = new TH2F((std::string("evt_Adj_tracksDz_max")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );


  }

} // branchesTracksAtECALadjustable()

std::vector<int> RecHitAnalyzer::findSubcrystal(const CaloGeometry* caloGeom, const float& eta, const float& phi, const int& granularityMultiEta, const int& granularityMultiPhi)
{

    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel )
    { 
      auto subDetGeometry = caloGeom->getSubdetectorGeometry(id);
      auto caloCellGeometry = subDetGeometry->getGeometry(id);
      auto corners = caloCellGeometry->getCornersREP();
      //float kappa= 4*TMath::Pi()/HBHE_IPHI_NUM;
      // if (phi>-kappa)
      //   phi=phi+kappa-TMath::Pi();
      // else  
      //   phi=phi+kappa+TMath::Pi();

      std::vector<float> eta_corners = {corners[0].eta(),corners[1].eta(),corners[2].eta(),corners[3].eta()};
      std::vector<float> phi_corners = {corners[0].phi(),corners[1].phi(),corners[2].phi(),corners[3].phi()};

      auto lowEta_lowPhi_index = 4;
      auto highEta_lowPhi_index = 4;
      auto lowEta_highPhi_index = 4;
      auto highEta_highPhi_index = 4;

      std::vector<size_t> eta_sorted_indices = {0,1,2,3};//(eta_corners.size());
      std::vector<size_t> phi_sorted_indices = {0,1,2,3};//(phi_corners.size());
      //std::iota(eta_sorted_indices.begin(), eta_sorted_indices.end(), 0);
      //std::iota(phi_sorted_indices.begin(), phi_sorted_indices.end(), 0);

      //index sort with lambdas
      std::sort(eta_sorted_indices.begin(), eta_sorted_indices.end(),[&eta_corners](size_t i1, size_t i2) {return eta_corners[i1] < eta_corners[i2];});
      std::sort(phi_sorted_indices.begin(), phi_sorted_indices.end(),[&phi_corners](size_t i1, size_t i2) {return phi_corners[i1] < phi_corners[i2];});

      for (unsigned int i =0; i<4; i++)
      {
        if      ((i==eta_sorted_indices[0] || i==eta_sorted_indices[1]) &&
                 (i==phi_sorted_indices[0] || i==phi_sorted_indices[1]) )
        { lowEta_lowPhi_index = i; }
        else if ((i==eta_sorted_indices[2] || i==eta_sorted_indices[3]) &&
                 (i==phi_sorted_indices[2] || i==phi_sorted_indices[3]) )
        { highEta_highPhi_index = i;  }
        else if ((i==eta_sorted_indices[0] || i==eta_sorted_indices[1]) &&
                 (i==phi_sorted_indices[2] || i==phi_sorted_indices[3]) )
        { lowEta_highPhi_index = i;  }
        else if ((i==eta_sorted_indices[2] || i==eta_sorted_indices[3]) &&
                 (i==phi_sorted_indices[0] || i==phi_sorted_indices[1]) )
        { highEta_lowPhi_index = i;  }
      }

      // if (lowEta_lowPhi_index  ==4 ||
      //     highEta_lowPhi_index ==4 ||
      //     lowEta_highPhi_index ==4 ||
      //     highEta_highPhi_index==4 ) std::cout<<"something went wrong\n";

      TVector2 lowEta_lowPhi_corner(corners[lowEta_lowPhi_index].eta(),corners[lowEta_lowPhi_index].phi());
      TVector2 highEta_lowPhi_corner(corners[highEta_lowPhi_index].eta(),corners[highEta_lowPhi_index].phi());
      TVector2 lowEta_highPhi_corner(corners[lowEta_highPhi_index].eta(),corners[lowEta_highPhi_index].phi());
      TVector2 highEta_highPhi_corner(corners[highEta_highPhi_index].eta(),corners[highEta_highPhi_index].phi());

      float subcrystal_eta_edges[granularityMultiEta+1][granularityMultiPhi+1];
      float subcrystal_phi_edges[granularityMultiEta+1][granularityMultiPhi+1];

      for (unsigned int etaIndex=0; etaIndex<granularityMultiEta+1; etaIndex++)
      {
        for (unsigned int phiIndex=0; phiIndex<granularityMultiPhi+1; phiIndex++)
        {
           TVector2 aveEta_lowPhi_corner (((granularityMultiEta-etaIndex)*lowEta_lowPhi_corner   + etaIndex*highEta_lowPhi_corner)/granularityMultiEta);
           TVector2 aveEta_highPhi_corner(((granularityMultiEta-etaIndex)*lowEta_highPhi_corner  + etaIndex*highEta_highPhi_corner)/granularityMultiEta);
           TVector2 aveEta_avePhi_corner (((granularityMultiPhi-phiIndex)*aveEta_lowPhi_corner   + phiIndex*aveEta_highPhi_corner )/granularityMultiPhi);
           subcrystal_eta_edges[etaIndex][phiIndex]=aveEta_avePhi_corner.X();
           subcrystal_phi_edges[etaIndex][phiIndex]=aveEta_avePhi_corner.Y();
        }
      }

      float subcrystal_eta_centers[granularityMultiEta][granularityMultiPhi];
      float subcrystal_phi_centers[granularityMultiEta][granularityMultiPhi];

      for (unsigned int etaIndex=0; etaIndex<granularityMultiEta; etaIndex++)
      {
        for (unsigned int phiIndex=0; phiIndex<granularityMultiPhi; phiIndex++)
        {
          float centerEta = (subcrystal_eta_edges[etaIndex][phiIndex]+
                             subcrystal_eta_edges[etaIndex+1][phiIndex+1]+
                             subcrystal_eta_edges[etaIndex+1][phiIndex]+
                             subcrystal_eta_edges[etaIndex][phiIndex+1])/4;

          float centerPhi = (subcrystal_phi_edges[etaIndex][phiIndex]+
                             subcrystal_phi_edges[etaIndex+1][phiIndex+1]+
                             subcrystal_phi_edges[etaIndex+1][phiIndex]+
                             subcrystal_phi_edges[etaIndex][phiIndex+1])/4;

          subcrystal_eta_centers[etaIndex][phiIndex]=centerEta;
          subcrystal_phi_centers[etaIndex][phiIndex]=centerPhi;          

        }
      }

      float minSubDist = 999;
      unsigned int subcrystal_phi_index=0;
      unsigned int subcrystal_eta_index=0;
      for (unsigned int etaIndex=0; etaIndex<granularityMultiEta; etaIndex++)
      {
        for (unsigned int phiIndex=0; phiIndex<granularityMultiPhi; phiIndex++)
        {
          float d=reco::deltaR(eta,phi,subcrystal_eta_centers[etaIndex][phiIndex],subcrystal_phi_centers[etaIndex][phiIndex]);
          if (d<minSubDist)
          {
            minSubDist = d;
            subcrystal_eta_index=etaIndex;
            subcrystal_phi_index=phiIndex;
          }
        }
      }

      //std::cout<<subcrystal_phi_index<<" "<<subcrystal_eta_index<<"\n";

      EBDetId ebId( id );
      // hEvt_Adj_tracksPt->SetBinContent(  ebId.iphi() - 1 +1, (ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta()) +141,
      //   hEvt_Adj_tracksPt->GetBinContent(ebId.iphi() - 1 +1, (ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta()) +141)+trackPt_ );
      int phi_base_coordinate = ebId.iphi() - 1 +1;
      int eta_base_coordinate =(ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta()) +141;
      phi_base_coordinate = (phi_base_coordinate - 1) * granularityMultiPhi + subcrystal_phi_index+1;
      eta_base_coordinate = (eta_base_coordinate - 1) * granularityMultiEta + subcrystal_eta_index+1;

      std::vector<int> return_vector;
      return_vector.push_back(phi_base_coordinate);
      return_vector.push_back(eta_base_coordinate);
      return return_vector;
}


void RecHitAnalyzer::fillByBinNumber(TH2F * histo, const std::vector<int>& phi_eta, const float& value)
{
  histo->SetBinContent(phi_eta[0],phi_eta[1],
  histo->GetBinContent(phi_eta[0],phi_eta[1])+value);
}


// Fill adjustable EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillTracksAtECALadjustable ( const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int proj ) {

  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, trackPt_, trackQPt_,trackD0_, trackDz_;
  GlobalPoint pos;

  edm::ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);

  vECAL_tracks_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  vECAL_tracksPt_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  vECAL_tracksQPt_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  vECAL_tracksD0_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  vECAL_tracksDz_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  vECAL_tracksPt_max_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  vECAL_tracksD0_max_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  vECAL_tracksDz_max_[proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
  hEvt_Adj_tracks[proj]->Reset();
  hEvt_Adj_tracksPt[proj]->Reset();
  hEvt_Adj_tracksQPt[proj]->Reset();
  hEvt_Adj_tracksD0[proj]->Reset();
  hEvt_Adj_tracksDz[proj]->Reset();
  hEvt_Adj_tracksPt_max[proj]->Reset();
  hEvt_Adj_tracksD0_max[proj]->Reset();
  hEvt_Adj_tracksDz_max[proj]->Reset();

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByLabel( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByLabel( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  edm::Handle<reco::TrackCollection> tracksH_;
  iEvent.getByLabel( trackCollectionT_, tracksH_ );

  bool isPVgood=false;
  edm::Handle<reco::VertexCollection> pvColl;
  iEvent.getByLabel(pvCollectionT_, pvColl);
  isPVgood = pvColl.product()->size()>0;
  reco::Vertex the_PV;
  if (isPVgood) the_PV = pvColl.product()->at(0);


  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  int bin;

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) { 
    if ( !(iTk->quality(tkQt_)) ) continue;
    if(iTk->pt()<=0.5)continue; 
    if(iTk->charge()==0) continue;

    bool isPropagationOk=false;
    eta = 0.;
    phi = 0.;

    auto propagatedECALTrack = spr::propagateTrackToECAL(&*iTk, magfield.product());
    isPropagationOk=propagatedECALTrack.ok;
    if (propagatedECALTrack.ok)
    {
      eta = propagatedECALTrack.direction.eta();
      phi = propagatedECALTrack.direction.phi();
    }
    
    if ( std::abs(eta) > 3. || !isPropagationOk ) continue;

    trackPt_ = iTk->pt();
    trackQPt_ = iTk->pt()*iTk->charge();

    trackD0_ = -iTk->dxy(the_PV.position());
    trackDz_ = iTk->dz(the_PV.position());
 
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalEndcap )
    {
      //int phi2;
      float kappa= 4*TMath::Pi()/HBHE_IPHI_NUM;
      if (phi>-kappa)
        phi=phi+kappa-TMath::Pi();
      else  
        phi=phi+kappa+TMath::Pi();

      hEvt_Adj_tracks[proj]->Fill( phi, eta );
      hEvt_Adj_tracksPt[proj]->Fill( phi, eta, trackPt_ );
      hEvt_Adj_tracksQPt[proj]->Fill( phi, eta, trackQPt_ );
      hEvt_Adj_tracksD0[proj]->Fill( phi, eta, trackD0_ );
      hEvt_Adj_tracksDz[proj]->Fill( phi, eta, trackDz_ );
      bin = hEvt_Adj_tracks[proj]->FindBin( phi, eta );
      if ( trackPt_ > hEvt_Adj_tracksPt_max[proj]->GetBinContent( bin ) ) {
        hEvt_Adj_tracksPt_max[proj]->SetBinContent( bin, trackPt_ );
        hEvt_Adj_tracksD0_max[proj]->SetBinContent( bin, trackD0_ );
        hEvt_Adj_tracksDz_max[proj]->SetBinContent( bin, trackDz_ );
      }

      hECALadj_tracks[proj]->Fill( phi, eta, 1.0 );
      hECALadj_tracksPt[proj]->Fill( phi, eta, trackPt_ );
      hECALadj_tracksQPt[proj]->Fill( phi, eta, trackQPt_ );
      hECALadj_tracksD0[proj]->Fill( phi, eta, trackD0_ );
      hECALadj_tracksDz[proj]->Fill( phi, eta, trackDz_ );

    }
    else if ( id.subdetId() == EcalBarrel ) { 

      std::vector<int> phi_eta = findSubcrystal(caloGeom, eta, phi, granularityMultiEta[proj], granularityMultiPhi[proj]);

      fillByBinNumber(hEvt_Adj_tracks[proj], phi_eta, 1.0);
      fillByBinNumber(hEvt_Adj_tracksPt[proj], phi_eta, trackPt_);
      fillByBinNumber(hEvt_Adj_tracksQPt[proj], phi_eta, trackQPt_);
      fillByBinNumber(hEvt_Adj_tracksD0[proj], phi_eta, trackD0_);
      fillByBinNumber(hEvt_Adj_tracksDz[proj], phi_eta, trackDz_);

      if ( trackPt_ > hEvt_Adj_tracksPt_max[proj]->GetBinContent( phi_eta[0],phi_eta[1] ) ) {
        hEvt_Adj_tracksPt_max[proj]->SetBinContent( phi_eta[0],phi_eta[1], trackPt_ );
        hEvt_Adj_tracksD0_max[proj]->SetBinContent( phi_eta[0],phi_eta[1], trackD0_ );
        hEvt_Adj_tracksDz_max[proj]->SetBinContent( phi_eta[0],phi_eta[1], trackDz_ );
      }


      fillByBinNumber(hECALadj_tracks[proj], phi_eta, 1.0);
      fillByBinNumber(hECALadj_tracksPt[proj], phi_eta, trackPt_);
      fillByBinNumber(hECALadj_tracksQPt[proj], phi_eta, trackQPt_);
      fillByBinNumber(hECALadj_tracksD0[proj], phi_eta, trackD0_);
      fillByBinNumber(hECALadj_tracksDz[proj], phi_eta, trackDz_);

  } //tracks loop

  int index1d=0;
  for ( int ieta=1; ieta<=totalEtaBins[proj]; ieta++ )
  {
    for ( int iphi=1; iphi<=totalPhiBins[proj]; iphi++ )
    {
      index1d= (ieta-1)*totalPhiBins[proj]+iphi-1;//ieta_global*EB_IPHI_MAX + iphi_; 

      vECALadj_tracksPt_[proj][index1d]=hEvt_Adj_tracksPt->GetBinContent(iphi,ieta);
      vECALadj_tracksQPt_[proj][index1d]=hEvt_Adj_tracksQPt->GetBinContent(iphi,ieta);
      vECALadj_tracksD0_[proj][index1d]=hEvt_Adj_tracksD0->GetBinContent(iphi,ieta);
      vECALadj_tracksDz_[proj][index1d]=hEvt_Adj_tracksDz->GetBinContent(iphi,ieta);
      vECALadj_tracks_[proj][index1d]=hEvt_Adj_tracks->GetBinContent(iphi,ieta);
      vECALadj_tracksPt_max_[proj][index1d]=hEvt_Adj_tracksPt_max->GetBinContent(iphi,ieta);
      vECALadj_tracksD0_max_[proj][index1d]=hEvt_Adj_tracksD0_max->GetBinContent(iphi,ieta);
      vECALadj_tracksDz_max_[proj][index1d]=hEvt_Adj_tracksDz_max->GetBinContent(iphi,ieta);

    }
  }

  // Get average D0 and Dz for each position
  for (unsigned int idx_=0;idx_<vECAL_tracks_[proj].size();idx_++) {
    if (vECAL_tracks_[proj][idx_] != 0) {
      vECAL_tracksD0_[proj][idx_] = vECAL_tracksD0_[proj][idx_] / vECAL_tracks_[proj][idx_];
      vECAL_tracksDz_[proj][idx_] = vECAL_tracksDz_[proj][idx_] / vECAL_tracks_[proj][idx_];
    }
  }

} // fillTracksAtECALadjustable()
