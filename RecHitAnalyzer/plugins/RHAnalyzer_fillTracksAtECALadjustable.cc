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

TProfile2D *hECALadj_tracks[Nadjproj];
TProfile2D *hECALadj_tracksPt[Nadjproj];
TProfile2D *hECALadj_tracksQPt[Nadjproj];
TProfile2D *hECALadj_tracksD0[Nadjproj];
TProfile2D *hECALadj_tracksDz[Nadjproj];

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
    hECALadj_tracks[proj] = fs->make<TProfile2D>((std::string("ECALadj_tracks")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksPt[proj] = fs->make<TProfile2D>((std::string("ECALadj_tracksPt")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksQPt[proj] = fs->make<TProfile2D>((std::string("ECALadj_tracksQPt")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksD0[proj] = fs->make<TProfile2D>((std::string("ECALadj_tracksD0")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
        totalPhiBins, -TMath::Pi(), TMath::Pi(),
        adjEtaBins.size()-1, &adjEtaBins[0] );
    hECALadj_tracksDz[proj] = fs->make<TProfile2D>((std::string("ECALadj_tracksDz")+adj_projections[proj]).c_str(), "E(#phi,#eta);#phi;#eta",
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

// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
//void fillTracksAtECAL_with_EEproj ( TH2F *hEvt_EE_tracks_, TH2F *hEvt_EE_tracksPt_, TH2F *hEvt_EE_tracksQPt_, TH2F *hEvt_EE_tracksD0_, TH2F *hEvt_EE_tracksDz_,  TH2F *hEvt_EE_tracksPt_max_, TH2F *hEvt_EE_tracksD0_max_, TH2F *hEvt_EE_tracksDz_max_, int ieta_global_offset, int ieta_signed_offset ) {
// void fillTracksAtECAL_with_EEproj ( int side, int ieta_global_offset, int ieta_signed_offset, int proj ) {

//   int ieta_global_, ieta_signed_;
//   int ieta_, iphi_, idx_;
//   float track_;
//   float trackPt_;
//   float trackQPt_;
//   float trackD0_; 
//   float trackDz_;
//   float trackPt_max_;
//   float trackD0_max_;
//   float trackDz_max_;
//   for (int ieta = 1; ieta < hEvt_EE_tracksPt[side]->GetNbinsY()+1; ieta++) {
//     ieta_ = ieta - 1;
//     ieta_global_ = ieta_ + ieta_global_offset;
//     ieta_signed_ = ieta_ + ieta_signed_offset;
//     for (int iphi = 1; iphi < hEvt_EE_tracksPt[side]->GetNbinsX()+1; iphi++) {

//       track_  = hEvt_EE_tracks[side]->GetBinContent( iphi, ieta );
//       trackPt_ = hEvt_EE_tracksPt[side]->GetBinContent( iphi, ieta );
//       trackQPt_ = hEvt_EE_tracksQPt[side]->GetBinContent( iphi, ieta );
//       trackD0_ = hEvt_EE_tracksD0[side]->GetBinContent( iphi, ieta );
//       trackDz_ = hEvt_EE_tracksDz[side]->GetBinContent( iphi, ieta );
//       trackPt_max_ = hEvt_EE_tracksPt_max[side]->GetBinContent( iphi, ieta );
//       trackD0_max_ = hEvt_EE_tracksD0_max[side]->GetBinContent( iphi, ieta );
//       trackDz_max_ = hEvt_EE_tracksDz_max[side]->GetBinContent( iphi, ieta );
//       if ( trackPt_ <= zs ) continue;
//       // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
//       iphi_ = iphi  + 5*38; // shift
//       iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
//       iphi_ = iphi_ - 1;
//       idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
//       // Fill vector for image
//       vECAL_tracks_[proj][idx_] = track_;
//       vECAL_tracksPt_[proj][idx_] = trackPt_;
//       vECAL_tracksQPt_[proj][idx_] = trackQPt_;
//       vECAL_tracksD0_[proj][idx_] = trackD0_;
//       vECAL_tracksDz_[proj][idx_] = trackDz_;
//       vECAL_tracksPt_max_[proj][idx_] = trackPt_max_;
//       vECAL_tracksD0_max_[proj][idx_] = trackD0_max_;
//       vECAL_tracksDz_max_[proj][idx_] = trackDz_max_;
//       // Fill histogram for monitoring
//       hECAL_tracks[proj]->Fill( iphi_, ieta_signed_, track_ );
//       hECAL_tracksPt[proj]->Fill( iphi_, ieta_signed_, trackPt_ );
//       hECAL_tracksQPt[proj]->Fill( iphi_, ieta_signed_, trackQPt_ );
//       hECAL_tracksD0[proj]->Fill( iphi_, ieta_signed_, trackD0_ );
//       hECAL_tracksDz[proj]->Fill( iphi_, ieta_signed_, trackDz_ );

//     } // iphi_
//   } // ieta_

// } // fillTracksAtECAL_with_EEproj



std::vector<int> findSubcrystal(const CaloGeometry* caloGeom, const float& eta, const float& phi, const int& granularityMultiEta, const int& granularityMultiPhi)
{

    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel )
    { 
      auto subDetGeometry = caloGeom->getSubdetectorGeometry(id);
      auto caloCellGeometry = subDetGeometry->getGeometry(id);
      auto corners = caloCellGeometry->getCornersREP();
      // auto reference= caloCellGeometry->getPosition();
      // auto reference_phi = reference.phi();
      // auto reference_eta = reference.eta();
      float kappa= 4*TMath::Pi()/HBHE_IPHI_NUM;
      // if (reference_phi>-kappa)
      //   reference_phi=reference_phi+kappa-TMath::Pi();
      // else  
      //   reference_phi=reference_phi+kappa+TMath::Pi();
      if (phi>-kappa)
        phi=phi+kappa-TMath::Pi();
      else  
        phi=phi+kappa+TMath::Pi();

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



// Fill adjustable EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillTracksAtECALadjustable ( const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int proj ) {

  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, trackPt_, trackQPt_,trackD0_, trackDz_;
  GlobalPoint pos;

  edm::ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);

  vECAL_tracks_[proj].assign( totalEtaBins*totalPhiBins, 0. );
  vECAL_tracksPt_[proj].assign( totalEtaBins*totalPhiBins, 0. );
  vECAL_tracksQPt_[proj].assign( totalEtaBins*totalPhiBins, 0. );
  vECAL_tracksD0_[proj].assign( totalEtaBins*totalPhiBins, 0. );
  vECAL_tracksDz_[proj].assign( totalEtaBins*totalPhiBins, 0. );
  vECAL_tracksPt_max_[proj].assign( totalEtaBins*totalPhiBins, 0. );
  vECAL_tracksD0_max_[proj].assign( totalEtaBins*totalPhiBins, 0. );
  vECAL_tracksDz_max_[proj].assign( totalEtaBins*totalPhiBins, 0. );
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
  // for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
  //       iTk != tracksH_->end(); ++iTk ) {

  //   if(iTk->pt()<=0.5)continue; 
  //   if(iTk->charge()==0) continue;// NO neutral objects
  //   if ( !(iTk->quality(tkQt_)) ) continue;

  //   if (proj==4)
  //   {
  //     bool pv_match=false;
  //     for ( reco::Vertex::trackRef_iterator iTkPV = the_PV.tracks_begin();
  //           iTkPV != the_PV.tracks_end(); ++iTkPV ) {
  //       if (fabs(iTk->pt()-iTkPV->get()->pt())<0.001 &&
  //           fabs(iTk->eta()-iTkPV->get()->eta())<0.001 &&
  //           fabs(iTk->phi()-iTkPV->get()->phi())<0.001)
  //       {
  //         pv_match=true;
  //         break;
  //       }
  //     }
  //     if (!pv_match)continue;
  //   }    


  //   bool isPropagationOk=false;
  //   eta = 0.;
  //   phi = 0.;
  //   switch (proj)
  //   {
  //     case 0:
  //     {
  //       eta = iTk->eta();
  //       phi = iTk->phi();
  //       isPropagationOk=true;
  //     }
  //     break;

  //     case 1: case 3: case 4:
  //     {
  //       auto propagatedECALTrack = spr::propagateTrackToECAL(&*iTk, magfield.product());
  //       isPropagationOk=propagatedECALTrack.ok;
  //       if (propagatedECALTrack.ok)
  //       {
  //         eta = propagatedECALTrack.direction.eta();
  //         phi = propagatedECALTrack.direction.phi();
  //       }
  //     }
  //     break;

  //     case 2:
  //     {
  //       auto propagatedHCALTrack = spr::propagateTrackToHCAL(&*iTk, magfield.product());
  //       isPropagationOk=propagatedHCALTrack.ok;
  //       if (propagatedHCALTrack.ok)
  //       {
  //         eta = propagatedHCALTrack.direction.eta();
  //         phi = propagatedHCALTrack.direction.phi();
  //       }
  //     }
  //     break;

  //     default:
  //     {
  //       isPropagationOk=false;
  //     }
  //     break;
  //   }
    
  //   if ( std::abs(eta) > 3. || !isPropagationOk) continue;
  //   DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
  //   if ( id.subdetId() == EcalBarrel ) continue;
  //   if ( id.subdetId() == EcalEndcap ) {
  //     iz_ = (eta > 0.) ? 1 : 0;
  //     switch(proj)
  //     {
  //       case 0: case 1: case 2:
  //       {
  //         trackD0_ = iTk->d0();
  //         trackDz_ = iTk->dz();
  //         break;
  //       }
  //       case 3: case 4: default:
  //       {
  //         trackD0_ = -iTk->dxy(the_PV.position());
  //         trackDz_ = iTk->dz(the_PV.position());
  //         break;
  //       }
  //     }
  //     // Fill intermediate helper histogram by eta,phi
  //     hEvt_EE_tracks[iz_]->Fill( phi, eta );
  //     hEvt_EE_tracksPt[iz_]->Fill( phi, eta, iTk->pt() );
  //     hEvt_EE_tracksQPt[iz_]->Fill( phi, eta, iTk->pt()*iTk->charge() );
  //     hEvt_EE_tracksD0[iz_]->Fill( phi, eta, trackD0_ );
  //     hEvt_EE_tracksDz[iz_]->Fill( phi, eta, trackDz_ );
  //     bin = hEvt_EE_tracks[iz_]->FindBin( phi, eta );
  //     if ( iTk->pt() > hEvt_EE_tracksPt_max[iz_]->GetBinContent( bin ) ) {
  //       hEvt_EE_tracksPt_max[iz_]->SetBinContent( bin, iTk->pt() );
  //       hEvt_EE_tracksD0_max[iz_]->SetBinContent( bin, trackD0_ );
  //       hEvt_EE_tracksDz_max[iz_]->SetBinContent( bin, trackDz_ );
  //     }
  //   }
  // } // tracks

  // // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  // ieta_global_offset = 0;
  // ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  // fillTracksAtECAL_with_EEproj( 0, ieta_global_offset, ieta_signed_offset, proj );

  // // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  // ieta_global_offset = 55;

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) { 
    if ( !(iTk->quality(tkQt_)) ) continue;
    if(iTk->pt()<=0.5)continue; 
    if(iTk->charge()==0) continue;


    if (proj==4)
    {
      bool pv_match=false;
      for ( reco::Vertex::trackRef_iterator iTkPV = the_PV.tracks_begin();
            iTkPV != the_PV.tracks_end(); ++iTkPV ) {
        if (fabs(iTk->pt()-iTkPV->get()->pt())<0.001 &&
            fabs(iTk->eta()-iTkPV->get()->eta())<0.001 &&
            fabs(iTk->phi()-iTkPV->get()->phi())<0.001)
        {
          pv_match=true;
          break;
        }
      }
      if (!pv_match)continue;
    }


    bool isPropagationOk=false;
    eta = 0.;
    phi = 0.;
    switch (proj)
    {
      case 0:
      {
        eta = iTk->eta();
        phi = iTk->phi();
        isPropagationOk=true;
      }
      break;

      case 1: case 3: case 4:
      {
        auto propagatedECALTrack = spr::propagateTrackToECAL(&*iTk, magfield.product());
        isPropagationOk=propagatedECALTrack.ok;
        if (propagatedECALTrack.ok)
        {
          eta = propagatedECALTrack.direction.eta();
          phi = propagatedECALTrack.direction.phi();
        }
      }
      break;

      case 2:
      {
        auto propagatedHCALTrack = spr::propagateTrackToHCAL(&*iTk, magfield.product());
        isPropagationOk=propagatedHCALTrack.ok;
        if (propagatedHCALTrack.ok)
        {
          eta = propagatedHCALTrack.direction.eta();
          phi = propagatedHCALTrack.direction.phi();
        }
      }
      break;

      default:
      {
        isPropagationOk=false;
      }
      break;
    }
    
    if ( std::abs(eta) > 3. || !isPropagationOk) continue;

    trackPt_ = iTk->pt();
    trackQPt_ = iTk->pt()*iTk->charge();
    switch(proj)
    {
      case 0: case 1: case 2:
      {
        trackD0_ = iTk->d0();
        trackDz_ = iTk->dz();
        break;
      }
      case 3: case 4: default:
      {
        trackD0_ = -iTk->dxy(the_PV.position());
        trackDz_ = iTk->dz(the_PV.position());
        break;
      }
    }
 //qui   
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalEndcap )
    {
      int phi2;
      float kappa= 4*TMath::Pi()/HBHE_IPHI_NUM;
      if (phi>-kappa)
        phi2=phi+kappa-TMath::Pi();
      else  
        phi2=phi+kappa+TMath::Pi();

      hEvt_Adj_tracksPt->Fill(  phi2, eta, trackPt_ );
      hEvt_Adj_tracksQPt->Fill( phi2, eta, trackPt_ );
      //hEvt_Adj_tracksQPt->Fill( phi, eta, trackQPt_ );
      hECALadj_tracksPt->Fill(  phi2, eta, trackPt_ );
      hECALadj_tracksQPt->Fill( phi2, eta, trackPt_ );


    }
    else if ( id.subdetId() == EcalBarrel ) { 
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      if ( trackPt_ <= zs ) continue;
      // Fill vector for image
      ieta_signed = ieta_;
      ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
      idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
      vECAL_tracks_[proj][idx_] += 1.;
      vECAL_tracksPt_[proj][idx_] += trackPt_;
      vECAL_tracksQPt_[proj][idx_] += trackQPt_;
      vECAL_tracksD0_[proj][idx_] += trackD0_;
      vECAL_tracksDz_[proj][idx_] += trackDz_;
      if ( trackPt_ > vECAL_tracksPt_max_[proj][idx_] ) {
        vECAL_tracksPt_max_[proj][idx_] = trackPt_;
        vECAL_tracksD0_max_[proj][idx_] = trackD0_;
        vECAL_tracksDz_max_[proj][idx_] = trackDz_;
      }
      // Fill histogram for monitoring
      hECAL_tracks[proj]->Fill( iphi_, ieta_signed, 1.0 );
      hECAL_tracksPt[proj]->Fill( iphi_, ieta_signed, trackPt_ );
      hECAL_tracksQPt[proj]->Fill( iphi_, ieta_signed, trackQPt_ );
      hECAL_tracksD0[proj]->Fill( iphi_, ieta_signed, trackD0_ );
      hECAL_tracksDz[proj]->Fill( iphi_, ieta_signed, trackDz_ );
    }

  } // EB
  // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  // ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
  // ieta_signed_offset = EB_IETA_MAX;
  // fillTracksAtECAL_with_EEproj( 1, ieta_global_offset, ieta_signed_offset, proj );

  // Get average D0 and Dz for each position
  for (unsigned int idx_=0;idx_<vECAL_tracks_[proj].size();idx_++) {
    if (vECAL_tracks_[proj][idx_] != 0) {
      vECAL_tracksD0_[proj][idx_] = vECAL_tracksD0_[proj][idx_] / vECAL_tracks_[proj][idx_];
      vECAL_tracksDz_[proj][idx_] = vECAL_tracksDz_[proj][idx_] / vECAL_tracks_[proj][idx_];
    }
  }

} // fillTracksAtECALadjustable()
