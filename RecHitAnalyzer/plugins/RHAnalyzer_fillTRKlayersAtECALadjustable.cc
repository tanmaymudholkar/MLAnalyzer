#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
//#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"


// Fill TRK rec hits ////////////////////////////////
// by layer at ECAL adjustable

TH2F *hTOB_ECALadj[nTOB][Nadjproj];
std::vector<float> vTOB_ECALadj_[nTOB][Nadjproj];
TH2F *hEvt_Adj_TOB[nTOB][Nadjproj];

TH2F *hTEC_ECALadj[nTEC][Nadjproj];
std::vector<float> vTEC_ECALadj_[nTEC][Nadjproj];
TH2F *hEvt_Adj_TEC[nTEC][Nadjproj];

TH2F *hTIB_ECALadj[nTIB][Nadjproj];
std::vector<float> vTIB_ECALadj_[nTIB][Nadjproj];
TH2F *hEvt_Adj_TIB[nTIB][Nadjproj];

TH2F *hTID_ECALadj[nTID][Nadjproj];
std::vector<float> vTID_ECALadj_[nTID][Nadjproj];
TH2F *hEvt_Adj_TID[nTID][Nadjproj];

TH2F *hBPIX_ECALadj[nBPIX][Nadjproj];
std::vector<float> vBPIX_ECALadj_[nBPIX][Nadjproj];
TH2F *hEvt_Adj_BPIX[nBPIX][Nadjproj];

TH2F *hFPIX_ECALadj[nFPIX][Nadjproj];
std::vector<float> vFPIX_ECALadj_[nFPIX][Nadjproj];
TH2F *hEvt_Adj_FPIX[nFPIX][Nadjproj];

// Initialize branches ____________________________________________________________//
void RecHitAnalyzer::branchesTRKlayersAtECALadjustable ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images

  int layer;
  char hname[50], htitle[50];
  //const double * eta_bins_EE[2] = {eta_bins_EEm,eta_bins_EEp};

  for ( unsigned int proj=0; proj < Nadjproj; proj++ ) {
  //TOB
  
  	for ( int iL(0); iL < nTOB; iL++ ) {
    	// Branches for images
    	layer = iL + 1;
    	sprintf(hname, "TOB_layer%d_ECALadj%s",layer,adj_projections[proj].c_str());
    	tree->Branch(hname,        &vTOB_ECALadj_[iL][proj]);

    	// Histograms for monitoring
    	sprintf(htitle,"N(#phi,#eta);#phi;#eta");
    	hTOB_ECALadj[iL][proj] = fs->make<TH2F>(hname, htitle,
        	totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

      sprintf(hname, "evt_TOB_layer%d_ECALAdj%s",layer,adj_projections[proj].c_str());
      hEvt_Adj_TOB[iL][proj] = new TH2F(hname, htitle,
      		totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );
  	} // iL
  

  //TEC
	 for ( int iL(0); iL < nTEC; iL++ ) {
	    // Branches for images
	    layer = iL + 1;
	    sprintf(hname, "TEC_layer%d_ECALadj%s",layer,adj_projections[proj].c_str());
	    tree->Branch(hname,        &vTEC_ECALadj_[iL][proj]);

	    // Histograms for monitoring
	    sprintf(htitle,"N(#phi,#eta);#phi;#eta");
	    hTEC_ECALadj[iL][proj] = fs->make<TH2F>(hname, htitle,
	        totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

		  sprintf(hname, "evt_TEC_layer%d_ECALAdj%s",layer,adj_projections[proj].c_str());
		  hEvt_Adj_TEC[iL][proj] = new TH2F(hname, htitle,
		      totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

	  } // iL

	  //TIB
	  for ( int iL(0); iL < nTIB; iL++ ) {
	    // Branches for images
	    layer = iL + 1;
	    sprintf(hname, "TIB_layer%d_ECALadj%s",layer,adj_projections[proj].c_str());
	    tree->Branch(hname,        &vTIB_ECALadj_[iL][proj]);

	    // Histograms for monitoring
	    sprintf(htitle,"N(#phi,#eta);#phi;#eta");
	    hTIB_ECALadj[iL][proj] = fs->make<TH2F>(hname, htitle,
	        totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

		  sprintf(hname, "evt_TIB_layer%d_ECALAdj%s",layer,adj_projections[proj].c_str());
		  hEvt_Adj_TIB[iL][proj] = new TH2F(hname, htitle,
		      totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

	  } // iL

  //TID
	  for ( int iL(0); iL < nTID; iL++ ) {
	    // Branches for images
	    layer = iL + 1;
	    sprintf(hname, "TID_layer%d_ECALadj%s",layer,adj_projections[proj].c_str());
	    tree->Branch(hname,        &vTID_ECALadj_[iL][proj]);

	    // Histograms for monitoring
	    sprintf(htitle,"N(#phi,#eta);#phi;#eta");
	    hTID_ECALadj[iL][proj] = fs->make<TH2F>(hname, htitle,
	        totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

		  sprintf(hname, "evt_TID_layer%d_ECALAdj%s",layer,adj_projections[proj].c_str());
		  hEvt_Adj_TID[iL][proj] = new TH2F(hname, htitle,
		      totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

	  } // iL

	  //BPIX
	  for ( int iL(0); iL < nBPIX; iL++ ) {
	    // Branches for images
	    layer = iL + 1;
	    sprintf(hname, "BPIX_layer%d_ECALadj%s",layer,adj_projections[proj].c_str());
	    tree->Branch(hname,        &vBPIX_ECALadj_[iL][proj]);

	    // Histograms for monitoring
	    sprintf(htitle,"N(#phi,#eta);#phi;#eta");
	    hBPIX_ECALadj[iL][proj] = fs->make<TH2F>(hname, htitle,
	        totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

		  sprintf(hname, "evt_BPIX_layer%d_ECALAdj%s",layer,adj_projections[proj].c_str());
		  hEvt_Adj_BPIX[iL][proj] = new TH2F(hname, htitle,
		      totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

	  } // iL

  //FPIX
	  for ( int iL(0); iL < nFPIX; iL++ ) {
	    // Branches for images
	    layer = iL + 1;
	    sprintf(hname, "FPIX_layer%d_ECALadj%s",layer,adj_projections[proj].c_str());
	    tree->Branch(hname,        &vFPIX_ECALadj_[iL][proj]);

	    // Histograms for monitoring
	    sprintf(htitle,"N(#phi,#eta);#phi;#eta");
	    hFPIX_ECALadj[iL][proj] = fs->make<TH2F>(hname, htitle,
	        totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

		  sprintf(hname, "evt_FPIX_layer%d_ECALAdj%s",layer,adj_projections[proj].c_str());
		  hEvt_Adj_FPIX[iL][proj] = new TH2F(hname, htitle,
		      totalPhiBins[proj], -TMath::Pi(), TMath::Pi(),
          adjEtaBins[proj].size()-1, &adjEtaBins[proj][0] );

	  } // iL
}//proj

} // branchesEB()

void RecHitAnalyzer::fillTRKlayerHelper (int layer_, unsigned int proj, TH2F *hSUBDET_ECAL[][Nadjproj], TH2F *hEvt_Adj_SUBDET[][Nadjproj], const CaloGeometry* caloGeom, const float& eta, const float& phi)
{
  DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
  if ( ecalId.subdetId() == EcalBarrel )
  {
    std::vector<int> phi_eta = findSubcrystal(caloGeom, eta, phi, granularityMultiEta[proj], granularityMultiPhi[proj]);
    fillByBinNumber(hSUBDET_ECAL[layer_-1][proj], phi_eta, 1.0);
    fillByBinNumber(hEvt_Adj_SUBDET[layer_-1][proj], phi_eta, 1.0);
  }
  else if ( ecalId.subdetId() == EcalEndcap )
  {
    float phi2;
    float kappa= 4*TMath::Pi()/HBHE_IPHI_NUM;
    if (phi>-kappa)
      phi2=phi+kappa-TMath::Pi();
    else  
      phi2=phi+kappa+TMath::Pi();

    //std::cout<<"lay:"<<layer_<<" proj:"<<proj<<std::endl;
    //std::cout<<"phi2:"<<phi2<<" eta:"<<eta<<std::endl;
    hSUBDET_ECAL[layer_-1][proj]->Fill( phi2, eta );
    hEvt_Adj_SUBDET[layer_-1][proj]->Fill( phi2, eta);
  }
}

unsigned int RecHitAnalyzer::getLayer(const DetId& detid)
{
  unsigned int subid=detid.subdetId();

          switch(subid){
            case 1:{//BPIX
              PXBDetId pdetId = PXBDetId(detid);
              return pdetId.layer();
            }break;

            case 2:{//FPIX
              PXFDetId pdetId = PXFDetId(detid.rawId());
              return pdetId.disk();
            }break;

            case 3:{//TIB
              TIBDetId pdetId = TIBDetId(detid);
              return pdetId.layer();
            }break;

            case 4:{//TID
              TIDDetId pdetId = TIDDetId(detid);
              return pdetId.wheel();
            }break;

            case 5:{//TOB
              TOBDetId pdetId = TOBDetId(detid);
              return pdetId.layer();
            }break;

            case 6:{//TEC
              TECDetId pdetId = TECDetId(detid);
              return pdetId.wheel();
            }break;
          }
          return 999;
}


// Fill TRK rechits at ECAL adjustable ______________________________________________________________//
void RecHitAnalyzer::fillTRKlayersAtECALadjustable ( const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int proj ) {

  float eta, phi;
  GlobalPoint pos;

  for ( int iL(0); iL < nTOB; iL++ ) {
    vTOB_ECALadj_[iL][proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
    hEvt_Adj_TOB[iL][proj]->Reset();
  }
  for ( int iL(0); iL < nTEC; iL++ ) {
    vTEC_ECALadj_[iL][proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
    hEvt_Adj_TEC[iL][proj]->Reset();
  }
  for ( int iL(0); iL < nTIB; iL++ ) {
    vTIB_ECALadj_[iL][proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
    hEvt_Adj_TIB[iL][proj]->Reset();
  }
  for ( int iL(0); iL < nTID; iL++ ) {
    vTID_ECALadj_[iL][proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
    hEvt_Adj_TID[iL][proj]->Reset();
  }
  for ( int iL(0); iL < nBPIX; iL++ ) {
    vBPIX_ECALadj_[iL][proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
    hEvt_Adj_BPIX[iL][proj]->Reset();
  }
  for ( int iL(0); iL < nFPIX; iL++ ) {
    vFPIX_ECALadj_[iL][proj].assign( totalEtaBins[proj]*totalPhiBins[proj], 0. );
    hEvt_Adj_FPIX[iL][proj]->Reset();
  }

  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  bool isPVgood=false;
  edm::Handle<reco::VertexCollection> pvColl;
  iEvent.getByLabel(pvCollectionT_, pvColl);
  isPVgood = pvColl.product()->size()>0;
  reco::Vertex the_PV;
  if (isPVgood) the_PV = pvColl.product()->at(0);
  TVector3 pv_v(the_PV.x(),the_PV.y(),the_PV.z());

//sipixel
  edm::ESHandle<TrackerGeometry> geom;
  iSetup.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);

  edm::Handle<SiPixelRecHitCollection> recHitColl;
  iEvent.getByLabel(siPixelRecHitCollectionT_, recHitColl);

  SiPixelRecHitCollection::const_iterator recHitIdIterator      = (recHitColl.product())->begin();
  SiPixelRecHitCollection::const_iterator recHitIdIteratorEnd   = (recHitColl.product())->end();

  for ( ; recHitIdIterator != recHitIdIteratorEnd; recHitIdIterator++)
  {
    SiPixelRecHitCollection::DetSet detset = *recHitIdIterator;
    DetId detId = DetId(detset.detId()); // Get the Detid object
    unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
    unsigned int layer = getLayer(detId);
    const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );

    SiPixelRecHitCollection::DetSet::const_iterator pixeliter=detset.begin();
    SiPixelRecHitCollection::DetSet::const_iterator rechitRangeIteratorEnd   = detset.end();
    for(;pixeliter!=rechitRangeIteratorEnd;++pixeliter)
    {//loop on the rechit
      if (pixeliter->isValid())
      {
        LocalPoint lp = pixeliter->localPosition();
        GlobalPoint GP = theGeomDet->surface().toGlobal(Local3DPoint(lp));

        TVector3 GP_v(GP.x(),GP.y(),GP.z());
        GP_v=GP_v-pv_v;
        phi=GP_v.Phi();
        eta=GP_v.Eta();

        //if ( std::abs(eta) > 3. ) continue;
        //DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
        if ( subid == PixelSubdetector::PixelBarrel ){
          fillTRKlayerHelper (layer, proj, hBPIX_ECALadj, hEvt_Adj_BPIX, caloGeom, eta, phi);
        }
        else if ( subid == PixelSubdetector::PixelEndcap ){
          fillTRKlayerHelper (layer, proj, hFPIX_ECALadj, hEvt_Adj_FPIX, caloGeom, eta, phi);
        }
      }
    }
  }

//sistrip
  for (const auto & itoken: siStripRecHitCollectionT_)
  {
    edm::Handle<SiStripRecHit2DCollection> stripRecHitColl;
    iEvent.getByLabel( itoken , stripRecHitColl);

    SiStripRecHit2DCollection::const_iterator stripRecHitIdIterator      = (stripRecHitColl.product())->begin();
    SiStripRecHit2DCollection::const_iterator stripRecHitIdIteratorEnd   = (stripRecHitColl.product())->end();

    for (; stripRecHitIdIterator != stripRecHitIdIteratorEnd; ++stripRecHitIdIterator)
    { 
     SiStripRecHit2DCollection::DetSet detset = *stripRecHitIdIterator;
     DetId detId = DetId(detset.detId()); // Get the Detid object
     unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
     unsigned int layer = getLayer(detId);
     const StripGeomDetUnit* theGeomDet = dynamic_cast<const StripGeomDetUnit*>( theTracker.idToDet( detId ) );

     SiStripRecHit2DCollection::DetSet::const_iterator stripiter=detset.begin();
     SiStripRecHit2DCollection::DetSet::const_iterator stripRechitRangeIteratorEnd   = detset.end();
     for(;stripiter!=stripRechitRangeIteratorEnd;++stripiter)
      {
        if (stripiter->isValid())
        {
          LocalPoint lp = stripiter->localPosition();
          GlobalPoint GP = theGeomDet->surface().toGlobal(Local3DPoint(lp));

          TVector3 GP_v(GP.x(),GP.y(),GP.z());
          GP_v=GP_v-pv_v;
          phi=GP_v.Phi();
          eta=GP_v.Eta();

          //if ( std::abs(eta) > 3. ) continue;
          DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
          if ( subid == StripSubdetector::TOB ) {
            fillTRKlayerHelper (layer, proj, hTOB_ECALadj, hEvt_Adj_TOB, caloGeom, eta, phi);
          }
          else if ( subid == StripSubdetector::TEC ) {
            fillTRKlayerHelper (layer, proj, hTEC_ECALadj, hEvt_Adj_TEC, caloGeom, eta, phi);
          }
          else if ( subid == StripSubdetector::TIB ) {
            fillTRKlayerHelper (layer, proj, hTIB_ECALadj, hEvt_Adj_TIB, caloGeom, eta, phi);
          }
          else if ( subid == StripSubdetector::TID ) {
            fillTRKlayerHelper (layer, proj, hTID_ECALadj, hEvt_Adj_TID, caloGeom, eta, phi);
          }
        }
      }
    }
  }

  int index1d=0;
  for ( int ieta=1; ieta<=totalEtaBins[proj]; ieta++ )
  {
    for ( int iphi=1; iphi<=totalPhiBins[proj]; iphi++ )
    {
      index1d= (ieta-1)*totalPhiBins[proj]+iphi-1;//ieta_global*EB_IPHI_MAX + iphi_; 

      for ( int iL(0); iL < nTOB; iL++ ) {
        vTOB_ECALadj_[iL][proj][index1d] = hEvt_Adj_TOB[iL][proj]->GetBinContent(iphi,ieta);
      }
      for ( int iL(0); iL < nTEC; iL++ ) {
        vTEC_ECALadj_[iL][proj][index1d] = hEvt_Adj_TEC[iL][proj]->GetBinContent(iphi,ieta);
      }
      for ( int iL(0); iL < nTIB; iL++ ) {
        vTIB_ECALadj_[iL][proj][index1d] = hEvt_Adj_TIB[iL][proj]->GetBinContent(iphi,ieta);
      }
      for ( int iL(0); iL < nTID; iL++ ) {
        vTID_ECALadj_[iL][proj][index1d] = hEvt_Adj_TID[iL][proj]->GetBinContent(iphi,ieta);
      }
      for ( int iL(0); iL < nBPIX; iL++ ) {
        vBPIX_ECALadj_[iL][proj][index1d] = hEvt_Adj_BPIX[iL][proj]->GetBinContent(iphi,ieta);
      }
      for ( int iL(0); iL < nFPIX; iL++ ) {
        vFPIX_ECALadj_[iL][proj][index1d] = hEvt_Adj_FPIX[iL][proj]->GetBinContent(iphi,ieta);
      }

    }
  }

} // fillEB()
