#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Fill EB rec hits ////////////////////////////////
// Store event rechits in a vector of length equal
// to number of crystals in EB (ieta:170 x iphi:360)

TProfile2D *hEB_simEnergy;
//TProfile2D *hEB_energy;
//TProfile2D *hEB_time;
std::vector<float> vEB_simEnergy_;
//std::vector<float> vEB_energy_;
//std::vector<float> vEB_time_;
unsigned int nRecHits_, nSimHits_, nDigis_;
std::vector<float> vEB_adc_[EcalDataFrame::MAXSAMPLES];
TProfile2D *hEB_adc[EcalDataFrame::MAXSAMPLES];

// Initialize branches _____________________________________________________//
void SCRegressor::branchesEB ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("EB_simEnergy", &vEB_simEnergy_);
  tree->Branch("EB_energy", &vEB_energy_);
  tree->Branch("EB_time",   &vEB_time_);
  tree->Branch("nRecHits",     &nRecHits_);
  tree->Branch("nDigis",       &nDigis_);
  tree->Branch("nSimHits",     &nSimHits_);

  // Histograms for monitoring
  hEB_energy = fs->make<TProfile2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hEB_time = fs->make<TProfile2D>("EB_time", "t(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );

  char hname[50], htitle[50];
  for ( int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS ) {

    // Branches for images
    sprintf(hname, "EB_adc%d",iS);
    tree->Branch(hname, &vEB_adc_[iS]);

    // Histograms for monitoring
    sprintf(htitle,"ADC(iphi, ieta);iphi;ieta");
    hEB_adc[iS] = fs->make<TProfile2D>(hname, htitle,
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );

  } // iS

  hEB_simEnergy = fs->make<TProfile2D>("EB_simEnergy", "E(#iphi,#ieta);#iphi;#ieta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );

} // branchesEB()

// Fill EB rechits _________________________________________________________________//
void SCRegressor::fillEB ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  float energy_;

  vEB_simEnergy_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vEB_energy_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vEB_time_.assign( EBDetId::kSizeForDenseIndexing, 0. );

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByToken( EBRecHitCollectionT_, EBRecHitsH_);
  nRecHits_ = EBRecHitsH_->size();

  // Fill EB rechits 
  for ( EcalRecHitCollection::const_iterator iRHit = EBRecHitsH_->begin();
        iRHit != EBRecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ <= zs ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi() - 1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    // Fill histograms for monitoring 
    hEB_energy->Fill( iphi_,ieta_,energy_ );
    hEB_time->Fill( iphi_,ieta_,iRHit->time() );
    // Get Hashed Index: provides convenient 
    // index mapping from [ieta][iphi] -> [idx]
    idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
    // Fill vectors for images
    vEB_energy_[idx_] = energy_;
    vEB_time_[idx_] = iRHit->time();

    //std::cout << "idx,ieta,iphi,E:" <<idx_<<","<< ieta_ << "," << iphi_ << "," << iRHit->energy() << std::endl;

  } // EB rechits

  /*
  edm::Handle<EBDigiCollection> EBDigisH;
  iEvent.getByToken(EBDigiCollectionT_, EBDigisH);
  nDigis_ = EBDigisH->size();

  int idx;
  // Initialize arrays
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS)
    vEB_adc_[iS].assign(EBDetId::kSizeForDenseIndexing,0);

  // Record signal-full entries
  for(EBDigiCollection::const_iterator iDigi = EBDigisH->begin();
      iDigi != EBDigisH->end();
      ++iDigi) {
    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iDigi->id() );
    //DetId id( iDigi->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    // Get Hashed Index & Cell Geometry
    // Hashed index provides a convenient index mapping
    // from [ieta][iphi] -> [idx]
    idx = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
    //cell  = caloGeom->getGeometry(ebId);
    // Unpack the digi into a dataframe
    EcalDataFrame df(*iDigi);
    for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {
      // Get the iS-th sample
      EcalMGPASample digiSample( df.sample(iS) );
      // Fill some histograms to monitor distributions
      // These will contain *cumulative* statistics and as such
      // should be used for monitoring purposes only
      hEB_adc[iS]->Fill( iphi_, ieta_, digiSample.adc() );
      // Fill event arrays
      // These are the actual inputs to the detector images
      vEB_adc_[iS][idx] += digiSample.adc();
      //vEB_adc_[iS][idx] += digiSample.adc()/TMath::CosH(cell->etaPos()); // pick out only transverse component
    } // sample

  } // EB digi

  edm::Handle<PCaloHitContainer> EBSimHitsH;
  iEvent.getByToken(EBSimHitCollectionT_, EBSimHitsH);
  nSimHits_ = EBSimHitsH->size();
  for ( PCaloHitContainer::const_iterator iRHit = EBSimHitsH->begin();
        iRHit != EBSimHitsH->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ <= zs ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi() - 1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    // Fill histograms for monitoring
    hEB_simEnergy->Fill( iphi_, ieta_, energy_ );
    // Get Hashed Index: provides convenient
    // index mapping from [ieta][iphi] -> [idx]
    idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
    // Fill vectors for images
    vEB_simEnergy_[idx_] = energy_;

  } // EB simhits
  */

} // fillEB()
