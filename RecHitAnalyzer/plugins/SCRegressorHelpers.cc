#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressorHelpers.h"

std::tuple<float, float>
SCRegressorHelpers::get_xy_at_given_z_from_eta_phi(const float & z, const float & eta, const float & phi) {
  assert(z > 0.);
  assert(eta > 0.);
  float theta = 2.0*std::atan(std::exp(-eta));
  // tan(theta) = r_cyl / z
  float r_cyl = z * std::tan(theta);
  assert(r_cyl > 0.);
  float x = r_cyl * std::cos(phi);
  float y = r_cyl * std::sin(phi);
  return std::make_tuple(x, y);
}

bool
SCRegressorHelpers::passesTriggerInspiredPreselection(const reco::PhotonRef & photon) {
  if (std::fabs(photon->pt()) <= 25.) return false;
  if (std::fabs(photon->eta()) >= 2.65) return false;

  if ( photon->full5x5_r9() <= 0.8 ) return false;
  if ( photon->hadTowOverEm() >= 0.1 ) return false;
  if ( photon->hasPixelSeed() ) return false;
  ///*
  //if ( photon->passElectronVeto() == true ) return false;
  //if ( photon->userFloat("phoChargedIsolation")/std::abs(photon->pt()) > 0.3 ) return false;

  ///*
  if ( photon->full5x5_r9() <= 0.9) {
    if ( photon->full5x5_sigmaIetaIeta() >= 0.035 ) return false;
    // if ( photon->userFloat("phoPhotonIsolation") >= 6.0 ) return false;
    if ( photon->ecalPFClusterIso() >= 6.0 ) return false;
    //if ( photon->photonIso() >= 4.0 ) return false;
    if ( photon->trkSumPtHollowConeDR03() >= 6.0 ) return false;
    //if ( photon->trackIso() >= 6. ) return false;
  }
  return true;
}
