#ifndef MLAnalyzer_RecHitAnalyzer_interface_SCRegressorHelpers_h
#define MLAnalyzer_RecHitAnalyzer_interface_SCRegressorHelpers_h

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include <tuple>

// Helper functions

namespace SCRegressorHelpers {
  std::tuple<float, float> get_xy_at_given_z_from_eta_phi(const float &, const float &, const float &);
  bool passesTriggerInspiredPreselection(const reco::PhotonRef &);
}

#endif
