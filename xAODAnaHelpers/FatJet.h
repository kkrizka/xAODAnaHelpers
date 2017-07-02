#ifndef xAODAnaHelpers_FatJet_H
#define xAODAnaHelpers_FatJet_H

#include <TLine.h>

#include "xAODAnaHelpers/Particle.h"
#include "xAODAnaHelpers/Jet.h"

namespace xAH {

  class FatJet : public Particle
  {
    ClassDef(FatJet, 1);

  public:

    FatJet() : Particle() {};
    virtual ~FatJet() {};

    // scale
    TLorentzVector JetConstitScaleMomentum;
    TLorentzVector JetEMScaleMomentum;

    // area
    float GhostArea;
    float ActiveArea;
    float VoronoiArea;
    float ActiveArea4vec_pt;
    float ActiveArea4vec_eta;
    float ActiveArea4vec_phi;
    float ActiveArea4vec_m;

    // substructure
    float Split12;
    float Split23;
    float Split34;
    float Tau1_wta;
    float Tau2_wta;
    float Tau3_wta;
    float Tau21_wta;
    float Tau32_wta;
    float ECF1;
    float ECF2;
    float ECF3;
    float C2;
    float D2;
    float NTrimSubjets;
    int   MyNClusters;
    int   GhostTrackCount;

    // constituent
    int   numConstituents;

    // constituentAll
    std::vector<float>  constituentWeights;
    std::vector<JetConstituent> constituents;

    // bosonCount
    int GhostTQuarksFinalCount;
    int GhostWBosonsCount;
    int GhostZBosonsCount;
    int GhostHBosonsCount;

    // VTags
    int Wtag_medium;
    int Ztag_medium;

    int Wtag_tight;
    int Ztag_tight;

    // trackJets
    std::vector<xAH::Jet> trkJets;
  };

}//xAH
#endif // xAODAnaHelpers_FatJet_H
