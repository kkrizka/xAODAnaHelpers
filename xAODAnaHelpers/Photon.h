#ifndef xAODAnaHelpers_Photon_H
#define xAODAnaHelpers_Photon_H

#include "xAODAnaHelpers/Particle.h"


namespace xAH {
  
  class Photon : public Particle
  {
    ClassDef(Photon, 1);

  public:
    Photon() : Particle() {};
    virtual ~Photon() {};

    // isolation
    char isIsolated_FixedCutTightCaloOnly;
    char isIsolated_FixedCutTight;
    char isIsolated_FixedCutLoose;

    float  ptcone20;
    float  ptcone30;
    float  ptcone40;
    float  ptvarcone20;
    float  ptvarcone30;
    float  ptvarcone40;
    float  topoetcone20;
    float  topoetcone30;
    float  topoetcone40;

    // PID
    bool   PhotonID_Loose;
    bool   PhotonID_Medium;
    bool   PhotonID_Tight;

    //Purity
    float  Rhad1;
    float  Rhad;
    float  e277;
    float  Reta;
    float  Rphi;
    float  weta2;
    float  f1;
    float  wtots1;
    float  DeltaE;
    float  Eratio;

    // effSF
    float PhotonID_Tight_EffSF;
    float PhotonID_Medium_EffSF;
    float PhotonID_Loose_EffSF;

    float PhotonID_Tight_EffSF_Error;
    float PhotonID_Medium_EffSF_Error;
    float PhotonID_Loose_EffSF_Error;

    // trigger
    std::vector<std::string> trigMatched;
  };

}//xAH
#endif // xAODAnaHelpers_Photon_H
