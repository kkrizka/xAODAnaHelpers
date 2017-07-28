#ifndef xAODAnaHelpers_TruthParticle_H
#define xAODAnaHelpers_TruthParticle_H

#include "xAODAnaHelpers/Particle.h"

namespace xAH {

  class TruthParticle : public Particle
  {
    ClassDef(TruthParticle, 1);

  public:    
    TruthParticle() : Particle() {};
    virtual ~TruthParticle() {};

    // all
    int pdgId;
    int status;
    int barcode;

    // type
    bool is_higgs;
    bool is_bhad;

    // BVtx
    float Bdecay_x;
    float Bdecay_y;
    float Bdecay_z;
      
    // Parents
    int nParents;
    std::vector<int> parent_pdgId;
    std::vector<int> parent_barcode;
    std::vector<int> parent_status;

    // Children
    int nChildren;
    std::vector<int> child_pdgId;
    std::vector<int> child_barcode;
    std::vector<int> child_status;
  };
}//xAH
#endif // xAODAnaHelpers_TruthParticle_H
