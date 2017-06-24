#ifndef xAODAnaHelpers_JetConstituent_H
#define xAODAnaHelpers_JetConstituent_H

#include "xAODAnaHelpers/Particle.h"
#include "xAODAnaHelpers/TrackParticle.h"
#include "xAODAnaHelpers/MuonContainer.h"


namespace xAH {

  class JetConstituent : public Particle
  {
    ClassDef(JetConstituent, 1);

  public:      
    JetConstituent() : Particle() {};
    virtual ~JetConstituent() {};
  };

} //xAH
#endif // xAODAnaHelpers_JetConstituent_H
