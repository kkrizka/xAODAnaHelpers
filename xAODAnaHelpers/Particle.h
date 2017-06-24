#ifndef xAODAnaHelpers_Particle_H
#define xAODAnaHelpers_Particle_H

#include <TLorentzVector.h>

namespace xAH {
  
  class Particle : public TObject
  {
    ClassDef(Particle, 1);
    
  public:

    Particle() : TObject() {};
    virtual ~Particle() {};

    TLorentzVector p4;
  };

}//xAH
#endif // xAODAnaHelpers_Particle_H
