#ifndef xAODAnaHelpers_TrackParticle_H
#define xAODAnaHelpers_TrackParticle_H

#include "xAODAnaHelpers/Particle.h"

namespace xAH {

  class TrackParticle : public Particle
  {
    ClassDef(TrackParticle, 1);

  public:
    TrackParticle() : Particle() {};
    virtual ~TrackParticle() {};

    float chiSquared;
    float d0;

    std::vector<float> definingParametersCovMatrix;
    unsigned char expectInnermostPixelLayerHit;
    unsigned char expectNextToInnermostPixelLayerHit;

    float numberDoF;

    unsigned char numberOfInnermostPixelLayerHits;
    unsigned char numberOfInnermostPixelLayerSharedHits;
    unsigned char numberOfInnermostPixelLayerSplitHits;
    unsigned char numberOfNextToInnermostPixelLayerHits;
    unsigned char numberOfNextToInnermostPixelLayerSharedHits;
    unsigned char numberOfNextToInnermostPixelLayerSplitHits;
    unsigned char numberOfPhiHoleLayers;
    unsigned char numberOfPhiLayers;
    unsigned char numberOfPixelDeadSensors;
    unsigned char numberOfPixelHits;
    unsigned char numberOfPixelHoles;
    unsigned char numberOfPixelSharedHits;
    unsigned char numberOfPixelSplitHits;
    unsigned char numberOfPrecisionHoleLayers;
    unsigned char numberOfPrecisionLayers;
    unsigned char numberOfSCTDeadSensors;
    unsigned char numberOfSCTHits;
    unsigned char numberOfSCTHoles;
    unsigned char numberOfSCTSharedHits;
    unsigned char numberOfTRTHits;
    unsigned char numberOfTRTOutliers;

    float phi;
    float qOverP;
    float theta;

    Int_t vertexLink;
    UInt_t vertexLink_persIndex;
    UInt_t vertexLink_persKey;

    float vz;
    float z0;

  };
} //xAH

#endif // xAODAnaHelpers_TrackParticle_H

