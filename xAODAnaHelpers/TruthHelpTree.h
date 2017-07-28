#ifndef xAODAnaHelpers_TruthHelpTree_H
#define xAODAnaHelpers_TruthHelpTree_H

#include <TTree.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

#include "xAODTruth/TruthParticle.h"

#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODAnaHelpers/TruthParticle.h>
#include <xAODAnaHelpers/ParticleHelpTree.h>

namespace xAH {

  class TruthHelpTree : public ParticleHelpTree<TruthParticle,HelperClasses::TruthInfoSwitch>
  {
  public:
    TruthHelpTree(const std::string& name = "truth", const std::string& detailStr="", float units = 1e3);
    virtual ~TruthHelpTree();
    
    virtual void createBranches(TTree *tree);
    virtual void fillTruth( const xAOD::TruthParticle* truth );
  };
}



#endif // xAODAnaHelpers_TruthHelpTree_H
