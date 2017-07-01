#ifndef xAODAnaHelpers_PhotonHelpTree_H
#define xAODAnaHelpers_PhotonHelpTree_H

#include <TTree.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

#include "xAODEgamma/PhotonContainer.h"

#include <xAODAnaHelpers/HelperClasses.h>

#include <xAODAnaHelpers/Photon.h>
#include <xAODAnaHelpers/ParticleHelpTree.h>

namespace xAH {

  class PhotonHelpTree : public ParticleHelpTree<Photon,HelperClasses::PhotonInfoSwitch>
  {
  public:
    PhotonHelpTree(const std::string& name = "ph", const std::string& detailStr="", float units = 1e3, bool mc = false);
    virtual ~PhotonHelpTree();
    
    virtual void createBranches(TTree *tree);
    virtual void clear();
    virtual void fillPhoton( const xAOD::Photon* photon );

  private:
    // PID
    int m_n_IsLoose;
    int m_n_IsMedium;
    int m_n_IsTight;
  };
}
#endif // xAODAnaHelpers_PhotonHelpTree_H
