#ifndef xAODAnaHelpers_FatJetContainer_H
#define xAODAnaHelpers_FatJetContainer_H

#include <TTree.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

#include <JetSubStructureUtils/BosonTag.h>

#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODAnaHelpers/FatJet.h>
#include <xAODAnaHelpers/ParticleHelpTree.h>
#include <xAODAnaHelpers/JetHelpTree.h>

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"


namespace xAH {

  class FatJetHelpTree : public ParticleHelpTree<FatJet,HelperClasses::JetInfoSwitch>
  {
  public:
    FatJetHelpTree(const std::string& name = "fatjet", const std::string& detailStr="", float units = 1e3, bool mc = false);
    virtual ~FatJetHelpTree();
    
    virtual void createBranches(TTree *tree);
    virtual void clear();
    virtual void fillFatJet( const xAOD::Jet* jet );

    float       m_trackJetPtCut;
    float       m_trackJetEtaCut;
    std::string m_trackJetName;

  private:
      
    JetSubStructureUtils::BosonTag*      m_WbosonTaggerMedium;
    JetSubStructureUtils::BosonTag*      m_ZbosonTaggerMedium;
    JetSubStructureUtils::BosonTag*      m_WbosonTaggerTight ;
    JetSubStructureUtils::BosonTag*      m_ZbosonTaggerTight ;

    xAH::JetHelpTree* m_trkJets;
    bool SelectTrackJet(const xAOD::Jet* TrackJet);
  };
}



#endif // xAODAnaHelpers_FatJetHelpTree_H
