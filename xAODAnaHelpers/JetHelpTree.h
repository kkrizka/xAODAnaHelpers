#ifndef xAODAnaHelpers_JetHelpTree_H
#define xAODAnaHelpers_JetHelpTree_H

#include <TTree.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODAnaHelpers/Jet.h>
#include <xAODAnaHelpers/ParticleHelpTree.h>

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"


namespace xAH {

  class btagOpPoint 
  {
  public:
    std::string m_name;

    bool m_mc;
    std::string m_accessorName;
    SG::AuxElement::ConstAccessor< char > m_isTag;
    SG::AuxElement::ConstAccessor< std::vector<float> > m_sf;

    int m_njets;
    std::vector<float> m_weight_sf;

    btagOpPoint(const std::string& name, bool mc, const std::string& accessorName);
    ~btagOpPoint();

    void setTree(TTree *tree, const std::string& jetName);
    void setBranch(TTree *tree, const std::string& jetName);
    void clear();
    void Fill( const xAOD::Jet* jet );
    void FillGlobalSF( const xAOD::EventInfo* eventInfo );
  
  };  //struct btagOpPoint

  class JetHelpTree : public ParticleHelpTree<Jet,HelperClasses::JetInfoSwitch>
  {
  public:
    JetHelpTree(const std::string& name = "jet", const std::string& detailStr="", float units = 1e3, bool mc = false);
    virtual ~JetHelpTree();
    
    virtual void createBranches(TTree *tree);
    virtual void clear();
    virtual void fillJet(const xAOD::Jet* jet, const xAOD::Vertex* pv, int pvLocation );
    virtual void fillGlobalBTagSF( const xAOD::EventInfo* eventInfo );

  private:
    /**
       @brief helper function to determien whether a list contains a number.
       
       Used to find if a b-tagging WP is requested.
       
       @param sfList List of b-tagging working points.
       @param workingPt Working point to search for
       @return true if workingPt is in sfList, otherwise false
    */
    bool haveBTagSF(const std::vector<int>& sfList, int workingPt);

    btagOpPoint* m_btag_Fix30;
    btagOpPoint* m_btag_Fix50;
    btagOpPoint* m_btag_Fix60;
    btagOpPoint* m_btag_Fix70;
    btagOpPoint* m_btag_Fix77;
    btagOpPoint* m_btag_Fix80;
    btagOpPoint* m_btag_Fix85;
    btagOpPoint* m_btag_Fix90;

    btagOpPoint* m_btag_Flt30;
    btagOpPoint* m_btag_Flt50;
    btagOpPoint* m_btag_Flt60;
    btagOpPoint* m_btag_Flt70;
    btagOpPoint* m_btag_Flt77;
    btagOpPoint* m_btag_Flt85;
    btagOpPoint* m_btag_Flt90;
  
    InDet::InDetTrackSelectionTool *m_trkSelTool;
  };
}



#endif // xAODAnaHelpers_JetHelpTree_H
