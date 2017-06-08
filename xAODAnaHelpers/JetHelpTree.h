#ifndef xAODAnaHelpers_JetHelpTree_H
#define xAODAnaHelpers_JetHelpTree_H

#include <TTree.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

#include <xAODJet/JetContainer.h>

#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODAnaHelpers/Jet.h>
#include <xAODAnaHelpers/ParticleHelpTree.h>

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"


namespace xAH {

  class btagOpPoint {
  public:
    std::string m_name;
    bool m_mc;
    std::string m_acessorName;
    std::string m_tagger;
    int m_njets;
    std::vector<int>*                  m_isTag;
    std::vector<float>                 m_weight_sf;
    std::vector< std::vector<float> >* m_sf;

    btagOpPoint(std::string name, bool mc, std::string acessorName, std::string tagger="mv2c10"): m_name(name), m_mc(mc), m_acessorName(acessorName), m_tagger(tagger) {
      m_isTag = new std::vector<int>();
      m_sf    = new std::vector< std::vector<float> >();
    }

    ~btagOpPoint(){
      delete m_isTag;
      delete m_sf;
    }

    void setTree(TTree *tree, std::string jetName){
      //tree->SetBranchStatus  (("n"+jetName+"s_"+m_tagger+"_"+m_name).c_str(), 1);
      //tree->SetBranchAddress (("n"+jetName+"s_"+m_tagger+"_"+m_name).c_str(), &m_njets);
      tree->SetBranchStatus  (("n"+jetName+"s_"+m_name).c_str(), 1);
      tree->SetBranchAddress (("n"+jetName+"s_"+m_name).c_str(), &m_njets);

      HelperFunctions::connectBranch<int>     (jetName, tree,"is"+m_name,      &m_isTag);
      if(m_mc) HelperFunctions::connectBranch<std::vector<float> >(jetName, tree,"SF"+m_name,       &m_sf);
    }

    void setBranch(TTree *tree, std::string jetName){
      tree->Branch(("n"+jetName+"s_"+m_name).c_str(), &m_njets, ("n"+jetName+"s_"+m_name+"/I").c_str());
      tree->Branch((jetName+"_is"+m_name).c_str(),        &m_isTag);

      if ( m_mc ) {
	tree->Branch((jetName+"_SF"+m_name).c_str(),        &m_sf);
	tree->Branch(("weight_"+jetName+"SF"+m_name).c_str(), &m_weight_sf);
      }
    }


    void clear(){
      m_njets = 0;
      m_isTag->clear();
      m_weight_sf.clear();
      m_sf->clear();
    }

    void Fill( const xAOD::Jet* jet ) {

      SG::AuxElement::ConstAccessor< char > isTag("BTag_"+m_acessorName);
      if( isTag.isAvailable( *jet ) ) {
	if ( isTag( *jet ) == 1 ) ++m_njets;
	m_isTag->push_back( isTag( *jet ) );
      } else { 
	m_isTag->push_back( -1 ); 
      }
    
      if(!m_mc) { return; }
      SG::AuxElement::ConstAccessor< std::vector<float> > sf("BTag_SF_"+m_acessorName);
      if ( sf.isAvailable( *jet ) ) {
	m_sf->push_back( sf( *jet ) );
      } else {
	std::vector<float> junk(1,-999);
	m_sf->push_back(junk);
      }

      return;
    } // Fill

    void FillGlobalSF( const xAOD::EventInfo* eventInfo ) {
      SG::AuxElement::ConstAccessor< std::vector<float> > sf_GLOBAL("BTag_SF_"+m_acessorName+"_GLOBAL");
      if ( sf_GLOBAL.isAvailable( *eventInfo ) ) { 
	m_weight_sf = sf_GLOBAL( *eventInfo ); 
      } else { 
	m_weight_sf.push_back(-999.0); 
      }

      return;
    }
  
  };  //struct btagOpPoint

  class JetHelpTree : public ParticleHelpTree<Jet,HelperClasses::JetInfoSwitch>
  {
  public:
    JetHelpTree(const std::string& name = "jet", const std::string& detailStr="", float units = 1e3, bool mc = false);
    virtual ~JetHelpTree();
    
    virtual void createBranches(TTree *tree);
    virtual void FillJet( const xAOD::Jet* jet,            const xAOD::Vertex* pv, int pvLocation );
    virtual void FillJet( const xAOD::IParticle* particle, const xAOD::Vertex* pv, int pvLocation );
    virtual void FillGlobalBTagSF( const xAOD::EventInfo* eventInfo );

    virtual void updateParticle(uint idx, Jet& jet);

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
