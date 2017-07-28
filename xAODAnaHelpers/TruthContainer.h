#ifndef xAODAnaHelpers_TruthContainer_H
#define xAODAnaHelpers_TruthContainer_H

#include <TTree.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

#include "xAODTruth/TruthParticle.h"

#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODAnaHelpers/TruthParticle.h>
#include <xAODAnaHelpers/ParticleContainer.h>

namespace xAH {

    class TruthContainer : public ParticleContainer<TruthParticle,HelperClasses::TruthInfoSwitch>
    {
    public:
      TruthContainer(const std::string& name = "truth", const std::string& detailStr="", float units = 1e3);
      virtual ~TruthContainer();
    
      virtual void setTree    (TTree *tree);
      virtual void setBranches(TTree *tree);
      virtual void clear();
      virtual void FillTruth( const xAOD::TruthParticle* truth );
      virtual void FillTruth( const xAOD::IParticle* particle );
      using ParticleContainer::setTree; // make other overloaded version of execute() to show up in subclass

    protected:

      virtual void updateParticle(uint idx, TruthParticle& truth);

    private:


      //
      // Vector branches
      //
      
      // All
      std::vector<int>* m_pdgId;
      std::vector<int>* m_status;
      std::vector<int>* m_barcode;
    
      // type
      std::vector<int>* m_is_higgs;
      std::vector<int>* m_is_bhad;

      // bVtx
      std::vector<float>* m_Bdecay_x;
      std::vector<float>* m_Bdecay_y;
      std::vector<float>* m_Bdecay_z;

      // parents
      std::vector<int>* m_nParents;
      std::vector< std::vector<int> >* m_parent_pdgId;
      std::vector< std::vector<int> >* m_parent_barcode;
      std::vector< std::vector<int> >* m_parent_status;

      // children
      std::vector<int>* m_nChildren;
      std::vector< std::vector<int> >* m_child_pdgId;
      std::vector< std::vector<int> >* m_child_barcode;
      std::vector< std::vector<int> >* m_child_status;

    };
}



#endif // xAODAnaHelpers_TruthContainer_H
