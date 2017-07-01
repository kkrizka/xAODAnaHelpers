#ifndef xAODAnaHelpers_PhotonContainer_H
#define xAODAnaHelpers_PhotonContainer_H

#include <TTree.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>

#include "xAODEgamma/PhotonContainer.h"

#include <xAODAnaHelpers/HelperClasses.h>

#include <xAODAnaHelpers/Photon.h>
#include <xAODAnaHelpers/ParticleContainer.h>

typedef SG::AuxElement::Accessor< std::vector< float > > floatAccessor ;

namespace xAH {

  class PhotonContainer : public ParticleContainer<Photon,HelperClasses::PhotonInfoSwitch>
    {
    public:
      PhotonContainer(const std::string& name = "ph", const std::string& detailStr="", float units = 1e3, bool mc = false);
      virtual ~PhotonContainer();
    
      virtual void setTree(TTree *tree);
      virtual void setBranches(TTree *tree);
      virtual void clear();
      virtual void FillPhoton( const xAOD::Photon* photon );
      virtual void FillPhoton( const xAOD::IParticle* particle );
      using ParticleContainer::setTree; // make other overloaded version of execute() to show up in subclass

    protected:
      virtual void updateParticle(uint idx, Photon& photon);
    
    private:

      // isolation
      std::vector<char>*  m_isIsolated_FixedCutTightCaloOnly;
      std::vector<char>*  m_isIsolated_FixedCutTight;
      std::vector<char>*  m_isIsolated_FixedCutLoose;
      //std::vector<float>* m_etcone20;
      std::vector<float>* m_ptcone20;
      std::vector<float>* m_ptcone30;
      std::vector<float>* m_ptcone40;
      std::vector<float>* m_ptvarcone20;
      std::vector<float>* m_ptvarcone30;
      std::vector<float>* m_ptvarcone40;
      std::vector<float>* m_topoetcone20;
      std::vector<float>* m_topoetcone30;
      std::vector<float>* m_topoetcone40;
    
      // PID
      int m_n_IsLoose;
      std::vector<bool>*   m_PhotonID_Loose;
      int m_n_IsMedium;
      std::vector<bool>*   m_PhotonID_Medium;
      int m_n_IsTight;
      std::vector<bool>*   m_PhotonID_Tight;
    
      //Purity
      std::vector<float>* m_Rhad1;
      std::vector<float>* m_Rhad;
      std::vector<float>* m_e277;
      std::vector<float>* m_Reta;
      std::vector<float>* m_Rphi;
      std::vector<float>* m_weta2;
      std::vector<float>* m_f1;
      std::vector<float>* m_wtots1;
      //std::vector<float>* m_w1;
      std::vector<float>* m_DeltaE;
      std::vector<float>* m_Eratio;

      // effSF
      std::vector<float> *m_PhotonID_Tight_EffSF;
      std::vector<float> *m_PhotonID_Medium_EffSF;
      std::vector<float> *m_PhotonID_Loose_EffSF;

      std::vector<float> *m_PhotonID_Tight_EffSF_Error;
      std::vector<float> *m_PhotonID_Medium_EffSF_Error;
      std::vector<float> *m_PhotonID_Loose_EffSF_Error;

      // trigger
      std::vector<std::vector<std::string> > *m_trigMatched;
    };
}
#endif // xAODAnaHelpers_PhotonContainer_H
