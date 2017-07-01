#include "xAODAnaHelpers/PhotonHelpTree.h"

#include <iostream>

using namespace xAH;

PhotonHelpTree::PhotonHelpTree(const std::string& name, const std::string& detailStr, float units, bool mc)
  : ParticleHelpTree("xAH::Photon",name, detailStr, units, mc)
{ }

PhotonHelpTree::~PhotonHelpTree()
{ }

void PhotonHelpTree::createBranches(TTree *tree)
{
  ParticleHelpTree::createBranches(tree);

  if(m_infoSwitch.m_isolation)
    {
      setBranchStatus(tree, "isIsolated_FixedCutTightCaloOnly", 1);
      setBranchStatus(tree, "isIsolated_FixedCutTight",         1);
      setBranchStatus(tree, "isIsolated_FixedCutLoose",         1);
      setBranchStatus(tree, "ptcone20",                         1);
      setBranchStatus(tree, "ptcone30",                         1);
      setBranchStatus(tree, "ptcone40",                         1);
      setBranchStatus(tree, "ptvarcone20",                      1);
      setBranchStatus(tree, "ptvarcone30",                      1);
      setBranchStatus(tree, "ptvarcone40",                      1);
      setBranchStatus(tree, "topoetcone20",                     1);
      setBranchStatus(tree, "topoetcone30",                     1);
      setBranchStatus(tree, "topoetcone40",                     1);
    }

  // PID
  if(m_infoSwitch.m_PID)
    {
      tree->Branch(("n"+m_name+"_IsLoose").c_str(), &m_n_IsLoose);
      setBranchStatus(tree,  "PhotonID_Loose"  , 1 );

      tree->Branch(("n"+m_name+"_IsMedium").c_str(), &m_n_IsMedium);
      setBranchStatus(tree,  "PhotonID_Medium" , 1 );

      tree->Branch(("n"+m_name+"_IsTight").c_str(), &m_n_IsTight);
      setBranchStatus(tree,  "PhotonID_Tight"  , 1 );
  }

  
  if(m_infoSwitch.m_purity)
    {
      setBranchStatus(tree,"Rhad1"  , 1);
      setBranchStatus(tree,"Rhad"   , 1);
      setBranchStatus(tree,"e277"   , 1);
      setBranchStatus(tree,"Reta"   , 1);
      setBranchStatus(tree,"Rphi"   , 1);
      setBranchStatus(tree,"weta2"  , 1);
      setBranchStatus(tree,"f1"     , 1);
      setBranchStatus(tree,"wtots1" , 1);
      setBranchStatus(tree,"DeltaE" , 1);
      setBranchStatus(tree,"Eratio" , 1);
    }

  if(m_infoSwitch.m_effSF && m_mc)
    {
      setBranchStatus(tree, "PhotonID_Tight_EffSF",  1);
      setBranchStatus(tree, "PhotonID_Medium_EffSF", 1);
      setBranchStatus(tree, "PhotonID_Loose_EffSF",  1);

      setBranchStatus(tree, "PhotonID_Tight_EffSF_Error",  1);
      setBranchStatus(tree, "PhotonID_Medium_EffSF_Error", 1);
      setBranchStatus(tree, "PhotonID_Loose_EffSF_Error",  1);
    }

  if(m_infoSwitch.m_trigger)
    {
      setBranchStatus(tree, "trigMatched", 1);
    }

}

void PhotonHelpTree::clear()
{  
  ParticleHelpTree::clear();

  // PID
  if(m_infoSwitch.m_PID)
    {
      m_n_IsLoose = 0;
      m_n_IsMedium= 0;
      m_n_IsTight = 0;
    }
}


void PhotonHelpTree::fillPhoton( const xAOD::Photon* photon )
{
  ParticleHelpTree::fillParticle(photon);
  xAH::Photon* myphoton=static_cast<xAH::Photon*>(m_particles->Last());

  if ( m_infoSwitch.m_isolation ) 
    {
      static SG::AuxElement::ConstAccessor<char> isIsolated_FixedCutTightCaloOnly ("isIsolated_FixedCutTightCaloOnly");
      SAFE_SET(myphoton,isIsolated_FixedCutTightCaloOnly,photon);

      static SG::AuxElement::ConstAccessor<char> isIsolated_FixedCutTight         ("isIsolated_FixedCutTight");
      SAFE_SET(myphoton,isIsolated_FixedCutTight,photon);

      static SG::AuxElement::ConstAccessor<char> isIsolated_FixedCutLoose         ("isIsolated_FixedCutLoose");
      SAFE_SET(myphoton,isIsolated_FixedCutLoose,photon);

      myphoton->ptcone20    = photon->isolation( xAOD::Iso::ptcone20    ) / m_units;
      myphoton->ptcone30    = photon->isolation( xAOD::Iso::ptcone30    ) / m_units;
      myphoton->ptcone40    = photon->isolation( xAOD::Iso::ptcone40    ) / m_units;
      myphoton->ptvarcone20 = photon->isolation( xAOD::Iso::ptvarcone20 ) / m_units;
      myphoton->ptvarcone30 = photon->isolation( xAOD::Iso::ptvarcone30 ) / m_units;
      myphoton->ptvarcone40 = photon->isolation( xAOD::Iso::ptvarcone40 ) / m_units;
      myphoton->topoetcone20= photon->isolation( xAOD::Iso::topoetcone20) / m_units;
      myphoton->topoetcone30= photon->isolation( xAOD::Iso::topoetcone30) / m_units;
      myphoton->topoetcone40= photon->isolation( xAOD::Iso::topoetcone40) / m_units;
    }

  if ( m_infoSwitch.m_PID )
    {
      static SG::AuxElement::ConstAccessor<bool> PhotonID_Loose  ("PhotonID_Loose");
      SAFE_SET(myphoton,PhotonID_Loose,photon);
      if(PhotonID_Loose.isAvailable(*photon) && PhotonID_Loose(*photon)) m_n_IsLoose++;

      static SG::AuxElement::ConstAccessor<bool> PhotonID_Medium ("PhotonID_Medium");
      SAFE_SET(myphoton,PhotonID_Medium,photon);
      if(PhotonID_Medium.isAvailable(*photon) && PhotonID_Medium(*photon)) m_n_IsMedium++;

      static SG::AuxElement::ConstAccessor<bool> PhotonID_Tight  ("PhotonID_Tight");
      SAFE_SET(myphoton,PhotonID_Tight,photon);
      if(PhotonID_Tight.isAvailable(*photon) && PhotonID_Tight(*photon)) m_n_IsTight++;
    }

  if (m_infoSwitch.m_purity) 
    {
      static SG::AuxElement::ConstAccessor<float> Rhad1  ("Rhad1"  );
      static SG::AuxElement::ConstAccessor<float> Rhad   ("Rhad"   );
      static SG::AuxElement::ConstAccessor<float> e277   ("e277"   );
      static SG::AuxElement::ConstAccessor<float> Reta   ("Reta"   );
      static SG::AuxElement::ConstAccessor<float> Rphi   ("Rphi"   );
      static SG::AuxElement::ConstAccessor<float> weta2  ("weta2"  );
      static SG::AuxElement::ConstAccessor<float> f1     ("f1"     );
      static SG::AuxElement::ConstAccessor<float> wtots1 ("wtots1" );
      static SG::AuxElement::ConstAccessor<float> w1     ("w1"     );
      static SG::AuxElement::ConstAccessor<float> DeltaE ("DeltaE" );
      static SG::AuxElement::ConstAccessor<float> Eratio ("Eratio" );

      SAFE_SET(myphoton,Rhad1  ,photon);
      SAFE_SET(myphoton,Rhad   ,photon);
      SAFE_SET(myphoton,e277   ,photon);
      SAFE_SET(myphoton,Reta   ,photon);
      SAFE_SET(myphoton,Rphi   ,photon);
      SAFE_SET(myphoton,weta2  ,photon);
      SAFE_SET(myphoton,f1     ,photon);
      SAFE_SET(myphoton,wtots1 ,photon);
      SAFE_SET(myphoton,DeltaE ,photon);
      SAFE_SET(myphoton,Eratio ,photon);
    }

  if (m_infoSwitch.m_effSF && m_mc)
    {
      static SG::AuxElement::ConstAccessor<float> PhotonID_Tight_EffSF  ("PhotonID_Tight_EffSF"  );
      static SG::AuxElement::ConstAccessor<float> PhotonID_Medium_EffSF ("PhotonID_Medium_EffSF" );
      static SG::AuxElement::ConstAccessor<float> PhotonID_Loose_EffSF  ("PhotonID_Loose_EffSF"  );

      static SG::AuxElement::ConstAccessor<float> PhotonID_Tight_EffSF_Error  ("PhotonID_Tight_EffSF_Error" );
      static SG::AuxElement::ConstAccessor<float> PhotonID_Medium_EffSF_Error ("PhotonID_Medium_EffSF_Error" );
      static SG::AuxElement::ConstAccessor<float> PhotonID_Loose_EffSF_Error  ("PhotonID_Loose_EffSF_Error" );


      SAFE_SET(myphoton,PhotonID_Tight_EffSF ,photon);
      SAFE_SET(myphoton,PhotonID_Medium_EffSF,photon);
      SAFE_SET(myphoton,PhotonID_Loose_EffSF ,photon);

      SAFE_SET(myphoton,PhotonID_Tight_EffSF_Error ,photon);
      SAFE_SET(myphoton,PhotonID_Medium_EffSF_Error,photon);
      SAFE_SET(myphoton,PhotonID_Loose_EffSF_Error ,photon);
  }

  if (m_infoSwitch.m_trigger) 
    {
      static SG::AuxElement::ConstAccessor< std::vector< std::string> > trigMatched("trigMatched");
      SAFE_SET(myphoton,trigMatched,photon);
    }

  return;
}
