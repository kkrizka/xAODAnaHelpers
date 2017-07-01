#include "xAODAnaHelpers/PhotonContainer.h"

#include <iostream>

using namespace xAH;

PhotonContainer::PhotonContainer(const std::string& name, const std::string& detailStr, float units, bool mc)
  : ParticleContainer(name, detailStr, units, mc, true)
{


  if(m_infoSwitch.m_isolation){
    m_isIsolated_FixedCutTightCaloOnly = new std::vector<char>   ();
    m_isIsolated_FixedCutTight         = new std::vector<char>   ();
    m_isIsolated_FixedCutLoose         = new std::vector<char>   ();
    m_ptcone20                         = new std::vector<float> ();
    m_ptcone30                         = new std::vector<float> ();
    m_ptcone40                         = new std::vector<float> ();
    m_ptvarcone20                      = new std::vector<float> ();
    m_ptvarcone30                      = new std::vector<float> ();
    m_ptvarcone40                      = new std::vector<float> ();
    m_topoetcone20                     = new std::vector<float> ();
    m_topoetcone30                     = new std::vector<float> ();
    m_topoetcone40                     = new std::vector<float> ();
  }    

      // PID
  if(m_infoSwitch.m_PID){
    m_n_IsLoose  = 0;
    m_n_IsMedium = 0;
    m_n_IsTight  = 0;  

    m_PhotonID_Loose    = new std::vector<bool>();
    m_PhotonID_Medium   = new std::vector<bool>();
    m_PhotonID_Tight    = new std::vector<bool>();
  }

  if(m_infoSwitch.m_purity){
      //Purity
      m_Rhad1    = new std::vector<float> ();
      m_Rhad     = new std::vector<float> ();
      m_e277     = new std::vector<float> ();
      m_Reta	 = new std::vector<float> ();
      m_Rphi	 = new std::vector<float> ();
      m_weta2    = new std::vector<float> ();
      m_f1	 = new std::vector<float> ();
      m_wtots1	 = new std::vector<float> ();
      m_DeltaE   = new std::vector<float> ();
      m_Eratio   = new std::vector<float> ();
      //std::vector<float> m_w1
  }
  
  if(m_infoSwitch.m_effSF && m_mc){
    m_PhotonID_Tight_EffSF =new std::vector<float>();
    m_PhotonID_Medium_EffSF=new std::vector<float>();
    m_PhotonID_Loose_EffSF =new std::vector<float>();

    m_PhotonID_Tight_EffSF_Error =new std::vector<float>();
    m_PhotonID_Medium_EffSF_Error=new std::vector<float>();
    m_PhotonID_Loose_EffSF_Error =new std::vector<float>();
  }

  if(m_infoSwitch.m_trigger){
    m_trigMatched=new std::vector<std::vector<std::string> >();
  }
}

PhotonContainer::~PhotonContainer()
{
  if(m_infoSwitch.m_isolation){
    delete m_isIsolated_FixedCutTightCaloOnly;
    delete m_isIsolated_FixedCutTight	   ;
    delete m_isIsolated_FixedCutLoose	   ;
    delete m_ptcone20		   ;
    delete m_ptcone30		   ;
    delete m_ptcone40		   ;
    delete m_ptvarcone20	   ;
    delete m_ptvarcone30	   ;
    delete m_ptvarcone40	   ;
    delete m_topoetcone20	   ;
    delete m_topoetcone30	   ;
    delete m_topoetcone40          ;
  }    

  // PID
  if(m_infoSwitch.m_PID){
    delete m_PhotonID_Loose;
    delete m_PhotonID_Medium;
    delete m_PhotonID_Tight;
  }

  if(m_infoSwitch.m_purity){
    delete m_Rhad1;
    delete m_Rhad;
    delete m_e277;
    delete m_Reta;
    delete m_Rphi;
    delete m_weta2;
    delete m_f1;
    delete m_wtots1;
    delete m_DeltaE;
    delete m_Eratio;
    //std::vector<float> m_w1
  }

  if(m_infoSwitch.m_effSF){
    delete m_PhotonID_Tight_EffSF;
    delete m_PhotonID_Medium_EffSF;
    delete m_PhotonID_Loose_EffSF;

    delete m_PhotonID_Tight_EffSF_Error;
    delete m_PhotonID_Medium_EffSF_Error;
    delete m_PhotonID_Loose_EffSF_Error;
  }

  if(m_infoSwitch.m_trigger){
    delete m_trigMatched;
  }
}

void PhotonContainer::setTree(TTree *tree)
{
  //
  // Connect branches
  ParticleContainer::setTree(tree);

  tree->SetBranchStatus  ("nph" , 1);
  tree->SetBranchAddress ("nph" , &m_n);

  if(m_infoSwitch.m_isolation){
    connectBranch<char> (tree, "isIsolated_FixedCutTightCaloOnly", &m_isIsolated_FixedCutTightCaloOnly );
    connectBranch<char> (tree, "isIsolated_FixedCutTight",         &m_isIsolated_FixedCutTight         );
    connectBranch<char> (tree, "isIsolated_FixedCutLoose",         &m_isIsolated_FixedCutLoose         );
    connectBranch<float>(tree, "ptcone20",                  &m_ptcone20                  );
    connectBranch<float>(tree, "ptcone30",                  &m_ptcone30                  );
    connectBranch<float>(tree, "ptcone40",                  &m_ptcone40                  );
    connectBranch<float>(tree, "ptvarcone20",               &m_ptvarcone20               );
    connectBranch<float>(tree, "ptvarcone30",               &m_ptvarcone30               );
    connectBranch<float>(tree, "ptvarcone40",               &m_ptvarcone40               );
    connectBranch<float>(tree, "topoetcone20",              &m_topoetcone20              );
    connectBranch<float>(tree, "topoetcone30",              &m_topoetcone30              );
    connectBranch<float>(tree, "topoetcone40",              &m_topoetcone40              );
  }    

  // PID
  if(m_infoSwitch.m_PID){
    tree->SetBranchStatus (("n"+m_name+"_IsLoose").c_str(),     1);
    tree->SetBranchAddress(("n"+m_name+"_IsLoose").c_str(),      &m_n_IsLoose);
    connectBranch<bool>(tree,  "PhotonID_Loose"  , &m_PhotonID_Loose );

    tree->SetBranchStatus (("n"+m_name+"_IsMedium").c_str(),     1);
    tree->SetBranchAddress(("n"+m_name+"_IsMedium").c_str(),      &m_n_IsMedium);
    connectBranch<bool>(tree,  "PhotonID_Medium" , &m_PhotonID_Medium);

    tree->SetBranchStatus (("n"+m_name+"_IsTight").c_str(),     1);
    tree->SetBranchAddress(("n"+m_name+"_IsTight").c_str(),      &m_n_IsTight);
    connectBranch<bool>(tree,  "PhotonID_Tight"  , &m_PhotonID_Tight );
  }


  if(m_infoSwitch.m_purity){
    connectBranch<float>(tree,"Rhad1"  , &m_Rhad1  );
    connectBranch<float>(tree,"Rhad"   , &m_Rhad   );
    connectBranch<float>(tree,"e277"   , &m_e277   );
    connectBranch<float>(tree,"Reta"   , &m_Reta   );
    connectBranch<float>(tree,"Rphi"   , &m_Rphi   );
    connectBranch<float>(tree,"weta2"  , &m_weta2  );
    connectBranch<float>(tree,"f1"     , &m_f1     );
    connectBranch<float>(tree,"wtots1" , &m_wtots1 );
    connectBranch<float>(tree,"DeltaE" , &m_DeltaE );
    connectBranch<float>(tree,"Eratio" , &m_Eratio );
  }

  if(m_infoSwitch.m_effSF && m_mc)
    {
      connectBranch<float>(tree, "PhotonID_Tight_EffSF", &m_PhotonID_Tight_EffSF);
      connectBranch<float>(tree, "PhotonID_Medium_EffSF",&m_PhotonID_Medium_EffSF);
      connectBranch<float>(tree, "PhotonID_Loose_EffSF", &m_PhotonID_Loose_EffSF);

      connectBranch<float>(tree, "PhotonID_Tight_EffSF_Error", &m_PhotonID_Tight_EffSF_Error);
      connectBranch<float>(tree, "PhotonID_Medium_EffSF_Error",&m_PhotonID_Medium_EffSF_Error);
      connectBranch<float>(tree, "PhotonID_Loose_EffSF_Error", &m_PhotonID_Loose_EffSF_Error);
    }

  if(m_infoSwitch.m_trigger)
    {
      connectBranch<std::vector<std::string> >(tree, "trigMatched", &m_trigMatched);
    }

}

void PhotonContainer::updateParticle(uint idx, Photon& photon)
{
  ParticleContainer::updateParticle(idx,photon);

  if(m_infoSwitch.m_isolation){
    photon.isIsolated_FixedCutTightCaloOnly =  m_isIsolated_FixedCutTightCaloOnly ->at(idx);
    photon.isIsolated_FixedCutTight =          m_isIsolated_FixedCutTight         ->at(idx);
    photon.isIsolated_FixedCutLoose =          m_isIsolated_FixedCutLoose         ->at(idx);
    photon.ptcone20 =                          m_ptcone20                         ->at(idx);
    photon.ptcone30 =                          m_ptcone30                         ->at(idx);
    photon.ptcone40 =                          m_ptcone40                         ->at(idx);
    photon.ptvarcone20 =                       m_ptvarcone20                      ->at(idx);
    photon.ptvarcone30 =                       m_ptvarcone30                      ->at(idx);
    photon.ptvarcone40 =                       m_ptvarcone40                      ->at(idx);
    photon.topoetcone20 =                      m_topoetcone20                     ->at(idx);
    photon.topoetcone30 =                      m_topoetcone30                     ->at(idx);
    photon.topoetcone40 =                      m_topoetcone40                     ->at(idx);
  }    

  // PID
  if(m_infoSwitch.m_PID){
    photon.PhotonID_Loose  = m_PhotonID_Loose ->at(idx);
    photon.PhotonID_Medium = m_PhotonID_Medium->at(idx);
    photon.PhotonID_Tight  = m_PhotonID_Tight ->at(idx);
  }

  
  if(m_infoSwitch.m_purity){
    photon.Rhad1 =  m_Rhad1  ->at(idx);
    photon.Rhad =   m_Rhad   ->at(idx);
    photon.e277 =   m_e277   ->at(idx);
    photon.Reta =   m_Reta   ->at(idx);
    photon.Rphi =   m_Rphi   ->at(idx);
    photon.weta2 =  m_weta2  ->at(idx);
    photon.f1 =     m_f1     ->at(idx);
    photon.wtots1 = m_wtots1 ->at(idx);
    photon.DeltaE = m_DeltaE ->at(idx);
    photon.Eratio = m_Eratio ->at(idx);
  }

  if(m_infoSwitch.m_effSF && m_mc){
    photon.PhotonID_Loose_EffSF =m_PhotonID_Loose_EffSF ->at(idx);
    photon.PhotonID_Medium_EffSF=m_PhotonID_Medium_EffSF->at(idx);
    photon.PhotonID_Tight_EffSF =m_PhotonID_Tight_EffSF ->at(idx);

    photon.PhotonID_Loose_EffSF_Error =m_PhotonID_Loose_EffSF_Error ->at(idx);
    photon.PhotonID_Medium_EffSF_Error=m_PhotonID_Medium_EffSF_Error->at(idx);
    photon.PhotonID_Tight_EffSF_Error =m_PhotonID_Tight_EffSF_Error ->at(idx);
  }

  if(m_infoSwitch.m_trigger){
    photon.trigMatched =m_trigMatched->at(idx);
  }

}


void PhotonContainer::setBranches(TTree *tree)
{
  ParticleContainer::setBranches(tree);


  if(m_infoSwitch.m_isolation){
    setBranch<char> (tree, "isIsolated_FixedCutTightCaloOnly", m_isIsolated_FixedCutTightCaloOnly );
    setBranch<char> (tree, "isIsolated_FixedCutTight",         m_isIsolated_FixedCutTight         );
    setBranch<char> (tree, "isIsolated_FixedCutLoose",         m_isIsolated_FixedCutLoose         );
    setBranch<float>(tree, "ptcone20",                         m_ptcone20                         );
    setBranch<float>(tree, "ptcone30",                         m_ptcone30                         );
    setBranch<float>(tree, "ptcone40",                         m_ptcone40                         );
    setBranch<float>(tree, "ptvarcone20",                      m_ptvarcone20                      );
    setBranch<float>(tree, "ptvarcone30",                      m_ptvarcone30                      );
    setBranch<float>(tree, "ptvarcone40",                      m_ptvarcone40                      );
    setBranch<float>(tree, "topoetcone20",                     m_topoetcone20                     );
    setBranch<float>(tree, "topoetcone30",                     m_topoetcone30                     );
    setBranch<float>(tree, "topoetcone40",                     m_topoetcone40                     );
  }    

  // PID
  if(m_infoSwitch.m_PID){
    tree->Branch(("n"+m_name+"_IsLoose").c_str(),      &m_n_IsLoose);
    setBranch<bool>(tree,  "PhotonID_Loose"  , m_PhotonID_Loose );

    tree->Branch(("n"+m_name+"_IsMedium").c_str(),      &m_n_IsMedium);
    setBranch<bool>(tree,  "PhotonID_Medium" , m_PhotonID_Medium);

    tree->Branch(("n"+m_name+"_IsTight").c_str(),      &m_n_IsTight);
    setBranch<bool>(tree,  "PhotonID_Tight"  , m_PhotonID_Tight );
  }

  // purity
  if(m_infoSwitch.m_purity){
    setBranch<float>(tree,"Rhad1"  , m_Rhad1  );
    setBranch<float>(tree,"Rhad"   , m_Rhad   );
    setBranch<float>(tree,"e277"   , m_e277   );
    setBranch<float>(tree,"Reta"   , m_Reta   );
    setBranch<float>(tree,"Rphi"   , m_Rphi   );
    setBranch<float>(tree,"weta2"  , m_weta2  );
    setBranch<float>(tree,"f1"     , m_f1     );
    setBranch<float>(tree,"wtots1" , m_wtots1 );
    setBranch<float>(tree,"Deltae" , m_DeltaE );
    setBranch<float>(tree,"Eratio" , m_Eratio );
  }

  // effSF
  if(m_infoSwitch.m_effSF && m_mc){
    setBranch<float>(tree, "PhotonID_Loose_EffSF" , m_PhotonID_Loose_EffSF);
    setBranch<float>(tree, "PhotonID_Medium_EffSF", m_PhotonID_Medium_EffSF);
    setBranch<float>(tree, "PhotonID_Tight_EffSF" , m_PhotonID_Tight_EffSF);

    setBranch<float>(tree, "PhotonID_Tight_EffSF_Error" , m_PhotonID_Tight_EffSF_Error);
    setBranch<float>(tree, "PhotonID_Medium_EffSF_Error", m_PhotonID_Medium_EffSF_Error);
    setBranch<float>(tree, "PhotonID_Loose_EffSF_Error" , m_PhotonID_Loose_EffSF_Error);
  }

  // trigger
  if(m_infoSwitch.m_trigger){
    setBranch<std::vector<std::string> >(tree, "trigMatched", m_trigMatched);
  }

  return;
}



void PhotonContainer::clear()
{
  
  ParticleContainer::clear();

  if(m_infoSwitch.m_isolation){
    m_isIsolated_FixedCutTightCaloOnly-> clear();
    m_isIsolated_FixedCutTight        -> clear();
    m_isIsolated_FixedCutLoose	      -> clear();
    m_ptcone20		              -> clear();
    m_ptcone30		              -> clear();
    m_ptcone40		              -> clear();
    m_ptvarcone20	              -> clear();
    m_ptvarcone30	              -> clear();
    m_ptvarcone40	              -> clear();
    m_topoetcone20	              -> clear();
    m_topoetcone30	              -> clear();
    m_topoetcone40                    -> clear();
  }

  // PID
  if(m_infoSwitch.m_PID){
    m_n_IsLoose = 0;
    m_PhotonID_Loose -> clear();

    m_n_IsMedium = 0;
    m_PhotonID_Medium-> clear();

    m_n_IsTight = 0;
    m_PhotonID_Tight -> clear();
  }

  // purity
  if(m_infoSwitch.m_purity){
    m_Rhad1  -> clear();
    m_Rhad   -> clear();
    m_e277   -> clear();
    m_Reta   -> clear();
    m_Rphi   -> clear();
    m_weta2  -> clear();
    m_f1     -> clear();
    m_wtots1 -> clear();
    m_DeltaE -> clear();
    m_Eratio -> clear();
    //std::vector<float> m_w1
  }

  // effSF
  if(m_infoSwitch.m_effSF && m_mc){
    m_PhotonID_Tight_EffSF ->clear();
    m_PhotonID_Medium_EffSF->clear();
    m_PhotonID_Loose_EffSF ->clear();

    m_PhotonID_Tight_EffSF_Error ->clear();
    m_PhotonID_Medium_EffSF_Error->clear();
    m_PhotonID_Loose_EffSF_Error ->clear();
  }

  // trigger
  if(m_infoSwitch.m_trigger){
    m_trigMatched->clear();
  }

}


void PhotonContainer::FillPhoton( const xAOD::Photon* photon ){
  return FillPhoton(static_cast<const xAOD::IParticle*>(photon));
}

void PhotonContainer::FillPhoton( const xAOD::IParticle* particle ) 
{

  ParticleContainer::FillParticle(particle);

  const xAOD::Photon* photon=dynamic_cast<const xAOD::Photon*>(particle);


  if ( m_infoSwitch.m_isolation ) {
    
    static SG::AuxElement::ConstAccessor<char> isIsolated_FixedCutTightCaloOnly    ("isIsolated_FixedCutTightCaloOnly");
    safeFill<char, char, xAOD::Photon>(photon, isIsolated_FixedCutTightCaloOnly, m_isIsolated_FixedCutTightCaloOnly, -1);
    
    static SG::AuxElement::ConstAccessor<char> isIsolated_FixedCutTight            ("isIsolated_FixedCutTight");
    safeFill<char, char, xAOD::Photon>(photon, isIsolated_FixedCutTight, m_isIsolated_FixedCutTight, -1);

    static SG::AuxElement::ConstAccessor<char> isIsolated_FixedCutLoose            ("isIsolated_FixedCutLoose");
    safeFill<char, char, xAOD::Photon>(photon, isIsolated_FixedCutLoose, m_isIsolated_FixedCutLoose, -1);

    m_ptcone20     -> push_back( photon->isolation( xAOD::Iso::ptcone20    ) / m_units  );
    m_ptcone30     -> push_back( photon->isolation( xAOD::Iso::ptcone30    ) / m_units  );
    m_ptcone40     -> push_back( photon->isolation( xAOD::Iso::ptcone40    ) / m_units  );
    m_ptvarcone20  -> push_back( photon->isolation( xAOD::Iso::ptvarcone20 ) / m_units  );
    m_ptvarcone30  -> push_back( photon->isolation( xAOD::Iso::ptvarcone30 ) / m_units  );
    m_ptvarcone40  -> push_back( photon->isolation( xAOD::Iso::ptvarcone40 ) / m_units  );
    m_topoetcone20 -> push_back( photon->isolation( xAOD::Iso::topoetcone20) / m_units  );
    m_topoetcone30 -> push_back( photon->isolation( xAOD::Iso::topoetcone30) / m_units  );
    m_topoetcone40 -> push_back( photon->isolation( xAOD::Iso::topoetcone40) / m_units  );
    
  }

  if ( m_infoSwitch.m_PID ) {
  
    static SG::AuxElement::ConstAccessor<bool> PhotonID_Loose  ("PhotonID_Loose");
    safeFill<bool, bool, xAOD::Photon>(photon, PhotonID_Loose, m_PhotonID_Loose, -1);

    static SG::AuxElement::ConstAccessor<bool> PhotonID_Medium ("PhotonID_Medium");
    safeFill<bool, bool, xAOD::Photon>(photon, PhotonID_Medium, m_PhotonID_Medium, -1);

    static SG::AuxElement::ConstAccessor<bool> PhotonID_Tight  ("PhotonID_Tight");
    safeFill<bool, bool, xAOD::Photon>(photon, PhotonID_Tight, m_PhotonID_Tight, -1);

  }

  if (m_infoSwitch.m_purity) {
    static SG::AuxElement::ConstAccessor<float> Rhad1  ("Rhad1"  );
    static SG::AuxElement::ConstAccessor<float> Rhad   ("Rhad"   );
    static SG::AuxElement::ConstAccessor<float> e277   ("e277"   );
    static SG::AuxElement::ConstAccessor<float> Reta   ("Reta"   );
    static SG::AuxElement::ConstAccessor<float> Rphi   ("Rphi"   );
    static SG::AuxElement::ConstAccessor<float> weta2  ("weta2"  );
    static SG::AuxElement::ConstAccessor<float> f1     ("f1"     );
    static SG::AuxElement::ConstAccessor<float> wtots1 ("wtots1" );
    //static SG::AuxElement::ConstAccessor<float> w1       ("w1"     );
    static SG::AuxElement::ConstAccessor<float> DeltaE ("DeltaE" );
    static SG::AuxElement::ConstAccessor<float> Eratio ("Eratio" );
    
    m_Rhad1  -> push_back( Rhad1  (*photon) );
    m_Rhad   -> push_back( Rhad   (*photon) );
    m_e277   -> push_back( e277   (*photon) );
    m_Reta   -> push_back( Reta   (*photon) );
    m_Rphi   -> push_back( Rphi   (*photon) );
    m_weta2  -> push_back( weta2  (*photon) );
    m_f1     -> push_back( f1     (*photon) );
    m_wtots1 -> push_back( wtots1 (*photon) );
    m_DeltaE -> push_back( DeltaE (*photon) );
    m_Eratio -> push_back( Eratio (*photon) );
  }

  if (m_infoSwitch.m_effSF && m_mc) {
    static SG::AuxElement::ConstAccessor<float> PhotonID_Tight_EffSF  ("PhotonID_Tight_EffSF"  );
    static SG::AuxElement::ConstAccessor<float> PhotonID_Medium_EffSF ("PhotonID_Medium_EffSF" );
    static SG::AuxElement::ConstAccessor<float> PhotonID_Loose_EffSF  ("PhotonID_Loose_EffSF"  );

    static SG::AuxElement::ConstAccessor<float> PhotonID_Tight_EffSF_Error  ("PhotonID_Tight_EffSF_Error" );
    static SG::AuxElement::ConstAccessor<float> PhotonID_Medium_EffSF_Error ("PhotonID_Medium_EffSF_Error" );
    static SG::AuxElement::ConstAccessor<float> PhotonID_Loose_EffSF_Error  ("PhotonID_Loose_EffSF_Error" );


    m_PhotonID_Tight_EffSF  ->push_back( PhotonID_Tight_EffSF (*photon) );
    m_PhotonID_Medium_EffSF ->push_back( PhotonID_Medium_EffSF(*photon) );
    m_PhotonID_Loose_EffSF  ->push_back( PhotonID_Loose_EffSF (*photon) );

    m_PhotonID_Tight_EffSF_Error  ->push_back( PhotonID_Tight_EffSF_Error (*photon) );
    m_PhotonID_Medium_EffSF_Error ->push_back( PhotonID_Medium_EffSF_Error(*photon) );
    m_PhotonID_Loose_EffSF_Error  ->push_back( PhotonID_Loose_EffSF_Error (*photon) );
  }

  if (m_infoSwitch.m_trigger) {
    static SG::AuxElement::ConstAccessor< std::vector< std::string> > trigMatched("trigMatched");

    m_trigMatched ->push_back( trigMatched(*photon) );
  }

  return;
}
