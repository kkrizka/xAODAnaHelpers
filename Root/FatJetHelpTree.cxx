#include "xAODAnaHelpers/FatJetHelpTree.h"
#include <xAODAnaHelpers/HelperFunctions.h>
#include <iostream>
#include "xAODTruth/TruthEventContainer.h"

using namespace xAH;

FatJetHelpTree::FatJetHelpTree(const std::string& name, const std::string& detailStr, float units, bool mc)
  : ParticleHelpTree("xAH::FatJet",name,detailStr,units,mc,false),
    m_trackJetPtCut(10e3),
    m_trackJetEtaCut(2.5)
{
  if ( m_infoSwitch.m_VTags ) 
    {
      m_WbosonTaggerMedium = new JetSubStructureUtils::BosonTag("medium", "smooth", "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_Wtagging_MC15_Prerecommendations_20150809.dat", false, false);
      m_ZbosonTaggerMedium = new JetSubStructureUtils::BosonTag("medium", "smooth", "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_Ztagging_MC15_Prerecommendations_20150809.dat", false, false);

      m_WbosonTaggerTight = new JetSubStructureUtils::BosonTag("tight", "smooth", "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_Wtagging_MC15_Prerecommendations_20150809.dat", false, false);
      m_ZbosonTaggerTight = new JetSubStructureUtils::BosonTag("tight", "smooth", "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_Ztagging_MC15_Prerecommendations_20150809.dat", false, false);
    }

  if( m_infoSwitch.m_trackJets )
    {
      std::string trkJetName = name;
      trkJetName += "_"+m_infoSwitch.m_trackJetName;
      m_trkJets = new xAH::JetHelpTree(trkJetName, "kinematic flavorTag constituent sfFTagFix77", m_units, m_mc);
  }
}

FatJetHelpTree::~FatJetHelpTree()
{ }

void FatJetHelpTree::createBranches(TTree *tree)
{
  //
  // Connect branches
  ParticleHelpTree::createBranches(tree);

  if( m_infoSwitch.m_scales )
    {
      setBranchStatus(tree, "JetConstitScaleMomentum", 1);
      setBranchStatus(tree, "JetEMScaleMomentum"     , 1);
    }

  if ( m_infoSwitch.m_area ) 
    {
      setBranchStatus(tree, "GhostArea"         , 1);
      setBranchStatus(tree, "ActiveArea"        , 1);
      setBranchStatus(tree, "VoronoiArea"       , 1);

      setBranchStatus(tree, "ActiveArea4vec_pt" , 1);
      setBranchStatus(tree, "ActiveArea4vec_eta", 1);
      setBranchStatus(tree, "ActiveArea4vec_phi", 1);
      setBranchStatus(tree, "ActiveArea4vec_m"  , 1);
  }

  if ( m_infoSwitch.m_substructure ) 
    {
      setBranchStatus(tree, "Split12",         1);
      setBranchStatus(tree, "Split23",         1);
      setBranchStatus(tree, "Split34",         1);
      setBranchStatus(tree, "Tau1_wta",        1);
      setBranchStatus(tree, "Tau2_wta",        1);
      setBranchStatus(tree, "Tau3_wta",        1);
      setBranchStatus(tree, "Tau21_wta",       1);
      setBranchStatus(tree, "Tau32_wta",       1);
      setBranchStatus(tree, "ECF1",            1);
      setBranchStatus(tree, "ECF2",            1);
      setBranchStatus(tree, "ECF3",            1);
      setBranchStatus(tree, "C2",              1);
      setBranchStatus(tree, "D2",              1);
      setBranchStatus(tree, "NTrimSubjets",    1);
      setBranchStatus(tree, "MyNClusters",     1);
      setBranchStatus(tree, "GhostTrackCount", 1);
    }

  if ( m_infoSwitch.m_constituent) 
    {
      setBranchStatus(tree, "numConstituents", 1);
    }

  if ( m_infoSwitch.m_constituentAll) 
    {
      setBranchStatus(tree, "constituentWeights", 1);
      setBranchStatus(tree, "constituents",       1);
    }

  if ( m_infoSwitch.m_bosonCount) 
    {
      setBranchStatus(tree, "GhostTQuarksFinalCount", 1);
      setBranchStatus(tree, "GhostWBosonsCount"     , 1);
      setBranchStatus(tree, "GhostZBosonsCount"     , 1);
      setBranchStatus(tree, "GhostHBosonsCount"     , 1);
    }

  if ( m_infoSwitch.m_VTags) 
    {
      setBranchStatus(tree, "Wtag_medium", 1);
      setBranchStatus(tree, "Ztag_medium", 1);

      setBranchStatus(tree, "Wtag_tight" , 1);
      setBranchStatus(tree, "Ztag_tight" , 1);
    }
  
  if( m_infoSwitch.m_trackJets )
    {
      setBranchStatus(tree, "trkJets" , 1);
    }
}

void FatJetHelpTree::clear()
{
  ParticleHelpTree::clear();
}

void FatJetHelpTree::fillFatJet( const xAOD::Jet* fatjet )
{
  ParticleHelpTree::fillParticle(fatjet);
  xAH::FatJet* myfatjet=static_cast<xAH::FatJet*>(m_particles->Last());

  if( m_infoSwitch.m_scales )
    {
      static SG::AuxElement::ConstAccessor<float> JetConstitScaleMomentum_pt ("JetConstitScaleMomentum_pt" );
      static SG::AuxElement::ConstAccessor<float> JetConstitScaleMomentum_eta("JetConstitScaleMomentum_eta");
      static SG::AuxElement::ConstAccessor<float> JetConstitScaleMomentum_phi("JetConstitScaleMomentum_phi");
      static SG::AuxElement::ConstAccessor<float> JetConstitScaleMomentum_m  ("JetConstitScaleMomentum_m"  );
      myfatjet->JetConstitScaleMomentum.SetPtEtaPhiM(JetConstitScaleMomentum_pt (*fatjet)/m_units,
						     JetConstitScaleMomentum_eta(*fatjet),
						     JetConstitScaleMomentum_phi(*fatjet),
						     JetConstitScaleMomentum_m  (*fatjet)/m_units);

      static SG::AuxElement::ConstAccessor<float> JetEMScaleMomentum_pt ("JetEMScaleMomentum_pt" );
      static SG::AuxElement::ConstAccessor<float> JetEMScaleMomentum_eta("JetEMScaleMomentum_eta");
      static SG::AuxElement::ConstAccessor<float> JetEMScaleMomentum_phi("JetEMScaleMomentum_phi");
      static SG::AuxElement::ConstAccessor<float> JetEMScaleMomentum_m  ("JetEMScaleMomentum_m"  );
      myfatjet->JetEMScaleMomentum.SetPtEtaPhiM(JetEMScaleMomentum_pt (*fatjet)/m_units,
						JetEMScaleMomentum_eta(*fatjet),
						JetEMScaleMomentum_phi(*fatjet),
						JetEMScaleMomentum_m  (*fatjet)/m_units);
    }

  if ( m_infoSwitch.m_area ) 
    {
      static SG::AuxElement::ConstAccessor<float> GhostArea("GhostArea");
      SAFE_SET(myfatjet,GhostArea,fatjet);
      static SG::AuxElement::ConstAccessor<float> ActiveArea("ActiveArea");
      SAFE_SET(myfatjet,ActiveArea,fatjet);
      static SG::AuxElement::ConstAccessor<float> VoronoiArea("VoronoiArea");
      SAFE_SET(myfatjet,VoronoiArea,fatjet);
      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_pt("ActiveArea4vec_pt");
      SAFE_SET(myfatjet,ActiveArea4vec_pt,fatjet);
      myfatjet->ActiveArea4vec_pt/=m_units;
      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_eta("ActiveArea4vec_eta");
      SAFE_SET(myfatjet,ActiveArea4vec_eta,fatjet);
      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_phi("ActiveArea4vec_phi");
      SAFE_SET(myfatjet,ActiveArea4vec_phi,fatjet);
      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_m("ActiveArea4vec_m");
      SAFE_SET(myfatjet,ActiveArea4vec_m, fatjet);
      myfatjet->ActiveArea4vec_m/=m_units;
    }

  if( m_infoSwitch.m_substructure )
    {
      static SG::AuxElement::ConstAccessor<float> Split12("Split12");
      SAFE_SET(myfatjet,Split12,fatjet);
      myfatjet->Split12/=m_units;

      static SG::AuxElement::ConstAccessor<float> Split23("Split23");
      SAFE_SET(myfatjet,Split23,fatjet);
      myfatjet->Split23/=m_units;

      static SG::AuxElement::ConstAccessor<float> Split34("Split34");
      SAFE_SET(myfatjet,Split34,fatjet);
      myfatjet->Split34/=m_units;

      static SG::AuxElement::ConstAccessor<float> Tau1_wta ("Tau1_wta");
      SAFE_SET(myfatjet,Tau1_wta,fatjet);

      static SG::AuxElement::ConstAccessor<float> Tau2_wta ("Tau2_wta");
      SAFE_SET(myfatjet,Tau2_wta,fatjet);

      static SG::AuxElement::ConstAccessor<float> Tau3_wta ("Tau3_wta");
      SAFE_SET(myfatjet,Tau3_wta,fatjet);

      static SG::AuxElement::ConstAccessor<float> Tau21_wta ("Tau21_wta");
      if(Tau21_wta.isAvailable( *fatjet ))
	myfatjet->Tau21_wta=Tau21_wta( *fatjet );
      else if ( Tau1_wta.isAvailable( *fatjet ) && Tau2_wta.isAvailable( *fatjet ) )
	myfatjet->Tau21_wta=Tau2_wta( *fatjet ) / Tau1_wta( *fatjet );

      static SG::AuxElement::ConstAccessor<float> Tau32_wta ("Tau32_wta");
      if(Tau32_wta.isAvailable( *fatjet ))
	myfatjet->Tau32_wta=Tau32_wta( *fatjet );
      else if ( Tau2_wta.isAvailable( *fatjet ) && Tau3_wta.isAvailable( *fatjet ) )
	myfatjet->Tau32_wta=Tau3_wta( *fatjet ) / Tau2_wta( *fatjet );

      static SG::AuxElement::ConstAccessor<int> MyNClusters ("MyNClusters");
      SAFE_SET(myfatjet,MyNClusters,fatjet);

      static SG::AuxElement::ConstAccessor<float> ECF1 ("ECF1");
      SAFE_SET(myfatjet,ECF1,fatjet);
      myfatjet->ECF1/=m_units;

      static SG::AuxElement::ConstAccessor<float> ECF2 ("ECF2");
      SAFE_SET(myfatjet,ECF2,fatjet);
      myfatjet->ECF2/=m_units;

      static SG::AuxElement::ConstAccessor<float> ECF3 ("ECF3");      
      SAFE_SET(myfatjet,ECF3,fatjet);
      myfatjet->ECF3/=m_units;

      static SG::AuxElement::ConstAccessor<int> NTrimSubjets("NTrimSubjets");
      SAFE_SET(myfatjet,NTrimSubjets,fatjet);

      static SG::AuxElement::ConstAccessor<float> D2 ("D2");
      if( D2.isAvailable( *fatjet ) )
	myfatjet->D2=D2( *fatjet );
      else if (ECF1.isAvailable( *fatjet ) && ECF2.isAvailable( *fatjet ) && ECF3.isAvailable( *fatjet ))
	{
	  float e2=(ECF2( *fatjet )/(ECF1( *fatjet )*ECF1( *fatjet )));
	  float e3=(ECF3( *fatjet )/(ECF1( *fatjet )*ECF1( *fatjet )*ECF1( *fatjet )));
	  myfatjet->D2= e3/(e2*e2*e2);
	}

      static SG::AuxElement::ConstAccessor<float> C2 ("C2");
      if(C2.isAvailable(*fatjet))
	myfatjet->C2=C2(*fatjet);
      else if( ECF1.isAvailable(*fatjet) && ECF2.isAvailable(*fatjet) && ECF3.isAvailable(*fatjet))
	myfatjet->C2= ECF3(*fatjet)*ECF1(*fatjet)/pow(ECF2(*fatjet),2.0);

      myfatjet->GhostTrackCount=fatjet->auxdata<int>("GhostTrackCount");

    }

  if( m_infoSwitch.m_constituent )
    myfatjet->numConstituents=fatjet->numConstituents();

  if( m_infoSwitch.m_constituentAll )
    {
      myfatjet->constituentWeights = fatjet->getAttribute< std::vector<float> >( "constituentWeights" );

      xAOD::JetConstituentVector consVec = fatjet->getConstituents();
      if( consVec.isValid() )
	{
	  // use the example provided in
	  // http://acode-browser.usatlas.bnl.gov/lxr/source/atlas/Event/xAOD/xAODJet/xAODJet/JetConstituentVector.h
	  JetConstituent myconstit;
	  xAOD::JetConstituentVector::const_iterator constit  = consVec.begin();
	  xAOD::JetConstituentVector::const_iterator constitE = consVec.end();
	  for( ; constit != constitE; constit++)
	    {
	      myconstit.p4.SetPtEtaPhiE(constit->pt() / m_units, constit->eta(), constit->phi(), constit->e() / m_units );
	      myfatjet->constituents.push_back(myconstit);
	    }
	}
    }

  if(m_infoSwitch.m_bosonCount)
    {
      const xAOD::Jet* fatJetParentJet = 0;

      try
	{
	  auto el = fatjet->auxdata<ElementLink<xAOD::JetContainer> >("Parent");
	  if(!el.isValid())
	    Warning("executeSingle()", "Invalid link to \"Parent\" from fat-jet.");
	  else
	    fatJetParentJet = (*el);
	}
      catch(...)
	{
	  Warning("executeSingle()", "Unable to get parent jet of fat-jet for truth labeling. Trimmed jet area would be used!");
	  fatJetParentJet = fatjet;
	}

      if(m_mc)
	{
	  static SG::AuxElement::ConstAccessor< int > GhostTQuarksFinalCount("GhostTQuarksFinalCount");
	  SAFE_SET(myfatjet,GhostTQuarksFinalCount,fatJetParentJet);

	  static SG::AuxElement::ConstAccessor< int > GhostWBosonsCount("GhostWBosonsCount");
	  SAFE_SET(myfatjet,GhostWBosonsCount,fatJetParentJet);

	  static SG::AuxElement::ConstAccessor< int > GhostZBosonsCount("GhostZBosonsCount");
	  SAFE_SET(myfatjet,GhostZBosonsCount,fatJetParentJet);

	  static SG::AuxElement::ConstAccessor< int > GhostHBosonsCount("GhostHBosonsCount");
	  SAFE_SET(myfatjet,GhostHBosonsCount,fatJetParentJet);
	}
    }
  
  if(m_infoSwitch.m_VTags)
    {
      myfatjet->Wtag_medium = m_WbosonTaggerMedium->result(*fatjet); 
      myfatjet->Ztag_medium = m_ZbosonTaggerMedium->result(*fatjet);

      myfatjet->Wtag_tight  = m_WbosonTaggerTight->result(*fatjet); 
      myfatjet->Ztag_tight  = m_ZbosonTaggerTight->result(*fatjet);
    }

  //
  // Associated track jets
  //
  if( m_infoSwitch.m_trackJets )
    {
      const xAOD::Jet* fatjet_parent = 0;
    
      try
	{
	  auto el = fatjet->auxdata<ElementLink<xAOD::JetContainer> >("Parent");
	  if(!el.isValid())
	    Warning("execute()", "Invalid link to \"Parent\" from leading calo-jet");
	  else
	    fatjet_parent = (*el);
	}
      catch(...)
	{
	  Warning("executeSingle()", "Unable to get parent jet of fat-jet for tracks labeling. Trimmed jet area would be used!");
	  fatjet_parent = fatjet;
	}

      std::vector<const xAOD::Jet*> assotrkjets;
      try
	{
	  assotrkjets = fatjet_parent->getAssociatedObjects<xAOD::Jet>(m_infoSwitch.m_trackJetName);
	}
      catch (...)
	{
	  Warning("execute()", "Unable to fetch \"%s\" link from leading calo-jet", m_infoSwitch.m_trackJetName.data());
	}

      for(const auto& TrackJet : assotrkjets)
	{
	  if(!SelectTrackJet(TrackJet)) continue;
	  m_trkJets->fillJet(TrackJet,0,0);
	  xAH::Jet *trkJet=static_cast<xAH::Jet*>(m_trkJets->particle(m_trkJets->size()-1));
	  myfatjet->trkJets.push_back(*trkJet);
	}
    }
}

bool FatJetHelpTree::SelectTrackJet(const xAOD::Jet* TrackJet)
{
  if( TrackJet->pt() < m_trackJetPtCut )            return false;
  if( fabs(TrackJet->eta()) > m_trackJetEtaCut )    return false;
  if( TrackJet->numConstituents() < 2 )             return false;

  return true;
}
