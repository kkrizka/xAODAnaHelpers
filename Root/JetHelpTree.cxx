#include "xAODAnaHelpers/JetHelpTree.h"

#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODTruth/TruthEventContainer.h>

#include <iostream>

using namespace xAH;

btagOpPoint::btagOpPoint(const std::string& name, bool mc, const std::string& accessorName)
  : m_name(name), m_mc(mc), m_accessorName(accessorName), m_isTag("BTag_"+accessorName), m_sf("BTag_SF_"+m_accessorName)
{ }

btagOpPoint::~btagOpPoint()
{ }

void btagOpPoint::setTree(TTree *tree, const std::string& jetName)
{
  tree->SetBranchStatus  (("n"+jetName+"s_"+m_name).c_str(), 1);
  tree->SetBranchAddress (("n"+jetName+"s_"+m_name).c_str(), &m_njets);
}

void btagOpPoint::setBranch(TTree *tree, const std::string& jetName)
{
  tree->Branch(("n"+jetName+"s_"+m_name).c_str(), &m_njets, ("n"+jetName+"s_"+m_name+"/I").c_str());

  if ( m_mc ) 
    tree->Branch(("weight_"+jetName+"SF"+m_name).c_str(), &m_weight_sf);
}

void btagOpPoint::clear()
{
  m_njets = 0;
  m_weight_sf.clear();
}

void btagOpPoint::Fill( const xAOD::Jet* jet ) 
{
  if( m_isTag.isAvailable( *jet ) ) 
    if ( m_isTag( *jet ) == 1 ) ++m_njets;
  return;
}

void btagOpPoint::FillGlobalSF( const xAOD::EventInfo* eventInfo )
{
  SG::AuxElement::ConstAccessor< std::vector<float> > sf_GLOBAL("BTag_SF_"+m_accessorName+"_GLOBAL");
  if ( sf_GLOBAL.isAvailable( *eventInfo ) )
    m_weight_sf = sf_GLOBAL( *eventInfo ); 
  else
    m_weight_sf.push_back(-999.0); 
}

JetHelpTree::JetHelpTree(const std::string& name, const std::string& detailStr, float units, bool mc)
  : ParticleHelpTree("xAH::Jet",name,detailStr,units,mc),
    m_trkSelTool(nullptr)

{
  if( !m_infoSwitch.m_sfFTagFix.empty() ) 
    {
      m_btag_Fix30   = new  btagOpPoint("Fix30",m_mc, "FixedCutBEff_30");
      m_btag_Fix50   = new  btagOpPoint("Fix50",m_mc, "FixedCutBEff_50");
      m_btag_Fix60   = new  btagOpPoint("Fix60",m_mc, "FixedCutBEff_60");
      m_btag_Fix70   = new  btagOpPoint("Fix70",m_mc, "FixedCutBEff_70");
      m_btag_Fix77   = new  btagOpPoint("Fix77",m_mc, "FixedCutBEff_77");
      m_btag_Fix80   = new  btagOpPoint("Fix80",m_mc, "FixedCutBEff_80");
      m_btag_Fix85   = new  btagOpPoint("Fix85",m_mc, "FixedCutBEff_85");
      m_btag_Fix90   = new  btagOpPoint("Fix90",m_mc, "FixedCutBEff_90");
    }

  if( !m_infoSwitch.m_sfFTagFlt.empty() ) 
    {
      m_btag_Flt30   = new  btagOpPoint("Flt30",m_mc, "FlatBEff_30");
      m_btag_Flt50   = new  btagOpPoint("Flt50",m_mc, "FlatBEff_50");
      m_btag_Flt60   = new  btagOpPoint("Flt60",m_mc, "FlatBEff_60");
      m_btag_Flt70   = new  btagOpPoint("Flt70",m_mc, "FlatBEff_70");
      m_btag_Flt77   = new  btagOpPoint("Flt77",m_mc, "FlatBEff_77");
      m_btag_Flt85   = new  btagOpPoint("Flt85",m_mc, "FlatBEff_85");
      m_btag_Flt90   = new  btagOpPoint("Flt90",m_mc, "FlatBEff_90");
    }
}

JetHelpTree::~JetHelpTree()
{
  if( !m_infoSwitch.m_sfFTagFix.empty() )
    {
      delete m_btag_Fix30;
      delete m_btag_Fix50;
      delete m_btag_Fix60;
      delete m_btag_Fix70;
      delete m_btag_Fix77;
      delete m_btag_Fix80;
      delete m_btag_Fix85;
      delete m_btag_Fix90;
    }

  if( !m_infoSwitch.m_sfFTagFlt.empty() ) 
    {
      delete m_btag_Flt30;
      delete m_btag_Flt50;
      delete m_btag_Flt60;
      delete m_btag_Flt70;
      delete m_btag_Flt77;
      delete m_btag_Flt85;
      delete m_btag_Flt90;
    }
}

void JetHelpTree::createBranches(TTree *tree)
{
  ParticleHelpTree::createBranches(tree);

  if( m_infoSwitch.m_rapidity ) 
    {
      setBranchStatus(tree, "rapidity", 1);
    }

  if( m_infoSwitch.m_clean ) 
    {
      setBranchStatus(tree, "Timing",                     1);
      setBranchStatus(tree, "LArQuality",                 1);
      setBranchStatus(tree, "HECQuality",                 1);
      setBranchStatus(tree, "NegativeE",                  1);
      setBranchStatus(tree, "AverageLArQF",               1);
      setBranchStatus(tree, "BchCorrCell",                1);
      setBranchStatus(tree, "N90Constituents",            1);
      setBranchStatus(tree, "LArBadHVEnergyFrac",         1);
      setBranchStatus(tree, "LArBadHVNCell",              1);
      setBranchStatus(tree, "OotFracClusters5",           1);
      setBranchStatus(tree, "OotFracClusters10",          1);
      setBranchStatus(tree, "LeadingClusterPt",           1);
      setBranchStatus(tree, "LeadingClusterSecondLambda", 1);
      setBranchStatus(tree, "LeadingClusterCenterLambda", 1);
      setBranchStatus(tree, "LeadingClusterSecondR",      1);
      setBranchStatus(tree, "clean_passLooseBad",         1);
      setBranchStatus(tree, "clean_passLooseBadUgly",     1);
      setBranchStatus(tree, "clean_passTightBad",         1);
      setBranchStatus(tree, "clean_passTightBadUgly",     1);
  }

  if ( m_infoSwitch.m_energy ) 
    {
      setBranchStatus(tree, "HECFrac",               1);
      setBranchStatus(tree, "EMFrac",                1);
      setBranchStatus(tree, "CentroidR",             1);
      setBranchStatus(tree, "FracSamplingMax",       1);
      setBranchStatus(tree, "FracSamplingMaxIndex",  1);
      setBranchStatus(tree, "LowEtConstituentsFrac", 1);
      setBranchStatus(tree, "GhostMuonSegmentCount", 1);
      setBranchStatus(tree, "Width",                 1);
    }

  if ( m_infoSwitch.m_scales ) 
    {
      setBranchStatus(tree, "emScalePt",              1);
      setBranchStatus(tree, "constScalePt",           1);
      setBranchStatus(tree, "pileupScalePt",          1);
      setBranchStatus(tree, "originConstitScalePt",   1);
      setBranchStatus(tree, "etaJESScalePt",          1);
      setBranchStatus(tree, "gscScalePt",             1);
      setBranchStatus(tree, "insituScalePt",          1);
    }

  if ( m_infoSwitch.m_constscaleEta ) 
    {
      setBranchStatus(tree, "constScaleEta", 1);
    }

  if ( m_infoSwitch.m_layer ) 
    {
      setBranchStatus(tree, "EnergyPerSampling", 1);
    }

  if ( m_infoSwitch.m_trackAll ) 
    {
      setBranchStatus(tree, "NumTrkPt1000",     1);
      setBranchStatus(tree, "SumPtTrkPt1000",   1);
      setBranchStatus(tree, "TrackWidthPt1000", 1);
      setBranchStatus(tree, "NumTrkPt500",      1);
      setBranchStatus(tree, "SumPtTrkPt500",    1);
      setBranchStatus(tree, "TrackWidthPt500",  1);
      setBranchStatus(tree, "JVF",              1);
    }

  if ( m_infoSwitch.m_trackPV ) 
    {
      setBranchStatus(tree, "NumTrkPt1000PV",      1);
      setBranchStatus(tree, "SumPtTrkPt1000PV",    1);
      setBranchStatus(tree, "TrackWidthPt1000PV",  1);
      setBranchStatus(tree, "NumTrkPt500PV",       1);
      setBranchStatus(tree, "SumPtTrkPt500PV",     1);
      setBranchStatus(tree, "TrackWidthPt500PV",   1);
      setBranchStatus(tree, "JVFPV",               1);
      setBranchStatus(tree, "Jvt",                 1);
      setBranchStatus(tree, "JvtJvfcorr",          1);
      setBranchStatus(tree, "JvtRpt",              1);
      if ( m_mc ) 
	{
	  setBranchStatus(tree, "JetJvtEfficiency_JVTSyst_JVT_Loose", 1);
	  setBranchStatus(tree, "JetJvtEfficiency_JVTSyst_JVT_Medium",1);
	  setBranchStatus(tree, "JetJvtEfficiency_JVTSyst_JVT_Tight", 1);
	}
    }

  if ( m_infoSwitch.m_allTrack ) 
    {
      // if want to apply the selection of the PV then need to setup track selection tool
      // this applies the JVF/JVT selection cuts
      // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JvtManualRecalculation
      if( m_infoSwitch.m_allTrackPVSel ) 
	{
	  m_trkSelTool = new InDet::InDetTrackSelectionTool( "JetTrackSelection", "Loose" );
	  m_trkSelTool->initialize();
	  // to do this need to have AddJets return a status code
	  //RETURN_CHECK( "HelpTreeBase::JetTrackSelection", m_trkSelTool->initialize(), "");
	}
      setBranchStatus(tree, "GhostTrackCount",   1);
      setBranchStatus(tree, "GhostTrackPt",      1);
      setBranchStatus(tree, "GhostTrack_pt",     1);
      setBranchStatus(tree, "GhostTrack_qOverP", 1);
      setBranchStatus(tree, "GhostTrack_eta",    1);
      setBranchStatus(tree, "GhostTrack_phi",    1);
      setBranchStatus(tree, "GhostTrack_e",      1);
      setBranchStatus(tree, "GhostTrack_d0",     1);
      setBranchStatus(tree, "GhostTrack_z0",     1);
      if ( m_infoSwitch.m_allTrackDetail ) 
	{
	  setBranchStatus(tree, "GhostTrack_nPixelHits",                           1);
	  setBranchStatus(tree, "GhostTrack_nSCTHits",                             1);
	  setBranchStatus(tree, "GhostTrack_nTRTHits",                             1);
	  setBranchStatus(tree, "GhostTrack_nPixelSharedHits",                     1);
	  setBranchStatus(tree, "GhostTrack_nPixelSplitHits",                      1);
	  setBranchStatus(tree, "GhostTrack_nInnermostPixelLayerHits",             1);
	  setBranchStatus(tree, "GhostTrack_nInnermostPixelLayerSharedHits",       1);
	  setBranchStatus(tree, "GhostTrack_nInnermostPixelLayerSplitHits",        1);
	  setBranchStatus(tree, "GhostTrack_nNextToInnermostPixelLayerHits",       1);
	  setBranchStatus(tree, "GhostTrack_nNextToInnermostPixelLayerSharedHits", 1);
	  setBranchStatus(tree, "GhostTrack_nNextToInnermostPixelLayerSplitHits",  1);
    }
  }

  if ( m_infoSwitch.m_constituent )
    {
      setBranchStatus(tree, "numConstituents",    1);
    }

  if ( m_infoSwitch.m_constituentAll )
    {
      setBranchStatus(tree, "constituentWeights", 1);
      setBranchStatus(tree, "constituent_pt",     1);
      setBranchStatus(tree, "constituent_eta",    1);
      setBranchStatus(tree, "constituent_phi",    1);
      setBranchStatus(tree, "constituent_e",      1);
    }

  if( m_infoSwitch.m_flavTag  || m_infoSwitch.m_flavTagHLT  ) 
    {
      setBranchStatus(tree, "MV2c00",  1);
      setBranchStatus(tree, "MV2c10",  1);
      setBranchStatus(tree, "MV2c20",  1);
      setBranchStatus(tree, "MV2c100", 1);

      setBranchStatus(tree, "HadronConeExclTruthLabelID",  1);

      if( m_infoSwitch.m_jetFitterDetails)
	{
	  setBranchStatus(tree, "JetFitter_nVTX"          ,1);
	  setBranchStatus(tree, "JetFitter_nSingleTracks" ,1);
	  setBranchStatus(tree, "JetFitter_nTracksAtVtx"  ,1);
	  setBranchStatus(tree, "JetFitter_mass"          ,1);
	  setBranchStatus(tree, "JetFitter_energyFraction",1);
	  setBranchStatus(tree, "JetFitter_significance3d",1);
	  setBranchStatus(tree, "JetFitter_deltaeta"      ,1);
	  setBranchStatus(tree, "JetFitter_deltaphi"      ,1);
	  setBranchStatus(tree, "JetFitter_N2Tpair"       ,1);
	}

      if( m_infoSwitch.m_svDetails)
	{
	  setBranchStatus(tree, "SV0",               1);
	  setBranchStatus(tree, "SV0_NGTinSvx",      1);
	  setBranchStatus(tree, "SV0_N2Tpair",       1);
	  setBranchStatus(tree, "SV0_massvx",        1);
	  setBranchStatus(tree, "SV0_efracsvx",      1);
	  setBranchStatus(tree, "SV0_normdist",      1);

	  setBranchStatus(tree, "SV1",               1);
	  setBranchStatus(tree, "SV1IP3D",           1);
	  setBranchStatus(tree, "SV1_pu",            1);
	  setBranchStatus(tree, "SV1_pb",            1);
	  setBranchStatus(tree, "SV1_pc",            1);
	  setBranchStatus(tree, "SV1_c",             1);
	  setBranchStatus(tree, "SV1_cu",            1);
	  setBranchStatus(tree, "SV1_NGTinSvx",      1);
	  setBranchStatus(tree, "SV1_N2Tpair",       1);
	  setBranchStatus(tree, "SV1_massvx",        1);
	  setBranchStatus(tree, "SV1_efracsvx",      1);
	  setBranchStatus(tree, "SV1_normdist",      1);
	  setBranchStatus(tree, "SV1_Lxy",           1);
	  setBranchStatus(tree, "SV1_L3d",           1);
	  setBranchStatus(tree, "SV1_distmatlay",    1);
	  setBranchStatus(tree, "SV1_dR",            1);
	}

    if( m_infoSwitch.m_ipDetails)
      {
	setBranchStatus(tree, "IP2D_pu",                   1);
	setBranchStatus(tree, "IP2D_pb",                   1);
	setBranchStatus(tree, "IP2D_pc",                   1);
	setBranchStatus(tree, "IP2D",                      1);
	setBranchStatus(tree, "IP2D_c",                    1);
	setBranchStatus(tree, "IP2D_cu",                   1);
	setBranchStatus(tree, "nIP2DTracks"              , 1);
	setBranchStatus(tree, "IP2D_gradeOfTracks"       , 1);
	setBranchStatus(tree, "IP2D_flagFromV0ofTracks"  , 1);
	setBranchStatus(tree, "IP2D_valD0wrtPVofTracks"  , 1);
	setBranchStatus(tree, "IP2D_sigD0wrtPVofTracks"  , 1);
	setBranchStatus(tree, "IP2D_weightBofTracks"     , 1);
	setBranchStatus(tree, "IP2D_weightCofTracks"     , 1);
	setBranchStatus(tree, "IP2D_weightUofTracks"     , 1);

	setBranchStatus(tree, "IP3D",                      1);
	setBranchStatus(tree, "IP3D_pu",                   1);
	setBranchStatus(tree, "IP3D_pb",                   1);
	setBranchStatus(tree, "IP3D_pc",                   1);
	setBranchStatus(tree, "IP3D",                      1);
	setBranchStatus(tree, "IP3D_c",                    1);
	setBranchStatus(tree, "IP3D_cu",                   1);
	setBranchStatus(tree, "nIP3DTracks"              , 1);
	setBranchStatus(tree, "IP3D_gradeOfTracks"       , 1);
	setBranchStatus(tree, "IP3D_flagFromV0ofTracks"  , 1);
	setBranchStatus(tree, "IP3D_valD0wrtPVofTracks"  , 1);
	setBranchStatus(tree, "IP3D_sigD0wrtPVofTracks"  , 1);
	setBranchStatus(tree, "IP3D_valZ0wrtPVofTracks"  , 1);
	setBranchStatus(tree, "IP3D_sigZ0wrtPVofTracks"  , 1);
	setBranchStatus(tree, "IP3D_weightBofTracks"     , 1);
	setBranchStatus(tree, "IP3D_weightCofTracks"     , 1);
	setBranchStatus(tree, "IP3D_weightUofTracks"     , 1);
      }

    if( m_infoSwitch.m_JVC ) 
      {
	setBranchStatus(tree, "JetVertexCharge_discriminant", 1);
      }
    }

  if( m_infoSwitch.m_flavTagHLT  ) 
    {
      setBranchStatus(tree, "vtxOnlineValid"     ,1);
      setBranchStatus(tree, "vtxHadDummy"        ,1);
      setBranchStatus(tree, "bs_online_vx"       ,1);
      setBranchStatus(tree, "bs_online_vy"       ,1);
      setBranchStatus(tree, "bs_online_vz"       ,1);

      setBranchStatus(tree, "vtx_offline_x0"     ,1);
      setBranchStatus(tree, "vtx_offline_y0"     ,1);
      setBranchStatus(tree, "vtx_offline_z0"     ,1);

      setBranchStatus(tree, "vtx_online_x0"      ,1);
      setBranchStatus(tree, "vtx_online_y0"      ,1);
      setBranchStatus(tree, "vtx_online_z0"      ,1);

      setBranchStatus(tree, "vtx_online_bkg_x0"  ,1);
      setBranchStatus(tree, "vtx_online_bkg_y0"  ,1);
      setBranchStatus(tree, "vtx_online_bkg_z0"  ,1);
    }

  if( !m_infoSwitch.m_sfFTagFix.empty() ) 
    {
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 30)) m_btag_Fix30->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 50)) m_btag_Fix50->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 60)) m_btag_Fix60->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 70)) m_btag_Fix70->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 77)) m_btag_Fix77->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 80)) m_btag_Fix80->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 85)) m_btag_Fix85->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFix, 90)) m_btag_Fix90->setBranch(tree, m_name);
    }

  if( !m_infoSwitch.m_sfFTagFlt.empty() ) 
    {
      if (haveBTagSF(m_infoSwitch.m_sfFTagFlt, 30)) m_btag_Flt30->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFlt, 50)) m_btag_Flt50->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFlt, 60)) m_btag_Flt60->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFlt, 70)) m_btag_Flt70->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFlt, 77)) m_btag_Flt77->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFlt, 85)) m_btag_Flt85->setBranch(tree, m_name);
      if (haveBTagSF(m_infoSwitch.m_sfFTagFlt, 90)) m_btag_Flt90->setBranch(tree, m_name);
    }// if sfFTagFlt

  if( m_infoSwitch.m_area ) 
    {
      setBranchStatus(tree, "GhostArea",          1);
      setBranchStatus(tree, "ActiveArea",         1);
      setBranchStatus(tree, "VoronoiArea",        1);
      setBranchStatus(tree, "ActiveArea4vec_pt",  1);
      setBranchStatus(tree, "ActiveArea4vec_eta", 1);
      setBranchStatus(tree, "ActiveArea4vec_phi", 1);
      setBranchStatus(tree, "ActiveArea4vec_m",   1);
  }

  if ( m_infoSwitch.m_truth && m_mc ) 
    {
      setBranchStatus(tree, "ConeTruthLabelID",   1);
      setBranchStatus(tree, "TruthCount",         1);
      setBranchStatus(tree, "TruthLabelDeltaR_B", 1);
      setBranchStatus(tree, "TruthLabelDeltaR_C", 1);
      setBranchStatus(tree, "TruthLabelDeltaR_T", 1);
      setBranchStatus(tree, "PartonTruthLabelID", 1);
      setBranchStatus(tree, "GhostTruthAssociationFraction", 1);
      setBranchStatus(tree, "truth_p4",           1);
      setBranchStatus(tree, "truth_pdgId",        1);
      setBranchStatus(tree, "truth_partonPt",     1);
      setBranchStatus(tree, "truth_partonDR",     1);
  }

  if ( m_infoSwitch.m_truthDetails )
    {
      setBranchStatus(tree, "GhostBHadronsFinalCount",   1);
      setBranchStatus(tree, "GhostBHadronsInitialCount", 1);
      setBranchStatus(tree, "GhostBQuarksFinalCount",    1);
      setBranchStatus(tree, "GhostBHadronsFinalPt",      1);
      setBranchStatus(tree, "GhostBHadronsInitialPt",    1);
      setBranchStatus(tree, "GhostBQuarksFinalPt",       1);

      setBranchStatus(tree, "GhostCHadronsFinalCount"  , 1);
      setBranchStatus(tree, "GhostCHadronsInitialCount", 1);
      setBranchStatus(tree, "GhostCQuarksFinalCount"   , 1);
      setBranchStatus(tree, "GhostCHadronsFinalPt"     , 1);
      setBranchStatus(tree, "GhostCHadronsInitialPt"   , 1);
      setBranchStatus(tree, "GhostCQuarksFinalPt"      , 1);

      setBranchStatus(tree, "GhostTausFinalCount",       1);
      setBranchStatus(tree, "GhostTausFinalPt"   ,       1);

      setBranchStatus(tree, "truth_pdgId"   , 1);
      setBranchStatus(tree, "truth_partonPt", 1);
      setBranchStatus(tree, "truth_partonDR", 1);
    }

  if ( m_infoSwitch.m_charge ) 
    {
      setBranchStatus(tree, "charge", 1);
    }
}

void JetHelpTree::fillJet( const xAOD::Jet* jet, const xAOD::Vertex* pv, int pvLocation )
{
  ParticleHelpTree::fillParticle(jet);
  xAH::Jet *myjet=static_cast<xAH::Jet*>(m_particles->Last());

  if ( m_infoSwitch.m_rapidity )
    {
      myjet->rapidity = jet->rapidity();
    }

  if ( m_infoSwitch.m_clean )
    {
      static SG::AuxElement::ConstAccessor<float> Timing ("Timing");
      myjet->Timing=Timing(*jet);

      static SG::AuxElement::ConstAccessor<float> LArQuality ("LArQuality");
      myjet->LArQuality=LArQuality(*jet);			      

      static SG::AuxElement::ConstAccessor<float> HECQuality ("HECQuality");
      myjet->HECQuality=HECQuality(*jet);

      static SG::AuxElement::ConstAccessor<float> NegativeE ("NegativeE");
      myjet->NegativeE=NegativeE(*jet)/m_units;

      static SG::AuxElement::ConstAccessor<float> AverageLArQF ("AverageLArQF");
      myjet->AverageLArQF=AverageLArQF(*jet);

      static SG::AuxElement::ConstAccessor<float> BchCorrCell ("BchCorrCell");
      if(BchCorrCell.isAvailable( *jet )) myjet->BchCorrCell=BchCorrCell(*jet);

      static SG::AuxElement::ConstAccessor<float> N90Constituents ("N90Constituents");
      myjet->N90Constituents=N90Constituents(*jet);

      static SG::AuxElement::ConstAccessor<float> LArBadHVEnergyFrac ("LArBadHVEnergyFrac");
      if(LArBadHVEnergyFrac.isAvailable( *jet )) myjet->LArBadHVEnergyFrac=LArBadHVEnergyFrac(*jet);

      static SG::AuxElement::ConstAccessor<int> LArBadHVNCell ("LArBadHVNCell");
      if(LArBadHVNCell.isAvailable( *jet )) myjet->LArBadHVNCell=LArBadHVNCell(*jet);

      static SG::AuxElement::ConstAccessor<float> OotFracClusters5 ("OotFracClusters5");
      myjet->OotFracClusters5=OotFracClusters5(*jet);

      static SG::AuxElement::ConstAccessor<float> OotFracClusters10 ("OotFracClusters10");
      myjet->OotFracClusters10=OotFracClusters10(*jet);

      static SG::AuxElement::ConstAccessor<float> LeadingClusterPt ("LeadingClusterPt");
      myjet->LeadingClusterPt=LeadingClusterPt(*jet);

      static SG::AuxElement::ConstAccessor<float> LeadingClusterSecondLambda ("LeadingClusterSecondLambda");
      myjet->LeadingClusterSecondLambda=LeadingClusterSecondLambda(*jet);

      static SG::AuxElement::ConstAccessor<float> LeadingClusterCenterLambda ("LeadingClusterCenterLambda");
      myjet->LeadingClusterCenterLambda=LeadingClusterCenterLambda(*jet);

      static SG::AuxElement::ConstAccessor<float> LeadingClusterSecondR ("LeadingClusterSecondR");
      myjet->LeadingClusterSecondR=LeadingClusterSecondR(*jet);

      static SG::AuxElement::ConstAccessor<char> clean_passLooseBad ("clean_passLooseBad");
      myjet->clean_passLooseBad=clean_passLooseBad(*jet);

      static SG::AuxElement::ConstAccessor<char> clean_passLooseBadUgly ("clean_passLooseBadUgly");
      myjet->clean_passLooseBadUgly=clean_passLooseBadUgly(*jet);

      static SG::AuxElement::ConstAccessor<char> clean_passTightBad ("clean_passTightBad");
      myjet->clean_passTightBad=clean_passTightBad(*jet);

      static SG::AuxElement::ConstAccessor<char> clean_passTightBadUgly ("clean_passTightBadUgly");
      myjet->clean_passTightBadUgly=clean_passTightBadUgly(*jet);

    } // clean


  if ( m_infoSwitch.m_energy ) 
    {
      static SG::AuxElement::ConstAccessor<float> HECFrac ("HECFrac");
      myjet->HECFrac=HECFrac(*jet);

      static SG::AuxElement::ConstAccessor<float> EMFrac ("EMFrac");
      myjet->EMFrac=EMFrac(*jet);

      static SG::AuxElement::ConstAccessor<float> CentroidR ("CentroidR");
      myjet->CentroidR=CentroidR(*jet);

      static SG::AuxElement::ConstAccessor<float> FracSamplingMax ("FracSamplingMax");
      myjet->FracSamplingMax=FracSamplingMax(*jet);

      static SG::AuxElement::ConstAccessor<int> FracSamplingMaxIndex ("FracSamplingMaxIndex");
      myjet->FracSamplingMaxIndex=FracSamplingMaxIndex(*jet);

      static SG::AuxElement::ConstAccessor<float> LowEtConstituentsFrac ("LowEtConstituentsFrac");
      if(LowEtConstituentsFrac.isAvailable( *jet )) myjet->LowEtConstituentsFrac=LowEtConstituentsFrac(*jet);

      static SG::AuxElement::ConstAccessor<int> GhostMuonSegmentCount ("GhostMuonSegmentCount");
      myjet->GhostMuonSegmentCount=GhostMuonSegmentCount(*jet);

      static SG::AuxElement::ConstAccessor<float> Width ("Width");
      if(Width.isAvailable( *jet )) myjet->Width=Width(*jet);

    } // energy


  // each step of the calibration sequence
  if ( m_infoSwitch.m_scales )
    {
      xAOD::JetFourMom_t fourVec;
      jet->getAttribute<xAOD::JetFourMom_t>( "JetEMScaleMomentum",            fourVec );
      myjet->JetEMScaleMomentum.SetPtEtaPhiM(fourVec.Pt()/m_units, fourVec.Eta(), fourVec.Phi(), fourVec.M()/m_units);
      jet->getAttribute<xAOD::JetFourMom_t>( "JetConstitScaleMomentum",       fourVec );
      myjet->JetConstitScaleMomentum.SetPtEtaPhiM(fourVec.Pt()/m_units, fourVec.Eta(), fourVec.Phi(), fourVec.M()/m_units);
      jet->getAttribute<xAOD::JetFourMom_t>( "JetPileupScaleMomentum",        fourVec );
      myjet->JetPileupScaleMomentum.SetPtEtaPhiM(fourVec.Pt()/m_units, fourVec.Eta(), fourVec.Phi(), fourVec.M()/m_units);
      jet->getAttribute<xAOD::JetFourMom_t>( "JetOriginConstitScaleMomentum", fourVec );
      myjet->JetOriginConstitScaleMomentum.SetPtEtaPhiM(fourVec.Pt()/m_units, fourVec.Eta(), fourVec.Phi(), fourVec.M()/m_units);
      jet->getAttribute<xAOD::JetFourMom_t>( "JetEtaJESScaleMomentum",        fourVec );
      myjet->JetEtaJESScaleMomentum.SetPtEtaPhiM(fourVec.Pt()/m_units, fourVec.Eta(), fourVec.Phi(), fourVec.M()/m_units);
      jet->getAttribute<xAOD::JetFourMom_t>( "JetGSCScaleMomentum",           fourVec );
      myjet->JetGSCScaleMomentum.SetPtEtaPhiM(fourVec.Pt()/m_units, fourVec.Eta(), fourVec.Phi(), fourVec.M()/m_units);
      jet->getAttribute<xAOD::JetFourMom_t>( "JetInsituScaleMomentum",        fourVec );
      myjet->JetInsituScaleMomentum.SetPtEtaPhiM(fourVec.Pt()/m_units, fourVec.Eta(), fourVec.Phi(), fourVec.M()/m_units);
    }

  if ( m_infoSwitch.m_layer ) 
    {
      static SG::AuxElement::ConstAccessor< std::vector<float> > EnergyPerSampling ("EnergyPerSampling");
      if ( EnergyPerSampling.isAvailable( *jet ) ) 
	{
	  myjet->EnergyPerSampling=EnergyPerSampling(*jet);
	  std::transform(myjet->EnergyPerSampling.begin(),
			 myjet->EnergyPerSampling.end(),
			 myjet->EnergyPerSampling.begin(),
			 std::bind2nd(std::divides<float>(), m_units));
	} 
      else
	{
	  // could push back a vector of 24...
	  // ... waste of space vs prevention of out of range down stream
	  static const std::vector<float> junk(1,-999);
	  myjet->EnergyPerSampling=junk;
	}
    }

  if ( m_infoSwitch.m_trackAll || m_infoSwitch.m_trackPV ) 
    {

      // several moments calculated from all verticies
      // one accessor for each and just use appropiately in the following
      static SG::AuxElement::ConstAccessor< std::vector<int> >   NumTrkPt1000("NumTrkPt1000");
      static SG::AuxElement::ConstAccessor< std::vector<float> > SumPtTrkPt1000("SumPtTrkPt1000");
      static SG::AuxElement::ConstAccessor< std::vector<float> > TrackWidthPt1000("TrackWidthPt1000");
      static SG::AuxElement::ConstAccessor< std::vector<int> >   NumTrkPt500("NumTrkPt500");
      static SG::AuxElement::ConstAccessor< std::vector<float> > SumPtTrkPt500("SumPtTrkPt500");
      static SG::AuxElement::ConstAccessor< std::vector<float> > TrackWidthPt500("TrackWidthPt500");
      static SG::AuxElement::ConstAccessor< std::vector<float> > JVF("JVF");
    
      if ( m_infoSwitch.m_trackAll ) 
	{
	  if ( NumTrkPt1000.isAvailable( *jet ) ) 
	    myjet->NumTrkPt1000=NumTrkPt1000( *jet );

	  if ( SumPtTrkPt1000.isAvailable( *jet ) )
	    {
	      myjet->SumPtTrkPt1000=SumPtTrkPt1000(*jet);
	      std::transform(myjet->SumPtTrkPt1000.begin(),
			     myjet->SumPtTrkPt1000.end(),
			     myjet->SumPtTrkPt1000.begin(),
			     std::bind2nd(std::divides<float>(), m_units));
	    }

	  if ( TrackWidthPt1000.isAvailable( *jet ) )
	    myjet->TrackWidthPt1000=TrackWidthPt1000( *jet );

	  if ( NumTrkPt500.isAvailable( *jet ) )
	    myjet->NumTrkPt500=NumTrkPt500( *jet );

	  if ( SumPtTrkPt500.isAvailable( *jet ) )
	    {
	      myjet->SumPtTrkPt500=SumPtTrkPt500( *jet );
	      std::transform(myjet->SumPtTrkPt500.begin(),
			     myjet->SumPtTrkPt500.end(),
			     myjet->SumPtTrkPt500.begin(),
			     std::bind2nd(std::divides<float>(), m_units));
	    }

	  if ( TrackWidthPt500.isAvailable( *jet ) )
	    myjet->TrackWidthPt500=TrackWidthPt500( *jet );

	  if ( JVF.isAvailable( *jet ) )
	    myjet->JVF=JVF( *jet );

	} // trackAll

      if ( m_infoSwitch.m_trackPV && pvLocation >= 0 )
	{

	  if ( NumTrkPt1000.isAvailable( *jet ) )
	    myjet->NumTrkPt1000PV=NumTrkPt1000( *jet )[pvLocation];

	  if ( SumPtTrkPt1000.isAvailable( *jet ) )
	    myjet->SumPtTrkPt1000PV=SumPtTrkPt1000( *jet )[pvLocation] / m_units;

	  if ( TrackWidthPt500.isAvailable( *jet ) )
	    myjet->TrackWidthPt500PV=TrackWidthPt500( *jet )[pvLocation];

	  if ( NumTrkPt500.isAvailable( *jet ) )
	    myjet->NumTrkPt500PV=NumTrkPt500( *jet )[pvLocation];

	  if ( SumPtTrkPt500.isAvailable( *jet ) )
	    myjet->SumPtTrkPt500PV=SumPtTrkPt500( *jet )[pvLocation] / m_units;

	  if ( TrackWidthPt500.isAvailable( *jet ) )
	    myjet->TrackWidthPt500PV=TrackWidthPt500( *jet )[pvLocation];


	  if ( JVF.isAvailable( *jet ) )
	    myjet->JVFPV=JVF( *jet )[pvLocation];

	  static SG::AuxElement::ConstAccessor< float > Jvt ("Jvt");
	  myjet->Jvt=Jvt(*jet);

	  static SG::AuxElement::ConstAccessor< float > JvtJvfcorr ("JvtJvfcorr");
	  myjet->JvtJvfcorr=JvtJvfcorr(*jet);

	  static SG::AuxElement::ConstAccessor< float > JvtRpt ("JvtRpt");
	  myjet->JvtRpt=JvtRpt(*jet);

	  if ( m_mc ) 
	    {
	      static SG::AuxElement::ConstAccessor< std::vector< float > > JetJvtEfficiency_JVTSyst_JVT_Loose("JetJvtEfficiency_JVTSyst_JVT_Loose");
	      static SG::AuxElement::ConstAccessor< std::vector< float > > JetJvtEfficiency_JVTSyst_JVT_Medium("JetJvtEfficiency_JVTSyst_JVT_Medium");
	      static SG::AuxElement::ConstAccessor< std::vector< float > > JetJvtEfficiency_JVTSyst_JVT_Tight("JetJvtEfficiency_JVTSyst_JVT_Tight");
	
	      if ( JetJvtEfficiency_JVTSyst_JVT_Loose .isAvailable( *jet ) ) myjet->JetJvtEfficiency_JVTSyst_JVT_Loose = JetJvtEfficiency_JVTSyst_JVT_Loose ( *jet );
	      if ( JetJvtEfficiency_JVTSyst_JVT_Medium.isAvailable( *jet ) ) myjet->JetJvtEfficiency_JVTSyst_JVT_Medium= JetJvtEfficiency_JVTSyst_JVT_Medium( *jet );
	      if ( JetJvtEfficiency_JVTSyst_JVT_Tight .isAvailable( *jet ) ) myjet->JetJvtEfficiency_JVTSyst_JVT_Tight = JetJvtEfficiency_JVTSyst_JVT_Tight ( *jet );
	    }

	} // trackPV

    }

  if ( m_infoSwitch.m_allTrack )
    {
      static SG::AuxElement::ConstAccessor< int > GhostTrackCount("GhostTrackCount");
      myjet->GhostTrackCount=GhostTrackCount( *jet );

      static SG::AuxElement::ConstAccessor< float > GhostTrackPt ("GhostTrackPt");
      myjet->GhostTrackPt=GhostTrackPt( *jet );

      static SG::AuxElement::ConstAccessor< std::vector<ElementLink<DataVector<xAOD::IParticle> > > > GhostTrack ("GhostTrack");
      if ( GhostTrack.isAvailable( *jet ) ) 
	{
	  std::vector<ElementLink<DataVector<xAOD::IParticle> > > trackLinks = GhostTrack( *jet );
	  //std::vector<float> pt(trackLinks.size(),-999);
	  for ( const auto& link_itr : trackLinks ) 
	    {
	      if( !link_itr.isValid() ) { continue; }
	      const xAOD::TrackParticle* track = dynamic_cast<const xAOD::TrackParticle*>( *link_itr );
	      // if asking for tracks passing PV selection ( i.e. JVF JVT tracks )
	      if( m_infoSwitch.m_allTrackPVSel ) 
		{
		  // PV selection from
		  // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JvtManualRecalculation
		  if( track->pt() < 500 )                { continue; } // pT cut
		  if( !m_trkSelTool->accept(*track,pv) ) { continue; } // ID quality cut
		  if( track->vertex() != pv )                          // if not in PV vertex fit
		    {
		      if( track->vertex() != 0 )         { continue; } // make sure in no vertex fits
		      if( fabs((track->z0()+track->vz()-pv->z())*sin(track->theta())) > 3.0 ) { continue; } // make sure close to PV in z
		    }
		}

	      TrackParticle mytrack;
	      mytrack.p4.SetPtEtaPhiE(track->pt() / m_units, track->eta(), track->phi(), track->e()  / m_units );
	      mytrack.qOverP = track->qOverP() * m_units;
	      mytrack.d0 = track->d0();
	      mytrack.z0 = track->z0() + track->vz() - pv->z(); // store z0 wrt PV...most useful
	      if( m_infoSwitch.m_allTrackDetail ) 
		{
		  // n pix, sct, trt
		  track->summaryValue( mytrack.numberOfPixelHits, xAOD::numberOfPixelHits );
		  track->summaryValue( mytrack.numberOfSCTHits  , xAOD::numberOfSCTHits );
		  track->summaryValue( mytrack.numberOfTRTHits  , xAOD::numberOfTRTHits );

		  // pixel split shared
		  track->summaryValue( mytrack.numberOfPixelSharedHits, xAOD::numberOfPixelSharedHits );
		  track->summaryValue( mytrack.numberOfPixelSplitHits , xAOD::numberOfPixelSplitHits );

		  // n ibl, split, shared
		  track->summaryValue( mytrack.numberOfInnermostPixelLayerHits      , xAOD::numberOfInnermostPixelLayerHits );
		  track->summaryValue( mytrack.numberOfInnermostPixelLayerSharedHits, xAOD::numberOfInnermostPixelLayerSharedHits );
		  track->summaryValue( mytrack.numberOfInnermostPixelLayerSplitHits , xAOD::numberOfInnermostPixelLayerSplitHits );

		  // n bl,  split, shared
		  track->summaryValue( mytrack.numberOfNextToInnermostPixelLayerHits      , xAOD::numberOfNextToInnermostPixelLayerHits );
		  track->summaryValue( mytrack.numberOfNextToInnermostPixelLayerSharedHits, xAOD::numberOfNextToInnermostPixelLayerSharedHits );
		  track->summaryValue( mytrack.numberOfNextToInnermostPixelLayerSplitHits , xAOD::numberOfNextToInnermostPixelLayerSplitHits );
		}
	      myjet->GhostTrack.push_back(mytrack);
	    }
	} // if GhostTrack available

    } // allTrack switch

  if( m_infoSwitch.m_constituent ) 
    {
      myjet->numConstituents=jet->numConstituents();
    }

  if( m_infoSwitch.m_constituentAll ) 
    {
      myjet->constituentWeights = jet->getAttribute< std::vector<float> >( "constituentWeights" );

      xAOD::JetConstituentVector consVec = jet->getConstituents();
      if( consVec.isValid() ) 
	{
	  // don't use auto since iterator can also set the scale ...
	  // not sure what that does with auto - probably default but just incase
	  // use the example provided in
	  // http://acode-browser.usatlas.bnl.gov/lxr/source/atlas/Event/xAOD/xAODJet/xAODJet/JetConstituentVector.h
	  JetConstituent myconstit;
	  xAOD::JetConstituentVector::const_iterator constit  = consVec.begin();
	  xAOD::JetConstituentVector::const_iterator constitE = consVec.end();
	  for( ; constit != constitE; constit++)
	    {
	      myconstit.p4.SetPtEtaPhiE(constit->pt() / m_units, constit->eta(), constit->phi(), constit->e() / m_units );
	      myjet->constituents.push_back(myconstit);
	    }
	}
    }

  if ( m_infoSwitch.m_flavTag || m_infoSwitch.m_flavTagHLT ) 
    {
      const xAOD::BTagging * myBTag(0);
    
      if(m_infoSwitch.m_flavTag)
	myBTag = jet->btagging();
      else if(m_infoSwitch.m_flavTagHLT)
	myBTag = jet->auxdata< const xAOD::BTagging* >("HLTBTag");

      if(m_infoSwitch.m_JVC ) 
	{
	  static SG::AuxElement::ConstAccessor<double> JetVertexCharge_discriminant("JetVertexCharge_discriminant");
	  myjet->JetVertexCharge_discriminant = JetVertexCharge_discriminant(*myBTag);
	}

      //MV2c00 MV2c20 MV2c10 MV2c100 MV2m
      myBTag->variable<double>("MV2c00" , "discriminant", myjet->MV2c00);
      myBTag->variable<double>("MV2c10" , "discriminant", myjet->MV2c10);
      myBTag->variable<double>("MV2c20" , "discriminant", myjet->MV2c20);
      myBTag->variable<double>("MV2c100", "discriminant", myjet->MV2c100);

      // flavor groups truth definition
      static SG::AuxElement::ConstAccessor<int> HadronConeExclTruthLabelID("HadronConeExclTruthLabelID");
      myjet->HadronConeExclTruthLabelID = HadronConeExclTruthLabelID( *jet );
      
      if(m_infoSwitch.m_jetFitterDetails )
	{

	  static SG::AuxElement::ConstAccessor< int   > JetFitter_nVTX          ("JetFitter_nVTX");
	  myjet->JetFitter_nVTX = JetFitter_nVTX( *myBTag );

	  static SG::AuxElement::ConstAccessor< int   > JetFitter_nSingleTracks ("JetFitter_nSingleTracks");
	  myjet->JetFitter_nSingleTracks = JetFitter_nSingleTracks( *myBTag );

	  static SG::AuxElement::ConstAccessor< int   > JetFitter_nTracksAtVtx  ("JetFitter_nTracksAtVtx");
	  myjet->JetFitter_nTracksAtVtx = JetFitter_nTracksAtVtx( *myBTag );

	  static SG::AuxElement::ConstAccessor< float > JetFitter_mass          ("JetFitter_mass");
	  myjet->JetFitter_mass = JetFitter_mass( *myBTag );

	  static SG::AuxElement::ConstAccessor< float > JetFitter_energyFraction("JetFitter_energyFraction");
	  myjet->JetFitter_energyFraction = JetFitter_energyFraction( *myBTag );

	  static SG::AuxElement::ConstAccessor< float > JetFitter_significance3d("JetFitter_significance3d");
	  myjet->JetFitter_significance3d = JetFitter_significance3d( *myBTag );

	  static SG::AuxElement::ConstAccessor< float > JetFitter_deltaeta      ("JetFitter_deltaeta");
	  myjet->JetFitter_deltaeta = JetFitter_deltaeta( *myBTag );

	  static SG::AuxElement::ConstAccessor< float > JetFitter_deltaphi      ("JetFitter_deltaphi");
	  myjet->JetFitter_deltaphi = JetFitter_deltaphi( *myBTag );

	  static SG::AuxElement::ConstAccessor< int   > JetFitter_N2Tpair        ("JetFitter_N2Tpair");
	  myjet->JetFitter_N2Tpar = JetFitter_N2Tpair( *myBTag );
	}

      if(m_infoSwitch.m_svDetails ) 
	{
	  /// @brief SV0 : Number of good tracks in vertex
	  static SG::AuxElement::ConstAccessor< int   >   SV0_NGTinSvx     ("SV0_NGTinSvx");
	  myjet->SV0_NGTinSvx=SV0_NGTinSvx( *myBTag );

	  // @brief SV0 : Number of 2-track pairs
	  static SG::AuxElement::ConstAccessor< int   >   SV0_N2Tpair      ("SV0_N2Tpair");
	  myjet->SV0_N2Tpair=SV0_N2Tpair( *myBTag );

	  /// @brief SV0 : vertex mass
	  static SG::AuxElement::ConstAccessor< float   > SV0_masssvx      ("SV0_masssvx");
	  myjet->SV0_masssvx=SV0_masssvx( *myBTag );

	  /// @brief SV0 : energy fraction
	  static SG::AuxElement::ConstAccessor< float   > SV0_efracsvx     ("SV0_efracsvx");           
	  myjet->SV0_efracsvx=SV0_efracsvx( *myBTag );

	  /// @brief SV0 : 3D vertex significance
	  static SG::AuxElement::ConstAccessor< float   > SV0_normdist     ("SV0_normdist");
	  myjet->SV0_normdist=SV0_normdist( *myBTag );

	  myBTag->variable<double>("SV0", "significance3D", myjet->SV0);

	  myjet->SV1plusIP3D_discriminant = myBTag -> SV1plusIP3D_discriminant();

	  /// @brief SV1 : Number of good tracks in vertex
	  static SG::AuxElement::ConstAccessor< int   >   SV1_NGTinSvx     ("SV1_NGTinSvx");
	  myjet->SV1_NGTinSvx=SV1_NGTinSvx( *myBTag );

	  // @brief SV1 : Number of 2-track pairs
	  static SG::AuxElement::ConstAccessor< int   >   SV1_N2Tpair      ("SV1_N2Tpair");
	  myjet->SV1_N2Tpair=SV1_N2Tpair( *myBTag );

	  /// @brief SV1 : vertex mass
	  static SG::AuxElement::ConstAccessor< float   > SV1_masssvx      ("SV1_masssvx");
	  myjet->SV1_masssvx=SV1_masssvx( *myBTag );

	  /// @brief SV1 : energy fraction
	  static SG::AuxElement::ConstAccessor< float   > SV1_efracsvx     ("SV1_efracsvx");
	  myjet->SV1_efracsvx=SV1_efracsvx( *myBTag );

	  /// @brief SV1 : 3D vertex significance
	  static SG::AuxElement::ConstAccessor< float   > SV1_normdist     ("SV1_normdist");
	  myjet->SV1_normdist=SV1_normdist( *myBTag );

	  myjet->SV1_pu = -30;  myBTag->variable<double>("SV1", "pu", myjet->SV1_pu);
	  myjet->SV1_pb = -30;  myBTag->variable<double>("SV1", "pb", myjet->SV1_pb);
	  myjet->SV1_pc = -30;  myBTag->variable<double>("SV1", "pc", myjet->SV1_pc);

	  myjet->SV1    = myBTag->calcLLR(myjet->SV1_pb,myjet->SV1_pu);
	  myjet->SV1_c  = myBTag->calcLLR(myjet->SV1_pb,myjet->SV1_pc);
	  myjet->SV1_cu = myBTag->calcLLR(myjet->SV1_pc,myjet->SV1_pu);

	  myBTag->variable<float>("SV1", "Lxy"         , myjet->SV1_Lxy);
	  myBTag->variable<float>("SV1", "L3d"         , myjet->SV1_L3d);
	  myBTag->variable<float>("SV1", "dstToMatLay" , myjet->SV1_distmatlay);
	  myBTag->variable<float>("SV1", "deltaR"      , myjet->SV1_dR );
	}

      if(m_infoSwitch.m_ipDetails ) 
	{
	  //
	  // IP2D
	  //

	  /// @brief IP2D: track grade
	  static SG::AuxElement::ConstAccessor< std::vector<int>   >   IP2D_gradeOfTracks     ("IP2D_gradeOfTracks");
	  myjet->IP2D_gradeOfTracks=IP2D_gradeOfTracks( *myBTag );

	  /// @brief IP2D : tracks from V0
	  static SG::AuxElement::ConstAccessor< std::vector<bool>   >  IP2D_flagFromV0ofTracks("IP2D_flagFromV0ofTracks");
	  myjet->IP2D_flagFromV0ofTracks=IP2D_flagFromV0ofTracks( *myBTag );

	  /// @brief IP2D : d0 value with respect to primary vertex
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP2D_valD0wrtPVofTracks("IP2D_valD0wrtPVofTracks");
	  myjet->IP2D_valD0wrtPVofTracks=IP2D_valD0wrtPVofTracks( *myBTag );

	  /// @brief IP2D : d0 significance with respect to primary vertex
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP2D_sigD0wrtPVofTracks("IP2D_sigD0wrtPVofTracks");
	  myjet->IP2D_sigD0wrtPVofTracks=IP2D_sigD0wrtPVofTracks( *myBTag );

	  /// @brief IP2D : track contribution to B likelihood
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP2D_weightBofTracks   ("IP2D_weightBofTracks");
	  myjet->IP2D_weightBofTracks=IP2D_weightBofTracks( *myBTag );

	  /// @brief IP2D : track contribution to C likelihood
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP2D_weightCofTracks   ("IP2D_weightCofTracks");
	  myjet->IP2D_weightCofTracks=IP2D_weightCofTracks( *myBTag );

	  /// @brief IP2D : track contribution to U likelihood
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP2D_weightUofTracks   ("IP2D_weightUofTracks");
	  myjet->IP2D_weightUofTracks=IP2D_weightUofTracks( *myBTag );

	  myjet->IP2D_pu = -99;  myBTag->variable<double>("IP2D", "pu", myjet->IP2D_pu);
	  myjet->IP2D_pb = -99;  myBTag->variable<double>("IP2D", "pb", myjet->IP2D_pb);
	  myjet->IP2D_pc = -99;  myBTag->variable<double>("IP2D", "pc", myjet->IP2D_pc);

	  myjet->IP2D    = myBTag->calcLLR(myjet->IP2D_pb,myjet->IP2D_pu);
	  myjet->IP2D_c  = myBTag->calcLLR(myjet->IP2D_pb,myjet->IP2D_pc);
	  myjet->IP2D_cu = myBTag->calcLLR(myjet->IP2D_pc,myjet->IP2D_pu);

	  //
	  // IP3D
	  //

	  /// @brief IP3D: track grade
	  static SG::AuxElement::ConstAccessor< std::vector<int>   >   IP3D_gradeOfTracks     ("IP3D_gradeOfTracks");
	  myjet->IP3D_gradeOfTracks=IP3D_gradeOfTracks( *myBTag );

	  /// @brief IP3D : tracks from V0
	  static SG::AuxElement::ConstAccessor< std::vector<bool>   >  IP3D_flagFromV0ofTracks("IP3D_flagFromV0ofTracks");
	  myjet->IP3D_flagFromV0ofTracks=IP3D_flagFromV0ofTracks( *myBTag );

	  /// @brief IP3D : d0 value with respect to primary vertex
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP3D_valD0wrtPVofTracks("IP3D_valD0wrtPVofTracks");
          myjet->IP3D_valD0wrtPVofTracks=IP3D_valD0wrtPVofTracks( *myBTag );

	  /// @brief IP3D : d0 significance with respect to primary vertex
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP3D_sigD0wrtPVofTracks("IP3D_sigD0wrtPVofTracks");
	  myjet->IP3D_sigD0wrtPVofTracks=IP3D_sigD0wrtPVofTracks( *myBTag );

	  /// @brief IP3D : z0 value with respect to primary vertex
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP3D_valZ0wrtPVofTracks("IP3D_valZ0wrtPVofTracks");
	  myjet->IP3D_valZ0wrtPVofTracks=IP3D_valZ0wrtPVofTracks( *myBTag );

	  /// @brief IP3D : z0 significance with respect to primary vertex
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP3D_sigZ0wrtPVofTracks("IP3D_sigZ0wrtPVofTracks");
	  myjet->IP3D_sigZ0wrtPVofTracks=IP3D_sigZ0wrtPVofTracks( *myBTag );

	  /// @brief IP3D : track contribution to B likelihood
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP3D_weightBofTracks   ("IP3D_weightBofTracks");
	  myjet->IP3D_weightBofTracks=IP3D_weightBofTracks( *myBTag );

	  /// @brief IP3D : track contribution to C likelihood
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP3D_weightCofTracks   ("IP3D_weightCofTracks");
	  myjet->IP3D_weightCofTracks=IP3D_weightCofTracks( *myBTag );

	  /// @brief IP3D : track contribution to U likelihood
	  static SG::AuxElement::ConstAccessor< std::vector<float>   > IP3D_weightUofTracks   ("IP3D_weightUofTracks");
	  myjet->IP3D_weightUofTracks=IP3D_weightUofTracks( *myBTag );

	  myjet->IP3D_pu = -30;  myBTag->variable<double>("IP3D", "pu", myjet->IP3D_pu);
	  myjet->IP3D_pb = -30;  myBTag->variable<double>("IP3D", "pb", myjet->IP3D_pb);
	  myjet->IP3D_pc = -30;  myBTag->variable<double>("IP3D", "pc", myjet->IP3D_pc);

	  myjet->IP3D_loglikelihoodratio=myBTag -> IP3D_loglikelihoodratio();

	  myjet->IP3D    = myBTag->calcLLR(myjet->IP3D_pb,myjet->IP3D_pu);
	  myjet->IP3D_c  = myBTag->calcLLR(myjet->IP3D_pb,myjet->IP3D_pc);
	  myjet->IP3D_cu = myBTag->calcLLR(myjet->IP3D_pc,myjet->IP3D_pu);

	}

      if(m_infoSwitch.m_flavTagHLT ) 
	{
	  const xAOD::Vertex *online_pvx       = jet->auxdata<const xAOD::Vertex*>("HLTBJetTracks_vtx");
	  const xAOD::Vertex *online_pvx_bkg   = jet->auxdata<const xAOD::Vertex*>("HLTBJetTracks_vtx_bkg");
	  const xAOD::Vertex *offline_pvx      = jet->auxdata<const xAOD::Vertex*>("offline_vtx");      

	  myjet->vtxOnlineValid = (online_pvx)?1.0:0.0;

	  char hadDummyPV = jet->auxdata< char >("hadDummyPV");
	       if( hadDummyPV == '0')  myjet->vtxHadDummy=0.0;
	  else if( hadDummyPV == '1')  myjet->vtxHadDummy=1.0;
	  else if( hadDummyPV == '2')  myjet->vtxHadDummy=2.0;

	  static SG::AuxElement::ConstAccessor< float > bs_online_vz ("bs_online_vz");
	  if(bs_online_vz.isAvailable( *jet) )
	    {
	      myjet->bs_online_vz = jet->auxdata< float >("bs_online_vz");
	      myjet->bs_online_vx = jet->auxdata< float >("bs_online_vx");
	      myjet->bs_online_vy = jet->auxdata< float >("bs_online_vy");
	    }

	  myjet->vtx_offline_x0 = offline_pvx->x();
	  myjet->vtx_offline_y0 = offline_pvx->y();
	  myjet->vtx_offline_z0 = offline_pvx->z();

	  if(online_pvx)
	    {
	      myjet->vtx_online_x0 = online_pvx->x();
	      myjet->vtx_online_y0 = online_pvx->y();
	      myjet->vtx_online_z0 = online_pvx->z();
	    }

	  if(online_pvx_bkg)
	    {
	      myjet->vtx_online_bkg_x0 = online_pvx_bkg->x();
	      myjet->vtx_online_bkg_y0 = online_pvx_bkg->y();
	      myjet->vtx_online_bkg_z0 = online_pvx_bkg->z();
	    }

	}// m_flavTagHLT
    }

  if( !m_infoSwitch.m_sfFTagFix.empty() ) 
    {
      for( unsigned int i=0; i<m_infoSwitch.m_sfFTagFix.size(); i++ ) 
	{
	  switch( m_infoSwitch.m_sfFTagFix.at(i) )
	    {
	    case 30 :
	      m_btag_Fix30->Fill( jet ); 
	      if(        m_btag_Fix30->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix30->m_isTag( *jet );
	      if(m_mc && m_btag_Fix30->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix30=m_btag_Fix30->m_sf   ( *jet );
	      break;
	    case 50 :
	      m_btag_Fix50->Fill( jet ); 
	      if(        m_btag_Fix50->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix50->m_isTag( *jet );
	      if(m_mc && m_btag_Fix50->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix50=m_btag_Fix50->m_sf   ( *jet );
	      break;
	    case 60 :
	      m_btag_Fix60->Fill( jet ); 
	      if(        m_btag_Fix60->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix60->m_isTag( *jet );
	      if(m_mc && m_btag_Fix60->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix60=m_btag_Fix60->m_sf   ( *jet );
	      break;
	    case 70 :
	      m_btag_Fix70->Fill( jet ); 
	      if(        m_btag_Fix70->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix70->m_isTag( *jet );
	      if(m_mc && m_btag_Fix70->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix70=m_btag_Fix70->m_sf   ( *jet );
	      break;
	    case 77 :
	      m_btag_Fix77->Fill( jet ); 
	      if(        m_btag_Fix77->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix77->m_isTag( *jet );
	      if(m_mc && m_btag_Fix77->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix77=m_btag_Fix77->m_sf   ( *jet );
	      break;
	    case 80 :
	      m_btag_Fix80->Fill( jet ); 
	      if(        m_btag_Fix80->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix80->m_isTag( *jet );
	      if(m_mc && m_btag_Fix80->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix80=m_btag_Fix80->m_sf   ( *jet );
	      break;
	    case 85 :
	      m_btag_Fix85->Fill( jet ); 
	      if(        m_btag_Fix85->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix85->m_isTag( *jet );
	      if(m_mc && m_btag_Fix85->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix85=m_btag_Fix85->m_sf   ( *jet );
	      break;
	    case 90 :
	      m_btag_Fix90->Fill( jet ); 
	      if(        m_btag_Fix90->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Fix90->m_isTag( *jet );
	      if(m_mc && m_btag_Fix90->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFix90=m_btag_Fix90->m_sf   ( *jet );
	      break;
	    }
	}
    } // sfFTagFix

  if( !m_infoSwitch.m_sfFTagFlt.empty() ) 
    {
      for( unsigned int i=0; i<m_infoSwitch.m_sfFTagFlt.size(); i++ ) 
	{
	  switch( m_infoSwitch.m_sfFTagFlt.at(i) ) 
	    {
	    case 30 :
	      m_btag_Flt30->Fill( jet ); 
	      if(        m_btag_Flt30->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Flt30->m_isTag( *jet );
	      if(m_mc && m_btag_Flt30->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFlt30=m_btag_Flt30->m_sf   ( *jet );
	      break;
	    case 50 :
	      m_btag_Flt50->Fill( jet ); 
	      if(        m_btag_Flt50->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Flt50->m_isTag( *jet );
	      if(m_mc && m_btag_Flt50->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFlt50=m_btag_Flt50->m_sf   ( *jet );
	      break;
	    case 60 :
	      m_btag_Flt60->Fill( jet ); 
	      if(        m_btag_Flt60->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Flt60->m_isTag( *jet );
	      if(m_mc && m_btag_Flt60->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFlt60=m_btag_Flt60->m_sf   ( *jet );
	      break;
	    case 70 :
	      m_btag_Flt70->Fill( jet ); 
	      if(        m_btag_Flt70->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Flt70->m_isTag( *jet );
	      if(m_mc && m_btag_Flt70->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFlt70=m_btag_Flt70->m_sf   ( *jet );
	      break;
	    case 77 :
	      m_btag_Flt77->Fill( jet ); 
	      if(        m_btag_Flt77->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Flt77->m_isTag( *jet );
	      if(m_mc && m_btag_Flt77->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFlt77=m_btag_Flt77->m_sf   ( *jet );
	      break;
	    case 85 :
	      m_btag_Flt85->Fill( jet ); 
	      if(        m_btag_Flt85->m_isTag.isAvailable( *jet )) myjet->MV2c20        =m_btag_Flt85->m_isTag( *jet );
	      if(m_mc && m_btag_Flt85->m_sf   .isAvailable( *jet )) myjet->MV2c20_sfFlt85=m_btag_Flt85->m_sf   ( *jet );
	      break;
	    }
	}
    } // sfFTagFlt

  if ( m_infoSwitch.m_area ) 
    {

      static SG::AuxElement::ConstAccessor<float> JetGhostArea("JetGhostArea");
      myjet->JetGhostArea=JetGhostArea( *jet );

      static SG::AuxElement::ConstAccessor<float> ActiveArea("ActiveArea");
      myjet->ActiveArea=ActiveArea( *jet );

      static SG::AuxElement::ConstAccessor<float> VoronoiArea("VoronoiArea");
      myjet->VoronoiArea=VoronoiArea( *jet );

      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_pt("ActiveArea4vec_pt");
      myjet->ActiveArea4vec_pt=ActiveArea4vec_pt( *jet );

      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_eta("ActiveArea4vec_eta");
      myjet->ActiveArea4vec_eta=ActiveArea4vec_eta( *jet );

      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_phi("ActiveArea4vec_phi");
      myjet->ActiveArea4vec_phi=ActiveArea4vec_phi( *jet );

      static SG::AuxElement::ConstAccessor<float> ActiveArea4vec_m("ActiveArea4vec_m");
      myjet->ActiveArea4vec_m=ActiveArea4vec_m( *jet );
    }

  if ( m_infoSwitch.m_truth && m_mc ) 
    {
      static SG::AuxElement::ConstAccessor<int> ConeTruthLabelID ("ConeTruthLabelID");
      myjet->ConeTruthLabelID=ConeTruthLabelID( *jet );

      static SG::AuxElement::ConstAccessor<int> TruthCount ("TruthCount");
      if(TruthCount.isAvailable( *jet )) myjet->TruthCount=TruthCount( *jet );

      //    seems to be empty
      //      static SG::AuxElement::ConstAccessor<float> TruthPt ("TruthPt");
      //      if ( TruthPt.isAvailable( *jet) ) {
      //        m_truthPt->push_back( TruthPt( *jet)/1000 );
      //      } else { m_truthPt->push_back( -999 ); }

      static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_B ("TruthLabelDeltaR_B");
      SAFE_SET(myjet,TruthLabelDeltaR_B,jet);

      static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_C ("TruthLabelDeltaR_C");
      SAFE_SET(myjet,TruthLabelDeltaR_C,jet);

      static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_T ("TruthLabelDeltaR_T");
      SAFE_SET(myjet,TruthLabelDeltaR_T,jet);

      static SG::AuxElement::ConstAccessor<int> PartonTruthLabelID("PartonTruthLabelID");
      SAFE_SET(myjet,PartonTruthLabelID,jet);

      static SG::AuxElement::ConstAccessor<float> GhostTruthAssociationFraction("GhostTruthAssociationFraction");
      SAFE_SET(myjet,GhostTruthAssociationFraction,jet);

      const xAOD::Jet* GhostTruthAssociationLink = HelperFunctions::getLink<xAOD::Jet>( jet, "GhostTruthAssociationLink" );
      if(GhostTruthAssociationLink) 
	myjet->truth_p4.SetPtEtaPhiE(GhostTruthAssociationLink->pt() / m_units,
				     GhostTruthAssociationLink->eta(),
				     GhostTruthAssociationLink->phi(),
				     GhostTruthAssociationLink->e() / m_units);
    }

  if ( m_infoSwitch.m_truthDetails ) 
    {
      //
      // B-Hadron Details
      //
      static SG::AuxElement::ConstAccessor<int> GhostBHadronsFinalCount ("GhostBHadronsFinalCount");
      myjet->GhostBHadronsFinalCount=GhostBHadronsFinalCount( *jet );

      static SG::AuxElement::ConstAccessor<int> GhostBHadronsInitialCount ("GhostBHadronsInitialCount");
      myjet->GhostBHadronsInitialCount=GhostBHadronsInitialCount( *jet );

      static SG::AuxElement::ConstAccessor<int> GhostBQuarksFinalCount ("GhostBQuarksFinalCount");
      myjet->GhostBQuarksFinalCount=GhostBQuarksFinalCount( *jet );

      static SG::AuxElement::ConstAccessor<float> GhostBHadronsFinalPt ("GhostBHadronsFinalPt");
      myjet->GhostBHadronsFinalPt=GhostBHadronsFinalPt( *jet );

      static SG::AuxElement::ConstAccessor<float> GhostBHadronsInitialPt ("GhostBHadronsInitialPt");
      myjet->GhostBHadronsInitialPt=GhostBHadronsInitialPt( *jet );

      static SG::AuxElement::ConstAccessor<float> GhostBQuarksFinalPt ("GhostBQuarksFinalPt");
      myjet->GhostBQuarksFinalPt=GhostBQuarksFinalPt( *jet );

      //
      // C-Hadron Details
      //
      static SG::AuxElement::ConstAccessor<int> GhostCHadronsFinalCount ("GhostCHadronsFinalCount");
      myjet->GhostCHadronsFinalCount=GhostCHadronsFinalCount( *jet );

      static SG::AuxElement::ConstAccessor<int> GhostCHadronsInitialCount ("GhostCHadronsInitialCount");
      myjet->GhostCHadronsInitialCount=GhostCHadronsInitialCount( *jet );

      static SG::AuxElement::ConstAccessor<int> GhostCQuarksFinalCount ("GhostCQuarksFinalCount");
      myjet->GhostCQuarksFinalCount=GhostCQuarksFinalCount( *jet );

      static SG::AuxElement::ConstAccessor<float> GhostCHadronsFinalPt ("GhostCHadronsFinalPt");
      myjet->GhostCHadronsFinalPt=GhostCHadronsFinalPt( *jet );

      static SG::AuxElement::ConstAccessor<float> GhostCHadronsInitialPt ("GhostCHadronsInitialPt");
      myjet->GhostCHadronsInitialPt=GhostCHadronsInitialPt( *jet );

      static SG::AuxElement::ConstAccessor<float> GhostCQuarksFinalPt ("GhostCQuarksFinalPt");
      myjet->GhostCQuarksFinalPt=GhostCQuarksFinalPt( *jet );

      //
      // Tau Details
      //
      static SG::AuxElement::ConstAccessor<int> GhostTausFinalCount ("GhostTausFinalCount");
      myjet->GhostTausFinalCount=GhostTausFinalCount( *jet );

      // THE ONLY UN-OFFICIAL PIECE OF CODE HERE USE WITH CAUTION
      static SG::AuxElement::ConstAccessor<float> GhostTausFinalPt ("GhostTausFinalPt");
      myjet->GhostTausFinalPt=GhostTausFinalPt( *jet );

      // light quark(1,2,3) , gluon (21 or 9), charm(4) and b(5)
      // GhostPartons should select for these pdgIds only
      //    static SG::AuxElement::ConstAccessor< std::vector<const xAOD::TruthParticle*> > ghostPartons("GhostPartons");
      //    if( ghostPartons.isAvailable( *jet )) {
      //    std::vector<const xAOD::TruthParticle*> truthPartons = ghostPartons( *jet );

      std::vector<const xAOD::TruthParticle*> GhostPartons = jet->getAssociatedObjects<xAOD::TruthParticle>("GhostPartons");

      if( GhostPartons.size() > 0)
	{
	  int iParent = 0;
	  for(unsigned int i=1; i < GhostPartons.size(); ++i)
	    {
	      if( (GhostPartons.at(i)->pt() > 0.001) && (GhostPartons.at(i)->e() > GhostPartons.at(iParent)->e()) )
		iParent = i;
	    }
	  myjet->truth_pdgId   =GhostPartons[iParent]->pdgId();
	  myjet->truth_partonPt=GhostPartons[iParent]->pt() / m_units;
	  myjet->truth_partonDR=GhostPartons[iParent]->p4().DeltaR( jet->p4() );
	}
    }

  if ( m_infoSwitch.m_charge ) 
    {
      xAOD::JetFourMom_t p4UsedInJetCharge;
      bool status = jet->getAttribute<xAOD::JetFourMom_t>( "JetPileupScaleMomentum", p4UsedInJetCharge );
      static SG::AuxElement::ConstAccessor<float> Charge ("Charge");

      if(status)
	{
	  float ptUsedInJetCharge   = p4UsedInJetCharge.Pt();
	  float calibratedJetCharge = jet->pt() ? (ptUsedInJetCharge * Charge(*jet) / jet->pt()) : -99;
	  myjet->charge=calibratedJetCharge;
	}
    }

  return;
}


void JetHelpTree::fillGlobalBTagSF( const xAOD::EventInfo* eventInfo ){

  if( !m_infoSwitch.m_sfFTagFix.empty() ) {
    for( unsigned int i=0; i<m_infoSwitch.m_sfFTagFix.size(); i++ ) {
      switch( m_infoSwitch.m_sfFTagFix.at(i) ) {
      case 30 :  m_btag_Fix30->FillGlobalSF(eventInfo); break;
      case 50 :  m_btag_Fix50->FillGlobalSF(eventInfo); break;
      case 60 :  m_btag_Fix60->FillGlobalSF(eventInfo); break;
      case 70 :  m_btag_Fix70->FillGlobalSF(eventInfo); break;
      case 77 :  m_btag_Fix77->FillGlobalSF(eventInfo); break;
      case 80 :  m_btag_Fix80->FillGlobalSF(eventInfo); break;
      case 85 :  m_btag_Fix85->FillGlobalSF(eventInfo); break;
      case 90 :  m_btag_Fix90->FillGlobalSF(eventInfo); break;
      }
    }
  } // sfFTagFix

  if( !m_infoSwitch.m_sfFTagFlt.empty() ) {
    for( unsigned int i=0; i<m_infoSwitch.m_sfFTagFlt.size(); i++ ) {
      switch( m_infoSwitch.m_sfFTagFlt.at(i) ) {
      case 30 :  m_btag_Flt30->FillGlobalSF(eventInfo); break;
      case 50 :	 m_btag_Flt50->FillGlobalSF(eventInfo); break;
      case 60 :	 m_btag_Flt60->FillGlobalSF(eventInfo); break;
      case 70 :	 m_btag_Flt70->FillGlobalSF(eventInfo); break;
      case 77 :	 m_btag_Flt77->FillGlobalSF(eventInfo); break;
      case 85 :	 m_btag_Flt85->FillGlobalSF(eventInfo); break;
      }	 
    }
  } // sfFTagFlt

  return;
}

bool JetHelpTree::haveBTagSF(const std::vector<int>& sfList, int workingPt){
  return (std::find(sfList.begin(), sfList.end(),workingPt ) != sfList.end());
 }

