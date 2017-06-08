#include "xAODAnaHelpers/JetHelpTree.h"
#include <xAODAnaHelpers/HelperFunctions.h>
#include <iostream>
#include "xAODTruth/TruthEventContainer.h"

using namespace xAH;
using std::vector;  using std::endl;  using std::cout;


JetHelpTree::JetHelpTree(const std::string& name, const std::string& detailStr, float units, bool mc)
  : ParticleHelpTree("xAH::Jet",name,detailStr,units,mc),
    m_trkSelTool(nullptr)

{ }

JetHelpTree::~JetHelpTree()
{ }

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
	  setBranchStatus(tree, "JvtEff_SF_Loose", 1);
	  setBranchStatus(tree, "JvtEff_SF_Medium",1);
	  setBranchStatus(tree, "JvtEff_SF_Tight", 1);
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
      setBranchStatus(tree, "coynstituent_pt",     1);
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
	  setBranchStatus(tree, "sv0_NGTinSvx",      1);
	  setBranchStatus(tree, "sv0_N2Tpair",       1);
	  setBranchStatus(tree, "sv0_massvx",        1);
	  setBranchStatus(tree, "sv0_efracsvx",      1);
	  setBranchStatus(tree, "sv0_normdist",      1);

	  setBranchStatus(tree, "SV1",               1);
	  setBranchStatus(tree, "SV1IP3D",           1);
	  setBranchStatus(tree, "sv1_pu",            1);
	  setBranchStatus(tree, "sv1_pb",            1);
	  setBranchStatus(tree, "sv1_pc",            1);
	  setBranchStatus(tree, "sv1_c",             1);
	  setBranchStatus(tree, "sv1_cu",            1);
	  setBranchStatus(tree, "sv1_NGTinSvx",      1);
	  setBranchStatus(tree, "sv1_N2Tpair",       1);
	  setBranchStatus(tree, "sv1_massvx",        1);
	  setBranchStatus(tree, "sv1_efracsvx",      1);
	  setBranchStatus(tree, "sv1_normdist",      1);
	  setBranchStatus(tree, "sv1_Lxy",           1);
	  setBranchStatus(tree, "sv1_L3d",           1);
	  setBranchStatus(tree, "sv1_distmatlay",    1);
	  setBranchStatus(tree, "sv1_dR",            1);
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
      setBranchStatus(tree, "truth_E",   1);
      setBranchStatus(tree, "truth_pt",  1);
      setBranchStatus(tree, "truth_phi", 1);
      setBranchStatus(tree, "truth_eta", 1);
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

void JetHelpTree::FillJet( const xAOD::Jet* jet, const xAOD::Vertex* pv, int pvLocation )
{
  return FillJet(static_cast<const xAOD::IParticle*>(jet), pv, pvLocation);
}

void JetHelpTree::FillJet( const xAOD::IParticle* particle, const xAOD::Vertex* pv, int pvLocation ){
  if(m_debug) cout << "In JetHelpTree::FillJet " << endl;

  ParticleHelpTree::FillParticle(particle);

  const xAOD::Jet* jet=dynamic_cast<const xAOD::Jet*>(particle);

  if( m_infoSwitch.m_rapidity ){
    m_rapidity->push_back( jet->rapidity() );
  }

  if (m_infoSwitch.m_clean) {
    static SG::AuxElement::ConstAccessor<float> jetTime ("Timing");
    safeFill<float, float, xAOD::Jet>(jet, jetTime, m_Timing, -999);

    static SG::AuxElement::ConstAccessor<float> LArQuality ("LArQuality");
    safeFill<float, float, xAOD::Jet>(jet, LArQuality, m_LArQuality, -999);

    static SG::AuxElement::ConstAccessor<float> hecq ("HECQuality");
    safeFill<float, float, xAOD::Jet>(jet, hecq, m_HECQuality, -999);

    static SG::AuxElement::ConstAccessor<float> negE ("NegativeE");
    safeFill<float, float, xAOD::Jet>(jet, negE, m_NegativeE, -999, m_units);

    static SG::AuxElement::ConstAccessor<float> avLArQF ("AverageLArQF");
    safeFill<float, float, xAOD::Jet>(jet, avLArQF, m_AverageLArQF, -999);

    static SG::AuxElement::ConstAccessor<float> bchCorrCell ("BchCorrCell");
    safeFill<float, float, xAOD::Jet>(jet, bchCorrCell, m_BchCorrCell, -999);

    static SG::AuxElement::ConstAccessor<float> N90Const ("N90Constituents");
    safeFill<float, float, xAOD::Jet>(jet, N90Const, m_N90Constituents, -999);

    static SG::AuxElement::ConstAccessor<float> LArBadHVEFrac ("LArBadHVEnergyFrac");
    safeFill<float, float, xAOD::Jet>(jet, LArBadHVEFrac, m_LArBadHVEnergyFrac, -999);

    static SG::AuxElement::ConstAccessor<int> LArBadHVNCell ("LArBadHVNCell");
    safeFill<int, int, xAOD::Jet>(jet, LArBadHVNCell, m_LArBadHVNCell, -999);

    static SG::AuxElement::ConstAccessor<float> OotFracClus5 ("OotFracClusters5");
    safeFill<float, float, xAOD::Jet>(jet, OotFracClus5, m_OotFracClusters5, -999);

    static SG::AuxElement::ConstAccessor<float> OotFracClus10 ("OotFracClusters10");
    safeFill<float, float, xAOD::Jet>(jet, OotFracClus10, m_OotFracClusters10, -999);

    static SG::AuxElement::ConstAccessor<float> leadClusPt ("LeadingClusterPt");
    safeFill<float, float, xAOD::Jet>(jet, leadClusPt, m_LeadingClusterPt, -999);

    static SG::AuxElement::ConstAccessor<float> leadClusSecondLambda ("LeadingClusterSecondLambda");
    safeFill<float, float, xAOD::Jet>(jet, leadClusSecondLambda, m_LeadingClusterSecondLambda, -999);

    static SG::AuxElement::ConstAccessor<float> leadClusCenterLambda ("LeadingClusterCenterLambda");
    safeFill<float, float, xAOD::Jet>(jet, leadClusCenterLambda, m_LeadingClusterCenterLambda, -999);

    static SG::AuxElement::ConstAccessor<float> leadClusSecondR ("LeadingClusterSecondR");
    safeFill<float, float, xAOD::Jet>(jet, leadClusSecondR, m_LeadingClusterSecondR, -999);

    static SG::AuxElement::ConstAccessor<char> clean_passLooseBad ("clean_passLooseBad");
    safeFill<char, int, xAOD::Jet>(jet, clean_passLooseBad, m_clean_passLooseBad, -999);

    static SG::AuxElement::ConstAccessor<char> clean_passLooseBadUgly ("clean_passLooseBadUgly");
    safeFill<char, int, xAOD::Jet>(jet, clean_passLooseBadUgly, m_clean_passLooseBadUgly, -999);

    static SG::AuxElement::ConstAccessor<char> clean_passTightBad ("clean_passTightBad");
    safeFill<char, int, xAOD::Jet>(jet, clean_passTightBad, m_clean_passTightBad, -999);

    static SG::AuxElement::ConstAccessor<char> clean_passTightBadUgly ("clean_passTightBadUgly");
    safeFill<char, int, xAOD::Jet>(jet, clean_passTightBadUgly, m_clean_passTightBadUgly, -999);

  } // clean


  if ( m_infoSwitch.m_energy ) {
    static SG::AuxElement::ConstAccessor<float> HECf ("HECFrac");
    safeFill<float, float, xAOD::Jet>(jet, HECf, m_HECFrac, -999);

    static SG::AuxElement::ConstAccessor<float> EMf ("EMFrac");
    safeFill<float, float, xAOD::Jet>(jet, EMf, m_EMFrac, -999);

    static SG::AuxElement::ConstAccessor<float> centroidR ("CentroidR");
    safeFill<float, float, xAOD::Jet>(jet, centroidR, m_CentroidR, -999);

    static SG::AuxElement::ConstAccessor<float> fracSampMax ("FracSamplingMax");
    safeFill<float, float, xAOD::Jet>(jet, fracSampMax, m_FracSamplingMax, -999);

    static SG::AuxElement::ConstAccessor<int> fracSampMaxIdx ("FracSamplingMaxIndex");
    safeFill<int, float, xAOD::Jet>(jet, fracSampMaxIdx, m_FracSamplingMaxIndex, -999);

    static SG::AuxElement::ConstAccessor<float> lowEtFrac ("LowEtConstituentsFrac");
    safeFill<float, float, xAOD::Jet>(jet, lowEtFrac, m_LowEtConstituentsFrac, -999);

    static SG::AuxElement::ConstAccessor<int> muonSegCount ("GhostMuonSegmentCount");
    safeFill<int, float, xAOD::Jet>(jet, muonSegCount, m_GhostMuonSegmentCount, -999);

    static SG::AuxElement::ConstAccessor<float> width ("Width");
    safeFill<float, float, xAOD::Jet>(jet, width, m_Width, -999);

  } // energy


  // each step of the calibration sequence
  if ( m_infoSwitch.m_scales ) {
    xAOD::JetFourMom_t fourVec;
    bool status(false);
    // EM Scale
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetEMScaleMomentum", fourVec );
    if( status ) { m_emScalePt->push_back( fourVec.Pt() / m_units ); }
    else { m_emScalePt->push_back( -999 ); }
    // Constit Scale
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetConstitScaleMomentum", fourVec );
    if( status ) { m_constScalePt->push_back( fourVec.Pt() / m_units ); }
    else { m_constScalePt->push_back( -999 ); }
    // Pileup Scale
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetPileupScaleMomentum", fourVec );
    if( status ) { m_pileupScalePt->push_back( fourVec.Pt() / m_units ); }
    else { m_pileupScalePt->push_back( -999 ); }
    // OriginConstit Scale
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetOriginConstitScaleMomentum", fourVec );
    if( status ) { m_originConstitScalePt->push_back( fourVec.Pt() / m_units ); }
    else { m_originConstitScalePt->push_back( -999 ); }
    // EtaJES Scale
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetEtaJESScaleMomentum", fourVec );
    if( status ) { m_etaJESScalePt->push_back( fourVec.Pt() / m_units ); }
    else { m_etaJESScalePt->push_back( -999 ); }
    // GSC Scale
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetGSCScaleMomentum", fourVec );
    if( status ) { m_gscScalePt->push_back( fourVec.Pt() / m_units ); }
    else { m_gscScalePt->push_back( -999 ); }
    // only available in data
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetInsituScaleMomentum", fourVec );
    if(status) { m_insituScalePt->push_back( fourVec.Pt() / m_units ); }
    else { m_insituScalePt->push_back( -999 ); }
  }

  if ( m_infoSwitch.m_constscaleEta ) {
    xAOD::JetFourMom_t fourVec;
    bool status(false);
    status = jet->getAttribute<xAOD::JetFourMom_t>( "JetConstitScaleMomentum", fourVec );
    if( status ) { m_constScaleEta->push_back( fourVec.Eta() ); }
    else { m_constScaleEta->push_back( -999 ); }
  }

  if ( m_infoSwitch.m_layer ) {
    static SG::AuxElement::ConstAccessor< std::vector<float> > ePerSamp ("EnergyPerSampling");
    if ( ePerSamp.isAvailable( *jet ) ) {
      m_EnergyPerSampling->push_back( ePerSamp( *jet ) );
      m_EnergyPerSampling->back();
      std::transform((m_EnergyPerSampling->back()).begin(),
                     (m_EnergyPerSampling->back()).end(),
                     (m_EnergyPerSampling->back()).begin(),
                     std::bind2nd(std::divides<float>(), m_units));
    } else {
      // could push back a vector of 24...
      // ... waste of space vs prevention of out of range down stream
      std::vector<float> junk(1,-999);
      m_EnergyPerSampling->push_back( junk );
    }
  }

  if ( m_infoSwitch.m_trackAll || m_infoSwitch.m_trackPV ) {

    // several moments calculated from all verticies
    // one accessor for each and just use appropiately in the following
    static SG::AuxElement::ConstAccessor< std::vector<int> >   nTrk1000("NumTrkPt1000");
    static SG::AuxElement::ConstAccessor< std::vector<float> > sumPt1000("SumPtTrkPt1000");
    static SG::AuxElement::ConstAccessor< std::vector<float> > trkWidth1000("TrackWidthPt1000");
    static SG::AuxElement::ConstAccessor< std::vector<int> >   nTrk500 ("NumTrkPt500");
    static SG::AuxElement::ConstAccessor< std::vector<float> > sumPt500 ("SumPtTrkPt500");
    static SG::AuxElement::ConstAccessor< std::vector<float> > trkWidth500 ("TrackWidthPt500");
    static SG::AuxElement::ConstAccessor< std::vector<float> > jvf("JVF");
    
    if ( m_infoSwitch.m_trackAll ) {

      std::vector<int> junkInt(1,-999);
      std::vector<float> junkFlt(1,-999);

      if ( nTrk1000.isAvailable( *jet ) ) {
        m_NumTrkPt1000->push_back( nTrk1000( *jet ) );
      } else { m_NumTrkPt1000->push_back( junkInt ); }

      if ( sumPt1000.isAvailable( *jet ) ) {
        m_SumPtTrkPt1000->push_back( sumPt1000( *jet ) );
        std::transform((m_SumPtTrkPt1000->back()).begin(),
                       (m_SumPtTrkPt1000->back()).end(),
                       (m_SumPtTrkPt1000->back()).begin(),
                       std::bind2nd(std::divides<float>(), m_units));
      } else { m_SumPtTrkPt1000->push_back( junkFlt ); }

      if ( trkWidth1000.isAvailable( *jet ) ) {
        m_TrackWidthPt1000->push_back( trkWidth1000( *jet ) );
      } else { m_TrackWidthPt1000->push_back( junkFlt ); }

      if ( nTrk500.isAvailable( *jet ) ) {
        m_NumTrkPt500->push_back( nTrk500( *jet ) );
      } else { m_NumTrkPt500->push_back( junkInt ); }

      if ( sumPt500.isAvailable( *jet ) ) {
        m_SumPtTrkPt500->push_back( sumPt500( *jet ) );
        std::transform((m_SumPtTrkPt500->back()).begin(),
                       (m_SumPtTrkPt500->back()).end(),
                       (m_SumPtTrkPt500->back()).begin(),
                       std::bind2nd(std::divides<float>(), m_units));
      } else { m_SumPtTrkPt500->push_back( junkFlt ); }

      if ( trkWidth500.isAvailable( *jet ) ) {
        m_TrackWidthPt500->push_back( trkWidth500( *jet ) );
      } else { m_TrackWidthPt500->push_back( junkFlt ); }

      if ( jvf.isAvailable( *jet ) ) {
        m_JVF->push_back( jvf( *jet ) );
      } else { m_JVF->push_back( junkFlt ); }

    } // trackAll

    if ( m_infoSwitch.m_trackPV && pvLocation >= 0 ) {

      if ( nTrk1000.isAvailable( *jet ) ) {
        m_NumTrkPt1000PV->push_back( nTrk1000( *jet )[pvLocation] );
      } else { m_NumTrkPt1000PV->push_back( -999 ); }

      if ( sumPt1000.isAvailable( *jet ) ) {
        m_SumPtTrkPt1000PV->push_back( sumPt1000( *jet )[pvLocation] / m_units );
      } else { m_SumPtTrkPt1000PV->push_back( -999 ); }

      if ( trkWidth1000.isAvailable( *jet ) ) {
        m_TrackWidthPt1000PV->push_back( trkWidth1000( *jet )[pvLocation] );
      } else { m_TrackWidthPt1000PV->push_back( -999 ); }

      if ( nTrk500.isAvailable( *jet ) ) {
        m_NumTrkPt500PV->push_back( nTrk500( *jet )[pvLocation] );
      } else { m_NumTrkPt500PV->push_back( -999 ); }

      if ( sumPt500.isAvailable( *jet ) ) {
        m_SumPtTrkPt500PV->push_back( sumPt500( *jet )[pvLocation] / m_units );
      } else { m_SumPtTrkPt500PV->push_back( -999 ); }

      if ( trkWidth500.isAvailable( *jet ) ) {
        m_TrackWidthPt500PV->push_back( trkWidth500( *jet )[pvLocation] );
      } else { m_TrackWidthPt500PV->push_back( -999 ); }

      if ( jvf.isAvailable( *jet ) ) {
        m_JVFPV->push_back( jvf( *jet )[pvLocation] );
      } else { m_JVFPV->push_back( -999 ); }

      static SG::AuxElement::ConstAccessor< float > jvt ("Jvt");
      safeFill<float, float, xAOD::Jet>(jet, jvt, m_Jvt, -999);

      static SG::AuxElement::ConstAccessor< float > jvtJvfcorr ("JvtJvfcorr");
      safeFill<float, float, xAOD::Jet>(jet, jvtJvfcorr, m_JvtJvfcorr, -999);

      static SG::AuxElement::ConstAccessor< float > jvtRpt ("JvtRpt");
      safeFill<float, float, xAOD::Jet>(jet, jvtRpt, m_JvtRpt, -999);

      if ( m_mc ) {
	static SG::AuxElement::ConstAccessor< std::vector< float > > jvtSF_Loose("JetJvtEfficiency_JVTSyst_JVT_Loose");
	static SG::AuxElement::ConstAccessor< std::vector< float > > jvtSF_Medium("JetJvtEfficiency_JVTSyst_JVT_Medium");
	static SG::AuxElement::ConstAccessor< std::vector< float > > jvtSF_Tight("JetJvtEfficiency_JVTSyst_JVT_Tight");
	std::vector<float> junkSF(1,1.0);
	
	if ( jvtSF_Loose.isAvailable( *jet ) )  { m_JvtEff_SF_Loose->push_back( jvtSF_Loose( *jet ) );   } else { m_JvtEff_SF_Loose->push_back( junkSF ); }
	if ( jvtSF_Medium.isAvailable( *jet ) ) { m_JvtEff_SF_Medium->push_back( jvtSF_Medium( *jet ) ); } else { m_JvtEff_SF_Medium->push_back( junkSF ); }
	if ( jvtSF_Tight.isAvailable( *jet ) )  { m_JvtEff_SF_Tight->push_back( jvtSF_Tight( *jet ) );   } else { m_JvtEff_SF_Tight->push_back( junkSF ); }
      }
      //      static SG::AuxElement::ConstAccessor<float> ghostTrackAssFrac("GhostTrackAssociationFraction");
      //      if ( ghostTrackAssFrac.isAvailable( *jet) ) {
      //        m_ghostTrackAssFrac->push_back( ghostTrackAssFrac( *jet) );
      //      } else { m_ghostTrackAssFrac->push_back( -999 ) ; }

    } // trackPV

  }

  if ( m_infoSwitch.m_allTrack ) {
    static SG::AuxElement::ConstAccessor< int > ghostTrackCount("GhostTrackCount");
    safeFill<int, int, xAOD::Jet>(jet, ghostTrackCount, m_GhostTrackCount, -999);

    static SG::AuxElement::ConstAccessor< float > ghostTrackPt ("GhostTrackPt");
    safeFill<float, float, xAOD::Jet>(jet, ghostTrackPt, m_GhostTrackPt, -999, m_units);

    std::vector<float> pt;
    std::vector<float> qOverP;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> e;
    std::vector<float> d0;
    std::vector<float> z0;
    std::vector<int> nPixHits;
    std::vector<int> nSCTHits;
    std::vector<int> nTRTHits;
    std::vector<int> nPixSharedHits;
    std::vector<int> nPixSplitHits;
    std::vector<int> nIMLPixHits;
    std::vector<int> nIMLPixSharedHits;
    std::vector<int> nIMLPixSplitHits;
    std::vector<int> nNIMLPixHits;
    std::vector<int> nNIMLPixSharedHits;
    std::vector<int> nNIMLPixSplitHits;
    static SG::AuxElement::ConstAccessor< std::vector<ElementLink<DataVector<xAOD::IParticle> > > >ghostTrack ("GhostTrack");
    if ( ghostTrack.isAvailable( *jet ) ) {
      std::vector<ElementLink<DataVector<xAOD::IParticle> > > trackLinks = ghostTrack( *jet );
      //std::vector<float> pt(trackLinks.size(),-999);
      for ( auto link_itr : trackLinks ) {
        if( !link_itr.isValid() ) { continue; }
        const xAOD::TrackParticle* track = dynamic_cast<const xAOD::TrackParticle*>( *link_itr );
        // if asking for tracks passing PV selection ( i.e. JVF JVT tracks )
        if( m_infoSwitch.m_allTrackPVSel ) {
          // PV selection from
          // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JvtManualRecalculation
          if( track->pt() < 500 )                { continue; } // pT cut
          if( !m_trkSelTool->accept(*track,pv) ) { continue; } // ID quality cut
          if( track->vertex() != pv ) {                        // if not in PV vertex fit
            if( track->vertex() != 0 )           { continue; } // make sure in no vertex fits
            if( fabs((track->z0()+track->vz()-pv->z())*sin(track->theta())) > 3.0 ) { continue; } // make sure close to PV in z
          }
        }
        pt. push_back( track->pt() / m_units );
        qOverP.push_back( track->qOverP() * m_units );
        eta.push_back( track->eta() );
        phi.push_back( track->phi() );
        e.  push_back( track->e()  / m_units );
        d0. push_back( track->d0() );
        z0. push_back( track->z0() + track->vz() - pv->z() ); // store z0 wrt PV...most useful
        if( m_infoSwitch.m_allTrackDetail ) {
          uint8_t getInt(0);
          // n pix, sct, trt
          track->summaryValue( getInt, xAOD::numberOfPixelHits );
          nPixHits.push_back( getInt );
          track->summaryValue( getInt, xAOD::numberOfSCTHits );
          nSCTHits.push_back( getInt );
          track->summaryValue( getInt, xAOD::numberOfTRTHits );
          nTRTHits.push_back( getInt );
          // pixel split shared
          track->summaryValue( getInt, xAOD::numberOfPixelSharedHits );
          nPixSharedHits.push_back( getInt );
          track->summaryValue( getInt, xAOD::numberOfPixelSplitHits );
          nPixSplitHits.push_back( getInt );
          // n ibl, split, shared
          track->summaryValue( getInt, xAOD::numberOfInnermostPixelLayerHits );
          nIMLPixHits.push_back( getInt );
          track->summaryValue( getInt, xAOD::numberOfInnermostPixelLayerSharedHits );
          nIMLPixSharedHits.push_back( getInt );
          track->summaryValue( getInt, xAOD::numberOfInnermostPixelLayerSplitHits );
          nIMLPixSplitHits.push_back( getInt );
          // n bl,  split, shared
          track->summaryValue( getInt, xAOD::numberOfNextToInnermostPixelLayerHits );
          nNIMLPixHits.push_back( getInt );
          track->summaryValue( getInt, xAOD::numberOfNextToInnermostPixelLayerSharedHits );
          nNIMLPixSharedHits.push_back( getInt );
          track->summaryValue( getInt, xAOD::numberOfNextToInnermostPixelLayerSplitHits );
          nNIMLPixSplitHits.push_back( getInt );
        }
      }
    } // if ghostTrack available
    m_GhostTrack_pt-> push_back( pt  );
    m_GhostTrack_qOverP-> push_back( qOverP );
    m_GhostTrack_eta->push_back( eta );
    m_GhostTrack_phi->push_back( phi );
    m_GhostTrack_e->  push_back( e   );
    m_GhostTrack_d0-> push_back( d0  );
    m_GhostTrack_z0-> push_back( z0  );
    if( m_infoSwitch.m_allTrackDetail ) {
      m_GhostTrack_nPixelHits->push_back( nPixHits );
      m_GhostTrack_nSCTHits->push_back( nSCTHits );
      m_GhostTrack_nTRTHits->push_back( nTRTHits );
      m_GhostTrack_nPixelSharedHits->push_back( nPixSharedHits );
      m_GhostTrack_nPixelSplitHits->push_back( nPixSplitHits );
      m_GhostTrack_nInnermostPixelLayerHits->push_back( nIMLPixHits );
      m_GhostTrack_nInnermostPixelLayerSharedHits->push_back( nIMLPixSharedHits );
      m_GhostTrack_nInnermostPixelLayerSplitHits->push_back( nIMLPixSplitHits );
      m_GhostTrack_nNextToInnermostPixelLayerHits->push_back( nNIMLPixHits );
      m_GhostTrack_nNextToInnermostPixelLayerSharedHits->push_back( nNIMLPixSharedHits );
      m_GhostTrack_nNextToInnermostPixelLayerSplitHits->push_back( nNIMLPixSplitHits );
    }
  } // allTrack switch

  if( m_infoSwitch.m_constituent ) {
    m_numConstituents->push_back( jet->numConstituents() );
  }

  if( m_infoSwitch.m_constituentAll ) {
    m_constituentWeights->push_back( jet->getAttribute< std::vector<float> >( "constituentWeights" ) );
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> e;
    xAOD::JetConstituentVector consVec = jet->getConstituents();
    if( consVec.isValid() ) {
      // don't use auto since iterator can also set the scale ...
      // not sure what that does with auto - probably default but just incase
      // use the example provided in
      // http://acode-browser.usatlas.bnl.gov/lxr/source/atlas/Event/xAOD/xAODJet/xAODJet/JetConstituentVector.h
      xAOD::JetConstituentVector::iterator constit = consVec.begin();
      xAOD::JetConstituentVector::iterator constitE = consVec.end();
      for( ; constit != constitE; constit++){
        pt. push_back( constit->pt() / m_units );
        eta.push_back( constit->eta() );
        phi.push_back( constit->phi() );
        e.  push_back( constit->e() / m_units  );
      }
    }
    m_constituent_pt-> push_back( pt  );
    m_constituent_eta->push_back( eta );
    m_constituent_phi->push_back( phi );
    m_constituent_e->  push_back( e   );
  }

  if ( m_infoSwitch.m_flavTag || m_infoSwitch.m_flavTagHLT ) {
    const xAOD::BTagging * myBTag(0);
    
    if(m_infoSwitch.m_flavTag){
      myBTag = jet->btagging();
    }else if(m_infoSwitch.m_flavTagHLT){
      myBTag = jet->auxdata< const xAOD::BTagging* >("HLTBTag");
    }

    if(m_infoSwitch.m_JVC ) {
      static SG::AuxElement::ConstAccessor<double> JetVertexCharge_discriminant("JetVertexCharge_discriminant");
      safeFill<double, double, xAOD::BTagging>(myBTag, JetVertexCharge_discriminant, m_JetVertexCharge_discriminant, -999);
    }

    //MV2c00 MV2c20 MV2c10 MV2c100 MV2m
    double val(-999);
    myBTag->variable<double>("MV2c00", "discriminant", val);
    m_MV2c00->push_back( val );
    myBTag->variable<double>("MV2c10", "discriminant", val);
    m_MV2c10->push_back( val );
    myBTag->variable<double>("MV2c20", "discriminant", val);
    m_MV2c20->push_back( val );
    myBTag->variable<double>("MV2c100", "discriminant", val);
    m_MV2c100->push_back( val );

    // flavor groups truth definition
    static SG::AuxElement::ConstAccessor<int> hadConeExclTruthLabel("HadronConeExclTruthLabelID");
    safeFill<int, int, xAOD::Jet>(jet, hadConeExclTruthLabel, m_HadronConeExclTruthLabelID, -999);

    if(m_infoSwitch.m_jetFitterDetails ) {

      static SG::AuxElement::ConstAccessor< int   > jf_nVTXAcc       ("JetFitter_nVTX");
      safeFill<int, float, xAOD::BTagging>(myBTag, jf_nVTXAcc, m_JetFitter_nVTX, -999);

      static SG::AuxElement::ConstAccessor< int   > jf_nSingleTracks ("JetFitter_nSingleTracks");
      safeFill<int, float, xAOD::BTagging>(myBTag, jf_nSingleTracks, m_JetFitter_nSingleTracks, -999);

      static SG::AuxElement::ConstAccessor< int   > jf_nTracksAtVtx  ("JetFitter_nTracksAtVtx");
      safeFill<int, float, xAOD::BTagging>(myBTag, jf_nTracksAtVtx, m_JetFitter_nTracksAtVtx, -999);

      static SG::AuxElement::ConstAccessor< float > jf_mass          ("JetFitter_mass");
      safeFill<float, float, xAOD::BTagging>(myBTag, jf_mass, m_JetFitter_mass, -999);

      static SG::AuxElement::ConstAccessor< float > jf_energyFraction("JetFitter_energyFraction");
      safeFill<float, float, xAOD::BTagging>(myBTag, jf_energyFraction, m_JetFitter_energyFraction, -999);

      static SG::AuxElement::ConstAccessor< float > jf_significance3d("JetFitter_significance3d");
      safeFill<float, float, xAOD::BTagging>(myBTag, jf_significance3d, m_JetFitter_significance3d, -999);

      static SG::AuxElement::ConstAccessor< float > jf_deltaeta      ("JetFitter_deltaeta");
      safeFill<float, float, xAOD::BTagging>(myBTag, jf_deltaeta, m_JetFitter_deltaeta, -999);

      static SG::AuxElement::ConstAccessor< float > jf_deltaphi      ("JetFitter_deltaphi");
      safeFill<float, float, xAOD::BTagging>(myBTag, jf_deltaphi, m_JetFitter_deltaphi, -999);

      static SG::AuxElement::ConstAccessor< int   > jf_N2Tpar        ("JetFitter_N2Tpair");
      safeFill<int, float, xAOD::BTagging>(myBTag, jf_N2Tpar, m_JetFitter_N2Tpar, -999);

      //static SG::AuxElement::ConstAccessor< double > jf_pb           ("JetFitterCombNN_pb");
      //safeFill<double, float, xAOD::BTagging>(myBTag, jf_pb, m_JetFitter_pb, -999);
      //
      //static SG::AuxElement::ConstAccessor< double > jf_pc           ("JetFitterCombNN_pc");
      //safeFill<double, float, xAOD::BTagging>(myBTag, jf_pc, m_JetFitter_pc, -999);
      //
      //static SG::AuxElement::ConstAccessor< double > jf_pu           ("JetFitterCombNN_pu");
      //safeFill<double, float, xAOD::BTagging>(myBTag, jf_pu, m_JetFitter_pu, -999);

    }

    if(m_infoSwitch.m_svDetails ) {
      if(m_debug) cout << "Filling m_svDetails " << endl;

      /// @brief SV0 : Number of good tracks in vertex
      static SG::AuxElement::ConstAccessor< int   >   sv0_NGTinSvxAcc     ("SV0_NGTinSvx");
      safeFill<int, float, xAOD::BTagging>(myBTag,    sv0_NGTinSvxAcc, m_sv0_NGTinSvx, -999);

      // @brief SV0 : Number of 2-track pairs
      static SG::AuxElement::ConstAccessor< int   >   sv0_N2TpairAcc      ("SV0_N2Tpair");
      safeFill<int, float, xAOD::BTagging>(myBTag, sv0_N2TpairAcc, m_sv0_N2Tpair, -999);

      /// @brief SV0 : vertex mass
      static SG::AuxElement::ConstAccessor< float   > sv0_masssvxAcc      ("SV0_masssvx");
      safeFill<float, float, xAOD::BTagging>(myBTag, sv0_masssvxAcc, m_sv0_massvx, -999);

      /// @brief SV0 : energy fraction
      static SG::AuxElement::ConstAccessor< float   > sv0_efracsvxAcc     ("SV0_efracsvx");                                    
      safeFill<float, float, xAOD::BTagging>(myBTag, sv0_efracsvxAcc, m_sv0_efracsvx, -999);
      
      /// @brief SV0 : 3D vertex significance
      static SG::AuxElement::ConstAccessor< float   > sv0_normdistAcc     ("SV0_normdist");
      safeFill<float, float, xAOD::BTagging>(myBTag, sv0_normdistAcc, m_sv0_normdist, -999);

      double sv0;
      myBTag->variable<double>("SV0", "significance3D", sv0);
      m_SV0->push_back(sv0);

      m_SV1IP3D->push_back( myBTag -> SV1plusIP3D_discriminant() );


      /// @brief SV1 : Number of good tracks in vertex
      static SG::AuxElement::ConstAccessor< int   >   sv1_NGTinSvxAcc     ("SV1_NGTinSvx");
      safeFill<int, float, xAOD::BTagging>(myBTag, sv1_NGTinSvxAcc, m_sv1_NGTinSvx, -999);

      // @brief SV1 : Number of 2-track pairs
      static SG::AuxElement::ConstAccessor< int   >   sv1_N2TpairAcc      ("SV1_N2Tpair");
      safeFill<int, float, xAOD::BTagging>(myBTag, sv1_N2TpairAcc, m_sv1_N2Tpair, -999);

      /// @brief SV1 : vertex mass
      static SG::AuxElement::ConstAccessor< float   > sv1_masssvxAcc      ("SV1_masssvx");
      safeFill<float, float, xAOD::BTagging>(myBTag, sv1_masssvxAcc, m_sv1_massvx, -999);

      /// @brief SV1 : energy fraction
      static SG::AuxElement::ConstAccessor< float   > sv1_efracsvxAcc     ("SV1_efracsvx");
      safeFill<float, float, xAOD::BTagging>(myBTag, sv1_efracsvxAcc, m_sv1_efracsvx, -999);

      /// @brief SV1 : 3D vertex significance
      static SG::AuxElement::ConstAccessor< float   > sv1_normdistAcc     ("SV1_normdist");
      safeFill<float, float, xAOD::BTagging>(myBTag, sv1_normdistAcc, m_sv1_normdist, -999);

      double sv1_pu = -30;  myBTag->variable<double>("SV1", "pu", sv1_pu);
      double sv1_pb = -30;  myBTag->variable<double>("SV1", "pb", sv1_pb);
      double sv1_pc = -30;  myBTag->variable<double>("SV1", "pc", sv1_pc);

      m_sv1_pu         ->push_back(sv1_pu);
      m_sv1_pb         ->push_back(sv1_pb);
      m_sv1_pc         ->push_back(sv1_pc);
      m_SV1            ->push_back( myBTag->calcLLR(sv1_pb,sv1_pu)  );
      m_sv1_c          ->push_back( myBTag->calcLLR(sv1_pb,sv1_pc)  );
      m_sv1_cu         ->push_back( myBTag->calcLLR(sv1_pc,sv1_pu)  );

      float sv1_Lxy;        myBTag->variable<float>("SV1", "Lxy"         , sv1_Lxy);
      float sv1_L3d;        myBTag->variable<float>("SV1", "L3d"         , sv1_L3d);
      float sv1_distmatlay; myBTag->variable<float>("SV1", "dstToMatLay" , sv1_distmatlay);
      float sv1_dR;         myBTag->variable<float>("SV1", "deltaR"      , sv1_dR );

      m_sv1_Lxy        ->push_back(sv1_Lxy        );
      m_sv1_L3d        ->push_back(sv1_L3d        );
      m_sv1_distmatlay ->push_back(sv1_distmatlay );
      m_sv1_dR         ->push_back(sv1_dR         );
  
      
    }

    if(m_infoSwitch.m_ipDetails ) {
      if(m_debug) cout << "Filling m_ipDetails " << endl;

      //
      // IP2D
      //

      /// @brief IP2D: track grade
      static SG::AuxElement::ConstAccessor< vector<int>   >   IP2D_gradeOfTracksAcc     ("IP2D_gradeOfTracks");
      safeVecFill<int, float, xAOD::BTagging>(myBTag, IP2D_gradeOfTracksAcc, m_IP2D_gradeOfTracks);

      /// @brief IP2D : tracks from V0
      static SG::AuxElement::ConstAccessor< vector<bool>   >  IP2D_flagFromV0ofTracksAcc("IP2D_flagFromV0ofTracks");
      safeVecFill<bool, float, xAOD::BTagging>(myBTag, IP2D_flagFromV0ofTracksAcc, m_IP2D_flagFromV0ofTracks);

      /// @brief IP2D : d0 value with respect to primary vertex
      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_valD0wrtPVofTracksAcc("IP2D_valD0wrtPVofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP2D_valD0wrtPVofTracksAcc, m_IP2D_valD0wrtPVofTracks);

      /// @brief IP2D : d0 significance with respect to primary vertex
      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_sigD0wrtPVofTracksAcc("IP2D_sigD0wrtPVofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP2D_sigD0wrtPVofTracksAcc, m_IP2D_sigD0wrtPVofTracks);

      /// @brief IP2D : track contribution to B likelihood
      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_weightBofTracksAcc   ("IP2D_weightBofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP2D_weightBofTracksAcc, m_IP2D_weightBofTracks);

      /// @brief IP2D : track contribution to C likelihood
      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_weightCofTracksAcc   ("IP2D_weightCofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP2D_weightCofTracksAcc, m_IP2D_weightCofTracks);

      /// @brief IP2D : track contribution to U likelihood
      static SG::AuxElement::ConstAccessor< vector<float>   > IP2D_weightUofTracksAcc   ("IP2D_weightUofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP2D_weightUofTracksAcc, m_IP2D_weightUofTracks);

      double ip2_pu = -99;  myBTag->variable<double>("IP2D", "pu", ip2_pu);
      double ip2_pb = -99;  myBTag->variable<double>("IP2D", "pb", ip2_pb);
      double ip2_pc = -99;  myBTag->variable<double>("IP2D", "pc", ip2_pc);

      m_IP2D_pu         ->push_back(ip2_pu);
      m_IP2D_pb         ->push_back(ip2_pb);
      m_IP2D_pc         ->push_back(ip2_pc);

      m_IP2D            ->push_back( myBTag->calcLLR(ip2_pb,ip2_pu)  );
      m_IP2D_c          ->push_back( myBTag->calcLLR(ip2_pb,ip2_pc)  );
      m_IP2D_cu         ->push_back( myBTag->calcLLR(ip2_pc,ip2_pu)  );


      //
      // IP3D
      //

      /// @brief IP3D: track grade
      static SG::AuxElement::ConstAccessor< vector<int>   >   IP3D_gradeOfTracksAcc     ("IP3D_gradeOfTracks");
      safeVecFill<int, float, xAOD::BTagging>(myBTag, IP3D_gradeOfTracksAcc, m_IP3D_gradeOfTracks);

      /// @brief IP3D : tracks from V0
      static SG::AuxElement::ConstAccessor< vector<bool>   >  IP3D_flagFromV0ofTracksAcc("IP3D_flagFromV0ofTracks");
      safeVecFill<bool, float, xAOD::BTagging>(myBTag, IP3D_flagFromV0ofTracksAcc, m_IP3D_flagFromV0ofTracks);

      /// @brief IP3D : d0 value with respect to primary vertex
      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_valD0wrtPVofTracksAcc("IP3D_valD0wrtPVofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP3D_valD0wrtPVofTracksAcc, m_IP3D_valD0wrtPVofTracks);

      /// @brief IP3D : d0 significance with respect to primary vertex
      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_sigD0wrtPVofTracksAcc("IP3D_sigD0wrtPVofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP3D_sigD0wrtPVofTracksAcc, m_IP3D_sigD0wrtPVofTracks);

      /// @brief IP3D : z0 value with respect to primary vertex
      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_valZ0wrtPVofTracksAcc("IP3D_valZ0wrtPVofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP3D_valZ0wrtPVofTracksAcc, m_IP3D_valZ0wrtPVofTracks);

      /// @brief IP3D : z0 significance with respect to primary vertex
      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_sigZ0wrtPVofTracksAcc("IP3D_sigZ0wrtPVofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP3D_sigZ0wrtPVofTracksAcc, m_IP3D_sigZ0wrtPVofTracks);

      /// @brief IP3D : track contribution to B likelihood
      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_weightBofTracksAcc   ("IP3D_weightBofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP3D_weightBofTracksAcc, m_IP3D_weightBofTracks);

      /// @brief IP3D : track contribution to C likelihood
      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_weightCofTracksAcc   ("IP3D_weightCofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP3D_weightCofTracksAcc, m_IP3D_weightCofTracks);

      /// @brief IP3D : track contribution to U likelihood
      static SG::AuxElement::ConstAccessor< vector<float>   > IP3D_weightUofTracksAcc   ("IP3D_weightUofTracks");
      safeVecFill<float, float, xAOD::BTagging>(myBTag, IP3D_weightUofTracksAcc, m_IP3D_weightUofTracks);

      double ip3_pu = -30;  myBTag->variable<double>("IP3D", "pu", ip3_pu);
      double ip3_pb = -30;  myBTag->variable<double>("IP3D", "pb", ip3_pb);
      double ip3_pc = -30;  myBTag->variable<double>("IP3D", "pc", ip3_pc);

      m_IP3D->push_back(    myBTag -> IP3D_loglikelihoodratio()  );

      m_IP3D_pu         ->push_back(ip3_pu  );
      m_IP3D_pb         ->push_back(ip3_pb  );
      m_IP3D_pc         ->push_back(ip3_pc  );

      m_IP3D            ->push_back( myBTag->calcLLR(ip3_pb,ip3_pu)  );
      m_IP3D_c          ->push_back( myBTag->calcLLR(ip3_pb,ip3_pc)  );
      m_IP3D_cu         ->push_back( myBTag->calcLLR(ip3_pc,ip3_pu)  );

    }



    if(m_infoSwitch.m_flavTagHLT ) {
      if(m_debug) cout << "Filling m_flavTagHLT " << endl;
      const xAOD::Vertex *online_pvx       = jet->auxdata<const xAOD::Vertex*>("HLTBJetTracks_vtx");
      const xAOD::Vertex *online_pvx_bkg   = jet->auxdata<const xAOD::Vertex*>("HLTBJetTracks_vtx_bkg");
      const xAOD::Vertex *offline_pvx      = jet->auxdata<const xAOD::Vertex*>("offline_vtx");      

      if(online_pvx)  m_vtxOnlineValid->push_back(1.0);
      else            m_vtxOnlineValid->push_back(0.0);
      
      char hadDummyPV = jet->auxdata< char >("hadDummyPV");
      if( hadDummyPV == '0')  m_vtxHadDummy->push_back(0.0);
      if( hadDummyPV == '1')  m_vtxHadDummy->push_back(1.0);
      if( hadDummyPV == '2')  m_vtxHadDummy->push_back(2.0);

      static SG::AuxElement::ConstAccessor< float > acc_bs_online_vs ("bs_online_vz");
      if(acc_bs_online_vs.isAvailable( *jet) ){
	if(m_debug) cout << "Have bs_online_vz " << endl;
	float bs_online_vz = jet->auxdata< float >("bs_online_vz");
	//std::cout << "**bs_online_vz " << bs_online_vz << std::endl;
	m_bs_online_vz->push_back( bs_online_vz );

	float bs_online_vx = jet->auxdata< float >("bs_online_vx");
	//std::cout << "**bs_online_vx " << bs_online_vx << std::endl;
	m_bs_online_vx->push_back( bs_online_vx );

	float bs_online_vy = jet->auxdata< float >("bs_online_vy");
	//std::cout << "**bs_online_vy " << bs_online_vy << std::endl;
	m_bs_online_vy->push_back( bs_online_vy );
      }else{
	m_bs_online_vz->push_back( -999 );
	m_bs_online_vx->push_back( -999 );
	m_bs_online_vy->push_back( -999 );
      }

      if(m_debug) cout << "Filling m_vtx_offline " << endl;
      m_vtx_offline_x0->push_back( offline_pvx->x() );
      m_vtx_offline_y0->push_back( offline_pvx->y() );
      m_vtx_offline_z0->push_back( offline_pvx->z() );
      if(m_debug) cout << "Done Filling m_vtx_offline " << endl;

      if(m_debug) cout << "Filling m_vtx_online... " << endl;
      if(online_pvx){
	if(m_debug) cout << " ... online_pvx valid " << endl;
        m_vtx_online_x0->push_back( online_pvx->x() );
        m_vtx_online_y0->push_back( online_pvx->y() );
        m_vtx_online_z0->push_back( online_pvx->z() );
      }else{           
        m_vtx_online_x0->push_back( -999 );
        m_vtx_online_y0->push_back( -999 );
        m_vtx_online_z0->push_back( -999 );
      }

      if(m_debug) cout << "Filling m_vtx_online... " << endl;
      if(online_pvx_bkg){
	if(m_debug) cout << " ...online_pvx_bkg valid " << endl;
        m_vtx_online_bkg_x0->push_back( online_pvx_bkg->x() );
        m_vtx_online_bkg_y0->push_back( online_pvx_bkg->y() );
        m_vtx_online_bkg_z0->push_back( online_pvx_bkg->z() );
      }else{           
        m_vtx_online_bkg_x0->push_back( -999 );
        m_vtx_online_bkg_y0->push_back( -999 );
        m_vtx_online_bkg_z0->push_back( -999 );
      }

    }// m_flavTagHLT
    if(m_debug) cout << "Done m_flavTagHLT " << endl;
  }


  if( !m_infoSwitch.m_sfFTagFix.empty() ) {
    for( unsigned int i=0; i<m_infoSwitch.m_sfFTagFix.size(); i++ ) {
      switch( m_infoSwitch.m_sfFTagFix.at(i) ) {
      case 30 : m_btag_Fix30->Fill( jet ); break;
      case 50 : m_btag_Fix50->Fill( jet ); break;
      case 60 : m_btag_Fix60->Fill( jet ); break;
      case 70 : m_btag_Fix70->Fill( jet ); break;
      case 77 : m_btag_Fix77->Fill( jet ); break;
      case 80 : m_btag_Fix80->Fill( jet ); break;
      case 85 : m_btag_Fix85->Fill( jet ); break;
      case 90 : m_btag_Fix90->Fill( jet ); break;
      }
    }
  } // sfFTagFix



  if( !m_infoSwitch.m_sfFTagFlt.empty() ) {
    for( unsigned int i=0; i<m_infoSwitch.m_sfFTagFlt.size(); i++ ) {
      switch( m_infoSwitch.m_sfFTagFlt.at(i) ) {
      case 30 : m_btag_Flt30->Fill( jet );  break;
      case 50 : m_btag_Flt50->Fill( jet );	break;
      case 60 : m_btag_Flt60->Fill( jet );	break;
      case 70 : m_btag_Flt70->Fill( jet );	break;
      case 77 : m_btag_Flt77->Fill( jet );	break;
      case 85 : m_btag_Flt85->Fill( jet );  break;
      }
    }
  } // sfFTagFlt


  if ( m_infoSwitch.m_area ) {

    static SG::AuxElement::ConstAccessor<float> ghostArea("JetGhostArea");
    safeFill<float, float, xAOD::Jet>(jet, ghostArea, m_GhostArea, -999);

    static SG::AuxElement::ConstAccessor<float> activeArea("ActiveArea");
    safeFill<float, float, xAOD::Jet>(jet, activeArea, m_ActiveArea, -999);

    static SG::AuxElement::ConstAccessor<float> voronoiArea("VoronoiArea");
    safeFill<float, float, xAOD::Jet>(jet, voronoiArea, m_VoronoiArea, -999);

    static SG::AuxElement::ConstAccessor<float> activeArea_pt("ActiveArea4vec_pt");
    safeFill<float, float, xAOD::Jet>(jet, activeArea_pt, m_ActiveArea4vec_pt, -999);

    static SG::AuxElement::ConstAccessor<float> activeArea_eta("ActiveArea4vec_eta");
    safeFill<float, float, xAOD::Jet>(jet, activeArea_eta, m_ActiveArea4vec_eta, -999);

    static SG::AuxElement::ConstAccessor<float> activeArea_phi("ActiveArea4vec_phi");
    safeFill<float, float, xAOD::Jet>(jet, activeArea_phi, m_ActiveArea4vec_phi, -999);

    static SG::AuxElement::ConstAccessor<float> activeArea_m("ActiveArea4vec_m");
    safeFill<float, float, xAOD::Jet>(jet, activeArea_m, m_ActiveArea4vec_m, -999);
  }


  if ( m_infoSwitch.m_truth && m_mc ) {

    static SG::AuxElement::ConstAccessor<int> ConeTruthLabelID ("ConeTruthLabelID");
    safeFill<int, int, xAOD::Jet>(jet, ConeTruthLabelID, m_ConeTruthLabelID, -999);

    static SG::AuxElement::ConstAccessor<int> TruthCount ("TruthCount");
    safeFill<int, int, xAOD::Jet>(jet, TruthCount, m_TruthCount, -999);

    //    seems to be empty
    //      static SG::AuxElement::ConstAccessor<float> TruthPt ("TruthPt");
    //      if ( TruthPt.isAvailable( *jet) ) {
    //        m_truthPt->push_back( TruthPt( *jet)/1000 );
    //      } else { m_truthPt->push_back( -999 ); }

    static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_B ("TruthLabelDeltaR_B");
    safeFill<float, float, xAOD::Jet>(jet, TruthLabelDeltaR_B, m_TruthLabelDeltaR_B, -999);

    static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_C ("TruthLabelDeltaR_C");
    safeFill<float, float, xAOD::Jet>(jet, TruthLabelDeltaR_C, m_TruthLabelDeltaR_C, -999);

    static SG::AuxElement::ConstAccessor<float> TruthLabelDeltaR_T ("TruthLabelDeltaR_T");
    safeFill<float, float, xAOD::Jet>(jet, TruthLabelDeltaR_T, m_TruthLabelDeltaR_T, -999);

    static SG::AuxElement::ConstAccessor<int> partonLabel("PartonTruthLabelID");
    safeFill<int, int, xAOD::Jet>(jet, partonLabel, m_PartonTruthLabelID, -999);

    static SG::AuxElement::ConstAccessor<float> ghostTruthAssFrac("GhostTruthAssociationFraction");
    safeFill<float, float, xAOD::Jet>(jet, ghostTruthAssFrac, m_GhostTruthAssociationFraction, -999);

    const xAOD::Jet* truthJet = HelperFunctions::getLink<xAOD::Jet>( jet, "GhostTruthAssociationLink" );
    if(truthJet) {
      m_truth_pt->push_back ( truthJet->pt() / m_units );
      m_truth_eta->push_back( truthJet->eta() );
      m_truth_phi->push_back( truthJet->phi() );
      m_truth_E->push_back  ( truthJet->e() / m_units );
    } else {
      m_truth_pt->push_back ( -999 );
      m_truth_eta->push_back( -999 );
      m_truth_phi->push_back( -999 );
      m_truth_E->push_back  ( -999 );
    }

  }

  if ( m_infoSwitch.m_truthDetails ) {

    //
    // B-Hadron Details
    //
    static SG::AuxElement::ConstAccessor<int> GhostBHadronsFinalCount ("GhostBHadronsFinalCount");
    safeFill<int, int, xAOD::Jet>(jet, GhostBHadronsFinalCount, m_GhostBHadronsFinalCount, -999);

    static SG::AuxElement::ConstAccessor<int> GhostBHadronsInitialCount ("GhostBHadronsInitialCount");
    safeFill<int, int, xAOD::Jet>(jet, GhostBHadronsInitialCount, m_GhostBHadronsInitialCount, -999);

    static SG::AuxElement::ConstAccessor<int> GhostBQuarksFinalCount ("GhostBQuarksFinalCount");
    safeFill<int, int, xAOD::Jet>(jet, GhostBQuarksFinalCount, m_GhostBQuarksFinalCount, -999);

    static SG::AuxElement::ConstAccessor<float> GhostBHadronsFinalPt ("GhostBHadronsFinalPt");
    safeFill<float, float, xAOD::Jet>(jet, GhostBHadronsFinalPt, m_GhostBHadronsFinalPt, -999);

    static SG::AuxElement::ConstAccessor<float> GhostBHadronsInitialPt ("GhostBHadronsInitialPt");
    safeFill<float, float, xAOD::Jet>(jet, GhostBHadronsInitialPt, m_GhostBHadronsInitialPt, -999);

    static SG::AuxElement::ConstAccessor<float> GhostBQuarksFinalPt ("GhostBQuarksFinalPt");
    safeFill<float, float, xAOD::Jet>(jet, GhostBQuarksFinalPt, m_GhostBQuarksFinalPt, -999);

    //
    // C-Hadron Details
    //
    static SG::AuxElement::ConstAccessor<int> GhostCHadronsFinalCount ("GhostCHadronsFinalCount");
    safeFill<int, int, xAOD::Jet>(jet, GhostCHadronsFinalCount, m_GhostCHadronsFinalCount, -999);

    static SG::AuxElement::ConstAccessor<int> GhostCHadronsInitialCount ("GhostCHadronsInitialCount");
    safeFill<int, int, xAOD::Jet>(jet, GhostCHadronsInitialCount, m_GhostCHadronsInitialCount, -999);

    static SG::AuxElement::ConstAccessor<int> GhostCQuarksFinalCount ("GhostCQuarksFinalCount");
    safeFill<int, int, xAOD::Jet>(jet, GhostCQuarksFinalCount, m_GhostCQuarksFinalCount, -999);

    static SG::AuxElement::ConstAccessor<float> GhostCHadronsFinalPt ("GhostCHadronsFinalPt");
    safeFill<float, float, xAOD::Jet>(jet, GhostCHadronsFinalPt, m_GhostCHadronsFinalPt, -999);

    static SG::AuxElement::ConstAccessor<float> GhostCHadronsInitialPt ("GhostCHadronsInitialPt");
    safeFill<float, float, xAOD::Jet>(jet, GhostCHadronsInitialPt, m_GhostCHadronsInitialPt, -999);

    static SG::AuxElement::ConstAccessor<float> GhostCQuarksFinalPt ("GhostCQuarksFinalPt");
    safeFill<float, float, xAOD::Jet>(jet, GhostCQuarksFinalPt, m_GhostCQuarksFinalPt, -999);

    //
    // Tau Details
    //
    static SG::AuxElement::ConstAccessor<int> GhostTausFinalCount ("GhostTausFinalCount");
    safeFill<int, int, xAOD::Jet>(jet, GhostTausFinalCount, m_GhostTausFinalCount, -999);

    // THE ONLY UN-OFFICIAL PIECE OF CODE HERE USE WITH CAUTION
    static SG::AuxElement::ConstAccessor<float> GhostTausFinalPt ("GhostTausFinalPt");
    safeFill<float, float, xAOD::Jet>(jet, GhostTausFinalPt, m_GhostTausFinalPt, -999);

    // light quark(1,2,3) , gluon (21 or 9), charm(4) and b(5)
    // GhostPartons should select for these pdgIds only
    //    static SG::AuxElement::ConstAccessor< std::vector<const xAOD::TruthParticle*> > ghostPartons("GhostPartons");
    //    if( ghostPartons.isAvailable( *jet )) {
    //    std::vector<const xAOD::TruthParticle*> truthPartons = ghostPartons( *jet );

    std::vector<const xAOD::TruthParticle*> truthPartons = jet->getAssociatedObjects<xAOD::TruthParticle>("GhostPartons");

    if( truthPartons.size() == 0){
      m_truth_pdgId->push_back(-999);
    } else {
      int iParent = 0;
      for(unsigned int i=1; i < truthPartons.size(); ++i){
        if( (truthPartons.at(i)->pt() > 0.001) && (truthPartons.at(i)->e() > truthPartons.at(iParent)->e()) )
          iParent = i;
      }
      m_truth_pdgId->push_back(truthPartons.at(iParent)->pdgId());
      m_truth_partonPt->push_back(truthPartons.at(iParent)->pt() / m_units);
      m_truth_partonDR->push_back(truthPartons.at(iParent)->p4().DeltaR( jet->p4() ));
    }

  }


  if ( m_infoSwitch.m_charge ) {
    xAOD::JetFourMom_t p4UsedInJetCharge;
    bool status = jet->getAttribute<xAOD::JetFourMom_t>( "JetPileupScaleMomentum", p4UsedInJetCharge );
    static SG::AuxElement::ConstAccessor<float>              uncalibratedJetCharge ("Charge");

    if(status){
      float ptUsedInJetCharge   = p4UsedInJetCharge.Pt();
      float calibratedJetCharge = jet->pt() ? (ptUsedInJetCharge * uncalibratedJetCharge(*jet) / jet->pt()) : -99;
      m_charge->push_back(calibratedJetCharge);
    }else{
      m_charge->push_back(-99);
    }

  }

  return;
}


void JetHelpTree::FillGlobalBTagSF( const xAOD::EventInfo* eventInfo ){

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

