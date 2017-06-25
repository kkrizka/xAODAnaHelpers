#ifndef xAODAnaHelpers_Jet_H
#define xAODAnaHelpers_Jet_H

#include "xAODAnaHelpers/Particle.h"
#include "xAODAnaHelpers/TrackParticle.h"
#include "xAODAnaHelpers/JetConstituent.h"
#include "xAODAnaHelpers/MuonContainer.h"


namespace xAH {

  class Jet : public Particle
  {
    ClassDef(Jet, 1);

  public:
      
    Jet();
    virtual ~Jet();

    float rapidity;

    // clean
    float Timing;
    float LArQuality;
    float HECQuality;
    float NegativeE;
    float AverageLArQF;
    float BchCorrCell;
    float N90Constituents;
    float LArBadHVEnergyFrac;
    int   LArBadHVNCell;
    float OotFracClusters5;
    float OotFracClusters10;
    float LeadingClusterPt;
    float LeadingClusterSecondLambda;
    float LeadingClusterCenterLambda;
    float LeadingClusterSecondR;
    int   clean_passLooseBad;
    int   clean_passLooseBadUgly;
    int   clean_passTightBad;
    int   clean_passTightBadUgly;

    // energy
    float HECFrac;
    float EMFrac;
    float CentroidR;
    float FracSamplingMax;
    float FracSamplingMaxIndex;
    float LowEtConstituentsFrac;
    float GhostMuonSegmentCount;
    float Width;

    // scales
    TLorentzVector JetEMScaleMomentum;
    TLorentzVector JetConstitScaleMomentum;
    TLorentzVector JetPileupScaleMomentum;
    TLorentzVector JetOriginConstitScaleMomentum;
    TLorentzVector JetEtaJESScaleMomentum;
    TLorentzVector JetGSCScaleMomentum;
    TLorentzVector JetInsituScaleMomentum;

    // layers
    std::vector<float> EnergyPerSampling;

    // trackPV
    float NumTrkPt1000PV;
    float SumPtTrkPt1000PV;
    float TrackWidthPt1000PV;
    float NumTrkPt500PV;
    float SumPtTrkPt500PV;
    float TrackWidthPt500PV;
    float JVFPV;

    // trackAll
    std::vector<int>   NumTrkPt1000;
    std::vector<float> SumPtTrkPt1000;
    std::vector<float> TrackWidthPt1000;
    std::vector<int>   NumTrkPt500;
    std::vector<float> SumPtTrkPt500;
    std::vector<float> TrackWidthPt500;
    std::vector<float> JVF;

    int GhostTrackCount;
    float GhostTrackPt;
    std::vector<xAH::TrackParticle> GhostTrack;

    // trackAll or trackPV
    float Jvt;
    float JvtJvfcorr;
    float JvtRpt;
    std::vector< float > JetJvtEfficiency_JVTSyst_JVT_Loose;
    std::vector< float > JetJvtEfficiency_JVTSyst_JVT_Medium;
    std::vector< float > JetJvtEfficiency_JVTSyst_JVT_Tight;

    // constituent
    int numConstituents;

    // constituentAll
    std::vector<float> constituentWeights;
    std::vector<JetConstituent> constituents;

    //JVC
    float JetVertexCharge_discriminant;

    // flavTag
    double SV0;
    float SV1;
    float IP3D;
    float MV1;
    double MV2c00;
    double MV2c10;
    double MV2c20;
    double MV2c100;
    float MV2;
    float SV1plusIP3D_discriminant;
    float IP3D_loglikelihoodratio;
    int  HadronConeExclTruthLabelID;

    float vtxOnlineValid;
    float vtxHadDummy;
  
    float bs_online_vx;
    float bs_online_vy;
    float bs_online_vz;

    float vtx_offline_x0;
    float vtx_offline_y0;
    float vtx_offline_z0;

    float vtx_online_x0;
    float vtx_online_y0;
    float vtx_online_z0;

    float vtx_online_bkg_x0;
    float vtx_online_bkg_y0;
    float vtx_online_bkg_z0;

    float JetFitter_nVTX           ;
    float JetFitter_nSingleTracks  ;
    float JetFitter_nTracksAtVtx   ;
    float JetFitter_mass           ;
    float JetFitter_energyFraction ;
    float JetFitter_significance3d ;
    float JetFitter_deltaeta       ;
    float JetFitter_deltaphi       ;
    float JetFitter_N2Tpar         ;

    float SV0_NGTinSvx  ;
    float SV0_N2Tpair   ;
    float SV0_masssvx   ;
    float SV0_efracsvx  ;
    float SV0_normdist  ;
    double SV1_pu       ;
    double SV1_pb       ;
    double SV1_pc       ;
    float SV1_c         ;
    float SV1_cu        ;
    float SV1_NGTinSvx  ;
    float SV1_N2Tpair   ;
    float SV1_masssvx   ;
    float SV1_efracsvx  ;
    float SV1_normdist  ;
    float SV1_Lxy       ;
    float SV1_L3d       ;
    float SV1_distmatlay;
    float SV1_dR        ;

    double IP2D_pu     ;
    double IP2D_pb     ;
    double IP2D_pc     ;
    float IP2D        ;
    float IP2D_c      ;
    float IP2D_cu     ;
    float nIP2DTracks ;

    std::vector<int  > IP2D_gradeOfTracks         ;
    std::vector<bool > IP2D_flagFromV0ofTracks    ;
    std::vector<float> IP2D_valD0wrtPVofTracks    ;
    std::vector<float> IP2D_sigD0wrtPVofTracks    ;
    std::vector<float> IP2D_weightBofTracks       ;
    std::vector<float> IP2D_weightCofTracks       ;
    std::vector<float> IP2D_weightUofTracks       ;

    double IP3D_pu     ;
    double IP3D_pb     ;
    double IP3D_pc     ;
    float IP3D_c      ;
    float IP3D_cu     ;
    float nIP3DTracks ;

    std::vector<int  > IP3D_gradeOfTracks      ;
    std::vector<bool > IP3D_flagFromV0ofTracks ;
    std::vector<float> IP3D_valD0wrtPVofTracks ;
    std::vector<float> IP3D_sigD0wrtPVofTracks ;
    std::vector<float> IP3D_valZ0wrtPVofTracks ;
    std::vector<float> IP3D_sigZ0wrtPVofTracks ;
    std::vector<float> IP3D_weightBofTracks    ;
    std::vector<float> IP3D_weightCofTracks    ;
    std::vector<float> IP3D_weightUofTracks    ;

    char isFix30;
    std::vector<float> sfFix30;

    char isFix50;
    std::vector<float> sfFix50;

    char isFix60;
    std::vector<float> sfFix60;

    char isFix70;
    std::vector<float> sfFix70;

    char isFix77;
    std::vector<float> sfFix77;

    char isFix80;
    std::vector<float> sfFix80;

    char isFix85;
    std::vector<float> sfFix85;

    char isFix90;
    std::vector<float> sfFix90;

    char isFlt30;
    std::vector<float> sfFlt30;

    char isFlt40;
    std::vector<float> sfFlt40;

    char isFlt50;
    std::vector<float> sfFlt50;

    char isFlt60;
    std::vector<float> sfFlt60;

    char isFlt70;
    std::vector<float> sfFlt70;

    char isFlt77;
    std::vector<float> sfFlt77;

    char isFlt85;
    std::vector<float> sfFlt85;

    // area
    float JetGhostArea;
    float ActiveArea;
    float VoronoiArea;
    float ActiveArea4vec_pt;
    float ActiveArea4vec_eta;
    float ActiveArea4vec_phi;
    float ActiveArea4vec_m;

    // truth
    int   ConeTruthLabelID;
    int   TruthCount;
    float TruthLabelDeltaR_B;
    float TruthLabelDeltaR_C;
    float TruthLabelDeltaR_T;
    int   PartonTruthLabelID;
    float GhostTruthAssociationFraction;

    int GhostBHadronsFinalCount;
    int GhostBHadronsInitialCount;
    int GhostBQuarksFinalCount;
    float GhostBHadronsFinalPt;
    float GhostBHadronsInitialPt;
    float GhostBQuarksFinalPt;

    int GhostCHadronsFinalCount;
    int GhostCHadronsInitialCount;
    int GhostCQuarksFinalCount;
    float GhostCHadronsFinalPt;
    float GhostCHadronsInitialPt;
    float GhostCQuarksFinalPt;

    int GhostTausFinalCount;
    float GhostTausFinalPt;

    TLorentzVector truth_p4;
    int truth_pdgId;
    double truth_partonPt;
    double truth_partonDR;

    // charge
    double charge;

    const Muon* matchedMuon;
    const Jet * matchedJet;

  public:

    void muonInJetCorrection(const xAH::MuonContainer* muons);

  };

} //xAH
#endif // xAODAnaHelpers_Jet_H
