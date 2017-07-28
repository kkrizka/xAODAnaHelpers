/********************************************************
 * HelpTreeBase:
 *
 * This class is meant to help the user write out a tree.
 * Some branches are included by default while others
 * need to be added by the user
 *
 * John Alison (john.alison@cern.ch)
 * Gabriel Facini (gabriel.facini@cern.ch)
 * Marco Milesi (marco.milesi@cern.ch)
 * Jeff Dandoy (jeff.dandoy@cern.ch)
 *
 ********************************************************/

// Dear emacs, this is -*-c++-*-
#ifndef xAODAnaHelpers_HelpTreeBase_H
#define xAODAnaHelpers_HelpTreeBase_H

#include <xAODRootAccess/TEvent.h>
#include <xAODRootAccess/TStore.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODEgamma/ElectronContainer.h>
#include <xAODEgamma/PhotonContainer.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODJet/JetContainer.h>
#include <xAODTrigger/JetRoIContainer.h>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTau/TauJetContainer.h>
#include <xAODMissingET/MissingETContainer.h>
#include <xAODTracking/TrackParticleContainer.h>

#include "xAODAnaHelpers/HelperClasses.h"
#include "xAODAnaHelpers/EventInfo.h"
#include "xAODAnaHelpers/MetContainer.h"
#include "xAODAnaHelpers/JetHelpTree.h"
#include "xAODAnaHelpers/ElectronContainer.h"
#include "xAODAnaHelpers/PhotonHelpTree.h"
#include "xAODAnaHelpers/FatJetHelpTree.h"
#include "xAODAnaHelpers/TruthHelpTree.h"
#include "xAODAnaHelpers/TrackContainer.h"
#include "xAODAnaHelpers/MuonContainer.h"
#include "xAODAnaHelpers/TauContainer.h"

#include "xAODAnaHelpers/TruthContainer.h"


#include <map>

// root includes
#include "TTree.h"
#include "TFile.h"

namespace TrigConf {
  class xAODConfigTool;
}

namespace Trig {
  class TrigDecisionTool;
}

typedef SG::AuxElement::Accessor< std::vector< float > > floatAccessor ;

class HelpTreeBase {

public:

  HelpTreeBase(xAOD::TEvent *event, TTree* tree, TFile* file, const float units = 1e3, bool debug = false, xAOD::TStore* store = nullptr );
  HelpTreeBase(TTree* tree, TFile* file, xAOD::TEvent *event = nullptr, xAOD::TStore* store = nullptr, const float units = 1e3, bool debug = false );
  virtual ~HelpTreeBase();

  void AddEvent       (const std::string& detailStr = "");
  void AddTrigger     (const std::string& detailStr = "");
  void AddJetTrigger  (const std::string& detailStr = "");
  void AddMuons       (const std::string& detailStr = "", const std::string& muonName = "muon");
  void AddElectrons   (const std::string& detailStr = "", const std::string& elecName = "el");
  void AddPhotons     (const std::string& detailStr = "", const std::string& photonName = "ph");
  void AddJets        (const std::string& detailStr = "", const std::string& jetName = "jet");
  void AddL1Jets      ();
  void AddTruth       (const std::string& detailStr = "", const std::string& truthName = "truth");
  void AddTrackParts  (const std::string& trackName,	 const std::string& detailStr = "");

  /**
   *  @brief  Declare a new collection of fatjets to be written to the output tree.
   *  @param  detailStr   A (space-separated) list of detail options. These keywords specify
   *                      exactly which information about each jet is written out. Current
   *                      influential options are: `kinematic` `substructure` `constituent`
   *                      `constituentAll`
   *  @param  fatjetName  The (prefix) name of the container. Default: `fatjet`.
   **/
  void AddFatJets     (const std::string& detailStr = "", const std::string& fatjetName = "fatjet");
  void AddTruthFatJets(const std::string& detailStr = "", const std::string& truthFatJetName = "truth_fatjet");

  void AddTaus        (const std::string& detailStr = "",  const std::string& tauName = "tau");
  void AddMET         (const std::string& detailStr = "");

  xAOD::TEvent* m_event;
  xAOD::TStore* m_store;

  // control which branches are filled
  HelperClasses::TriggerInfoSwitch*    m_trigInfoSwitch;
  HelperClasses::JetTriggerInfoSwitch* m_jetTrigInfoSwitch;
  HelperClasses::MuonInfoSwitch*       m_muonInfoSwitch;

  std::string                  m_triggerSelection;
  TrigConf::xAODConfigTool*    m_trigConfTool;
  Trig::TrigDecisionTool*      m_trigDecTool;

  void FillEvent( const xAOD::EventInfo* eventInfo, xAOD::TEvent* event = nullptr );

  void FillTrigger( const xAOD::EventInfo* eventInfo );
  void FillJetTrigger();

  void FillMuons( const xAOD::MuonContainer* muons, const xAOD::Vertex* primaryVertex, const std::string& muonName = "muon" );
  void FillMuon( const xAOD::Muon* muon, const xAOD::Vertex* primaryVertex, const std::string& muonName = "muon" );

  void FillElectrons( const xAOD::ElectronContainer* electrons, const xAOD::Vertex* primaryVertex, const std::string& elecName = "el" );
  void FillElectron ( const xAOD::Electron* elec, const xAOD::Vertex* primaryVertex, const std::string& elecName = "el" );

  void FillPhotons( const xAOD::PhotonContainer* photons, const std::string& photonName = "ph" );
  void FillPhoton ( const xAOD::Photon*          photon,  const std::string& photonName = "ph" );

  void FillJets( const xAOD::JetContainer* jets, int pvLocation = -1, const std::string& jetName = "jet" );
  void FillJet( const xAOD::Jet* jet_itr, const xAOD::Vertex* pv, int pvLocation, const std::string& jetName = "jet" );
  void FillL1Jets( const xAOD::JetRoIContainer* jets );

  void FillTruth( const xAOD::TruthParticleContainer* truth, const std::string& truthName = "truth" );
  void FillTruth( const xAOD::TruthParticle* truthPart, const std::string& truthName );

  void FillTracks( const std::string& trackName, const xAOD::TrackParticleContainer* tracks);
  void FillTrack( const xAOD::TrackParticle* trackPart, const std::string& trackName );

  /**
   *  @brief  Write a container of jets to the specified container name. The
   *          container name should be declared beforehand using `AddFatJets()`.
   *          This clears the current branch state for the collection so it only makes sense to
   *          call once per call to `Fill()`.
   *  @param  fatJets     A container of jets to be written out.
   *  @param  fatjetName  The name of the output collection to write to.
   */
  void FillFatJets( const xAOD::JetContainer* fatJets , const std::string& fatjetName = "fatjet");
  void FillFatJet ( const xAOD::Jet* fatjet_itr,        const std::string& fatjetName = "fatjet");

  void FillTruthFatJets( const xAOD::JetContainer* truthFatJets,     const std::string& truthFatJetName="truth_fatjet");
  void FillTruthFatJet ( const xAOD::Jet*          truth_fatjet_itr, const std::string& truthFatJetName="truth_fatjet");

  void FillTaus( const xAOD::TauJetContainer* taus, const std::string& tauName = "tau" );
  void FillTau ( const xAOD::TauJet* tau,           const std::string& tauName = "tau" );
  void FillMET( const xAOD::MissingETContainer* met );

  void Fill();
  void ClearEvent();
  void ClearTrigger();
  void ClearJetTrigger();
  void ClearMuons       (const std::string& muonName = "muon");
  void ClearElectrons   (const std::string& elecName = "el");
  void ClearPhotons     (const std::string& photonName = "ph");
  void ClearJets        (const std::string& jetName = "jet");
  void ClearL1Jets      ();
  void ClearTruth       (const std::string& truthName = "truth" );
  void ClearTracks	(const std::string& trackName);
  void ClearFatJets     (const std::string& fatjetName = "fatjet" );
  void ClearTruthFatJets(const std::string& truthFatJetName = "truth_fatjet" );
  void ClearTaus        (const std::string& tauName = "tau" );
  void ClearMET();

  bool writeTo( TFile *file );

  virtual void AddEventUser(const std::string& detailStr = "")      {
    if(m_debug) Info("AddEventUser","Empty function called from HelpTreeBase %s",detailStr.c_str());
    return;
  };

  virtual void AddTriggerUser(const std::string& detailStr = "")      {
    if(m_debug) Info("AddTriggerUser","Empty function called from HelpTreeBase %s",detailStr.c_str());
    return;
  };

  virtual void AddJetTriggerUser(const std::string& detailStr = "")      {
    if(m_debug) Info("AddJetTriggerUser","Empty function called from HelpTreeBase %s",detailStr.c_str());
    return;
  };

  virtual void AddMuonsUser(const std::string& detailStr = "")      {
    if(m_debug) Info("AddMuonsUser","Empty function called from HelpTreeBase %s",detailStr.c_str());
    return;
  };

  virtual void AddElectronsUser(const std::string& detailStr = "")  {
    if(m_debug) Info("AddElectronsUser","Empty function called from HelpTreeBase %s",detailStr.c_str());
    return;
  };

  virtual void AddTracksUser(const std::string& trackName, const std::string& detailStr = "")       {
    if(m_debug) Info("AddTracksUser","Empty function called from HelpTreeBase %s %s",trackName.c_str(), detailStr.c_str());
    return;
  };

  virtual void AddTausUser(const std::string& detailStr = "")       {
    if(m_debug) Info("AddTausUser","Empty function called from HelpTreeBase %s",detailStr.c_str());
    return;
  };

  virtual void AddMETUser(const std::string& detailStr = "")       {
    if(m_debug) Info("AddMETUser","Empty function called from HelpTreeBase %s",detailStr.c_str());
    return;
  };

  virtual void ClearEventUser       ()     { return; };
  virtual void ClearTriggerUser     ()   { return; };
  virtual void ClearMuonsUser       (const std::string& /*muonName = muon"*/)     { return; };
  virtual void ClearElectronsUser   (const std::string& /*elecName = "el"*/) { return; };
  virtual void ClearTracksUser       (const std::string& /*trackName*/)       { return; };
  virtual void ClearTausUser        (const std::string& /*tauName = "tau"*/) 	    { return; };
  virtual void ClearMETUser         ()       { return; };

  virtual void FillEventUser    ( const xAOD::EventInfo*  )        { return; };
  virtual void FillMuonsUser    ( const xAOD::Muon*,     const std::string& /*muonName = "muon"*/  )             { return; };
  virtual void FillElectronsUser( const xAOD::Electron*, const std::string& /*elecName = "el"*/ )     { return; };
  virtual void FillTracksUser   ( const std::string& /*trackName*/, const xAOD::TrackParticle*  )               { return; };
  virtual void FillTausUser( const xAOD::TauJet*,           const std::string& /*tauName = "tau"*/  )            { return; };
  virtual void FillMETUser( const xAOD::MissingETContainer*  ) { return; };
  virtual void FillTriggerUser( const xAOD::EventInfo*  )      { return; };
  virtual void FillJetTriggerUser()                            { return; };

 protected:

  template<typename T, typename U, typename V>
    void safeFill(const V* xAODObj, SG::AuxElement::ConstAccessor<T>& accessor, std::vector<U>& destination, U defaultValue, int m_units = 1);

  template<typename T, typename U, typename V>
    void safeVecFill(const V* xAODObj, SG::AuxElement::ConstAccessor<std::vector<T> >& accessor, std::vector<std::vector<U> >& destination, int m_units = 1);

  template<typename T>
    void setBranch(const std::string& prefix, const std::string& varName, std::vector<T>* localVectorPtr);

protected:

  TTree* m_tree;

  int m_units; //For MeV to GeV conversion in output

  bool m_debug;
  bool m_isMC;

  // event
  xAH::EventInfo*      m_eventInfo;

  // trigger
  int m_passL1;
  int m_passHLT;
  unsigned int m_masterKey;
  unsigned int m_L1PSKey;
  unsigned int m_HLTPSKey;
  std::vector<std::string>  m_elTrigForMatching;   /* each event can have a list of electron trigger chains to which each electron could be matched.
  						   / This list is created when configuring ElectronSelector.cxx, where the electron trigger matching is actually performed
						   */

  // jet trigger
  std::vector<std::string> m_passTriggers;
  std::vector<float> m_triggerPrescales;
  std::vector<std::string>  m_isPassBitsNames;
  std::vector<unsigned int> m_isPassBits;

  //
  //  Jets
  //
  std::map<std::string, xAH::JetHelpTree*> m_jets;

  //
  // L1 Jets
  //
  int m_nL1Jet;
  std::vector<float> m_l1Jet_et8x8;
  std::vector<float> m_l1Jet_eta;
  std::vector<float> m_l1Jet_phi;

  //
  // Truth
  //
  std::map<std::string, xAH::TruthHelpTree*> m_truth;

  //
  // Tracks
  //
  std::map<std::string, xAH::TrackContainer*> m_tracks;

  //
  // fat jets
  //
  std::map<std::string, xAH::FatJetHelpTree*> m_fatjets;

  //
  // truth fat jets
  //
  std::map<std::string, xAH::FatJetHelpTree*> m_truth_fatjets;

  //
  // muons
  //
  std::map<std::string, xAH::MuonContainer*> m_muons;
  std::map<std::string, std::vector<std::string> > m_RecoEff_SF_sysNames;
  std::map<std::string, std::vector<std::string> > m_IsoEff_SF_sysNames;
  std::map<std::string, std::vector<std::string> > m_TrigEff_SF_sysNames;
  std::vector<std::string>  m_TTVAEff_SF_sysNames;
  
  //
  // electrons
  //
  std::map<std::string, xAH::ElectronContainer*> m_elecs;
  
  //
  // photons
  //
  std::map<std::string, xAH::PhotonHelpTree*> m_photons;

  //
  // taus
  //
  std::map<std::string, xAH::TauContainer*> m_taus;

  // met
  xAH::MetContainer*      m_met;

};


template<typename T, typename U, typename V>
void HelpTreeBase::safeFill(const V* xAODObj, SG::AuxElement::ConstAccessor<T>& accessor, std::vector<U>& destination, U defaultValue, int m_units){
  if ( accessor.isAvailable( *xAODObj ) ) {
    destination.push_back( accessor( *xAODObj ) / m_units );
  } else {
    destination.push_back( defaultValue );
  }
  return;
}

template<typename T, typename U, typename V>
void HelpTreeBase::safeVecFill(const V* xAODObj, SG::AuxElement::ConstAccessor<std::vector<T> >& accessor, std::vector<std::vector<U> >& destination, int m_units){
  destination.push_back( std::vector<U>() );

  if ( accessor.isAvailable( *xAODObj ) ) {
    for(U itemInVec : accessor(*xAODObj))        destination.back().push_back(itemInVec / m_units);
  }
  return;
}


template<typename T>
void HelpTreeBase::setBranch(const std::string& prefix, const std::string& varName, std::vector<T>* localVectorPtr){
  m_tree->Branch((prefix+"_"+varName).c_str(),        localVectorPtr);
}


#endif

