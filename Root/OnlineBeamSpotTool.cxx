#include <xAODAnaHelpers/OnlineBeamSpotTool.h>
#include "PathResolver/PathResolver.h"
#include <iostream>

// ROOT include(s):
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

using namespace xAH;

OnlineBeamSpotTool::OnlineBeamSpotTool() :
  m_cachedRunNum(-1),
  m_cachedLB(-1),
  m_cachedRunInfo(nullptr),
  m_cachedLBData(nullptr),
  m_mcLBData(nullptr)
{
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.A.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.B.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.C.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.D.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.E.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.F.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.G.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.H.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.I.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.K.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2016.L.root");
  readFile("xAODAnaHelpers/OnlineBSInfo/OnlineBSInfo.2017.root");

  m_mcLBData = new LBData(0,999999,0,0,0);
}

OnlineBeamSpotTool::~OnlineBeamSpotTool()
{
  //std::cout << "In ~OnlineBeamSpotTool" << std::endl;
}


void OnlineBeamSpotTool::setRunInfo(int runNumber){
  RunToLBDataMapItr it = m_runList.find(runNumber);

  if(it != m_runList.end()){
    m_cachedRunInfo = &(it->second);
  } else {
    m_cachedRunInfo = nullptr;
  }
  m_cachedRunNum = runNumber;
  return;
}

const OnlineBeamSpotTool::LBData* OnlineBeamSpotTool::getLBData(int lumiBlock){
  if(!m_cachedRunInfo) return nullptr;

  for(const LBData& thisLBData : *m_cachedRunInfo){
    if((lumiBlock >= thisLBData.m_LBStart) && (lumiBlock <= thisLBData.m_LBEnd)){
      return &thisLBData;
    }
  }

  return nullptr;
}

const OnlineBeamSpotTool::LBData* OnlineBeamSpotTool::getLBData(int runNumber, int lumiBlock, bool isMC){

  //
  // Check MC
  //
  if(isMC)
    return m_mcLBData;

  //
  // Check cached data
  //
  if((runNumber == m_cachedRunNum) && (lumiBlock == m_cachedLB))
    return m_cachedLBData;

  //
  // GetRun info
  //
  if(runNumber != m_cachedRunNum)
    setRunInfo(runNumber);

  return getLBData(lumiBlock);
}

float OnlineBeamSpotTool::getOnlineBSInfo(const xAOD::EventInfo* eventInfo, OnlineBeamSpotTool::BSData datakey){
  return getOnlineBSInfo(eventInfo->runNumber(), eventInfo->lumiBlock(), eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ), datakey);
}

float OnlineBeamSpotTool::getOnlineBSInfo(const xAH::EventInfo* eventInfo, OnlineBeamSpotTool::BSData datakey){
  return getOnlineBSInfo(eventInfo->m_runNumber, eventInfo->m_lumiBlock, eventInfo->m_mc, datakey);
}

float OnlineBeamSpotTool::getOnlineBSInfo(int runNumber, int lumiBlock, bool isMC, OnlineBeamSpotTool::BSData datakey){
  //std::cout << "In OnlineBeamSpotTool (" << runNumber << " , " << lumiBlock<< ")" << std::endl;

  const LBData* thisLBInfo = getLBData(runNumber, lumiBlock, isMC);

  if(!thisLBInfo){
    std::cout << "OnlineBeamSpotTool::ERROR no LB data for run: " << runNumber << " LB: " << lumiBlock << std::endl;
    return -999;
  }

  if(datakey == OnlineBeamSpotTool::BSx)
    return thisLBInfo->m_BSx;

  if(datakey == OnlineBeamSpotTool::BSy)
    return thisLBInfo->m_BSy;

  return thisLBInfo->m_BSz;
}

void OnlineBeamSpotTool::readFile(std::string rootFileName){

  std::string fullRootFileName = PathResolverFindCalibFile( rootFileName );

  TFile* thisFile = new TFile(fullRootFileName.c_str(),"READ");
  TTree* tree = (TTree*)thisFile->Get("LBInfo");

  int RunNumber;
  std::vector<int>*   LBStart  = new std::vector<int>();
  std::vector<int>*   LBEnd    = new std::vector<int>();

  std::vector<float>* BSx      = new std::vector<float>();
  std::vector<float>* BSy      = new std::vector<float>();
  std::vector<float>* BSz      = new std::vector<float>();

  tree->SetBranchAddress("RunNumber",&RunNumber);
  tree->SetBranchAddress("LBStart",  &LBStart);
  tree->SetBranchAddress("LBEnd",    &LBEnd);
  tree->SetBranchAddress("BSx",      &BSx);
  tree->SetBranchAddress("BSy",      &BSy);
  tree->SetBranchAddress("BSz",      &BSz);

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i=0;i<nentries;i++) {
    tree->GetEntry(i);
    RunInfo thisRunInfo;

    for(unsigned int LBIt = 0; LBIt < LBStart->size(); ++LBIt){
      thisRunInfo.push_back(LBData(LBStart ->at(LBIt),
				   LBEnd   ->at(LBIt),
				   BSx     ->at(LBIt),
				   BSy     ->at(LBIt),
				   BSz     ->at(LBIt)
				   ));
    }

    m_runList.insert( std::make_pair(RunNumber, thisRunInfo) );
  }

  thisFile->Close();
}
