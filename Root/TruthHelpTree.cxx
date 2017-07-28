#include "xAODAnaHelpers/TruthHelpTree.h"

#include <xAODAnaHelpers/HelperFunctions.h>
#include <iostream>
#include <xAODTruth/TruthVertex.h>
#include <xAODTruth/TruthParticle.h>
//#include "xAODTruth/TruthEventHelpTree.h"

using namespace xAH;

TruthHelpTree::TruthHelpTree(const std::string& name, const std::string& detailStr, float units)
  : ParticleHelpTree("xAH::TruthParticle", name,detailStr,units,true)

{ }

TruthHelpTree::~TruthHelpTree()
{ }

void TruthHelpTree::createBranches(TTree *tree)
{
  ParticleHelpTree::createBranches(tree);

  setBranchStatus(tree,"pdgId",   1);
  setBranchStatus(tree,"status",  1);
  setBranchStatus(tree,"barcode", 1);

  if(m_infoSwitch.m_type)
    {
      setBranchStatus(tree,"is_higgs", 1);
      setBranchStatus(tree,"is_bhad",  1);
    }

  if(m_infoSwitch.m_bVtx)
    {
      setBranchStatus(tree,"Bdecay_x", 1);
      setBranchStatus(tree,"Bdecay_y", 1);
      setBranchStatus(tree,"Bdecay_z", 1);
    }

  if(m_infoSwitch.m_parents)
    {
      setBranchStatus(tree,"nParents",       1);
      setBranchStatus(tree,"parent_pdgId",   1);
      setBranchStatus(tree,"parent_barcode", 1);
      setBranchStatus(tree,"parent_status",  1);
    }

  if(m_infoSwitch.m_children)
    {
      setBranchStatus(tree,"nChildren",     1);
      setBranchStatus(tree,"child_pdgId",   1);
      setBranchStatus(tree,"child_barcode", 1);
      setBranchStatus(tree,"child_status",  1);
  }
}

void TruthHelpTree::fillTruth( const xAOD::TruthParticle* truth )
{
  ParticleHelpTree::fillParticle(truth);
  xAH::TruthParticle* mytruth=static_cast<xAH::TruthParticle*>(m_particles->Last());

  mytruth->pdgId  = truth->pdgId  ();
  mytruth->status = truth->status ();
  mytruth->barcode= truth->barcode();

  if(m_infoSwitch.m_type)
    {
      mytruth->is_higgs= truth->isHiggs       ();
      mytruth->is_bhad = truth->isBottomHadron();
    }



  if(m_infoSwitch.m_bVtx)
    {
      if(truth->isBottomHadron() && truth->hasDecayVtx())
	{
	  const xAOD::TruthVertex* vtx = truth->decayVtx();
	  mytruth->Bdecay_x=vtx->x();
	  mytruth->Bdecay_y=vtx->y();
	  mytruth->Bdecay_z=vtx->z();
	}
    }

  if(m_infoSwitch.m_parents)
    {
      int nParents = truth->nParents();
      mytruth->nParents=nParents;

      for(int iparent = 0; iparent < nParents; ++iparent)
	{
	  const xAOD::TruthParticle* parent = truth->parent(iparent);
	  if(parent)
	    {
	      mytruth->parent_pdgId  .push_back(parent->pdgId  ());
	      mytruth->parent_barcode.push_back(parent->barcode());
	      mytruth->parent_status .push_back(parent->status ());
	    }
	}
    }

  if(m_infoSwitch.m_children)
    {
      int nChildren = truth->nChildren();
      mytruth->nChildren=nChildren;

      for(int ichild = 0; ichild < nChildren; ++ichild)
	{
	  const xAOD::TruthParticle* child = truth->child(ichild);
	  if(child)
	    {
	      mytruth->child_pdgId  .push_back(child->pdgId  ());
	      mytruth->child_barcode.push_back(child->barcode());
	      mytruth->child_status .push_back(child->status ());
	    }
	}
    }

  return;
}

