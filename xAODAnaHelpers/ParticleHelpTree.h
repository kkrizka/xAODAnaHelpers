#ifndef xAODAnaHelpers_ParticleHelpTree_H
#define xAODAnaHelpers_ParticleHelpTree_H

#include <TTree.h>
#include <TClonesArray.h>

#include <vector>
#include <string>

#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODAnaHelpers/Particle.h>

#define SAFE_SET(obj, variable, xaodobj) if(variable.isAvailable( *xaodobj )) obj->variable=variable( *xaodobj );

namespace xAH {

  template <class T_PARTICLE, class T_INFOSWITCH>
  class ParticleHelpTree
  {
  public:
    ParticleHelpTree(const std::string& className,
		     const std::string& name, 
		     const std::string& detailStr="", 
		     float units = 1e3, 
		     bool mc = false, 
		     bool useMass=false)
      : m_name(name),
	m_infoSwitch(detailStr), 
	m_mc(mc),
	m_debug(false), 
	m_units(units), 
	m_useMass(useMass)
    {
      m_particles=new TClonesArray(className.c_str());
    }

    virtual ~ParticleHelpTree()
    { }

    int size()
    { return m_particles->GetEntries(); }

    TObject* particle(uint idx) const
    { return m_particles->At(idx); }

    virtual void createBranches(TTree *tree)
    {
      tree->Branch(branchName().c_str(),&m_particles);

      setBranchStatus(tree, "*",    false);
      if(m_infoSwitch.m_kinematic)
	setBranchStatus(tree, "p4", true);
    }

    virtual void clear()
    {
      m_particles->Clear("C");
    }

    virtual void fillParticle(const xAOD::IParticle* particle)
    {
      new((*m_particles)[m_particles->GetEntries()]) T_PARTICLE();
      xAH::Particle* myparticle=static_cast<xAH::Particle*>(m_particles->Last());

      if(m_infoSwitch.m_kinematic)
	{
	  if(m_useMass)
	    myparticle->p4.SetPtEtaPhiE(particle->pt() / m_units,
					particle->eta(),
					particle->phi(),
					particle->e() / m_units);
	  else
	    myparticle->p4.SetPtEtaPhiM(particle->pt() / m_units,
					particle->eta(),
					particle->phi(),
					particle->m() / m_units);
	}
    }

  protected:
    std::string branchName()
    {
      return m_name;
    }

    void setBranchStatus(TTree *tree, const std::string& varName, bool status)
    {
      tree->SetBranchStatus((branchName()+"."+varName).c_str(), status); 
    }

  public:
    std::string m_name;
    T_INFOSWITCH m_infoSwitch;
    bool  m_mc;
    bool  m_debug;
    float m_units;

  protected:
    TClonesArray *m_particles;

  private:
    bool          m_useMass;
    std::string   m_suffix;
  };

}//xAH
#endif // xAODAnaHelpers_ParticleHelpTree_H
