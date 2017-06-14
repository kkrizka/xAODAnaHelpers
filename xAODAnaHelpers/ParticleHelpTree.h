#ifndef xAODAnaHelpers_ParticleHelpTree_H
#define xAODAnaHelpers_ParticleHelpTree_H

#include <TTree.h>
#include <TClonesArray.h>

#include <vector>
#include <string>

#include <xAODAnaHelpers/HelperClasses.h>
#include <xAODAnaHelpers/HelperFunctions.h>

#include <xAODAnaHelpers/Particle.h>
#include <xAODBase/IParticle.h>

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
		     bool useMass=false, 
		     const std::string& suffix="")
      : m_name(name),
	m_infoSwitch(detailStr), 
	m_mc(mc),
	m_debug(false), 
	m_units(units), 
	m_useMass(useMass), 
	m_suffix(suffix)
    {
      m_particles=new TClonesArray(className.c_str());
    }

    virtual ~ParticleHelpTree()
    { }

    virtual void createBranches(TTree *tree)
    {
      tree->Branch(branchName().c_str(),&m_particles);

      setBranchStatus(tree, "*",    false);
      if(m_infoSwitch.m_kinematic)
	setBranchStatus(tree, "p4", true);
    }

    virtual void clear()
    {
      m_particles->Clear();
    }

    virtual void fillParticle(const xAOD::IParticle* particle, xAH::Particle* myparticle)
    {
      if(m_infoSwitch.m_kinematic)
	{
	  if(m_useMass)
	    myparticle->p4.SetPtEtaPhiE(particle->pt() / m_units,
					particle->eta(),
					particle->phi(),
					particle->e() / m_units);
	  else
	    myparticle->p4.SetPtEtaPhiE(particle->pt() / m_units,
					particle->eta(),
					particle->phi(),
					particle->m() / m_units);
	}
    }

  protected:
    std::string branchName()
    {
      std::string name = m_name;
      if (! m_suffix.empty()) { name += "_" + m_suffix; }
      return name;
    }

    void setBranchStatus(TTree *tree, const std::string& varName, bool status)
    {
      tree->SetBranchStatus((branchName()+"."+varName).c_str(), status); 
    }

    template<typename T, typename U, typename V> void safeFill(const V* xAODObj, SG::AuxElement::ConstAccessor<T>& accessor, std::vector<U>* destination, U defaultValue, int units = 1){
      if ( accessor.isAvailable( *xAODObj ) ) {
	destination->push_back( accessor( *xAODObj ) / units );
      } else {
	destination->push_back( defaultValue );
      }
    }

    template<typename T, typename U, typename V> void safeVecFill(const V* xAODObj, SG::AuxElement::ConstAccessor<std::vector<T> >& accessor, std::vector<std::vector<U> >* destination, int units = 1){
      destination->push_back( std::vector<U>() );

      if ( accessor.isAvailable( *xAODObj ) ) {
	for(U itemInVec : accessor(*xAODObj))        destination->back().push_back(itemInVec / units);
      } 
      return;
    }

  public:
    std::string m_name;
    T_INFOSWITCH m_infoSwitch;
    bool  m_mc;
    bool  m_debug;
    float m_units;

  private:
    bool          m_useMass;
    std::string   m_suffix;
    TClonesArray *m_particles;
  };

}//xAH
#endif // xAODAnaHelpers_ParticleHelpTree_H
