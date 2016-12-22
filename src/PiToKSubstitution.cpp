// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "PiToKSubstitution.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PiToKSubstitution
//
// 2013-10-29 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( PiToKSubstitution )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PiToKSubstitution::PiToKSubstitution( const std::string& type,
                                      const std::string& name,
                                      const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IParticleManipulator>(this);

}
//=============================================================================
// Destructor
//=============================================================================
PiToKSubstitution::~PiToKSubstitution() {}
//=============================================================================
StatusCode PiToKSubstitution::initialize()
{
  m_reFit = tool<IParticleReFitter>("OfflineVertexFitter",this);
  m_ppSvc = svc<LHCb::IParticlePropertySvc>("LHCb::ParticlePropertySvc",true);

  const LHCb::ParticleProperty* piInfo = m_ppSvc->find( "pi+" );
  const LHCb::ParticleProperty* KInfo  = m_ppSvc->find( "K+" );

  m_piMass = piInfo->mass();
  m_KMass = KInfo->mass();

  m_pimID = piInfo->antiParticle()->particleID();
  m_pipID = piInfo->particleID();
  m_KmID =  KInfo->antiParticle()->particleID();
  m_KpID = KInfo->particleID();

  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode PiToKSubstitution::doCorrection( LHCb::Particle* Part) 
{
  LHCb::Particle::ConstVector PartDau = Part->daughtersVector();
  LHCb::Particle::ConstVector::const_iterator i;
  for( i = PartDau.begin(); i != PartDau.end(); ++i){
    if ((*i)->particleID().abspid() == m_pipID.abspid() ){
      const Gaudi::LorentzVector& oldMom1 = (*i)->momentum() ;
      const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_KMass*m_KMass ) ;
      const_cast<LHCb::Particle*>((*i)) -> setParticleID( (*i)->charge() > 0 ? m_KpID : m_KmID); 
      const_cast<LHCb::Particle*>((*i)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
									      oldMom1.Py() ,
									      oldMom1.Pz() , 
									      newEnergy1  ) );
      const_cast<LHCb::Particle*>((*i)) -> setMeasuredMass( m_KMass );
    }
   if ((*i)->particleID().abspid() == m_KpID.abspid() ){
      const Gaudi::LorentzVector& oldMom1 = (*i)->momentum() ;
      const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_piMass*m_piMass ) ;
      const_cast<LHCb::Particle*>((*i)) -> setParticleID( (*i)->charge() > 0 ? m_pipID : m_pimID ); 
      const_cast<LHCb::Particle*>((*i)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
									      oldMom1.Py() ,
									      oldMom1.Pz() , 
									      newEnergy1  ) );
      const_cast<LHCb::Particle*>((*i)) -> setMeasuredMass( m_piMass ); 
   }
  }
  
  debug() << "Refit of " << Part->particleID().pid() << endmsg;
  StatusCode sc = m_reFit -> reFit(*Part);
  return sc;
}
//=============================================================================
StatusCode PiToKSubstitution::undoCorrection( LHCb::Particle* Part) 
{
  LHCb::Particle::ConstVector PartDau = Part->daughtersVector();
  LHCb::Particle::ConstVector::const_iterator i;
  for( i = PartDau.begin(); i != PartDau.end(); ++i){
    if ((*i)->particleID().abspid() == m_pipID.abspid() ){
      const Gaudi::LorentzVector& oldMom1 = (*i)->momentum() ;
      const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_KMass*m_KMass ) ;
      const_cast<LHCb::Particle*>((*i)) -> setParticleID( (*i)->charge() > 0 ? m_KpID : m_KmID); 
      const_cast<LHCb::Particle*>((*i)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
									      oldMom1.Py() ,
									      oldMom1.Pz() , 
									      newEnergy1  ) );
      const_cast<LHCb::Particle*>((*i)) -> setMeasuredMass( m_KMass );
    }
   if ((*i)->particleID().abspid() == m_KpID.abspid() ){
      const Gaudi::LorentzVector& oldMom1 = (*i)->momentum() ;
      const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_piMass*m_piMass ) ;
      const_cast<LHCb::Particle*>((*i)) -> setParticleID( (*i)->charge() > 0 ? m_pipID : m_pimID ); 
      const_cast<LHCb::Particle*>((*i)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
									      oldMom1.Py() ,
									      oldMom1.Pz() , 
									      newEnergy1  ) );
      const_cast<LHCb::Particle*>((*i)) -> setMeasuredMass( m_piMass ); 
   }
  }
  
  debug() << "Refit of " << Part->particleID().pid() << endmsg;
  StatusCode sc = m_reFit -> reFit(*Part);
  return sc;
}
