// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "PiToPSubstitution.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PiToPSubstitution
//
// 2013-11-18 : Robert Zillmer (LHCb)
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( PiToPSubstitution )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PiToPSubstitution::PiToPSubstitution( const std::string& type,
                                      const std::string& name,
                                      const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IParticleManipulator>(this);

}
//=============================================================================
// Destructor
//=============================================================================
PiToPSubstitution::~PiToPSubstitution() {}
//=============================================================================

StatusCode PiToPSubstitution::initialize()
{
  _reFit = tool<IParticleReFitter>( "OfflineVertexFitter", this );
  _ppSvc = svc<LHCb::IParticlePropertySvc>( "LHCb::ParticlePropertySvc", true );
  
  const LHCb::ParticleProperty* piInfo = _ppSvc->find( "pi+" );
  const LHCb::ParticleProperty* pInfo  = _ppSvc->find( "p+" );
  
  _piMass = piInfo->mass();
  _pMass  = pInfo->mass();
  
  _pipID  = piInfo->particleID();
  _pimID  = piInfo->antiParticle()->particleID();
  _pID    = pInfo->particleID();
  _pbarID = pInfo->antiParticle()->particleID();
  
  return StatusCode::SUCCESS;
}

StatusCode  PiToPSubstitution::doCorrection( LHCb::Particle* particle )
{
   // daughter particles stored in [pip, pim, L] order
  LHCb::Particle::ConstVector daughter = particle->daughtersVector();
 
  debug() << "Taking a " << (daughter.at(0))->particleID().pid() << " at first" << endmsg;  

  // pip -> p
  const Gaudi::LorentzVector& oldMomentum1 = (daughter.at(0))->momentum();
  const double newEnergy1 = std::sqrt( oldMomentum1.P2() + _pMass*_pMass );
  const_cast<LHCb::Particle*> (daughter.at(0))->setParticleID( _pID );
  const_cast<LHCb::Particle*> (daughter.at(0))->setMomentum( Gaudi::LorentzVector( oldMomentum1.Px(),
                                                                                   oldMomentum1.Py(), 
                                                                                   oldMomentum1.Pz(), 
                                                                                   newEnergy1 ) );
  const_cast<LHCb::Particle*> (daughter.at(0))->setMeasuredMass( _pMass );


  debug() << "Taking a " << (daughter.at(1))->particleID().pid() << endmsg;

  // pim -> pbar
  const Gaudi::LorentzVector& oldMomentum2 = (daughter.at(1))->momentum();
  const double newEnergy2 = std::sqrt( oldMomentum2.P2() + _pMass*_pMass );
  const_cast<LHCb::Particle*> (daughter.at(1))->setParticleID( _pbarID );
  const_cast<LHCb::Particle*> (daughter.at(1))->setMomentum( Gaudi::LorentzVector( oldMomentum2.Px(), 
                                                                       oldMomentum2.Py(), 
                                                                       oldMomentum2.Pz(), 
                                                                       newEnergy2 ) );
  const_cast<LHCb::Particle*> (daughter.at(1))->setMeasuredMass( _pMass );

  debug() << "Refit of " << particle->particleID().pid() << endmsg;
  StatusCode sc = _reFit->reFit( *particle );
  return sc;
}

StatusCode  PiToPSubstitution::undoCorrection( LHCb::Particle* particle )
{
  // daughter particles stored in [p, pbar, L] order
  LHCb::Particle::ConstVector daughter = particle->daughtersVector();
  LHCb::Particle::ConstVector::const_iterator i = daughter.begin();

  // p -> pip
  const Gaudi::LorentzVector& oldMomentum1 = (*i)->momentum();
  const double newEnergy1 = std::sqrt( oldMomentum1.P2() + _piMass*_piMass );
  const_cast<LHCb::Particle*> (*i)->setParticleID( _pipID );
  const_cast<LHCb::Particle*> (*i)->setMomentum( Gaudi::LorentzVector( oldMomentum1.Px(),
                                                                       oldMomentum1.Py(), 
                                                                       oldMomentum1.Pz(), 
                                                                       newEnergy1 ) );
  const_cast<LHCb::Particle*> (*i)->setMeasuredMass( _piMass );
  
  i++;
  
  // pbar -> pim
  const Gaudi::LorentzVector& oldMomentum2 = (*i)->momentum();
  const double newEnergy2 = std::sqrt( oldMomentum2.P2() + _piMass*_piMass );
  const_cast<LHCb::Particle*> (*i)->setParticleID( _pimID );
  const_cast<LHCb::Particle*> (*i)->setMomentum( Gaudi::LorentzVector( oldMomentum2.Px(),
                                           oldMomentum2.Py(), 
                                           oldMomentum2.Pz(), 
                                           newEnergy2 ) );
  const_cast<LHCb::Particle*> (*i)->setMeasuredMass( _piMass );

  debug() << "Refit of " << particle->particleID().pid() << endmsg;
  StatusCode sc = _reFit->reFit( *particle );  
  return sc;
}
