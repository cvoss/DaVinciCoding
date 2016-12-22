// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "LambdaToKsSubstitution.h"

//-----------------------------------------------------------------------------
// Implementation file for class : LambdaToKsSubstitution
//
// 2014-06-23 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( LambdaToKsSubstitution )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
LambdaToKsSubstitution::LambdaToKsSubstitution( const std::string& type,
                                                const std::string& name,
                                                const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IParticleManipulator>(this);

}
//=============================================================================
// Destructor
//=============================================================================
LambdaToKsSubstitution::~LambdaToKsSubstitution() {} 

//=============================================================================

StatusCode LambdaToKsSubstitution::initialize()
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

StatusCode  LambdaToKsSubstitution::doCorrection( LHCb::Particle* particle )
{
   // daughter particles stored in [pip, pim, L] order
  LHCb::Particle::ConstVector daughters      = particle->daughtersVector();
  LHCb::Particle::ConstVector grandDaughters = daughters.at(2)->daughtersVector();

  // pip -> p
  const Gaudi::LorentzVector& pMom = (grandDaughters.at(0))->momentum();
  const double newEnergy = std::sqrt( pMom.P2() + _piMass*_piMass );
  const_cast<LHCb::Particle*> (grandDaughters.at(0))->setParticleID( grandDaughters.at(0)->particleID() == _pID ? _pipID : _pimID );
  const_cast<LHCb::Particle*> (grandDaughters.at(0))->setMomentum( Gaudi::LorentzVector( pMom.Px(),
                                                                                   pMom.Py(), 
                                                                                   pMom.Pz(), 
                                                                                   newEnergy ) );
  const_cast<LHCb::Particle*> (grandDaughters.at(0))->setMeasuredMass( _piMass );
  info() << " pi Mass " << daughters.at(2)->daughtersVector().at(0)->momentum().M()  << endmsg;
  LHCb::Particle *L = const_cast<LHCb::Particle*>(daughters.at(2));
  
  StatusCode sc = _reFit->reFit( *L );
  info() << " Ks Mass " << daughters.at(2)->momentum().M()  << endmsg;
  return sc;
}

StatusCode  LambdaToKsSubstitution::undoCorrection( LHCb::Particle* particle )
{ 
  LHCb::Particle::ConstVector daughters      = particle->daughtersVector();
  LHCb::Particle::ConstVector grandDaughters = daughters.at(2)->daughtersVector();

  const Gaudi::LorentzVector& pMom = (grandDaughters.at(0))->momentum();
  const double newEnergy = std::sqrt( pMom.P2() + _pMass*_pMass );
  
  const_cast<LHCb::Particle*> (grandDaughters.at(0))->setParticleID( grandDaughters.at(0)->particleID() == _pipID ? _pID : _pbarID );
  const_cast<LHCb::Particle*> (grandDaughters.at(0))->setMomentum( Gaudi::LorentzVector( pMom.Px(),
                                                                                         pMom.Py(), 
                                                                                         pMom.Pz(), 
                                                                                         newEnergy ) );
  const_cast<LHCb::Particle*> (grandDaughters.at(0))->setMeasuredMass( _pMass );
  info() << " p Mass " << daughters.at(2)->daughtersVector().at(0)->momentum().M()  << endmsg;
  LHCb::Particle *L = const_cast<LHCb::Particle*>(daughters.at(2));
    
  StatusCode sc = _reFit->reFit( *L );  
  info() << " Lambda Mass " << daughters.at(2)->momentum().M()  << endmsg;
  return sc;
}
