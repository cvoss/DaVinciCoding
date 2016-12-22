// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "BLppiToolNoCorrect.h"

//-----------------------------------------------------------------------------
// Implementation file for class : BLppiToolNoCorrect
//
// 2014-07-09 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( BLppiToolNoCorrect )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
BLppiToolNoCorrect::BLppiToolNoCorrect( const std::string& type,
                                        const std::string& name,
                                        const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IParticleManipulator>(this);

}
//=============================================================================
// Destructor
//=============================================================================
BLppiToolNoCorrect::~BLppiToolNoCorrect() {} 

//=============================================================================
StatusCode BLppiToolNoCorrect::initialize()
{
  m_reFit = tool<IParticleReFitter>("OfflineVertexFitter",this);
  m_ppSvc = svc<LHCb::IParticlePropertySvc>("LHCb::ParticlePropertySvc",true);

  const LHCb::ParticleProperty* BInfo  = m_ppSvc->find( "B0" ) ;
  const LHCb::ParticleProperty* LbInfo = m_ppSvc->find( "Lambda_b0" );
  const LHCb::ParticleProperty* pInfo  = m_ppSvc->find( "p+" );
  const LHCb::ParticleProperty* piInfo = m_ppSvc->find( "pi+" );

  m_pMass = pInfo->mass();
  m_piMass = piInfo->mass();

  m_BID = BInfo->particleID();
  m_BbarID = BInfo->antiParticle()->particleID();
  m_LbID = LbInfo->particleID();
  m_LbbarID = LbInfo->antiParticle()->particleID();
  m_pID = pInfo->particleID();
  m_pbarID = pInfo->antiParticle()->particleID();
  m_pimID = piInfo->antiParticle()->particleID();
  m_pipID = piInfo->particleID();

  return StatusCode::SUCCESS;
}

//=============================================================================
StatusCode BLppiToolNoCorrect::doCorrection( LHCb::Particle* Part) 
{
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode BLppiToolNoCorrect::undoCorrection( LHCb::Particle* Part) 
{
  debug() << "Start Re-Substitution" << endmsg;
  LHCb::Particle::ConstVector PartDau = Part->daughtersVector();
  if(msgLevel(MSG::DEBUG)){
    LHCb::Particle::ConstVector::const_iterator i;
    for( i = PartDau.begin(); i != PartDau.end(); ++i)
      debug() << "PID of daughter " << (*i)->particleID().pid() << endmsg;
  }
  if(Part->particleID() == m_BbarID){
    const Gaudi::LorentzVector& oldMom = PartDau.at(1)->momentum() ;
    const double newEnergy = std::sqrt ( oldMom.P2() + m_piMass*m_piMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setParticleID(m_pimID); 
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMomentum( Gaudi::LorentzVector( oldMom.Px() ,
										     oldMom.Py() ,
										     oldMom.Pz() , 
										     newEnergy  ) );
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMeasuredMass( m_piMass );
  }
  if(Part->particleID() == m_BID){
    const Gaudi::LorentzVector& oldMom = PartDau.at(0)->momentum() ;
    const double newEnergy = std::sqrt ( oldMom.P2() + m_piMass*m_piMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setParticleID(m_pipID); 
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMomentum( Gaudi::LorentzVector( oldMom.Px() ,
										     oldMom.Py() ,
										     oldMom.Pz() , 
										     newEnergy  ) );
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMeasuredMass( m_piMass );
  }
  debug() << "Substitution of " << Part->particleID().pid() << endmsg;
  Part->setParticleID( ( Part->particleID().pid() < 0 ? m_LbID : m_LbbarID ) );
  debug() << "Refit of " << Part->particleID().pid() << endmsg;
  StatusCode sc = m_reFit -> reFit(*Part);
  return sc;
}
