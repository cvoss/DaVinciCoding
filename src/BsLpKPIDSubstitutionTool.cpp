// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "BsLpKPIDSubstitutionTool.h"

//-----------------------------------------------------------------------------
// Implementation file for class : BsLpKPIDSubstitutionTool
//
// 2013-07-29 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( BsLpKPIDSubstitutionTool )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
BsLpKPIDSubstitutionTool::BsLpKPIDSubstitutionTool( const std::string& type,
                                                    const std::string& name,
                                                    const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IParticleManipulator>(this);

}
//=============================================================================
// Destructor
//=============================================================================
BsLpKPIDSubstitutionTool::~BsLpKPIDSubstitutionTool() {} 

//=============================================================================
StatusCode BsLpKPIDSubstitutionTool::initialize()
{
  m_reFit = tool<IParticleReFitter>("OfflineVertexFitter",this);
  m_ppSvc = svc<LHCb::IParticlePropertySvc>("LHCb::ParticlePropertySvc",true);

  const LHCb::ParticleProperty* BInfo  = m_ppSvc->find( "B_s0" ) ;
  const LHCb::ParticleProperty* LbInfo = m_ppSvc->find( "Lambda_b0" );
  const LHCb::ParticleProperty* pInfo  = m_ppSvc->find( "p+" );
  const LHCb::ParticleProperty* piInfo = m_ppSvc->find( "pi+" );
  const LHCb::ParticleProperty* KInfo  = m_ppSvc->find( "K+" );

  m_pMass = pInfo->mass();
  m_piMass = piInfo->mass();
  m_KMass = KInfo->mass();

  m_BID = BInfo->particleID();
  m_BbarID = BInfo->antiParticle()->particleID();
  m_LbID = LbInfo->particleID();
  m_LbbarID = LbInfo->antiParticle()->particleID();
  m_pID = pInfo->particleID();
  m_pbarID = pInfo->antiParticle()->particleID();
  m_pimID = piInfo->antiParticle()->particleID();
  m_pipID = piInfo->particleID();
  m_KmID =  KInfo->antiParticle()->particleID();
  m_KpID = KInfo->particleID();

  return StatusCode::SUCCESS;
}
//=============================================================================
// StatusCode TrackStateProvider::finalize()
// {
//   m_reFit.release().ignore() ;
//   m_ppSvc.release().ignore() ;
//   return GaudiTool::finalize();
// }
//=============================================================================
StatusCode BsLpKPIDSubstitutionTool::doCorrection( LHCb::Particle* Part) 
{
  debug() << "Start Substitution" << endmsg;
  LHCb::Particle::ConstVector PartDau = Part->daughtersVector();
  if(msgLevel(MSG::DEBUG)){
    LHCb::Particle::ConstVector::const_iterator i;
      for( i = PartDau.begin(); i != PartDau.end(); ++i)
	debug() << "PID of daughter " << (*i)->particleID().pid() << endmsg;
  }
  if(Part->particleID() == m_LbID){
    debug() << "Substitution of Lb-Dau " << PartDau.at(1)->particleID().pid() << endmsg;
    const Gaudi::LorentzVector& oldMom1 = PartDau.at(1)->momentum() ;
    const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_pMass*m_pMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setParticleID(m_pbarID); 
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
										     oldMom1.Py() ,
										     oldMom1.Pz() , 
										     newEnergy1  ) );
    debug() << "Set Mass " << endmsg;
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMeasuredMass( m_pMass );

    const Gaudi::LorentzVector& oldMom0 = PartDau.at(0)->momentum() ;
    const double newEnergy0 = std::sqrt ( oldMom0.P2() + m_KMass*m_KMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setParticleID(m_KpID); 
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMomentum( Gaudi::LorentzVector( oldMom0.Px() ,
										     oldMom0.Py() ,
										     oldMom0.Pz() , 
										     newEnergy0  ) );
    debug() << "Set Mass " << endmsg;
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMeasuredMass( m_KMass );

  }
  if(Part->particleID() == m_LbbarID){
    debug() << "Substitution of Lbbar-Dau " << PartDau.at(0)->particleID().pid() << endmsg;
    const Gaudi::LorentzVector& oldMom0 = PartDau.at(0)->momentum() ;
    debug() << "Momentum vector " << &oldMom0 << endmsg;
    const double newEnergy0 = std::sqrt ( oldMom0.P2() + m_pMass*m_pMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setParticleID(m_pID); 
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMomentum( Gaudi::LorentzVector( oldMom0.Px() ,
										     oldMom0.Py() ,
										     oldMom0.Pz() , 
										     newEnergy0  ) );
    debug() << "Set Mass" << endmsg;
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMeasuredMass( m_pMass );

    const Gaudi::LorentzVector& oldMom1 = PartDau.at(1)->momentum() ;
    debug() << "Momentum vector " << &oldMom1 << endmsg;
    const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_KMass*m_KMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setParticleID(m_KmID); 
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
										     oldMom1.Py() ,
										     oldMom1.Pz() , 
										     newEnergy1  ) );
    debug() << "Set Mass" << endmsg;
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMeasuredMass( m_KMass );
  }
  debug() << "Substitution of " << Part->particleID().pid() << endmsg;
  Part->setParticleID( ( Part->particleID().pid() < 0 ? m_BID : m_BbarID ) );
  debug() << "Refit of " << Part->particleID().pid() << endmsg;
  StatusCode sc = m_reFit -> reFit(*Part);
  return sc;
}
//=============================================================================
StatusCode BsLpKPIDSubstitutionTool::undoCorrection( LHCb::Particle* Part) 
{
  debug() << "Start Re-Substitution" << endmsg;
  LHCb::Particle::ConstVector PartDau = Part->daughtersVector();
  if(msgLevel(MSG::DEBUG)){
    LHCb::Particle::ConstVector::const_iterator i;
    for( i = PartDau.begin(); i != PartDau.end(); ++i)
      debug() << "PID of daughter " << (*i)->particleID().pid() << endmsg;
  }
  if(Part->particleID() == m_BbarID){
    const Gaudi::LorentzVector& oldMom1 = PartDau.at(1)->momentum() ;
    const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_piMass*m_piMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setParticleID(m_pimID); 
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
										     oldMom1.Py() ,
										     oldMom1.Pz() , 
										     newEnergy1  ) );
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMeasuredMass( m_piMass );

    const Gaudi::LorentzVector& oldMom0 = PartDau.at(0)->momentum() ;
    const double newEnergy0 = std::sqrt ( oldMom0.P2() + m_piMass*m_piMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setParticleID(m_pipID); 
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMomentum( Gaudi::LorentzVector( oldMom0.Px() ,
										     oldMom0.Py() ,
										     oldMom0.Pz() , 
										     newEnergy0  ) );
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMeasuredMass( m_piMass );
  }
  if(Part->particleID() == m_BID){
    const Gaudi::LorentzVector& oldMom0 = PartDau.at(0)->momentum() ;
    const double newEnergy0 = std::sqrt ( oldMom0.P2() + m_piMass*m_piMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setParticleID(m_pipID); 
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMomentum( Gaudi::LorentzVector( oldMom0.Px() ,
										     oldMom0.Py() ,
										     oldMom0.Pz() , 
										     newEnergy0  ) );
    const_cast<LHCb::Particle*>(PartDau.at(0)) -> setMeasuredMass( m_piMass );

    const Gaudi::LorentzVector& oldMom1 = PartDau.at(1)->momentum() ;
    const double newEnergy1 = std::sqrt ( oldMom1.P2() + m_piMass*m_piMass ) ;
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setParticleID(m_pimID); 
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMomentum( Gaudi::LorentzVector( oldMom1.Px() ,
										     oldMom1.Py() ,
										     oldMom1.Pz() , 
										     newEnergy1  ) );
    const_cast<LHCb::Particle*>(PartDau.at(1)) -> setMeasuredMass( m_piMass );
  }
  debug() << "Substitution of " << Part->particleID().pid() << endmsg;
  Part->setParticleID( ( Part->particleID().pid() < 0 ? m_LbID : m_LbbarID ) );
  debug() << "Refit of " << Part->particleID().pid() << endmsg;
  StatusCode sc = m_reFit -> reFit(*Part);
  return sc;
}
