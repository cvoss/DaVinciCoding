// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "BdToBsSubstitution.h"

//-----------------------------------------------------------------------------
// Implementation file for class : BdToBsSubstitution
//
// 2013-10-28 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( BdToBsSubstitution )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
BdToBsSubstitution::BdToBsSubstitution( const std::string& type,
                                        const std::string& name,
                                        const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IParticleManipulator>(this);

}
//=============================================================================
// Destructor
//=============================================================================
BdToBsSubstitution::~BdToBsSubstitution() {} 

//=============================================================================
//=============================================================================
StatusCode BdToBsSubstitution::initialize()
{
  m_ppSvc = svc<LHCb::IParticlePropertySvc>("LHCb::ParticlePropertySvc",true);

  const LHCb::ParticleProperty* BInfo  = m_ppSvc->find( "B0" ) ;
  const LHCb::ParticleProperty* BsInfo  = m_ppSvc->find( "B_s0" ) ;

  m_deltaMass = BsInfo->mass() - BInfo->mass();

  m_BID     = BInfo->particleID();
  m_BbarID  = BInfo->antiParticle()->particleID();
  m_BsID    = BsInfo->particleID();
  m_BsbarID = BsInfo->antiParticle()->particleID();

  return StatusCode::SUCCESS;
}
//=============================================================================
// StatusCode BdToBsSubstitution::finalize()
// {
//   m_reFit.release().ignore() ;
//   m_ppSvc.release().ignore() ;
//   return GaudiTool::finalize();
// }
//=============================================================================
StatusCode BdToBsSubstitution::doCorrection( LHCb::Particle* Part) 
{  
  Part->setMeasuredMass(Part->measuredMass() + m_deltaMass );
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode BdToBsSubstitution::undoCorrection( LHCb::Particle* Part) 
{
  Part->setMeasuredMass(Part->measuredMass() - m_deltaMass );
  return StatusCode::SUCCESS;
}
