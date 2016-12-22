// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "NoPIDSubstitution.h"

//-----------------------------------------------------------------------------
// Implementation file for class : NoPIDSubstitution
//
// 2013-10-25 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( NoPIDSubstitution )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
NoPIDSubstitution::NoPIDSubstitution( const std::string& type,
                                      const std::string& name,
                                      const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IParticleManipulator>(this);

}
//=============================================================================
// Destructor
//=============================================================================
NoPIDSubstitution::~NoPIDSubstitution() {} 

//=============================================================================
StatusCode NoPIDSubstitution::initialize(){
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode NoPIDSubstitution::doCorrection( LHCb::Particle* Part) {
  debug() << "ParticleID is " << Part->particleID().pid() << endmsg;
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode NoPIDSubstitution::undoCorrection( LHCb::Particle* Part) {
  debug() << "ParticleID is " << Part->particleID().pid() << endmsg;
  return StatusCode::SUCCESS;
}
