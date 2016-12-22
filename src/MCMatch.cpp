// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MCMatch.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MCMatch
//
// 2013-05-21 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MCMatch )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MCMatch::MCMatch( const std::string& name,
                  ISvcLocator* pSvcLocator)
  : DaVinciAlgorithm ( name , pSvcLocator )
{
  declareProperty("DoMatching", m_match = true);
  declareProperty("deepMatch", m_deepMatch = false);
}
//=============================================================================
// Destructor
//=============================================================================
MCMatch::~MCMatch() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MCMatch::initialize() {
  StatusCode sc = DaVinciAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  _assoc = tool<IParticle2MCAssociator>("MCMatchObjP2MCRelator",this);
  _PID = tool<IParticleManipulator>("PIDSubstitutionTool",this);

  info() << "Output written to " << outputLocation() << endmsg;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MCMatch::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  if (!m_match && m_deepMatch)
    warning() << "Configuration Error, Match needs to be true for DeepMatch to work" << endmsg;

  LHCb::Particle::Range Bs = this->particles();
  debug() << "Size of Inputs " << Bs.size() << endmsg;
  LHCb::Particle::ConstVector::const_iterator IterB;
  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    StatusCode s1 = _PID -> doCorrection(const_cast<LHCb::Particle*>(*IterB));
    if(m_match && m_deepMatch)
      if( MCMatch::_truthMatch(*(*IterB)) ) 
      {
        StatusCode s2 = _PID -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
        info() << "True B Found" << endmsg;
        cloneAndMarkTree(*IterB);
        counter("BCand")++;
        continue;
      }
    if(m_match && !m_deepMatch)
      if(_assoc->relatedMCP( (*IterB), LHCb::MCParticleLocation::Default))
        if( _assoc->relatedMCP( (*IterB), LHCb::MCParticleLocation::Default)->particleID().hasBottom() )
        {
          StatusCode s2 = _PID -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
          info() << "True B Found" << endmsg;
          cloneAndMarkTree(*IterB);
          counter("BCand")++;
          continue;
        }
    StatusCode s2 = _PID -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
    info() << "B Found" << endmsg;
    cloneAndMarkTree(*IterB);
    counter("BCand")++;
  }  
  
  //this->_saveInTES();
  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MCMatch::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciAlgorithm::finalize();  // must be called after all other actions
}

//=============================================================================

bool MCMatch::_truthMatch  (const LHCb::Particle& TopPart)
{
  const LHCb::MCParticle* MCTopPart = _assoc->relatedMCP(&TopPart, LHCb::MCParticleLocation::Default);
  if( !MCTopPart ) return false;
  if( MCTopPart->particleID() != TopPart.particleID() ) return false;
  debug() << "==> TopPartID " << MCTopPart->particleID().pid() << endmsg;
  LHCb::Particle::ConstVector daus = TopPart.daughtersVector();
  LHCb::Particle::ConstVector::const_iterator iDau;
  for( iDau = daus.begin(); iDau != daus.end(); ++iDau){
    if( !_truthMatch( *(*iDau) ) ) return false;
    if( _assoc->relatedMCP( (*iDau), LHCb::MCParticleLocation::Default)->mother()->particleID() != TopPart.particleID() )
      return false;
  }
  debug() << "==> True Decay found" << endmsg;
  return true;
} 
