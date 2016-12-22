// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MultCand.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MultCand
//
// 2015-02-05 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MultCand )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MultCand::MultCand( const std::string& name,
                    ISvcLocator* pSvcLocator)
  : DaVinciHistoAlgorithm ( name , pSvcLocator )
{
  declareProperty("PIDToolnames", m_PIDToolnames, "List of PID tools");
}
//=============================================================================
// Destructor
//=============================================================================
MultCand::~MultCand() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MultCand::initialize() {
  StatusCode sc = DaVinciHistoAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  for (unsigned int i = 0; i < m_PIDToolnames.size(); i++)
    m_PIDTools.push_back( tool<IParticleManipulator>( m_PIDToolnames.at(i),this) );

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MultCand::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  const LHCb::Particle::Range Bs = this->particles();

  LHCb::Particle::ConstVector::const_iterator IterB;
  
  std::vector<IParticleManipulator*>::iterator iTool;
  std::vector<IParticleManipulator*>::reverse_iterator riTool;

  LHCb::Particle::ConstVector::const_iterator Mult = select_randomly(Bs.begin(), Bs.end());
  warning() << "==> Size of Input " << Bs.size() << endmsg;
  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    warning() << "==> Iter Mass " << (*IterB)->measuredMass() << endmsg;
    warning() << "==> Iter Type " << (int)(*IterB)->daughtersVector().at(2)->daughtersVector().at(0)->proto()->track()->type() << endmsg;
  }
  
  warning() << "==> Chosen index " << (*Mult)->measuredMass() << endmsg;
  
  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    for ( iTool = m_PIDTools.begin(); iTool != m_PIDTools.end(); ++iTool)
      (*iTool) -> doCorrection(const_cast<LHCb::Particle*>(*IterB));

    if ( (*IterB)->measuredMass() > 5230. && (*IterB)->measuredMass() < 5340.)
      plot(Bs.size(), "CandsPerEvent", "", 0., 100., 100);
    plot((*IterB)->measuredMass(), "Mass", "", 5000., 5600., 100);
    
    for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
      (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
  }
  
  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MultCand::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciHistoAlgorithm::finalize();  // must be called after all other actions
}

//=============================================================================
