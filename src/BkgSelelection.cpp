// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "BkgSelelection.h"

//-----------------------------------------------------------------------------
// Implementation file for class : BkgSelelection
//
// 2014-11-17 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( BkgSelelection )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
BkgSelelection::BkgSelelection( const std::string& name,
                                ISvcLocator* pSvcLocator)
  : MySelectionAnalysis ( name , pSvcLocator )
{
  declareProperty("BMassLowEdge", m_BMassLowEdge = 5420.*Gaudi::Units::MeV);
  declareProperty("BMassHighEdge", m_BMassHighEdge = 5800.*Gaudi::Units::MeV);
  declareProperty("LMassLowEdge", m_LMassLowEdge = 1105.*Gaudi::Units::MeV);
  declareProperty("LMassHighEdge", m_LMassHighEdge = 1125.*Gaudi::Units::MeV);
  declareProperty("ProtonPID", m_p_pNNPIDCut = 0.1);
  declareProperty("PIDToolnames", m_PIDToolnames, "List of PID tools");
  declareProperty("FilterTrigger", m_TriggerFilter = true);
  declareProperty("SelectLL", m_SelectLL = true);
}
//=============================================================================
// Destructor
//=============================================================================
BkgSelelection::~BkgSelelection() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode BkgSelelection::initialize() {
  StatusCode sc = MySelectionAnalysis::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  for (unsigned int i = 0; i < m_PIDToolnames.size(); i++)
    m_PIDTools.push_back( tool<IParticleManipulator>( m_PIDToolnames.at(i),this) );

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode BkgSelelection::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.

  LHCb::Particle::Range Bs = this->particles();
  
  debug() << "Size of Inputs " << Bs.size() << endmsg;
  if ( 0 == Bs.size() )
    warning() << "Empty particle container" << endmsg;
  
  LHCb::Particle::ConstVector::const_iterator IterB;
  
  std::vector<IParticleManipulator*>::iterator iTool;
  std::vector<IParticleManipulator*>::reverse_iterator riTool;
  
  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    for ( iTool = m_PIDTools.begin(); iTool != m_PIDTools.end(); ++iTool)
      (*iTool) -> doCorrection(const_cast<LHCb::Particle*>(*IterB));
    
    if ( FilterTrigger( *IterB, m_TriggerFilter, "2012b" ) == StatusCode::SUCCESS )
      if (PreSelection( *IterB, m_BMassLowEdge, m_BMassHighEdge, m_LMassLowEdge, m_LMassHighEdge) == StatusCode::SUCCESS)
        if ( PerformFit( *IterB ) == StatusCode::SUCCESS )
          if ( ( m_SelectLL ? LLCleanUp() : DDCleanUp() ) == StatusCode::SUCCESS )
            if ( m_p_pNNPID > m_p_pNNPIDCut
                 && m_p->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5 &&
                 m_pi->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5)
            {
              counter("Candidates")++;
              info() << " Canidate found " << endmsg;
              cloneAndMarkTree(*IterB);
            }
    
    for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
      (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode BkgSelelection::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return MySelectionAnalysis::finalize();  // must be called after all other actions
}

//=============================================================================
