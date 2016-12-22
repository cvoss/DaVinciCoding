// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "WritePureKinematics.h"

//-----------------------------------------------------------------------------
// Implementation file for class : WritePureKinematics
//
// 2015-02-05 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( WritePureKinematics )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
WritePureKinematics::WritePureKinematics( const std::string& name,
                                          ISvcLocator* pSvcLocator)
  : DaVinciTupleAlgorithm ( name , pSvcLocator )
{
  declareProperty("PIDToolnames", m_PIDToolnames, "List of PID tools");
}
//=============================================================================
// Destructor
//=============================================================================
WritePureKinematics::~WritePureKinematics() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode WritePureKinematics::initialize() {
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;
  
  for (unsigned int i = 0; i < m_PIDToolnames.size(); i++)
    m_PIDTools.push_back( tool<IParticleManipulator>( m_PIDToolnames.at(i),this) );

  const LHCb::ParticleProperty* pInfo = ppSvc()->find( "p+" );

  m_pID = pInfo->particleID();
  m_pbarID = pInfo->antiParticle()->particleID();

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode WritePureKinematics::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  LHCb::Particle::Range Bs = this->particles();
  
  debug() << "Size of Cands" << Bs.size() << endmsg;
  
  LHCb::Particle::ConstVector::const_iterator IterB;
  
  std::vector<IParticleManipulator*>::iterator iTool;
  std::vector<IParticleManipulator*>::reverse_iterator riTool;

  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    for ( iTool = m_PIDTools.begin(); iTool != m_PIDTools.end(); ++iTool)
      (*iTool) -> doCorrection(const_cast<LHCb::Particle*>(*IterB));
    
    Tuple tuple = nTuple( "BLppiDPKin" ) ;
   
    debug() << "Get Tuple" << endmsg;
    

    WriteParticleKinematic("B",*(*IterB),tuple);
    WriteParticleKinematic("L",*((*IterB)->daughtersVector().at(2)),tuple);
    if((*IterB)->daughtersVector().at(0)->particleID() == m_pID )
    {
      WriteParticleKinematic("p",*((*IterB)->daughtersVector().at(0)),tuple);
      WriteParticleKinematic("pi",*((*IterB)->daughtersVector().at(1)),tuple);
    }
    
    if ((*IterB)->daughtersVector().at(1)->particleID() == m_pbarID )
    {
      WriteParticleKinematic("p",*((*IterB)->daughtersVector().at(1)),tuple);
      WriteParticleKinematic("pi",*((*IterB)->daughtersVector().at(0)),tuple);
    }
    
    for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
      (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
  }

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode WritePureKinematics::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize();  // must be called after all other actions
}
//============================================================================
StatusCode WritePureKinematics::setColumn(Tuple& tuple, std::string candname,
                                          std::string data, double value)
  {
    candname.append(data);
    tuple->column(candname,value);
    return StatusCode::SUCCESS;
  }
//=============================================================================
StatusCode WritePureKinematics::WriteParticleKinematic( std::string name, 
                                                        const LHCb::Particle& part, 
                                                        Tuple& tuple)
{
  setColumn( tuple, name, "PDG", (double)part.particleID().pid() );
  setColumn( tuple, name, "Charge", (double)part.charge() );
  setColumn( tuple, name, "pMag",part.p() );
  setColumn( tuple, name, "pt",part.pt() );  
  setColumn( tuple, name, "pX", part.momentum().px() );
  setColumn( tuple, name, "pY", part.momentum().py() );
  setColumn( tuple, name, "pZ", part.momentum().pz() );
  setColumn( tuple, name, "pE", part.momentum().E() );
  
  return StatusCode::SUCCESS;
} 
//=============================================================================
