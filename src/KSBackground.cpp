// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "KSBackground.h"

//-----------------------------------------------------------------------------
// Implementation file for class : KSBackground
//
// 2014-06-23 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( KSBackground )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
KSBackground::KSBackground( const std::string& name,
                            ISvcLocator* pSvcLocator)
  : DaVinciTupleAlgorithm ( name , pSvcLocator )
  , NKS(0),NLZ(0)
{
  
}
//=============================================================================
// Destructor
//=============================================================================
KSBackground::~KSBackground() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode KSBackground::initialize() {
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  ArmenterosPodolanski = tool<IPhysicsComputation>("AmenterosPodolanski",this);


  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode KSBackground::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.

  LHCb::Particle::Range Bs = this->particles();
  
  if ( 0 == Bs.size() )
    warning() << "Empty particle container" << endmsg;

  LHCb::Particle::ConstVector::const_iterator IterB;
  
  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    const LHCb::Vertex* m_BPV = (const LHCb::Vertex*) bestPV((*IterB));

    double flength(0);
    distanceCalculator()->distance( *IterB, m_BPV, flength);

    TVector3 FlL, Lmom;
    FlL.SetX( - m_BPV->position().x() + (*IterB)->endVertex()->position().x() );
    FlL.SetY( - m_BPV->position().y() + (*IterB)->endVertex()->position().y() );
    FlL.SetZ( - m_BPV->position().z() + (*IterB)->endVertex()->position().z() );
    Lmom.SetX( (*IterB)->momentum().px() );
    Lmom.SetY( (*IterB)->momentum().py() );
    Lmom.SetZ( (*IterB)->momentum().pz() );
 
    LHCb::Particle::ConstVector APCands;
    std::vector<double> AP_vars(2);
    APCands.push_back(*IterB);
    APCands.push_back( (*IterB)->daughtersVector().at(0) );
    APCands.push_back( (*IterB)->daughtersVector().at(1) );
    
    ArmenterosPodolanski->initPhysics( APCands );
    ArmenterosPodolanski->computePhysics( &AP_vars );
    
    Tuple tuple = nTuple("KSTest");
    tuple->column("V0ID", (*IterB)->particleID().pid() );
    tuple->column("AP_pt", AP_vars.at(0));
    tuple->column("AP_alpha", AP_vars.at(1));
    tuple->column("Dira",Lmom.Angle(FlL));
    tuple->column("Mass", (*IterB)->measuredMass() );
    tuple->column("flength",flength);
    
    tuple->write();
      
      
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode KSBackground::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize();  // must be called after all other actions
}

//=============================================================================
