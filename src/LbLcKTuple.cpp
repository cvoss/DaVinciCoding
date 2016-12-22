// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "LbLcKTuple.h"

//-----------------------------------------------------------------------------
// Implementation file for class : LbLcKTuple
//
// 2016-07-11 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( LbLcKTuple )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
LbLcKTuple::LbLcKTuple( const std::string& name,
                        ISvcLocator* pSvcLocator)
  : MyDVAlgorithm ( name , pSvcLocator )
{

}
//=============================================================================
// Destructor
//=============================================================================
LbLcKTuple::~LbLcKTuple() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode LbLcKTuple::initialize() {
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode LbLcKTuple::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;
  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.

  debug() << "Getting Particles" << endmsg;
  const LHCb::Particle::Range Bs = this->particles(); 
  const LHCb::RecVertex::Range _pvs = this->primaryVertices();
  
  for ( LHCb::Particle::ConstVector::const_iterator iB = Bs.begin() ; iB != Bs.end() ; ++iB ){

    LHCb::Particle::ConstVector daus = (*iB)->daughtersVector();
    Tuple tuple = nTuple( "LbLcK" ) ; 
    _CandInfo("Xb", *(*iB), _pvs, tuple, false);      
    _CandInfo("LC", *(daus.at(0)), _pvs, tuple, false);
    _CandInfo("K", *(daus.at(1)), _pvs, tuple, false);
    tuple->write();
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode LbLcKTuple::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize();  // must be called after all other actions
}

//=============================================================================
