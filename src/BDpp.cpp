// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "BDpp.h"

//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/bind.hpp>

//#include "Kernel/ParticleFilters.h"

//using namespace boost::lambda;

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( BDpp )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
BDpp::BDpp( const std::string& name,
                        ISvcLocator* pSvcLocator)
  : DaVinciHistoAlgorithm ( name , pSvcLocator )
{
  declareProperty("useMC", _useMC = false);
  declareProperty("nBins", m_Bins = 150);
  declareProperty("PIDToolname", m_PIDToolname = "NoPIDSubstitution");
}
//=============================================================================
// Destructor
//=============================================================================
BDpp::~BDpp() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode BDpp::initialize() {
  StatusCode sc = DaVinciHistoAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  _ppsvc = ppSvc();
  //_assoc = tool<IParticle2MCAssociator>("MCMatchObjP2MCRelator");
  _PID = tool<IParticleManipulator>(m_PIDToolname,this);
 
  const LHCb::ParticleProperty* BInfo = _ppsvc->find( "B0" ) ;
  const LHCb::ParticleProperty* LbInfo = _ppsvc->find( "Lambda_b0" );
  const LHCb::ParticleProperty* LInfo = _ppsvc->find( "Lambda0" );
  const LHCb::ParticleProperty* pInfo = _ppsvc->find( "p+" );
  const LHCb::ParticleProperty* piInfo = _ppsvc->find( "pi+" );

  _BMass = BInfo->mass();
  _LbMass = LbInfo->mass();
  _piMass = piInfo->mass();
  _pMass = pInfo->mass();

  _BID = BInfo->particleID();
  _BbarID = BInfo->antiParticle()->particleID();
  _LID = LInfo->particleID();
  _LbID = LbInfo->particleID();
  _LbbarID = LbInfo->antiParticle()->particleID();
  _pID = pInfo->particleID();
  _pbarID = pInfo->antiParticle()->particleID();
  _pimID = piInfo->antiParticle()->particleID();
  _pipID = piInfo->particleID();

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode BDpp::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;
  StatusCode sc = StatusCode::SUCCESS ;

  debug() << "Getting Particles" << endmsg;
  const LHCb::Particle::Range Bs = this->particles(); 
  const LHCb::RecVertex::Range _pvs = this->primaryVertices();

  //LHCb::Particle *Lambda(0), *p(0), *pi(0);
  info() << Bs.size() << endmsg;

  for ( LHCb::Particle::ConstVector::const_iterator iB = Bs.begin() ; iB != Bs.end() ; ++iB ){

    const LHCb::Vertex* BPV = (const LHCb::Vertex*) bestPV((*iB));
    if(!BPV)
      return StatusCode::FAILURE;

    plot((*iB)->measuredMass(), "M_B", "", 4500, 6000*Gaudi::Units::MeV, m_Bins ); // MM
    // info() << (*iB)->particleID().abspid() << endmsg;
    _PID -> doCorrection(const_cast<LHCb::Particle*>(*iB));
    
    plot((*iB)->measuredMass(), "M_B_Sub", "", 4500*Gaudi::Units::MeV, 6000*Gaudi::Units::MeV, m_Bins); // MM
    _PID -> undoCorrection(const_cast<LHCb::Particle*>(*iB));
    plot((*iB)->measuredMass(), "M_B_reSub", "", 4500, 600*Gaudi::Units::MeV, m_Bins); // MM

  }
  
  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.

  return StatusCode::SUCCESS;
  
}

  

//=============================================================================
//  Finalize
//=============================================================================
StatusCode BDpp::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;
  return DaVinciHistoAlgorithm::finalize();  // must be called after all other actions
}

