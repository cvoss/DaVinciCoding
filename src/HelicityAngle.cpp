// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "HelicityAngle.h"

//-----------------------------------------------------------------------------
// Implementation file for class : HelicityAngle
//
// 2016-01-25 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( HelicityAngle )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
HelicityAngle::HelicityAngle( const std::string& type,
                              const std::string& name,
                              const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IPhysicsComputation>(this);

}
//=============================================================================
// Destructor
//=============================================================================
HelicityAngle::~HelicityAngle() {} 

//=============================================================================
StatusCode HelicityAngle::initPhysics(std::vector<TLorentzVector> cands)
{
  //Per definition
  //cands.at(0) -> Grandmother
  //cands.at(1) -> mother
  //cands.at(2) -> daughter 

  m_grandparent = TLorentzVector( cands.at(0).Px(), cands.at(0).Py(), cands.at(0).Pz(), cands.at(0).E() );
  m_parent      = TLorentzVector( cands.at(1).Px(), cands.at(1).Py(), cands.at(1).Pz(), cands.at(1).E() );
  m_particle    = TLorentzVector( cands.at(2).Px(), cands.at(2).Py(), cands.at(2).Pz(), cands.at(2).E() );

  // m_grandparent.Print("v");
  // m_parent.Print("v");
  // m_particle.Print("v");
  debug() << "Grandparent " <<m_grandparent.M() << " Parent " << m_parent.M() << " Daughter " << m_particle.M() << endmsg;

  if ( msgLevel(MSG::INFO) )
  {
    info() << "Parent 4Vector lab frame" << endmsg;
    m_parent.Print("");
    info() << "-----------------------" << endmsg;
    info() << "Grandparent 4Vector lab frame" << endmsg;
    m_grandparent.Print("");
    info() << "-----------------------" << endmsg;
    info() << "Daughter 4Vector lab frame" << endmsg;
    m_particle.Print("");
  }
  
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode HelicityAngle::initPhysics(const LHCb::Particle::ConstVector cands)
{
  //Per definition
  //cands.at(0) -> Grandmother
  //cands.at(1) -> mother
  //cands.at(2) -> daughter
  debug() << cands.size() << endmsg;
    
  return StatusCode::SUCCESS;
} 
//=============================================================================
StatusCode HelicityAngle::computePhysics(std::vector<double> *vars)
{
  debug() << "Parent 4Vector" << endmsg;
  if ( msgLevel(MSG::DEBUG) )  m_parent.Print("v");
  debug() << "Parent Boost Vector" << endmsg;
  if ( msgLevel(MSG::DEBUG) ) m_parent.BoostVector().Print("v");

  TVector3 LzBoost(0.,0.,0.);
  LzBoost.SetX(m_parent.BoostVector().X());
  LzBoost.SetY(m_parent.BoostVector().Y());
  LzBoost.SetZ(m_parent.BoostVector().Z());
  
  m_particle.Boost(-LzBoost);
  m_parent.Boost(-LzBoost);
  m_grandparent.Boost(-LzBoost);
  
  if ( msgLevel(MSG::INFO) )
  {
    info() << "Parent 4Vector parent frame" << endmsg;
    m_parent.Print("");
    info() << "-----------------------" << endmsg;
    info() << "Grandparent 4Vector parent frame" << endmsg;
    m_grandparent.Print("");
    info() << "-----------------------" << endmsg;
    info() << "Daughter 4Vector parent frame" << endmsg;
    m_particle.Print("");
  }
  info() << "BoostedMass " << m_particle.M() << " BoostedMomentum " << m_particle.Vect().Mag() << endmsg;
  info() << "BoostedMass " << m_parent.M() << " BoostedMomentum " << m_parent.Vect().Mag() << endmsg;
  info() << "BoostedMass " << m_grandparent.M() << " BoostedMomentum " << m_grandparent.Vect().Mag() << endmsg;
  
  TVector3 particle3 = m_particle.Vect();
  TVector3 grandparent3 = -m_grandparent.Vect();

  debug() << "PartMom " << particle3.Mag() << " GrandMotherMom " << grandparent3.Mag() << endmsg;
  
  vars->at(0) = particle3.Angle(grandparent3);
  vars->at(1) = cos(particle3.Angle(grandparent3));
  
  return StatusCode::SUCCESS;
}
