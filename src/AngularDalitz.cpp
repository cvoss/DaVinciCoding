// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "AngularDalitz.h"

//-----------------------------------------------------------------------------
// Implementation file for class : AngularDalitz
//
// 2014-07-31 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( AngularDalitz )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
AngularDalitz::AngularDalitz( const std::string& type,
                              const std::string& name,
                              const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IPhysicsComputation>(this);

}
//=============================================================================
// Destructor
//=============================================================================
AngularDalitz::~AngularDalitz() {} 

//=============================================================================
StatusCode AngularDalitz::initialize()
{
   
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode AngularDalitz::computePhysics(std::vector<double> *vars)
{
  vars->at(0) = m_Lambda.Dot(m_z);
  vars->at(1) = atan2( m_Lambda.Dot(m_y) , m_Lambda.Dot(m_x) );
  vars->at(2) = atan2( m_v.Dot(m_w), m_v.Dot(m_u) ); 

  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode AngularDalitz::initPhysics(const LHCb::Particle::ConstVector cands)
{
  TLorentzVector B(cands.at(0)->momentum().Px(), cands.at(0)->momentum().Py(),
                   cands.at(0)->momentum().Pz(), cands.at(0)->momentum().E());
  TLorentzVector beamAxis(0.,0.,4.*Gaudi::Units::TeV,4.*Gaudi::Units::TeV);
  beamAxis.Boost(-B.BoostVector() );
  TVector3 beamVect = beamAxis.Vect();
  
  m_z = B.Vect().Unit();
  m_x = ( beamVect.Cross(m_z) ).Unit();
  m_y = ( m_z.Cross(m_x) ).Unit();
  
  TLorentzVector L (cands.at(1)->momentum().Px(), cands.at(1)->momentum().Py(),
                    cands.at(1)->momentum().Pz(), cands.at(1)->momentum().E()); 
  L.Boost(-B.BoostVector() );
  TLorentzVector p (cands.at(2)->momentum().Px(), cands.at(2)->momentum().Py(),
                    cands.at(2)->momentum().Pz(), cands.at(2)->momentum().E());
  p.Boost(-B.BoostVector() );
  TLorentzVector pi(cands.at(3)->momentum().Px(), cands.at(3)->momentum().Py(),
                    cands.at(3)->momentum().Pz(), cands.at(3)->momentum().E());
  pi.Boost(-B.BoostVector() ); 
  
  m_Lambda = L.Vect().Unit();
  m_u = ( m_z.Cross(m_Lambda) ).Unit();
  m_v = ( p.Vect().Cross( pi.Vect() ) ).Unit();
  m_w = ( m_Lambda.Cross( m_u ) ).Unit();
  
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode AngularDalitz::initPhysics(const std::vector<TLorentzVector> cands)
{
  TLorentzVector B(cands.at(0));
  TLorentzVector beamAxis(0.,0.,4.*Gaudi::Units::TeV,4.*Gaudi::Units::TeV);
  beamAxis.Boost(-B.BoostVector() );
  TVector3 beamVect = beamAxis.Vect();
  
  m_x = B.Vect().Unit();
  m_z = ( beamVect.Cross(m_x) ).Unit();
  m_y = ( m_z.Cross(m_x) ).Unit();
  
  TLorentzVector L (cands.at(1)); 
  L.Boost(-B.BoostVector() );
  TLorentzVector p (cands.at(2));
  p.Boost(-B.BoostVector() );
  TLorentzVector pi(cands.at(3));
  pi.Boost(-B.BoostVector() );

  m_Lambda = L.Vect().Unit();
  m_u = ( m_x.Cross(m_Lambda) ).Unit();
  m_v = ( p.Vect().Cross( pi.Vect() ) ).Unit();
  m_w = ( m_Lambda.Cross( m_u ) ).Unit();
  
  return StatusCode::SUCCESS;
}
