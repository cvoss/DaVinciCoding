// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "AmenterosPodolanski.h"

//-----------------------------------------------------------------------------
// Implementation file for class : AmenterosPodolanski
//
// 2015-01-16 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( AmenterosPodolanski )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
AmenterosPodolanski::AmenterosPodolanski( const std::string& type,
                                          const std::string& name,
                                          const IInterface* parent )
: GaudiTool ( type, name , parent )
{
  declareInterface<IPhysicsComputation>(this);
  
}
//=============================================================================
// Destructor
//=============================================================================
AmenterosPodolanski::~AmenterosPodolanski() {} 

//=============================================================================
StatusCode AmenterosPodolanski::initialize()
{
   
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode AmenterosPodolanski::initPhysics(const LHCb::Particle::ConstVector cands)
{
  TLorentzVector L (cands.at(0)->momentum().Px(), cands.at(0)->momentum().Py(),
                    cands.at(0)->momentum().Pz(), cands.at(0)->momentum().E()); 
  TLorentzVector p (cands.at(1)->momentum().Px(), cands.at(1)->momentum().Py(),
                    cands.at(1)->momentum().Pz(), cands.at(1)->momentum().E());
  TLorentzVector pi(cands.at(2)->momentum().Px(), cands.at(2)->momentum().Py(),
                    cands.at(2)->momentum().Pz(), cands.at(2)->momentum().E());

  p0     = L.Vect();
  if ( cands.at(1)->charge() > 0)
  {
    pplus = p.Vect();
    pminus = pi.Vect();
  }
  if( cands.at(1)->charge() < 0. )
  {
    pminus = p.Vect();
    pplus = pi.Vect();
  }

  return StatusCode::SUCCESS;
} 
//=============================================================================
StatusCode AmenterosPodolanski::initPhysics(std::vector<TLorentzVector> cands)
{
  // TLorentzVector L (cands.at(0).Px(), cands.at(0).Py(),
  //                   cands.at(0).Pz(), cands.at(0).E()); 
  // TLorentzVector p (cands.at(1).Px(), cands.at(1).Py(),
  //                   cands.at(1).Pz(), cands.at(1).E());
  // TLorentzVector pi(cands.at(2).Px(), cands.at(2).Py(),
  //                   cands.at(2).Pz(), cands.at(2).E());

  // p0     = L.Vect();
  // if ( charges.at(1) > 0)
  // {
  //   pplus = p.Vect();
  //   pminus = pi.Vect();
  // }
  // if( charges.at(1) < 0. )
  // {
  //   pminus = p.Vect();
  //   pplus = pi.Vect();
  // }
  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode AmenterosPodolanski::computePhysics(std::vector<double> *vars)
{
  double phiplus  = pplus.Angle(p0);
  double plus_L   = pplus.Mag() * cos(phiplus);
  
  double phiminus = pminus.Angle(p0);
  double minus_L  = pminus.Mag() * cos(phiminus);

  double pt = pplus.Mag() * sin(phiplus);
  double aL = ( plus_L - minus_L ) / ( plus_L + minus_L );

  vars->at(0) = pt;
  vars->at(1) = aL;
  
  return StatusCode::SUCCESS;
}
