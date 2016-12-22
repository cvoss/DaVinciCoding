// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "SquareDalitz.h"

//-----------------------------------------------------------------------------
// Implementation file for class : SquareDalitz
//
// 2014-07-28 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( SquareDalitz )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
SquareDalitz::SquareDalitz( const std::string& name,
                            ISvcLocator* pSvcLocator)
  : DaVinciHistoAlgorithm ( name , pSvcLocator )
{

}
//=============================================================================
// Destructor
//=============================================================================
SquareDalitz::~SquareDalitz() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode SquareDalitz::initialize() {
  StatusCode sc = DaVinciHistoAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  AngularDalitz = tool<IPhysicsComputation>("AngularDalitz",this);
  
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  const LHCb::ParticleProperty* LambdaInfo = ppSvc()->find( "Lambda0" ) ;
  const LHCb::ParticleProperty* BInfo = ppSvc()->find( "B0" ) ;
  const LHCb::ParticleProperty* pInfo = ppSvc()->find( "p+" );
  const LHCb::ParticleProperty* pipInfo = ppSvc()->find( "pi+" );
  const LHCb::ParticleProperty* pimInfo = ppSvc()->find( "pi-" );
    
  m_BID = BInfo->particleID();
  m_BbarID = BInfo->antiParticle()->particleID();
  m_LID = LambdaInfo->particleID();
  m_LbarID = LambdaInfo->antiParticle()->particleID();
  m_pID = pInfo->particleID();
  m_pbarID = pInfo->antiParticle()->particleID();
  m_pipID = pipInfo->particleID();
  m_pimID = pimInfo->particleID();
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode SquareDalitz::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  const LHCb::MCParticles* parts = getIfExists<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );

  info() << parts->size() << endmsg;
  
  // const Gaudi::LorentzVector B(0), L(0), p(0), pi(0);
   LHCb::MCParticles::const_iterator iMC;
  for(iMC = parts->begin(); iMC != parts->end(); ++iMC)
  {
    if ( (*iMC)->particleID().abspid() == 511 ) debug() << "MC Top Part ID " << (*iMC)->particleID().pid() 
                                                       << " Oscilated " << (*iMC)->hasOscillated()  << endmsg;
    
    TLorentzVector B(0,0,0,0), L(0,0,0,0), p(0,0,0,0), pi(0,0,0,0);
  
    if( (*iMC)->particleID().abspid() == m_BbarID.abspid() )
    {
      debug() << "##################################################" << endmsg;
      B.SetPxPyPzE( (*(*iMC)).momentum().px(),(*(*iMC)).momentum().py(), (*(*iMC)).momentum().pz(),(*(*iMC)).momentum().E() );
      SmartRefVector<LHCb::MCVertex>::const_iterator vi;
      for( vi = (*iMC)->endVertices().begin(); vi != (*iMC)->endVertices().end(); vi++ )
      {
        SmartRefVector<LHCb::MCParticle>::const_iterator ip;
        for( ip = (*vi)->products().begin(); ip != (*vi)->products().end(); ip++ )
        {
          debug() << "ID " << (*ip)->particleID().pid() << endmsg; 
          if( (*ip)->particleID() == m_LID )    L.SetPxPyPzE( (*(*ip)).momentum().px(),
                                                              (*(*ip)).momentum().py(), 
                                                              (*(*ip)).momentum().pz(),
                                                              (*(*ip)).momentum().E() );
          if( (*ip)->particleID() == m_pbarID ) p.SetPxPyPzE( (*(*ip)).momentum().px(),
                                                              (*(*ip)).momentum().py(), 
                                                              (*(*ip)).momentum().pz(),
                                                              (*(*ip)).momentum().E() );
          if( (*ip)->particleID() == m_pipID ) pi.SetPxPyPzE( (*(*ip)).momentum().px(),
                                                              (*(*ip)).momentum().py(), 
                                                              (*(*ip)).momentum().pz(),
                                                              (*(*ip)).momentum().E() );
        }
        debug() << "Vtx done" << endmsg;
        debug() << "##################################################" << endmsg;
      }
    }
    if( (*iMC)->particleID().abspid() == m_BID.abspid() )
    {
      debug() << "##################################################" << endmsg;
      B.SetPxPyPzE( (*(*iMC)).momentum().px(),(*(*iMC)).momentum().py(), (*(*iMC)).momentum().pz(),(*(*iMC)).momentum().E() );
      SmartRefVector<LHCb::MCVertex>::const_iterator vi;
      for( vi = (*iMC)->endVertices().begin(); vi != (*iMC)->endVertices().end(); vi++ )
      {
        SmartRefVector<LHCb::MCParticle>::const_iterator ip;
        for( ip = (*vi)->products().begin(); ip != (*vi)->products().end(); ip++ )
        {
          debug() << "ID " << (*ip)->particleID().pid() << endmsg;   
          if( (*ip)->particleID() == m_LbarID )    L.SetPxPyPzE( (*(*ip)).momentum().px(),
                                                              (*(*ip)).momentum().py(), 
                                                              (*(*ip)).momentum().pz(),
                                                              (*(*ip)).momentum().E() );
          if( (*ip)->particleID() == m_pID )       p.SetPxPyPzE( (*(*ip)).momentum().px(),
                                                              (*(*ip)).momentum().py(), 
                                                              (*(*ip)).momentum().pz(),
                                                              (*(*ip)).momentum().E() );
          if( (*ip)->particleID() == m_pimID )    pi.SetPxPyPzE( (*(*ip)).momentum().px(),
                                                              (*(*ip)).momentum().py(), 
                                                              (*(*ip)).momentum().pz(),
                                                              (*(*ip)).momentum().E() );
        }
        debug() << "Vtx done" << endmsg;
        debug() << "##################################################" << endmsg;
      }
    }
    if ( B.M()>0 && L.M()>0 && p.M()>0 && pi.M()>0)
    {
      std::vector<TLorentzVector> Cands;
      Cands.push_back(B);
      Cands.push_back(L);
      Cands.push_back(p);
      Cands.push_back(pi);
      
      std::vector<double> LDAngleVars(3);
      
      AngularDalitz->initPhysics( Cands );
      AngularDalitz->computePhysics( &LDAngleVars );
      
      plot(LDAngleVars.at(0),"cosTheta","",-1.1,1.1);
      plot(LDAngleVars.at(1),"PhiL","PhiL",-3.2,3.2);
      plot(LDAngleVars.at(2),"Phihh","Phihh",-3.2,3.2);
      if( L.Theta() > 0.1 && L.Theta() < 0.2 
          && p.Theta() > 0.1 && p.Theta() < 0.2 
          && pi.Theta() > 0.1 && pi.Theta() < 0.2 
          )
      {
        plot(LDAngleVars.at(0),"cosThetaCut","",-1.1,1.1);
        plot(LDAngleVars.at(1),"PhiLCut","PhiL",-3.2,3.2);
        plot(LDAngleVars.at(2),"PhihhCut","Phihh",-3.2,3.2);
      }
        
    }
    
  }
  
  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode SquareDalitz::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciHistoAlgorithm::finalize();  // must be called after all other actions
}

//=============================================================================
// Calculate Helicity angles
//=============================================================================
Double_t SquareDalitz::coshel(TLorentzVector particle,
                              TLorentzVector parent, 
                              TLorentzVector grandparent) 
{
  TVector3 boosttoparent = -(parent.BoostVector());
  
  particle.Boost(boosttoparent);
  grandparent.Boost(boosttoparent);
  
  TVector3 particle3 = particle.Vect();
  TVector3 grandparent3 = grandparent.Vect();
  
  Float_t numerator = particle3.Dot(grandparent3);
  Float_t denominator = (particle3.Mag())*(grandparent3.Mag());
  Float_t temp = numerator/denominator;
  
  return temp;

}
//=============================================================================
