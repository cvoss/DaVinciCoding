// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "WriteGenMCKinematics.h"

//-----------------------------------------------------------------------------
// Implementation file for class : WriteGenMCKinematics
//
// 2014-07-29 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( WriteGenMCKinematics )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
WriteGenMCKinematics::WriteGenMCKinematics( const std::string& name,
                                            ISvcLocator* pSvcLocator)
  : DaVinciTupleAlgorithm ( name , pSvcLocator )
{
  
}
//=============================================================================
// Destructor
//=============================================================================
WriteGenMCKinematics::~WriteGenMCKinematics() {
 
} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode WriteGenMCKinematics::initialize() {
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  const LHCb::ParticleProperty* LambdaInfo = ppSvc()->find( "Lambda0" ) ;
  const LHCb::ParticleProperty* BInfo = ppSvc()->find( "B0" ) ;
  const LHCb::ParticleProperty* pInfo = ppSvc()->find( "p+" );
  const LHCb::ParticleProperty* pipInfo = ppSvc()->find( "pi+" );
  const LHCb::ParticleProperty* pimInfo = ppSvc()->find( "pi-" );
  const LHCb::ParticleProperty* gammaInfo = ppSvc()->find( LHCb::ParticleID(22) );
    
  m_BID = BInfo->particleID();
  m_BbarID = BInfo->antiParticle()->particleID();
  m_LID = LambdaInfo->particleID();
  m_LbarID = LambdaInfo->antiParticle()->particleID();
  m_pID = pInfo->particleID();
  m_pbarID = pInfo->antiParticle()->particleID();
  m_pipID = pipInfo->particleID();
  m_pimID = pimInfo->particleID();
  m_gammaID = gammaInfo->particleID();
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode WriteGenMCKinematics::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  const LHCb::MCParticles* parts = getIfExists<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );
  LHCb::MCParticles::const_iterator iMC;
  
  for(iMC = parts->begin(); iMC != parts->end(); ++iMC)
  {
    Tuple tuple = nTuple( "BLppiTruth" ) ;
    TLorentzVector B(0,0,0,0), L(0,0,0,0), p(0,0,0,0), pi(0,0,0,0), Lp(0,0,0,0), Lpi(0,0,0,0);
    TVector3 PV(0,0,0), BV(0,0,0),LV(0,0,0);
    int BzTrue(0), BzbTrue(0), LzTrue(0), LzbarTrue(0), pTrue(0), pbarTrue(0), pimTrue(0), pipTrue(0);
    int NBzVtx(0), NBzbVtx(0);
    
    if( (*iMC)->particleID().abspid() == m_BID.abspid() )
    {
      BzbTrue = 1;
      B.SetPxPyPzE( (*(*iMC)).momentum().px(),(*(*iMC)).momentum().py(), (*(*iMC)).momentum().pz(),(*(*iMC)).momentum().E() );
      PV.SetXYZ((*(*iMC)).primaryVertex()->position().x(),(*(*iMC)).primaryVertex()->position().y(),(*(*iMC)).primaryVertex()->position().z());
      Gaudi::LorentzVector pi_Temp(0.,0.,0.,0.);   
      SmartRefVector<LHCb::MCVertex>::const_iterator vi;
      NBzbVtx =  (*iMC)->endVertices().size() ;

      for( vi = (*iMC)->endVertices().begin(); vi != (*iMC)->endVertices().end(); vi++ )
      {
        SmartRefVector<LHCb::MCParticle>::const_iterator ip;
        info() << "Number of B daugters " << (*vi)->products().size() << endmsg;
        for( ip = (*vi)->products().begin(); ip != (*vi)->products().end(); ip++ )
        {
          debug() << "Particle ID " << (*ip)->particleID().pid() << endmsg;
                    
          if( (*ip)->particleID() == m_LID ) 
          {
            LzTrue = 1;
            BV.SetXYZ((*(*ip)).originVertex()->position().x(),(*(*ip)).originVertex()->position().y(),(*(*ip)).originVertex()->position().z());
            L.SetPxPyPzE( (*(*ip)).momentum().px(), (*(*ip)).momentum().py(), (*(*ip)).momentum().pz(), (*(*ip)).momentum().E() );
            SmartRefVector<LHCb::MCVertex>::const_iterator vL;
            for( vL = (*ip)->endVertices().begin(); vL != (*ip)->endVertices().end(); vL++ )
            {
              SmartRefVector<LHCb::MCParticle>::const_iterator iL;
              for( iL = (*vL)->products().begin(); iL != (*vL)->products().end(); iL++ )
              {
                if( (*iL)->particleID() == m_pID ){
                  Lp.SetPxPyPzE( (*(*iL)).momentum().px(), (*(*iL)).momentum().py(), 
                                 (*(*iL)).momentum().pz(), (*(*iL)).momentum().E() );
                  LV.SetXYZ((*(*iL)).originVertex()->position().x(),(*(*iL)).originVertex()->position().y(),(*(*iL)).originVertex()->position().z());
                }
                
                if( (*iL)->particleID() == m_pimID )    Lpi.SetPxPyPzE( (*(*iL)).momentum().px(), (*(*iL)).momentum().py(), 
                                                                        (*(*iL)).momentum().pz(), (*(*iL)).momentum().E() );
              }
            }
          }
          if( (*ip)->particleID() == m_pbarID )
          {
            pbarTrue = 1;
            p.SetPxPyPzE( (*(*ip)).momentum().px(), (*(*ip)).momentum().py(), (*(*ip)).momentum().pz(), (*(*ip)).momentum().E() );
          }
          if( (*ip)->particleID() == m_pipID )
          {
            pipTrue = 1;
            pi_Temp = pi_Temp + (*ip)->momentum();
          }
          if( (*ip)->particleID() == m_gammaID )
            pi_Temp = pi_Temp + (*ip)->momentum();
        }
      }
      pi.SetPxPyPzE( pi_Temp.px(), pi_Temp.py(), pi_Temp.pz(), pi_Temp.E() );
    }

    if( (*iMC)->particleID().abspid() == m_BbarID.abspid() )
    {
      BzTrue = 1;
      B.SetPxPyPzE( (*(*iMC)).momentum().px(),(*(*iMC)).momentum().py(), (*(*iMC)).momentum().pz(),(*(*iMC)).momentum().E() );
      PV.SetXYZ((*(*iMC)).primaryVertex()->position().x(),(*(*iMC)).primaryVertex()->position().y(),(*(*iMC)).primaryVertex()->position().z());
      Gaudi::LorentzVector pi_Temp(0.,0.,0.,0.);
      SmartRefVector<LHCb::MCVertex>::const_iterator vi;
      NBzVtx =  (*iMC)->endVertices().size() ;
      
      for( vi = (*iMC)->endVertices().begin(); vi != (*iMC)->endVertices().end(); vi++ )
      {
        SmartRefVector<LHCb::MCParticle>::const_iterator ip;
        for( ip = (*vi)->products().begin(); ip != (*vi)->products().end(); ip++ )
        {
          if( (*ip)->particleID() == m_LbarID ) 
          {
            LzbarTrue = 1;
            BV.SetXYZ((*(*ip)).originVertex()->position().x(),(*(*ip)).originVertex()->position().y(),(*(*ip)).originVertex()->position().z());
            L.SetPxPyPzE( (*(*ip)).momentum().px(), (*(*ip)).momentum().py(), 
                            (*(*ip)).momentum().pz(), (*(*ip)).momentum().E() );
            SmartRefVector<LHCb::MCVertex>::const_iterator vL;
            for( vL = (*ip)->endVertices().begin(); vL != (*ip)->endVertices().end(); vL++ )
            {
              SmartRefVector<LHCb::MCParticle>::const_iterator iL;
              for( iL = (*vL)->products().begin(); iL != (*vL)->products().end(); iL++ )
              {
                if( (*iL)->particleID() == m_pbarID )
                {
                  Lp.SetPxPyPzE( (*(*iL)).momentum().px(), (*(*iL)).momentum().py(), 
                                   (*(*iL)).momentum().pz(), (*(*iL)).momentum().E() );
                  LV.SetXYZ((*(*iL)).originVertex()->position().x(),(*(*iL)).originVertex()->position().y(),(*(*iL)).originVertex()->position().z()); 
                }
                
                if( (*iL)->particleID() == m_pipID )    Lpi.SetPxPyPzE( (*(*iL)).momentum().px(), (*(*iL)).momentum().py(), 
                                                                          (*(*iL)).momentum().pz(), (*(*iL)).momentum().E() );
              }
            }
          }
          if( (*ip)->particleID() == m_pID )
          {
            pTrue = 1;
            p.SetPxPyPzE( (*(*ip)).momentum().px(), (*(*ip)).momentum().py(), (*(*ip)).momentum().pz(), (*(*ip)).momentum().E() );
          }
          if( (*ip)->particleID() == m_pimID )
          {
            pimTrue = 1;
            pi_Temp = pi_Temp + (*ip)->momentum();
          }
          if( (*ip)->particleID() == m_gammaID )
            pi_Temp = pi_Temp + (*ip)->momentum();

        }
      }
      pi.SetPxPyPzE( pi_Temp.px(), pi_Temp.py(), pi_Temp.pz(), pi_Temp.E() );
    }    
    
    
    
    //if ( B.M()>0 && L.M()>0 && p.M()>0 && pi.M()>0 && Lp.M()>0 && Lpi.M()>0 )
    {
      WriteParticleTruth("B",B,tuple);
      WriteParticleTruth("L",L,tuple);
      WriteParticleTruth("p",p,tuple);
      WriteParticleTruth("pi",pi,tuple);
      WriteVertexTruth("PV",PV,tuple);
      WriteVertexTruth("BV",BV,tuple);
      WriteVertexTruth("LV",LV,tuple);
      tuple->column("BzTrue",BzTrue);
      tuple->column("BzbTrue",BzbTrue);
      tuple->column("LzTrue",LzTrue);
      tuple->column("LzbarTrue",LzbarTrue);
      tuple->column("pTrue",pTrue);
      tuple->column("pbarTrue",pbarTrue);
      tuple->column("pimTrue",pimTrue);
      tuple->column("pipTrue",pipTrue);

      tuple->write();
    }
          
  }
  
  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode WriteGenMCKinematics::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;
  
  return DaVinciTupleAlgorithm::finalize();  // must be called after all other actions
}
//============================================================================
StatusCode   WriteGenMCKinematics::setColumn(Tuple& tuple, std::string candname,
                                              std::string data, double value)
  {
    candname.append(data);
    tuple->column(candname,value);
    return StatusCode::SUCCESS;
  }
//=============================================================================
StatusCode WriteGenMCKinematics::WriteParticleTruth( std::string name, TLorentzVector& p4, Tuple& tuple)
{
  setColumn( tuple, name, "_Px", p4.Px() );
  setColumn( tuple, name, "_Py", p4.Py() );
  setColumn( tuple, name, "_Pz", p4.Pz() );
  setColumn( tuple, name, "_E", p4.E() );
  return StatusCode::SUCCESS;
} 
//=============================================================================
StatusCode WriteGenMCKinematics::WriteVertexTruth( std::string name, TVector3& r, Tuple& tuple)
{
  setColumn( tuple, name, "_X", r.X() );
  setColumn( tuple, name, "_Y", r.Y() );
  setColumn( tuple, name, "_Z", r.Z() );
  return StatusCode::SUCCESS;
}
