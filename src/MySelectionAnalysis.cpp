// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MySelectionAnalysis.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MySelectionAnalysis
//
// 2014-05-19 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MySelectionAnalysis )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MySelectionAnalysis::MySelectionAnalysis( const std::string& name,
                                          ISvcLocator* pSvcLocator)
  : DaVinciAlgorithm ( name , pSvcLocator )
  , m_log_pi_IP(0), m_log_p_IP(0), m_log_Lpi_IP(0), m_log_Lp_IP(0)
  , m_log_BIP(0), m_log_LIP(0)
  , m_log_BVtxchi2(0) 
  , m_log_B_DauSumchi2(0), m_log_B_FLBchi2(0)
  , m_log_Lambda_FDchi2(0), m_B_Eta(0)
  , m_B_PT(0) , m_L_PT(0) , m_p_PT(0) , m_pi_PT(0) , m_Lp_PT(0) , m_Lpi_PT(0)
  , m_BWerner(0), m_LWerner(0), m_ppiDist(0)
  , m_LAngle(0), m_LMassVal(0), m_LP_PT(0), m_pointing(0)
  , m_BMassHighEdge(0), m_BMassLowEdge(0)
  , m_LMassHighEdge(0), m_LMassLowEdge(0)
  , m_B(0), m_Lambda(0), m_p(0), m_pi(0), m_Lp(0), m_Lpi(0)
{

}
//=============================================================================
// Destructor
//=============================================================================
MySelectionAnalysis::~MySelectionAnalysis() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MySelectionAnalysis::initialize() {
  StatusCode sc = DaVinciAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  m_ppsvc = ppSvc();
  m_states = tool<ITrackStateProvider>("TrackStateProvider",this);

  _L0TriggerTool     = tool<ITriggerTisTos>( "L0TriggerTisTos","L0TriggerTisTos",this);
  _Hlt1TriggerTool   = tool<ITriggerTisTos>( "Hlt1TriggerTisTos","Hlt1TriggerTisTos",this );
  _Hlt2TriggerTool   = tool<ITriggerTisTos>( "Hlt2TriggerTisTos","Hlt2TriggerTisTos",this );
  _HltTriggerTool    = tool<ITriggerTisTos>( "TriggerTisTos","TriggerTisTos",this);

  // copy setting of TriggerTisTos to settings of Hlt1, Hlt2 
  _Hlt1TriggerTool->setTOSFrac( _HltTriggerTool->getTOSFrac() );
  _Hlt1TriggerTool->setTISFrac( _HltTriggerTool->getTISFrac() );
  _Hlt2TriggerTool->setTOSFrac( _HltTriggerTool->getTOSFrac() );
  _Hlt2TriggerTool->setTISFrac( _HltTriggerTool->getTISFrac() );
   
  const LHCb::ParticleProperty* LambdaInfo = m_ppsvc->find( "Lambda0" ) ;
  const LHCb::ParticleProperty* BInfo      = m_ppsvc->find( "B0" ) ;
  const LHCb::ParticleProperty* BsInfo     = m_ppsvc->find( "B_s0" ) ;
  const LHCb::ParticleProperty* pInfo      = m_ppsvc->find( "p+" );
  const LHCb::ParticleProperty* pipInfo    = m_ppsvc->find( "pi+" );
  const LHCb::ParticleProperty* pimInfo    = m_ppsvc->find( "pi-" );
  const LHCb::ParticleProperty* KpInfo     = m_ppsvc->find( "K+" );  
  const LHCb::ParticleProperty* KmInfo     = m_ppsvc->find( "K-" ); 

  m_BID         = BInfo->particleID();
  m_BbarID      = BInfo->antiParticle()->particleID();
  m_BsbarID     = BsInfo->antiParticle()->particleID();
  m_LambdaID    = LambdaInfo->particleID();
  m_LambdabarID = LambdaInfo->antiParticle()->particleID();
  m_protonID    = pInfo->particleID();
  m_protonbarID = pInfo->antiParticle()->particleID();
  m_pimID       = pimInfo->particleID();
  m_pipID       = pipInfo->particleID();
  m_KpID        = KpInfo->particleID();
  m_KmID        = KmInfo->particleID();

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MySelectionAnalysis::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MySelectionAnalysis::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciAlgorithm::finalize();  // must be called after all other actions
}

//=============================================================================
// Select Lambda and BMass window
//=============================================================================
StatusCode MySelectionAnalysis::PreSelection(const LHCb::Particle* BCand,
                                             double BMassLowEdge, double BMassHighEdge,
                                             double LMassLowEdge, double LMassHighEdge)
{
  if( (BCand)->measuredMass() < BMassLowEdge || (BCand)->measuredMass() > BMassHighEdge )
  {
    debug() << "BMass Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
  if( (BCand)->daughtersVector().at(2) )
    if( (BCand)->daughtersVector().at(2)->measuredMass() < LMassLowEdge ||
        (BCand)->daughtersVector().at(2)->measuredMass() > LMassHighEdge )
    {
      debug() << "LMass Failure" << endmsg;
      return StatusCode::FAILURE;
    }
  
  m_LMassVal = (BCand)->daughtersVector().at(2)->measuredMass();
  
  m_BPV = (const LHCb::Vertex*) bestPV((BCand));

  if(!m_BPV)
  {
    debug() << "No PV Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
  TVector3 FlB0, Bp;
  FlB0.SetX( - m_BPV->position().x() + (BCand)->endVertex()->position().x() );
  FlB0.SetY( - m_BPV->position().y() + (BCand)->endVertex()->position().y() );
  FlB0.SetZ( - m_BPV->position().z() + (BCand)->endVertex()->position().z() );
  Bp.SetX( (BCand)->momentum().px() );
  Bp.SetY( (BCand)->momentum().py() );
  Bp.SetZ( (BCand)->momentum().pz() );

  m_BWerner = Bp.Angle(FlB0);
  
  if ( isnan(m_BWerner) )
  {
    info() << "NAN for B-Werner" << endmsg;
    return StatusCode::FAILURE;
  }
  if ( !MySelectionAnalysis::ThetaAngleInB(*BCand, *((BCand)->daughtersVector().at(2)), m_LAngle) )
  {
    debug() << "Calculate HeliAngle Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
  TVector3 FlL, Lmom;
  FlL.SetX( - (BCand)->endVertex()->position().x() + (BCand)->daughtersVector().at(2)->endVertex()->position().x());
  FlL.SetY( - (BCand)->endVertex()->position().y() + (BCand)->daughtersVector().at(2)->endVertex()->position().y());
  FlL.SetZ( - (BCand)->endVertex()->position().z() + (BCand)->daughtersVector().at(2)->endVertex()->position().z());
  Lmom.SetX( (BCand)->daughtersVector().at(2)->momentum().px() );
  Lmom.SetY( (BCand)->daughtersVector().at(2)->momentum().py() );
  Lmom.SetZ( (BCand)->daughtersVector().at(2)->momentum().pz() );
  
  double LWerner(-1e99);
  if ( BoostLWerner(*BCand, *(BCand)->daughtersVector().at(2), FlL, LWerner) ==  StatusCode::FAILURE)
  {
    debug() << "LWerner Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
  m_LWerner = log10(LWerner);

  debug() << "B Mass before DTF " << BCand->measuredMass() 
          << " and Lambda before DTF " << (BCand)->daughtersVector().at(2)->measuredMass() << endmsg;

  distanceCalculator()->distance( (BCand)->daughtersVector().at(0), (BCand)->daughtersVector().at(1), m_ppiDist);

  return StatusCode::SUCCESS;
}
//=============================================================================
// Perform DecayTreeFit
//=============================================================================
StatusCode MySelectionAnalysis::PerformFit(const LHCb::Particle* BCand)
{
  DecayTreeFitter::Fitter decayTreeFitter((*BCand), *m_BPV, m_states);
  decayTreeFitter.setMassConstraint ( m_LambdaID );
  decayTreeFitter.updateCand(*const_cast<LHCb::Particle*>(BCand));
  decayTreeFitter.fit();
  
  if (decayTreeFitter.status() != 0 ) 
  {
    debug() << "DTF Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
  if(decayTreeFitter.chiSquare()<= 0.)
  {
    debug() << "DTF-II Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
  LHCb::DecayTree fitted = decayTreeFitter.getFittedTree() ;
    
  m_B = fitted.head();if( !m_B ) return StatusCode::FAILURE;
  m_Lambda = fitted.find( m_B->particleID() == m_BbarID ? m_LambdaID    : m_LambdabarID );
  m_p      = fitted.find( m_B->particleID() == m_BbarID ? m_protonbarID : m_protonID );
  m_pi     = fitted.find( m_B->particleID() == m_BbarID ? m_pipID       : m_pimID );
  m_Lpi    = fitted.find( m_B->particleID() == m_BbarID ? m_pimID       : m_pipID );
  m_Lp     = fitted.find( m_B->particleID() == m_BbarID ? m_protonID    : m_protonbarID );
  if(m_B->particleID().abspid() == m_BsbarID.abspid())
  {
    m_Lambda = fitted.find( m_B->particleID() == m_BsbarID ? m_LambdaID    : m_LambdabarID );
    m_Lpi    = fitted.find( m_B->particleID() == m_BsbarID ? m_pimID       : m_pipID );
    m_Lp     = fitted.find( m_B->particleID() == m_BsbarID ? m_protonID    : m_protonbarID );
    m_p      = fitted.find( m_B->particleID() == m_BsbarID ? m_protonbarID : m_protonID );
    m_pi     = fitted.find( m_B->particleID() == m_BsbarID ? m_KpID       : m_KmID );
  }
  
  debug() << m_B << " " << m_Lambda << " " << m_p << " " << m_pi << " " << m_Lp << " " << m_Lpi << endmsg;

  debug() << "B Mass" << m_B->measuredMass() <<endmsg;
  debug() << "L Mass" << m_Lambda->measuredMass() <<endmsg;
  debug() << "Mass p+pi " << (m_p->momentum()+m_pi->momentum()).M();
  debug() << "Mass Lp+Lpi " << (m_Lp->momentum()+m_Lpi->momentum()).M() << endmsg;

  double BLF(0), BLFErr(0), LLF(0), LLFErr(0),Bctau(0),Lctau(0);
  const Gaudi::Math::ParticleParams* BParams = decayTreeFitter.fitParams(BCand);
  debug() << BParams << endmsg;
  if (BParams){
    Bctau = BParams->ctau().value();
    BLF = BParams->decayLength().value();
    BLFErr = BParams->decayLength().error();
  }
  const Gaudi::Math::ParticleParams* LParams = decayTreeFitter.fitParams(BCand->daughtersVector().at(2));
  debug() << LParams << endmsg;
  if (LParams){
    Lctau = LParams->ctau().value();
    LLF = LParams->decayLength().value();
    LLFErr = LParams->decayLength().error();
  }

  if( BLF <= 0. || LLF <= 0. || BLFErr<=0. || LLFErr<=0. || Bctau <= 0. || Lctau <= 0.)
  {
    debug() << "Calculate ctau Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
 m_log_B_ctau = log10(Bctau);
 m_log_L_ctau = log10(Lctau);
 m_log_BVtxchi2 =((decayTreeFitter.chiSquare()/decayTreeFitter.nDof()));
 m_FitProb =(TMath::Prob(decayTreeFitter.chiSquare(),decayTreeFitter.nDof()));
 m_log_B_FLBchi2 =(log10(BLF/BLFErr));
 m_log_Lambda_FDchi2 =(log10(LLF/LLFErr));

 debug() << "Finished Fit procedure" << endmsg;
 
 return StatusCode::SUCCESS;
}
//=============================================================================
// Calculate Common Variables
//=============================================================================
StatusCode MySelectionAnalysis::Calculate()
{
  debug() << "Begin Calculations" << endmsg;
  
  double p_pNNPID = m_p->proto()->info(LHCb::ProtoParticle::ProbNNp, -1e9) ;
  double p_piNNPID = m_p->proto()->info(LHCb::ProtoParticle::ProbNNpi, -1e9) ;
  double p_kNNPID = m_p->proto()->info(LHCb::ProtoParticle::ProbNNk, -1e9) ;

  double pi_pNNPID = m_pi->proto()->info(LHCb::ProtoParticle::ProbNNp, -1e9) ;
  double pi_piNNPID = m_pi->proto()->info(LHCb::ProtoParticle::ProbNNpi, -1e9) ;
  double pi_kNNPID = m_pi->proto()->info(LHCb::ProtoParticle::ProbNNk, -1e9) ;

  double Lp_pNNPID = m_Lp->proto()->info(LHCb::ProtoParticle::ProbNNp, -1e9) ;
  double Lp_pDLLPID = m_Lp->proto()->info(LHCb::ProtoParticle::RichDLLp, -1e9) ;
  double Lp_piNNPID = m_Lp->proto()->info(LHCb::ProtoParticle::ProbNNpi, -1e9) ;
  double Lp_kNNPID = m_Lp->proto()->info(LHCb::ProtoParticle::ProbNNk, -1e9) ;

  double Lpi_pNNPID = m_Lpi->proto()->info(LHCb::ProtoParticle::ProbNNp, -1e9) ;
  double Lpi_piNNPID = m_Lpi->proto()->info(LHCb::ProtoParticle::ProbNNpi, -1e9) ;
  double Lpi_kNNPID = m_Lpi->proto()->info(LHCb::ProtoParticle::ProbNNk, -1e9) ;

  m_p_pNNPID = ( p_pNNPID / ( p_pNNPID + p_piNNPID + p_kNNPID) );
  m_p_KNNPID = ( p_kNNPID / ( p_kNNPID + p_piNNPID + p_pNNPID) );
  m_p_piNNPID = ( p_piNNPID );

  m_pi_pNNPID = pi_pNNPID / ( pi_pNNPID + pi_piNNPID + pi_kNNPID ) ;
  m_pi_KNNPID = pi_kNNPID / ( pi_kNNPID + pi_piNNPID + pi_pNNPID ) ;
  m_pi_piNNPID = pi_piNNPID ;

  m_Lp_pNNPID = Lp_pNNPID / ( Lp_pNNPID + Lp_piNNPID + Lp_kNNPID) ;
  m_Lp_pDLLPID = Lp_pDLLPID;
  m_Lp_KNNPID = Lp_kNNPID / ( Lp_kNNPID + Lp_piNNPID + Lp_piNNPID) ;
  m_Lp_piNNPID = Lp_piNNPID ;
  
  m_Lpi_pNNPID = Lpi_pNNPID / ( Lpi_pNNPID + Lpi_piNNPID + Lpi_kNNPID) ;
  m_Lpi_KNNPID = Lpi_kNNPID / ( Lpi_kNNPID + Lpi_piNNPID + Lpi_pNNPID) ;
  m_Lpi_piNNPID = Lpi_piNNPID ;

  m_pi_KDLLPID = m_pi->proto()->info(LHCb::ProtoParticle::CombDLLk, -1e9);
  
  // if ( m_p_pNNPID < 0.05)
  //   return StatusCode::FAILURE;

  double BPVIP(0), LPVIP(0), pPVIP(0), piPVIP(0), LpPVIP(0), LpiPVIP(0);
  double BPVIPchi2(0), LPVIPchi2(0), pPVIPchi2(0), piPVIPchi2(0), LpPVIPchi2(0), LpiPVIPchi2(0);
  
  distanceCalculator()->distance( m_B, m_BPV, BPVIP, BPVIPchi2);
  distanceCalculator()->distance( (m_Lambda), m_BPV, LPVIP, LPVIPchi2);
  distanceCalculator()->distance( (m_p), m_BPV, pPVIP, pPVIPchi2);
  distanceCalculator()->distance( (m_pi), m_BPV, piPVIP, piPVIPchi2);
  distanceCalculator()->distance( (m_Lp), m_BPV, LpPVIP, LpPVIPchi2);
  distanceCalculator()->distance( (m_Lpi), m_BPV, LpiPVIP, LpiPVIPchi2);
  
  if( LpiPVIP<=0. || pPVIP<=0. || LpPVIP<=0. || piPVIP<=0. || LPVIP<=0.)
  {
    debug() << "Calculate IP Failure" << endmsg;
    return StatusCode::FAILURE;
  }
  
  m_log_pi_IP =(log10(piPVIP));
  m_log_p_IP =(log10(pPVIP));
  m_log_Lpi_IP =(log10(LpiPVIP));
  m_log_Lp_IP =(log10(LpPVIP));
  m_log_BIP = log10(BPVIP);
  m_log_LIP = log10(LPVIP);
  m_log_B_DauSumchi2 = log10(pPVIPchi2 + piPVIPchi2);
  
  m_B_Eta =( 0.5 * log( ( m_B->momentum().e() + m_B->p() ) / ( m_B->momentum().e() - m_B->p() ) ) );
  m_B_PT =((m_B->pt()));
  m_L_PT =((m_Lambda->pt()));
  m_p_PT =((m_p->pt()));
  m_pi_PT =((m_pi->pt()));
  m_Lp_PT =((m_Lp->pt()));
  m_Lpi_PT =((m_Lpi->pt()));

  m_LP_PT = ( m_Lambda->momentum() + m_p->momentum() ).Pt();
  
  double sump(0), sumpt(0);
  sumpt = m_Lambda->pt() + m_p->pt() + m_pi->pt() ;
  sump  = m_Lambda->p() + m_p->p() + m_pi->p()  ;
  m_pointing = ( sump * sin( m_BWerner ) ) / ( sump * sin( m_BWerner ) + sumpt );

  if ( m_log_BVtxchi2 > 10. ) //|| m_LWerner > 0. )
  {
    debug() << "TMVA Pesel" << endmsg;
    return StatusCode::FAILURE;
  }  

  return StatusCode::SUCCESS;
} 
//=============================================================================
StatusCode MySelectionAnalysis::BoostLWerner(LHCb::Particle B, LHCb::Particle Lambda, TVector3 rLambda, double &angle){
 
  TLorentzVector pB( B.momentum().Px(), B.momentum().Py(), B.momentum().Pz(), B.momentum().E() );
  TLorentzVector pLambda( Lambda.momentum().Px(), Lambda.momentum().Py(), Lambda.momentum().Pz(), Lambda.momentum().E() );  
  TLorentzVector xLambda( rLambda.x(), rLambda.y(), rLambda.z(), sqrt(rLambda.Mag2()) / ( Lambda.p() / Lambda.momentum().E()) );
  
  pLambda.Boost( -pB.BoostVector() );
  xLambda.Boost( -pB.BoostVector() ); 

  angle = pLambda.Vect().Angle( xLambda.Vect());

  if ( isnan(angle) )
    return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
  
}
//==============================================================================
StatusCode MySelectionAnalysis::ThetaAngleInB(LHCb::Particle B, LHCb::Particle Part, double &angle)
{
  TLorentzVector pB( B.momentum().Px(), B.momentum().Py(), B.momentum().Pz(), B.momentum().E() );
  TLorentzVector pPart( Part.momentum().Px(), Part.momentum().Py(), Part.momentum().Pz(), Part.momentum().E() );
  pPart.Boost( -pB.BoostVector() );

  angle = (pPart.Vect().Angle( pB.Vect() ));
  if ( isnan(angle) )
    return StatusCode::FAILURE;
  
  return StatusCode::SUCCESS;
} 
//=============================================================================
// Calculate Common Variables customised for the LL-Selection
//=============================================================================
StatusCode MySelectionAnalysis::LLCleanUp()
{
  debug() << "Begin LL cleanup" << endmsg;
  StatusCode sc = MySelectionAnalysis::Calculate();
  if ( sc = StatusCode::FAILURE)
    return sc;
  
  if ( m_ppiDist > 0.15)
    return StatusCode::FAILURE;
  
  if ( m_BWerner > 0.04 )
    return StatusCode::FAILURE;

  if ( m_log_BVtxchi2 > 5.)
    return StatusCode::FAILURE;

  if ( m_log_B_FLBchi2 < 0.75 )
    return StatusCode::FAILURE;

  if ( m_log_Lambda_FDchi2 < 0.)
    return StatusCode::FAILURE;

  if ( m_log_L_ctau < -1. || m_log_L_ctau > 2.5 )
    return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}
//=============================================================================
// Calculate Common Variables customised for the DD-Selection
//=============================================================================
StatusCode MySelectionAnalysis::DDCleanUp()
{
  debug() << "Begin DD cleanup" << endmsg;
  StatusCode sc = MySelectionAnalysis::Calculate();
  if ( sc = StatusCode::FAILURE)
    return sc;
  
  if ( m_ppiDist > 0.13 )
    return StatusCode::FAILURE;
  
  if ( m_BWerner > 0.03 )
    return StatusCode::FAILURE;

  if ( m_log_BVtxchi2 > 5.)
    return StatusCode::FAILURE;

  if ( m_log_Lambda_FDchi2 < 0.75)
    return StatusCode::FAILURE;

  if ( m_log_L_ctau < 0.8 || m_log_L_ctau > 2.5 )
    return StatusCode::FAILURE;

  if ( m_log_B_ctau > 0.5 || m_log_B_ctau < -1.25)
    return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Filter Events that triggered
//=============================================================================
StatusCode MySelectionAnalysis::FilterTrigger( const LHCb::Particle* IterB, 
                                               bool doFilter = false, 
                                               std::string DataPeriod = "2012b" )
{
  ITriggerTisTos::TisTosTob classifiedDec;
  bool bL0HadronTos(false), bL0GlobalTis(false);
  bool bHlt1TrackAllL0Tos(false);

  bool bHlt2Topo2BodyBBDTTos(false);
  bool bHlt2Topo3BodyBBDTTos(false);
  bool bHlt2Topo4BodyBBDTTos(false);

  bool bHlt2Topo2BodySimpleTos(false);
  bool bHlt2Topo3BodySimpleTos(false);
  bool bHlt2Topo4BodySimpleTos(false);

  _L0TriggerTool->setOfflineInput(*IterB);

  // Check for Hadron L0
  _L0TriggerTool->setTriggerInput("L0HadronDecision");
  classifiedDec = _L0TriggerTool->tisTosTobTrigger();
  bL0HadronTos = classifiedDec.tos();
  // Check for Global TIS
  _L0TriggerTool->setTriggerInput("L0Global");
  classifiedDec = _L0TriggerTool->tisTosTobTrigger();
  bL0GlobalTis = classifiedDec.tis();

  // Discover if split Hlt1,Hlt2 or old style
  bool split(false);
  // split = exist<LHCb::HltDecReports>("Hlt1/DecReports",false);
  // if( !split ) split = exist<LHCb::HltDecReports>("Hlt1/DecReports");
  
  // Write Hlt1 information
  ITriggerTisTos * triggerTisTosTool = _HltTriggerTool;
  _HltTriggerTool->setOfflineInput(*IterB);
  if( split ){
    triggerTisTosTool = _Hlt1TriggerTool;
    triggerTisTosTool->setOfflineInput(*IterB);
  }
  triggerTisTosTool->setTriggerInput("Hlt1TrackAllL0Decision");
  classifiedDec = triggerTisTosTool->tisTosTobTrigger();
  bHlt1TrackAllL0Tos = classifiedDec.tos();

  // Write Hlt2 information
  triggerTisTosTool = _HltTriggerTool;
  _HltTriggerTool->setOfflineInput(*IterB);
  if( split ){
    triggerTisTosTool = _Hlt2TriggerTool;
    triggerTisTosTool->setOfflineInput(*IterB);
  }
  triggerTisTosTool->setTriggerInput("Hlt2Topo2BodyBBDTDecision");
  classifiedDec = triggerTisTosTool->tisTosTobTrigger();
  bHlt2Topo2BodyBBDTTos = classifiedDec.tos();
  triggerTisTosTool->setTriggerInput("Hlt2Topo3BodyBBDTDecision");
  classifiedDec = triggerTisTosTool->tisTosTobTrigger();
  bHlt2Topo3BodyBBDTTos = classifiedDec.tos();
  classifiedDec = triggerTisTosTool->tisTosTobTrigger();
  triggerTisTosTool->setTriggerInput("Hlt2Topo4BodyBBDTDecision");
  classifiedDec = triggerTisTosTool->tisTosTobTrigger();
  bHlt2Topo4BodyBBDTTos = classifiedDec.tos();

  if (DataPeriod == "2011")
  {
    triggerTisTosTool->setTriggerInput("Hlt2Topo2BodySimpleDecision");
    classifiedDec = triggerTisTosTool->tisTosTobTrigger();
    bHlt2Topo2BodySimpleTos = classifiedDec.tos();
    triggerTisTosTool->setTriggerInput("Hlt2Topo3BodySimpleDecision");
    classifiedDec = triggerTisTosTool->tisTosTobTrigger();
    bHlt2Topo3BodySimpleTos = classifiedDec.tos();
    triggerTisTosTool->setTriggerInput("Hlt2Topo4BodySimpleDecision");
    classifiedDec = triggerTisTosTool->tisTosTobTrigger();
    bHlt2Topo4BodySimpleTos = classifiedDec.tos();
  }
  
  // m_L0HadronTos = bL0HadronTos; 
  // m_L0GlobalTis = bL0GlobalTis ;
  // m_Hlt1TrackAllL0Tos = bHlt1TrackAllL0Tos;
  // m_Hlt2TopoBBDTTos = bHlt2TopoBBDTTos; 

  if ( !doFilter ) return StatusCode::SUCCESS;

  bool triggered(false);

  if (DataPeriod != "2011")
    if(  bL0HadronTos == 1 || bL0GlobalTis == 1 ) 
      if ( bHlt1TrackAllL0Tos == 1 )
        if ( bHlt2Topo2BodyBBDTTos  == 1 || 
             bHlt2Topo3BodyBBDTTos  == 1 || 
             bHlt2Topo4BodyBBDTTos  == 1 )
          triggered = true;

  if (DataPeriod == "2011")
    if(  bL0HadronTos == 1 || bL0GlobalTis == 1 ) 
      if ( bHlt1TrackAllL0Tos == 1 )
        if ( bHlt2Topo2BodyBBDTTos  == 1 || 
             bHlt2Topo3BodyBBDTTos  == 1 || 
             bHlt2Topo4BodyBBDTTos  == 1 || 
             bHlt2Topo2BodySimpleTos == 1 ||
             bHlt2Topo3BodySimpleTos == 1 ||
             bHlt2Topo4BodySimpleTos == 1 )
          triggered = true;
  
  if( triggered )
    return StatusCode::SUCCESS;
  else
  {
    debug() << "Trigger-Filter Failure" << endmsg;
    return StatusCode::FAILURE;  
  }
}
