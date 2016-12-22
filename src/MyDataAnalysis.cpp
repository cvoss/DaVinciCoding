// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MyDataAnalysis.h"

#include "ReadMLPLL2012Uni.h"
#include "ReadMLPDD2012.h"
#include "ReadMLPLL2011Uni.h"
#include "ReadMLPDD2011Uni.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MyDataAnalysis
//
// 2014-05-19 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MyDataAnalysis )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MyDataAnalysis::MyDataAnalysis( const std::string& name,
                                ISvcLocator* pSvcLocator)
  : MySelectionAnalysis ( name , pSvcLocator )
  , LMass("LMass","LMass",1100,1130)
  , BMass("BMass","BMass",4900,5600)
  , BMassDTF("BMassDTF","BMassDTF",4900,5600)
  , BMassDTF_P("BMassDTF_P","BMassDTF_P", 0,500e3)
  , BMassDTF_PT("BMassDTF_PT","BMassDTF_PT", 0,500e3)
  , BMassDTF_Eta("BMassDTF_Eta","BMassDTF_Eta", 0,500e3)
  , p_pNNPID("p_pNNPID","p_pNNPID",-0.5,1.5)
  , p_KNNPID("p_KNNPID","p_KNNPID",-0.5,1.5)
  , p_piNNPID("p_piNNPID","p_piNNPID",-0.5,1.5)
  , pi_pNNPID("pi_pNNPID","pi_pNNPID",-0.5,1.5)
  , pi_KNNPID("pi_KNNPID","pi_KNNPID",-0.5,1.5)
  , pi_piNNPID("pi_piNNPID","pi_piNNPID",-0.5,1.5)
  , Lp_pNNPID("Lp_pNNPID","Lp_pNNPID",-0.5,1.5)
  , Lp_pDLLPID("Lp_pDLLPID","Lp_pDLLPID",-1e99,1e99)
  , Lp_KNNPID("Lp_KNNPID","Lp_KNNPID",-0.5,1.5)
  , Lp_piNNPID("Lp_piNNPID","Lp_piNNPID",-0.5,1.5)
  , Lpi_pNNPID("Lpi_pNNPID","Lpi_pNNPID",-0.5,1.5)
  , Lpi_KNNPID("Lpi_KNNPID","Lpi_KNNPID",-0.5,1.5)
  , Lpi_piNNPID("Lpi_piNNPID","Lpi_piNNPID",-0.5,1.5)
  , mLP("mLP","mLP",2053,5140)
  , mPpi("mPpi","mPpi",1053,5140)
  , mLpi("mLpi","mLpi",1053,5140)
  , mPLp("mPLp","mPLp",1870.,5010.)
  , mpiLpi("mpiLpi","mpiLpi",275,3410.)
  , Lp_TrkType("Lp_TrkType","Lp_TrkType",0,10)
  , Lpi_TrkType("Lpi_TrkType","Lpi_TrkType",0,10)
  , MVALL("MVALL","MVALL",-3,2)
  , MVADD("MVADD","MVADD",-3,2)
  , Polarity("Polarity","Polarity", -2,2)
  , log_pi_IP("log_pi_IP", "log_pi_IP", -10, 10)
  , log_p_IP("log_p_IP","log_p_IP", -10, 10)
  , log_Lp_IP("log_Lp_IP","log_Lp_IP", -10, 10)
  , log_Lpi_IP("log_Lpi_IP","log_Lpi_IP", -10, 10)
  , log_LIP("log_LIP","log_LIP",-10,10)
  , log_BIP("log_BIP","log_BIP",-10,10)
  , log_BVtxchi2("log_BVtxchi2","log_BVtxchi2", -10, 10)
  , FitProb("FitProb", "FitProb", -1, 2)
  , log_B_DauSumchi2("log_B_DauSumchi2", "log_B_DauSumchi2",-10, 10)
  , log_B_FLBchi2("log_B_FLBchi2", "log_B_FLBchi2", 0, 10)
  , log_Lambda_FDchi2("log_Lambda_FDchi2","log_LambdaFDchi2", -5, 10)
  , B_PT("B_PT","B_PT", 0,500e3)
  , B_Eta("B_Eta", "B_Eta", 2,6)
  , log_B_ctau("log_B_ctau","log_B_ctau",-10,10)
  , log_L_ctau("log_L_ctau","log_L_ctau",-10,10)
  , WernerB("BWerner","BWerner",0,1)
  , WernerL("LWerner","LWerner",-10,5)
  , ppi_Dist("ppiDist","ppiDist",-0.1,0.5)
  , LAngle("LAngle", "LAngle",-1,3.2)
  , Pointing("Pointing","Pointing",-1.,2.)
  , LP_PT("LP_PT","LP_PT",0.,15000.)
  , L_PT("L_PT","L_PT", 0,500e3)
  , p_PT("p_PT","p_PT", 0,500e3)
  , pi_PT("pi_PT","pi_PT", 0,500e3)
  , Lp_PT("Lp_PT","Lp_PT", 0,500e3)
  , Lpi_PT("Lpi_PT","Lpi_PT", 0,500e3)
  , DLcosTheta("DLcosTheta", "DLcosTheta",-1.1,1.1)
  , DLPhiL("DLPhiL","DLPhiL",-3.2,3.2)
  , DLPhihh("DLPhihh","DLPhihh",-3.2,3.2)
  , DpcosTheta("DpcosTheta", "DpcosTheta",-1.1,1.1)
  , DpPhiL("DpPhiL","DpPhiL",-3.2,3.2)
  , DpPhihh("DpPhihh","DpPhihh",-3.2,3.2)
  , Lambda_Phi("Lambda_Phi","Lambda_Phi",-3.2,3.2)
  , BMassMC("BMassMC","BMassMC",5000,5600)
  , mLPMC("mLPMC","mLPMC",2053,5140)
  , mPpiMC("mPpiMC","mPpiMC",1053,5140)
  , mLpiMC("mLpiMC","mLpiMC",1053,5140)
  , mPPMC("mPPMC","mPPMC",1800,5300)
  , mpipiMC("mpipiMC","mpipiMC",100,5300)
  , AP_pt("AP_pt","AP_pt", 0, 3000)
  , AP_alpha("AP_alpha","AP_alpha", -1.1,1.1)
  , DTF2Status("DTF2Status","DTF2Status",-1000,1000)
  , P_P("P_P","P_P",0,500e3)
  , P_TRACK_Eta("P_TRACK_Eta","P_TRACK_Eta",-1,6)
  , P_nTracks("P_nTracks","P_nTracks",-1,500)
  , Pi_P("Pi_P","Pi_P",0,500e3)
  , Pi_TRACK_Eta("Pi_TRACK_Eta","Pi_TRACK_Eta",-1,6)
  , Pi_nTracks("Pi_nTracks","Pi_nTracks",-1,500)
  , Lp_P("Lp_P","Lp_P",0,500e3)
  , Lp_TRACK_Eta("Lp_TRACK_Eta","Lp_TRACK_Eta",-1,6)
  , Lp_nTracks("Lp_nTracks","Lp_nTracks",-1,500)
  , nTracks("nTracks","nTracks",-1,500)
  , RunNumber("RunNumber","RunNumber",0,1000000.)
  , EventNumber("EventNumber","EventNumber",0,10e99)
  , OdinTCK("OdinTCK","OdinTCK",0,10e99)
  , DataPeriod("DataPeriod","DataPeriod", 0 ,3)
  , pi_KDLLPID("pi_KDLLPID","pi_KDLLPID",-1e99,1e99)
  , L0HadronTos("L0HadronTos","L0HadronTos",-1,2)
  , L0GlobalTis("L0GlobalTis","L0GlobalTis",-1,2)
  , Hlt1TrackAllL0Tos("Hlt1TrackAllL0Tos","Hlt1TrackAllL0Tos",-1,2)
  , Hlt2Topo2BodyBBDTTos("Hlt2Topo2BodyBBDTTos","Hlt2Topo2BodyBBDTTos",-1,2)
  , Hlt2Topo3BodyBBDTTos("Hlt2Topo3BodyBBDTTos","Hlt2Topo3BodyBBDTTos",-1,2)
  , Hlt2Topo4BodyBBDTTos("Hlt2Topo4BodyBBDTTos","Hlt2Topo4BodyBBDTTos",-1,2)
  , Hlt2Topo2BodySimpleTos("Hlt2Topo2BodySimpleTos","Hlt2Topo2BodySimpleTos",-1,2)
  , Hlt2Topo3BodySimpleTos("Hlt2Topo3BodySimpleTos","Hlt2Topo3BodySimpleTos",-1,2)
  , Hlt2Topo4BodySimpleTos("Hlt2Topo4BodySimpleTos","Hlt2Topo4BodySimpleTos",-1,2)
  , SanityB("SanityB","SanityB",-1e99,1e99)
  , SanityL("SanityL","SanityL",-1e99,1e99)
  , SpinTPCUnit("SpinTPCUnit","SpinTPCUnit",-1e99,1e99)
  , SpinTPC("SpinTPC","SpinTPC",-1e99,1e99)
  , SpinTPCMC("SpinTPCMC","SpinTPCMC",-1e99,1e99)
  , BTag("BTag","BTag",-2.,2.)
  , pHelicity("pHelicity","pHelicity",-0.5,3.5)
  , cospHelicity("cospHelicity","cospHelicity",-1.1,1.1)
  , cosLpHelicity_Enh("cosLpHelicity_Enh","cosLpHelicity_Enh",-1.1,1.1)
  , cosLpHelicity_3b("cosLpHelicity_3b","cosLpHelicity_3b",-1.1,1.1)
  , E_L_Brest("E_L_Brest","E_L_Brest",0.,5000.)
  , cosPoleBaryHelicity("cosPoleBaryHelicity","cosPoleBaryHelicity",-1.1,1.1)
  , cosPBPoleAngle("cosDauPPiAngle","cosDauPpiAngle",-1.1,1.1)
  , nCand("nCand","nCand",-1,1e99)
{
  declareProperty("BMassLowEdge", m_BMassLowEdge = 4900.*Gaudi::Units::MeV);
  declareProperty("BMassHighEdge", m_BMassHighEdge = 5650.*Gaudi::Units::MeV);
  declareProperty("LambdaMassLowEdge", m_LMassLowEdge = 1105.*Gaudi::Units::MeV);
  declareProperty("LambdaMassHighEdge", m_LMassHighEdge = 1125.*Gaudi::Units::MeV);
  declareProperty("RootFile", RootFile = "sPlot.root");
  declareProperty("DoSplot", m_sPlot = false);
  declareProperty("PIDToolnames", m_PIDToolnames, "List of PID tools");
  declareProperty("TriggerFiltering", m_TriggerFilter = true);
  declareProperty("DataPeriod", m_dataPeriod = "2012b" );
  declareProperty("PerformSanityChecks", m_checks = false );
  declareProperty("MCTuples", m_mctuple = false);
  
}
//=============================================================================
// Destructor
//=============================================================================
MyDataAnalysis::~MyDataAnalysis() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MyDataAnalysis::initialize() {
  StatusCode sc = MySelectionAnalysis::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  warning() << "Using the following configuration \n"
            << "BMassLowEdge " << m_BMassLowEdge << "\n"
            << "BMassHighEdge" <<  m_BMassHighEdge << "\n"
            << "LambdaMassLowEdge" <<  m_LMassLowEdge << "\n"
            << "LambdaMassHighEdge" << m_LMassHighEdge << "\n"
            << "RootFile" << RootFile << "\n"
            << "DoSplot" << m_sPlot << "\n"
            << "TrigerFiltering" << m_TriggerFilter << "\n"
            << "DataPeriod " << m_dataPeriod << "\n"
            << endmsg;
  for (unsigned int i = 0; i < m_PIDToolnames.size(); i++)
    warning() << i <<"th tool" << m_PIDToolnames.at(i) << endmsg;

  if ( (m_dataPeriod != "2011") && (m_dataPeriod != "2012a") && (m_dataPeriod != "2012b") )
    fatal() << "Data Period not correctly configured " << m_dataPeriod << endmsg;
  
  m_magfieldsvc = svc<ILHCbMagnetSvc> ( "MagneticFieldSvc", true );
  
  // _L0TriggerTool = tool<ITriggerTisTos>("L0TriggerTisTos",this);
  // _HltTriggerTool = tool<ITriggerTisTos>("TriggerTisTos",this );
  // _Hlt1TriggerTool = tool<ITriggerTisTos>("Hlt1TriggerTisTos",this);
  // _Hlt2TriggerTool = tool<ITriggerTisTos>("Hlt2TriggerTisTos",this);
  
  for (unsigned int i = 0; i < m_PIDToolnames.size(); i++)
    m_PIDTools.push_back( tool<IParticleManipulator>( m_PIDToolnames.at(i),this) );

  AngularDalitz            = tool<IPhysicsComputation>("AngularDalitz",this);
  ArmenterosPodolanski     = tool<IPhysicsComputation>("AmenterosPodolanski",this);
  HelicityAngle            = tool<IPhysicsComputation>("HelicityAngle",this);
  PolarityAngleEnhancement = tool<IPhysicsComputation>("HelicityAngle",this);
  PolarityAngleBMeson      = tool<IPhysicsComputation>("HelicityAngle",this);
  AnglePoleBMeson          = tool<IPhysicsComputation>("HelicityAngle",this);
  AngleBP                  = tool<IPhysicsComputation>("HelicityAngle",this);
  
  _assoc               = tool<IParticle2MCAssociator>("MCMatchObjP2MCRelator",this);

  debug() << "Get properties" << endmsg;
    
  const LHCb::ParticleProperty* LambdaInfo = m_ppsvc->find( "Lambda0" ) ;
  const LHCb::ParticleProperty* BInfo      = m_ppsvc->find( "B0" ) ;
  const LHCb::ParticleProperty* pInfo      = m_ppsvc->find( "p+" );
  const LHCb::ParticleProperty* piInfo     = m_ppsvc->find( "pi+" );

  m_LMass = LambdaInfo->mass();
  m_BMass = BInfo->mass();
  m_pMass = pInfo->mass();
  m_piMass = piInfo->mass();

  BMassDTF.setMin( m_BMassLowEdge);
  BMassDTF.setMax( m_BMassHighEdge);
  BMass.setMin( m_BMassLowEdge);
  BMass.setMax( m_BMassHighEdge);

  LMass.setMin( m_LMassLowEdge );
  LMass.setMax( m_LMassHighEdge );
  
  mLP.setMin( m_LMass + m_pMass );
  mLP.setMax( m_BMass - m_piMass );
  mPpi.setMin( m_pMass + m_piMass );
  mPpi.setMax( m_BMass - m_LMass );
  mLpi.setMin( m_LMass + m_piMass );
  mLpi.setMax( m_BMass - m_pMass );

  mLPMC.setMin( m_LMass + m_pMass );
  mLPMC.setMax( m_BMass - m_piMass );
  mPpiMC.setMin( m_pMass + m_piMass );
  mPpiMC.setMax( m_BMass - m_LMass );
  mLpiMC.setMin( m_LMass + m_piMass );
  mLpiMC.setMax( m_BMass - m_pMass );
  mPPMC.setMin( 2*m_pMass );
  mPPMC.setMax( m_BMass - 2*m_piMass );
  mpipiMC.setMin( 2*m_piMass );
  mpipiMC.setMax( m_BMass - 2*m_pMass );

  debug() << "Define Output ArgSets" << endmsg;
    
  Set = new RooArgSet(LMass,BMassDTF,p_pNNPID,pi_KNNPID,Lp_TrkType,Lpi_TrkType,MVALL,MVADD,"ArgSet");
  Set->add(RooArgSet(p_KNNPID,p_piNNPID,pi_pNNPID,pi_piNNPID));
  Set->add(RooArgSet(Lp_pNNPID,Lp_pDLLPID, Lp_KNNPID,Lp_piNNPID,Lpi_pNNPID,Lpi_KNNPID,Lpi_piNNPID));
  Set->add(RooArgSet(mLP,mPpi,mLpi,mPLp,mpiLpi,Polarity));
  
  Set->add(RooArgList(log_p_IP, log_pi_IP, log_Lp_IP, log_Lpi_IP, log_BIP, log_LIP, log_BVtxchi2, FitProb));
  Set->add(RooArgSet(log_B_DauSumchi2, log_B_FLBchi2, log_Lambda_FDchi2, B_PT, B_Eta, log_B_ctau, log_L_ctau));
  Set->add(RooArgSet(WernerB, WernerL, ppi_Dist,LAngle,LP_PT,Pointing,BMass));
  Set->add(RooArgSet(L_PT,p_PT,pi_PT,Lp_PT,Lpi_PT,DLcosTheta, DLPhiL, DLPhihh));
  Set->add(RooArgSet(DpcosTheta, DpPhiL, DpPhihh, Lambda_Phi));
  Set->add(RooArgSet(mLPMC,mPpiMC,mLpiMC,mPPMC,mpipiMC,BMassMC,AP_pt,AP_alpha,DTF2Status));
  Set->add(RooArgSet(P_P, P_TRACK_Eta, P_nTracks,Pi_P, Pi_TRACK_Eta, Pi_nTracks));
  Set->add(RooArgSet(Lp_P, Lp_TRACK_Eta, Lp_nTracks,nTracks,BMassDTF_P,BMassDTF_PT,BMassDTF_Eta));
  Set->add(RooArgSet(pi_KDLLPID));
  Set->add(RooArgSet(RunNumber, EventNumber, OdinTCK));
  Set->add(RooArgSet(DataPeriod));
  Set->add(RooArgSet(L0HadronTos,L0GlobalTis,Hlt1TrackAllL0Tos));
  Set->add(RooArgSet(Hlt2Topo2BodyBBDTTos,Hlt2Topo3BodyBBDTTos));
  Set->add(RooArgSet(Hlt2Topo4BodyBBDTTos,Hlt2Topo2BodySimpleTos));
  Set->add(RooArgSet(Hlt2Topo3BodySimpleTos,Hlt2Topo4BodySimpleTos));
  Set->add(RooArgSet(SanityB,SanityL,SpinTPCUnit,SpinTPC,SpinTPCMC,BTag));
  Set->add(RooArgSet(pHelicity, cospHelicity));
  Set->add(RooArgSet(cosLpHelicity_Enh, cosLpHelicity_3b, E_L_Brest,cosPoleBaryHelicity,cosPBPoleAngle));
  Set->add(RooArgSet(nCand));
  
  debug() << "Initialise Output RooDataSet" << endmsg;

  rawdata = new RooDataSet("rawdata","rawdata", *Set);
  
  debug() << "Setting Data Periods" << endmsg;

  if ( m_dataPeriod == "2011" ) DataPeriod.setVal(1.);
  if ( m_dataPeriod == "2012a" ) DataPeriod.setVal(2.);
  if ( m_dataPeriod == "2012b" ) DataPeriod.setVal(3.);

  debug() << "Set Data Periods" << endmsg;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MyDataAnalysis::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.
  
  LHCb::Particle::Range Bs = this->particles();
  
  const LHCb::RecSummary * tracks = getIfExists<LHCb::RecSummary>("/Event/Rec/Summary");
  if( tracks )
    m_nTracks = tracks->info(LHCb::RecSummary::nTracks,-1);
  else m_nTracks = -1;
  
  //LHCb::MCParticles* particles = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default ); 
  
  debug() << "Size of Inputs " << Bs.size() << endmsg;
  if ( 0 == Bs.size() )
  {
    warning() << "Empty particle container" << endmsg;
    return StatusCode::SUCCESS;
  }
  nCand.setVal(Bs.size());
  // debug() << "Initialise the Particle iterator" << endmsg;
  LHCb::Particle::ConstVector::const_iterator IterB;

  // if ( 1 == Bs.size() )
  // {
  //   debug() << "Initialise the Particle iterator for n=1" << endmsg;
  //   IterB = Bs.begin();
  //   debug() << "Initialise done" << endmsg;
  //       }
  
  // else {
  //   debug() << "Initialise the Particle iterator for n!=1" << endmsg;
  //   IterB = select_randomly(Bs.begin(), Bs.end());
  //   debug() << "Initialise done" << endmsg;
  // }
  
  // if (!(*IterB))
  // {
  //   warning() << "No Candidate chosen" << endmsg;
  //   return StatusCode::SUCCESS;
  // }
  
  std::vector<IParticleManipulator*>::iterator iTool;
  std::vector<IParticleManipulator*>::reverse_iterator riTool;

  debug() << "Begin Event Selection" << endmsg;
  
  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    debug() << "###################################################################" << endmsg;
    debug() << "Taking a " << (*IterB)->particleID().pid() << " and do substitution" << endmsg;
    for ( iTool = m_PIDTools.begin(); iTool != m_PIDTools.end(); ++iTool)
      (*iTool) -> doCorrection(const_cast<LHCb::Particle*>(*IterB));

    if ( FilterTrigger( *IterB, m_TriggerFilter, m_dataPeriod) == StatusCode::SUCCESS )
      if (PreSelection( *IterB, m_BMassLowEdge, m_BMassHighEdge, m_LMassLowEdge, m_LMassHighEdge) == StatusCode::SUCCESS)
        if ( PerformFit( *IterB ) == StatusCode::SUCCESS )
        {
          bool LLCase(false);
          if ( (int)m_Lp->proto()->track()->type() == 3 || (int)m_Lpi->proto()->track()->type() == 3 )
            LLCase = true;
          debug() << "LL or DD Event" << endmsg;
          if ( ( LLCase ? LLCleanUp() : DDCleanUp() ) == StatusCode::SUCCESS )
            if ( m_p->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5 &&
                 m_pi->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5)
            {
              if ( m_mctuple )
                if( truthMatch(*(*IterB)) ) 
                {
                  MyDataAnalysis::FillDataSet();
                  MyDataAnalysis::TriggerWriteOut( *IterB );
                  DTF2Status.setVal(MyDataAnalysis::constrainBMass(*IterB));
                  MyDataAnalysis::EventInfo();
                  MyDataAnalysis::CalculateMCTPC(*IterB);
                  
                  counter("Found_Bs")++;
                  rawdata->add(*Set,1.0);
                }
              if ( !m_mctuple )
                {
                  MyDataAnalysis::FillDataSet();
                  MyDataAnalysis::TriggerWriteOut( *IterB );
                  DTF2Status.setVal(MyDataAnalysis::constrainBMass(*IterB));
                  MyDataAnalysis::EventInfo();
                                   
                  counter("Found_Bs")++;
                  rawdata->add(*Set,1.0);
                }
            }
        }
    
    for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
      (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MyDataAnalysis::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;
  MyDataAnalysis::WriteDataSet();
  return MySelectionAnalysis::finalize();  // must be called after all other actions
}
//============================================================================
// Fill the Data Set
//============================================================================
StatusCode MyDataAnalysis::FillDataSet()
{
  debug() << "Filling Data sets" << endmsg;
  
  std::vector<double> m_MVAinputValues;
  m_MVAinputValues.push_back( m_log_pi_IP );
  m_MVAinputValues.push_back( m_log_p_IP );
  m_MVAinputValues.push_back( m_log_Lpi_IP );
  m_MVAinputValues.push_back( m_log_Lp_IP );
  m_MVAinputValues.push_back( m_log_BIP );
  m_MVAinputValues.push_back( m_log_LIP );
  m_MVAinputValues.push_back( m_log_BVtxchi2 );
  m_MVAinputValues.push_back( m_FitProb );
  m_MVAinputValues.push_back( m_log_B_DauSumchi2 );
  m_MVAinputValues.push_back( m_log_B_FLBchi2 );
  m_MVAinputValues.push_back( m_log_Lambda_FDchi2 );
  m_MVAinputValues.push_back( m_B_PT );
  m_MVAinputValues.push_back( m_B_Eta );
  m_MVAinputValues.push_back( m_log_B_ctau );
  m_MVAinputValues.push_back( m_log_L_ctau );
  m_MVAinputValues.push_back( m_BWerner );
  m_MVAinputValues.push_back( m_LWerner );
  m_MVAinputValues.push_back( m_ppiDist );
  m_MVAinputValues.push_back( m_LAngle );  
  m_MVAinputValues.push_back( m_pointing );  

  log_pi_IP.setVal(m_MVAinputValues.at(0));
  log_p_IP.setVal(m_MVAinputValues.at(1));
  log_Lpi_IP.setVal(m_MVAinputValues.at(2));
  log_Lp_IP.setVal(m_MVAinputValues.at(3));
  log_BIP.setVal(m_MVAinputValues.at(4));
  log_LIP.setVal(m_MVAinputValues.at(5));
  log_BVtxchi2.setVal(m_MVAinputValues.at(6));
  FitProb.setVal(m_MVAinputValues.at(7));
  log_B_DauSumchi2.setVal(m_MVAinputValues.at(8));
  log_B_FLBchi2.setVal(m_MVAinputValues.at(9));
  log_Lambda_FDchi2.setVal(m_MVAinputValues.at(10));
  B_PT.setVal(m_MVAinputValues.at(11));
  B_Eta.setVal(m_MVAinputValues.at(12));
  log_B_ctau.setVal(m_MVAinputValues.at(13));
  log_L_ctau.setVal(m_MVAinputValues.at(14));
  WernerB.setVal(m_MVAinputValues.at(15));
  WernerL.setVal(m_MVAinputValues.at(16));
  ppi_Dist.setVal(m_MVAinputValues.at(17));
  LAngle.setVal(m_MVAinputValues.at(18));
  Pointing.setVal(m_MVAinputValues.at(19));
  
  const char* MVAinputVars[] = 
    {
      "log_pi_IP",
      "log_p_IP",
      "log_Lpi_IP",
      "log_Lp_IP", 
      "log_BIP", 
      "log_LIP", 
      "log_BVtxchi2", 
      "FitProb",
      "log_B_DauSumchi2",
      "log_B_FLBchi2",
      "log_Lambda_FDchi2",
      "B_PT",
      "B_Eta", 
      "log_B_ctau", 
      "log_L_ctau", 
      "BWerner", 
      "LWerner",
      "ppiDist", 
      "LAngle",
      "Pointing"
    };

  std::vector<std::string> theMVAInputVars;
  for ( unsigned int j = 0; j < m_MVAinputValues.size() ; ++j ) theMVAInputVars.push_back(MVAinputVars[j]);      

  ReadMLPLL2011Uni MLPSelLL2011     (theMVAInputVars);
  ReadMLPDD2011Uni MLPSelDD2011     (theMVAInputVars);

  ReadMLPLL2012Uni MLPSelLL2012     (theMVAInputVars);
  ReadMLPDD2012    MLPSelDD2012     (theMVAInputVars);

  double MVALLVal(-2), MVADDVal(-2);

  if ( (int)m_Lp->proto()->track()->type() == 3 || (int)m_Lpi->proto()->track()->type() == 3 )
  {
    if ( m_dataPeriod == "2012b" || m_dataPeriod == "2012a" )
    {
      MVALLVal = MLPSelLL2012.GetMvaValue(m_MVAinputValues);
    }
    if ( m_dataPeriod == "2011" )
      MVALLVal = MLPSelLL2011.GetMvaValue(m_MVAinputValues);
  }
  
  if ( (int)m_Lp->proto()->track()->type() == 5 || (int)m_Lpi->proto()->track()->type() == 5 )
  {
    if ( m_dataPeriod == "2012b" || m_dataPeriod == "2012a" )
    {
      MVADDVal = MLPSelDD2012.GetMvaValue(m_MVAinputValues); 
    }
    if ( m_dataPeriod == "2011" )
      MVADDVal = MLPSelDD2011.GetMvaValue(m_MVAinputValues); 
  }
 
  MVALL.setVal(MVALLVal);
  MVADD.setVal(MVADDVal);
  //Pointing.setVal(m_pointing);
  
  BMassDTF.setVal( m_B->momentum().M() );
  BMassDTF_P.setVal( m_B->momentum().P() );
  BMassDTF_PT.setVal( m_B->momentum().Pt() );
  BMassDTF_Eta.setVal( m_B->momentum().Eta() );

  BMass.setVal( m_B->measuredMass() );
  LMass.setVal(m_LMassVal);
  
  mLP.setVal( (m_Lambda->momentum() + m_p->momentum()).M() );
  mPpi.setVal( (m_p->momentum() + m_pi->momentum()).M() );
  mLpi.setVal( (m_pi->momentum() + m_Lambda->momentum()).M() );
  mPLp.setVal( (m_p->momentum() + m_Lp->momentum()).M() );
  mpiLpi.setVal( (m_pi->momentum() + m_Lpi->momentum()).M() );
  
  Lp_TrkType.setVal((double)m_Lp->proto()->track()->type());
  Lpi_TrkType.setVal((double)m_Lpi->proto()->track()->type());

  L_PT.setVal(m_L_PT);
  p_PT.setVal(m_p_PT);
  pi_PT.setVal(m_pi_PT);
  Lp_PT.setVal(m_Lp_PT);
  Lpi_PT.setVal(m_Lpi_PT);
  
  p_pNNPID.setVal   ( m_p_pNNPID );
  p_KNNPID.setVal   ( m_p_KNNPID );
  p_piNNPID.setVal  ( m_p_piNNPID );
  pi_pNNPID.setVal  ( m_pi_pNNPID );
  pi_KNNPID.setVal  ( m_pi_KNNPID );
  pi_piNNPID.setVal ( m_pi_piNNPID );
  Lp_pNNPID.setVal  ( m_Lp_pNNPID );
  Lp_pDLLPID.setVal ( m_Lp_pDLLPID );
  Lp_KNNPID.setVal  ( m_Lp_KNNPID );
  Lp_piNNPID.setVal ( m_Lp_piNNPID );
  Lpi_pNNPID.setVal ( m_Lpi_pNNPID );
  Lpi_KNNPID.setVal ( m_Lpi_KNNPID );
  Lpi_piNNPID.setVal( m_Lpi_piNNPID );
  pi_KDLLPID.setVal( m_pi_KDLLPID );

  Polarity.setVal( m_magfieldsvc->isDown() ? -1. : 1.);

  LHCb::Particle::ConstVector LCands;
  LCands.push_back(m_B);
  LCands.push_back(m_Lambda);
  LCands.push_back(m_p);
  LCands.push_back(m_pi);

  LHCb::Particle::ConstVector pCands;
  pCands.push_back(m_B);
  pCands.push_back(m_p);
  pCands.push_back(m_Lambda);
  pCands.push_back(m_pi);

  std::vector<double> LDAngleVars(3);
  std::vector<double> pDAngleVars(3);

  AngularDalitz->initPhysics( LCands );
  AngularDalitz->computePhysics( &LDAngleVars );

  AngularDalitz->initPhysics( pCands );
  AngularDalitz->computePhysics( &pDAngleVars );

  DLcosTheta.setVal(LDAngleVars.at(0));
  DLPhiL.setVal(LDAngleVars.at(1));
  DLPhihh.setVal(LDAngleVars.at(2));

  DpcosTheta.setVal(LDAngleVars.at(0));
  DpPhiL.setVal(LDAngleVars.at(1));
  DpPhihh.setVal(LDAngleVars.at(2));

  Lambda_Phi.setVal(m_Lambda->momentum().Phi());
  
  LP_PT.setVal(m_LP_PT);

  std::vector<double> AP_vars(2);
  
  LHCb::Particle::ConstVector APCands;
  APCands.push_back(m_Lambda);
  APCands.push_back(m_Lp);
  APCands.push_back(m_Lpi);

  ArmenterosPodolanski->initPhysics( APCands );
  ArmenterosPodolanski->computePhysics( &AP_vars );

  AP_pt.setVal( AP_vars.at(0) );
  AP_alpha.setVal( AP_vars.at(1) );

  P_P.setVal(m_p->p());
  Pi_P.setVal(m_pi->p());
  Lp_P.setVal(m_Lp->p());
  P_TRACK_Eta.setVal(m_p->proto()->track()->pseudoRapidity());
  Pi_TRACK_Eta.setVal(m_pi->proto()->track()->pseudoRapidity());
  Lp_TRACK_Eta.setVal(m_Lp->proto()->track()->pseudoRapidity());
  
  debug() << " and mDTF = " << BMassDTF.getVal() << " found" << endmsg;
  
  return StatusCode::SUCCESS;
}
 
//=============================================================================
// Write Trigger Output
//=============================================================================
StatusCode MyDataAnalysis::TriggerWriteOut( const LHCb::Particle *IterB ){

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
  _HltTriggerTool->setOfflineInput(*IterB);
 
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
  split = exist<LHCb::HltDecReports>("Hlt1/DecReports",false);
  if( !split ) split = exist<LHCb::HltDecReports>("Hlt1/DecReports");
  
  // Write Hlt1 information
  {
    ITriggerTisTos *triggerTisTosTool = _HltTriggerTool;
    if( split ){
      triggerTisTosTool = _Hlt1TriggerTool;
      triggerTisTosTool->setOfflineInput(_HltTriggerTool->offlineLHCbIDs() );
    }
    triggerTisTosTool->setTriggerInput("Hlt1TrackAllL0Decision");
    classifiedDec = triggerTisTosTool->tisTosTobTrigger();
    bHlt1TrackAllL0Tos = classifiedDec.tos();
  }
  
  // Write Hlt2 information
  {
    ITriggerTisTos *triggerTisTosTool = _HltTriggerTool;
    if( split ){
      triggerTisTosTool = _Hlt2TriggerTool;
      triggerTisTosTool->setOfflineInput(_HltTriggerTool->offlineLHCbIDs() );
    }
    triggerTisTosTool->setTriggerInput("Hlt2Topo2BodyBBDTDecision");
    classifiedDec = triggerTisTosTool->tisTosTobTrigger();
    bHlt2Topo2BodyBBDTTos = classifiedDec.tos();
    triggerTisTosTool->setTriggerInput("Hlt2Topo3BodyBBDTDecision");
    classifiedDec = triggerTisTosTool->tisTosTobTrigger();
    bHlt2Topo3BodyBBDTTos = classifiedDec.tos();
    triggerTisTosTool->setTriggerInput("Hlt2Topo4BodyBBDTDecision");
    classifiedDec = triggerTisTosTool->tisTosTobTrigger();
    bHlt2Topo4BodyBBDTTos = classifiedDec.tos();
  }
  
  if (m_dataPeriod == "2011")
  {
    ITriggerTisTos *triggerTisTosTool = _HltTriggerTool;
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
  
  L0HadronTos.setVal((double)bL0HadronTos);
  L0GlobalTis.setVal((double)bL0GlobalTis);
  Hlt1TrackAllL0Tos.setVal((double)bHlt1TrackAllL0Tos);
  Hlt2Topo2BodyBBDTTos.setVal((double)bHlt2Topo2BodyBBDTTos);
  Hlt2Topo3BodyBBDTTos.setVal((double)bHlt2Topo3BodyBBDTTos);
  Hlt2Topo4BodyBBDTTos.setVal((double)bHlt2Topo4BodyBBDTTos);
  Hlt2Topo2BodySimpleTos.setVal((double)bHlt2Topo2BodySimpleTos);
  Hlt2Topo3BodySimpleTos.setVal((double)bHlt2Topo3BodySimpleTos);
  Hlt2Topo4BodySimpleTos.setVal((double)bHlt2Topo4BodySimpleTos);
  
  return StatusCode::SUCCESS;
}
//=============================================================================
// Write DataSets to file
//=============================================================================
StatusCode MyDataAnalysis::WriteDataSet(){

  RooDataSet *LL_Lambda      = (RooDataSet*)rawdata        ->reduce("Lp_TrkType==3&&Lpi_TrkType==3");
  RooDataSet *DD_Lambda      = (RooDataSet*)rawdata        ->reduce("Lp_TrkType==5&&Lpi_TrkType==5");

  // RooDataSet *LL_Lambda_Up   = (RooDataSet*)rawdata        ->reduce("Lp_TrkType==3&&Lpi_TrkType==3&&Polarity>0.");
  // RooDataSet *LL_Lambda_Down = (RooDataSet*)rawdata        ->reduce("Lp_TrkType==3&&Lpi_TrkType==3&&Polarity<0.");

  // RooDataSet *DD_Lambda_Up   = (RooDataSet*)rawdata        ->reduce("Lp_TrkType==5&&Lpi_TrkType==5&&Polarity>0.");
  // RooDataSet *DD_Lambda_Down = (RooDataSet*)rawdata        ->reduce("Lp_TrkType==5&&Lpi_TrkType==5&&Polarity<0.");

  RooDataSet *LL_B           = (RooDataSet*)LL_Lambda      ->reduce("LMass>1113&&LMass<1118");
  RooDataSet *DD_B           = (RooDataSet*)DD_Lambda      ->reduce("LMass>1111&&LMass<1122");

  LL_B                       ->SetName("ReducedLL");
  DD_B                       ->SetName("ReducedDD");

  warning() << "Constructing Datasets" << endmsg;
  RooDataSet* data = new RooDataSet(*LL_B, "MergedData");
  data->append(*DD_B);
  data->SetName("ReducedData");
  warning() << "# of Entries " << data->numEntries() << endmsg;

  TFile Output(RootFile.c_str(),"recreate");
  Output.cd();
  rawdata->Write();
  data->Write();
  LL_B->Write();
  DD_B->Write();

  // LL_B_Up->Write();
  // DD_B_Up->Write();
  
  // LL_B_Down->Write();
  // DD_B_Down->Write();
  Output.Close();
  return StatusCode::SUCCESS;
}
//=============================================================================
// Constrain to BMass and refit
//=============================================================================
StatusCode MyDataAnalysis::constrainBMass(const LHCb::Particle *BCand)
{
  DecayTreeFitter::Fitter decayTreeFitter((*BCand), *m_BPV, m_states);
  decayTreeFitter.setMassConstraint ( m_LambdaID );
  decayTreeFitter.setMassConstraint ( m_BbarID );
  decayTreeFitter.updateCand(*const_cast<LHCb::Particle*>(BCand));
  decayTreeFitter.fit();
  
  LHCb::DecayTree fitted = decayTreeFitter.getFittedTree() ;
  const LHCb::Particle *B, *Lambda, *p, *pi, *Lp, *Lpi;
  
  debug() << "2nd DTF Status: " << decayTreeFitter.status() << endmsg;
  
  B = fitted.head();
  if( !B ) return StatusCode::FAILURE;

  Lambda = fitted.find( B->particleID() == m_BbarID ? m_LambdaID    : m_LambdabarID );
  p      = fitted.find( B->particleID() == m_BbarID ? m_protonbarID : m_protonID );
  pi     = fitted.find( B->particleID() == m_BbarID ? m_pipID       : m_pimID );
  Lpi    = fitted.find( m_B->particleID() == m_BbarID ? m_pimID       : m_pipID );
  Lp     = fitted.find( m_B->particleID() == m_BbarID ? m_protonID    : m_protonbarID );
  if(B->particleID().abspid() == m_BsbarID.abspid())
  {
    Lambda = fitted.find( m_B->particleID() == m_BsbarID ? m_LambdaID    : m_LambdabarID );
    p      = fitted.find( m_B->particleID() == m_BsbarID ? m_protonbarID : m_protonID );
    pi     = fitted.find( m_B->particleID() == m_BsbarID ? m_KpID       : m_KmID );
    Lpi    = fitted.find( m_B->particleID() == m_BsbarID ? m_pimID       : m_pipID );
    Lp     = fitted.find( m_B->particleID() == m_BsbarID ? m_protonID    : m_protonbarID );
  }
  debug() << B << " " << Lambda << " " << p << " " << pi << endmsg;

  BMassMC.setVal(B->momentum().M());
  mLPMC.setVal( (Lambda->momentum() + p->momentum()).M() );
  mPpiMC.setVal( (p->momentum() + pi->momentum()).M() );
  mLpiMC.setVal( (pi->momentum() + Lambda->momentum()).M() );
  mPPMC.setVal( (p->momentum() + Lp->momentum()).M() );
  mpipiMC.setVal( (pi->momentum() + Lpi->momentum()).M() );

  //===========================================================================
  // Calculation of the TPC and the helicity angle
  //===========================================================================
  TLorentzVector BMeson    ( (Double_t)(B->momentum().Px()),
                             (Double_t)(B->momentum().Py()),
                             (Double_t)(B->momentum().Pz()),
                             (Double_t)(B->momentum().E()) );
  
  TLorentzVector LambdaBary( (Double_t)(Lambda->momentum().Px()),
                             (Double_t)(Lambda->momentum().Py()),
                             (Double_t)(Lambda->momentum().Pz()),
                             (Double_t)(Lambda->momentum().E()) );

  TLorentzVector Proton    ( (Double_t)(p->momentum().Px()),
                             (Double_t)(p->momentum().Py()),
                             (Double_t)(p->momentum().Pz()),
                             (Double_t)(p->momentum().E()) );

  TLorentzVector Pion      ( (Double_t)(pi->momentum().Px()),
                             (Double_t)(pi->momentum().Py()),
                             (Double_t)(pi->momentum().Pz()),
                             (Double_t)(pi->momentum().E()) );
  
  TLorentzVector DauProton ( (Double_t)(Lp->momentum().Px()),
                             (Double_t)(Lp->momentum().Py()),
                             (Double_t)(Lp->momentum().Pz()),
                             (Double_t)(Lp->momentum().E()) );

  TLorentzVector Baryonium = Proton + LambdaBary;
  
  std::vector<TLorentzVector> cands;

  info()<< "Applying first Initial Helicity Tool"<< endmsg;
  cands.push_back( BMeson );
  cands.push_back( Baryonium );
  cands.push_back( Proton );
  
  HelicityAngle->initPhysics(cands);
  std::vector<double> HeliAngle(2);
  HelicityAngle->computePhysics( &HeliAngle);
  pHelicity.setVal(  HeliAngle.at(0) );
  cospHelicity.setVal(  HeliAngle.at(1) );

  info()<< "Applying first Helicity Tool"<< endmsg;
  std::vector<TLorentzVector> PolarityPolecands;
  PolarityPolecands.push_back( Baryonium );
  PolarityPolecands.push_back( LambdaBary );
  PolarityPolecands.push_back( DauProton );
  std::vector<double> PolarityAngle(2);
  PolarityAngleEnhancement->initPhysics(PolarityPolecands);
  PolarityAngleEnhancement->computePhysics( &PolarityAngle);
  cosLpHelicity_Enh.setVal(  PolarityAngle.at(1) );
  
  info()<< "Applying second Helicity Tool"<< endmsg;
  std::vector<TLorentzVector> PolarityBcands;
  PolarityBcands.push_back( BMeson );
  PolarityBcands.push_back( LambdaBary );
  PolarityBcands.push_back( DauProton );
  std::vector<double> PolarityFullAngle(2);
  PolarityAngleBMeson->initPhysics(PolarityBcands);
  PolarityAngleBMeson->computePhysics( &PolarityFullAngle);
  cosLpHelicity_3b.setVal(  PolarityFullAngle.at(1) );

  info()<< "Applying third Helicity Tool"<< endmsg;
  std::vector<TLorentzVector> PolarityPoleBcands;
  PolarityPoleBcands.push_back( BMeson );
  PolarityPoleBcands.push_back( LambdaBary );
  PolarityPoleBcands.push_back( Baryonium );
  AnglePoleBMeson->initPhysics(PolarityPoleBcands);
  std::vector<double> PoleFullAngle(2);
  AnglePoleBMeson->computePhysics( &PoleFullAngle);
  cosPoleBaryHelicity.setVal(  PoleFullAngle.at(1) );

  //Plot B zu Pole
  info()<< "Applying four Helicity Tool"<< endmsg;
  std::vector<TLorentzVector> PolarityBProton;
  PolarityBProton.push_back( Pion );
  PolarityBProton.push_back( LambdaBary );
  PolarityBProton.push_back( DauProton );
  AngleBP->initPhysics(PolarityBProton);
  std::vector<double> PBPoleAngle(2);
  AngleBP->computePhysics( &PBPoleAngle);
  cosPBPoleAngle.setVal( PBPoleAngle.at(1) );
  
  DauProton.Boost( -LambdaBary.BoostVector() );
  TLorentzVector LSpin( 0., 0., 0., 0.);
  
  if ( Lambda->particleID() == m_LambdaID ) LSpin.SetVect( - DauProton.Vect().Unit() );
  if ( Lambda->particleID() == m_LambdabarID ) LSpin.SetVect( DauProton.Vect().Unit() );
  
  BTag.setVal( p->particleID() == m_protonbarID ? -1. : 1. );

  TLorentzVector LzCopy(LambdaBary);
  
  LzCopy.Boost( -LambdaBary.BoostVector() );
  SanityL = LSpin.Dot( LzCopy );
  LSpin.Boost( -BMeson.BoostVector() );

  LambdaBary.Boost( -BMeson.BoostVector() );
  Proton.Boost( -BMeson.BoostVector() );
  
  E_L_Brest.setVal(LambdaBary.E());
    
  DauProton.Boost( -BMeson.BoostVector() );
  Pion.Boost( -BMeson.BoostVector() );
  BMeson.Boost( -BMeson.BoostVector() );
  
  SanityL.setVal( LSpin.Dot( LambdaBary ) );
  SanityB.setVal( BMeson.E() ) ;

  info() << " L Brestframe energy " << LambdaBary.E()
         << " P Brestframe energy " << Proton.E()
         << " Pi Brestframe energy " << Pion.E()
         << " B Brestframe energy " << BMeson.E() << endmsg;

  SpinTPCUnit.setVal( ( LSpin.Vect().Unit() ).Dot( LambdaBary.Vect().Unit().Cross( Pion.Vect().Unit() ) ) );
  SpinTPC.setVal( ( LSpin.Vect() ).Dot( LambdaBary.Vect().Cross( Pion.Vect() ) ) );

 
  return decayTreeFitter.status();
  
}
//=============================================================================
// EventInfo (to be implemented)
//=============================================================================
StatusCode MyDataAnalysis::EventInfo()
{
  const LHCb::ODIN* odin = getIfExists<LHCb::ODIN>(evtSvc(),LHCb::ODINLocation::Default);
  if ( !odin ) { odin = getIfExists<LHCb::ODIN>(evtSvc(),LHCb::ODINLocation::Default,false); }
  if ( !odin )
  {
    // should always be available ...
    return Error( "Cannot load the ODIN data object", StatusCode::SUCCESS );
  }
  const LHCb::RecSummary * tracks = getIfExists<LHCb::RecSummary>("/Event/Rec/Summary");
  int m_nTracks;
  if( tracks )
    m_nTracks = tracks->info(LHCb::RecSummary::nTracks,-1);
  else {
    m_nTracks = -1;
    return Error( "Cannot load the Rec data object", StatusCode::SUCCESS );
  }
  
  nTracks.setVal(m_nTracks);
  P_nTracks.setVal(m_nTracks);
  Pi_nTracks.setVal(m_nTracks);
  Lp_nTracks.setVal(m_nTracks);

  RunNumber.setVal(odin->runNumber());
  EventNumber.setVal(odin->eventNumber()) ;
  OdinTCK.setVal( odin->triggerConfigurationKey() );

  return StatusCode::SUCCESS;
}
//=============================================================================
// TruthMachting for sanity checks
//=============================================================================
bool MyDataAnalysis::truthMatch  (const LHCb::Particle& TopPart)
{
  debug() << "==> Begin Matching" << endmsg;
  const LHCb::MCParticle* MCTopPart = _assoc->relatedMCP(&TopPart, LHCb::MCParticleLocation::Default);
  if( !MCTopPart ) return false;
  if( MCTopPart->particleID() != TopPart.particleID() ) return false;
  debug() << "==> TopPartID " << MCTopPart->particleID().pid() << endmsg;
  LHCb::Particle::ConstVector daus = TopPart.daughtersVector();
  LHCb::Particle::ConstVector::const_iterator iDau;
  for( iDau = daus.begin(); iDau != daus.end(); ++iDau){
    if( !truthMatch( *(*iDau) ) ) return false;
    if( _assoc->relatedMCP( (*iDau), LHCb::MCParticleLocation::Default)->mother()->particleID() != TopPart.particleID() )
      return false;
  }
  info() << "==> True Decay found" << endmsg;
  return true;
} 
//===========================================================================
// Calculation of the TPC
//===========================================================================
StatusCode MyDataAnalysis::CalculateMCTPC(const LHCb::Particle *BCand)
{
  const LHCb::MCParticle *B(0), *Lambda(0), *pi(0), *Lp(0);
    
  B      = _assoc->relatedMCP(BCand,
                              LHCb::MCParticleLocation::Default);
  Lambda = _assoc->relatedMCP((BCand->daughtersVector().at(2)),
                              LHCb::MCParticleLocation::Default);
  Lp     = _assoc->relatedMCP((BCand->daughtersVector().at(2)->daughtersVector().at(0)),
                              LHCb::MCParticleLocation::Default);
  pi     = BCand->particleID() == m_BbarID ? _assoc->relatedMCP((BCand->daughtersVector().at(0)), LHCb::MCParticleLocation::Default)
                                           : _assoc->relatedMCP((BCand->daughtersVector().at(1)), LHCb::MCParticleLocation::Default);
  
  TLorentzVector BMeson    ( (Double_t)(B->momentum().Px()/1000.),
                             (Double_t)(B->momentum().Py()/1000.),
                             (Double_t)(B->momentum().Pz()/1000.),
                             (Double_t)(B->momentum().E()/1000.) );
  
  TLorentzVector LambdaBary( (Double_t)(Lambda->momentum().Px()/1000.),
                             (Double_t)(Lambda->momentum().Py()/1000.),
                             (Double_t)(Lambda->momentum().Pz()/1000.),
                             (Double_t)(Lambda->momentum().E()/1000.) );
  
  TLorentzVector Pion      ( (Double_t)(pi->momentum().Px()/1000.),
                             (Double_t)(pi->momentum().Py()/1000.),
                             (Double_t)(pi->momentum().Pz()/1000.),
                             (Double_t)(pi->momentum().E()/1000.) );
  
  TLorentzVector DauProton ( (Double_t)(Lp->momentum().Px()/1000.),
                             (Double_t)(Lp->momentum().Py()/1000.),
                             (Double_t)(Lp->momentum().Pz()/1000.),
                             (Double_t)(Lp->momentum().E()/1000.) );                 
  
  DauProton.Boost( -LambdaBary.BoostVector() );
  TLorentzVector LSpin( 0., 0., 0., 0.);
  
  if ( Lambda->particleID() == m_LambdaID ) LSpin.SetVect( - DauProton.Vect().Unit() );
  if ( Lambda->particleID() == m_LambdabarID ) LSpin.SetVect( DauProton.Vect().Unit() );
  
  LambdaBary.Boost( -LambdaBary.BoostVector() );
  LSpin.Boost( -BMeson.BoostVector() );
  LambdaBary.Boost( -BMeson.BoostVector() );
  DauProton.Boost( -BMeson.BoostVector() );
  Pion.Boost( -BMeson.BoostVector() );
  BMeson.Boost( -BMeson.BoostVector() );

  SpinTPCMC.setVal( ( LSpin.Vect().Unit() ).Dot( LambdaBary.Vect().Unit().Cross( Pion.Vect().Unit() ) ) );

  return StatusCode::SUCCESS;
}
