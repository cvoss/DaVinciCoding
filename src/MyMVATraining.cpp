// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MyMVATraining.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MyMVATraining
//
// 2014-05-19 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MyMVATraining )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MyMVATraining::MyMVATraining( const std::string& name,
                              ISvcLocator* pSvcLocator)
  : MySelectionAnalysis ( name , pSvcLocator )
  , Trees()
  , DatanTrain(0), DatanTest(0), MCnTrain(0),MCnTest(0)
  , m_L0HadronTos(0), m_L0GlobalTis(0), m_Hlt1TrackAllL0Tos(0), m_Hlt2TopoBBDTTos(0)
{
  declareProperty("MCMatchLocation", m_MCMatchLocation = "");
  declareProperty("MCPVLocation", m_MCPVLocation = "");
  declareProperty("DataLocation", m_DataLocation = "");
  declareProperty("DataPVLocation", m_DataPVLocation = "");
  declareProperty("BMassLowEdge", m_BMassLowEdge = 5400.*Gaudi::Units::MeV);
  declareProperty("BMassHighEdge", m_BMassHighEdge = 5600.*Gaudi::Units::MeV);
  declareProperty("LambdaMassLowEdge", m_LMassLowEdge = 1100.*Gaudi::Units::MeV);
  declareProperty("LambdaMassHighEdge", m_LMassHighEdge = 1125.*Gaudi::Units::MeV);
  declareProperty("NSignalTrainEvents", m_NTrain = 1000);
  declareProperty("NSignalTestEvents", m_NTest = 1000);
  declareProperty("NBkgTrainEvents", m_BkgNTrain = 1000);
  declareProperty("NBkgTestEvents", m_BkgNTest = 1000);
  declareProperty("ROOTFile", m_outputFile = "TMVA.root");
  declareProperty("BDTClassifierName", m_ClassifierBDT = "BDT");
  declareProperty("MLPClassifierName", m_ClassifierMLP = "MLP");
  declareProperty("PIDToolnames", m_PIDToolnames, "List of PID tools");
  declareProperty("DoTraining", m_DoTraining = true);
  declareProperty("SaveRootTrees", m_save = false);
  declareProperty("TreeLocation", m_Location = "Void/of/time/and/Space");
  declareProperty("LLTraining", m_LLTraining = true);
  
  Trees.push_back(new TTree("tSTrainLL","tSTrainLL"));
  Trees.push_back(new TTree("tSTestLL","tSTestLL"));
  Trees.push_back(new TTree("tBTrainLL","tBTrainLL"));
  Trees.push_back(new TTree("tBTestLL","tBTestLL"));

  for(std::vector<TTree*>::iterator iT = Trees.begin(); iT != Trees.end(); ++iT)
  {
    (*iT)->SetDirectory(0);
    (*iT)->Branch("log_pi_IP",&m_log_pi_IP,"log_pi_IP/D");
    (*iT)->Branch("log_p_IP",&m_log_p_IP,"log_p_IP/D");
    (*iT)->Branch("log_Lpi_IP",&m_log_Lpi_IP,"log_Lpi_IP/D");
    (*iT)->Branch("log_Lp_IP",&m_log_Lp_IP,"log_Lp_IP/D");
    (*iT)->Branch("log_BIP", &m_log_BIP, "log_BIP/D");
    (*iT)->Branch("log_LIP", &m_log_LIP, "log_LIP/D");
    (*iT)->Branch("log_BVtxchi2",&m_log_BVtxchi2,"log_BVtxchi2/D");
    (*iT)->Branch("FitProb",&m_FitProb,"FitProb/D");
    (*iT)->Branch("log_B_DauSumchi2",&m_log_B_DauSumchi2,"log_B_DauSumchi2/D");
    (*iT)->Branch("log_B_FLBchi2",&m_log_B_FLBchi2,"log_B_FLBchi2/D");
    (*iT)->Branch("log_Lambda_FDchi2",&m_log_Lambda_FDchi2,"log_Lambda_FDchi2/D");
    (*iT)->Branch("B_Eta",&m_B_Eta,"B_Eta/D");
    (*iT)->Branch("B_PT",&m_B_PT,"B_PT/D");
    (*iT)->Branch("log_B_ctau",&m_log_B_ctau,"log_B_ctau/D");
    (*iT)->Branch("log_L_ctau",&m_log_L_ctau,"log_L_ctau/D");
    (*iT)->Branch("BWerner", &m_BWerner, "BWerner/D");
    (*iT)->Branch("LWerner", &m_LWerner, "LWerner/D");
    (*iT)->Branch("ppiDist", &m_ppiDist, "ppiDist/D");
    (*iT)->Branch("LAngle", &m_LAngle,"LAngle/D");
    (*iT)->Branch("LP_PT", &m_LP_PT,"LP_PT/D");
    (*iT)->Branch("Pointing", &m_pointing,"Pointing/D");
    (*iT)->Branch("L0HadronTos", &m_L0HadronTos , "L0HadronTos/I");
    (*iT)->Branch("L0GlobalTis", &m_L0GlobalTis , "L0GlobalTis/I");
    (*iT)->Branch("Hlt1TrackAllL0Tos", &m_Hlt1TrackAllL0Tos , "Hlt1TrackAllL0Tos/I");
    (*iT)->Branch("Hlt2TopoBBDTTos", &m_Hlt2TopoBBDTTos , "Hlt2TopoBBDTTos/I");   
  }
}
//=============================================================================
// Destructor
//=============================================================================
MyMVATraining::~MyMVATraining() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MyMVATraining::initialize() {
  StatusCode sc = MySelectionAnalysis::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  for (unsigned int i = 0; i < m_PIDToolnames.size(); i++)
    m_PIDTools.push_back( tool<IParticleManipulator>( m_PIDToolnames.at(i),this) );
 
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MyMVATraining::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.

  if(DatanTrain == m_BkgNTrain - 1 
     && DatanTest == m_BkgNTest - 1
     && MCnTrain == m_NTrain - 1
     && MCnTest == m_NTest - 1 )
  {
    debug() << "Skipping further Events" << endmsg;
    return StatusCode::SUCCESS;
  }
  
  LHCb::Particle::Range MCMatchparticles = getIfExists<LHCb::Particle::Range>( m_MCMatchLocation ); 
  LHCb::RecVertex::Range MCPvs = getIfExists<LHCb::RecVertex::Range>(m_MCPVLocation);

  LHCb::Particle::Range Dataparticles    = getIfExists<LHCb::Particle::Range>( m_DataLocation ); 
  LHCb::RecVertex::Range DataPvs = getIfExists<LHCb::RecVertex::Range>(m_DataPVLocation);

  LHCb::Particle::ConstVector::const_iterator IterB;
  std::vector<IParticleManipulator*>::iterator iTool;
  std::vector<IParticleManipulator*>::reverse_iterator riTool;
  
  debug() << "# of MC Events: " << MCMatchparticles.size() << endmsg;
  
  if(MCnTrain != m_NTrain && MCnTest != m_NTest)
    if(MCMatchparticles.size())
      if(MCPvs.size()){
        debug() << "Numbers " << MCnTrain << " " << MCnTest << endmsg;      
        for (IterB = MCMatchparticles.begin(); IterB != MCMatchparticles.end(); IterB++){
          if(MCnTrain == m_NTrain -1 && MCnTest == m_NTest - 1 )
            return StatusCode::SUCCESS;
          for ( iTool = m_PIDTools.begin(); iTool != m_PIDTools.end(); ++iTool)
            (*iTool) -> doCorrection(const_cast<LHCb::Particle*>(*IterB));
          if(MCnTrain < m_NTrain - 1){
            debug() << "Filering" << endmsg;
            if (FilterTrigger( *IterB, true, "2012b") == StatusCode::SUCCESS ) 
              if ( PreSelection( *IterB, 5000., 5800., m_LMassLowEdge, m_LMassHighEdge) == StatusCode::SUCCESS )
                if ( PerformFit( *IterB ) == StatusCode::SUCCESS )
                  if ( ( m_LLTraining ? LLCleanUp() : DDCleanUp() ) == StatusCode::SUCCESS )
                {
                  Trees[0]->Fill();
                  if( MCnTrain%100==0 )
                    info() << "Fill MC Train Trees " << MCnTrain << endmsg;  
                  MCnTrain++;
                }
            for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
              (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
            continue;
          }
          debug() << "Numbers " << MCnTrain << " " << MCnTest << endmsg;
          if(MCnTest < m_NTest - 1){
            debug() << "Filering" << endmsg;
            if (FilterTrigger( *IterB , true, "2012b") == StatusCode::SUCCESS ) 
              if ( PreSelection( *IterB, 5000., 5800., m_LMassLowEdge, m_LMassHighEdge) == StatusCode::SUCCESS )
                if ( PerformFit( *IterB ) == StatusCode::SUCCESS )
                  if ( ( m_LLTraining ? LLCleanUp() : DDCleanUp() ) == StatusCode::SUCCESS )
                  {
                    Trees[1]->Fill();
                    if( MCnTest%100==0 )
                      info() << "Fill MC Test Trees " << MCnTest << endmsg;
                    MCnTest++;
                  }
           }
           for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
             (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));  
        }
      }

  debug() << "# of Data Events: " << Dataparticles.size() << endmsg;
  
  if(DatanTrain != m_BkgNTrain && DatanTest != m_BkgNTest)
    if(Dataparticles.size())
      if(DataPvs.size()){
        for (IterB = Dataparticles.begin(); IterB != Dataparticles.end(); IterB++){
          if(DatanTrain == m_NTrain -1 && DatanTest == m_NTest - 1 ) return StatusCode::SUCCESS;
          
          for ( iTool = m_PIDTools.begin(); iTool != m_PIDTools.end(); ++iTool)
            (*iTool) -> doCorrection(const_cast<LHCb::Particle*>(*IterB));
          if(DatanTrain < m_BkgNTrain - 1){
            if (FilterTrigger( *IterB, true, "2012b") == StatusCode::SUCCESS ) 
              if (PreSelection( *IterB, m_BMassLowEdge, m_BMassHighEdge, m_LMassLowEdge, m_LMassHighEdge) == StatusCode::SUCCESS)
                if ( PerformFit( *IterB ) == StatusCode::SUCCESS )
                  if ( ( m_LLTraining ? LLCleanUp() : DDCleanUp() ) == StatusCode::SUCCESS )
                    if ( m_p_pNNPID > 0.2
                         && m_p->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5 &&
                         m_pi->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5)
                    {
                      Trees[2]->Fill();
                      if( DatanTrain%100==0 )
                        info() << "Fill Data Train Trees " << DatanTrain << endmsg;  
                      DatanTrain++;
                    }
            for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
              (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB)); 
            continue;
          }
          if(DatanTest < m_BkgNTest - 1){
             if (FilterTrigger( *IterB, true, "2012b") == StatusCode::SUCCESS ) 
              if (PreSelection( *IterB, m_BMassLowEdge, m_BMassHighEdge, m_LMassLowEdge, m_LMassHighEdge) == StatusCode::SUCCESS)
                if ( PerformFit( *IterB ) == StatusCode::SUCCESS )
                  if ( ( m_LLTraining ? LLCleanUp() : DDCleanUp() ) == StatusCode::SUCCESS )
                    if ( m_p_pNNPID > 0.2 
                         && m_p->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5 &&
                         m_pi->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5)
                    {
                      Trees[3]->Fill();
                      if( DatanTest%100==0 )
                        info() << "Fill Data Test Trees " << DatanTest << endmsg;
                      DatanTest++;
                    }
          }
          for ( riTool = m_PIDTools.rbegin(); riTool != m_PIDTools.rend(); ++riTool)
            (*riTool) -> undoCorrection(const_cast<LHCb::Particle*>(*IterB));
          
        } 
      }
  
  debug() << "Event Voll" << endmsg; 
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MyMVATraining::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  if( m_DoTraining ) MyMVATraining::DoTraining();
  if( m_save ) MyMVATraining::SaveRootTrees();
  return MySelectionAnalysis::finalize();  // must be called after all other actions
}
//=============================================================================
// Configure Training
//=============================================================================
StatusCode MyMVATraining::DoTraining()
{
  info() << "==> Start TMVAClassification" << endmsg;
  
  TFile* outputFile = TFile::Open( m_outputFile.c_str(), "RECREATE" );

  factory = new TMVA::Factory( "TMVAClassification", outputFile,
                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;:AnalysisType=Classification" );

  info() << "define Variables" << endmsg;
  
  factory->AddVariable( "log_pi_IP", 'F' );
  factory->AddVariable( "log_p_IP", 'F' );
  // if ( m_LLTraining ) 
  factory->AddVariable( "log_Lpi_IP", 'F' );
  // if ( m_LLTraining ) 
  factory->AddVariable( "log_Lp_IP", 'F' );
  factory->AddVariable( "log_BIP", 'F');
  // if ( m_LLTraining ) 
  factory->AddVariable( "log_LIP", 'F');
  factory->AddVariable( "log_BVtxchi2", 'F');
  factory->AddVariable( "FitProb", 'F');
  factory->AddVariable( "log_B_DauSumchi2", 'F');
  factory->AddVariable( "log_B_FLBchi2", 'F');
  // if ( m_LLTraining ) 
  factory->AddVariable( "log_Lambda_FDchi2", 'F');
  factory->AddVariable( "B_PT", 'F');
  factory->AddVariable( "B_Eta", 'F' );
  factory->AddVariable( "log_B_ctau", 'F');
  // if ( m_LLTraining ) 
  factory->AddVariable( "log_L_ctau", 'F');
  factory->AddVariable( "BWerner", 'F');
  // if ( m_LLTraining ) 
  factory->AddVariable( "LWerner", 'F');
  factory->AddVariable( "ppiDist", 'F');
  // if ( m_LLTraining ) 
  factory->AddVariable( "LAngle" , 'F');
  // factory->AddVariable( "LP_PT" , 'F');
  factory->AddVariable( "Pointing", 'F');
  
  info() << "Add trees" << endmsg;
  
  factory->AddSignalTree( Trees[0], 2.9e-6, "Training" );
  factory->AddSignalTree( Trees[1], 2.9e-6,  "Test" );

  factory->AddBackgroundTree( Trees[2], 1, "Training" );
  factory->AddBackgroundTree( Trees[3], 1,  "Test" );

  info() << "Book BDT " << endmsg;
  
  // ---- Select methods
  factory->BookMethod( TMVA::Types::kBDT, m_ClassifierBDT.c_str(),
                     "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=30:NNodesMax=3");

  std::string test2(m_ClassifierBDT);

  factory->BookMethod( TMVA::Types::kBDT, (test2.append(std::string("Ada"))).c_str(),
                     "!H:!V:NTrees=850:MinNodeSize=10.0%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30" );

  test2.clear();
  test2.append(m_ClassifierBDT);

  factory->BookMethod( TMVA::Types::kBDT, (test2.append(std::string("Bagging"))).c_str(),
                           "!H:!V:NTrees=800:MaxDepth=3:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

  test2.clear();
  test2.append(m_ClassifierBDT);
  
  factory->BookMethod( TMVA::Types::kBDT, (test2.append(std::string("Fisher"))).c_str(),
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );
 
  test2.clear();
  test2.append(m_ClassifierBDT);
  
  factory->BookMethod( TMVA::Types::kBDT,(test2.append(std::string("Forest"))).c_str(),
                     "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=30:NNodesMax=3:UseRandomisedTrees");
  info() << "Book MLP" << endmsg;

  std::string test1(m_ClassifierMLP);
  factory->BookMethod( TMVA::Types::kMLP, test1.c_str(),
                       "H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N+5:TestRate=5:UseRegulator" );
  test1.clear();
  test1.append(m_ClassifierMLP);

  factory->BookMethod( TMVA::Types::kMLP, (test1.append(std::string("Uni"))).c_str(),
                       "H:!V:NeuronType=tanh:VarTransform=U:NCycles=500:HiddenLayers=N+1:TestRate=5:UseRegulator" );
  test1.clear();
  test1.append(m_ClassifierMLP);
  factory->BookMethod( TMVA::Types::kMLP, (test1.append(std::string("UniLayers"))).c_str(),
                        "H:!V:NeuronType=tanh:VarTransform=Norm:NCycles=500:HiddenLayers=N-10,N-11,N-12:TestRate=5:UseRegulator" );


  // ---- Now you can tell the factory to train, test, and evaluate the MVAs
  
  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  
  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  info() << "==> TMVAClassification is done!" << endmsg;

  outputFile->Close();

  return StatusCode::SUCCESS;
}
//=============================================================================
// Save Trees
//=============================================================================
StatusCode MyMVATraining::SaveRootTrees(){

  TFile Output( m_Location.c_str(), "recreate");
  Output.cd();
  for(std::vector<TTree*>::iterator iT = Trees.begin(); iT != Trees.end(); ++iT)
    (*iT)->Write();
  Output.Close();
  
  return StatusCode::SUCCESS;
}
