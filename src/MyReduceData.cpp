// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MyReduceData.h"

#include "ReadMLPLL2012Uni.h"
#include "ReadMLPDD2012.h"
#include "ReadMLPLL2011Uni.h"
#include "ReadMLPDD2011Uni.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MyReduceData
//
// 2014-05-26 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MyReduceData )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MyReduceData::MyReduceData( const std::string& name,
                            ISvcLocator* pSvcLocator)
  : MySelectionAnalysis ( name , pSvcLocator )
{
  declareProperty("BMassLowEdge", m_BMassLowEdge = 4950.*Gaudi::Units::MeV);
  declareProperty("BMassHighEdge", m_BMassHighEdge = 5650.*Gaudi::Units::MeV);
  declareProperty("MassLowEdge", m_LMassLowEdge = 1105.*Gaudi::Units::MeV);
  declareProperty("LMassHighEdge", m_LMassHighEdge = 1125.*Gaudi::Units::MeV);
  declareProperty("MVALLSelection", m_MVALLValCut = -2. );
  declareProperty("MVADDSelection", m_MVADDValCut = -2. );
  declareProperty("ProtonPID", m_p_pNNPIDCut = -1.);
  declareProperty("PIDToolnames", m_PIDToolnames, "List of PID tools");
  declareProperty("FilterTrigger", m_TriggerFilter = true);
  declareProperty("DataPeriod", m_dataPeriod = "2012b" );
  declareProperty("doLcVeto", m_doLcVeto = false);
}
//=============================================================================
// Destructor
//=============================================================================
MyReduceData::~MyReduceData() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MyReduceData::initialize() {
  StatusCode sc = MySelectionAnalysis::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;
  
  if ( (m_dataPeriod != "2011") && (m_dataPeriod != "2012a") && (m_dataPeriod != "2012b") )
    fatal() << "Data Period not correctly configured " << m_dataPeriod << endmsg;

  // _L0TriggerTool = tool<ITriggerTisTos>("L0TriggerTisTos",this);
  // _HltTriggerTool = tool<ITriggerTisTos>("TriggerTisTos",this);
  for (unsigned int i = 0; i < m_PIDToolnames.size(); i++)
    m_PIDTools.push_back( tool<IParticleManipulator>( m_PIDToolnames.at(i),this) );

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MyReduceData::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.  

  LHCb::Particle::Range Bs = this->particles();
  
  debug() << "Size of Inputs " << Bs.size() << endmsg;
  if ( 0 == Bs.size() )
    warning() << "Empty particle container" << endmsg;
  
  LHCb::Particle::ConstVector::const_iterator IterB;
  
  std::vector<IParticleManipulator*>::iterator iTool;
  std::vector<IParticleManipulator*>::reverse_iterator riTool;
  
  for (IterB = Bs.begin(); IterB != Bs.end(); IterB++)
  {
    for ( iTool = m_PIDTools.begin(); iTool != m_PIDTools.end(); ++iTool)
      (*iTool) -> doCorrection(const_cast<LHCb::Particle*>(*IterB));
    
    if ( FilterTrigger( *IterB, m_TriggerFilter, m_dataPeriod ) == StatusCode::SUCCESS )
      if (PreSelection( *IterB, m_BMassLowEdge, m_BMassHighEdge, m_LMassLowEdge, m_LMassHighEdge) == StatusCode::SUCCESS)
        if ( PerformFit( *IterB ) == StatusCode::SUCCESS )
        {
          bool LLCase(false);
          if ( (int)m_Lp->proto()->track()->type() == 3 || (int)m_Lpi->proto()->track()->type() == 3 )
            LLCase = true;
          if ( ( LLCase ? LLCleanUp() : DDCleanUp() ) == StatusCode::SUCCESS )
            if ( m_p->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5 &&
                 m_pi->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) < 0.5)
              if ( FilterMVA() == StatusCode::SUCCESS )
                if ( LcVeto( m_doLcVeto ) == StatusCode::SUCCESS )
                {
                  counter("Candidates")++;
                  info() << " Canidate found " << endmsg;
                  cloneAndMarkTree(*IterB);
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
StatusCode MyReduceData::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return MySelectionAnalysis::finalize();  // must be called after all other actions
}

//=================================================================================
// Calculate and filter the MVA response
//=================================================================================
StatusCode MyReduceData::FilterMVA( )
{
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
  
  ReadMLPLL2012Uni    MLPSelLL2012    (theMVAInputVars);
  ReadMLPDD2012       MLPSelDD2012    (theMVAInputVars);

  double MVALLVal(-2), MVADDVal(-2);
  
  if ( (int)m_Lp->proto()->track()->type() == 3 || (int)m_Lpi->proto()->track()->type() == 3 )
  {
    if ( m_dataPeriod == "2012b" )
      MVALLVal = MLPSelLL2012.GetMvaValue(m_MVAinputValues); 
    if ( m_dataPeriod == "2011" )
      MVALLVal = MLPSelLL2011.GetMvaValue(m_MVAinputValues);
  }

  if ( (int)m_Lp->proto()->track()->type() == 5 || (int)m_Lpi->proto()->track()->type() == 5 )
  {
    if ( m_dataPeriod == "2012b" )
      MVADDVal = MLPSelDD2012.GetMvaValue(m_MVAinputValues); 
    if ( m_dataPeriod == "2011" )
      MVADDVal = MLPSelDD2011.GetMvaValue(m_MVAinputValues); 
  }

  if( m_p_pNNPID <  m_p_pNNPIDCut )
    return StatusCode::FAILURE;

  if( MVADDVal == -2 && MVALLVal < m_MVALLValCut )
    return StatusCode::FAILURE;
  if( MVALLVal == -2 && MVADDVal < m_MVADDValCut )
    return StatusCode::FAILURE;
  
  return StatusCode::SUCCESS;
}
//=================================================================================
// LC Mass veto
//=================================================================================
StatusCode MyReduceData::LcVeto( bool doLcVeto = false )
{
  if ( !doLcVeto )
    return StatusCode::SUCCESS;

  double MassLpi( (m_pi->momentum() + m_Lambda->momentum()).M() );
  if ( (int)m_Lp->proto()->track()->type() == 3 || (int)m_Lpi->proto()->track()->type() == 3 )
    if ( (MassLpi<2266.)|| (MassLpi>2306.) )
      return StatusCode::SUCCESS;
  if ( (int)m_Lp->proto()->track()->type() == 5 || (int)m_Lpi->proto()->track()->type() == 5 )
    if ( (MassLpi<2271.)|| (MassLpi>2302.) )
      return StatusCode::SUCCESS;
  
  return StatusCode::FAILURE; 
}
