#ifndef MYSELECTIONANALYSIS_H 
#define MYSELECTIONANALYSIS_H 1

// Include files 
// from DaVinci
#include "Kernel/DaVinciAlgorithm.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "DecayTreeFitter/Fitter.h"
#include "TrackInterfaces/ITrackStateProvider.h"
#include "Kernel/ITriggerTisTos.h"

/** @class MySelectionAnalysis MySelectionAnalysis.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-05-19
 */

class MySelectionAnalysis : public DaVinciAlgorithm {

public: 
/// Standard constructor
  MySelectionAnalysis( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MySelectionAnalysis( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

  StatusCode BoostLWerner  (LHCb::Particle B, LHCb::Particle Lambda, TVector3 rLambda, double &angle);
  StatusCode ThetaAngleInB (LHCb::Particle B, LHCb::Particle Part, double &angle);
  StatusCode PreSelection  (const LHCb::Particle* BCand, 
                            double BMassLowEdge, double BMassHighEdge,
                            double LMassLowEdge, double LMassHighEdge);

  StatusCode Calculate     ();
  StatusCode PerformFit    ( const LHCb::Particle* BCand);
  StatusCode FilterTrigger ( const LHCb::Particle* IterB, bool doFilter, std::string DataPeriod );
  StatusCode LLCleanUp     ();
  StatusCode DDCleanUp     ();

  bool m_LLEvent;
  double m_log_pi_IP, m_log_p_IP, m_log_Lpi_IP, m_log_Lp_IP, m_log_BIP, m_log_LIP;
  double m_log_BVtxchi2;
  double m_log_B_DauSumchi2;
  double m_log_B_FLBchi2, m_log_Lambda_FDchi2, m_B_Eta;
  double m_B_PT, m_L_PT, m_p_PT, m_pi_PT, m_Lp_PT, m_Lpi_PT;
  double m_log_B_ctau, m_log_L_ctau, m_FitProb;
  double m_BWerner, m_LWerner;
  double m_ppiDist;
  double m_LAngle;
  double m_LMassVal;
  double m_LP_PT;
  double m_pointing;
  
  double m_BMassHighEdge, m_BMassLowEdge;
  double m_LMassHighEdge, m_LMassLowEdge;
  
  double m_p_pNNPID, m_p_KNNPID, m_p_piNNPID, m_pi_pNNPID, m_pi_KNNPID, m_pi_piNNPID,
    m_Lp_pNNPID, m_Lp_pDLLPID, m_Lp_KNNPID, m_Lp_piNNPID, m_Lpi_pNNPID, m_Lpi_KNNPID, m_Lpi_piNNPID;

  double m_pi_KDLLPID;
    
  const LHCb::Particle *m_B, *m_Lambda, *m_p, *m_pi, *m_Lp, *m_Lpi;
  LHCb::ParticleID m_BID, m_BbarID, m_BsbarID, 
    m_LambdaID, m_LambdabarID,
    m_protonID, m_protonbarID,
    m_pimID, m_pipID,
    m_KpID, m_KmID ;
  const LHCb::Vertex* m_BPV;
  
  const LHCb::IParticlePropertySvc *m_ppsvc;
  ITrackStateProvider *m_states;

  ITriggerTisTos *_L0TriggerTool;
  ITriggerTisTos *_HltTriggerTool;
  ITriggerTisTos *_Hlt1TriggerTool;
  ITriggerTisTos *_Hlt2TriggerTool;

  bool m_TriggerFilter;

private:

};
#endif // MYSELECTIONANALYSIS_H
