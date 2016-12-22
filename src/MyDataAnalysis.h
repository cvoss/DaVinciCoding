#ifndef MYDATAANALYSIS_H 
#define MYDATAANALYSIS_H 1

// Include files 
// from DaVinci.
#include "MySelectionAnalysis.h"
#include  <random>
#include  <iterator>

#include "Event/HltDecReports.h"
#include "IParticleManipulator.h"
#include "Event/RecSummary.h"
#include "DecayTreeFitter/Fitter.h"
#include "TrackInterfaces/ITrackStateProvider.h"
#include "Kernel/ILHCbMagnetSvc.h"
#include "TFile.h"
#include "TROOT.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "IPhysicsComputation.h"
#include "Event/ODIN.h"
#include "Kernel/IParticle2MCAssociator.h"

/** @class MyDataAnalysis MyDataAnalysis.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-05-19
 */

class MyDataAnalysis : public MySelectionAnalysis {
public: 
  /// Standard constructor
  MyDataAnalysis( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MyDataAnalysis( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  StatusCode TriggerWriteOut(const LHCb::Particle *IterB );
  StatusCode WriteDataSet();
  StatusCode FillDataSet();
  StatusCode constrainBMass(const LHCb::Particle *BCand);
  StatusCode EventInfo();
  StatusCode CalculateMCTPC(const LHCb::Particle *BCand);
  bool       truthMatch  (const LHCb::Particle& TopPart);
      
  template<typename Iter, typename RandomGenerator>
  Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
  }
  
  template<typename Iter>
  Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
  }
  
  ILHCbMagnetSvc* m_magfieldsvc;
  
  std::vector<IParticleManipulator*> m_PIDTools;
  IParticle2MCAssociator *_assoc;
  
  IPhysicsComputation *AngularDalitz, *ArmenterosPodolanski,
    *HelicityAngle, *PolarityAngleEnhancement, *PolarityAngleBMeson, *AnglePoleBMeson, *AngleBP;
  
  RooRealVar LMass,BMass,BMassDTF,BMassDTF_P,BMassDTF_PT,BMassDTF_Eta,p_pNNPID,p_KNNPID,
    p_piNNPID,pi_pNNPID,pi_KNNPID,pi_piNNPID,Lp_pNNPID, Lp_pDLLPID, Lp_KNNPID,Lp_piNNPID,Lpi_pNNPID,
    Lpi_KNNPID,Lpi_piNNPID,mLP,mPpi,mLpi,mPLp,mpiLpi,Lp_TrkType,Lpi_TrkType,MVALL,MVADD,Polarity;

  RooRealVar log_pi_IP, log_p_IP, log_Lp_IP, log_Lpi_IP,log_LIP, log_BIP, log_BVtxchi2, FitProb, log_B_DauSumchi2,
    log_B_FLBchi2, log_Lambda_FDchi2, B_PT, B_Eta, log_B_ctau, log_L_ctau,
    WernerB, WernerL, ppi_Dist, LAngle,LP_PT, Pointing;

  RooRealVar L_PT, p_PT, pi_PT, Lp_PT, Lpi_PT;
  RooRealVar DLcosTheta, DLPhiL, DLPhihh;
  RooRealVar DpcosTheta, DpPhiL, DpPhihh;
  RooRealVar Lambda_Phi;

  RooRealVar BMassMC,mLPMC,mPpiMC,mLpiMC, mPPMC, mpipiMC, AP_pt, AP_alpha, DTF2Status;
  
  RooRealVar P_P, P_TRACK_Eta, P_nTracks;
  RooRealVar Pi_P, Pi_TRACK_Eta, Pi_nTracks;
  RooRealVar Lp_P, Lp_TRACK_Eta, Lp_nTracks;
 	RooRealVar nTracks;

  RooRealVar RunNumber, EventNumber, OdinTCK;
  RooRealVar DataPeriod;

  RooRealVar pi_KDLLPID;
  
  RooRealVar L0HadronTos; 
  RooRealVar L0GlobalTis;
  RooRealVar Hlt1TrackAllL0Tos;
  RooRealVar Hlt2Topo2BodyBBDTTos;
  RooRealVar Hlt2Topo3BodyBBDTTos;
  RooRealVar Hlt2Topo4BodyBBDTTos;
  RooRealVar Hlt2Topo2BodySimpleTos;
  RooRealVar Hlt2Topo3BodySimpleTos;
  RooRealVar Hlt2Topo4BodySimpleTos;

  RooRealVar SanityB, SanityL, SpinTPCUnit, SpinTPC, SpinTPCMC, BTag,
    pHelicity, cospHelicity, cosLpHelicity_Enh, cosLpHelicity_3b, cosPBPoleAngle;
  RooRealVar E_L_Brest,cosPoleBaryHelicity;
  RooRealVar nCand;
  
  RooArgSet *Set;
  RooDataSet *rawdata;
  
  std::string RootFile, m_dataPeriod;
  std::vector<std::string> m_PIDToolnames;
  //std::vector<double> m_MVAinputValues;

  bool m_sPlot, m_checks, m_mctuple;
  int m_nTracks;
  double m_LMass, m_BMass, m_pMass, m_piMass;
  
};
#endif // MYDATAANALYSIS_H
