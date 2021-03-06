#ifndef MYA2MUMUANALYSIS_H 
#define MYA2MUMUANALYSIS_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciAlgorithm.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "IPhysicsComputation.h"
#include "Event/ODIN.h"
#include "Kernel/IParticle2MCAssociator.h"
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "TrackInterfaces/ITrackStateProvider.h"


#include "GaudiKernel/AlgFactory.h" 
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ObjectVector.h"

#include "Kernel/ParticleFilters.h"

/** @class MyA2MuMuAnalysis MyA2MuMuAnalysis.h
 *  
 *
 *  @author Christian Voss
 *  @date   2016-12-16
 */
using namespace boost::lambda;

class MyA2MuMuAnalysis : public DaVinciAlgorithm {
public: 
  /// Standard constructor
  MyA2MuMuAnalysis( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MyA2MuMuAnalysis( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  std::size_t fetchParts           ( LHCb::ParticleID& type,
                                     LHCb::Particle::Range& alltracks,
                                     LHCb::Particle::ConstVector& parts1,
                                     LHCb::Particle::ConstVector& parts2);
  StatusCode checkforBMeson        ( LHCb::Particle::ConstVector& B,
                                     LHCb::ParticleID& BHypo,
                                     LHCb::Particle& BMeson);
  
  StatusCode checkforMuonOverlap   ( const LHCb::Particle *A,
                                     StatusCode BzTag,
                                     StatusCode BzbTag,
                                     StatusCode BplTag,
                                     StatusCode BmTag );
  
  StatusCode calculateMomError     ( const LHCb::Particle& part, double& Error);
  
  const LHCb::IParticlePropertySvc *m_ppsvc;
  ITrackStateProvider *m_states;

  RooRealVar AMass, AMMass, nA, nBz, nBc, mu_mu_NNPID, K_KNNPID, BVertexChi2, DVertexChi2, BDIRA, K_IPchi2,
    Mu1_probChi2, Mu2_probChi2, Mu1_pErr, Mu2_pErr, Mu1_Ghost, Mu2_Ghost, Mu_Overlap, BTagMap;
  RooDataSet *rawdata;
  RooArgSet *Set;

  std::string m_B_zLocation, m_Bpl_Location, m_A_Location, m_DataPVLocation, RootFile;
  std::vector<double> m_DMassWindow;
  double m_Mu_PID_Cut, m_D0_VertexCut, m_K_PID_Cut, m_B_VertexCut, m_BDIRA, m_KIPchi2;

  LHCb::ParticleID  m_BzID, m_BzbarID, m_BplusID, m_BminusID, m_DzID, m_DzbarID,
    m_DplusID, m_DminusID, m_pimID, m_pipID, m_KpID, m_KmID, m_mumID, m_mupID;
  LHCb::Particle TagBz, TagBzbar, TagBplus, TagBminus;
  
  std::size_t NeutralBs, ChargedBs;
  
};
#endif // MYA2MUMUANALYSIS_H
