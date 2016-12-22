#ifndef KSBACKGROUND_H 
#define KSBACKGROUND_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciTupleAlgorithm.h"
#include "IPhysicsComputation.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
/** @class KSBackground KSBackground.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-06-23
 */
class KSBackground : public DaVinciTupleAlgorithm {
public: 
  /// Standard constructor
  KSBackground( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~KSBackground( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:
   const LHCb::IParticlePropertySvc *m_ppsvc;
  IPhysicsComputation *ArmenterosPodolanski;

  int NKS, NLZ;
  
};
#endif // KSBACKGROUND_H
