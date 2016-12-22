#ifndef WRITEGENMCKINEMATICS_H 
#define WRITEGENMCKINEMATICS_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciTupleAlgorithm.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "Kernel/IParticle2MCAssociator.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"

/** @class WriteGenMCKinematics WriteGenMCKinematics.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-07-29
 */
class WriteGenMCKinematics : public DaVinciTupleAlgorithm {
public: 
  /// Standard constructor
  WriteGenMCKinematics( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~WriteGenMCKinematics( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  StatusCode setColumn(Tuple& tuple, std::string candname, std::string data, double value);
  StatusCode WriteParticleTruth( std::string name, TLorentzVector& p4, Tuple& tuple);
  StatusCode WriteVertexTruth( std::string name, TVector3& r, Tuple& tuple);
  

  LHCb::ParticleID m_BID, m_BbarID, m_LID, m_LbarID, m_pID, m_pbarID, m_pipID, m_pimID, m_gammaID;
  
};
#endif // WRITEGENMCKINEMATICS_H
