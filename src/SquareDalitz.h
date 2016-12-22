#ifndef SQUAREDALITZ_H 
#define SQUAREDALITZ_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciHistoAlgorithm.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "Kernel/IParticle2MCAssociator.h"
#include "IPhysicsComputation.h"

/** @class SquareDalitz SquareDalitz.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-07-28
 */
class SquareDalitz : public DaVinciHistoAlgorithm {
public: 
  /// Standard constructor
  SquareDalitz( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~SquareDalitz( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  Double_t coshel(TLorentzVector particle,
                  TLorentzVector parent, 
                  TLorentzVector grandparent);

  LHCb::ParticleID m_BID, m_BbarID, m_LID, m_LbarID, m_pID, m_pbarID, m_pipID, m_pimID;

  IPhysicsComputation *AngularDalitz;
};
#endif // SQUAREDALITZ_H
