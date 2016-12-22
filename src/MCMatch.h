#ifndef MCMATCH_H 
#define MCMATCH_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciAlgorithm.h"
#include "Kernel/IParticle2MCAssociator.h"
#include "IParticleManipulator.h"

/** @class MCMatch MCMatch.h
 *
 *  @author Christian Voss
 *  @date   2013-05-21
 */
class MCMatch : public DaVinciAlgorithm {
public: 
  /// Standard constructor
  MCMatch( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MCMatch( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:
  IParticle2MCAssociator *_assoc;
  IParticleManipulator   *_PID;
  bool                   _truthMatch  (const LHCb::Particle& TopPart);
  bool                   m_deepMatch, m_match;
};
#endif // MCMATCH_H
