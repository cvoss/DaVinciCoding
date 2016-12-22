#ifndef IMOMEMTUMCORRECTIONTOOL_H 
#define IMOMEMTUMCORRECTIONTOOL_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"

namespace LHCb { 
  class Particle ; 
  class Track ;
}

static const InterfaceID IID_iMomentumCorrectionTool ( "iMomentumCorrectionTool", 1, 0 );

/** @class iMomentumCorrectionTool iMomentumCorrectionTool.h
 *  
 *
 *  @author Christian Voss
 *  @date   2012-10-03
 */
class iMomentumCorrectionTool : virtual public IAlgTool {
public: 

  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_iMomentumCorrectionTool; }

  virtual void ParticleCorrectMomError( const LHCb::Particle* Cand) = 0;
  virtual void TrackCorrectMomError( const LHCb::Track* Trk) = 0;

protected:

private:

};
#endif // IMOMEMTUMCORRECTIONTOOL_H
