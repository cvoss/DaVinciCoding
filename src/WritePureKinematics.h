#ifndef WRITEPUREKINEMATICS_H 
#define WRITEPUREKINEMATICS_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciTupleAlgorithm.h"
#include "IParticleManipulator.h"

/** @class WritePureKinematics WritePureKinematics.h
 *  
 *
 *  @author Christian Voss
 *  @date   2015-02-05
 */
class WritePureKinematics : public DaVinciTupleAlgorithm {
public: 
  /// Standard constructor
  WritePureKinematics( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~WritePureKinematics( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:
  std::vector<IParticleManipulator*> m_PIDTools;
  
  StatusCode setColumn              (Tuple& tuple, std::string candname, std::string data, double value);
  StatusCode WriteParticleKinematic ( std::string name, const LHCb::Particle& part, Tuple& tuple);

  std::vector<std::string> m_PIDToolnames;
  LHCb::ParticleID m_pID, m_pbarID;
};
#endif // WRITEPUREKINEMATICS_H
