#ifndef IPARTICLEMANIPULATOR_H 
#define IPARTICLEMANIPULATOR_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"
#include "Event/Particle.h"
#include "Kernel/ParticleID.h"
#include "Kernel/ParticleProperty.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/IParticleReFitter.h"

// namespace LHCb { class Particle; }

static const InterfaceID IID_IParticleManipulator ( "IParticleManipulator", 1, 0 );

/** @class IParticleManipulator IParticleManipulator.h
 *  
 *
 *  @author Christian Voss
 *  @date   2012-11-26
 */
class IParticleManipulator : virtual public IAlgTool {
public: 

  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IParticleManipulator; }

  virtual StatusCode   doCorrection( LHCb::Particle* Part) = 0;
  virtual StatusCode undoCorrection( LHCb::Particle* Part) = 0;

protected:

private:

};
#endif // IPARTICLEMANIPULATOR_H
