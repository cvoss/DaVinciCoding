#ifndef IPHYSICSCOMPUTATION_H 
#define IPHYSICSCOMPUTATION_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"
#include "Event/Particle.h"

static const InterfaceID IID_IPhysicsComputation ( "IPhysicsComputation", 1, 0 );

/** @class IPhysicsComputation IPhysicsComputation.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-07-31
 */
class TLorentzVector;

class IPhysicsComputation : virtual public IAlgTool {
public: 

  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IPhysicsComputation; }

  virtual StatusCode computePhysics(std::vector<double> *vars) = 0 ;
  virtual StatusCode initPhysics(const LHCb::Particle::ConstVector cands) = 0;
  virtual StatusCode initPhysics(const std::vector<TLorentzVector> cands) = 0;

protected:

private:

};
#endif // IPHYSICSCOMPUTATION_H
