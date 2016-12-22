#ifndef HELICITYANGLE_H 
#define HELICITYANGLE_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IPhysicsComputation.h"            // Interface
#include "Event/Particle.h"
#include "TLorentzVector.h"
#include "TVector3.h"


/** @class HelicityAngle HelicityAngle.h
 *  
 *
 *  @author Christian Voss
 *  @date   2016-01-25
 */
class HelicityAngle : public GaudiTool, virtual public IPhysicsComputation {
public: 
  /// Standard constructor
  HelicityAngle( const std::string& type, 
                 const std::string& name,
                 const IInterface* parent);

  virtual ~HelicityAngle( ); ///< Destructor

  virtual StatusCode computePhysics(std::vector<double> *vars);
  virtual StatusCode initPhysics(const LHCb::Particle::ConstVector cands);
  virtual StatusCode initPhysics(const std::vector<TLorentzVector> cands);

protected:

private:

  TLorentzVector m_grandparent, m_parent, m_particle;
    
};
#endif // HELICITYANGLE_H
