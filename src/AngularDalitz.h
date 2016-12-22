#ifndef ANGULARDALITZ_H 
#define ANGULARDALITZ_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IPhysicsComputation.h"            // Interface
#include "Event/Particle.h"
#include "TLorentzVector.h"
#include "TVector3.h"

/** @class AngularDalitz AngularDalitz.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-07-31
 */
class AngularDalitz : public GaudiTool, virtual public IPhysicsComputation {
public: 
  /// Standard constructor
  AngularDalitz( const std::string& type, 
                 const std::string& name,
                 const IInterface* parent);

  ~AngularDalitz( ); ///< Destructor

  virtual StatusCode initialize();
  virtual StatusCode computePhysics(std::vector<double> *vars);
  virtual StatusCode initPhysics(const LHCb::Particle::ConstVector cands);
  virtual StatusCode initPhysics(const std::vector<TLorentzVector> cands);
  

protected:

private:

  TVector3 m_x,m_y,m_z,m_u,m_v,m_w,m_Lambda;
  
  
};
#endif // ANGULARDALITZ_H
