#ifndef AMENTEROSPODOLANSKI_H 
#define AMENTEROSPODOLANSKI_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IPhysicsComputation.h"            // Interface
#include "Event/Particle.h"
#include "TLorentzVector.h"
#include "TVector3.h"


/** @class AmenterosPodolanski AmenterosPodolanski.h
 *  
 *
 *  @author Christian Voss
 *  @date   2015-01-16
 */
class AmenterosPodolanski : public GaudiTool, virtual public IPhysicsComputation {
public: 
  /// Standard constructor
  AmenterosPodolanski( const std::string& type, 
                       const std::string& name,
                       const IInterface* parent);

  virtual ~AmenterosPodolanski( ); ///< Destructor

  virtual StatusCode initialize();
  virtual StatusCode computePhysics(std::vector<double> *vars);
  virtual StatusCode initPhysics(const LHCb::Particle::ConstVector cands);
  virtual StatusCode initPhysics(const std::vector<TLorentzVector> cands);
  
protected:

private:

  TVector3 p0,pplus,pminus;
  
};
#endif // AMENTEROSPODOLANSKI_H
