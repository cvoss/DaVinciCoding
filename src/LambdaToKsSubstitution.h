#ifndef LAMBDATOKSSUBSTITUTION_H 
#define LAMBDATOKSSUBSTITUTION_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IParticleManipulator.h"            // Interface


/** @class LambdaToKsSubstitution LambdaToKsSubstitution.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-06-23
 */
class LambdaToKsSubstitution : public GaudiTool, virtual public IParticleManipulator {
public: 
  /// Standard constructor
  LambdaToKsSubstitution( const std::string& type, 
                          const std::string& name,
                          const IInterface* parent);

  virtual ~LambdaToKsSubstitution( ); ///< Destructor

  virtual StatusCode initialize();

  virtual StatusCode   doCorrection( LHCb::Particle* particle );
  virtual StatusCode undoCorrection( LHCb::Particle* particle );

protected:

private:

  double _piMass, _pMass;
  IParticleReFitter* _reFit;
  LHCb::IParticlePropertySvc* _ppSvc;
  LHCb::ParticleID _pipID, _pimID, _pID, _pbarID;

};
#endif // LAMBDATOKSSUBSTITUTION_H
