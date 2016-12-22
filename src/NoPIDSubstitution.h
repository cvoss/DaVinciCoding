#ifndef NOPIDSUBSTITUTION_H 
#define NOPIDSUBSTITUTION_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IParticleManipulator.h"            // Interface


/** @class NoPIDSubstitution NoPIDSubstitution.h
 *
 *  @author Christian Voss
 *  @date   2013-10-24
 */
class NoPIDSubstitution : public GaudiTool, virtual public IParticleManipulator {
public: 
  /// Standard constructor
  NoPIDSubstitution( const std::string& type, 
                     const std::string& name,
                     const IInterface* parent) ;

  ~NoPIDSubstitution( ) ; ///< Destructor

  StatusCode initialize();

  virtual StatusCode   doCorrection( LHCb::Particle* Part) ;

  virtual StatusCode undoCorrection( LHCb::Particle* Part) ;

protected:

private:

};
#endif // NOPIDSUBSTITUTION_H
