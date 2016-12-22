#ifndef PITOKSUBSTITUTION_H 
#define PITOKSUBSTITUTION_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IParticleManipulator.h"            // Interface


/** @class PiToKSubstitution PiToKSubstitution.h
 *  
 *
 *  @author Christian Voss
 *  @date   2013-10-29
 */
class PiToKSubstitution : public GaudiTool, virtual public IParticleManipulator {
public: 
  /// Standard constructor
  PiToKSubstitution( const std::string& type, 
                     const std::string& name,
                     const IInterface* parent);
  ~PiToKSubstitution( ); ///< Destructor

  virtual StatusCode initialize();
  
  virtual StatusCode   doCorrection( LHCb::Particle* Part) ;
  virtual StatusCode undoCorrection( LHCb::Particle* Part) ;

protected:

private:
  double m_piMass, m_KMass;
  IParticleReFitter *m_reFit;
  LHCb::IParticlePropertySvc *m_ppSvc;
  LHCb::ParticleID m_pimID, m_pipID, m_KmID, m_KpID;
};
#endif // PITOKSUBSTITUTION_H
