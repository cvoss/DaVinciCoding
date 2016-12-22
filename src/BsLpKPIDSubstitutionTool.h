#ifndef BSLPKPIDSUBSTITUTIONTOOL_H 
#define BSLPKPIDSUBSTITUTIONTOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IParticleManipulator.h"            // Interface
#include "Event/Particle.h"
#include "Kernel/ParticleID.h"
#include "Kernel/ParticleProperty.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/IParticleReFitter.h"

/** @class BsLpKPIDSubstitutionTool BsLpKPIDSubstitutionTool.h
 *  
 *
 *  @author Christian Voss
 *  @date   2013-07-29
 */
class BsLpKPIDSubstitutionTool : public GaudiTool, virtual public IParticleManipulator {
public: 
  /// Standard constructor
  BsLpKPIDSubstitutionTool( const std::string& type, 
                            const std::string& name,
                            const IInterface* parent);

  ~BsLpKPIDSubstitutionTool( ); ///< Destructor

  virtual StatusCode initialize();

  virtual StatusCode   doCorrection( LHCb::Particle* Part) ;
  virtual StatusCode undoCorrection( LHCb::Particle* Part) ;

protected:

private:

  double m_pMass, m_piMass, m_KMass;
  IParticleReFitter *m_reFit;
  LHCb::IParticlePropertySvc *m_ppSvc;
  LHCb::ParticleID m_BID, m_BbarID, m_LbID, m_LbbarID, m_pID, m_pbarID, m_pimID, m_pipID, m_KmID, m_KpID;

};
#endif // BSLPKPIDSUBSTITUTIONTOOL_H
