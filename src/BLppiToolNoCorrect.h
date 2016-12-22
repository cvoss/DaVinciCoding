#ifndef BLPPITOOLNOCORRECT_H 
#define BLPPITOOLNOCORRECT_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IParticleManipulator.h"            // Interface


/** @class BLppiToolNoCorrect BLppiToolNoCorrect.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-07-09
 */
class BLppiToolNoCorrect : public GaudiTool, virtual public IParticleManipulator {
public: 
  /// Standard constructor
  BLppiToolNoCorrect( const std::string& type, 
                      const std::string& name,
                      const IInterface* parent);

  virtual ~BLppiToolNoCorrect( ); ///< Destructor
  
  virtual StatusCode initialize();

  virtual StatusCode   doCorrection( LHCb::Particle* Part) ;
  virtual StatusCode undoCorrection( LHCb::Particle* Part) ;

protected:

private:

  double m_pMass, m_piMass;
  IParticleReFitter *m_reFit;
  LHCb::IParticlePropertySvc *m_ppSvc;
  LHCb::ParticleID m_BID, m_BbarID, m_LbID, m_LbbarID, m_pID, m_pbarID, m_pimID, m_pipID;

};
#endif // BLPPITOOLNOCORRECT_H
