#ifndef BDTOBSSUBSTITUTION_H 
#define BDTOBSSUBSTITUTION_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "IParticleManipulator.h"            // Interface

/** @class BdToBsSubstitution BdToBsSubstitution.h
 *  
 *  @author Christian Voss
 *  @date   2013-10-28
 */
class BdToBsSubstitution : public GaudiTool, virtual public IParticleManipulator {
public: 
  /// Standard constructor
  BdToBsSubstitution( const std::string& type, 
                      const std::string& name,
                      const IInterface* parent);

  ~BdToBsSubstitution( ); ///< Destructor

  virtual StatusCode initialize();
  
  virtual StatusCode   doCorrection( LHCb::Particle* Part) ;
  virtual StatusCode undoCorrection( LHCb::Particle* Part) ;
  

protected:

private:

  double m_deltaMass;
  LHCb::IParticlePropertySvc *m_ppSvc;
  LHCb::ParticleID m_BID, m_BbarID, m_BsID, m_BsbarID;
};
#endif // BDTOBSSUBSTITUTION_H
