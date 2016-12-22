#ifndef MUMUMOMENTUM_H 
#define MUMUMOMENTUM_H 1

// Include files
// from DaVinci, this is a specialized GaudiAlgorithm
#include "MyDVAlgorithm.h"
#include "iMomentumCorrectionTool.h"

/** @class MuMuMomentum MuMuMomentum.h
 *  
 *
 *  @author Christian Voss
 *  @date   2012-08-29
 */
class MuMuMomentum : public MyDVAlgorithm {
public: 
  /// Standard constructor
  MuMuMomentum( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MuMuMomentum( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  StatusCode   _JpsiReco    (LHCb::Particle::ConstVector& mup,
                             LHCb::Particle::ConstVector& mun);

  std::size_t  _fetchMuons  (LHCb::ParticleID& type,
                             LHCb::Particle::Range& alltracks,
                             LHCb::Particle::ConstVector& mup,
                             LHCb::Particle::ConstVector& mun);
  StatusCode   _DocaPull    (LHCb::Particle& cand1,
			     LHCb::Particle& cand2,
			     LHCb::Vertex& Vtx, Tuple& tuple);

  double       _JpsiMassWindow, _JpsiChi2;
  double       _UpsilonMassWindow, _UpsilonChi2;
  double       _ZMassWindow, _ZChi2;
  double       _PhiMassWindow, _PhiChi2;
  double       _JpsiMass, _UpsilonMass, _ZMass, _PhiMass;
  bool         _useMC, _doJpsi, _doZ0, _doUpsilon, _doPhi;

  LHCb::ParticleID  _muonID, _amuonID, _JpsiID, _UpsilonID, _ZID, _PhiID;

  iMomentumCorrectionTool *_corr;
};
#endif // MUMUMOMENTUM_H
