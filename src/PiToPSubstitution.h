#ifndef PITOPSUBSTITUTION_H
#define PITOPSUBSTITUTION_H

#include "GaudiAlg/GaudiTool.h"
#include "IParticleManipulator.h"


class PiToPSubstitution : public GaudiTool, virtual public IParticleManipulator
{
public:

  PiToPSubstitution( const std::string& type,
		      const std::string& name,
		      const IInterface* parent );
  ~PiToPSubstitution();

  virtual StatusCode initialize();

  virtual StatusCode   doCorrection( LHCb::Particle* particle );
  virtual StatusCode undoCorrection( LHCb::Particle* particle );
  
private:

  double _piMass, _pMass;
  IParticleReFitter* _reFit;
  LHCb::IParticlePropertySvc* _ppSvc;
  LHCb::ParticleID _pipID, _pimID, _pID, _pbarID;
};

#endif // PITOPSUBSTITUTION_H
