#ifndef BDPP_H 
#define BDPP_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciHistoAlgorithm.h"
#include "Kernel/IParticle2MCAssociator.h"
#include "Kernel/IEventTupleTool.h"
#include "IParticleManipulator.h"

//#include <boost/lambda/lambda.hpp>
//#include <boost/lambda/bind.hpp>

#include "GaudiKernel/MsgStream.h"

#include "Kernel/ITriggerTisTos.h"

#include "TFile.h"
#include "TROOT.h"
#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooArgusBG.h"
#include "RooChebychev.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "TH1.h"

#include "RooStats/SPlot.h"

class BDpp : public DaVinciHistoAlgorithm {
public: 
  /// Standard constructor
  BDpp( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~BDpp( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  std::string m_PIDToolname;
  int m_Bins;
  bool _useMC;
  double _BMass, _LbMass, _pMass, _piMass;
  LHCb::ParticleID _BID, _BbarID, _LbID, _LbbarID, _pID, _pbarID, _pimID, _pipID, _LID;

  const LHCb::IParticlePropertySvc *_ppsvc;
  IParticle2MCAssociator           *_assoc;
  IParticleManipulator             *_PID;
};
#endif // BDPP_H
