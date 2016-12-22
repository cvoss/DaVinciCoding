#ifndef MYREDUCEDATA_H 
#define MYREDUCEDATA_H 1

// Include files 
// from DaVinci.
#include "MySelectionAnalysis.h"

#include "Kernel/ITriggerTisTos.h"
#include "Event/HltDecReports.h"
#include "IParticleManipulator.h"
#include "Event/RecSummary.h"
#include "DecayTreeFitter/Fitter.h"
#include "TrackInterfaces/ITrackStateProvider.h"
#include "Kernel/ILHCbMagnetSvc.h"
#include "TMath.h"

/** @class MyReduceData MyReduceData.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-05-26
 */
class MyReduceData : public MySelectionAnalysis {
public: 
  /// Standard constructor
  MyReduceData( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MyReduceData( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:
  StatusCode FilterMVA    ( );
  StatusCode LcVeto       ( bool doLcVeto );
  
  std::vector<IParticleManipulator*> m_PIDTools;

  std::vector<std::string> m_PIDToolnames;
  std::string m_dataPeriod;
  double m_p_pNNPIDCut, m_MVALLValCut, m_MVADDValCut;
  bool m_doLcVeto;
  
};
#endif // MYREDUCEDATA_H
