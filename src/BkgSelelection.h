#ifndef BKGSELELECTION_H 
#define BKGSELELECTION_H 1

// Include files 
// from DaVinci.
#include "MySelectionAnalysis.h"
#include "IParticleManipulator.h"

/** @class BkgSelelection BkgSelelection.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-11-17
 */
class BkgSelelection : public MySelectionAnalysis {
public: 
  /// Standard constructor
  BkgSelelection( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~BkgSelelection( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  bool m_SelectLL;
  std::vector<std::string> m_PIDToolnames;
  std::vector<IParticleManipulator*> m_PIDTools;
  double m_p_pNNPIDCut;

};
#endif // BKGSELELECTION_H
