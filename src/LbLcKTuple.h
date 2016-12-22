#ifndef LBLCKTUPLE_H 
#define LBLCKTUPLE_H 1

// Include files 
// from DaVinci.
#include "MyDVAlgorithm.h"

/** @class LbLcKTuple LbLcKTuple.h
 *  
 *
 *  @author Christian Voss
 *  @date   2016-07-11
 */
class LbLcKTuple : public MyDVAlgorithm {
public: 
  /// Standard constructor
  LbLcKTuple( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~LbLcKTuple( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

};
#endif // LBLCKTUPLE_H
