#ifndef MULTCAND_H 
#define MULTCAND_H 1

// Include files 
// from DaVinci.
#include "Kernel/DaVinciHistoAlgorithm.h"
#include "IParticleManipulator.h"
#include  <random>
#include  <iterator>

/** @class MultCand MultCand.h
 *  
 *
 *  @author Christian Voss
 *  @date   2015-02-03
 */

class MultCand : public DaVinciHistoAlgorithm {
public: 
  /// Standard constructor
  MultCand( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MultCand( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  std::vector<IParticleManipulator*> m_PIDTools;
  std::vector<std::string> m_PIDToolnames;

  template<typename Iter, typename RandomGenerator>
  Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
  }
  
  template<typename Iter>
  Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
  }
  
};
#endif // MULTCAND_H
