#ifndef MYMVATRAINING_H 
#define MYMVATRAINING_H 1

// Include files 
// from DaVinci.
#include "MySelectionAnalysis.h"

#include "IParticleManipulator.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "DecayTreeFitter/Fitter.h"
#include "TrackInterfaces/ITrackStateProvider.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

/** @class MyMVATraining MyMVATraining.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-05-19
 */
class MyMVATraining : public  MySelectionAnalysis{
public: 
  /// Standard constructor
  MyMVATraining( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~MyMVATraining( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

private:

  StatusCode DoTraining();
  StatusCode SaveRootTrees();
  
  TMVA::Factory *factory;
    
  std::vector<TTree*> Trees;
  
  int m_NTrain, m_NTest, m_BkgNTrain, m_BkgNTest;
  int DatanTrain, DatanTest, MCnTrain, MCnTest;
  
  Int_t m_L0HadronTos; 
  Int_t m_L0GlobalTis;
  Int_t m_Hlt1TrackAllL0Tos;
  Int_t m_Hlt2TopoBBDTTos;

  std::string   m_MCMatchLocation, m_MCPVLocation;
  std::string   m_DataLocation, m_DataPVLocation;
  std::string   m_outputFile;
  std::string   m_ClassifierBDT;
  std::string   m_ClassifierMLP;
  std::vector<std::string> m_PIDToolnames;
  bool m_DoTraining, m_save;
  std::string m_Location;

  double m_LambdaMass, m_BMass,m_BMassLowEdge, m_BMassHighEdge, m_LambdaMassWindow;

  std::vector< IParticleManipulator*> m_PIDTools;

  bool m_LLTraining;
  
};
#endif // MYMVATRAINING_H
