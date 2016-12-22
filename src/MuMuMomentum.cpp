// Include files 

// from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MuMuMomentum.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MuMuMomentum
//
// 2012-08-29 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MuMuMomentum )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MuMuMomentum::MuMuMomentum( const std::string& name,
                            ISvcLocator* pSvcLocator)
  : MyDVAlgorithm ( name , pSvcLocator )
{
  declareProperty("JpsiMassWindow", _JpsiMassWindow = 400.*Gaudi::Units::MeV);  
  declareProperty("UpsilonMassWindow", _UpsilonMassWindow = 400.*Gaudi::Units::MeV);  
  declareProperty("ZMassWindow", _ZMassWindow = 30000.*Gaudi::Units::MeV);  
  declareProperty("PhiMassWindow", _PhiMassWindow = 400*Gaudi::Units::MeV);
  declareProperty("JpsiChi2", _JpsiChi2 = 50.);
  declareProperty("useMC", _useMC = true);
  declareProperty("DoJpsi", _doJpsi = false);
  declareProperty("DoZ0", _doZ0 = false);
  declareProperty("DoUpsilin", _doUpsilon = false);
  declareProperty("DoPhi", _doPhi = false);
}
//=============================================================================
// Destructor
//=============================================================================
MuMuMomentum::~MuMuMomentum() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MuMuMomentum::initialize() {
  StatusCode sc = MyDVAlgorithm::initialize(); 
  if ( sc.isFailure() ) return sc;

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;
  _ppsvc = ppSvc();
  //_assoc = tool<IParticle2MCAssociator>("MCMatchObjP2MCRelator");
  _EvtInfo = tool<IEventTupleTool>("TupleToolEventInfo");
  //_TriggerInfo = tool<IEventTupleTool>("TupleToolTrigger");
  //_TriggerInfo = tool<IEventTupleTool>("MyTupleTriggerTool");
  //_StrippingInfo = tool<IEventTupleTool>("TupleToolStripping");

  _corr = tool<iMomentumCorrectionTool>("MomentumCorrectionTool");

  const LHCb::ParticleProperty* mupInfo = _ppsvc->find( "mu+" ) ;
  const LHCb::ParticleProperty* munInfo = _ppsvc->find( "mu-" ) ;
  const LHCb::ParticleProperty* JpsiInfo = _ppsvc->find( "J/psi(1S)" ) ;
  const LHCb::ParticleProperty* UpsilonInfo = _ppsvc->find( "Upsilon(1S)" ) ;
  const LHCb::ParticleProperty* Z0Info = _ppsvc->find( "Z0" ) ;
  const LHCb::ParticleProperty* PhiInfo = _ppsvc->find("phi(1020)");
  
  _JpsiMass = JpsiInfo->mass();
  _JpsiID = JpsiInfo->particleID();
  _UpsilonID = UpsilonInfo->particleID();
  _UpsilonMass = UpsilonInfo->mass();
  _ZMass = Z0Info->mass();
  _ZID = Z0Info->particleID();
  _PhiMass = PhiInfo->mass();
  _PhiID = PhiInfo->particleID();
  _muonID = munInfo->particleID();
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MuMuMomentum::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  // code goes here  
  LHCb::Particle::Range tracks = this->particles();
  LHCb::Particle::ConstVector mup, mun;

  if ( msgLevel(MSG::DEBUG) )
    warning() << "==> Number of Muobs" << tracks.size() << endmsg;
  
  std::size_t Npion = _fetchMuons (_muonID, tracks, mup, mun);

  _JpsiReco(mup, mun);

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted. 
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MuMuMomentum::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return MyDVAlgorithm::finalize();
}

//=============================================================================

StatusCode MuMuMomentum::_JpsiReco(LHCb::Particle::ConstVector& mup,
                        LHCb::Particle::ConstVector& mun)
{
  LHCb::Particle *munCand(0), *mupCand(0);
    
  LHCb::Particle::ConstVector::const_iterator ip, in;
  for (ip = mup.begin(); ip != mup.end();++ip)
    for (in = mun.begin(); in != mun.end();++in)
    {
      if ( _overlap(*(*ip), *(*in)) ) continue;
      if ( ( (*ip)->charge() * (*in)->charge() ) > 0 ) continue;
      
      //distanceCalculator()->distance((*ip), (*in), mumudist);
      //if ( mumudist*Gaudi::Units::mm > _mumudist ) continue;
      
      Gaudi::LorentzVector Two = (*ip)->momentum() + (*in)->momentum();
      double TwoBodyMass = Two.M();

      LHCb::Vertex VtxTwoBody;
      LHCb::Particle Comp;

      if(_doJpsi){
	if(fabs(TwoBodyMass*Gaudi::Units::MeV - _JpsiMass*Gaudi::Units::MeV) < _JpsiMassWindow)
	  Comp.setParticleID(_JpsiID );
	else continue;
      }
      if(_doUpsilon){
        if (fabs(TwoBodyMass*Gaudi::Units::MeV - _UpsilonMass*Gaudi::Units::MeV) < _UpsilonMassWindow)
          Comp.setParticleID(_UpsilonID );
	else continue;
      }
      if(_doZ0){
	if (fabs(TwoBodyMass*Gaudi::Units::MeV - _ZMass*Gaudi::Units::MeV) < _ZMassWindow) 
	  Comp.setParticleID(_ZID );
	else continue;
      }
      if(_doPhi){
	if (fabs(TwoBodyMass*Gaudi::Units::MeV - _PhiMass*Gaudi::Units::MeV) < _PhiMassWindow) 
	  Comp.setParticleID(_PhiID );
	else continue;
      }
      if(!_doJpsi && !_doUpsilon && !_doZ0 && !_doPhi)
	return StatusCode::SUCCESS;

      mupCand = new LHCb::Particle(*(*ip));
      munCand = new LHCb::Particle(*(*in));      
      LHCb::Particle::ConstVector daus;
      daus.push_back(mupCand);
      daus.push_back(munCand);

      //particleCombiner("MomentumCombiner")->combine(daus,Comp,VtxTwoBody); 

      StatusCode scFit = vertexFitter()->fit(*(*ip), *(*in), VtxTwoBody, Comp);
      if (!scFit) continue;
      //if ( VtxTwoBody.chi2() > _JpsiChi2 ) continue;

      // Comp.setMomentum(Two);


      // Comp.addToDaughters(mupCand);
      // Comp.addToDaughters(munCand);
      
      const LHCb::RecVertex::Range _pvs = this->primaryVertices();
      Tuple tuple = nTuple( "JpsiData" ) ;
      _EvtInfo->fill(tuple);
      _DocaPull(*munCand, *mupCand, VtxTwoBody, tuple);
      _CandInfo("Comp", Comp, _pvs, tuple, _useMC);      
      _CandInfo("mup", *mupCand, _pvs, tuple, _useMC);
      _CandInfo("mun", *munCand, _pvs, tuple, _useMC);

      _corr->ParticleCorrectMomError(mupCand);
      _corr->ParticleCorrectMomError(munCand);
      _CandInfo("mupCorr", *mupCand, _pvs, tuple, _useMC);
      _CandInfo("munCorr", *munCand, _pvs, tuple, _useMC);

      if(_useMC) tuple->column("Truth", _truthMatch(Comp) );
      counter("J/psi Candidates")++;
        
      const LHCb::Track* Trk = mupCand->proto()->track();
      tuple->column("TrkTest_before", Trk->firstState().covariance()[4][4]);
      _corr->TrackCorrectMomError( Trk );
      tuple->column("TrkTest_after", Trk->firstState().covariance()[4][4]);

      tuple->write();
        
    }  
  delete munCand; delete mupCand;
  
  return StatusCode::SUCCESS;
  
}

//====================================================================================

std::size_t MuMuMomentum::_fetchMuons (LHCb::ParticleID& type,
                                 LHCb::Particle::Range& alltracks,
                                 LHCb::Particle::ConstVector& mup,
                                 LHCb::Particle::ConstVector& mun)
{
  LHCb::ParticleID atype = _ppsvc->find( type )->antiParticle()->particleID();
  std::size_t listlenpart = DaVinci::filter ( alltracks ,
                                              bind(&LHCb::Particle::particleID,_1) == type ,
                                              mun ) ;
  std::size_t listlenapart = DaVinci::filter ( alltracks ,
                                               bind(&LHCb::Particle::particleID,_1) == atype ,
                                               mup ) ;
  
  return listlenpart + listlenapart;
  
}

//====================================================================================

StatusCode MuMuMomentum::_DocaPull(LHCb::Particle& cand1, LHCb::Particle& cand2,
				  LHCb::Vertex& Vtx, Tuple& tuple)

{
  Gaudi::XYZVector dx3 = cand1.proto()->track()->firstState().position() -
    cand2.proto()->track()->firstState().position() ;
  Gaudi::XYZVector n3  = cand1.proto()->track()->firstState().slopes().Cross(
			     cand2.proto()->track()->firstState().slopes() 
									     ) ;
  double doca = dx3.Dot(n3) / n3.R() ;

  double doca_pull = std::sqrt(Vtx.chi2()) * (doca>0 ? 1 : -1) ; 

  _setColumn( tuple, "Doca", "Value", doca );
  _setColumn( tuple, "Doca", "Pull", doca_pull );

  return StatusCode::SUCCESS;
}

//=====================================================================================
