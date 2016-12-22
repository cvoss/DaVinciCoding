#ifndef MYDVALGORITHM_H 
#define MYDVALGORITHM_H 1

// Include files
// from DaVinci, this is a specialized GaudiAlgorithm
#include "Kernel/DaVinciTupleAlgorithm.h"

#include "GaudiKernel/GenericAddress.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

#include "Kernel/IParticle2MCAssociator.h"
#include "Kernel/IEventTupleTool.h"
#include "Kernel/IParticleTupleTool.h"

#include "GaudiKernel/AlgFactory.h" 
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ObjectVector.h"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "Kernel/ParticleFilters.h"
#include "Kernel/IDecayTreeFit.h"
//#include "DecayTreeFitter/Fitter.h"

#include "LoKi/LoKi.h"
#include "LoKi/Child.h"
//#include "DecayTreeFitter/Fitter.h"
//#include "DecayTreeFitter/Tree.h"

#include "Kernel/IANNSvc.h"
#include "Event/L0DUReport.h"
#include "Event/HltDecReports.h"

/** @class MyDVAlgorithm MyDVAlgorithm.h
 *  
 *
 *  @author Christian Voss
 *  @date   2012-07-27
 */

using namespace boost::lambda;

class MyDVAlgorithm : public DaVinciTupleAlgorithm {
public: 

  MyDVAlgorithm ( const std::string& name, ISvcLocator* pSvcLocator )
    : DaVinciTupleAlgorithm ( name , pSvcLocator )
  { }
  
  //virtual ~MyDVAlgorithm( ); ///< Destructor

  //virtual StatusCode initialize();    ///< Algorithm initialization
  //virtual StatusCode execute   ();    ///< Algorithm execution
  //virtual StatusCode finalize  ();    ///< Algorithm finalization

protected:

  StatusCode   _setColumn   (Tuple& tuple,
                             std::string candname,
                             std::string data,
                             double value)
  {
    candname.append(data);
    tuple->column(candname,value);
    return StatusCode::SUCCESS;
  }

  StatusCode   _setColumn   (Tuple& tuple,
                             std::string candname,
                             std::string data,
                             int value)
  {
    candname.append(data);
    tuple->column(candname,value);
    return StatusCode::SUCCESS;
  }

  StatusCode _writeGeneratorInfo( std::string name, const LHCb::MCParticle* MCcand, Tuple& tuple )
  {
    if(MCcand)
      {
	_setColumn( tuple, name, "PDG", MCcand->particleID().pid() );
	debug() << "Mother Check" << endmsg;
        if(MCcand->mother() )_setColumn( tuple, name, "MotherID", MCcand->mother()->particleID().pid() );
        else _setColumn( tuple, name, "MotherID", 0 );
        _setColumn( tuple, name, "pMag",MCcand->p() );
        _setColumn( tuple, name, "pX", MCcand->momentum().px() );
        _setColumn( tuple, name, "pY", MCcand->momentum().py() );
        _setColumn( tuple, name, "pZ", MCcand->momentum().pz() );
        _setColumn( tuple, name, "pE", MCcand->momentum().E() );
        _setColumn( tuple, name, "GenMass", MCcand->virtualMass() );
	_setColumn( tuple, name, "Beta", MCcand->beta() );
	_setColumn( tuple, name, "Gamma", MCcand->gamma() );
	_setColumn( tuple, name, "BetaGamma", MCcand->betaGamma() );
	_setColumn( tuple, name, "Eta", MCcand->pseudoRapidity() );
	//_setColumn( tuple, name, "PrimVtxX", MCcand->primaryVertex()->position().x() );
	//_setColumn( tuple, name, "PrimVtxY", MCcand->primaryVertex()->position().y() );
	//_setColumn( tuple, name, "PrimVtxZ", MCcand->primaryVertex()->position().z() );
	debug() << "OriginVtx Check" << endmsg;
	if (MCcand->originVertex()){ 
	  _setColumn( tuple, name, "OriginVtxX", MCcand->originVertex()->position().x() );
	  _setColumn( tuple, name, "OriginVtxY", MCcand->originVertex()->position().y() );
	  _setColumn( tuple, name, "OriginVtxZ", MCcand->originVertex()->position().z() );
	}
	else {
	  _setColumn( tuple, name, "OriginVtxX", -1e99 );
	  _setColumn( tuple, name, "OriginVtxY", -1e99 );
	  _setColumn( tuple, name, "OriginVtxZ", -1e99 );
	}
      }
    else
      {
	_setColumn( tuple, name, "PDG", 0 );
        _setColumn( tuple, name, "MotherID", 0 );
        _setColumn( tuple, name, "pMag", -1e9 );
        _setColumn( tuple, name, "pX", -1e9 );
        _setColumn( tuple, name, "pY", -1e9 );
        _setColumn( tuple, name, "pZ", -1e9 );
        _setColumn( tuple, name, "pE", -1e9 );
        _setColumn( tuple, name, "GenMass", -1e9 );
	_setColumn( tuple, name, "Beta", -1e99 );
	_setColumn( tuple, name, "Gamma", -1e99 );
	_setColumn( tuple, name, "BetaGamma", -1e99 );
	_setColumn( tuple, name, "Eta", -1e99 );
	//_setColumn( tuple, name, "PrimVtxX", -1e99 );
	//_setColumn( tuple, name, "PrimVtxY", -1e99 );
	//_setColumn( tuple, name, "PrimVtxZ", -1e99 );
	_setColumn( tuple, name, "OriginVtxX", -1e99 );
	_setColumn( tuple, name, "OriginVtxY", -1e99 );
	_setColumn( tuple, name, "OriginVtxZ", -1e99 );
      }
    return StatusCode::SUCCESS;
  }
 
  StatusCode   _CandInfo    (std::string name,
                             const LHCb::Particle& part,
                             const LHCb::RecVertex::Range& _pvs,
                             Tuple& tuple,
                             bool useMC = false)
  {
    _setColumn( tuple, name, "PDG", (double)part.particleID().pid() );
    _setColumn( tuple, name, "Charge", (double)part.charge() );
    _setColumn( tuple, name, "pMag",part.p() );
    _setColumn( tuple, name, "pt",part.pt() );  
    _setColumn( tuple, name, "pX", part.momentum().px() );
    _setColumn( tuple, name, "pY", part.momentum().py() );
    _setColumn( tuple, name, "pZ", part.momentum().pz() );
    _setColumn( tuple, name, "pE", part.momentum().E() );
    if(part.proto())
    {
      _setColumn( tuple, name, "PIDe", part.proto()->info(LHCb::ProtoParticle::CombDLLe, -1e9) );
      _setColumn( tuple, name, "PIDp", part.proto()->info(LHCb::ProtoParticle::CombDLLp, -1e9) );
      _setColumn( tuple, name, "PIDk", part.proto()->info(LHCb::ProtoParticle::CombDLLk, -1e9) );
      _setColumn( tuple, name, "PIDmu", part.proto()->info(LHCb::ProtoParticle::CombDLLmu, -1e9) );
      _setColumn( tuple, name, "PIDpi", part.proto()->info(LHCb::ProtoParticle::CombDLLpi, -1e9) );
      _setColumn( tuple, name, "PIDNNe", part.proto()->info(LHCb::ProtoParticle::ProbNNe, -1e9) );
      _setColumn( tuple, name, "PIDNNp", part.proto()->info(LHCb::ProtoParticle::ProbNNp, -1e9) );
      _setColumn( tuple, name, "PIDNNk", part.proto()->info(LHCb::ProtoParticle::ProbNNk, -1e9) );
      _setColumn( tuple, name, "PIDNNmu", part.proto()->info(LHCb::ProtoParticle::ProbNNmu, -1e9) );
      _setColumn( tuple, name, "PIDNNpi", part.proto()->info(LHCb::ProtoParticle::ProbNNpi, -1e9) );
      //      _setColumn( tuple, name, "TrkChi2", LoKi::Cuts::TRPCHI2 ( &part ) );
      _setColumn( tuple, name, "TrkEta", part.proto()->track()->pseudoRapidity() );
      _setColumn( tuple, name, "TrkType", (double)part.proto()->track()->type() );
      _setColumn( tuple, name, "TrkX", part.proto()->track()->firstState().x() );
      _setColumn( tuple, name, "TrkY", part.proto()->track()->firstState().y() );
      _setColumn( tuple, name, "TrkZ", part.proto()->track()->firstState().z() );
      _setColumn( tuple, name, "TrkTX", part.proto()->track()->firstState().tx() );
      _setColumn( tuple, name, "TrkTY", part.proto()->track()->firstState().ty() );
      _setColumn( tuple, name, "TrkQoverP", part.proto()->track()->firstState().qOverP() );
      std::vector<LHCb::LHCbID> IDs = part.proto()->track()->lhcbIDs();
      std::vector<LHCb::LHCbID>::const_iterator iID;

      int isVelo(0);
      int isVeloR(0); 
      int isVeloPhi(0); 
      int isVeloPileUp(0); 
      int isVL(0); 
      int isVLR(0); 
      int isVLPhi(0); 
      int isIT(0); 
      int isOT(0); 
      int isTT(0); 
      int isST(0); 
      int isRich(0); 
      int isCalo(0); 
      int isMuon(0);
      for ( iID = IDs.begin(); iID != IDs.end(); ++iID)
	  {
	    if( iID->isVelo()) isVelo++; 
	    if( iID->isVeloR()) isVeloR++; 
	    if( iID->isVeloPhi()) isVeloPhi++; 
	    if( iID->isVeloPileUp()) isVeloPileUp++; 
	    // if( iID->isVL()) isVL++; 
	    // if( iID->isVLR()) isVLR++; 
	    // if( iID->isVLPhi()) isVLPhi++; 
	    if( iID->isIT()) isIT++; 
	    if( iID->isOT()) isOT++; 
	    if( iID->isTT()) isTT++; 
	    if( iID->isST()) isST++; 
	    if( iID->isRich()) isRich++; 
	    if( iID->isCalo()) isCalo++; 
	    if( iID->isMuon()) isMuon++; 
      
	  }
      _setColumn( tuple, name, "TrkisVelo", isVelo); 
      _setColumn( tuple, name, "TrkisVeloR", isVeloR); 
      _setColumn( tuple, name, "TrkisVeloPhi", isVeloPhi); 
      _setColumn( tuple, name, "TrkisVeloPileUp", isVeloPileUp); 
      // _setColumn( tuple, name, "TrkisVL", isVL); 
      // _setColumn( tuple, name, "TrkisVLR", isVLR); 
      // _setColumn( tuple, name, "TrkisVLPhi", isVLPhi); 
      _setColumn( tuple, name, "TrkisIT", isIT); 
      _setColumn( tuple, name, "TrkisOT", isOT); 
      _setColumn( tuple, name, "TrkisTT", isTT); 
      _setColumn( tuple, name, "TrkisST", isST); 
      _setColumn( tuple, name, "TrkisRich", isRich); 
      _setColumn( tuple, name, "TrkisCalo", isCalo); 
      _setColumn( tuple, name, "TrkisMuon", isMuon); 

  }
  
    debug() << "==> CandInfo: Momentum Cov Matrix" << endmsg;
    for(int i = 0 ; i < 4; ++i)
      for(int j = 0; j < 4; ++j)
      {
        std::stringstream out;
        out << "PxPyPzCov_" << i << "_" << j;
        _setColumn( tuple, name, out.str(), part.momCovMatrix()[i][j] );
      }
    if(part.proto())
      for(int i = 0 ; i < 5; ++i)
        for(int j = 0; j < 5; ++j)
        {
          std::stringstream out;
          out << "TrkCov_" << i << "_" << j;
          _setColumn( tuple, name, out.str(), part.proto()->track()->firstState().covariance()[i][j] );
        }
    else
      for(int i = 0 ; i < 5; ++i)
        for(int j = 0; j < 5; ++j)
        {
          std::stringstream out;
          out << "TrkCov_" << i << "_" << j; 
          _setColumn( tuple, name, out.str(), -1e9);
        }
    
    debug() << "==> CandInfo: IP Block" << endmsg;
    double IP, IPchi2;
    debug() << "==> CandInfo: Calculating min IP" << endmsg;
    _minIP( part, _pvs, IP, IPchi2);
    debug() << "==> CandInfo: Filling IP" << endmsg;
    _setColumn( tuple, name, "IP", IP) ;
    _setColumn( tuple, name, "IPchi2", IPchi2);
    
    debug() << "==> CandInfo: Mass Block" << endmsg;
    _setColumn( tuple, name, "Mass", part.measuredMass() );
    
    debug() << "==> CandInfo: Mass Error" << endmsg;
    _setColumn( tuple, name, "MassErr", part.measuredMassErr() );
    
    debug() << "==> CandInfo: Vertex Block" << endmsg;
    if( part.endVertex() ){
      _setColumn( tuple, name, "reducedchi2", part.endVertex()->chi2PerDoF() );
      _setColumn( tuple, name, "chi2", part.endVertex()->chi2() );
      _setColumn( tuple, name, "nDoF", (double)part.endVertex()->nDoF() );
      _setColumn( tuple, name, "VtxX", part.endVertex()->position().x() );
      _setColumn( tuple, name, "VtxY", part.endVertex()->position().y() );
      _setColumn( tuple, name, "VtxZ", part.endVertex()->position().z() );
      debug() << "==> CandInfo: Vertex Block Outgoing Tracks" << endmsg;
      //_setColumn( tuple, name, "VtxNTracks", (double)part.endVertex()->outgoingParticlesVector().size() );
      debug() << "==> CandInfo: Vertex Cov Matrix" << endmsg;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          std::stringstream out;
          out << "VtxCov_" << i << "_" << j;
          _setColumn( tuple, name, out.str(), part.endVertex()->covMatrix()[i][j] );
        }
    }
    
    name.append("MC");
    debug() << "==> CandInfo: MC Block" << endmsg;
    if(useMC){
      const LHCb::MCParticle* MCpart = _assoc->relatedMCP(&part, LHCb::MCParticleLocation::Default);
      _writeGeneratorInfo( name, MCpart, tuple );    
    }
    return StatusCode::SUCCESS;
  }
  
  std::size_t  _fetchList   (LHCb::ParticleID& type, 
                             LHCb::Particle::Range& alltracks, 
                             LHCb::Particle::ConstVector& seltracks)
  {
    LHCb::ParticleID atype = _ppsvc->find( type )->antiParticle()->particleID();
    
    std::size_t listlenpart = DaVinci::filter ( alltracks , 
                                                bind(&LHCb::Particle::particleID,_1) == type , 
                                                seltracks ) ;
    
    std::size_t listlenapart = DaVinci::filter ( alltracks , 
                                                 bind(&LHCb::Particle::particleID,_1) == atype , 
                                                 seltracks ) ;
    return listlenpart + listlenapart;
  }
  
  bool         _overlap     (const LHCb::Particle& part1,
                             const LHCb::Particle& part2)
  {
    if( !part1.endVertex() )
    {
      if( !part2.endVertex() )
        if( part1.proto() == part2.proto() ) return true;
      if( part2.endVertex() )
      {
        LHCb::Particle::ConstVector daughters = part2.daughtersVector() ;
        LHCb::Particle::ConstVector::const_iterator iDau;
        for ( iDau = daughters.begin(); iDau != daughters.end(); ++iDau)
        {
          if( (*iDau)->proto() == part1.proto() ) return true;
        }
      }
    }
    if( part1.endVertex() )
    {
      if( !part2.endVertex() )
      {
        LHCb::Particle::ConstVector daughters = part1.daughtersVector() ;
        LHCb::Particle::ConstVector::const_iterator iDau;
        for ( iDau = daughters.begin(); iDau != daughters.end(); ++iDau)
        {
          if( (*iDau)->proto() == part2.proto() ) return true;
        }
      }
      if( part2.endVertex() )
      {
        LHCb::Particle::ConstVector daughters1 = part1.daughtersVector() ;
        LHCb::Particle::ConstVector daughters2 = part2.daughtersVector() ;
        LHCb::Particle::ConstVector::const_iterator iDau1, iDau2;
        for ( iDau1 = daughters1.begin(); iDau1 != daughters1.end(); ++iDau1)
          for ( iDau2 = daughters2.begin(); iDau2 != daughters2.end(); ++iDau2)
          {
            if( (*iDau1)->proto() == (*iDau1)->proto() ) return true;
          }
      } 
    }
    return false;
  }
  
  StatusCode   _minIP       (const LHCb::Particle& part, 
                             const LHCb::RecVertex::Range& _pvs, 
                             double& IP, 
                             double& IPchi2)
  {
    LHCb::RecVertex::Range::iterator ivtx;
    double IP_temp, IPchi2_temp;
    double IP_save(1e100), IPchi2_save(1e100);
    debug() << "==> minIP: Starting Iteration for " << (int)_pvs.size() << " Vertices" << endmsg;

    for( ivtx = _pvs.begin() ; ivtx!=_pvs.end() ; ++ivtx)
    {
      debug() << "==> minIP: Calculation for " << (*ivtx) << " and " << &part << endmsg;
      distanceCalculator()->distance( &part, (*ivtx), IP_temp, IPchi2_temp);
      debug() << "==> minIP: Finished calculation for " << (*ivtx) << endmsg;
      if(IPchi2_temp < IPchi2_save ) 
      { 
        IPchi2_save = IPchi2_temp; 
        IP_save = IP_temp; 
      }
    }
    IP = IP_save;
    IPchi2 = IPchi2_save;
    debug() << "==> minIP: Finished" << endmsg;

    return StatusCode::SUCCESS;
  }
  
  int         _truthMatch  (const LHCb::Particle& TopPart)
  {
    const LHCb::MCParticle* MCTopPart = _assoc->relatedMCP(&TopPart, LHCb::MCParticleLocation::Default);
    if( !MCTopPart ) return 0;
    if( MCTopPart->particleID() != TopPart.particleID() ) return 0;
    LHCb::Particle::ConstVector daus = TopPart.daughtersVector();
    LHCb::Particle::ConstVector::const_iterator iDau;
    for( iDau = daus.begin(); iDau != daus.end(); ++iDau){
      if( !_truthMatch( *(*iDau) ) ) return 0;
      if( _assoc->relatedMCP( (*iDau), LHCb::MCParticleLocation::Default)->mother()->particleID() != TopPart.particleID() )
        return 0;
    }
    debug() << "==> True Decay found" << endmsg;
    return 1;
  } 
  
  StatusCode   _writePVTXs  (const LHCb::RecVertex::Range& _pvs,
                             LHCb::Particle& Comp,
                             Tuple& tuple)
  {
    std::ostringstream Name;
    uint i(1);
    _setColumn( tuple, "N", "PrimVtx",  (int)_pvs.size() );
    for( int ivtx = 0 ; ivtx < 15 ; ++ivtx)
    {
      Name << "PrimVtx" << i ;
      if( i <_pvs.size() )
      {
        double dist(0);
        distanceCalculator()->distance(_pvs.at(ivtx), Comp.endVertex(), dist);
        _setColumn( tuple, Name.str(), "reducedchi2",  _pvs.at(ivtx)->chi2PerDoF() );
        _setColumn( tuple, Name.str(), "chi2",         _pvs.at(ivtx)->chi2() );
        _setColumn( tuple, Name.str(), "nDoF", (double)_pvs.at(ivtx)->nDoF() );
        _setColumn( tuple, Name.str(), "X",            _pvs.at(ivtx)->position().x() );
        _setColumn( tuple, Name.str(), "Y",            _pvs.at(ivtx)->position().y() );
        _setColumn( tuple, Name.str(), "Z",            _pvs.at(ivtx)->position().z() );
        _setColumn( tuple, Name.str(), "NTracks",      (int)_pvs.at(ivtx)->tracks().size() );
        _setColumn( tuple, Name.str(), "DistToComp",      dist );
        
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
          {
            std::stringstream out;
            out << "VtxCov_" << i << "_" << j;
            _setColumn( tuple, Name.str(), out.str(), _pvs.at(ivtx)->covMatrix()[i][j] );
          }
      }
    else
    {
      _setColumn( tuple, Name.str(), "reducedchi2", -1e99 );
      _setColumn( tuple, Name.str(), "chi2", -1e99 );
      _setColumn( tuple, Name.str(), "nDoF", -1e99 );
      _setColumn( tuple, Name.str(), "X", -1e99 );
      _setColumn( tuple, Name.str(), "Y", -1e99 );
      _setColumn( tuple, Name.str(), "Z", -1e99 );
      _setColumn( tuple, Name.str(), "NTracks",-1e99);
      _setColumn( tuple, Name.str(), "DistToComp", -1e99 );
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          std::stringstream out;
          out << "VtxCov_" << i << "_" << j;
          _setColumn( tuple, Name.str(), out.str(), -1e99 );
        }
    }
      
      Name.str("");
      i++;
    }
    return StatusCode::SUCCESS;
  }

  const LHCb::IParticlePropertySvc *_ppsvc;
  IParticle2MCAssociator           *_assoc;
  IEventTupleTool                  *_EvtInfo;
  IEventTupleTool                  *_TriggerInfo;
  IEventTupleTool                  *_StrippingInfo;

private:
  
};
#endif // MYDVALGORITHM_H
