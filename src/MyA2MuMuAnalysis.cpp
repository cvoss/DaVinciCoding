// Include files 

 // from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "MyA2MuMuAnalysis.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MyA2MuMuAnalysis
//
// 2016-12-16 : Christian Voss
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( MyA2MuMuAnalysis )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MyA2MuMuAnalysis::MyA2MuMuAnalysis( const std::string& name,
                                    ISvcLocator* pSvcLocator)
  : DaVinciAlgorithm ( name , pSvcLocator )
  , AMass("AMass","AMass",14000,56000)
  , AMMass("AMMass","AMMass",14000,56000)
  , nA("nA","nA",1.0,1000.)
  , nBz("nBz","nBz",1.0,1000.)
  , nBc("nBc","nBc",1.0,1000.)
  , mu_mu_NNPID("mu_muNNPID","mu_muNNPID",-0.5,1.5)
  , K_KNNPID("K_KNNPID","K_KNNPID",-0.5,1.5)
  , BVertexChi2("BVertexChi2","BVertexChi2",-1.,1e99)
  , DVertexChi2("DVertexChi2","DVertexChi2",-1.,1e99)
  , BDIRA("BDIRA","BDIRA",-1.,4.)
  , K_IPchi2("K_IPchi2","K_IPchi2",-1,1e99)
  , Mu1_probChi2("Mu1_probChi2","Mu1_probChi2",-0.1,2.)
  , Mu2_probChi2("Mu2_probChi2","Mu2_probChi2",-0.1,2.)
  , Mu1_pErr("Mu1_pErr","Mu1_pErr",-0.1,2.)
  , Mu2_pErr("Mu2_pErr","Mu2_pErr",-0.1,2.)
  , Mu1_Ghost("Mu1_Ghost","Mu1_Ghost",-0.1,2.)
  , Mu2_Ghost("Mu2_Ghost","Mu2_Ghost",-0.1,2.)
  , Mu_Overlap("Mu_Overlap","Mu_Overlap",-2.,2.)
  , BTagMap("BTagMap","BTagMap",0.,50.)
{
  declareProperty("B_zLocation", m_B_zLocation = "");
  declareProperty("Bpl_Location", m_Bpl_Location = "");
  declareProperty("A_Location", m_A_Location = "");
  declareProperty("DataPVLocation", m_DataPVLocation = "");
  declareProperty("DMassCut", m_DMassWindow, "D Meson Mass Window");
  declareProperty("MuonPID", m_Mu_PID_Cut);
  declareProperty("KaonPID", m_K_PID_Cut);
  declareProperty("D0Vertex", m_D0_VertexCut );
  declareProperty("BVertex", m_B_VertexCut );
  declareProperty("BDIRA", m_BDIRA = 0.998 );
  declareProperty("KIPchi2", m_KIPchi2 = 12.);
  declareProperty("RootFile", RootFile = "A2MuMuData.root");
}
//=============================================================================
// Destructor
//=============================================================================
MyA2MuMuAnalysis::~MyA2MuMuAnalysis() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode MyA2MuMuAnalysis::initialize() {
  StatusCode sc = DaVinciAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  m_ppsvc = ppSvc();
  m_states = tool<ITrackStateProvider>("TrackStateProvider",this);
  const LHCb::ParticleProperty* BzInfo     = m_ppsvc->find( "B0" ) ;
  const LHCb::ParticleProperty* BplInfo    = m_ppsvc->find( "B+" ) ;
  const LHCb::ParticleProperty* DzInfo     = m_ppsvc->find( "D0" ) ;
  const LHCb::ParticleProperty* DplInfo    = m_ppsvc->find( "D+" ) ;
  const LHCb::ParticleProperty* pipInfo    = m_ppsvc->find( "pi+" );
  const LHCb::ParticleProperty* pimInfo    = m_ppsvc->find( "pi-" );
  const LHCb::ParticleProperty* KpInfo     = m_ppsvc->find( "K+" );  
  const LHCb::ParticleProperty* KmInfo     = m_ppsvc->find( "K-" );
  const LHCb::ParticleProperty* mupInfo    = m_ppsvc->find( "mu+" ) ;
  const LHCb::ParticleProperty* munInfo    = m_ppsvc->find( "mu-" ) ;

  m_BzID        = BzInfo->particleID();
  m_BzbarID     = BzInfo->antiParticle()->particleID();
  m_BplusID     = BplInfo->particleID();
  m_BminusID    = BplInfo->antiParticle()->particleID();

  m_DzID        = DzInfo->particleID();
  m_DzbarID     = DzInfo->antiParticle()->particleID();
  m_DplusID     = DplInfo->particleID();
  m_DminusID    = DplInfo->antiParticle()->particleID();
  
  m_pimID       = pimInfo->particleID();
  m_pipID       = pipInfo->particleID();
  m_KpID        = KpInfo->particleID();
  m_KmID        = KmInfo->particleID();

  m_mumID       = munInfo->particleID();
  m_mupID       = mupInfo->particleID();

  Set = new RooArgSet(AMass, AMMass, nA, nBz, nBc, mu_mu_NNPID, K_KNNPID,"ArgSet");
  Set->add( RooArgSet ( BVertexChi2, DVertexChi2, BDIRA, K_IPchi2, Mu_Overlap, BTagMap ) );
  Set->add( RooArgSet ( Mu1_probChi2, Mu2_probChi2, Mu1_pErr, Mu2_pErr, Mu1_Ghost, Mu2_Ghost ) );
  rawdata = new RooDataSet("rawdata","rawdata", *Set);
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MyA2MuMuAnalysis::execute() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

  setFilterPassed(true);  // Mandatory. Set to true if event is accepted.

  // TagBz = NULL; TagBzbar = NULL; TagBplus = NULL; TagBminus = NULL;
  
  LHCb::Particle::Range BzParticles   = getIfExists<LHCb::Particle::Range>( m_B_zLocation ); 
  LHCb::Particle::Range BplsParticles = getIfExists<LHCb::Particle::Range>( m_Bpl_Location ); 
  LHCb::Particle::Range AParticles    = getIfExists<LHCb::Particle::Range>( m_A_Location );  
  
  LHCb::RecVertex::Range DataPvs = getIfExists<LHCb::RecVertex::Range>(m_DataPVLocation);

  LHCb::Particle::ConstVector Bz, Bzb, Bpl, Bm;
  
  NeutralBs = MyA2MuMuAnalysis::fetchParts ( m_BzID, BzParticles, Bz, Bzb );
  ChargedBs = MyA2MuMuAnalysis::fetchParts ( m_BplusID, BplsParticles, Bpl, Bm );
  
  StatusCode BzTag  = MyA2MuMuAnalysis::checkforBMeson( Bz, m_BzID, TagBz);
  StatusCode BzbTag = MyA2MuMuAnalysis::checkforBMeson( Bzb, m_BzbarID, TagBzbar);
  StatusCode BplTag = MyA2MuMuAnalysis::checkforBMeson( Bpl, m_BplusID, TagBplus);
  StatusCode BmTag  = MyA2MuMuAnalysis::checkforBMeson( Bm, m_BminusID, TagBminus);  

  //if ( TagBz )
  debug() << &TagBz << endmsg;
  //if ( TagBzbar )
  debug() << &TagBzbar << endmsg;
  //if ( TagBplus )
  debug() << &TagBplus << endmsg;
  //if ( TagBminus )
  debug() << &TagBminus << endmsg;

  double BTagMapping(0);

  if ( BzTag.isFailure() && BzbTag.isFailure() && BplTag.isFailure() && BmTag.isFailure() )
    BTagMapping += TMath::Power( 2.,0 ) ;

  if ( BzTag.isSuccess() )
    BTagMapping += TMath::Power( 2.,1 ) ;
  if ( BzbTag.isSuccess() )
    BTagMapping += TMath::Power( 2.,2 ) ;
  if ( BplTag.isSuccess() )
    BTagMapping += TMath::Power( 2.,3 ) ;
  if ( BmTag.isSuccess() )
    BTagMapping += TMath::Power( 2.,4 ) ;
  
  // bool BTag(false);
  
  // if ( (BzTag == StatusCode::SUCCESS && BzbTag == StatusCode::SUCCESS)  ||
  //      (BplTag == StatusCode::SUCCESS && BmTag == StatusCode::SUCCESS)  ||
  //      (BzTag == StatusCode::SUCCESS && BplTag == StatusCode::SUCCESS)  ||
  //      (BzTag == StatusCode::SUCCESS && BmTag == StatusCode::SUCCESS)   ||
  //      (BzbTag == StatusCode::SUCCESS && BplTag == StatusCode::SUCCESS) ||
  //      (BzbTag == StatusCode::SUCCESS && BmTag == StatusCode::SUCCESS)
  //      )
  // {
  //   debug() << "B Tag sucessful with" << endmsg;
  //   debug() << "B Candidate of type " <<  m_ppsvc->find( m_BzID )->name() << " " << BzTag.isSuccess() << endmsg;
  //   debug() << "B Candidate of type " <<  m_ppsvc->find( m_BzbarID )->name() << " " << BzbTag.isSuccess() << endmsg;
  //   debug() << "B Candidate of type " <<  m_ppsvc->find( m_BplusID )->name() << " " << BplTag.isSuccess() << endmsg;
  //   debug() << "B Candidate of type " <<  m_ppsvc->find( m_BminusID )->name() << " " << BmTag.isSuccess() << endmsg;

  //   BTag = true;
  // }

  BTagMap.setVal( BTagMapping ) ;  

  LHCb::Particle::ConstVector::const_iterator IterB;
  for (IterB = AParticles.begin(); IterB != AParticles.end(); IterB++)
  {
    // if ( BTag )
    {
      debug() << "Search for A->mu mu" << endmsg;
      StatusCode CheckOverlap = MyA2MuMuAnalysis::checkforMuonOverlap( (*IterB), BzTag, BzbTag, BplTag, BmTag );
      Mu_Overlap.setVal( (double) CheckOverlap.isSuccess() );

      Mu1_probChi2.setVal( (*IterB)->daughtersVector().at(0)->proto()->track()->probChi2() );
      Mu2_probChi2.setVal( (*IterB)->daughtersVector().at(1)->proto()->track()->probChi2() );
      
      double Mu1_pErrVal(0), Mu2_pErrVal(0);
      MyA2MuMuAnalysis::calculateMomError( *(*IterB)->daughtersVector().at(0), Mu1_pErrVal);
      MyA2MuMuAnalysis::calculateMomError( *(*IterB)->daughtersVector().at(1), Mu2_pErrVal);
      
      Mu1_pErr.setVal( Mu1_pErrVal );
      Mu2_pErr.setVal( Mu2_pErrVal );
      
      Mu1_Ghost.setVal( (*IterB)->daughtersVector().at(0)->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) );
      Mu2_Ghost.setVal( (*IterB)->daughtersVector().at(1)->proto()->info(LHCb::ProtoParticle::ProbNNghost,1e9) );
      
      nBz.setVal((double)NeutralBs);
      nBc.setVal((double)ChargedBs);
      nA.setVal((double)AParticles.size());
      
      AMass.setVal( (*IterB)->momentum().M() );
      AMMass.setVal( (*IterB)->measuredMass() );
      
      debug() << "A Candidate found" << endmsg;
      counter("Found_As")++;
      rawdata->add(*Set,1.0);
      
    }
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MyA2MuMuAnalysis::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;
  warning() << "Constructing Datasets" << endmsg;
  TFile Output(RootFile.c_str(),"recreate");
  Output.cd();
  rawdata->Write();
  Output.Close();

  return DaVinciAlgorithm::finalize();  // must be called after all other actions
}
//=============================================================================
std::size_t MyA2MuMuAnalysis::fetchParts (LHCb::ParticleID& type,
                                          LHCb::Particle::Range& alltracks,
                                          LHCb::Particle::ConstVector& parts1,
                                          LHCb::Particle::ConstVector& parts2)
{
  LHCb::ParticleID atype   = m_ppsvc->find( type )->antiParticle()->particleID();
  std::size_t listlenpart  = DaVinci::filter ( alltracks ,
                                               bind(&LHCb::Particle::particleID,_1) == type ,
                                               parts1 ) ;
  std::size_t listlenapart = DaVinci::filter ( alltracks ,
                                               bind(&LHCb::Particle::particleID,_1) == atype ,
                                               parts2 ) ;
  
  return listlenpart + listlenapart;
  
}
//=============================================================================
StatusCode MyA2MuMuAnalysis::checkforBMeson(LHCb::Particle::ConstVector& B,
                                            LHCb::ParticleID& BHypo,
                                            LHCb::Particle& BMeson)
{
  LHCb::Particle::ConstVector::const_iterator iB;
  const LHCb::Vertex *BPV = NULL;
  const LHCb::Particle *TagB = NULL;
  bool BTag(false);
  double K_PID(1e-9), Mu_PID(1e-9), BVtx(1e-9), DVtx(1e-9), BWerner(1e-9), KIPchi2(-1e9);
  
  for (iB = B.begin(); iB != B.end();++iB)
    {
      BPV = (const LHCb::Vertex*) bestPV((*iB));
      if ( BPV ) {
        TVector3 FlB0, Bp;
        FlB0.SetX( - BPV->position().x() + (*iB)->endVertex()->position().x() );
        FlB0.SetY( - BPV->position().y() + (*iB)->endVertex()->position().y() );
        FlB0.SetZ( - BPV->position().z() + (*iB)->endVertex()->position().z() );
        Bp.SetX( (*iB)->momentum().px() );
        Bp.SetY( (*iB)->momentum().py() );
        Bp.SetZ( (*iB)->momentum().pz() );
        
        double tmpBWerner = Bp.Angle(FlB0);
        if ( cos(tmpBWerner) > m_BDIRA ) 
        {
          double tmpBVtx = (*iB)->endVertex()->chi2PerDoF();
          if ( tmpBVtx < m_B_VertexCut )
          {
            LHCb::Particle::ConstVector PartDau = (*iB)->daughtersVector();
            const LHCb::Particle *DMeson = PartDau.at(0);
            const LHCb::Particle *Muon   = PartDau.at(1);
            double tmpMu_PID = Muon->proto()->info(LHCb::ProtoParticle::ProbNNmu, -1e9);
            if ( tmpMu_PID > m_Mu_PID_Cut)
            {
              if ( DMeson->momentum().M() > m_DMassWindow.at(0) || DMeson->momentum().M() < m_DMassWindow.at(1) )
              {
                double tmpDVtx = DMeson->endVertex()->chi2PerDoF();
                if ( tmpDVtx < m_D0_VertexCut )
                {
                  LHCb::Particle::ConstVector DDau = DMeson->daughtersVector();
                  const LHCb::Particle *Kaon = DDau.at(0);
                  double tmpK_PID = Kaon->proto()->info(LHCb::ProtoParticle::ProbNNk, -1e9);
                  if ( tmpK_PID > m_K_PID_Cut)
                  {
                    BPV = NULL;
                    BPV = (const LHCb::Vertex*) bestPV((Kaon));
                    if ( BPV ) {
                      double tmpKIP(-1e99), tmpKIPchi2(-1e99);
                      distanceCalculator()->distance( Kaon, BPV, tmpKIP, tmpKIPchi2);
                      if ( tmpKIPchi2 > m_KIPchi2 )
                      {
                        BVtx = tmpBVtx;
                        DVtx = tmpDVtx;
                        BWerner = tmpBWerner;
                        Mu_PID = tmpMu_PID;
                        K_PID = tmpK_PID;
                        KIPchi2 = tmpKIPchi2;
                        BTag = true;
                        if ( BHypo.hasDown() ) NeutralBs--;
                        if ( BHypo.hasUp () ) ChargedBs--;
                        TagB = ( (*iB) );
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  if ( BTag )
  {
    K_KNNPID.setVal( K_PID );
    mu_mu_NNPID.setVal( Mu_PID );
    BVertexChi2.setVal( BVtx );
    DVertexChi2.setVal( DVtx );
    BDIRA.setVal(BWerner);
    K_IPchi2.setVal(KIPchi2);

    BMeson = *TagB ;
    debug() << BMeson << endmsg;

    debug() << "B Candidate of type " <<  m_ppsvc->find( BHypo )->name() << " found" << endmsg;
    // debug() << "Deleting Tag B Particle" << endmsg;
    // delete TagB;
    return StatusCode::SUCCESS;
  }
  
  return StatusCode::FAILURE;
}
//=============================================================================
StatusCode MyA2MuMuAnalysis::checkforMuonOverlap( const LHCb::Particle *A,
                                                  StatusCode BzTag,
                                                  StatusCode BzbTag,
                                                  StatusCode BplTag,
                                                  StatusCode BmTag )
{
  const LHCb::Particle *BMuon = NULL;
  const LHCb::Particle *A_mu1 = NULL;
  const LHCb::Particle *A_mu2 = NULL;
  debug() << "Checking overlap for B0 " << BzTag.isSuccess() << endmsg;
  //Check overlap A B0 if tagged
  if ( BzTag.isSuccess() )
  {
    debug() << TagBz << endmsg;
    BMuon = TagBz.daughtersVector().at(1);
    A_mu1 = A->daughtersVector().at(0);
    A_mu2 = A->daughtersVector().at(1);

    if ( checkOverlap()->foundOverlap( BMuon, A_mu1 ) )
      return StatusCode::FAILURE;
    if ( checkOverlap()->foundOverlap( BMuon, A_mu2 ) )
      return StatusCode::FAILURE;
  }
  
  debug() << "Checking overlap for B~0 " << BzTag.isSuccess() << endmsg;
  //Check overlap A B~0 if tagged
  if ( BzbTag.isSuccess() )
  {
    debug() << TagBzbar << endmsg;
    BMuon = TagBzbar.daughtersVector().at(1);
    A_mu1 = A->daughtersVector().at(0);
    A_mu2 = A->daughtersVector().at(1);

    if ( checkOverlap()->foundOverlap( BMuon, A_mu1 ) )
      return StatusCode::FAILURE;
    if ( checkOverlap()->foundOverlap( BMuon, A_mu2 ) )
      return StatusCode::FAILURE;
  }

  debug() << "Checking overlap for B+ " << BzTag.isSuccess() << endmsg;
  //Check overlap A B+ if tagged
  if ( BplTag.isSuccess() )
  {
    debug() << TagBplus << endmsg;
    BMuon = TagBplus.daughtersVector().at(1);
    A_mu1 = A->daughtersVector().at(0);
    A_mu2 = A->daughtersVector().at(1);

    if ( checkOverlap()->foundOverlap( BMuon, A_mu1 ) )
      return StatusCode::FAILURE;
    if ( checkOverlap()->foundOverlap( BMuon, A_mu2 ) )
      return StatusCode::FAILURE;
  }

  debug() << "Checking overlap for B- " << BzTag.isSuccess() << endmsg;
  //Check overlap A B- if tagged
  if ( BmTag.isSuccess() )
  {
    debug() << TagBminus << endmsg;
    BMuon = TagBminus.daughtersVector().at(1);
    A_mu1 = A->daughtersVector().at(0);
    A_mu2 = A->daughtersVector().at(1);

    if ( checkOverlap()->foundOverlap( BMuon, A_mu1 ) )
      return StatusCode::FAILURE;
    if ( checkOverlap()->foundOverlap( BMuon, A_mu2 ) )
      return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}
//=============================================================================
StatusCode MyA2MuMuAnalysis::calculateMomError( const LHCb::Particle& part, double& Error)
{
  double qOverP = part.proto()->track()->firstState().qOverP();
  double u_qOverP = sqrt(part.proto()->track()->firstState().covariance()[4][4] );

  Error = u_qOverP / qOverP;
  
  return StatusCode::SUCCESS;
}
//=============================================================================
