// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/NTupler
// Class:      MiniFromPat
// 
/**\class MiniFromPat MiniFromPat.cc PhaseTwoAnalysis/NTupler/plugins/MiniFromPat.cc

Description: produces flat ntuples from PAT collections
   - storing gen, reco, and pf leptons with pT > 10 GeV and |eta| < 3
   - storing gen and reco jets with pT > 20 GeV and |eta| < 5

Implementation:
   - muon isolation comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_isolation
   - muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_identification
   - electron isolation might need to be refined
   - electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
      /!\ no ID is implemented for forward electrons as:
      - PFClusterProducer does not run on miniAOD
      - jurassic isolation needs tracks
   - PF jet ID comes from Run-2 https://github.com/cms-sw/cmssw/blob/CMSSW_9_1_1_patch1/PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
   - no JEC applied
   - b-tagging WPs come from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#B_tagging 
*/

//
// Original Author:  Elvire Bouvier
//         Created:  Tue, 20 Jun 2017 11:27:12 GMT
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"//
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MiniFromPat : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
  public:
    explicit MiniFromPat(const edm::ParameterSet&);
    ~MiniFromPat();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0,
      TRUE_PROMPT_ELECTRON,
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON};  


  private:
    virtual void beginJob() override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    void genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    void recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endJob() override;

    bool isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot); 
    bool isME0MuonSel(reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);
    bool isME0MuonSelNew(reco::Muon, double, double, double);

    // ----------member data ---------------------------
    edm::Service<TFileService> fs_;

    unsigned int pileup_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
    edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
    PFJetIDSelectionFunctor jetIDLoose_;
    PFJetIDSelectionFunctor jetIDTight_;
    edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
    edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
    edm::EDGetTokenT<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>> genVerticesToken_;
    edm::EDGetTokenT<std::vector<pat::Photon>> photonsToken_;
    const ME0Geometry* ME0Geometry_; 
    double mvaThres_[3];
    double deepThres_[3];

    TTree *t_event_, *t_genParts_, *t_vertices_, *t_genJets_, *t_genPhotons_, *t_looseElecs_, *t_tightElecs_, *t_looseMuons_, *t_tightMuons_, *t_puppiJets_, *t_puppiMET_, *t_loosePhotons_, *t_tightPhotons_;

    MiniEvent_t ev_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MiniFromPat::MiniFromPat(const edm::ParameterSet& iConfig):
  pileup_(iConfig.getParameter<unsigned int>("pileup")),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
  elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),  
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
  jetIDTight_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT), 
  metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
  genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
  genPartsToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
  genVerticesToken_(consumes<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>>(iConfig.getParameter<edm::InputTag>("genVertices"))),
  photonsToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons")))
{
  //now do what ever initialization is needed
  if (pileup_ == 0) {
    mvaThres_[0] = -0.694;
    mvaThres_[1] = 0.128;
    mvaThres_[2] = 0.822;
    deepThres_[0] = 0.131;
    deepThres_[1] = 0.432;
    deepThres_[2] = 0.741;
  } else if (pileup_ == 140) {
    mvaThres_[0] = -0.654;
    mvaThres_[1] = 0.214;
    mvaThres_[2] = 0.864;
    deepThres_[0] = 0.159;
    deepThres_[1] = 0.507;
    deepThres_[2] = 0.799;
  } else if (pileup_ == 200) {
    mvaThres_[0] = -0.642;
    mvaThres_[1] = 0.236;
    mvaThres_[2] = 0.878;
    deepThres_[0] = 0.170;
    deepThres_[1] = 0.527;
    deepThres_[2] = 0.821;
  } else {
    mvaThres_[0] = -1.;
    mvaThres_[1] = -1.;
    mvaThres_[2] = -1.;
    deepThres_[0] = 0.;
    deepThres_[1] = 0.;
    deepThres_[2] = 0.;
  }  

  usesResource("TFileService");

  t_event_      = fs_->make<TTree>("Event","Event");
  t_genParts_   = fs_->make<TTree>("Particle","Particle");
  t_genPhotons_ = fs_->make<TTree>("GenPhoton","GenPhoton"); 
  t_vertices_   = fs_->make<TTree>("Vertex","Vertex");
  t_genJets_    = fs_->make<TTree>("GenJet","GenJet");
  t_looseElecs_ = fs_->make<TTree>("ElectronLoose","ElectronLoose");
  t_tightElecs_ = fs_->make<TTree>("ElectronTight","ElectronTight");
  t_looseMuons_ = fs_->make<TTree>("MuonLoose","MuonLoose");
  t_tightMuons_ = fs_->make<TTree>("MuonTight","MuonTight");
  t_puppiJets_  = fs_->make<TTree>("JetPUPPI","JetPUPPI");
  t_puppiMET_   = fs_->make<TTree>("PuppiMissingET","PuppiMissingET");
  t_loosePhotons_ = fs_->make<TTree>("PhotonLoose","PhotonLoose");
  t_tightPhotons_ = fs_->make<TTree>("PhotonTight","PhotonTight");
  createMiniEventTree(t_event_, t_genParts_, t_vertices_, t_genJets_, t_genPhotons_, t_looseElecs_, t_tightElecs_, t_looseMuons_, t_tightMuons_, t_puppiJets_, t_puppiMET_, t_loosePhotons_, t_tightPhotons_, ev_);

}


MiniFromPat::~MiniFromPat()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method to fill gen level pat -------------
  void
MiniFromPat::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<std::vector<pat::PackedGenParticle>> genParts;
  iEvent.getByToken(genPartsToken_, genParts);

  Handle<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>> genVertices;
  iEvent.getByToken(genVerticesToken_, genVertices);

  Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);

  // Jets
  std::vector<size_t> jGenJets;
  ev_.ngj = 0;
  for (size_t i = 0; i < genJets->size(); i++) {
    if (genJets->at(i).pt() < 20.) continue;
    if (fabs(genJets->at(i).eta()) > 5) continue;

    bool overlaps = false;
    for (size_t j = 0; j < genParts->size(); j++) {
      if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
      if (fabs(genJets->at(i).pt()-genParts->at(j).pt()) < 0.01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    jGenJets.push_back(i);

    ev_.gj_pt[ev_.ngj]   = genJets->at(i).pt();
    ev_.gj_phi[ev_.ngj]  = genJets->at(i).phi();
    ev_.gj_eta[ev_.ngj]  = genJets->at(i).eta();
    ev_.gj_mass[ev_.ngj] = genJets->at(i).mass();
    ev_.ngj++;
  }

  // Leptons
  ev_.ngl = 0;
  for (size_t i = 0; i < genParts->size(); i++) {
    if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    if (genParts->at(i).pt() < 10.) continue;
    if (fabs(genParts->at(i).eta()) > 3.) continue;
    double genIso = 0.;
    for (size_t j = 0; j < jGenJets.size(); j++) {
      if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue; 
      std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
      for (size_t k = 0; k < jconst.size(); k++) {
        double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
        if (deltaR < 0.01 || deltaR > 0.4) continue;
        genIso = genIso + jconst[k]->pt();
      }
    }
    genIso = genIso / genParts->at(i).pt();
    ev_.gl_pid[ev_.ngl]    = genParts->at(i).pdgId();
    ev_.gl_ch[ev_.ngl]     = genParts->at(i).charge();
    ev_.gl_st[ev_.ngl]     = genParts->at(i).status();
    ev_.gl_p[ev_.ngl]      = genParts->at(i).p();
    ev_.gl_px[ev_.ngl]     = genParts->at(i).px();
    ev_.gl_py[ev_.ngl]     = genParts->at(i).py();
    ev_.gl_pz[ev_.ngl]     = genParts->at(i).pz();
    ev_.gl_nrj[ev_.ngl]    = genParts->at(i).energy();
    ev_.gl_pt[ev_.ngl]     = genParts->at(i).pt();
    ev_.gl_phi[ev_.ngl]    = genParts->at(i).phi();
    ev_.gl_eta[ev_.ngl]    = genParts->at(i).eta();
    ev_.gl_mass[ev_.ngl]   = genParts->at(i).mass();
    ev_.gl_relIso[ev_.ngl] = genIso; 
    ev_.ngl++;
  }

  // Vertex
  float vtxZ = genVertices->z();
  std::cout << "vertex z is s = " << genVertices->z() << std::endl;

  // Photons
  ev_.ngp = 0;
  for (size_t i = 0; i < genParts->size(); i++) {
    if (abs(genParts->at(i).pdgId()) != 22) continue;
    if (genParts->at(i).pt() < 10.) continue;
    if (fabs(genParts->at(i).eta()) > 3.) continue;

    ev_.gp_st[ev_.ngp]     = genParts->at(i).status();
    ev_.gp_p[ev_.ngp]      = genParts->at(i).p();
    ev_.gp_px[ev_.ngp]     = genParts->at(i).px();
    ev_.gp_py[ev_.ngp]     = genParts->at(i).py();
    ev_.gp_pz[ev_.ngp]     = genParts->at(i).pz();
    ev_.gp_nrj[ev_.ngp]    = genParts->at(i).energy();
    ev_.gp_pt[ev_.ngp]     = genParts->at(i).pt();
    ev_.gp_phi[ev_.ngp]    = genParts->at(i).phi();
    ev_.gp_eta[ev_.ngp]    = genParts->at(i).eta();
    ev_.gp_vtxz[ev_.ngp]   = vtxZ;
    ev_.ngp++;
  }
}

// ------------ method to fill reco level pat -------------
  void
MiniFromPat::recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);

  Handle<std::vector<pat::Electron>> elecs;
  iEvent.getByToken(elecsToken_, elecs);
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convToken_, conversions);
  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();  

  Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);

  Handle<std::vector<pat::MET>> mets;
  iEvent.getByToken(metsToken_, mets);

  Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);

  Handle<std::vector<pat::Photon>> photons;
  iEvent.getByToken(photonsToken_, photons);
   
  // Vertices
  int prVtx = -1;
  ev_.nvtx = 0;
  for (size_t i = 0; i < vertices->size(); i++) {
    if (vertices->at(i).isFake()) continue;
    if (vertices->at(i).ndof() <= 4) continue;
    if (prVtx < 0) prVtx = i;
    ev_.v_pt2[ev_.nvtx] = vertices->at(i).p4().pt();
    ev_.nvtx++;
  }
  if (prVtx < 0) return;

  // Muons
  ev_.nlm = 0;
  ev_.ntm = 0;

  for (size_t i = 0; i < muons->size(); i++) {
    if (muons->at(i).pt() < 2.) continue;
    if (fabs(muons->at(i).eta()) > 2.8) continue;

    // Loose ID
    double dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.056);
    double dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0096);    
    bool isLoose = (fabs(muons->at(i).eta()) < 2.4 && muon::isLooseMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut));

    // Medium ID -- needs to be updated
    bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    if (muons->at(i).innerTrack().isNonnull()){
    	ipxy = std::abs(muons->at(i).muonBestTrack()->dxy(vertices->at(prVtx).position())) < 0.2;
    	ipz = std::abs(muons->at(i).muonBestTrack()->dz(vertices->at(prVtx).position())) < 0.5;
    	validPxlHit = muons->at(i).innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    	highPurity = muons->at(i).innerTrack()->quality(reco::Track::highPurity);
    }    
    // bool isMedium = (fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);

    // Tight ID
    dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.032);
    dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0041);
    bool isTight = (fabs(muons->at(i).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(i),vertices->at(prVtx))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.048, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);

    if (!isLoose) continue;

    ev_.lm_ch[ev_.nlm]     = muons->at(i).charge();
    ev_.lm_pt[ev_.nlm]     = muons->at(i).pt();
    ev_.lm_phi[ev_.nlm]    = muons->at(i).phi();
    ev_.lm_eta[ev_.nlm]    = muons->at(i).eta();
    ev_.lm_mass[ev_.nlm]   = muons->at(i).mass();
    ev_.lm_relIso[ev_.nlm] = (muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt();
    ev_.lm_g[ev_.nlm] = -1;
    for (int ig = 0; ig < ev_.ngl; ig++) {
      if (abs(ev_.gl_pid[ig]) != 13) continue;
      if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.lm_eta[ev_.nlm],ev_.lm_phi[ev_.nlm]) > 0.4) continue;
      ev_.lm_g[ev_.nlm]    = ig;
    }
    ev_.nlm++;

    if (!isTight) continue;

    ev_.tm_ch[ev_.ntm]     = muons->at(i).charge();
    ev_.tm_pt[ev_.ntm]     = muons->at(i).pt();
    ev_.tm_phi[ev_.ntm]    = muons->at(i).phi();
    ev_.tm_eta[ev_.ntm]    = muons->at(i).eta();
    ev_.tm_mass[ev_.ntm]   = muons->at(i).mass();
    ev_.tm_relIso[ev_.ntm] = (muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt();
    ev_.tm_g[ev_.ntm] = -1;
    for (int ig = 0; ig < ev_.ngl; ig++) {
      if (abs(ev_.gl_pid[ig]) != 13) continue;
      if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.tm_eta[ev_.ntm],ev_.tm_phi[ev_.ntm]) > 0.4) continue;
      ev_.tm_g[ev_.ntm]    = ig;
    }
    ev_.ntm++;
  }

  // Electrons

  ev_.nle = 0;
  ev_.nte = 0;

  for (size_t i = 0; i < elecs->size(); i++) {
    if (elecs->at(i).pt() < 10.) continue;
    if (fabs(elecs->at(i).eta()) > 3.) continue;

    bool isLoose = isLooseElec(elecs->at(i),conversions,beamspot);    
    // bool isMedium = isMediumElec(elecs->at(i),conversions,beamspot);    
    bool isTight = isTightElec(elecs->at(i),conversions,beamspot);    

    if (!isLoose) continue;

    ev_.le_ch[ev_.nle]     = elecs->at(i).charge();
    ev_.le_pt[ev_.nle]     = elecs->at(i).pt();
    ev_.le_phi[ev_.nle]    = elecs->at(i).phi();
    ev_.le_eta[ev_.nle]    = elecs->at(i).eta();
    ev_.le_mass[ev_.nle]   = elecs->at(i).mass();
    ev_.le_relIso[ev_.nle] = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();
    ev_.le_g[ev_.nle] = -1;
    for (int ig = 0; ig < ev_.ngl; ig++) {
      if (abs(ev_.gl_pid[ig]) != 11) continue;
      if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.le_eta[ev_.nle],ev_.le_phi[ev_.nle]) > 0.4) continue;
      ev_.le_g[ev_.nle]    = ig;
    }
    ev_.nle++;

    if (!isTight) continue;

    ev_.te_ch[ev_.nte]     = elecs->at(i).charge();
    ev_.te_pt[ev_.nte]     = elecs->at(i).pt();
    ev_.te_phi[ev_.nte]    = elecs->at(i).phi();
    ev_.te_eta[ev_.nte]    = elecs->at(i).eta();
    ev_.te_mass[ev_.nte]   = elecs->at(i).mass();
    ev_.te_relIso[ev_.nte] = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();
    ev_.te_g[ev_.nte] = -1;
    for (int ig = 0; ig < ev_.ngl; ig++) {
      if (abs(ev_.gl_pid[ig]) != 11) continue;
      if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.te_eta[ev_.nte],ev_.te_phi[ev_.nte]) > 0.4) continue;
      ev_.te_g[ev_.nte]    = ig;
    }
    ev_.nte++;
  }

  // Jets
  ev_.nj = 0;
  for (size_t i =0; i < jets->size(); i++) {
    if (jets->at(i).pt() < 20.) continue;
    if (fabs(jets->at(i).eta()) > 5) continue;

    //    std::cout << "SCZ DEBUG JET " << i << " pt=" << jets->at(i).pt() << " eta=" << jets->at(i).eta() << std::endl;

    float sumCandPt = 0.;
    float sumCandPtSq = 0;
    float rmscand_n = 0.;
    float rmscand_d = 0.;
    float M11 = 0.;
    float M12 = 0.;
    float M21 = 0.;
    float M22 = 0.;
    float tot_wt = 0.;

    ev_.j_chargedSumConst[ev_.nj] = 0.;
    ev_.j_neutralSumConst[ev_.nj] = 0.;
    ev_.j_hfemSumConst[ev_.nj] = 0.;
    ev_.j_hfhadSumConst[ev_.nj] = 0.;
    ev_.j_chargedNConst[ev_.nj] = 0;
    ev_.j_neutralNConst[ev_.nj] = 0;
    ev_.j_hfemNConst[ev_.nj] = 0;
    ev_.j_hfhadNConst[ev_.nj]= 0;
    ev_.j_eSumConst[ev_.nj] = 0.;
    ev_.j_eNConst[ev_.nj] = 0;
    ev_.j_muSumConst[ev_.nj] = 0.;
    ev_.j_muNConst[ev_.nj] = 0;
    ev_.j_photonSumConst[ev_.nj] = 0.;
    ev_.j_photonNConst[ev_.nj] = 0;

    for ( unsigned k = 0; k < jets->at(i).numberOfSourceCandidatePtrs(); ++k ) {
      reco::CandidatePtr pfJetConstituent = jets->at(i).sourceCandidatePtr(k);
                    
      const reco::Candidate* kcand = pfJetConstituent.get();
      const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>( kcand );
      if ( !lPack ) throw cms::Exception( "NoPackedConstituent" ) << " For jet " << i << " failed to get constituent " << k << std::endl;
      float candPt = kcand->pt();
      float candDr   = reco::deltaR(*kcand,jets->at(i));

      sumCandPt += candPt;
      sumCandPtSq += candPt*candPt;

      float wt = pow(kcand->pt(),2);
      tot_wt += wt;
          
      //! RMSCand
      float deta = kcand->eta()- jets->at(i).eta();
      float dphi = deltaPhi(kcand->phi(), jets->at(i).phi());
      float dr = deltaR(kcand->eta(), kcand->phi(), jets->at(i).eta(),  jets->at(i).phi());
      rmscand_n += wt*dr*dr;
      rmscand_d += wt; 
          
      M11 += wt*deta*deta;
      M22 += wt*dphi*dphi;
      M12 += wt*deta*dphi;
      M21 += wt*deta*dphi;

      //      std::cout << "  SCZ DEBUG CONSTITUENT" << k << " " << candPt << " "  << candDr << std::endl;

      int pfid2     = kcand->pdgId();

      switch(std::abs(pfid2)){
      case 211:  //PFCandidate::h charged hadron
      case 321:
      case 2212:
	ev_.j_chargedSumConst[ev_.nj] += kcand->pt();
	ev_.j_chargedNConst[ev_.nj]++;
	break;
      case 11:  //PFCandidate::e // electron
	ev_.j_eSumConst[ev_.nj] += kcand->pt();
	ev_.j_eNConst[ev_.nj]++;
	break;
      case 13: //PFCandidate::mu // muon
	ev_.j_muSumConst[ev_.nj] += kcand->pt();
	ev_.j_muNConst[ev_.nj]++;
	break;
      case 22: //PFCandidate:gamma // gamma
	ev_.j_photonSumConst[ev_.nj] += kcand->pt();
	ev_.j_photonNConst[ev_.nj]++;
	break;
      case 130:  //PFCandidate::h0 //Neutral hadron
      case 2112:
	ev_.j_neutralSumConst[ev_.nj] += kcand->pt();
	ev_.j_neutralNConst[ev_.nj]++;
	break;
      case 1:  //PFCandidate::h_HF //hadron in HF
	ev_.j_hfhadSumConst[ev_.nj] += kcand->pt(); 
	ev_.j_hfhadNConst[ev_.nj]++;
	break;
      case 2:  //PFCandidate::egamma_HF //electromagnetic in HF
	ev_.j_hfemSumConst[ev_.nj] += kcand->pt(); 
	ev_.j_hfemNConst[ev_.nj]++;
	break;
      default:
	std::cout << "   SCZ DEBUG UNEXPECTED PF ID " << pfid2 << std::endl;
	break;
      }

    }

    ev_.j_RMSCand[ev_.nj] = (rmscand_d > 0 ) ? sqrt ( rmscand_n / rmscand_d ) : -999;
      
    M12 = -1.*M12;
    M21 = -1.*M21;
      
    //! eign values
    float trace = M11 + M22;
    float detrm = (M11*M22) - (M12*M21);
      
    float lam1 = trace/2. + sqrt( pow(trace,2)/4. - detrm );
    float lam2 = trace/2. - sqrt( pow(trace,2)/4. - detrm );
      
    ev_.j_Axis1[ev_.nj] = (tot_wt > 0 && lam1 >= 0) ? sqrt( lam1 / tot_wt ) : -999;
    ev_.j_Axis2[ev_.nj] = (tot_wt > 0 && lam2 >=0 ) ? sqrt( lam2 / tot_wt ) : -999;
  
    ev_.j_Sigma[ev_.nj] = (tot_wt > 0 ) ? sqrt( pow(ev_.j_Axis1[ev_.nj],2) + pow(ev_.j_Axis2[ev_.nj],2) ) : -999;
      
    ev_.j_ptD[ev_.nj] = (sumCandPt > 0.) ? std::sqrt(sumCandPtSq)/sumCandPt : -999.;
    
    
    //    std::cout << " SCZ DEBUG JET " << i << " ptD=" << ev_.j_ptD[ev_.nj] << " RMS=" << ev_.j_RMSCand[ev_.nj] << " Axis1=" << ev_.j_Axis1[ev_.nj] << " Axis2=" << ev_.j_Axis2[ev_.nj] << std::endl;
    //    std::cout << "   SCZ DEBUG chargedSumConst=" << ev_.j_chargedSumConst[ev_.nj] << " neutralSumConst=" << ev_.j_neutralSumConst[ev_.nj]
    //	      << " hfemSumConst=" << ev_.j_hfemSumConst[ev_.nj] << " hfhadSumConst=" << ev_.j_hfhadSumConst[ev_.nj] << std::endl;
    //    std::cout << "   SCZ DEBUG eSumConst=" << ev_.j_eSumConst[ev_.nj] << " muSumConst=" << ev_.j_muSumConst[ev_.nj] << " photonSumConst=" << ev_.j_photonSumConst[ev_.nj] << std::endl;
    //    std::cout << "   SCZ DEBUG chargedNConst=" << ev_.j_chargedNConst[ev_.nj] << " neutralNConst=" << ev_.j_neutralNConst[ev_.nj]
    //              << " hfemNConst=" << ev_.j_hfemNConst[ev_.nj] << " hfhadNConst=" << ev_.j_hfhadNConst[ev_.nj] << std::endl;
    //    std::cout << "   SCZ DEBUG eNConst=" << ev_.j_eNConst[ev_.nj] << " muNConst=" << ev_.j_muNConst[ev_.nj] << " photonNConst=" << ev_.j_photonNConst[ev_.nj] << std::endl;


    bool overlaps = false;
    for (size_t j = 0; j < elecs->size(); j++) {
      if (fabs(jets->at(i).pt()-elecs->at(j).pt()) < 0.01*elecs->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(elecs->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) continue;
    for (size_t j = 0; j < muons->size(); j++) {
      if (fabs(jets->at(i).pt()-muons->at(j).pt()) < 0.01*muons->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(muons->at(j).p4(),jets->at(i).p4()) < 0.01) {
        overlaps = true;
        break;
      }
    }

    //    std::cout << " SCZ DEBUG OVERLAPS=" << overlaps << std::endl;

    if (overlaps) continue;

    pat::strbitset retLoose = jetIDLoose_.getBitTemplate();
    retLoose.set(false);
    bool isLoose = jetIDLoose_(jets->at(i), retLoose);
    pat::strbitset retTight = jetIDTight_.getBitTemplate();
    retTight.set(false);
    bool isTight = jetIDTight_(jets->at(i), retTight);

    double mvav2   = jets->at(i).bDiscriminator("pfCombinedMVAV2BJetTags"); 
    bool isLooseMVAv2  = mvav2 > mvaThres_[0];
    bool isMediumMVAv2 = mvav2 > mvaThres_[1];
    bool isTightMVAv2  = mvav2 > mvaThres_[2];
    double deepcsv = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
                            jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
    bool isLooseDeepCSV  = deepcsv > deepThres_[0];
    bool isMediumDeepCSV = deepcsv > deepThres_[1];
    bool isTightDeepCSV  = deepcsv > deepThres_[2];

    ev_.j_id[ev_.nj]      = (isTight | (isLoose<<1));
    ev_.j_pt[ev_.nj]      = jets->at(i).pt();
    ev_.j_phi[ev_.nj]     = jets->at(i).phi();
    ev_.j_eta[ev_.nj]     = jets->at(i).eta();
    ev_.j_mass[ev_.nj]    = jets->at(i).mass();
    ev_.j_mvav2[ev_.nj]   = (isTightMVAv2 | (isMediumMVAv2<<1) | (isLooseMVAv2<<2)); 
    ev_.j_deepcsv[ev_.nj] = (isTightDeepCSV | (isMediumDeepCSV<<1) | (isLooseDeepCSV<<2));
    ev_.j_flav[ev_.nj]    = jets->at(i).partonFlavour();
    ev_.j_hadflav[ev_.nj] = jets->at(i).hadronFlavour();
    ev_.j_pid[ev_.nj]     = (jets->at(i).genParton() ? jets->at(i).genParton()->pdgId() : 0);
    ev_.j_g[ev_.nj] = -1;
    for (int ig = 0; ig < ev_.ngj; ig++) {
      if (reco::deltaR(ev_.gj_eta[ig],ev_.gj_phi[ig],ev_.j_eta[ev_.nj],ev_.j_phi[ev_.nj]) > 0.4) continue;
      ev_.j_g[ev_.nj]     = ig;
      break;
    }	

    std::cout << " SCZ DEBUG isLoose=" << isLoose << " isTight=" << isTight << " flav=" << ev_.j_flav[ev_.nj] << " j_pid=" << ev_.j_pid[ev_.nj] << std::endl;

    ev_.nj++;

  }
  
  // MET
  ev_.nmet = 0;
  if (mets->size() > 0) {
    ev_.met_pt[ev_.nmet]  = mets->at(0).pt();
    ev_.met_eta[ev_.nmet] = mets->at(0).eta();
    ev_.met_phi[ev_.nmet] = mets->at(0).phi();
    ev_.nmet++;
  }

  // Photons

  ev_.nlp = 0;
  ev_.ntp = 0;

  for (size_t i = 0; i < photons->size(); i++) {
    if (photons->at(i).pt() < 10.) continue;
    if (fabs(photons->at(i).eta()) > 3.) continue;

    float mvaValue = photons->at(i).userFloat("mvaValue");
    bool isEB = photons->at(i).isEB();
     
    bool isLoose = 0;
    bool isTight = 0;

    if( isEB )
      {
	 isLoose = (mvaValue > 0.00);
	 isTight = (mvaValue > 0.56);
      }     
     else
      {
	 isLoose = (mvaValue > 0.20);
	 isTight = (mvaValue > 0.68);
      }          

    if (!isLoose) continue;

    ev_.lp_pt[ev_.nlp]     = photons->at(i).pt();
    ev_.lp_phi[ev_.nlp]    = photons->at(i).phi();
    ev_.lp_eta[ev_.nlp]    = photons->at(i).eta();
    ev_.lp_nrj[ev_.nlp]    = photons->at(i).energy();
    ev_.lp_g[ev_.nlp] = -1;
    // add multicluster quantities too
    ev_.lp_pt_multi[ev_.nlp]     = photons->at(i).superCluster()->seed()->energy() / cosh(photons->at(i).superCluster()->seed()->eta());
    ev_.lp_phi_multi[ev_.nlp]    = photons->at(i).superCluster()->seed()->phi();
    ev_.lp_eta_multi[ev_.nlp]    = photons->at(i).superCluster()->seed()->eta();
    ev_.lp_nrj_multi[ev_.nlp]    = photons->at(i).superCluster()->seed()->energy();
    for (int ig = 0; ig < ev_.ngp; ig++) {
      if (reco::deltaR(ev_.gp_eta[ig],ev_.gp_phi[ig],ev_.lp_eta[ev_.nlp],ev_.lp_phi[ev_.nlp]) > 0.4) continue;
      ev_.lp_g[ev_.nlp]    = ig;
    }
    ev_.nlp++;

    if (!isTight) continue;

    ev_.tp_pt[ev_.ntp]     = photons->at(i).pt();
    ev_.tp_phi[ev_.ntp]    = photons->at(i).phi();
    ev_.tp_eta[ev_.ntp]    = photons->at(i).eta();
    ev_.tp_nrj[ev_.ntp]    = photons->at(i).energy();
    ev_.tp_g[ev_.ntp] = -1;
    // add multicluster quantities too
    ev_.tp_pt_multi[ev_.ntp]     = photons->at(i).superCluster()->seed()->energy() / cosh(photons->at(i).superCluster()->seed()->eta());
    ev_.tp_phi_multi[ev_.ntp]    = photons->at(i).superCluster()->seed()->phi();
    ev_.tp_eta_multi[ev_.ntp]    = photons->at(i).superCluster()->seed()->eta();
    ev_.tp_nrj_multi[ev_.ntp]    = photons->at(i).superCluster()->seed()->energy();
    for (int ig = 0; ig < ev_.ngp; ig++) {
      if (reco::deltaR(ev_.gp_eta[ig],ev_.gp_phi[ig],ev_.tp_eta[ev_.ntp],ev_.tp_phi[ev_.ntp]) > 0.4) continue;
      ev_.tp_g[ev_.ntp]    = ig;
    }
    ev_.ntp++;
  }
   
}

// ------------ method called for each event  ------------
  void
MiniFromPat::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //analyze the event
  if(!iEvent.isRealData()) genAnalysis(iEvent, iSetup);
  recoAnalysis(iEvent, iSetup);
  
  //save event if at least one lepton at gen or reco level
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  t_event_->Fill();
  t_genParts_->Fill();
  t_genPhotons_->Fill();
  t_vertices_->Fill();
  t_genJets_->Fill();
  t_looseElecs_->Fill();
  t_tightElecs_->Fill();
  t_looseMuons_->Fill();
  t_tightMuons_->Fill();
  t_puppiJets_->Fill();
  t_puppiMET_->Fill();
  t_loosePhotons_->Fill();
  t_tightPhotons_->Fill();

}


// ------------ method check that an e passes loose ID ----------------------------------
  bool
MiniFromPat::isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.02992) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.004119) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.05176) return false;
  if (patEl.hcalOverEcal() > 6.741) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 2.5) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 73.76) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes medium ID ----------------------------------
  bool
MiniFromPat::isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01609) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001766) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.03130) return false;
  if (patEl.hcalOverEcal() > 7.371) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.325) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 22.6) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method check that an e passes tight ID ----------------------------------
  bool
MiniFromPat::isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot) 
{
  if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
  if (patEl.full5x5_sigmaIetaIeta() > 0.01614) return false;
  if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001322) return false;
  if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.06129) return false;
  if (patEl.hcalOverEcal() > 4.492) return false;
  if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.255) return false;
  double Ooemoop = 999.;
  if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
  else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
  else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
  if (Ooemoop > 18.26) return false;
  if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
  return true;
}

// ------------ method to improve ME0 muon ID ----------------
  bool 
MiniFromPat::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaX = 999;
    double deltaY = 999;
    double pullX = 999;
    double pullY = 999;
    double deltaPhi = 999;

    bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for(std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber){

      for (std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment){

        if (chamber->detector() == 5){

          deltaX   = std::abs(chamber->x - segment->x);
          deltaY   = std::abs(chamber->y - segment->y);
          pullX    = std::abs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
          pullY    = std::abs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
          deltaPhi = std::abs(atan(chamber->dXdZ) - atan(segment->dXdZ));

        }
      }
    }

    if ((pullX < pullXCut) || (deltaX < dXCut)) X_MatchFound = true;
    if ((pullY < pullYCut) || (deltaY < dYCut)) Y_MatchFound = true;
    if (deltaPhi < dPhi) Dir_MatchFound = true;

    result = X_MatchFound && Y_MatchFound && Dir_MatchFound;

  }

  return result;

}

bool 
MiniFromPat::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut)
{

  bool result = false;
  bool isME0 = muon.isME0Muon();

  if(isME0){

    double deltaEta = 999;
    double deltaPhi = 999;
    double deltaPhiBend = 999;

    const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
    for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

      if (chamber->detector() == 5){

        for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

          LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
          LocalPoint seg_loc_coord(segment->x, segment->y, 0);
          LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
          LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);

          const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);

          GlobalPoint trk_glb_coord = me0chamber->toGlobal(trk_loc_coord);
          GlobalPoint seg_glb_coord = me0chamber->toGlobal(seg_loc_coord);

          //double segDPhi = segment->me0SegmentRef->deltaPhi();
          // need to check if this works
          double segDPhi = me0chamber->computeDeltaPhi(seg_loc_coord, seg_loc_vec);
          double trackDPhi = me0chamber->computeDeltaPhi(trk_loc_coord, trk_loc_vec);

          deltaEta = std::abs(trk_glb_coord.eta() - seg_glb_coord.eta() );
          deltaPhi = std::abs(trk_glb_coord.phi() - seg_glb_coord.phi() );
          deltaPhiBend = std::abs(segDPhi - trackDPhi);

          if (deltaEta < dEtaCut && deltaPhi < dPhiCut && deltaPhiBend < dPhiBendCut) result = true;

        }
      }
    }

  }

  return result;

}

// ------------ method called once each job just before starting event loop  ------------
  void 
MiniFromPat::beginJob()
{
}

// ------------ method called once each run ----------------
void
MiniFromPat::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  ME0Geometry_ =( &*hGeom);
}

// ------------ method called when ending the processing of a run  ------------
  void
MiniFromPat::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
MiniFromPat::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniFromPat::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniFromPat);
