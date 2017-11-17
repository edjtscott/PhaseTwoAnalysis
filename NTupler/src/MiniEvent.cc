#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

void createMiniEventTree(TTree *t_event_, TTree *t_genParts_, TTree *t_vertices_, TTree *t_genJets_, TTree *t_genPhotons_, TTree *t_looseElecs_, TTree *t_tightElecs_, TTree *t_looseMuons_, TTree *t_tightMuons_, TTree *t_puppiJets_, TTree *t_puppiMET_, TTree *t_loosePhotons_, TTree *t_tightPhotons_, MiniEvent_t &ev)
{
  //event header
  t_event_->Branch("Run",               &ev.run,        "Run/I");
  t_event_->Branch("Event",             &ev.event,      "Event/I");
  t_event_->Branch("Lumi",              &ev.lumi,       "Lumi/I");

  //gen level event
  t_genParts_->Branch("Particle_size",  &ev.ngl,        "Particle_size/I");
  t_genParts_->Branch("PID",            ev.gl_pid,      "PID[Particle_size]/I");
  t_genParts_->Branch("Charge",         ev.gl_ch,       "Charge[Particle_size]/I");
  t_genParts_->Branch("Status",         ev.gl_st,       "Status[Particle_size]/I");
  t_genParts_->Branch("P",              ev.gl_p,        "P[Particle_size]/F");
  t_genParts_->Branch("Px",             ev.gl_px,       "Px[Particle_size]/F");
  t_genParts_->Branch("Py",             ev.gl_py,       "Py[Particle_size]/F");
  t_genParts_->Branch("Pz",             ev.gl_pz,       "Pz[Particle_size]/F");
  t_genParts_->Branch("E",              ev.gl_nrj,      "E[Particle_size]/F");
  t_genParts_->Branch("PT",             ev.gl_pt,       "PT[Particle_size]/F");
  t_genParts_->Branch("Eta",            ev.gl_eta,      "Eta[Particle_size]/F");
  t_genParts_->Branch("Phi",            ev.gl_phi,      "Phi[Particle_size]/F");
  t_genParts_->Branch("Mass",           ev.gl_mass,     "Mass[Particle_size]/F");
  t_genParts_->Branch("IsolationVar",   ev.gl_relIso,   "IsolationVar/F");

  t_genJets_->Branch("GenJet_size",     &ev.ngj,        "GenJet_size/I");
  t_genJets_->Branch("PT",              ev.gj_pt,       "PT[GenJet_size]/F");
  t_genJets_->Branch("Eta",             ev.gj_eta,      "Eta[GenJet_size]/F");
  t_genJets_->Branch("Phi",             ev.gj_phi,      "Phi[GenJet_size]/F");
  t_genJets_->Branch("Mass",            ev.gj_mass,     "Mass[GenJet_size]/F");

  t_genPhotons_->Branch("GenPhoton_size",  &ev.ngp,        "GenPhoton_size/I");
  t_genPhotons_->Branch("Status",          ev.gp_st,       "Status[GenPhoton_size]/I");
  t_genPhotons_->Branch("P",               ev.gp_p,        "P[GenPhoton_size]/F");
  t_genPhotons_->Branch("Px",              ev.gp_px,       "Px[GenPhoton_size]/F");
  t_genPhotons_->Branch("Py",              ev.gp_py,       "Py[GenPhoton_size]/F");
  t_genPhotons_->Branch("Pz",              ev.gp_pz,       "Pz[GenPhoton_size]/F");
  t_genPhotons_->Branch("E",               ev.gp_nrj,      "E[GenPhoton_size]/F");
  t_genPhotons_->Branch("PT",              ev.gp_pt,       "PT[GenPhoton_size]/F");
  t_genPhotons_->Branch("Eta",             ev.gp_eta,      "Eta[GenPhoton_size]/F");
  t_genPhotons_->Branch("Phi",             ev.gp_phi,      "Phi[GenPhoton_size]/F");
  t_genPhotons_->Branch("Vtxz",            ev.gp_vtxz,     "Vtxz[GenPhoton_size]/F");
   
  //reco level event
  t_vertices_->Branch("Vertex_size",    &ev.nvtx,       "Vertex_size/I");
  t_vertices_->Branch("SumPT2",         &ev.v_pt2,      "SumPT2[Vertex_size]/F");

  t_looseElecs_->Branch("ElectronLoose_size", &ev.nle,  "ElectronLoose_size/I");
  t_looseElecs_->Branch("Charge",       ev.le_ch,       "Charge[ElectronLoose_size]/I");
  t_looseElecs_->Branch("Particle",     ev.le_g,        "Particle[ElectronLoose_size]/I");
  t_looseElecs_->Branch("PT",           ev.le_pt,       "PT[ElectronLoose_size]/F");
  t_looseElecs_->Branch("Eta",          ev.le_eta,      "Eta[ElectronLoose_size]/F");
  t_looseElecs_->Branch("Phi",          ev.le_phi,      "Phi[ElectronLoose_size]/F");
  t_looseElecs_->Branch("Mass",         ev.le_mass,     "Mass[ElectronLoose_size]/F");
  t_looseElecs_->Branch("IsolationVar", ev.le_relIso,   "IsolationVar[ElectronLoose_size]/F");

  t_tightElecs_->Branch("ElectronTight_size", &ev.nte,  "ElectronTight_size/I");
  t_tightElecs_->Branch("Charge",       ev.te_ch,       "Charge[ElectronTight_size]/I");
  t_tightElecs_->Branch("Particle",     ev.te_g,        "Particle[ElectronTight_size]/I");
  t_tightElecs_->Branch("PT",           ev.te_pt,       "PT[ElectronTight_size]/F");
  t_tightElecs_->Branch("Eta",          ev.te_eta,      "Eta[ElectronTight_size]/F");
  t_tightElecs_->Branch("Phi",          ev.te_phi,      "Phi[ElectronTight_size]/F");
  t_tightElecs_->Branch("Mass",         ev.te_mass,     "Mass[ElectronTight_size]/F");
  t_tightElecs_->Branch("IsolationVar", ev.te_relIso,   "IsolationVar[ElectronTight_size]/F");

  t_looseMuons_->Branch("MuonLoose_size", &ev.nlm,      "MuonLoose_size/I");
  t_looseMuons_->Branch("Charge",       ev.lm_ch,       "Charge[MuonLoose_size]/I");
  t_looseMuons_->Branch("Particle",     ev.lm_g,        "Particle[MuonLoose_size]/I");
  t_looseMuons_->Branch("PT",           ev.lm_pt,       "PT[MuonLoose_size]/F");
  t_looseMuons_->Branch("Eta",          ev.lm_eta,      "Eta[MuonLoose_size]/F");
  t_looseMuons_->Branch("Phi",          ev.lm_phi,      "Phi[MuonLoose_size]/F");
  t_looseMuons_->Branch("Mass",         ev.lm_mass,     "Mass[MuonLoose_size]/F");
  t_looseMuons_->Branch("IsolationVar", ev.lm_relIso,   "IsolationVar[MuonLoose_size]/F");

  t_tightMuons_->Branch("MuonTight_size", &ev.ntm,      "MuonTight_size/I");
  t_tightMuons_->Branch("Charge",       ev.tm_ch,       "Charge[MuonTight_size]/I");
  t_tightMuons_->Branch("Particle",     ev.tm_g,        "Particle[MuonTight_size]/I");
  t_tightMuons_->Branch("PT",           ev.tm_pt,       "PT[MuonTight_size]/F");
  t_tightMuons_->Branch("Eta",          ev.tm_eta,      "Eta[MuonTight_size]/F");
  t_tightMuons_->Branch("Phi",          ev.tm_phi,      "Phi[MuonTight_size]/F");
  t_tightMuons_->Branch("Mass",         ev.tm_mass,     "Mass[MuonTight_size]/F");
  t_tightMuons_->Branch("IsolationVar", ev.tm_relIso,   "IsolationVar[MuonTight_size]/F");

  t_puppiJets_->Branch("JetPUPPI_size", &ev.nj,         "JetPUPPI_size/I");
  t_puppiJets_->Branch("ID",            ev.j_id,        "ID[JetPUPPI_size]/I");
  t_puppiJets_->Branch("GenJet",        ev.j_g,         "GenJet[JetPUPPI_size]/I");
  t_puppiJets_->Branch("PT",            ev.j_pt,        "PT[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Eta",           ev.j_eta,       "Eta[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Phi",           ev.j_phi,       "Phi[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Mass",          ev.j_mass,      "Mass[JetPUPPI_size]/F");
  t_puppiJets_->Branch("MVAv2",         ev.j_mvav2,     "MVAv2[JetPUPPI_size]/I");
  t_puppiJets_->Branch("DeepCSV",       ev.j_deepcsv,   "DeepCSV[JetPUPPI_size]/I");
  t_puppiJets_->Branch("PartonFlavor",  ev.j_flav,      "PartonFlavor[JetPUPPI_size]/I");
  t_puppiJets_->Branch("HadronFlavor",  ev.j_hadflav,   "HadronFlavor[JetPUPPI_size]/I");
  t_puppiJets_->Branch("GenPartonPID",  ev.j_pid,       "GenPartonPID[JetPUPPI_size]/I");
  t_puppiJets_->Branch("chargedSumPtConst", ev.j_chargedSumPtConst, "chargedSumPtConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("neutralSumPtConst", ev.j_neutralSumPtConst, "neutralSumPtConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("hfemSumPtConst", ev.j_hfemSumPtConst, "hfemSumPtConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("hfhadSumPtConst", ev.j_hfhadSumPtConst, "hfhadSumPtConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("chargedNConst", ev.j_chargedNConst, "chargedNConst[JetPUPPI_size]/I");
  t_puppiJets_->Branch("neutralNConst", ev.j_neutralNConst, "neutralNConst[JetPUPPI_size]/I");
  t_puppiJets_->Branch("hfemNConst", ev.j_hfemNConst, "hfemNConst[JetPUPPI_size]/I");
  t_puppiJets_->Branch("hfhadNConst", ev.j_hfhadNConst, "hfhadNConst[JetPUPPI_size]/I");
  t_puppiJets_->Branch("eSumPtConst", ev.j_eSumPtConst, "eSumPtConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("eNConst", ev.j_eNConst, "eNConst[JetPUPPI_size]/I");
  t_puppiJets_->Branch("muSumPtConst", ev.j_muSumPtConst, "muSumPtConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("muNConst", ev.j_muNConst, "muNConst[JetPUPPI_size]/I");
  t_puppiJets_->Branch("photonSumPtConst", ev.j_photonSumPtConst, "photonSumPtConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("photonNConst", ev.j_photonNConst, "photonNConst[JetPUPPI_size]/I");
  t_puppiJets_->Branch("chargedSumPuppiWConst", ev.j_chargedSumPuppiWConst, "chargedSumPuppiWConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("neutralSumPuppiWConst", ev.j_neutralSumPuppiWConst, "neutralSumPuppiWConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("hfemSumPuppiWConst", ev.j_hfemSumPuppiWConst, "hfemSumPuppiWConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("hfhadSumPuppiWConst", ev.j_hfhadSumPuppiWConst, "hfhadSumPuppiWConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("eSumPuppiWConst", ev.j_eSumPuppiWConst, "eSumPuppiWConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("muSumPuppiWConst", ev.j_muSumPuppiWConst, "muSumPuppiWConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("photonSumPuppiWConst", ev.j_photonSumPuppiWConst, "photonSumPuppiWConst[JetPUPPI_size]/F");
  t_puppiJets_->Branch("RMSCand", ev.j_RMSCand, "RMSCand[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Axis1", ev.j_Axis1, "Axis1[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Axis2", ev.j_Axis2, "Axis2[JetPUPPI_size]/F");
  t_puppiJets_->Branch("Sigma", ev.j_Sigma, "Sigma[JetPUPPI_size]/F");
  t_puppiJets_->Branch("ptD", ev.j_ptD, "ptD[JetPUPPI_size]/F");

  t_puppiMET_->Branch("PuppiMissingET_size", &ev.nmet,  "PuppiMissingET_size/I");
  t_puppiMET_->Branch("MET",            ev.met_pt,      "MET[PuppiMissingET_size]/F");
  t_puppiMET_->Branch("Phi",            ev.met_phi,     "Phi[PuppiMissingET_size]/F");
  t_puppiMET_->Branch("Eta",            ev.met_eta,     "Eta[PuppiMissingET_size]/F");

  t_loosePhotons_->Branch("PhotonLoose_size", &ev.nlp,     "PhotonLoose_size/I");
  t_loosePhotons_->Branch("Particle",     ev.lp_g,         "Particle[PhotonLoose_size]/I");
  t_loosePhotons_->Branch("PT",           ev.lp_pt,        "PT[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Eta",          ev.lp_eta,       "Eta[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Phi",          ev.lp_phi,       "Phi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("E",            ev.lp_nrj,       "E[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("PT_multi",     ev.lp_pt_multi,  "PT_multi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Eta_multi",    ev.lp_eta_multi, "Eta_multi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("Phi_multi",    ev.lp_phi_multi, "Phi_multi[PhotonLoose_size]/F");
  t_loosePhotons_->Branch("E_multi",      ev.lp_nrj_multi, "E_multi[PhotonLoose_size]/F");

  t_tightPhotons_->Branch("PhotonTight_size", &ev.ntp,     "PhotonTight_size/I");
  t_tightPhotons_->Branch("Particle",     ev.tp_g,         "Particle[PhotonTight_size]/I");
  t_tightPhotons_->Branch("PT",           ev.tp_pt,        "PT[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Eta",          ev.tp_eta,       "Eta[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Phi",          ev.tp_phi,       "Phi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("E",            ev.tp_nrj,       "E[PhotonTight_size]/F");
  t_tightPhotons_->Branch("PT_multi",     ev.tp_pt_multi,  "PT_multi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Eta_multi",    ev.tp_eta_multi, "Eta_multi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("Phi_multi",    ev.tp_phi_multi, "Phi_multi[PhotonTight_size]/F");
  t_tightPhotons_->Branch("E_multi",      ev.tp_nrj_multi, "E_multi[PhotonTight_size]/F");
}

