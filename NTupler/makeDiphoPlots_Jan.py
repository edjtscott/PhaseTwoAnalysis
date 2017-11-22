#!/usr/bin/env python
# code to loop over HGCAL photon ntuples and make some plots

import os
import ROOT as r
from diphoHelpers import initHistByEta, fillHistByEta, getEffSigma, printHists, initHistByJetType, fillHistByJetType, drawJetHist

r.gROOT.SetBatch(True)

def main():
  #setup histos with several eta breakdown scenarios
  etaHists = {}
  initHistByEta(etaHists, 'mggTightID', 40, 115, 135)
  initHistByEta(etaHists, 'mggPtCutsTightID', 40, 115, 135)
  initHistByEta(etaHists, 'mggVBFPhaseSpace', 40, 115, 135)
  #vtx histos
  vtxHists = {}
  vtxHists['choseCorrectVtx'] = r.TH1F('choseCorrectVtx', 'choseCorrectVtx', 2, -0.5, 1.5)
  vtxHists['choseCorrectVtxVBFPhaseSpace'] = r.TH1F('choseCorrectVtxVBFPhaseSpace', 'choseCorrectVtxVBFPhaseSpace', 2, -0.5, 1.5)
  #histos in VBF phase space, for jet shape comparisons
  jetHists     = {}
  jetVals    = {}
  jetVariables = ['chargedSumPtConst', 'neutralSumPtConst', 'hfemSumPtConst', 'hfhadSumPtConst', 'chargedNConst', 'neutralNConst', 'hfemNConst', 'hfhadNConst', 'eSumPtConst', 'eNConst', 'muSumPtConst', 'muNConst', 'photonSumPtConst', 'photonNConst', 'RMSCand', 'Axis1', 'Axis2', 'Sigma', 'ptD']
  initHistByJetType(jetHists, 'jetPt')
  initHistByJetType(jetHists, 'jetEta')
  for jetVar in jetVariables:
    initHistByJetType(jetHists, jetVar)
    initHistByJetType(jetHists, jetVar+'_limited')
  #the 2D energy correction hist. FIXME: Binning could be changed
  correctionHists = {}
  correctionHists['photonECorrection'] = r.TH2F('photonECorrection','photonECorrection',6,0.,3.,25,30.,530.)
  correctionHists['photonECorrectionEntries'] = r.TH2F('photonECorrectionEntries','photonECorrectionEntries',6,0.,3.,25,30.,530.)
  correctionHists['photonECorrection'].Sumw2()
  correctionHists['photonECorrectionEntries'].Sumw2()

  #file lists
  nFiles = {}
  nFiles['VBF_PU0']   = 5
  nFiles['VBF_PU200'] = 23
  nFiles['ggH_PU0']   = 3
  nFiles['ggH_PU200'] = 8
  baseNames = {}
  baseNames['VBF_PU0']   = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/VBF_PU0_17Nov17/VBFHToGG_M125_14TeV_amcatnlo_pythia8/crab_VBFHToGG_M125_14TeV_amcatnlo_pythia8-PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/171117_164919/0000/MiniEvents_'
  baseNames['VBF_PU200'] = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/VBF_PU200_20Nov17/VBFHToGG_M125_14TeV_amcatnlo_pythia8/crab_VBFHToGG_M125_14TeV_amcatnlo_pythia8--PU200_93X_upgrade2023_realistic_v2-v2_2ndAttempt/171120_121557/0000/MiniEvents_'
  baseNames['ggH_PU0']   = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/ggH_PU0_20Nov17/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/crab_ggHToGG_M125_14TeV_amcatnloFXFX_pythia8--noPU_93X_upgrade2023_realistic_v2-v1_2ndAttempt/171120_121021/0000/MiniEvents_'
  baseNames['ggH_PU200'] = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/ggH_PU200_20Nov17/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/crab_ggHToGG_M125_14TeV_amcatnloFXFX_pythia8--PU200_93X_upgrade2023_realistic_v2-v2_2ndAttempt/171120_121647/0000/MiniEvents_'
  fileVetos = {}
  fileVetos['VBF_PU0']   = []
  fileVetos['VBF_PU200'] = []
  fileVetos['ggH_PU0']   = []
  fileVetos['ggH_PU200'] = [3,5] #FIXME these might need updating

  #FIXME make configurable
  #theKey = 'VBF_PU0'
  #theKey = 'VBF_PU200'
  theKey = 'ggH_PU0'
  #theKey = 'ggH_PU200'
  theTree       = r.TChain('ntuple/PhotonTight')
  genPhotonTree = r.TChain('ntuple/GenPhoton')
  genVtxTree    = r.TChain('ntuple/GenVertex')
  recoVtxTree   = r.TChain('ntuple/Vertex')
  jetTree       = r.TChain('ntuple/JetPUPPI')
  for fileNum in range(1,nFiles[theKey]+1):
    if fileNum in fileVetos[theKey]: continue
    fileName = 'root://eoscms.cern.ch//%s%g.root'%(baseNames[theKey],fileNum)
    theTree.Add(fileName)
    genPhotonTree.Add(fileName)
    genVtxTree.Add(fileName)
    recoVtxTree.Add(fileName)
    jetTree.Add(fileName)

  print 'got trees from files'
  for i,ev in enumerate(theTree):
    #if i==1000: break #FIXME: just for testing
    if i==0: print 'Processing first event'
    elif i%10000==0: print 'Processing event %g'%i
    #check two photons
    nPhotons = getattr(ev,'PhotonTight_size')
    if not nPhotons==2: 
      continue
    #setup collections
    #gen and reco vertex z
    genVtxTree.GetEntry(i)
    genVtxZ  = getattr(genVtxTree,'Z')
    recoVtxTree.GetEntry(i)
    recoVtxZ = getattr(recoVtxTree,'Z')
    if len(recoVtxZ) < 1:
      continue
    recoVtxZ = recoVtxZ[0]
    correctVtx = abs(genVtxZ-recoVtxZ) < 1.
    vtxHists['choseCorrectVtx'].Fill(correctVtx)
    #get gen photon energies
    genPhotonTree.GetEntry(i)
    genPhotonE = getattr(genPhotonTree,'E')
    #check if photon is barrel or endcap
    photonIsEB = getattr(ev,'IsEB')
    #get gen photon index
    photonGenIndex = getattr(ev,'Particle')
    #default pat quantities
    photonPt  = getattr(ev,'PT')
    photonE   = getattr(ev,'E')
    photonEta = getattr(ev,'Eta')
    photonPhi = getattr(ev,'Phi')
    #multicluster quantities
    photonPtMulti  = getattr(ev,'PT_multi')
    photonEMulti   = getattr(ev,'E_multi')
    photonEtaMulti = getattr(ev,'Eta_multi')
    photonPhiMulti = getattr(ev,'Phi_multi')
    #jet quantities
    jetTree.GetEntry(i)
    nJets       = getattr(jetTree,'JetPUPPI_size')
    jetPt       = getattr(jetTree,'PT')
    jetEta      = getattr(jetTree,'Eta')
    jetPhi      = getattr(jetTree,'Phi')
    jetMass     = getattr(jetTree,'Mass')
    jetGenIndex = getattr(jetTree,'GenJet')
    jetGenPID   = getattr(jetTree,'GenPartonPID')
    jetID       = getattr(jetTree,"ID")
    for jetVar in jetVariables:
      jetVals[jetVar] = getattr(jetTree,jetVar)

    #processing (use default quantities for barrel, multi ones for endcap)
    leadIndex       = int(photonPt[1] > photonPt[0])
    leadGenIndex    = photonGenIndex[leadIndex]
    subleadIndex    = 1 - leadIndex
    subleadGenIndex = photonGenIndex[subleadIndex]
    if leadGenIndex < 0 or subleadGenIndex < 0:
      continue
    leadPhoton = r.TLorentzVector()
    if photonIsEB[leadIndex]:
      leadPhoton.SetPtEtaPhiM(photonPt[leadIndex], photonEta[leadIndex], photonPhi[leadIndex], 0)
    else:
      leadPhoton.SetPtEtaPhiM(photonPtMulti[leadIndex], photonEtaMulti[leadIndex], photonPhiMulti[leadIndex], 0)
    subleadPhoton = r.TLorentzVector()
    if photonIsEB[subleadIndex]:
      subleadPhoton.SetPtEtaPhiM(photonPt[subleadIndex], photonEta[subleadIndex], photonPhi[subleadIndex], 0)
    else:
      subleadPhoton.SetPtEtaPhiM(photonPtMulti[subleadIndex], photonEtaMulti[subleadIndex], photonPhiMulti[subleadIndex], 0)

    #FIXME recalculate four-vectors if reco vertex further than 1cm from gen vertex (skipping for now)
    if not correctVtx: 
      continue

    #2D energy correction hist stuff can come first
    leadPhotonCorrection    = genPhotonE[leadGenIndex] / leadPhoton.E()
    subleadPhotonCorrection = genPhotonE[subleadGenIndex] / subleadPhoton.E()
    #print 'lead, sublead correction factors = %1.3f, %1.3f'%(leadPhotonCorrection,subleadPhotonCorrection)
    correctionHists['photonECorrection'].Fill(abs(leadPhoton.Eta()),    leadPhoton.E(),    leadPhotonCorrection)
    correctionHists['photonECorrection'].Fill(abs(subleadPhoton.Eta()), subleadPhoton.E(), subleadPhotonCorrection)
    correctionHists['photonECorrectionEntries'].Fill(abs(leadPhoton.Eta()),    leadPhoton.E())
    correctionHists['photonECorrectionEntries'].Fill(abs(subleadPhoton.Eta()), subleadPhoton.E())

    #then diphoton mass plots
    theDiphoton = leadPhoton + subleadPhoton
    diphoMass = theDiphoton.M()
    fillHistByEta(etaHists, 'mggTightID', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)
    if 3*leadPhoton.Pt()>diphoMass and 4*subleadPhoton.Pt()>diphoMass:
      fillHistByEta(etaHists, 'mggPtCutsTightID', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)

    #find lead and sublead jets
    if nJets < 2:
      continue
    leadJetIndex = -1
    leadJetPt    = -1.
    subleadJetIndex = -1
    subleadJetPt    = -1.
    for iJet in range(nJets):
      if not jetID[iJet]:
        continue
      tempPt = jetPt[iJet]
      if tempPt>leadJetPt:
        leadJetPt = jetPt[iJet]
        leadJetIndex = iJet
    for iJet2 in range(nJets):
      if iJet2 == leadJetIndex:
        continue
      if not jetID[iJet2]:
        continue
      tempPt = jetPt[iJet2]
      if tempPt>subleadJetPt:
        subleadJetPt = jetPt[iJet2]
        subleadJetIndex = iJet2
    if not (leadJetIndex>-1 and subleadJetIndex>-1):
      continue
    leadJet = r.TLorentzVector()
    leadJet.SetPtEtaPhiM(leadJetPt, jetEta[leadJetIndex], jetPhi[leadJetIndex], jetMass[leadJetIndex])
    subleadJet = r.TLorentzVector()
    subleadJet.SetPtEtaPhiM(subleadJetPt, jetEta[subleadJetIndex], jetPhi[subleadJetIndex], jetMass[subleadJetIndex])
    theDijet = leadJet + subleadJet
    dijetMass = theDijet.M()
    leadJetGenPID      = jetGenPID[leadJetIndex]
    subleadJetGenPID   = jetGenPID[subleadJetIndex]
    leadJetGenIndex    = jetGenIndex[leadJetIndex]
    subleadJetGenIndex = jetGenIndex[subleadJetIndex]
    #jetDeta = abs(jetEta[leadJetIndex]l - jetEta[subleadJetIndex]) #FIXME: no dEta cut in flashgg?
    #then project to VBF phase space
    if not (3*leadPhoton.Pt()>diphoMass and 4*subleadPhoton.Pt()>diphoMass):
      continue
    if not (abs(jetEta[leadJetIndex])<4.7 and abs(jetEta[subleadJetIndex])<4.7):
      continue
    if not (dijetMass > 250.):
      continue
    if not (jetID[leadJetIndex] and jetID[subleadJetIndex]): #FIXME need to investigate jetID
      continue
    #make jet, other plots here
    vtxHists['choseCorrectVtxVBFPhaseSpace'].Fill(correctVtx)
    fillHistByEta(etaHists, 'mggVBFPhaseSpace', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)
    fillHistByJetType(jetHists, 'jetPt', leadJetGenIndex, leadJetGenPID, 0, leadJet.Pt())
    fillHistByJetType(jetHists, 'jetPt', subleadJetGenIndex, subleadJetGenPID, 1, subleadJet.Pt())
    fillHistByJetType(jetHists, 'jetEta', leadJetGenIndex, leadJetGenPID, 0, leadJet.Eta())
    fillHistByJetType(jetHists, 'jetEta', subleadJetGenIndex, subleadJetGenPID, 1, subleadJet.Eta())
    for jetVar in jetVariables:
      if 'Axis' in jetVar:
        if jetVals[jetVar][leadJetIndex]>-1.: 
          fillHistByJetType(jetHists, jetVar, leadJetGenIndex,    leadJetGenPID,    0, jetVals[jetVar][leadJetIndex])
        if jetVals[jetVar][subleadJetIndex]>-1.: 
          fillHistByJetType(jetHists, jetVar, subleadJetGenIndex, subleadJetGenPID, 1, jetVals[jetVar][subleadJetIndex])
      else:
        fillHistByJetType(jetHists, jetVar, leadJetGenIndex,    leadJetGenPID,    0, jetVals[jetVar][leadJetIndex])
        fillHistByJetType(jetHists, jetVar, subleadJetGenIndex, subleadJetGenPID, 1, jetVals[jetVar][subleadJetIndex])
      if leadJetPt>50. and leadJetPt<100.:
        fillHistByJetType(jetHists, jetVar+'_limited', leadJetGenIndex,    leadJetGenPID,    0, jetVals[jetVar][leadJetIndex])
      if subleadJetPt>50. and subleadJetPt<100.:
        fillHistByJetType(jetHists, jetVar+'_limited', subleadJetGenIndex, subleadJetGenPID, 1, jetVals[jetVar][subleadJetIndex])

  #post-processing
  correctionHists['photonECorrection'].Divide(correctionHists['photonECorrectionEntries'])
  correctionHists['photonECorrection'].SetMinimum(0.)
  correctionHists['photonECorrection'].SetMaximum(3.)

  #draw hists, send to web
  canv = r.TCanvas('canv','canv')
  outdirName = 'DiphoPlots_%s/'%theKey
  os.system('mkdir -p %s'%outdirName)
  webDir = '/afs/cern.ch/user/e/escott/www/HFuture/Pass1/%s'%theKey
  os.system('mkdir -p %s'%webDir)
  os.system('cp /afs/cern.ch/user/e/escott/www/HFuture/Pass1/index.php %s'%webDir)
  printHists(canv, etaHists, outdirName)
  printHists(canv, correctionHists, outdirName)
  printHists(canv, vtxHists, outdirName)
  printHists(canv, jetHists, outdirName)
  drawJetHist(canv, jetHists, 'jetPt', outdirName)
  drawJetHist(canv, jetHists, 'jetEta', outdirName)
  for jetVar in jetVariables:
    drawJetHist(canv, jetHists, jetVar, outdirName)
    drawJetHist(canv, jetHists, jetVar+'_limited', outdirName)
  os.system('cp %s* %s'%(outdirName,webDir))
  print 'plots moved to %s'%webDir
  #save the correction hist for combination
  outFile = r.TFile('CorrectionHists/corrHist_%s.root'%theKey,'RECREATE')
  correctionHists['photonECorrection'].Write()
  outFile.Close()


if __name__ == '__main__':
  main()
