#!/usr/bin/env python
# code to loop over HGCAL photon ntuples and make some plots

import os
from numpy import sin, cos, tan, sqrt
import ROOT as r
from diphoHelpers import initHistByEta, fillHistByEta, getEffSigma, printHists, initHistByJetType, fillHistByJetType, drawJetHist, deltaR

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-k', '--key', default='VBF_PU200', help='choose the sample to run on')
parser.add_option('-d', '--doLoose', default=False, action='store_true', help='use loose photons (default false, ie use only tight photons)')
parser.add_option('-w', '--writePlots', default=False, action='store_true', help='send plots to web etc')
parser.add_option('-c', '--applyCorrections', default=True, action='store_false', help='apply photon corrections (default true)')
parser.add_option('-m', '--maxEvents', type='int', default=-1, help='specify number of events on which to run')
parser.add_option('-l', '--lumi', type='float', default=35.9, help='specify luminosity (in /fb)')
(opts,args) = parser.parse_args()

r.gROOT.SetBatch(True)

def main():
  #setup histos with several eta breakdown scenarios
  etaHists = {}
  initHistByEta(etaHists, 'mggTightID', 40, 115, 135)
  initHistByEta(etaHists, 'mggPtCutsTightID', 40, 115, 135)
  initHistByEta(etaHists, 'mggVBFPhaseSpace', 40, 115, 135)
  #misc histos
  miscHists = {}
  miscHists['choseCorrectVtx'] = r.TH1F('choseCorrectVtx', 'choseCorrectVtx', 2, -0.5, 1.5)
  miscHists['choseCorrectVtxVBFPhaseSpace'] = r.TH1F('choseCorrectVtxVBFPhaseSpace', 'choseCorrectVtxVBFPhaseSpace', 2, -0.5, 1.5)
  miscHists['cutFlow0_TotEvts'] = r.TH1F('cutFlow0_TotEvts', 'cutFlow0_TotEvts', 1, 0.5, 1.5)
  miscHists['cutFlow1_GtTwoPhotons'] = r.TH1F('cutFlow1_GtTwoPhotons', 'cutFlow1_GtTwoPhotons', 1, 0.5, 1.5)
  miscHists['cutFlow1a_RecoVtxExists'] = r.TH1F('cutFlow1a_RecoVtxExists', 'cutFlow1a_RecoVtxExists', 1, 0.5, 1.5)
  miscHists['cutFlow2_TwoMatchingPhotons'] = r.TH1F('cutFlow2_TwoMatchingPhotons', 'cutFlow2_TwoMatchingPhotons', 1, 0.5, 1.5)
  miscHists['cutFlow2a_Mgg100to180'] = r.TH1F('cutFlow2a_cutFlow2a_Mgg100to180', 'cutFlow2a_cutFlow2a_Mgg100to180', 1, 0.5, 1.5)
  miscHists['cutFlow3_ScaledPtCuts'] = r.TH1F('cutFlow3_ScaledPtCuts', 'cutFlow3_ScaledPtCuts', 1, 0.5, 1.5)
  miscHists['cutFlow4_TwoValidJetsIDd'] = r.TH1F('cutFlow4_TwoValidJetsIDd', 'cutFlow4_TwoValidJetsIDd', 1, 0.5, 1.5)
  miscHists['cutFlow5_PtCutJets'] = r.TH1F('cutFlow5_PtCutJets', 'cutFlow5_PtCutJets', 1, 0.5, 1.5)
  miscHists['cutFlow6_DijetMassCut'] = r.TH1F('cutFlow6_DijetMassCut', 'cutFlow6_DijetMassCut', 1, 0.5, 1.5)
  miscHists['cutFlow7_EffSigmaWindow'] = r.TH1F('cutFlow7_EffSigmaWindow', 'cutFlow7_EffSigmaWindow', 1, 0.5, 1.5)
  #histos in VBF phase space, for jet shape comparisons
  jetHists     = {}
  jetVals    = {}
  jetVariables = ['chargedSumPtConst', 'neutralSumPtConst', 'hfemSumPtConst', 'hfhadSumPtConst', 'chargedNConst', 'neutralNConst', 'hfemNConst', 'hfhadNConst', 'eSumPtConst', 'eNConst', 'muSumPtConst', 'muNConst', 'photonSumPtConst', 'photonNConst', 'RMSCand', 'Axis1', 'Axis2', 'Sigma', 'ptD']
  initHistByJetType(jetHists, 'jetPt')
  initHistByJetType(jetHists, 'jetEta')
  for jetVar in jetVariables:
    initHistByJetType(jetHists, jetVar)
    initHistByJetType(jetHists, jetVar+'_limited')
  #the 2D energy correction hist. FIXME: Binning could be changed, ideally variable
  correctionHists = {}
  correctionHists['photonECorrection'] = r.TH2F('photonECorrection','photonECorrection',6,0.,3.,25,30.,530.)
  correctionHists['photonECorrectionEntries'] = r.TH2F('photonECorrectionEntries','photonECorrectionEntries',6,0.,3.,25,30.,530.)
  correctionHists['photonECorrection'].Sumw2()
  correctionHists['photonECorrectionEntries'].Sumw2()
  correctionHists['photonECorrectionClosure'] = r.TH2F('photonECorrectionClosure','photonECorrectionClosure',6,0.,3.,25,30.,530.)
  correctionHists['photonECorrectionClosureEntries'] = r.TH2F('photonECorrectionClosureEntries','photonECorrectionClosureEntries',6,0.,3.,25,30.,530.)
  correctionHists['photonECorrectionClosure'].Sumw2()
  correctionHists['photonECorrectionClosureEntries'].Sumw2()

  if opts.applyCorrections:
    corrFile = r.TFile('combinedCorrectionHist.root','READ')
    corrHist = corrFile.Get('combinedECorrection')

  #file lists, cross-sections etc
  nFiles = {}
  nFiles['VBF_PU0']   = 5
  nFiles['VBF_PU200'] = 23
  nFiles['ggH_PU0']   = 3
  nFiles['ggH_PU200'] = 8
  nFiles['GJet_PU200'] = 345
  nFiles['DiPhoJetsBox_PU200'] = 72
  baseNames = {}
  baseNames['VBF_PU0']   = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/VBF_PU0_17Nov17/VBFHToGG_M125_14TeV_amcatnlo_pythia8/crab_VBFHToGG_M125_14TeV_amcatnlo_pythia8-PhaseIITDRFall17MiniAOD-noPU_93X_upgrade2023_realistic_v2-v1/171117_164919/0000/MiniEvents_'
  baseNames['VBF_PU200'] = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/VBF_PU200_20Nov17/VBFHToGG_M125_14TeV_amcatnlo_pythia8/crab_VBFHToGG_M125_14TeV_amcatnlo_pythia8--PU200_93X_upgrade2023_realistic_v2-v2_2ndAttempt/171120_121557/0000/MiniEvents_'
  baseNames['ggH_PU0']   = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/ggH_PU0_20Nov17/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/crab_ggHToGG_M125_14TeV_amcatnloFXFX_pythia8--noPU_93X_upgrade2023_realistic_v2-v1_2ndAttempt/171120_121021/0000/MiniEvents_'
  baseNames['ggH_PU200'] = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/ggH_PU200_20Nov17/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/crab_ggHToGG_M125_14TeV_amcatnloFXFX_pythia8--PU200_93X_upgrade2023_realistic_v2-v2_2ndAttempt/171120_121647/0000/MiniEvents_'
  baseNames['GJet_PU200'] = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/GJet_PU200_23Nov17/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/crab_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_14TeV--PU200_93X_upgrade2023_realistic_v2-v2/171123_232640/0000/MiniEvents_'
  baseNames['DiPhoJetsBox_PU200'] = 'eos/cms/store/group/dpg_hgcal/comm_hgcal/escott/DiPhotonJetsBox_PU200_23Nov17/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/crab_DiPhotonJetsBox_MGG-80toInf_14TeV--PU200_93X_upgrade2023_realistic_v2-v2/171123_232953/0000/MiniEvents_'
  fileVetos = {}
  fileVetos['VBF_PU0']   = []
  fileVetos['VBF_PU200'] = []
  fileVetos['ggH_PU0']   = []
  fileVetos['ggH_PU200'] = [3] #FIXME these might need updating
  fileVetos['GJet_PU200'] = []
  fileVetos['DiPhoJetsBox_PU200'] = []
  xSecs = {} #in pb
  brHgg = 0.00227
  xSecs['VBF_PU0']   = 4.278*brHgg
  xSecs['VBF_PU200'] = 4.278*brHgg
  xSecs['ggH_PU0']   = 54.67*brHgg
  xSecs['ggH_PU200'] = 54.67*brHgg
  xSecs['GJet_PU200']        = 950.3
  xSecs['DiPhoJetsBox_PU200'] = 95.81
 

  theKey = opts.key
  if theKey not in nFiles.keys():
    raise Exception('invalid key: use one of %s'%(nFiles.keys()))
  print 'Running on %s events'%theKey
  if opts.maxEvents>0: print 'Up to a maximum of %g events'%opts.maxEvents
  print 'Write plots is set to',opts.writePlots
  print 'and applyCorrections is',opts.applyCorrections
  doLoose = opts.doLoose
  if not doLoose: 
    theTree     = r.TChain('ntuple/PhotonTight')
  else: 
    theTree     = r.TChain('ntuple/PhotonLoose')
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

  #flatTree = r.TTree('flatTree_%s'%theKey,'flatTree_%s'%theKey)
  #leadPtOvMgg = -9999.
  #flatTree.Branch('leadPtOvMgg',leadPtOvMgg,'leadPtOvMgg/D')

  print 'got trees from files'
  nEvents = theTree.GetEntries()
  print 'found %g events'%nEvents
  if opts.maxEvents > 0 and nEvents > opts.maxEvents:
    nEvents = opts.maxEvents
    print 'but max is %g so will only use that many'%opts.maxEvents
  nExpEvents = 1000. * opts.lumi * xSecs[theKey]
  evtWeight = nExpEvents / float(nEvents)
  print 'expecting %1.0f events'%nExpEvents
  print 'applying event weight of %1.3f to cutflow hists'%evtWeight
  for i,ev in enumerate(theTree):
    if i==0: print 'Processing first event'
    elif i%10000==0: print 'Processing event %g'%i
    if i==opts.maxEvents: break
    miscHists['cutFlow0_TotEvts'].Fill(1,evtWeight)
    #check two photons
    if not doLoose:
      nPhotons = getattr(ev,'PhotonTight_size')
    else:
      nPhotons = getattr(ev,'PhotonLoose_size')
    if nPhotons < 2:
      continue
    miscHists['cutFlow1_GtTwoPhotons'].Fill(1,evtWeight)
    #setup collections
    #gen and reco vertex z
    genVtxTree.GetEntry(i)
    genVtxZ  = getattr(genVtxTree,'Z')
    recoVtxTree.GetEntry(i)
    recoVtxZ = getattr(recoVtxTree,'Z')
    if len(recoVtxZ) < 1:
      continue
    miscHists['cutFlow1a_RecoVtxExists'].Fill(1,evtWeight)
    recoVtxZ = recoVtxZ[0]
    correctVtx = abs(genVtxZ-recoVtxZ) < 1.
    miscHists['choseCorrectVtx'].Fill(correctVtx)
    #get gen photon energies
    genPhotonTree.GetEntry(i)
    genPhotonE = getattr(genPhotonTree,'E')
    genPhotonEta = getattr(genPhotonTree,'Eta')
    genPhotonPhi = getattr(genPhotonTree,'Phi')
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
    photonZMulti = getattr(ev,'Z_multi')
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
    #leadIndex       = int(photonPt[1] > photonPt[0])
    #leadGenIndex    = photonGenIndex[leadIndex]
    #subleadIndex    = 1 - leadIndex
    #subleadGenIndex = photonGenIndex[subleadIndex]
    leadIndex    = -1
    leadPt       = -1.
    subleadIndex = -1.
    subleadPt    = -1.
    for iPho in range(nPhotons):
      tempPt = photonPt[iPho]
      if abs(photonEta[iPho]) > 3.:
        continue
      elif abs(photonEta[iPho]) > 1.5:
        tempPt = photonPtMulti[iPho]
      if tempPt > leadPt:
        leadPt = tempPt
        leadIndex = iPho
    for iPho2 in range(nPhotons):
      if iPho2 == leadIndex:
        continue
      tempPt = photonPt[iPho2]
      if abs(photonEta[iPho2]) > 3.:
        continue
      elif abs(photonEta[iPho2]) > 1.5:
        tempPt = photonPtMulti[iPho2]
      if tempPt > subleadPt:
        subleadPt = tempPt
        subleadIndex = iPho2
    if leadIndex < 0 or subleadIndex < 0:
      continue
    leadGenIndex    = photonGenIndex[leadIndex]
    subleadGenIndex = photonGenIndex[subleadIndex]
    if leadGenIndex < 0 or subleadGenIndex < 0:
      continue
    if leadPt < 30. or subleadPt < 30.:
      continue
    miscHists['cutFlow2_TwoMatchingPhotons'].Fill(1,evtWeight)
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
    theDiphoton = leadPhoton + subleadPhoton

    #recalculate four-vectors if reco vertex further than 1cm from gen vertex
    if not correctVtx: 
      #print 'Original lead photon:    ',leadPhoton.Print()
      #print 'Original sublead photon: ',subleadPhoton.Print()
      #print 'Original diphoton mass:  ',theDiphoton.M()
      leadZ = photonZMulti[leadIndex]
      leadX = leadZ * tan(leadPhoton.Theta()) * cos(leadPhoton.Phi())
      leadY = leadZ * tan(leadPhoton.Theta()) * sin(leadPhoton.Phi())
      leadDirection = r.Math.XYZVector(leadX, leadY, leadZ - genVtxZ)
      leadFourVector = r.Math.XYZVector.Unit(leadDirection) * leadPhoton.E()
      leadPhoton.SetPxPyPzE(leadFourVector.X(), leadFourVector.Y(), leadFourVector.Z(), leadPhoton.E())
      subleadZ = photonZMulti[subleadIndex]
      subleadX = subleadZ * tan(subleadPhoton.Theta()) * cos(subleadPhoton.Phi())
      subleadY = subleadZ * tan(subleadPhoton.Theta()) * sin(subleadPhoton.Phi())
      subleadDirection = r.Math.XYZVector(subleadX, subleadY, subleadZ - genVtxZ)
      subleadFourVector = r.Math.XYZVector.Unit(subleadDirection) * subleadPhoton.E()
      subleadPhoton.SetPxPyPzE(subleadFourVector.X(), subleadFourVector.Y(), subleadFourVector.Z(), subleadPhoton.E())
      theDiphoton = leadPhoton + subleadPhoton
      #print 'New lead photon:    ',leadPhoton.Print()
      #print 'New sublead photon: ',subleadPhoton.Print()
      #print 'New diphoton mass:  ',theDiphoton.M()
      #print '\n\n'

    #2D energy correction hist stuff can come first
    leadPhotonCorrection    = genPhotonE[leadGenIndex] / leadPhoton.E()
    subleadPhotonCorrection = genPhotonE[subleadGenIndex] / subleadPhoton.E()
    #print 'lead, sublead correction factors = %1.3f, %1.3f'%(leadPhotonCorrection,subleadPhotonCorrection)
    correctionHists['photonECorrection'].Fill(abs(leadPhoton.Eta()),    leadPhoton.E(),    leadPhotonCorrection)
    correctionHists['photonECorrection'].Fill(abs(subleadPhoton.Eta()), subleadPhoton.E(), subleadPhotonCorrection)
    correctionHists['photonECorrectionEntries'].Fill(abs(leadPhoton.Eta()),    leadPhoton.E())
    correctionHists['photonECorrectionEntries'].Fill(abs(subleadPhoton.Eta()), subleadPhoton.E())

    #then diphoton mass plots
    diphoMass = theDiphoton.M()
    if opts.applyCorrections:
      corrBinsX, corrBinsY = corrHist.GetNbinsX(), corrHist.GetNbinsY()
      corrXlow, corrXhigh = corrHist.GetXaxis().GetXmin(), corrHist.GetXaxis().GetXmax()
      corrYlow, corrYhigh = corrHist.GetYaxis().GetXmin(), corrHist.GetYaxis().GetXmax()
      leadCorrBinX = int( (abs(leadPhoton.Eta()) - corrXlow) * corrBinsX / (corrXhigh - corrXlow) ) + 1
      leadCorrBinY = int( (leadPhoton.E() - corrYlow) * corrBinsY / (corrYhigh - corrYlow) ) + 1
      leadCorrValue = corrHist.GetBinContent( leadCorrBinX, leadCorrBinY )
      subleadCorrBinX = int( (abs(subleadPhoton.Eta()) - corrXlow) * corrBinsX / (corrXhigh - corrXlow) ) + 1
      subleadCorrBinY = int( (subleadPhoton.E() - corrYlow) * corrBinsY / (corrYhigh - corrYlow) ) + 1
      subleadCorrValue = corrHist.GetBinContent( subleadCorrBinX, subleadCorrBinY )
      #print 'eta, E, binX, binY, corr: %1.3f, %1.3f, %g, %g, %1.3f'%(leadPhoton.Eta(),leadPhoton.E(),leadCorrBinX,leadCorrBinY,leadCorrValue)
      #print 'eta, E, binX, binY, corr: %1.3f, %1.3f, %g, %g, %1.3f\n'%(subleadPhoton.Eta(),subleadPhoton.E(),subleadCorrBinX,subleadCorrBinY,subleadCorrValue)
      diphoMass *= sqrt( leadCorrValue * subleadCorrValue )
      if leadCorrValue > 0.:
        correctionHists['photonECorrectionClosure'].Fill(abs(leadPhoton.Eta()),       leadPhoton.E(),       genPhotonE[leadGenIndex]    / (leadPhoton.E() *    leadCorrValue))
        correctionHists['photonECorrectionClosureEntries'].Fill(abs(leadPhoton.Eta()),       leadPhoton.E())
      if subleadCorrValue > 0.:
        correctionHists['photonECorrectionClosure'].Fill(abs(subleadPhoton.Eta()),    subleadPhoton.E(),    genPhotonE[subleadGenIndex] / (subleadPhoton.E() * subleadCorrValue))
        correctionHists['photonECorrectionClosureEntries'].Fill(abs(subleadPhoton.Eta()),    subleadPhoton.E())
    if diphoMass < 100. or diphoMass > 180.:
      continue
    miscHists['cutFlow2a_Mgg100to180'].Fill(1,evtWeight)
    fillHistByEta(etaHists, 'mggTightID', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)
    if not (3*leadPhoton.Pt()>diphoMass and 4*subleadPhoton.Pt()>diphoMass):
      continue
    miscHists['cutFlow3_ScaledPtCuts'].Fill(1,evtWeight)
    fillHistByEta(etaHists, 'mggPtCutsTightID', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)

    #find lead and sublead jets
    if nJets < 2:
      continue
    leadJetIndex = -1
    leadJetPt    = -1.
    subleadJetIndex = -1
    subleadJetPt    = -1.
    for iJet in range(nJets):
      #if not jetID[iJet]:
      if not jetID[iJet] > 2:
        continue
      if not abs(jetEta[iJet]) < 4.7:
        continue
      photonMatch = False
      for iGenPho in range(len(genPhotonEta)):
        tempDeltaR = deltaR(jetEta[iJet], genPhotonEta[iGenPho], jetPhi[iJet], genPhotonPhi[iGenPho])
        if tempDeltaR<0.4: 
          photonMatch = True
      if photonMatch:
        continue
      tempPt = jetPt[iJet]
      if tempPt>leadJetPt:
        leadJetPt = jetPt[iJet]
        leadJetIndex = iJet
    for iJet2 in range(nJets):
      if iJet2 == leadJetIndex:
        continue
      #if not jetID[iJet2]:
      if not jetID[iJet2] > 2:
        continue
      if not abs(jetEta[iJet2]) < 4.7:
        continue
      photonMatch = False
      for iGenPho2 in range(len(genPhotonEta)):
        tempDeltaR = deltaR(jetEta[iJet2], genPhotonEta[iGenPho2], jetPhi[iJet2], genPhotonPhi[iGenPho2])
        if tempDeltaR<0.4: 
          photonMatch = True
      if photonMatch:
        continue
      tempPt = jetPt[iJet2]
      if tempPt>subleadJetPt:
        subleadJetPt = jetPt[iJet2]
        subleadJetIndex = iJet2
    if not (leadJetIndex>-1 and subleadJetIndex>-1):
      continue
    miscHists['cutFlow4_TwoValidJetsIDd'].Fill(1,evtWeight)
    if not (leadJetPt>40. and subleadJetPt>30.):
      continue
    miscHists['cutFlow5_PtCutJets'].Fill(1,evtWeight)
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
    #jetDeta = abs(jetEta[leadJetIndex] - jetEta[subleadJetIndex]) #FIXME: no dEta cut in flashgg?
    #then project to VBF phase space
    if not (dijetMass > 250.):
      continue
    miscHists['cutFlow6_DijetMassCut'].Fill(1,evtWeight)
    #make jet, other plots here
    miscHists['choseCorrectVtxVBFPhaseSpace'].Fill(correctVtx)
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
    if diphoMass < 123 or diphoMass > 127: 
      continue
    miscHists['cutFlow7_EffSigmaWindow'].Fill(1,evtWeight)

  #post-processing
  print 'Processed last event'
  correctionHists['photonECorrection'].Divide(correctionHists['photonECorrectionEntries'])
  correctionHists['photonECorrection'].SetMinimum(0.)
  correctionHists['photonECorrection'].SetMaximum(3.)
  correctionHists['photonECorrectionClosure'].Divide(correctionHists['photonECorrectionClosureEntries'])
  correctionHists['photonECorrectionClosure'].SetMinimum(0.)
  correctionHists['photonECorrectionClosure'].SetMaximum(3.)

  #draw hists, send to web
  if opts.writePlots:
    canv = r.TCanvas('canv','canv')
    outdirName = 'DiphoPlots_%s/'%theKey
    if doLoose:
      outdirName = 'DiphoPlots_Loose_%s/'%theKey
    os.system('mkdir -p %s'%outdirName)
    webDir = '/afs/cern.ch/user/e/escott/www/HFuture/Pass1d/%s'%theKey
    if doLoose:
      webDir = '/afs/cern.ch/user/e/escott/www/HFuture/Pass1d/Loose/%s'%theKey
    os.system('mkdir -p %s'%webDir)
    os.system('cp /afs/cern.ch/user/e/escott/www/HFuture/Pass1/index.php %s'%webDir)
    printHists(canv, etaHists, outdirName)
    printHists(canv, correctionHists, outdirName)
    printHists(canv, miscHists, outdirName)
    #printHists(canv, jetHists, outdirName)
    drawJetHist(canv, jetHists, 'jetPt', outdirName)
    drawJetHist(canv, jetHists, 'jetEta', outdirName)
    for jetVar in jetVariables:
      drawJetHist(canv, jetHists, jetVar, outdirName)
      #drawJetHist(canv, jetHists, jetVar+'_limited', outdirName)
    os.system('cp %s* %s'%(outdirName,webDir))
    print 'plots moved to %s'%webDir
    os.system('mkdir -p EventYields/')
    yieldFile = open('EventYields/%s.txt'%theKey,'w')
    for cfName,cfHist in miscHists.iteritems():
      if not 'cutFlow' in cfName:
        continue
      yieldFile.write('nEvents at %s is %1.3f \n'%(cfName,cfHist.Integral()))
    yieldFile.close()
    #save the correction hists for combination
    os.system('mkdir -p CorrectionHists/')
    if not opts.applyCorrections:
      if not doLoose:
        outFile = r.TFile('CorrectionHists/corrHist_%s.root'%theKey,'RECREATE')
      else:
        outFile = r.TFile('CorrectionHistsLoose/corrHist_%s.root'%theKey,'RECREATE')
      correctionHists['photonECorrection'].Write()
      correctionHists['photonECorrectionEntries'].Write()
      correctionHists['photonECorrectionClosure'].Write()
      correctionHists['photonECorrectionClosureEntries'].Write()
      outFile.Close()
    else:
      corrFile.Close()


if __name__ == '__main__':
  main()
