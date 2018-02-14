#!/usr/bin/env python
# code to loop over HGCAL photon ntuples and make flat trees
# selection is loosened for training

import os
from numpy import sin, cos, tan, sqrt, array
import ROOT as r
from diphoHelpers import initHistByEta, fillHistByEta, getEffSigma, printHists, initHistByJetType, fillHistByJetType, drawJetHist, deltaR, deltaPhi

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-k', '--key', default='VBF_PU200', help='choose the sample to run on')
parser.add_option('-d', '--doLoose', default=False, action='store_true', help='use loose photons (default false, ie use only tight photons)')
parser.add_option('-w', '--writePlots', default=False, action='store_true', help='send plots to web etc')
parser.add_option('-c', '--applyCorrections', default=True, action='store_false', help='apply photon corrections (default true)')
parser.add_option('-m', '--maxEvents', type='int', default=-1, help='specify number of events on which to run')
parser.add_option('-t', '--maxTreeEvents', type='int', default=-1, help='max number of events to write to tree')
parser.add_option('-l', '--lumi', type='float', default=35.9, help='specify luminosity (in /fb)')
parser.add_option('--forTesting', default=False, action='store_true', help='apply normal VBF preselection, rather than loosened')
(opts,args) = parser.parse_args()

r.gROOT.SetBatch(True)

def main():
  jetVals    = {}
  jetVariables = ['chargedSumPtConst', 'neutralSumPtConst', 'hfemSumPtConst', 'hfhadSumPtConst', 'chargedNConst', 'neutralNConst', 'hfemNConst', 'hfhadNConst', 'eSumPtConst', 'eNConst', 'muSumPtConst', 'muNConst', 'photonSumPtConst', 'photonNConst', 'RMSCand', 'Axis1', 'Axis2', 'Sigma', 'ptD']

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
  xSecs['GJet_PU200']         = 950.3
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

  #attempt to write a flat tree with desired variables in VBF phase space...
  if not opts.forTesting:
    flatFile = r.TFile('/afs/cern.ch/work/e/escott/public/HFuture/flatTreeForTraining_%s.root'%theKey,'RECREATE')
    flatTree = r.TTree('flatTreeForTraining','flatTreeForTraining')
  else:
    flatFile = r.TFile('/afs/cern.ch/work/e/escott/public/HFuture/flatTreeForTesting_%s.root'%theKey,'RECREATE')
    flatTree = r.TTree('flatTreeForTesting','flatTreeForTesting')
  forTree_leadPtOvMgg             = [-9999.]
  forTree_leadPtOvMgg             = array(forTree_leadPtOvMgg, dtype=float)
  forTree_subleadPtOvMgg          = [-9999.]
  forTree_subleadPtOvMgg          = array(forTree_subleadPtOvMgg, dtype=float)
  flatTree.Branch('leadPtOvMgg',forTree_leadPtOvMgg,'leadPtOvMgg/D')
  flatTree.Branch('subleadPtOvMgg',forTree_subleadPtOvMgg,'subleadPtOvMgg/D')
  forTree_leadEta              = [-9999.]
  forTree_leadEta = array(forTree_leadEta, dtype=float)
  forTree_subleadEta           = [-9999.]
  forTree_subleadEta = array(forTree_subleadEta, dtype=float)
  flatTree.Branch('leadEta',forTree_leadEta,'leadEta/D')
  flatTree.Branch('subleadEta',forTree_subleadEta,'subleadEta/D')
  forTree_cosDeltaPhi               = [-9999.]
  forTree_cosDeltaPhi = array(forTree_cosDeltaPhi, dtype=float)
  flatTree.Branch('cosDeltaPhi',forTree_cosDeltaPhi,'cosDeltaPhi/D')
  forTree_dijetDiphoDeltaPhi               = [-9999.]
  forTree_dijetDiphoDeltaPhi = array(forTree_dijetDiphoDeltaPhi, dtype=float)
  flatTree.Branch('dijetDiphoDeltaPhi',forTree_dijetDiphoDeltaPhi,'dijetDiphoDeltaPhi/D')
  forTree_leadJetPt               = [-9999.]
  forTree_leadJetPt = array(forTree_leadJetPt, dtype=float)
  forTree_subleadJetPt            = [-9999.]
  forTree_subleadJetPt = array(forTree_subleadJetPt, dtype=float)
  flatTree.Branch('leadJetPt',forTree_leadJetPt,'leadJetPt/D')
  flatTree.Branch('subleadJetPt',forTree_subleadJetPt,'subleadJetPt/D')
  forTree_jetDeltaEta               = [-9999.]
  forTree_jetDeltaEta = array(forTree_jetDeltaEta, dtype=float)
  flatTree.Branch('jetDeltaEta',forTree_jetDeltaEta,'jetDeltaEta/D')
  forTree_jetDeltaPhi               = [-9999.]
  forTree_jetDeltaPhi = array(forTree_jetDeltaPhi, dtype=float)
  flatTree.Branch('jetDeltaPhi',forTree_jetDeltaPhi,'jetDeltaPhi/D')
  forTree_leadJetPtD              = [-9999.]
  forTree_leadJetPtD = array(forTree_leadJetPtD, dtype=float)
  forTree_subleadJetPtD           = [-9999.]
  forTree_subleadJetPtD = array(forTree_subleadJetPtD, dtype=float)
  flatTree.Branch('leadJetPtD',forTree_leadJetPtD,'leadJetPtD/D')
  flatTree.Branch('subleadJetPtD',forTree_subleadJetPtD,'subleadJetPtD/D')
  forTree_leadJetChargedNConst    = [-9999]
  forTree_leadJetChargedNConst = array(forTree_leadJetChargedNConst, dtype=int)
  forTree_subleadJetChargedNConst = [-9999]
  forTree_subleadJetChargedNConst = array(forTree_subleadJetChargedNConst, dtype=int)
  flatTree.Branch('leadJetChargedNConst',forTree_leadJetChargedNConst,'leadJetChargedNConst/I')
  flatTree.Branch('subleadJetChargedNConst',forTree_subleadJetChargedNConst,'subleadJetChargedNConst/I')
  forTree_leadJetAxis2            = [-9999.]
  forTree_leadJetAxis2 = array(forTree_leadJetAxis2, dtype=float)
  forTree_subleadJetAxis2         = [-9999.]
  forTree_subleadJetAxis2 = array(forTree_subleadJetAxis2, dtype=float)
  flatTree.Branch('leadJetAxis2',forTree_leadJetAxis2,'leadJetAxis2/D')
  flatTree.Branch('subleadJetAxis2',forTree_subleadJetAxis2,'subleadJetAxis2/D')
  forTree_dijetMass               = [-9999.]
  forTree_dijetMass = array(forTree_dijetMass, dtype=float)
  flatTree.Branch('dijetMass',forTree_dijetMass,'dijetMass/D')

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
  nTreeEvents = 0
  for i,ev in enumerate(theTree):
    if i==0: print 'Processing first event'
    elif i%10000==0: print 'Processing event %g'%i
    if i==opts.maxEvents: 
      break
    if nTreeEvents==opts.maxTreeEvents: 
      print 'terminated on event %g, therefore total efficiency is %1.3f' % (i,float(nTreeEvents)/i)
      break
    #check two photons
    if not doLoose:
      nPhotons = getattr(ev,'PhotonTight_size')
    else:
      nPhotons = getattr(ev,'PhotonLoose_size')
    if nPhotons < 2:
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

    #recalculate four-vectors if reco vertex further than 1cm from gen vertex
    if not correctVtx: 
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

    #2D energy correction hist stuff can come first
    leadPhotonCorrection    = genPhotonE[leadGenIndex] / leadPhoton.E()
    subleadPhotonCorrection = genPhotonE[subleadGenIndex] / subleadPhoton.E()

    #then diphoton mass plots
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
      leadPhoton.SetE( leadPhoton.E() * leadCorrValue )
      subleadPhoton.SetE( subleadPhoton.E() * subleadCorrValue )
    theDiphoton = leadPhoton + subleadPhoton
    diphoMass = theDiphoton.M()
    if diphoMass < 100. or diphoMass > 180.:
      continue
    if not opts.forTesting:
      if not (4*leadPhoton.Pt()>diphoMass and 5*subleadPhoton.Pt()>diphoMass):
        continue
    else:
      if not (3*leadPhoton.Pt()>diphoMass and 4*subleadPhoton.Pt()>diphoMass):
        continue

    #find lead and sublead jets
    if nJets < 2:
      continue
    leadJetIndex = -1
    leadJetPt    = -1.
    subleadJetIndex = -1
    subleadJetPt    = -1.
    for iJet in range(nJets):
      #if not jetID[iJet] > 2:
      #if not jetID[iJet]:
      #  continue
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
      #if not jetID[iJet2] > 2:
      #if not jetID[iJet2]:
      #  continue
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
    if not opts.forTesting:
      if not (leadJetPt>30. and subleadJetPt>20.):
        continue
    else:
      if not (leadJetPt>40. and subleadJetPt>30.):
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
    #then project to VBF phase space
    if not opts.forTesting:
      if not (dijetMass > 100.):
        continue
    else:
      if not (dijetMass > 250.):
        continue
      if not (diphoMass > 123. and diphoMass < 127.):
        if 'GJet' in theKey or 'DiPho' in theKey:
          continue
    #print diphoMass
    #assign vars, write tree
    forTree_leadPtOvMgg[0]             = leadPhoton.Pt() / diphoMass
    forTree_subleadPtOvMgg[0]          = subleadPhoton.Pt() / diphoMass
    forTree_leadEta[0]                 = leadPhoton.Eta()
    forTree_subleadEta[0]              = subleadPhoton.Eta()
    forTree_cosDeltaPhi[0]             = cos( deltaPhi(leadPhoton.Phi(),subleadPhoton.Phi()) )
    forTree_dijetDiphoDeltaPhi[0]      = deltaPhi(theDiphoton.Phi(),theDijet.Phi())
    forTree_leadJetPt[0]               = leadJet.Pt()
    forTree_subleadJetPt[0]            = subleadJet.Pt()
    forTree_dijetMass[0]               = theDijet.M()
    forTree_jetDeltaEta[0]             = abs(leadJet.Eta()-subleadJet.Eta())
    forTree_jetDeltaPhi[0]             = deltaPhi(leadJet.Phi(),subleadJet.Phi())
    forTree_leadJetPtD[0]              = jetVals['ptD'][leadJetIndex]
    forTree_subleadJetPtD[0]           = jetVals['ptD'][subleadJetIndex]
    forTree_leadJetChargedNConst[0]    = jetVals['chargedNConst'][leadJetIndex]
    forTree_subleadJetChargedNConst[0] = jetVals['chargedNConst'][subleadJetIndex]
    forTree_leadJetAxis2[0]            = jetVals['Axis2'][leadJetIndex]
    forTree_subleadJetAxis2[0]         = jetVals['Axis2'][subleadJetIndex]
    flatTree.Fill()
    nTreeEvents += 1

  #post-processing
  flatFile.Write()
  flatFile.Close()
  corrFile.Close()


if __name__ == '__main__':
  main()
