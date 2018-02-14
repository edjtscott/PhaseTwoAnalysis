#!/usr/bin/env python
# code to loop over HGCAL photon ntuples and make plots of pT(H) in different eta regions

import os
from numpy import sin, cos, tan, sqrt
from collections import OrderedDict as odict
import ROOT as r
from diphoHelpers import initHistByEta, fillHistByEta, getEffSigma, printHists, initHistByJetType, fillHistByJetType, drawJetHist, deltaR

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-k', '--key', default='VBF_PU200', help='choose the sample to run on')
parser.add_option('-m', '--maxEvents', type='int', default=-1, help='specify number of events on which to run')
(opts,args) = parser.parse_args()

r.gROOT.SetBatch(True)


def main():
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
  #setup trees
  genPhotonTree = r.TChain('ntuple/GenPhoton')
  genVtxTree    = r.TChain('ntuple/GenVertex')
  recoVtxTree   = r.TChain('ntuple/Vertex')
  for fileNum in range(1,nFiles[theKey]+1):
    if fileNum in fileVetos[theKey]: continue
    fileName = 'root://eoscms.cern.ch//%s%g.root'%(baseNames[theKey],fileNum)
    genPhotonTree.Add(fileName)
    genVtxTree.Add(fileName)
    recoVtxTree.Add(fileName)

  print 'got trees from files'
  nEvents = genPhotonTree.GetEntries()
  print 'found %g events'%nEvents
  if opts.maxEvents > 0 and nEvents > opts.maxEvents:
    print 'but max is %g so will only use that many'%opts.maxEvents

  #setup the hists for different eta combos
  ptHists = odict()
  #pT(H) distributions
  ptHists['BarrelBarrel'] = r.TH1F('ptH_BB', 'ptH_BB', 50, 0., 250.)
  ptHists['BarrelEndcap'] = r.TH1F('ptH_EB', 'ptH_EB', 50, 0., 250.)
  ptHists['EndcapEndcap'] = r.TH1F('ptH_EE', 'ptH_EE', 50, 0., 250.)
  ptHists['SameEndcap']   = r.TH1F('ptH_sameE', 'ptH_sameE', 50, 0., 250.)
  #deltaEta distributions
  ptHists['deltaEta']   = r.TH1F('deltaEta', 'deltaEta', 30, 0., 3.)
  ptHists['dEtaLow']   = r.TH1F('dEtaLow', 'dEtaLow', 30, 0., 3.)
  ptHists['dEtaLowPtLow']   = r.TH1F('dEtaLowPtLow', 'dEtaLowPtLow', 30, 0., 3.)
  ptHists['dEtaLowPtHigh']   = r.TH1F('dEtaLowPtHigh', 'dEtaLowPtHigh', 30, 0., 3.)
  ptHists['dEtaHigh']   = r.TH1F('dEtaHigh', 'dEtaHigh', 30, 0., 3.)
  ptHists['dEtaHighPtLow']   = r.TH1F('dEtaHighPtLow', 'dEtaHighPtLow', 30, 0., 3.)
  ptHists['dEtaHighPtHigh']   = r.TH1F('dEtaHighPtHigh', 'dEtaHighPtHigh', 30, 0., 3.)
  ptHists['dEtaForLowPtH']   = r.TH1F('dEtaForLowPtH', 'dEtaForLowPtH', 30, 0., 3.)
  ptHists['dEtaForHighPtH']   = r.TH1F('dEtaForHighPtH', 'dEtaForHighPtH', 30, 0., 3.)
  #2D / using both
  ptHists['BB_dEtaLow'] = r.TH1F('BB_dEtaLow', 'BB_dEtaLow', 1, 0.5, 1.5)
  ptHists['BB_dEtaHigh'] = r.TH1F('BB_dEtaHigh', 'BB_dEtaHigh', 1, 0.5, 1.5)
  ptHists['EB_dEtaLow'] = r.TH1F('EB_dEtaLow', 'EB_dEtaLow', 1, 0.5, 1.5)
  ptHists['EB_dEtaHigh'] = r.TH1F('EB_dEtaHigh', 'EB_dEtaHigh', 1, 0.5, 1.5)
  ptHists['EE_dEtaLow'] = r.TH1F('EE_dEtaLow', 'EE_dEtaLow', 1, 0.5, 1.5)
  ptHists['EE_dEtaHigh'] = r.TH1F('EE_dEtaHigh', 'EE_dEtaHigh', 1, 0.5, 1.5)
  ptHists['PtVsDeltaEta']   = r.TH2F('PtVsDeltaEta', 'PtVsDeltaEta', 30, 0., 3., 25, 0., 250.)

  #main loop starts here
  for i,ev in enumerate(genPhotonTree):
    if i==0: print 'Processing first event'
    elif i%10000==0: print 'Processing event %g'%i
    if i==opts.maxEvents: break
    #check two photons
    #setup collections
    #gen and reco vertex z
    genVtxTree.GetEntry(i)
    genVtxZ  = getattr(genVtxTree,'Z')
    recoVtxTree.GetEntry(i)
    recoVtxZ = getattr(recoVtxTree,'Z')
    if len(recoVtxZ) < 1:
      recoVtxZ = recoVtxZ[0]
    else: recoVtxZ = 0.
    correctVtx = abs(genVtxZ-recoVtxZ) < 1.
    #get gen photon energies
    genPhotonTree.GetEntry(i)
    genPhotonE = getattr(genPhotonTree,'E')
    genPhotonPt = getattr(genPhotonTree,'PT')
    genPhotonEta = getattr(genPhotonTree,'Eta')
    genPhotonPhi = getattr(genPhotonTree,'Phi')
    nPhotons = len(genPhotonE)
    if nPhotons < 2:
      continue
    #get leading and subleading photon
    leadIndex = -1
    leadPt = -9999.
    subleadIndex = -1
    subleadPt = -9999.
    for iPho in range(len(genPhotonPt)):
      tempPt = genPhotonPt[iPho]
      if tempPt > leadPt: 
        leadPt = tempPt
        leadIndex = iPho
    for iPho in range(len(genPhotonPt)):
      if iPho == leadIndex: continue
      tempPt = genPhotonPt[iPho]
      if tempPt > subleadPt: 
        subleadPt = tempPt
        subleadIndex = iPho
    if leadPt < 40.: continue
    if subleadPt < 40.: continue
    leadPhoton = r.TLorentzVector()
    leadPhoton.SetPtEtaPhiM(genPhotonPt[leadIndex], genPhotonEta[leadIndex], genPhotonPhi[leadIndex], 0)
    leadEnergy = leadPhoton.Energy()
    subleadPhoton = r.TLorentzVector()
    subleadPhoton.SetPtEtaPhiM(genPhotonPt[subleadIndex], genPhotonEta[subleadIndex], genPhotonPhi[subleadIndex], 0)
    subleadEnergy = subleadPhoton.Energy()
    theAngle = leadPhoton.Angle(subleadPhoton.Vect())
    theDiphoton = leadPhoton + subleadPhoton
    diphoMass = theDiphoton.M()
    #print diphoMass
    leadIsBarrel = abs(leadPhoton.Eta()) < 1.44
    leadIsEndcap = abs(leadPhoton.Eta()) > 1.6 and abs(leadPhoton.Eta()) < 2.8
    subleadIsBarrel = abs(subleadPhoton.Eta()) < 1.44
    subleadIsEndcap = abs(subleadPhoton.Eta()) > 1.6 and abs(subleadPhoton.Eta()) < 2.8
    deltaEta = leadPhoton.Eta() - subleadPhoton.Eta()
    deltaEta = abs(deltaEta)
    if leadIsBarrel and subleadIsBarrel: 
      ptHists['BarrelBarrel'].Fill( theDiphoton.Pt() )
      if abs(deltaEta) < 0.8:   ptHists['BB_dEtaLow'].Fill(1)
      elif abs(deltaEta) > 0.8: ptHists['BB_dEtaHigh'].Fill(1)
    elif leadIsEndcap and subleadIsBarrel: 
      ptHists['BarrelEndcap'].Fill( theDiphoton.Pt() )
      if abs(deltaEta) < 0.8:   ptHists['EB_dEtaLow'].Fill(1)
      elif abs(deltaEta) > 0.8: ptHists['EB_dEtaHigh'].Fill(1)
    elif leadIsBarrel and subleadIsEndcap: 
      ptHists['BarrelEndcap'].Fill( theDiphoton.Pt() )
      if abs(deltaEta) < 0.8:   ptHists['EB_dEtaLow'].Fill(1)
      elif abs(deltaEta) > 0.8: ptHists['EB_dEtaHigh'].Fill(1)
    elif leadIsEndcap and subleadIsEndcap: 
      ptHists['EndcapEndcap'].Fill( theDiphoton.Pt() )
      if abs(deltaEta) < 0.8:   ptHists['EE_dEtaLow'].Fill(1)
      elif abs(deltaEta) > 0.8: ptHists['EE_dEtaHigh'].Fill(1)
      if leadPhoton.Eta() * subleadPhoton.Eta() > 0.: 
        ptHists['SameEndcap'].Fill( theDiphoton.Pt() )
    ptHists['deltaEta'].Fill( deltaEta )
    if abs(deltaEta) < 0.8: 
      ptHists['dEtaLow'].Fill( deltaEta )
      if abs(theDiphoton.Pt()) < 50:  ptHists['dEtaLowPtLow'].Fill( deltaEta )
      elif abs(theDiphoton.Pt()) > 50:  ptHists['dEtaLowPtHigh'].Fill( deltaEta )
    elif abs(deltaEta) > 0.8: 
      ptHists['dEtaHigh'].Fill( deltaEta )
      if abs(theDiphoton.Pt()) < 50:  ptHists['dEtaHighPtLow'].Fill( deltaEta )
      elif abs(theDiphoton.Pt()) > 50:  ptHists['dEtaHighPtHigh'].Fill( deltaEta )
    if theDiphoton.Pt() < 50.: 
      ptHists['dEtaForLowPtH'].Fill( deltaEta )
    elif theDiphoton.Pt() > 50.: 
      ptHists['dEtaForHighPtH'].Fill( deltaEta )
    ptHists['deltaEta'].Fill( deltaEta )
    ptHists['PtVsDeltaEta'].Fill( deltaEta, theDiphoton.Pt() )

  #end of event loop
  print 'Processed last event'

  #print fractions etc
  nBB = ptHists['BarrelBarrel'].GetEntries()
  nEB = ptHists['BarrelEndcap'].GetEntries()
  nEE = ptHists['EndcapEndcap'].GetEntries()
  nSameE = ptHists['SameEndcap'].GetEntries()
  nLow = ptHists['dEtaLow'].GetEntries()
  nHigh = ptHists['dEtaHigh'].GetEntries()
  nLowLow = ptHists['dEtaLowPtLow'].GetEntries()
  nLowHigh = ptHists['dEtaLowPtHigh'].GetEntries()
  nHighLow = ptHists['dEtaHighPtLow'].GetEntries()
  nHighHigh= ptHists['dEtaHighPtHigh'].GetEntries()
  nBBlow  = ptHists['BB_dEtaLow'].GetEntries()
  nBBhigh = ptHists['BB_dEtaHigh'].GetEntries()
  nEBlow  = ptHists['EB_dEtaLow'].GetEntries()
  nEBhigh = ptHists['EB_dEtaHigh'].GetEntries()
  nEElow  = ptHists['EE_dEtaLow'].GetEntries()
  nEEhigh = ptHists['EE_dEtaHigh'].GetEntries()
  nTotPtLow = float(nLowLow + nHighLow)
  nTotPtHigh = float(nLowHigh+ nHighHigh)
  nTotEta = float(nLow + nHigh)
  nTot = float(nBB + nEB + nEE)
  print 'fraction of events in BB, EB and EE: %2.1f%%, %2.1f%%,, %2.1f%%,'%(100*nBB/nTot, 100*nEB/nTot, 100*nEE/nTot)
  print 'and fraction of events SAME endcap: %2.1f%%'%(100*nSameE/nTot)
  print 'fraction of low, high delta eta is %2.1f%%, %2.1f%%'%(100*nLow/nTotEta,100*nHigh/nTotEta)
  print 'fraction of low, high delta eta for ptH < 50 is %2.1f%%, %2.1f%%'%(100*nLowLow/nTotPtLow,100*nHighLow/nTotPtLow)
  print 'fraction of low, high delta eta for ptH > 50 is %2.1f%%, %2.1f%%'%(100*nLowHigh/nTotPtHigh,100*nHighHigh/nTotPtHigh)
  print 'fraction of low, high delta eta for BB events is %2.1f%%, %2.1f%%'%(100.*nBBlow/nBB,100.*nBBhigh/nBB)
  print 'fraction of low, high delta eta for EB events is %2.1f%%, %2.1f%%'%(100.*nEBlow/nEB,100.*nEBhigh/nEB)
  print 'fraction of low, high delta eta for EE events is %2.1f%%, %2.1f%%'%(100.*nEElow/nEE,100.*nEEhigh/nEE)
  print 'with total events passing selection being %g'%(int(nTot))
  #write out mass histos
  outFile = r.TFile('PtHiggsHists_%s.root'%theKey,'RECREATE')
  for key,hist in ptHists.iteritems():
    hist.Write()
  outFile.Close()
  canv = r.TCanvas()
  for key,hist in ptHists.iteritems():
    if 'Vs' in hist.GetName(): hist.Draw('colz')
    else: hist.Draw()
    canv.Print('%s_%s.pdf'%(hist.GetName(),theKey))
    canv.Print('%s_%s.png'%(hist.GetName(),theKey))

if __name__ == '__main__':
  main()
