#!/usr/bin/env python
# code to loop over HGCAL photon ntuples and make some plots

import os
import sys
import ROOT as r
from diphoHelpers import initHistByEta, fillHistByEta, getEffSigma, printHists

r.gROOT.SetBatch(True)

def main():
  #setup histos with several eta breakdown scenarios
  etaHists = {}
  initHistByEta(etaHists, 'mggTightID', 40, 115, 135)
  initHistByEta(etaHists, 'mggTightIDMulti', 40, 115, 135)
  initHistByEta(etaHists, 'mggPtCutsTightID', 40, 115, 135)
  initHistByEta(etaHists, 'mggPtCutsTightIDMulti', 40, 115, 135)
  
  #PUname = 'PU0'
  PUname = 'PU200_v2'
  fileName = 'new_VBF_%s.root'%PUname
  webDir = '/afs/cern.ch/user/e/escott/www/HFuture/Pass1/VBF_%s'%PUname
  
  theFile = r.TFile(fileName)
  print 'got file %s'%fileName
  theTree = theFile.Get('ntuple/PhotonTight')
  print 'got tree ntuple/PhotonTight from file %s'%fileName
  for i,ev in enumerate(theTree):
    if i==0: print 'Processing first event'
    elif i%1000==0: print 'Processing event %g'%i
    #check two photons
    nPhotons = getattr(ev,'PhotonTight_size')
    if not nPhotons==2: 
      continue
    #setup collections
    #gen and reco vertex z
    #genVtxZ  = getattr(ev,'')
    #recoVtxZ = getattr(ev,'')
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

    #process normal values
    leadIndex    = int(photonPt[1] > photonPt[0])
    subleadIndex = 1 - leadIndex
    leadPhoton = r.TLorentzVector()
    leadPhoton.SetPtEtaPhiM(photonPt[leadIndex], photonEta[leadIndex], photonPhi[leadIndex], 0)
    subleadPhoton = r.TLorentzVector()
    subleadPhoton.SetPtEtaPhiM(photonPt[subleadIndex], photonEta[subleadIndex], photonPhi[subleadIndex], 0)
    theDiphoton = leadPhoton+subleadPhoton
    diphoMass = theDiphoton.M()
    fillHistByEta(etaHists, 'mggTightID', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)
    if 3*leadPhoton.Pt()>diphoMass and 4*subleadPhoton.Pt()>diphoMass:
      fillHistByEta(etaHists, 'mggPtCutsTightID', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)

    #again for multi values
    leadIndex    = int(photonPtMulti[1] > photonPtMulti[0])
    subleadIndex = 1 - leadIndex
    leadPhoton = r.TLorentzVector()
    leadPhoton.SetPtEtaPhiM(photonPtMulti[leadIndex], photonEtaMulti[leadIndex], photonPhiMulti[leadIndex], 0)
    subleadPhoton = r.TLorentzVector()
    subleadPhoton.SetPtEtaPhiM(photonPtMulti[subleadIndex], photonEtaMulti[subleadIndex], photonPhiMulti[subleadIndex], 0)
    theDiphoton = leadPhoton+subleadPhoton
    diphoMass = theDiphoton.M()
    fillHistByEta(etaHists, 'mggTightIDMulti', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)
    if 3*leadPhoton.Pt()>diphoMass and 4*subleadPhoton.Pt()>diphoMass:
      fillHistByEta(etaHists, 'mggPtCutsTightIDMulti', leadPhoton.Eta(), subleadPhoton.Eta(), diphoMass)

  canv = r.TCanvas('canv','canv')
  outdirName = 'DiphoPlots_VBF_%s/'%PUname
  os.system('mkdir -p %s'%outdirName)
  printHists(canv, etaHists, outdirName)
  os.system('cp %s* %s'%(outdirName,webDir))
  print 'plots moved to %s'%webDir


if __name__ == '__main__':
  main()
