#!/usr/bin/env python
# code to loop over HGCAL photon ntuples and make some plots

import os
import sys
import ROOT as r
from diphoHelpers import initHistByEta, fillHistByEta, getEffSigma, printHists

r.gROOT.SetBatch(True)

#const values for MVA ID
barrelLooseWP = 0.
barrelTightWP = 0.56
endcapLooseWP = 0.2
endcapTightWP = 0.68


def main():
  #setup histos
  etaHists = {}
  initHistByEta(etaHists, 'mggAllPrompt', 40, 115, 135)
  initHistByEta(etaHists, 'mggLooseID', 40, 115, 135)
  initHistByEta(etaHists, 'mggTightID', 40, 115, 135)
  initHistByEta(etaHists, 'mggPtCuts', 40, 115, 135)
  initHistByEta(etaHists, 'mggPtCutsTightID', 40, 115, 135)
  
  #fileName = '/vols/cms/es811/HFuture/200PU-06A536FA-E5B3-E711-B827-24BE05CE1E51.root'
  fileName = '/vols/cms/es811/HFuture/tempPU200.root'
  #fileName = '/vols/cms/es811/HFuture/0PU-4files.root'
  #webDir = '/afs/cern.ch/user/e/escott/www/FinalFits/HFuture/VBF_PU200' 
  
  f = r.TFile(fileName)
  print 'got file %s'%fileName,f
  t = f.Get('ntupler/photons')
  print 'got tree ntupler/photons',f
  for i,ev in enumerate(t):
    if i==0: print 'Processing first event'
    elif i%1000==0: print 'Processing event %g'%i
    #setup collections
    gen_pt  = getattr(ev,'gen_pt')
    gen_eta = getattr(ev,'gen_eta')
    gen_phi = getattr(ev,'gen_phi')
    gen_conversionZ = getattr(ev,'gen_conversionZ')
    gen_isPromptFinalState = getattr(ev,'gen_isPromptFinalState')
    gen_iLocalReco = getattr(ev,'gen_iLocalReco')
    gen_iGedReco = getattr(ev,'gen_iGedReco')

    gedReco_pt  = getattr(ev,'gedReco_pt')
    gedReco_eta = getattr(ev,'gedReco_eta')
    gedReco_phi = getattr(ev,'gedReco_phi')
    gedReco_mvaValue = getattr(ev,'gedReco_mvaValue')
    localReco_pt  = getattr(ev,'localReco_pt')
    localReco_eta = getattr(ev,'localReco_eta')
    localReco_phi = getattr(ev,'localReco_phi')
    localReco_mvaValue = getattr(ev,'localReco_mvaValue')
  
    #loop over gen photons, make prompt diphotons
    nGenPhotons = len(gen_pt)
    #print 'nGenPhotons=%g'%nGenPhotons
    for iGen in range(nGenPhotons):
      genIsPrompt = gen_isPromptFinalState[iGen]
      if not genIsPrompt: continue
      genPt  = gen_pt[iGen]
      genEta = gen_eta[iGen]
      genPhi = gen_phi[iGen]
      #print 'gen photon 1 pt, eta, phi:   %1.2f, %1.2f, %1.2f'%(genPt,genEta,genPhi)
      genGedReco = gen_iGedReco[iGen]
      genLocalReco = gen_iLocalReco[iGen]
      recoPt  = -999.
      recoEta = -999.
      recoPhi = -999.
      if abs(genEta) < 1.5 and genGedReco>-1:
        recoPt  = gedReco_pt[genGedReco]
        recoEta = gedReco_eta[genGedReco]
        recoPhi = gedReco_phi[genGedReco]
        recoMVA = gedReco_mvaValue[genGedReco]
        recoIsEB = True
      elif abs(genEta) < 3.0 and genLocalReco>-1:
        recoPt  = localReco_pt[genLocalReco]
        recoEta = localReco_eta[genLocalReco]
        recoPhi = localReco_phi[genLocalReco]
        recoMVA = localReco_mvaValue[genLocalReco]
        recoIsEB = False
      else: continue
      recoP4 = r.TLorentzVector()
      recoP4.SetPtEtaPhiM(recoPt, recoEta, recoPhi, 0)
      #print 'reco photon 1 pt, eta, phi:   %1.2f, %1.2f, %1.2f'%(recoP4.Pt(), recoP4.Eta(), recoP4.Phi())
  
      #loop over possible pairings
      for iGen2 in range(iGen+1,nGenPhotons):
        gen2IsPrompt = gen_isPromptFinalState[iGen2]
        if not gen2IsPrompt: continue
        gen2Pt  = gen_pt[iGen2]
        gen2Eta = gen_eta[iGen2]
        gen2Phi = gen_phi[iGen2]
        #print 'gen photon 2 pt, eta, phi:   %1.2f, %1.2f, %1.2f'%(gen2Pt,gen2Eta,gen2Phi)
        gen2GedReco = gen_iGedReco[iGen2]
        gen2LocalReco = gen_iLocalReco[iGen2]
        reco2Pt  = -999.
        reco2Eta = -999.
        reco2Phi = -999.
        reco2IsEB = None
        if abs(gen2Eta) < 1.5 and gen2GedReco>-1:
          reco2Pt  = gedReco_pt[gen2GedReco]
          reco2Eta = gedReco_eta[gen2GedReco]
          reco2Phi = gedReco_phi[gen2GedReco]
          reco2MVA = gedReco_mvaValue[gen2GedReco]
          reco2IsEB = True
        elif abs(gen2Eta) < 3.0 and gen2LocalReco>-1:
          reco2Pt  = localReco_pt[gen2LocalReco]
          reco2Eta = localReco_eta[gen2LocalReco]
          reco2Phi = localReco_phi[gen2LocalReco]
          reco2MVA = localReco_mvaValue[gen2LocalReco]
          reco2IsEB = False
        else: continue
        reco2P4 = r.TLorentzVector()
        reco2P4.SetPtEtaPhiM(reco2Pt, reco2Eta, reco2Phi, 0)
        #print 'reco photon 2 pt, eta, phi:   %1.2f, %1.2f, %1.2f'%(recoP4.Pt(), recoP4.Eta(), recoP4.Phi())
        diphoP4 = recoP4+reco2P4
        diphoMass = diphoP4.M()
        #print 'just created diphoton with mass %1.2f GeV'%diphoP4.M()
        fillHistByEta(etaHists, 'mggAllPrompt', recoP4.Eta(), reco2P4.Eta(), diphoMass)
        recoPassLoose = (recoIsEB and recoMVA>barrelLooseWP) or (not recoIsEB and recoMVA>endcapLooseWP)
        recoPassTight = (recoIsEB and recoMVA>barrelTightWP) or (not recoIsEB and recoMVA>endcapTightWP)
        reco2PassLoose = (reco2IsEB and reco2MVA>barrelLooseWP) or (not reco2IsEB and reco2MVA>endcapLooseWP)
        reco2PassTight = (reco2IsEB and reco2MVA>barrelTightWP) or (not reco2IsEB and reco2MVA>endcapTightWP)
        if recoPassLoose and reco2PassLoose:
          fillHistByEta(etaHists, 'mggLooseID', recoP4.Eta(), reco2P4.Eta(), diphoMass)
        if recoPassTight and reco2PassTight:
          fillHistByEta(etaHists, 'mggTightID', recoP4.Eta(), reco2P4.Eta(), diphoMass)
        leadPt = max(recoP4.Pt(),reco2P4.Pt())
        subPt  = min(recoP4.Pt(),reco2P4.Pt())
        if 3*leadPt>diphoMass and 4*subPt>diphoMass:
          fillHistByEta(etaHists, 'mggPtCuts', recoP4.Eta(), reco2P4.Eta(), diphoMass)
          if recoPassTight and reco2PassTight:
            fillHistByEta(etaHists, 'mggPtCutsTightID', recoP4.Eta(), reco2P4.Eta(), diphoMass)
        
      #end second loop over gen photons
    #end first loop over gen photons
  canv = r.TCanvas('canv','canv')
  outdirName = 'DiphoPlots_%s/'%fileName.split('/')[-1].split('.root')[0]
  os.system('mkdir -p %s'%outdirName)
  printHists(canv, etaHists, outdirName)


if __name__ == '__main__':
  main()
