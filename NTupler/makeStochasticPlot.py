#!/usr/bin/env python
# code to loop over HGCAL photon ntuples and make plot of mass resolution vs stochastic term

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

  #setup the stochastic values, and their gaussians and histograms
  # the sigma depends on energy, so have to set that in each event
  stoVals = odict()
  massHists = odict()
  for i in range(0,21):
    stoVals[i] = 0.2 + 0.01*i
    massHists[i] = r.TH1F('mass_%g'%i, 'mass_%g'%i, 40, 120, 130)
  #get random number generator to smear out the gen energies
  #according to formula sigma/E = s/sqrt(E) (+) 1%, for s between 0.2 and 0.4
  theGaus = r.TF1('smearGaus','gaus',0.,2.)

  #for debugging
  leadEnergySum = 0.
  subleadEnergySum = 0.
  photonCounter = 0
  leadEnergySumEE = 0.
  subleadEnergySumEE = 0.
  photonCounterEE = 0

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
      if abs(genPhotonEta[iPho]) > 3.: continue
      tempPt = genPhotonPt[iPho]
      if tempPt > leadPt: 
        leadPt = tempPt
        leadIndex = iPho
    for iPho in range(len(genPhotonPt)):
      if iPho == leadIndex: continue
      if abs(genPhotonEta[iPho]) > 3.: continue
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
    leadIsBarrel = abs(leadPhoton.Eta()) < 1.44
    leadIsEndcap = abs(leadPhoton.Eta()) > 1.6 and abs(leadPhoton.Eta()) < 2.8
    subleadIsBarrel = abs(subleadPhoton.Eta()) < 1.44
    subleadIsEndcap = abs(subleadPhoton.Eta()) > 1.6 and abs(subleadPhoton.Eta()) < 2.8
    eventIsBB = leadIsBarrel and subleadIsBarrel
    eventIsEB = (leadIsBarrel and subleadIsEndcap) or (leadIsEndcap and subleadIsBarrel)
    eventIsEE = leadIsEndcap and subleadIsEndcap
    #FIXME why does this have such a huge effect?! 
    #if not eventIsEE: continue
    leadEnergySum += leadEnergy
    subleadEnergySum += subleadEnergy
    photonCounter += 1
    if eventIsEE: 
      leadEnergySumEE += leadEnergy
      subleadEnergySumEE += subleadEnergy
      photonCounterEE += 1
    theAngle = leadPhoton.Angle(subleadPhoton.Vect())
    theDiphoton = leadPhoton + subleadPhoton
    diphoMass = theDiphoton.M()

    #print diphoMass
    #print sqrt( 2 * leadPhoton.Energy() * subleadPhoton.Energy() * (1. - cos(leadPhoton.Angle(subleadPhoton.Vect()))))
    #TODO: get random number generator to smear out the gen energies
    #according to formula sigma/E = s/sqrt(E) (+) 1%, for s between 0.2 and 0.4
    for key,val in stoVals.items(): 
      #leadSigma = sqrt( (val/sqrt(leadEnergy))**2 + 0.01*0.01 )
      leadSigma = sqrt( (val/sqrt(leadEnergy))**2 + 0.008*0.008 )
      #FIXME sigmas set here
      theGaus.SetParameters(1., 1., leadSigma)
      tempLeadEnergy = leadEnergy * theGaus.GetRandom()
      #subleadSigma = sqrt( (val/sqrt(subleadEnergy))**2 + 0.01*0.01 )
      subleadSigma = sqrt( (val/sqrt(subleadEnergy))**2 + 0.008*0.008 )
      theGaus.SetParameters(1., 1., subleadSigma)
      tempSubleadEnergy = subleadEnergy * theGaus.GetRandom()
      tempMass = sqrt( 2 * tempLeadEnergy * tempSubleadEnergy * (1. - cos(theAngle)))
      massHists[key].Fill( tempMass )

  #end of event loop
  print 'Processed last event'

  #fit the hists
  canv = r.TCanvas()
  sigmaFile = open('massHistSigmas.txt','w')
  theVals = []
  for key,hist in massHists.items():
    hist.Draw()
    tempFitter = r.TF1('tempFit','gaus',120.,130.)
    hist.Fit(tempFitter)
    massVal = tempFitter.GetParameter(1)
    sigmaVal = tempFitter.GetParameter(2)
    sigmaErr = tempFitter.GetParError(2)
    stochTerm = 20 + key
    theVals.append( (stochTerm,sigmaVal/125.) )
    doPretty = False
    if doPretty:
      sigmaFile.write('Sigma/M for s = %g%% is %1.3f%%, error %1.3f%%\n'%(stochTerm, 100.*sigmaVal/125., 100.*sigmaErr/125.))
    else :
      sigmaFile.write('%g,%1.3f,%1.3f\n'%(stochTerm, 100.*sigmaVal/125., 100.*sigmaErr/125.))
  sigmaFile.close()
  canv.Clear()

  #now (attempt to) make a pretty plot... 
  #theGraph = r.TH1F('resoVsStoch','',21,19.5,40.5)
  #for pair in theVals:
  #  theGraph.Fill(pair[0], 100*pair[1])
  #theGraph.GetXaxis().SetTitle('Stochastic term (%)')
  #theGraph.GetYaxis().SetTitle('Mass resolution (%)')
  #theGraph.SetMarkerStyle(33)
  #theGraph.SetMinimum(1.4)
  #theGraph.SetMaximum(2.8)
  #theGraph.SetStats(0)
  #lat = r.TLatex()
  #lat.SetTextFont(42)
  #lat.SetTextSize(0.025)
  #lat.SetLineWidth(2)
  #lat.SetTextSize(0.05)
  #lat.SetTextAlign(11)
  #lat.DrawLatex(0.05,0.95,"#bf{CMS} #it{Simulation}     H#rightarrow#gamma#gamma")
  #theGraph.Draw('P')
  ##theGraph.Draw('hist,c')
  #canv.Print('MassVsStochasticTerm.pdf')
  #canv.Print('MassVsStochasticTerm.png')

  ##write out mass histos
  #outFile = r.TFile('StochasticHists.root','RECREATE')
  #for key,hist in massHists.items():
  #  hist.Write()
  #outFile.Close()

  print 'inclusive mean lead energy is %3.1f'%(leadEnergySum/photonCounter)
  print 'endcap mean lead energy is    %3.1f'%(leadEnergySumEE/photonCounterEE)
  print 'inclusive mean sublead energy is %3.1f'%(subleadEnergySum/photonCounter)
  print 'endcap mean sublead energy is    %3.1f'%(subleadEnergySumEE/photonCounterEE)

if __name__ == '__main__':
  main()
