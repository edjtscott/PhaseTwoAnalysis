#!/usr/bin/env python
# takes the corrction hists, sum them up to get overall correction

import os
import ROOT as r

r.gROOT.SetBatch(True)

doLoose = False
#doLoose = True

def main():
  combinedECorrection        = r.TH2F('combinedECorrection',       'combinedECorrection',        6, 0., 3., 25, 30., 530.)
  combinedECorrectionEntries = r.TH2F('combinedECorrectionEntries','combinedECorrectionEntries', 6, 0., 3., 25, 30., 530.)

  theKeys = ['ggH_PU0', 'ggH_PU200', 'VBF_PU0', 'VBF_PU200']

  entriesMap = {}
  histMap    = {}
  inFiles = {}
  for key in theKeys:
    if not doLoose: inFiles[key] = r.TFile('CorrectionHists/corrHist_%s.root'%key, 'READ')
    else: inFiles[key] = r.TFile('CorrectionHistsLoose/corrHist_%s.root'%key, 'READ')
    tempCorrectionHist = inFiles[key].Get('photonECorrection')
    histMap[key] = tempCorrectionHist.Clone('photonECorrection_%s'%key)
    tempCorrectionHistEntries = inFiles[key].Get('photonECorrectionEntries')
    entriesMap[key] = tempCorrectionHistEntries.Clone('photonECorrectionEntries_%s'%key)
    combinedECorrectionEntries.Add(entriesMap[key].Clone('contrib_%s'%key))

  for key in theKeys:
    entriesMap[key].Divide(combinedECorrectionEntries.Clone('denom_%s'%key))
    histMap[key].Multiply(entriesMap[key].Clone('mult_%s'%key))
    combinedECorrection.Add(histMap[key].Clone('add_%s'%key))


  combinedECorrection.SetMinimum(0.)
  combinedECorrection.SetMaximum(3.)

  for key in theKeys: inFiles[key].Close()

  #draw hists, send to web
  canv = r.TCanvas('canv','canv')
  outdirName = 'DiphoPlots_Combined/'
  if doLoose: outdirName = 'DiphoPlots_Combined_Loose/'
  os.system('mkdir -p %s'%outdirName)
  webDir = '/afs/cern.ch/user/e/escott/www/HFuture/Pass1/Combined'
  if doLoose: webDir = '/afs/cern.ch/user/e/escott/www/HFuture/Pass1/Loose/Combined'
  os.system('mkdir -p %s'%webDir)
  os.system('cp /afs/cern.ch/user/e/escott/www/HFuture/Pass1/index.php %s'%webDir)
  #canv.cd()
  combinedECorrection.Draw('colz,text')
  canv.Print(outdirName+combinedECorrection.GetName()+'.pdf')
  canv.Print(outdirName+combinedECorrection.GetName()+'.png')
  os.system('cp %s* %s'%(outdirName,webDir))
  print 'plots moved to %s'%webDir
  #save the combined correction hist
  outFile = r.TFile('CorrectionHists/combinedCorrectionHist.root','RECREATE')
  if doLoose: outFile = r.TFile('CorrectionHistsLoose/combinedCorrectionHist.root','RECREATE')
  #outFile.cd()
  combinedECorrection.Write()
  outFile.Close()


if __name__ == '__main__':
  main()
