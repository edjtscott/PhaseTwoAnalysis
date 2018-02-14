#!/usr/bin/env python
# code to loop over HGCAL jet ntuples and make some plots

import os
import sys
import ROOT as r
from diphoHelpers import initHistByEta, fillHistByEta, getEffSigma, printHists

r.gROOT.SetBatch(True)

variables = ['chargedSumConst', 'neutralSumConst', 'hfemSumConst', 'hfhadSumConst', 'chargedNConst', 'neutralNConst', 'hfemNConst', 'hfhadNConst', 'eSumConst', 'eNConst', 'muSumConst', 'muNConst', 'photonSumConst', 'photonNConst', 'RMSCand', 'Axis1', 'Axis2', 'Sigma', 'ptD']

jetcats = [ "PU","VBF","Rad" ]

maxes = [50,50,50,50,20,20,20,20,50,20,50,20,50,20,0.4,0.4,0.4,0.4,1.0]

print len(variables),len(maxes)

bins = []
for maxval in maxes:
  if maxval == 20 or maxval == 50:
    bins.append(maxval)
  else:
    bins.append(40)

h = {}
b = {}

def main():
  #setup histos
  
  for i in range(len(variables)):
    var = variables[i]
    hbin = bins[i]
    hmax = maxes[i]
    for jetcat in jetcats:
      hname = '%s%s' % (var,jetcat)
      if hmax == 20:
        h[hname] =  r.TH1F(hname,hname,int(hbin), -0.5, float(hmax)-0.5)
      else:
        h[hname] =  r.TH1F(hname,hname,int(hbin), 0.0, float(hmax))
      h[hname].SetLineWidth(2)

#  fileName = "MiniTest.root"
  fileName = "/afs/cern.ch/work/e/escott/public/FinalFits/ForSeth/new_VBF_PU0.root"
#  fileName = "/afs/cern.ch/work/e/escott/public/FinalFits/ForSeth/new_VBF_PU200.root"

# PU = GenJet < 0
# radiation = GenJet >= 0  && (GenPartonPID == 0 ||  abs(GenPartonPID) > 6)
# vbfquark = GenJet >= 0 && GenPartonPID != 0 && abs(GenPartonPID) <= 6
 
  f = r.TFile(fileName)
  print 'got file %s'%fileName,f
  t = f.Get('ntuple/JetPUPPI')
  print 'got tree ntuple/JetPUPPI',t
  for i,ev in enumerate(t):
    if i==0: print 'Processing first event'
    elif i%1000==0: print 'Processing event %g'%i
    #setup collections
    PT = getattr(ev,'PT')
    Eta = getattr(ev,'Eta')
    GenJet = getattr(ev,'GenJet')
    GenPartonPID = getattr(ev,'GenPartonPID')
    ID = getattr(ev,"ID")
    for var in variables:
      b[var] = getattr(ev,var)
    for j in range(len(PT)):
      pt = PT[j]
      eta = Eta[j]
      gen = GenJet[j]
      parton = GenPartonPID[j]
      jetid = ID[j]
      if gen < 0: 
        strlabel = "PU"
      elif parton == 0 or abs(parton) > 6:
        strlabel = "Rad"
      else:
        strlabel = "VBF"
      print j,pt,eta,gen,parton,strlabel,jetid
      if pt < 30 or abs(eta) < 1.5 or abs(eta) > 3.0 or not jetid:
        continue
      for var in variables:
        h["%s%s" % (var,strlabel)].Fill(b[var][j])
        
      #end second loop over gen photons
    #end first loop over gen photons
  canv = r.TCanvas('canv','canv')
#  outdirName = 'JetPlots_%s/'%fileName.split('/')[-1].split('.root')[0]
  outdirName = "~/www/cms/JetPlots_15Nov_PU0"
  os.system('mkdir -p %s'%outdirName)
  for var in variables:
    h["%s%s" % (var,"VBF")].SetLineColor(1)
    h["%s%s" % (var,"Rad")].SetLineColor(2)
    h["%s%s" % (var,"PU")].SetLineColor(4)
    h["%s%s" % (var,"VBF")].Draw()
    h["%s%s" % (var,"Rad")].Draw("same")
    h["%s%s" % (var,"PU")].Draw("same")
    vpu = h["%s%s" % (var,"PU")].GetMaximum()
    vvbf = h["%s%s" % (var,"VBF")].GetMaximum()
    vrad = h["%s%s" % (var,"Rad")].GetMaximum()
    h["%s%s" % (var,"VBF")].SetMaximum(1.1*max(vpu,vvbf,vrad))
    h["%s%s" % (var,"VBF")].Draw()
    h["%s%s" % (var,"Rad")].Draw("same")
    h["%s%s" % (var,"PU")].Draw("same")
    for ext in ["png","pdf"]:
      canv.SaveAs("%s/%s.%s" % (outdirName,var,ext))

if __name__ == '__main__':
  main()
