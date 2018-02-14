import ROOT as r
from numpy import sqrt

#helper functions
def initHistByEta(histColl, histName, bins=50, low=1., high=0.):
  for etaName in ['_inclusive','_BB','_B1','_B2','_B3','_EB','_EE','_E1','_E2','_E3']:
    name =  histName + etaName
    histColl[name] = r.TH1F(name,name,bins,low,high)

def fillHistByEta(histColl, histName, eta1, eta2, val):
  minEta = abs(min(abs(eta1),abs(eta2)))
  maxEta = abs(max(abs(eta1),abs(eta2)))
  histColl[histName+'_inclusive'].Fill(val)
  if minEta<1.5:
    if maxEta<1.5: histColl[histName+'_BB'].Fill(val)
    elif maxEta<2.0: histColl[histName+'_B1'].Fill(val)
    elif maxEta<2.5: histColl[histName+'_B2'].Fill(val)
    elif maxEta<3.0: histColl[histName+'_B3'].Fill(val)
    if maxEta>1.5: histColl[histName+'_EB'].Fill(val)
  elif minEta<3.0:
    histColl[histName+'_EE'].Fill(val)
    if maxEta<2.0: histColl[histName+'_E1'].Fill(val)
    elif maxEta<2.5: histColl[histName+'_E2'].Fill(val)
    elif maxEta<3.0: histColl[histName+'_E3'].Fill(val)

def initHistByJetType(histColl, histName, bins=50, low=1., high=0.):
  for order in ['_all', '_lead', '_sub']:
    for jetName in ['_Quark','_PU','_NoParton','_OtherParton']:
      name =  histName + order + jetName
      if 'chargedNConst' in name:
        histColl[name] = r.TH1F(name,name,26,-0.5,25.5)
      else:
        histColl[name] = r.TH1F(name,name,bins,low,high)

def fillHistByJetType(histColl, histName, genIndex, partonID, order, val):
  if genIndex < 0:
    histColl[histName+'_all_PU'].Fill(val)
    if order==0: histColl[histName+'_lead_PU'].Fill(val)
    if order==1: histColl[histName+'_sub_PU'].Fill(val)
  elif partonID==0:
    histColl[histName+'_all_NoParton'].Fill(val)
    if order==0: histColl[histName+'_lead_NoParton'].Fill(val)
    if order==1: histColl[histName+'_sub_NoParton'].Fill(val)
  elif abs(partonID)>6:
    histColl[histName+'_all_OtherParton'].Fill(val)
    if order==0: histColl[histName+'_lead_OtherParton'].Fill(val)
    if order==1: histColl[histName+'_sub_OtherParton'].Fill(val)
  else:
    histColl[histName+'_all_Quark'].Fill(val)
    if order==0: histColl[histName+'_lead_Quark'].Fill(val)
    if order==1: histColl[histName+'_sub_Quark'].Fill(val)

def drawJetHist(canv, histColl, jetVar, outdir):
  canv.cd() 
  for order in ['_all', '_lead', '_sub']:
    axisMax = 1.1 * max( [ histColl[jetVar+order+'_Quark'].GetMaximum(), histColl[jetVar+order+'_NoParton'].GetMaximum(), histColl[jetVar+order+'_OtherParton'].GetMaximum(), histColl[jetVar+order+'_PU'].GetMaximum()] )
    histColl[jetVar+order+'_Quark'].SetLineColor(1)
    histColl[jetVar+order+'_OtherParton'].SetLineColor(2)
    histColl[jetVar+order+'_NoParton'].SetLineColor(3)
    histColl[jetVar+order+'_PU'].SetLineColor(4)
    histColl[jetVar+order+'_Quark'].SetMaximum(axisMax)
    histColl[jetVar+order+'_OtherParton'].SetMaximum(axisMax)
    histColl[jetVar+order+'_NoParton'].SetMaximum(axisMax)
    histColl[jetVar+order+'_PU'].SetMaximum(axisMax)
    histColl[jetVar+order+'_Quark'].Draw()
    histColl[jetVar+order+'_OtherParton'].Draw('same')
    histColl[jetVar+order+'_NoParton'].Draw('same')
    histColl[jetVar+order+'_PU'].Draw('same')
    canv.Print(outdir+jetVar+order+'_comb.pdf')
    canv.Print(outdir+jetVar+order+'_comb.png')

def printHists(canv, hists, outdir, option='hist'):
  for key,hist in hists.iteritems():
    canv.cd() 
    setXTitle(key,hist)
    if isinstance(hist,r.TH2):
      hist.Draw('colz,text')
    elif isinstance(hist,r.TH1):
      hist.SetMarkerStyle(r.kFullCircle)
      hist.SetMarkerColor(r.kBlack)
      hist.Draw(option)
    canv.Print(outdir+hist.GetName()+'.pdf')
    canv.Print(outdir+hist.GetName()+'.png')

def getEffSigma(theHist, wmin=115., wmax=135., step=0.001, epsilon=0.005):
  point = wmin
  weight = 0.
  points = [] #vector<pair<double,double> > 
  thesum = theHist.Integral()
  if thesum <= 0: return 0.
  for i in range(theHist.GetNbinsX()):
    weight += theHist.GetBinContent(i)
    if weight > epsilon:
      points.append( [theHist.GetBinCenter(i),weight/thesum] )
  low = wmin
  high = wmax
  width = wmax-wmin
  for i in range(len(points)):
    for j in range(i,len(points)):
      wy = points[j][1] - points[i][1]
      if abs(wy-0.683) < epsilon:
        wx = points[j][0] - points[i][0]
        if wx < width:
          low = points[i][0]
          high = points[j][0]
          width=wx
  return 0.5*(high-low)

def setXTitle(key,hist):
  if 'mgg' in hist.GetName():
    hist.GetXaxis().SetTitle('m_{#gamma#gamma}')

def deltaPhi( phi1, phi2):
  dPhi = phi1 - phi2;
  pi = 3.14159265;
  if   dPhi <= -1.*pi: dPhi += 2.*pi;
  elif dPhi > pi:      dPhi -= 2.*pi;
  return dPhi;

def deltaR( eta1, eta2, phi1, phi2 ):
  dPhi = deltaPhi(phi1,phi2)
  dEta = eta1 - eta2
  return sqrt( dPhi*dPhi + dEta*dEta )

#FIXME
class DiPhoton():
  pass