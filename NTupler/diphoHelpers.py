import ROOT as r

#helper functions
def initHistByEta(histColl, histName, bins=100, low=1., high=0.):
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

def printHists(canv, hists, outdir):
  for key,hist in hists.iteritems():
    canv.cd() 
    setXTitle(key,hist)
    hist.Draw("hist")
    canv.Print(outdir+hist.GetName()+".pdf")
    canv.Print(outdir+hist.GetName()+".png")

def getEffSigma(theHist, wmin=120., wmax=130., step=0.001, epsilon=0.005):
  point = wmin
  weight = 0.
  points = [] #vector<pair<double,double> > 
  thesum = theHist.Integral()
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
  hist.GetXaxis().SetTitle('m_{#gamma#gamma}')


#FIXME
class DiPhoton():
  pass
