import phase2tdrStyle
from ROOT import TH1F, TH2F, TLegend, gStyle, TGraph, kGreen, kBlue, gROOT, TGraphErrors, kBlack
gROOT.SetBatch(True)

doFudge = False
#doFudge = True

tgs = []
#tgs.append(TGraph())
tgs.append(TGraphErrors())
tgs[-1].SetName("tg%i"%len(tgs))

fileName = 'massHistSigmas_0pt8'
#fileName = 'massHistSigmas_1pt0'
#fileName = 'massHistSigmas_0pt8_debug'
for line in open(fileName+'_basic.txt','r'):
  x = float(line.strip().split(",")[0])
  y = float(line.strip().split(",")[1])
  e = float(line.strip().split(",")[2])
  if doFudge: 
    y *= 1.6/1.77
    e *= 1.6/1.77
  num = tgs[-1].GetN()
  tgs[-1].SetPoint(tgs[-1].GetN(),x,y)
  tgs[-1].SetPointError(tgs[-1].GetN()-1,0.,e)

can = phase2tdrStyle.setCanvas()
tgs[0].SetTitle('H#rightarrow#gamma#gamma, p_{T}^{#gamma} > 40 GeV, 1.6 < |#eta^{#gamma}| < 2.8')
tgs[0].SetLineColor(kGreen+2)
tgs[0].SetLineWidth(2)
tgs[0].SetMarkerStyle(20)
tgs[0].SetMarkerColor(kGreen+2)
#tgs[0].Draw("AL")
tgs[0].Draw("AP")
copyGraph = TGraph(tgs[0])
copyGraph.SetMarkerStyle(24)
copyGraph.SetMarkerColor(kBlack)
copyGraph.Draw("P")
#tgs[0].SetMarkerStyle(24)
#tgs[0].SetMarkerSize(0.5)
#tgs[0].SetMarkerColor(kGreen+2)
#tgs[0].Draw("P")
tgs[0].GetHistogram().GetXaxis().SetTitle("Stochastic term (%)")
tgs[0].GetHistogram().GetYaxis().SetTitle("Mass resolution (%)")
phase2tdrStyle.formatHisto(tgs[0].GetHistogram())
tgs[0].GetHistogram().GetXaxis().SetTitleOffset(1.3)
phase2tdrStyle.drawCMS()
phase2tdrStyle.drawEnPu()
leg = TLegend(0.45, 0.35-0.05*(len(tgs)+1), 0.9, 0.35)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetFillColor(0)
for h in tgs:
    leg.AddEntry(h,h.GetTitle(),'l')
#tgs.append(TGraph())
#tgs[-1].SetTitle("Run 2, 13TeV, inclusive WP")
#tgs[-1].SetPoint(0,0.297,0.57)
#tgs[-1].SetMarkerSize(2)
#tgs[-1].SetMarkerStyle(33)
#tgs[-1].Draw("p")
#leg.AddEntry(tgs[-1],tgs[-1].GetTitle(),'p')
leg.Draw()
#fileName = 'MassResoVsStochastic_'+fileName.split('_')[1]
fileName = 'MassResoVsStochastic_'+fileName.split('_')[1]+'_alt'
if doFudge: fileName += 'Fudged'
#FIXME
fileName += '_bugFixed'
can.SaveAs(fileName+'.pdf')
can.SaveAs(fileName+'.png')
