#!/usr/bin/env python
import ROOT 
from ROOT import TH1F 

ROOT.gStyle.SetOptStat(0)

tag="0.4_0.1_1"
samples=["decayAll"]#,"decay0","decay1","decay10"]
labels=["all taus"]#, "#pi", "#rho","a_{1} (3#pi)"]
colors=[ROOT.kBlack]# ,ROOT.kRed,ROOT.kBlue,ROOT.kGreen+2]

#variabs=["histoRecoTauType","histoRecoTauP","histoRecoTauMass","histoRecoTauTheta","histoGenTauP","histoGenTauType","histoGenTauVisP","histoGenTauTheta","histoGenTauVisMass"]
#xLabels=["Tau Type, Reco","Reco Tau P (GeV)","Reco Tau Mass (GeV)","Tau Theta","Gen Tau P","Gen Tau Type","Gen Tau Visible P (GeV)","Gen Tau Theta","Gen Tau Visible Mass"]

variabs=["histoRecoTauP"]
xLabels=["Reco Tau P (GeV)"]


files={}
for i in range(0,len(samples)):
  files[i]=ROOT.TFile("firstTest_"+samples[i]+"_"+tag+".root")

#Format for the rate histograms:
def formatHisto(file,variab,rename,titleX,color=ROOT.kBlack):
        histo = file.Get(variab)
        histo.SetName(rename)
        histo.SetXTitle(titleX)
        histo.SetLineColor(color)
        histo.SetLineWidth(2)
        #histo.SetMarkerColor(color)
        #histo.SetMarkerStyle(20)
        histo.Sumw2()
        return histo

# One variable per canvas 
iv=0
for var in variabs:
   c=ROOT.TCanvas("c"+var)
   leg=ROOT.TLegend(0.75,0.89,0.95,0.75)
   leg.SetFillStyle(0)
   leg.SetFillColor(0)
   leg.SetLineColor(0)

   histo={}
   for i in range(0,len(samples)):
     histo[i]=formatHisto(files[i],var,samples[i]+var,xLabels[iv],colors[i])
     leg.AddEntry(histo[i],labels[i],"l")

   # hack loop get the maximum right
   histo[0].Draw("hist")
   max=histo[0].GetMaximum()

   for i in range(1,len(samples)):
       histo[i].Draw("hist,sames")
       if histo[i].GetMaximum()>max:
            max=histo[i].GetMaximum()

   # some style tricks 
   #if "Mass" in var:
   #   histo[0].GetXaxis().SetRangeUser(0,3)

   #if "Type" in var:
   #   histo[0].GetXaxis().SetRangeUser(-1,20)
   #   histo[0].SetMinimum(1)

   histo[0].SetMaximum(max*1.4)
   histo[0].GetYaxis().SetTitle("Events (not normalized to xsec)")

   iv=iv+1
   leg.Draw()

   c.SaveAs(var+".png")
   #c.SetLogy()
   #c.SaveAs(var+"_LOG.png")
