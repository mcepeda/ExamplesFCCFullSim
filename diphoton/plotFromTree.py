# Small (ugly) example of running over the tree, looking at gen and reco photons
# Many things to understand, this is just a fast check


import ROOT
from ROOT import TFile, TTree, TH1F, TH2F, TLorentzVector
import numpy as np
import math

file=TFile("miniTree.root")
tree=file.Get("outtree")

def dRAngle(p1,p2):
   dphi=p1.Phi()-p2.Phi()
   if (dphi>math.pi) : dphi=2*math.pi-dphi
   dtheta=p1.Theta()-p2.Theta()
   dR=math.sqrt(dtheta*dtheta+dphi*dphi)
   return dR

MAXDR=0.1

massPlot=ROOT.TH1F("massPi0Cands","Mass Pi0",60,0.01,0.3)
dRPlot=ROOT.TH1F("DRDiPhoton","",100,0,0.2)
massDRPlot=ROOT.TH2F("massDR","",60,0.01,0.3,100,0,0.2)

massGenPlot=ROOT.TH1F("massPi0CandsGen","Mass Pi0",60,0.01,0.3)
dRGenPlot=ROOT.TH1F("DRDiPhotonGen","",100,0,0.2)
massDRGenPlot=ROOT.TH2F("massDRGen","",60,0.01,0.3,100,0,0.2)

numberOfPi0s=ROOT.TH1F("Pi0sfound","",5,0,5)
numberOfPi0sGen=ROOT.TH1F("Pi0sfoundGen","",5,0,5)
numberOfPi0sGenResolved=ROOT.TH1F("Pi0sfoundGenResolved","",5,0,5)


# run over all events 
for event in tree:

    genPi0=event.nGenPi0s
    recoPi0=0
    genPi0Resolved=0

    # Gen Loop 
    for i in range(0,int(event.nGenPhotons)):

        v1=ROOT.TLorentzVector()
        v1.SetXYZM(event.genPhotonPx[i],event.genPhotonPy[i],event.genPhotonPz[i],0)

        best=-1
        bestDR=0.2
        #bestMassDiff=0.05

        for j in range(i+1,int(event.nGenPhotons)):

          vj=ROOT.TLorentzVector()
          vj.SetXYZM(event.genPhotonPx[j],event.genPhotonPy[j],event.genPhotonPz[j],0)
        
          dRij=dRAngle(v1,vj)
          #Pi0Cand=v1+vj

          if dRij<bestDR:
             best=j
             bestDR=dRij
          #if abs(Pi0.M()-0.13498)>bestMassDiff:
          #   best=j

        if best!=-1:
   
         v2=ROOT.TLorentzVector()
         v2.SetXYZM(event.genPhotonPx[best],event.genPhotonPy[best],event.genPhotonPz[best],0)
         dR12=dRAngle(v1,v2)
         Pi0=v1+v2

         dRGenPlot.Fill(dR12)
         massGenPlot.Fill(Pi0.M())
         massDRGenPlot.Fill(Pi0.M(),dR12)
 
         if dR12<0.1 and abs(Pi0.M()-0.13498)<0.02:
             genPi0Resolved+=1 


    # Reco Loop 
    for i in range(0,int(event.nPhotons)):

        #print (i,event.photonP[i])
        v1=ROOT.TLorentzVector()
        v1.SetXYZM(event.photonPx[i],event.photonPy[i],event.photonPz[i],0)

        best=-1
        bestDR=0.2
        #bestMassDiff=0.05

        for j in range(i+1,int(event.nPhotons)):

          vj=ROOT.TLorentzVector()
          vj.SetXYZM(event.photonPx[j],event.photonPy[j],event.photonPz[j],0)
        
          dRij=dRAngle(v1,vj)
          #Pi0Cand=v1+vj


          if dRij<bestDR:
             best=j
             bestDR=dRij
          #if abs(Pi0.M()-0.13498)>bestMassDiff:
          #   best=j

        if best!=-1:
   
         v2=ROOT.TLorentzVector()
         v2.SetXYZM(event.photonPx[best],event.photonPy[best],event.photonPz[best],0)
         dR12=dRAngle(v1,v2)
         Pi0=v1+v2

         dRPlot.Fill(dR12)
         massPlot.Fill(Pi0.M())
         massDRPlot.Fill(Pi0.M(),dR12)
 
         if dR12<0.1 and abs(Pi0.M()-0.13498)<0.02:
             recoPi0+=1 

    #print ("Pi0s? %d %d" %(genPi0,recoPi0))
    numberOfPi0s.Fill(recoPi0) 
    numberOfPi0sGen.Fill(genPi0)
    numberOfPi0sGenResolved.Fill(genPi0Resolved)
    


canvas=ROOT.TCanvas("c")
canvas.cd()
massPlot.Draw()
#massGenPlot.Draw()
#massPlot.Draw("sames")
#massGenPlot.SetLineColor(ROOT.kRed)
massPlot.SetXTitle("mass, GeV")
canvas.SaveAs("testMass.png")

#numberOfPi0s.Draw("")
#numberOfPi0sGen.Draw("sames")
#numberOfPi0sGen.SetLineColor(ROOT.kGreen+2)
#numberOfPi0sGenResolved.Draw("sames")
#numberOfPi0sGenResolved.SetLineColor(ROOT.kRed)
#numberOfPi0s.SetXTitle("number of Pi0s")
#canvas.SaveAs("numberOfPi0s.png")

massDRPlot.Draw("colz")
massDRPlot.SetXTitle("mass, GeV")
massDRPlot.SetYTitle("dR between photons")
canvas.SaveAs("testMassDR.png")

#massDRGenPlot.Draw("colz")
#massDRGenPlot.SetXTitle("mass, GeV")
#massDRGenPlot.SetYTitle("dR between photons")
#canvas.SaveAs("testMassDRGen.png")

dRGenPlot.Draw()
dRPlot.Draw("sames")
dRGenPlot.SetLineColor(ROOT.kRed)
dRPlot.SetXTitle("dR between photons (gen vs reco)")
canvas.SetLogy()
canvas.SaveAs("testDR.png")


